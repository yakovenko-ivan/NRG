module data_save_class

    use kind_parameters
    use global_data
    use computational_domain_class
    use computational_mesh_class
    use boundary_conditions_class
    use field_pointers
    use data_manager_class
    use field_output_types
    use field_writer_class
    use tecplot_writer_class
    use cgns_writer_class

    implicit none

    private
    public :: data_save, data_save_c

    type :: data_save
        private
        type(computational_domain) :: domain
        type(computational_mesh_pointer) :: mesh
        type(boundary_conditions_pointer) :: boundaries

        type(field_scalar_cons_pointer), dimension(:), allocatable :: visible_fields
        character(len=40), dimension(:), allocatable :: visible_fields_names

        real(dp) :: save_time
        character(len=20) :: save_time_units
        real(sp) :: save_time_coefficient
        character(len=25) :: save_format
        character(len=10) :: save_time_units_abbreviation
        character(len=32) :: dataset_name

        integer :: output_counter
        character(len=256) :: data_save_folder
        logical :: debug_flag

        class(field_writer), allocatable :: writer
    contains
        procedure, private :: read_properties
        procedure, private :: write_properties
        procedure, private :: set_properties
        procedure, private :: initialize_writer
        procedure, private :: write_output_counter
        procedure, private :: build_snapshot

        procedure :: save_all_data

        procedure :: get_data_save_folder
        procedure :: get_save_time
        procedure :: get_save_time_units_abbreviation
        procedure :: get_save_time_coefficient
        procedure :: get_visible_fields_number

        procedure :: set_save_time
        procedure :: write_log
    end type data_save

    interface data_save_c
        module procedure constructor
        module procedure constructor_file
    end interface

contains

    type(data_save) function constructor(manager, visible_fields_names, save_time, &
            save_time_units, save_format, data_save_folder, debug_flag, dataset_name)
        type(data_manager), intent(in) :: manager
        character(len=40), dimension(:), intent(in) :: visible_fields_names
        real(dp), intent(in) :: save_time
        character(len=*), intent(in) :: save_time_units
        character(len=*), intent(in) :: save_format
        character(len=*), intent(in) :: data_save_folder
        logical, intent(in) :: debug_flag
        character(len=*), intent(in), optional :: dataset_name

        integer :: output_counter, io_unit
        character(len=32) :: dataset_name_local

        output_counter = 0
        dataset_name_local = 'NRG'
        if (present(dataset_name)) then
            if (len_trim(dataset_name) == 0) then
                error stop 'Data save: dataset_name must not be empty'
            end if
            if (len_trim(dataset_name) > len(dataset_name_local)) then
                error stop 'Data save: dataset_name exceeds the 32-character CGNS name limit'
            end if
            dataset_name_local = trim(adjustl(dataset_name))
        end if

        call constructor%set_properties(manager, visible_fields_names, save_time, &
            save_time_units, save_format, data_save_folder, debug_flag, output_counter, &
            dataset_name_local)

        open(newunit=io_unit, file=data_save_data_file_name, status='replace', &
            form='formatted', delim='quote')
        call constructor%write_properties(io_unit)
        close(io_unit)
    end function constructor

    type(data_save) function constructor_file(manager)
        type(data_manager), intent(in) :: manager
        integer :: io_unit

        open(newunit=io_unit, file=data_save_data_file_name, status='old', &
            form='formatted', delim='quote')
        call constructor_file%read_properties(manager, io_unit)
        close(io_unit)
    end function constructor_file

    subroutine read_properties(this, manager, data_save_file_unit)
        class(data_save), intent(inout) :: this
        type(data_manager), intent(in) :: manager
        integer, intent(in) :: data_save_file_unit

        character(len=40), dimension(:), allocatable :: visible_fields_names
        integer :: visible_fields_number
        real(dp) :: save_time
        character(len=20) :: save_time_units
        character(len=25) :: save_format
        character(len=256) :: data_save_folder
        logical :: debug_flag
        character(len=32) :: dataset_name
        integer :: output_counter, ierr

        namelist /data_save_parameters/ visible_fields_number, save_time, &
            save_time_units, save_format, data_save_folder, debug_flag, dataset_name
        namelist /visible_fields/ visible_fields_names
        namelist /data_save_output_counter/ output_counter

        ! Keep older data_save namelists readable: if dataset_name is absent,
        ! namelist input leaves this default value unchanged.
        dataset_name = 'NRG'
        read(unit=data_save_file_unit, nml=data_save_parameters)
        allocate(visible_fields_names(visible_fields_number))
        read(unit=data_save_file_unit, nml=visible_fields)
        read(unit=data_save_file_unit, nml=data_save_output_counter, iostat=ierr)

        if (ierr /= 0) output_counter = 0
        call this%set_properties(manager, visible_fields_names, save_time, &
            save_time_units, save_format, data_save_folder, debug_flag, output_counter, &
            dataset_name)
    end subroutine read_properties

    subroutine write_properties(this, data_save_file_unit)
        class(data_save), intent(in) :: this
        integer, intent(in) :: data_save_file_unit

        character(len=40), dimension(:), allocatable :: visible_fields_names
        integer :: visible_fields_number
        real(dp) :: save_time
        character(len=20) :: save_time_units
        character(len=25) :: save_format
        character(len=256) :: data_save_folder
        logical :: debug_flag
        character(len=32) :: dataset_name

        namelist /data_save_parameters/ visible_fields_number, save_time, &
            save_time_units, save_format, data_save_folder, debug_flag, dataset_name
        namelist /visible_fields/ visible_fields_names

        visible_fields_number = size(this%visible_fields_names)
        allocate(visible_fields_names(visible_fields_number))
        visible_fields_names = this%visible_fields_names

        save_time = this%save_time
        save_time_units = this%save_time_units
        save_format = this%save_format
        data_save_folder = this%data_save_folder
        debug_flag = this%debug_flag
        dataset_name = this%dataset_name

        write(unit=data_save_file_unit, nml=data_save_parameters)
        write(unit=data_save_file_unit, nml=visible_fields)
    end subroutine write_properties

    subroutine set_properties(this, manager, visible_fields_names, save_time, &
            save_time_units, save_format, data_save_folder, debug_flag, output_counter, &
            dataset_name)
        class(data_save), intent(inout) :: this
        type(data_manager), intent(in) :: manager
        character(len=*), dimension(:), intent(in) :: visible_fields_names
        real(dp), intent(in) :: save_time
        character(len=*), intent(in) :: save_time_units
        character(len=*), intent(in) :: save_format
        character(len=*), intent(in) :: data_save_folder
        logical, intent(in) :: debug_flag
        integer, intent(in) :: output_counter
        character(len=*), intent(in) :: dataset_name

        integer :: vector_projections_number
        integer :: number_of_cons_scalar_fields, number_of_cons_vector_fields
        integer :: processor_rank
        integer :: number_of_visible_fields
        integer :: list_count, scal_count, vect_count, proj_count
        type(field_scalar_cons_pointer) :: scal_ptr
        type(field_vector_cons_pointer) :: vect_ptr
        type(field_tensor_cons_pointer) :: tens_ptr

        this%domain = manager%domain
        this%save_time = save_time
        this%save_time_units = trim(save_time_units)
        this%save_format = trim(adjustl(save_format))
        this%data_save_folder = trim(data_save_folder)
        this%output_counter = output_counter
        this%debug_flag = debug_flag

        if (len_trim(dataset_name) == 0) then
            error stop 'Data save: dataset_name must not be empty'
        end if
        if (len_trim(dataset_name) > len(this%dataset_name)) then
            error stop 'Data save: dataset_name exceeds the 32-character CGNS name limit'
        end if
        this%dataset_name = trim(adjustl(dataset_name))

        select case(trim(this%save_time_units))
            case('minutes')
                this%save_time_units_abbreviation = 'm'
                this%save_time_coefficient = 1.0_dp/60.0_dp
            case('seconds')
                this%save_time_units_abbreviation = 's'
                this%save_time_coefficient = 1.0_dp
            case('milliseconds')
                this%save_time_units_abbreviation = 'ms'
                this%save_time_coefficient = 1.0e3_dp
            case('microseconds')
                this%save_time_units_abbreviation = 'us'
                this%save_time_coefficient = 1.0e6_dp
            case('nanoseconds')
                this%save_time_units_abbreviation = 'ns'
                this%save_time_coefficient = 1.0e9_dp
            case default
                error stop 'Data save: unsupported save_time_units'
        end select

        this%mesh = manager%computational_mesh_pointer
        this%boundaries = manager%boundary_conditions_pointer

        allocate(this%visible_fields_names(size(visible_fields_names)))
        do list_count = 1, size(visible_fields_names)
            this%visible_fields_names(list_count) = visible_fields_names(list_count)
        end do

        number_of_cons_scalar_fields = manager%get_number_of_cons_scalar_fields()
        number_of_cons_vector_fields = manager%get_number_of_cons_vector_fields()

        number_of_visible_fields = 0
        if (this%debug_flag) then
            number_of_visible_fields = number_of_cons_scalar_fields
            do vect_count = 1, number_of_cons_vector_fields
                call manager%get_cons_field_pointer_by_number( &
                    scal_ptr, vect_ptr, tens_ptr, 'vector', vect_count)
                number_of_visible_fields = number_of_visible_fields + &
                    vect_ptr%v_ptr%get_projections_number()
            end do
        else
            do list_count = 1, size(visible_fields_names)
                call manager%get_cons_field_pointer_by_name( &
                    scal_ptr, vect_ptr, tens_ptr, visible_fields_names(list_count))
                if (associated(scal_ptr%s_ptr)) number_of_visible_fields = number_of_visible_fields + 1
                if (associated(vect_ptr%v_ptr)) then
                    number_of_visible_fields = number_of_visible_fields + &
                        vect_ptr%v_ptr%get_projections_number()
                end if
            end do
        end if

        allocate(this%visible_fields(number_of_visible_fields))
        number_of_visible_fields = 0

        if (this%debug_flag) then
            do scal_count = 1, number_of_cons_scalar_fields
                number_of_visible_fields = number_of_visible_fields + 1
                call manager%get_cons_field_pointer_by_number( &
                    scal_ptr, vect_ptr, tens_ptr, 'scalar', scal_count)
                this%visible_fields(number_of_visible_fields)%s_ptr => scal_ptr%s_ptr
            end do
            do vect_count = 1, number_of_cons_vector_fields
                call manager%get_cons_field_pointer_by_number( &
                    scal_ptr, vect_ptr, tens_ptr, 'vector', vect_count)
                vector_projections_number = vect_ptr%v_ptr%get_projections_number()
                do proj_count = 1, vector_projections_number
                    number_of_visible_fields = number_of_visible_fields + 1
                    this%visible_fields(number_of_visible_fields)%s_ptr => &
                        vect_ptr%v_ptr%pr(proj_count)
                end do
            end do
        else
            do list_count = 1, size(visible_fields_names)
                call manager%get_cons_field_pointer_by_name( &
                    scal_ptr, vect_ptr, tens_ptr, visible_fields_names(list_count))
                if (associated(scal_ptr%s_ptr)) then
                    number_of_visible_fields = number_of_visible_fields + 1
                    this%visible_fields(number_of_visible_fields)%s_ptr => scal_ptr%s_ptr
                end if
                if (associated(vect_ptr%v_ptr)) then
                    vector_projections_number = vect_ptr%v_ptr%get_projections_number()
                    do proj_count = 1, vector_projections_number
                        number_of_visible_fields = number_of_visible_fields + 1
                        this%visible_fields(number_of_visible_fields)%s_ptr => &
                            vect_ptr%v_ptr%pr(proj_count)
                    end do
                end if
            end do
        end if

        call this%initialize_writer()

        processor_rank = this%domain%get_processor_rank()
        if (processor_rank == 0) call make_directory(this%data_save_folder)
    end subroutine set_properties

    subroutine initialize_writer(this)
        class(data_save), intent(inout) :: this

        if (allocated(this%writer)) deallocate(this%writer)

        select case(trim(this%save_format))
            case('tecplot')
                allocate(tecplot_writer :: this%writer)
            case('cgns', 'cgns_hdf5', 'cgns/hdf5')
                if (.not. cgns_backend_enabled) then
                    error stop 'Data save: CGNS requested but NRG was built with NRG_ENABLE_CGNS=OFF'
                end if
                allocate(cgns_writer :: this%writer)
            case default
                write(*,'(A)') 'Unknown data save format: '//trim(this%save_format)
                error stop 'Data save: unsupported output format'
        end select
    end subroutine initialize_writer

    subroutine build_snapshot(this, time_seconds, display_time, debug, snapshot)
        class(data_save), intent(in) :: this
        real(dp), intent(in) :: time_seconds, display_time
        logical, intent(in) :: debug
        type(field_snapshot), intent(out) :: snapshot

        integer :: field_index
        character(len=20) :: time_text

        snapshot%domain = this%domain
        snapshot%mesh = this%mesh
        snapshot%boundaries = this%boundaries
        snapshot%time_seconds = time_seconds
        snapshot%display_time = display_time
        snapshot%output_index = this%output_counter
        snapshot%include_ghost_cells = debug
        snapshot%time_units_abbreviation = this%save_time_units_abbreviation
        snapshot%dataset_name = this%dataset_name

        ! Match the legacy Tecplot scientific-notation time label, but remove
        ! both the old "zone" prefix and the format field's leading blank.
        write(time_text,'(E11.4)') display_time
        snapshot%zone_name = trim(adjustl(time_text)) // &
            trim(this%save_time_units_abbreviation)

        allocate(snapshot%fields(size(this%visible_fields)))
        do field_index = 1, size(this%visible_fields)
            snapshot%fields(field_index)%values%s_ptr => this%visible_fields(field_index)%s_ptr
            snapshot%fields(field_index)%nrg_name = this%visible_fields(field_index)%s_ptr%name_long
            snapshot%fields(field_index)%short_name = this%visible_fields(field_index)%s_ptr%name_short
            snapshot%fields(field_index)%grid_location = FIELD_LOCATION_CELL_CENTER
        end do
    end subroutine build_snapshot

    subroutine save_all_data(this, time, stop_flag, make_save)
        class(data_save), intent(inout) :: this
        real(dp), intent(in) :: time
        logical, intent(in), optional :: make_save
        logical, intent(in) :: stop_flag

        real(dp) :: display_time
        character(len=100) :: file_stem
        logical :: debug, make_flag
        integer :: processor_rank
        type(field_snapshot) :: snapshot

        processor_rank = this%domain%get_processor_rank()

        make_flag = .false.
        if (stop_flag) make_flag = .true.
        if (present(make_save)) make_flag = make_save

        if ((time*this%save_time_coefficient >= &
             this%save_time*(this%output_counter+1)) .or. make_flag) then

            display_time = time*this%save_time_coefficient
            debug = make_flag .or. this%debug_flag
            write(file_stem,'(I6.6,A)') int(display_time), &
                trim(this%save_time_units_abbreviation)

            ! Preserve legacy debug naming/time semantics for the Tecplot backend.
            if (this%debug_flag) then
                display_time = real(this%output_counter, dp)
                write(file_stem,'(I6.6)') this%output_counter
            end if

            call this%build_snapshot(time, display_time, debug, snapshot)
            call this%writer%write_snapshot(snapshot, this%data_save_folder, file_stem)

            this%output_counter = this%output_counter + 1
        end if

        if (stop_flag .and. processor_rank == 0) call this%write_output_counter()
    end subroutine save_all_data

    subroutine write_output_counter(this)
        class(data_save), intent(in) :: this
        integer :: output_counter, io_unit, ierr
        namelist /data_save_output_counter/ output_counter

        open(newunit=io_unit, file=data_save_data_file_name, status='old', form='formatted')
        read(unit=io_unit, nml=data_save_output_counter, iostat=ierr)

        if (ierr == 0) then
            backspace(io_unit)
            backspace(io_unit)
            backspace(io_unit)
        end if

        output_counter = this%output_counter
        write(unit=io_unit, nml=data_save_output_counter)
        close(io_unit)
    end subroutine write_output_counter

    pure function get_data_save_folder(this) result(folder)
        class(data_save), intent(in) :: this
        character(len=256) :: folder
        folder = this%data_save_folder
    end function get_data_save_folder

    pure function get_save_time(this) result(value)
        class(data_save), intent(in) :: this
        real(dp) :: value
        value = this%save_time
    end function get_save_time

    pure function get_save_time_units_abbreviation(this) result(units)
        class(data_save), intent(in) :: this
        character(len=10) :: units
        units = this%save_time_units_abbreviation
    end function get_save_time_units_abbreviation

    pure function get_save_time_coefficient(this) result(value)
        class(data_save), intent(in) :: this
        real(sp) :: value
        value = this%save_time_coefficient
    end function get_save_time_coefficient

    pure function get_visible_fields_number(this) result(number)
        class(data_save), intent(in) :: this
        integer :: number
        number = size(this%visible_fields)
    end function get_visible_fields_number

    pure subroutine set_save_time(this, new_save_time)
        class(data_save), intent(inout) :: this
        real(dp), intent(in) :: new_save_time

        this%output_counter = this%output_counter * int(this%save_time/new_save_time)
        this%save_time = new_save_time
    end subroutine set_save_time

    subroutine write_log(this, log_unit)
        class(data_save), intent(in) :: this
        integer, intent(in) :: log_unit

        write(log_unit,'(A)') '************************************************************************************* '
        write(log_unit,'(A)') ' Data save setup:'
        write(log_unit,'(A,E14.7,A)') ' Data save time: ', this%save_time, &
            this%save_time_units_abbreviation
        write(log_unit,'(A,A)') ' Data set name: ', trim(this%dataset_name)
        write(log_unit,'(A,A)') ' Data save format: ', trim(this%save_format)
        write(log_unit,'(A,A)') ' Data save folder: ', trim(this%data_save_folder)
        write(log_unit,'(A)') '************************************************************************************* '
    end subroutine write_log

    subroutine make_directory(path)
        character(len=*), intent(in) :: path
        integer :: exitstat

#ifdef WIN
        call execute_command_line('mkdir "'//trim(path)//'"', &
            wait=.true., exitstat=exitstat)
#else
        call execute_command_line('mkdir -p "'//trim(path)//'"', &
            wait=.true., exitstat=exitstat)
#endif
        if (exitstat /= 0) then
            write(*,'(A)') 'WARNING: unable to create output directory: '//trim(path)
        end if
    end subroutine make_directory

end module data_save_class
