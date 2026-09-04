module post_processor_class

    use kind_parameters
    use global_data
    use computational_domain_class
    use data_manager_class
    use field_pointers
    use post_processor_operation_class
    use boundary_conditions_class

    implicit none

    private
    public :: post_processor, post_processor_c

    type post_processor
        private
        type(post_processor_operation), dimension(:), allocatable :: post_processor_operations
        real(dp), dimension(:), allocatable :: values
        type(computational_domain) :: domain
        type(boundary_conditions_pointer) :: boundaries
        character(len=30) :: post_processor_output_file
        character(len=100) :: post_processor_title
        real(sp) :: save_time
        character(len=20) :: save_time_units
        real(sp) :: save_time_coefficient
        character(len=10) :: save_time_units_abbreviation
        integer :: operations_number
        integer :: operation_counter
        integer :: output_counter
    contains
        procedure, private :: read_properties
        procedure, private :: write_properties
        procedure, private :: set_properties
        procedure, private :: write_output_header
        procedure, private :: leading_indexes_to_coordinates

        procedure :: write_log
        procedure :: write_output_counter
        procedure :: create_post_processor_operation
        procedure :: process_data
    end type post_processor

    interface post_processor_c
        module procedure constructor
        module procedure constructor_file
    end interface

contains

    type(post_processor) function constructor( &
        manager, post_processor_output_file, post_processor_title, operations_number, &
        save_time, save_time_units, post_processor_data_file_name)

        type(data_manager), intent(in) :: manager
        character(len=*), intent(in) :: post_processor_output_file
        character(len=*), intent(in) :: post_processor_title
        integer, intent(in) :: operations_number
        real(dp), intent(in) :: save_time
        character(len=*), intent(in) :: save_time_units
        character(len=*), intent(in) :: post_processor_data_file_name

        integer :: output_counter = 0
        integer :: io_unit

        call constructor%set_properties( &
            manager, post_processor_output_file, post_processor_title, &
            operations_number, save_time, save_time_units, output_counter)

        open(newunit=io_unit, file=post_processor_data_file_name, &
            status='replace', form='formatted', delim='quote')
        call constructor%write_properties(io_unit)
        close(io_unit)
    end function constructor

    type(post_processor) function constructor_file(manager, post_processor_data_file_name)
        type(data_manager), intent(in) :: manager
        character(len=50), intent(in) :: post_processor_data_file_name

        integer :: io_unit
        integer :: operations_counter

        open(newunit=io_unit, file=post_processor_data_file_name, &
            status='old', form='formatted')
        call constructor_file%read_properties(manager, io_unit)

        do operations_counter = 1, size(constructor_file%post_processor_operations)
            constructor_file%post_processor_operations(operations_counter) = &
                post_processor_operation_c(manager, io_unit)
        end do
        close(io_unit)
    end function constructor_file

    subroutine read_properties(this, manager, pproc_unit)
        class(post_processor), intent(inout) :: this
        type(data_manager), intent(in) :: manager
        integer, intent(in) :: pproc_unit

        character(len=30) :: post_processor_output_file
        character(len=100) :: post_processor_title
        integer :: operations_number
        real(dp) :: save_time
        character(len=20) :: save_time_units
        integer :: output_counter
        integer :: ierr

        namelist /pproc_parameters/ post_processor_output_file, post_processor_title, &
            operations_number, save_time, save_time_units
        namelist /pproc_output_counter/ output_counter

        post_processor_title = 'NRG post-processor'

        read(unit=pproc_unit, nml=pproc_parameters)
        read(unit=pproc_unit, nml=pproc_output_counter, iostat=ierr)

        if (ierr /= 0) output_counter = 0
        rewind(pproc_unit)
        call this%set_properties( &
            manager, post_processor_output_file, post_processor_title, &
            operations_number, save_time, save_time_units, output_counter)
    end subroutine read_properties

    subroutine write_properties(this, pproc_unit)
        class(post_processor), intent(in) :: this
        integer, intent(in) :: pproc_unit

        character(len=30) :: post_processor_output_file
        character(len=100) :: post_processor_title
        integer :: operations_number
        real(dp) :: save_time
        character(len=20) :: save_time_units

        namelist /pproc_parameters/ post_processor_output_file, post_processor_title, &
            operations_number, save_time, save_time_units

        post_processor_output_file = this%post_processor_output_file
        post_processor_title = this%post_processor_title
        operations_number = size(this%post_processor_operations)
        save_time = this%save_time
        save_time_units = this%save_time_units

        write(unit=pproc_unit, nml=pproc_parameters)
    end subroutine write_properties

    subroutine set_properties( &
        this, manager, post_processor_output_file, post_processor_title, &
        operations_number, save_time, save_time_units, output_counter)

        class(post_processor), intent(inout) :: this
        character(len=*), intent(in) :: post_processor_output_file
        character(len=*), intent(in) :: post_processor_title
        type(data_manager), intent(in) :: manager
        integer, intent(in) :: operations_number
        real(dp), intent(in) :: save_time
        character(len=*), intent(in) :: save_time_units
        integer, intent(in) :: output_counter

        this%domain = manager%domain
        this%boundaries%bc_ptr => manager%boundary_conditions_pointer%bc_ptr

        allocate(this%post_processor_operations(operations_number))
        allocate(this%values(operations_number + 1 + this%domain%get_domain_dimensions()))

        if (len_trim(post_processor_title) == 0) then
            error stop 'Postprocessor: post_processor_title must not be empty'
        end if
        if (index(post_processor_title, '"') /= 0) then
            error stop 'Postprocessor: post_processor_title must not contain quotation marks'
        end if

        this%post_processor_output_file = trim(post_processor_output_file)
        this%post_processor_title = trim(post_processor_title)
        this%save_time = save_time
        this%save_time_units = trim(save_time_units)
        this%output_counter = output_counter
        this%operations_number = operations_number
        this%operation_counter = 0

        select case(trim(this%save_time_units))
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
            error stop 'Postprocessor: unsupported save_time_units'
        end select
    end subroutine set_properties

    subroutine write_log(this, log_unit)
        class(post_processor), intent(in) :: this
        integer, intent(in) :: log_unit
        integer :: operations_counter

        write(log_unit,'(A,A)') ' Title: ', trim(this%post_processor_title)
        write(log_unit,'(A,E14.7,A)') &
            ' Write time: ', this%save_time, this%save_time_units_abbreviation
        write(log_unit,'(A)') ' Operations: '

        do operations_counter = 1, size(this%post_processor_operations)
            write(log_unit,'(A,I2,A)') achar(9)//'#', operations_counter, ':'
            call this%post_processor_operations(operations_counter)%write_log(log_unit)
        end do
    end subroutine write_log

    subroutine write_output_counter(this, post_processor_output_file)
        class(post_processor), intent(in) :: this
        character(len=*), intent(in) :: post_processor_output_file

        integer :: io_unit
        integer :: ierr
        integer :: output_counter

        namelist /pproc_output_counter/ output_counter

        open(newunit=io_unit, file=post_processor_output_file, &
            status='old', form='formatted')
        read(io_unit, nml=pproc_output_counter, iostat=ierr)

        if (ierr == 0) then
            backspace(io_unit)
            backspace(io_unit)
            backspace(io_unit)
        end if

        output_counter = this%output_counter
        write(unit=io_unit, nml=pproc_output_counter)
        close(io_unit)
    end subroutine write_output_counter

    subroutine create_post_processor_operation( &
        this, manager, field_name, operation_type, operation_area, &
        operation_area_distance, grad_projection, post_processor_data_file_name, &
        location_mode)

        class(post_processor), intent(inout) :: this
        type(data_manager), intent(in) :: manager
        character(len=*), intent(in) :: field_name
        character(len=*), intent(in) :: operation_type
        integer, dimension(3,2), intent(in) :: operation_area
        integer, dimension(3), intent(in) :: operation_area_distance
        integer, intent(in) :: grad_projection
        character(len=50), intent(in) :: post_processor_data_file_name
        character(len=*), intent(in), optional :: location_mode

        integer :: io_unit

        this%operation_counter = this%operation_counter + 1

        if (this%operation_counter > this%operations_number) then
            error stop &
                'Postprocessor: trying to create more operations than operations_number'
        end if

        open(newunit=io_unit, file=post_processor_data_file_name, status='old', &
            form='formatted', position='append', delim='quote')

        if (present(location_mode)) then
            this%post_processor_operations(this%operation_counter) = &
                post_processor_operation_c( &
                    manager=manager, &
                    field_name=field_name, &
                    operation_type=operation_type, &
                    operation_area=operation_area, &
                    operation_area_distance=operation_area_distance, &
                    grad_projection=grad_projection, &
                    post_processor_data_file_unit=io_unit, &
                    location_mode=location_mode)
        else
            this%post_processor_operations(this%operation_counter) = &
                post_processor_operation_c( &
                    manager=manager, &
                    field_name=field_name, &
                    operation_type=operation_type, &
                    operation_area=operation_area, &
                    operation_area_distance=operation_area_distance, &
                    grad_projection=grad_projection, &
                    post_processor_data_file_unit=io_unit)
        end if

        close(io_unit)
    end subroutine create_post_processor_operation

    subroutine write_output_header(this, output_file_unit)
        class(post_processor), intent(in) :: this
        integer, intent(in) :: output_file_unit

        character(len=5), dimension(3) :: axis_names
        character(len=20) :: coordinate_system
        character(len=10) :: coordinate_units
        character(len=128) :: variable_name
        integer :: dimensions
        integer :: dim
        integer :: operations_counter

        dimensions = this%domain%get_domain_dimensions()
        axis_names = this%domain%get_axis_names()
        coordinate_system = this%domain%get_coordinate_system_name()

        write(output_file_unit,'(A)') &
            'TITLE="' // trim(this%post_processor_title) // '"'

        write(output_file_unit,'(A)',advance='no') 'VARIABLES='
        write(output_file_unit,'(A)',advance='no') &
            '"time [' // trim(this%save_time_units_abbreviation) // ']"'

        variable_name = this%post_processor_operations(1)%get_output_variable_name()
        write(output_file_unit,'(A)',advance='no') &
            ' "' // trim(variable_name) // '"'

        do dim = 1, dimensions
            coordinate_units = 'm'
            select case(trim(coordinate_system))
            case('cylindrical')
                if (dim == 2) coordinate_units = 'rad'
            case('spherical')
                if (dim >= 2) coordinate_units = 'rad'
            end select

            write(output_file_unit,'(A)',advance='no') &
                ' "' // trim(axis_names(dim)) // ' [' // trim(coordinate_units) // ']"'
        end do

        do operations_counter = 2, size(this%post_processor_operations)
            variable_name = &
                this%post_processor_operations(operations_counter)%get_output_variable_name()
            write(output_file_unit,'(A)',advance='no') &
                ' "' // trim(variable_name) // '"'
        end do

        write(output_file_unit,*)
    end subroutine write_output_header

    subroutine leading_indexes_to_coordinates(this, indexes, coordinates)
        class(post_processor), intent(in) :: this
        real(dp), dimension(3), intent(in) :: indexes
        real(dp), dimension(3), intent(out) :: coordinates

        real(dp), dimension(:,:), allocatable :: domain_lengths
        integer, dimension(3) :: cells_number
        real(dp) :: cell_size
        integer :: dimensions
        integer :: dim

        dimensions = this%domain%get_domain_dimensions()
        domain_lengths = this%domain%get_domain_lengths()
        cells_number = this%domain%get_global_cells_number()
        cells_number(1:dimensions) = cells_number(1:dimensions) - 2

        coordinates = 0.0_dp

        do dim = 1, dimensions
            if (indexes(dim) <= 0.0_dp) cycle

            cell_size = &
                (domain_lengths(dim,2) - domain_lengths(dim,1)) / &
                real(cells_number(dim), dp)

            coordinates(dim) = &
                domain_lengths(dim,1) + (indexes(dim) - 0.5_dp)*cell_size
        end do
    end subroutine leading_indexes_to_coordinates

    subroutine process_data(this, time)
        class(post_processor), intent(inout) :: this
        real(dp), intent(in) :: time

        integer :: output_file_unit
        integer, dimension(3) :: leading_point_indexes, point_indexes
        real(dp), dimension(3) :: leading_point_indexes_real
        real(dp), dimension(3) :: leading_point_coordinates
        real(dp) :: current_output_time
        real(dp) :: target_output_time
        real(dp) :: save_interval
        real(dp) :: schedule_tolerance

        logical :: file_exists
        logical :: write_header
        integer :: output_file_size
        integer :: i, dim, operations_counter
        integer :: processor_rank

        processor_rank = this%domain%get_processor_rank()

        save_interval       = real(this%save_time, dp)
        current_output_time = time*real(this%save_time_coefficient, dp)
        target_output_time  = save_interval*real(this%output_counter, dp)

        schedule_tolerance = 100.0_dp*epsilon(1.0_dp)* &
            max(1.0_dp, abs(target_output_time))

        if ((current_output_time + schedule_tolerance >= target_output_time) .or. &
            (time == 0.0_dp)) then

            if (processor_rank == 0) then
                output_file_size = 0
                inquire(file=this%post_processor_output_file, exist=file_exists, &
                    size=output_file_size)

                write_header = (.not. file_exists)
                if (file_exists) write_header = (output_file_size == 0)

                if (file_exists) then
                    open(newunit=output_file_unit, &
                        file=this%post_processor_output_file, status='old', &
                        form='formatted', position='append')
                else
                    open(newunit=output_file_unit, &
                        file=this%post_processor_output_file, status='new', &
                        form='formatted')
                end if

                if (write_header) call this%write_output_header(output_file_unit)
            end if

            this%values(1) = current_output_time

            call this%post_processor_operations(1)%process_operation( &
                this%boundaries, leading_point_indexes, this%values(2), &
                leading_point_indexes_real)

            call this%leading_indexes_to_coordinates( &
                leading_point_indexes_real, leading_point_coordinates)

            do dim = 1, this%domain%get_domain_dimensions()
                this%values(2 + dim) = leading_point_coordinates(dim)
            end do

            do operations_counter = 2, size(this%post_processor_operations)
                call this%post_processor_operations(operations_counter)%set_point( &
                    leading_point_indexes)

                if (associated( &
                    this%post_processor_operations(operations_counter)% &
                    operation_field_scal%s_ptr)) then

                    call this%post_processor_operations(operations_counter)% &
                        process_operation( &
                            this%boundaries, point_indexes, &
                            this%values(operations_counter + &
                            this%domain%get_domain_dimensions() + 1))
                end if

                if (associated( &
                    this%post_processor_operations(operations_counter)% &
                    operation_field_vect%v_ptr)) then

                    call this%post_processor_operations(operations_counter)% &
                        process_operation( &
                            this%boundaries, point_indexes, &
                            this%values(operations_counter + &
                            this%domain%get_domain_dimensions() + 1))
                end if
            end do

            if (processor_rank == 0) then
                write(output_file_unit,'(*(E20.12))') &
                    (this%values(i), i=1,size(this%values))
                close(output_file_unit)
            end if

            this%output_counter = max( &
                this%output_counter + 1, &
                floor(current_output_time/save_interval) + 1)
        end if
    end subroutine process_data

end module post_processor_class
