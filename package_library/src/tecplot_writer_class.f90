module tecplot_writer_class

    use kind_parameters
    use global_data
    use field_writer_class
    use field_output_types

#ifdef mpi
    use MPI
#endif

    implicit none

    private
    public :: tecplot_writer

    type, extends(field_writer) :: tecplot_writer
        real(sp), dimension(:,:,:), allocatable :: io_buffer
    contains
        procedure :: write_snapshot => tecplot_write_snapshot
        procedure, private :: ensure_buffer
        procedure, private :: generate_header
        procedure, private :: save_mesh
        procedure, private :: save_fields
        procedure, private :: save_bounds
    end type tecplot_writer

contains

    subroutine tecplot_write_snapshot(this, snapshot, output_folder, file_stem)
        class(tecplot_writer), intent(inout) :: this
        type(field_snapshot), intent(in) :: snapshot
        character(len=*), intent(in) :: output_folder
        character(len=*), intent(in) :: file_stem

        character(len=256) :: file_name, proc_rank, snapshot_folder
        integer :: unit_io, processor_rank
#ifdef mpi
        integer :: error
#endif

        processor_rank = snapshot%domain%get_processor_rank()
        call this%ensure_buffer(snapshot)

#ifdef mpi
        snapshot_folder = trim(output_folder)//trim(fold_sep)//trim(file_stem)

        if (processor_rank == 0) then
            call make_directory(snapshot_folder)
            file_name = trim(snapshot_folder)//trim(fold_sep)//trim(file_stem)//'_header.plt'
            open(newunit=unit_io, file=trim(file_name), status='replace', &
                 access='stream', form='unformatted')
            call this%generate_header(snapshot, unit_io)
            close(unit_io)
        end if

        call MPI_BARRIER(MPI_COMM_WORLD, error)

        write(proc_rank,'(A,I4.4)') '_proc_', processor_rank
        file_name = trim(snapshot_folder)//trim(fold_sep)//trim(file_stem)// &
                    trim(proc_rank)//'_data.plt'
        open(newunit=unit_io, file=trim(file_name), status='replace', &
             access='stream', form='unformatted')
        call this%save_mesh(snapshot, unit_io)
        call this%save_fields(snapshot, unit_io)
        if (snapshot%include_ghost_cells) call this%save_bounds(snapshot, unit_io)
        close(unit_io)
#else
        file_name = trim(output_folder)//trim(fold_sep)//trim(file_stem)//'.plt'
        open(newunit=unit_io, file=trim(file_name), status='replace', &
             access='stream', form='unformatted')
        call this%generate_header(snapshot, unit_io)
        call this%save_mesh(snapshot, unit_io)
        call this%save_fields(snapshot, unit_io)
        if (snapshot%include_ghost_cells) call this%save_bounds(snapshot, unit_io)
        close(unit_io)
#endif
    end subroutine tecplot_write_snapshot

    subroutine ensure_buffer(this, snapshot)
        class(tecplot_writer), intent(inout) :: this
        type(field_snapshot), intent(in) :: snapshot
        integer, dimension(3,2) :: allocation_bounds

        if (allocated(this%io_buffer)) return

        allocation_bounds = snapshot%domain%get_local_utter_cells_bounds()
        allocate(this%io_buffer( &
            allocation_bounds(1,1):allocation_bounds(1,2), &
            allocation_bounds(2,1):allocation_bounds(2,2), &
            allocation_bounds(3,1):allocation_bounds(3,2)))
    end subroutine ensure_buffer

    subroutine generate_header(this, snapshot, tecplot_io_unit)
        class(tecplot_writer), intent(inout) :: this
        type(field_snapshot), intent(in) :: snapshot
        integer, intent(in) :: tecplot_io_unit

        integer :: file_type
        integer, dimension(60) :: converted_string
        character(len=40) :: zone_name
        integer :: domain_dimensions
        integer, dimension(3) :: utter_cells_number, inner_cells_number
        character(len=5), dimension(3) :: axis_names
        real(dp), dimension(:,:), allocatable :: domain_lengths
        real(dp) :: min_value, max_value
        integer :: variables_number, variables_counter
        integer :: dim, counter

        file_type = 1
        domain_dimensions = snapshot%domain%get_domain_dimensions()
        utter_cells_number = snapshot%domain%get_global_cells_number()
        inner_cells_number = 1
        inner_cells_number(:domain_dimensions) = &
            utter_cells_number(:domain_dimensions) - 2

        axis_names = snapshot%domain%get_axis_names()
        domain_lengths = snapshot%domain%get_domain_lengths()

        variables_number = domain_dimensions + size(snapshot%fields)
        if (snapshot%include_ghost_cells) variables_number = variables_number + 1

        ! Preserve the legacy TDV112 byte layout exactly.
        write(tecplot_io_unit) '#!TDV112'
        write(tecplot_io_unit) 1
        write(tecplot_io_unit) file_type

        write(tecplot_io_unit) (converted_string(counter), counter=1, &
            convert_string(converted_string, snapshot%dataset_name))

        write(tecplot_io_unit) variables_number

        do dim = 1, domain_dimensions
            write(tecplot_io_unit) (converted_string(counter), counter=1, &
                convert_string(converted_string, axis_names(dim)))
        end do
        do variables_counter = 1, size(snapshot%fields)
            write(tecplot_io_unit) (converted_string(counter), counter=1, &
                convert_string(converted_string, &
                snapshot%fields(variables_counter)%short_name))
        end do
        if (snapshot%include_ghost_cells) then
            write(tecplot_io_unit) (converted_string(counter), counter=1, &
                convert_string(converted_string, 'bc_markers'))
        end if

        write(tecplot_io_unit) 299.0
        zone_name = snapshot%zone_name
        write(tecplot_io_unit) (converted_string(counter), counter=1, &
            convert_string(converted_string, zone_name))
        write(tecplot_io_unit) -1
        write(tecplot_io_unit) -1
        write(tecplot_io_unit) snapshot%display_time
        write(tecplot_io_unit) -1
        write(tecplot_io_unit) 0
        write(tecplot_io_unit) 1
        write(tecplot_io_unit) 0
        write(tecplot_io_unit) 0

        do variables_counter = 1, variables_number
            write(tecplot_io_unit) 0
        end do

        if (snapshot%include_ghost_cells) then
            write(tecplot_io_unit) utter_cells_number
        else
            write(tecplot_io_unit) inner_cells_number
        end if

        write(tecplot_io_unit) 0
        write(tecplot_io_unit) 357.0
        write(tecplot_io_unit) 299.0

        do variables_counter = 1, variables_number
            write(tecplot_io_unit) 1
        end do

        write(tecplot_io_unit) 0
        write(tecplot_io_unit) 0
        write(tecplot_io_unit) -1

        do dim = 1, domain_dimensions
            write(tecplot_io_unit) domain_lengths(dim,1)
            write(tecplot_io_unit) domain_lengths(dim,2)
        end do

        do variables_counter = 1, size(snapshot%fields)
            call snapshot%fields(variables_counter)%values%s_ptr%get_field_min( &
                bc_ptr=snapshot%boundaries, min_value=min_value)
            write(tecplot_io_unit) min_value
            call snapshot%fields(variables_counter)%values%s_ptr%get_field_max( &
                bc_ptr=snapshot%boundaries, max_value=max_value)
            write(tecplot_io_unit) max_value
        end do

        if (snapshot%include_ghost_cells) then
            write(tecplot_io_unit) 0.0_dp
            write(tecplot_io_unit) dble(size(snapshot%boundaries%bc_ptr%boundary_types))
        end if
    end subroutine generate_header

    subroutine save_mesh(this, snapshot, tecplot_io_unit)
        class(tecplot_writer), intent(inout) :: this
        type(field_snapshot), intent(in) :: snapshot
        integer, intent(in) :: tecplot_io_unit

        integer :: dimensions, i, j, k, dim
        integer, dimension(3,2) :: loop, inner_loop, utter_loop

        dimensions = snapshot%domain%get_domain_dimensions()
        utter_loop = snapshot%domain%get_local_utter_cells_bounds()
        inner_loop = snapshot%domain%get_local_inner_cells_bounds()

        if (snapshot%include_ghost_cells) then
            loop = utter_loop
        else
            loop = inner_loop
        end if

        do dim = 1, dimensions
            this%io_buffer = snapshot%mesh%mesh_ptr%mesh(dim,:,:,:)
            write(tecplot_io_unit) (((this%io_buffer(i,j,k), &
                i=loop(1,1),loop(1,2)), j=loop(2,1),loop(2,2)), &
                k=loop(3,1),loop(3,2))
        end do
    end subroutine save_mesh

    subroutine save_fields(this, snapshot, tecplot_io_unit)
        class(tecplot_writer), intent(inout) :: this
        type(field_snapshot), intent(in) :: snapshot
        integer, intent(in) :: tecplot_io_unit

        integer, dimension(3,2) :: loop, inner_loop, utter_loop
        integer :: i, j, k, fields_count

        utter_loop = snapshot%domain%get_local_utter_cells_bounds()
        inner_loop = snapshot%domain%get_local_inner_cells_bounds()

        if (snapshot%include_ghost_cells) then
            loop = utter_loop
        else
            loop = inner_loop
        end if

        do fields_count = 1, size(snapshot%fields)
            this%io_buffer = snapshot%fields(fields_count)%values%s_ptr%cells(:,:,:)
            write(tecplot_io_unit) (((this%io_buffer(i,j,k), &
                i=loop(1,1),loop(1,2)), j=loop(2,1),loop(2,2)), &
                k=loop(3,1),loop(3,2))
        end do
    end subroutine save_fields

    subroutine save_bounds(this, snapshot, tecplot_io_unit)
        class(tecplot_writer), intent(inout) :: this
        type(field_snapshot), intent(in) :: snapshot
        integer, intent(in) :: tecplot_io_unit

        integer, dimension(3,2) :: loop
        integer :: i, j, k

        loop = snapshot%domain%get_local_utter_cells_bounds()
        this%io_buffer = snapshot%boundaries%bc_ptr%bc_markers(:,:,:)
        write(tecplot_io_unit) (((this%io_buffer(i,j,k), &
            i=loop(1,1),loop(1,2)), j=loop(2,1),loop(2,2)), &
            k=loop(3,1),loop(3,2))
    end subroutine save_bounds

    integer function convert_string(int_str, char_str) result(num)
        integer, dimension(:), intent(out) :: int_str
        character(len=*), intent(in) :: char_str

        int_str = 0
        do num = 1, len(char_str)
            int_str(num) = ichar(char_str(num:num))
        end do
        int_str(num) = 0
    end function convert_string

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

end module tecplot_writer_class
