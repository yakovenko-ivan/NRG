module cgns_writer_class

    use, intrinsic :: iso_c_binding, only: c_char, c_double, c_int, c_int64_t, c_null_char

    use kind_parameters
    use global_data
    use field_writer_class
    use field_output_types

    implicit none

    private
    public :: cgns_writer, cgns_backend_enabled

#ifdef NRG_CGNS
    logical, parameter :: cgns_backend_enabled = .true.
#else
    logical, parameter :: cgns_backend_enabled = .false.
#endif

    type, extends(field_writer) :: cgns_writer
    contains
        procedure :: write_snapshot => cgns_write_snapshot
    end type cgns_writer

#ifdef NRG_CGNS
    interface
        integer(c_int) function nrg_cgns_begin(filename, dataset_name, zone_name, &
                cell_dim, phys_dim, cell_counts, time_seconds, output_index, &
                coordinate_system, fn, base, zone, solution) bind(C, name='nrg_cgns_begin')
            import :: c_char, c_double, c_int, c_int64_t
            character(kind=c_char), intent(in) :: filename(*)
            character(kind=c_char), intent(in) :: dataset_name(*)
            character(kind=c_char), intent(in) :: zone_name(*)
            integer(c_int), value :: cell_dim, phys_dim
            integer(c_int64_t), intent(in) :: cell_counts(3)
            real(c_double), value :: time_seconds
            integer(c_int), value :: output_index
            character(kind=c_char), intent(in) :: coordinate_system(*)
            integer(c_int), intent(out) :: fn, base, zone, solution
        end function nrg_cgns_begin

        integer(c_int) function nrg_cgns_write_coordinate(fn, base, zone, name, data) &
                bind(C, name='nrg_cgns_write_coordinate')
            import :: c_char, c_double, c_int
            integer(c_int), value :: fn, base, zone
            character(kind=c_char), intent(in) :: name(*)
            real(c_double), intent(in) :: data(*)
        end function nrg_cgns_write_coordinate

        integer(c_int) function nrg_cgns_set_solution_rind(fn, base, zone, solution, &
                index_dim, rind_data) bind(C, name='nrg_cgns_set_solution_rind')
            import :: c_int
            integer(c_int), value :: fn, base, zone, solution, index_dim
            integer(c_int), intent(in) :: rind_data(6)
        end function nrg_cgns_set_solution_rind

        integer(c_int) function nrg_cgns_write_field_double(fn, base, zone, solution, &
                name, data, nrg_long_name, nrg_short_name, units) &
                bind(C, name='nrg_cgns_write_field_double')
            import :: c_char, c_double, c_int
            integer(c_int), value :: fn, base, zone, solution
            character(kind=c_char), intent(in) :: name(*)
            real(c_double), intent(in) :: data(*)
            character(kind=c_char), intent(in) :: nrg_long_name(*)
            character(kind=c_char), intent(in) :: nrg_short_name(*)
            character(kind=c_char), intent(in) :: units(*)
        end function nrg_cgns_write_field_double

        integer(c_int) function nrg_cgns_write_field_int(fn, base, zone, solution, &
                name, data) bind(C, name='nrg_cgns_write_field_int')
            import :: c_char, c_int
            integer(c_int), value :: fn, base, zone, solution
            character(kind=c_char), intent(in) :: name(*)
            integer(c_int), intent(in) :: data(*)
        end function nrg_cgns_write_field_int

        integer(c_int) function nrg_cgns_close(fn) bind(C, name='nrg_cgns_close')
            import :: c_int
            integer(c_int), value :: fn
        end function nrg_cgns_close
    end interface
#endif

contains

    subroutine cgns_write_snapshot(this, snapshot, output_folder, file_stem)
        class(cgns_writer), intent(inout) :: this
        type(field_snapshot), intent(in) :: snapshot
        character(len=*), intent(in) :: output_folder
        character(len=*), intent(in) :: file_stem

#ifndef NRG_CGNS
        error stop 'CGNS output requested, but NRG was built with NRG_ENABLE_CGNS=OFF'
#else
        integer :: dimensions, dim, field_index
        integer :: processor_count
        integer, dimension(3) :: cell_counts, vertex_counts
        integer, dimension(3,2) :: loop, storage_loop
        integer(c_int64_t) :: cell_counts_c(3)
        integer(c_int) :: fn, base, zone, solution, status
        integer(c_int) :: rind_data(6)
        real(c_double), allocatable :: coordinate_data(:)
        real(c_double), allocatable :: field_data(:)
        integer(c_int), allocatable :: marker_data(:)
        character(len=256) :: file_name
        character(len=32) :: name
        character(len=32) :: units
        character(len=20) :: coordinate_system
        character(len=5), dimension(3) :: axis_names
        character(kind=c_char), allocatable :: c_filename(:), c_name(:), c_coord_system(:)
        character(kind=c_char), allocatable :: c_dataset_name(:), c_zone_name(:)
        character(kind=c_char), allocatable :: c_long_name(:), c_short_name(:), c_units(:)

        processor_count = snapshot%domain%get_mpi_communicator_size()
        if (processor_count /= 1) then
            error stop 'CGNS writer currently supports one MPI rank; use Tecplot for parallel output until PCGNS is enabled'
        end if

        dimensions = snapshot%domain%get_domain_dimensions()
        coordinate_system = trim(snapshot%domain%get_coordinate_system_name())
        axis_names = snapshot%domain%get_axis_names()

        cell_counts = snapshot%domain%get_global_cells_number()
        cell_counts(1:dimensions) = cell_counts(1:dimensions) - 2
        cell_counts(dimensions+1:3) = 1
        vertex_counts = 1
        vertex_counts(1:dimensions) = cell_counts(1:dimensions) + 1
        cell_counts_c = int(cell_counts, c_int64_t)

        call make_directory(output_folder)
        file_name = trim(output_folder) // trim(fold_sep) // trim(file_stem) // '.cgns'

        c_filename = to_c_string(trim(file_name))
        c_dataset_name = to_c_string(trim(snapshot%dataset_name))
        c_zone_name = to_c_string(trim(snapshot%zone_name))
        c_coord_system = to_c_string(trim(coordinate_system))

        status = nrg_cgns_begin(c_filename, c_dataset_name, c_zone_name, &
            int(dimensions,c_int), int(dimensions,c_int), cell_counts_c, &
            real(snapshot%time_seconds,c_double), int(snapshot%output_index,c_int), &
            c_coord_system, fn, base, zone, solution)
        call check_status(status, 'creating CGNS snapshot')

        do dim = 1, dimensions
            call build_vertex_coordinate(snapshot, dim, vertex_counts, coordinate_data)
            name = trim(axis_names(dim))
            c_name = to_c_string(trim(name))
            status = nrg_cgns_write_coordinate(fn, base, zone, c_name, coordinate_data)
            call check_status(status, 'writing coordinate '//trim(name))
            deallocate(coordinate_data)
        end do

        if (snapshot%include_ghost_cells) then
            rind_data = 0_c_int
            do dim = 1, dimensions
                rind_data(2*dim-1) = 1_c_int
                rind_data(2*dim) = 1_c_int
            end do
            status = nrg_cgns_set_solution_rind(fn, base, zone, solution, &
                int(dimensions,c_int), rind_data)
            call check_status(status, 'writing FlowSolution rind metadata')
        end if

        storage_loop = snapshot%domain%get_local_utter_cells_bounds()
        if (snapshot%include_ghost_cells) then
            loop = storage_loop
        else
            loop = snapshot%domain%get_local_inner_cells_bounds()
        end if

        do field_index = 1, size(snapshot%fields)
            call pack_real_field(snapshot%fields(field_index)%values%s_ptr%cells, storage_loop(:,1), loop, field_data)
            name = cgns_field_name(snapshot, field_index)
            units = cgns_field_units(snapshot, field_index)
            c_name = to_c_string(trim(name))
            c_long_name = to_c_string(trim(snapshot%fields(field_index)%nrg_name))
            c_short_name = to_c_string(trim(snapshot%fields(field_index)%short_name))
            c_units = to_c_string(trim(units))
            status = nrg_cgns_write_field_double(fn, base, zone, solution, c_name, field_data, &
                c_long_name, c_short_name, c_units)
            call check_status(status, 'writing field '//trim(name))
            deallocate(field_data)
        end do

        if (snapshot%include_ghost_cells) then
            call pack_integer_field(snapshot%boundaries%bc_ptr%bc_markers, storage_loop(:,1), loop, marker_data)
            c_name = to_c_string('bc_markers')
            status = nrg_cgns_write_field_int(fn, base, zone, solution, c_name, marker_data)
            call check_status(status, 'writing NRG boundary markers')
            deallocate(marker_data)
        end if

        status = nrg_cgns_close(fn)
        call check_status(status, 'closing CGNS snapshot')
#endif
    end subroutine cgns_write_snapshot

#ifdef NRG_CGNS
    subroutine build_vertex_coordinate(snapshot, coordinate_index, vertex_counts, data)
        type(field_snapshot), intent(in) :: snapshot
        integer, intent(in) :: coordinate_index
        integer, dimension(3), intent(in) :: vertex_counts
        real(c_double), allocatable, intent(out) :: data(:)

        real(dp), dimension(:,:), allocatable :: lengths
        real(dp) :: dx
        integer :: i, j, k, n

        lengths = snapshot%domain%get_domain_lengths()
        dx = (lengths(coordinate_index,2) - lengths(coordinate_index,1)) / &
             real(vertex_counts(coordinate_index)-1, dp)

        allocate(data(product(vertex_counts(1:snapshot%domain%get_domain_dimensions()))))
        n = 0
        do k = 1, vertex_counts(3)
        do j = 1, vertex_counts(2)
        do i = 1, vertex_counts(1)
            n = n + 1
            select case(coordinate_index)
                case(1)
                    data(n) = real(lengths(1,1) + real(i-1,dp)*dx, c_double)
                case(2)
                    data(n) = real(lengths(2,1) + real(j-1,dp)*dx, c_double)
                case(3)
                    data(n) = real(lengths(3,1) + real(k-1,dp)*dx, c_double)
            end select
        end do
        end do
        end do
    end subroutine build_vertex_coordinate

    subroutine pack_real_field(source, source_lower, loop, data)
        real(dp), dimension(:,:,:), intent(in) :: source
        integer, dimension(3), intent(in) :: source_lower
        integer, dimension(3,2), intent(in) :: loop
        real(c_double), allocatable, intent(out) :: data(:)

        integer :: i, j, k, n, total

        total = product(loop(:,2) - loop(:,1) + 1)
        allocate(data(total))
        n = 0
        do k = loop(3,1), loop(3,2)
        do j = loop(2,1), loop(2,2)
        do i = loop(1,1), loop(1,2)
            n = n + 1
            data(n) = real(source(i-source_lower(1)+1, &
                                  j-source_lower(2)+1, &
                                  k-source_lower(3)+1), c_double)
        end do
        end do
        end do
    end subroutine pack_real_field

    subroutine pack_integer_field(source, source_lower, loop, data)
        integer(i1), dimension(:,:,:), intent(in) :: source
        integer, dimension(3), intent(in) :: source_lower
        integer, dimension(3,2), intent(in) :: loop
        integer(c_int), allocatable, intent(out) :: data(:)

        integer :: i, j, k, n, total

        total = product(loop(:,2) - loop(:,1) + 1)
        allocate(data(total))
        n = 0
        do k = loop(3,1), loop(3,2)
        do j = loop(2,1), loop(2,2)
        do i = loop(1,1), loop(1,2)
            n = n + 1
            data(n) = int(source(i-source_lower(1)+1, &
                                j-source_lower(2)+1, &
                                k-source_lower(3)+1), c_int)
        end do
        end do
        end do
    end subroutine pack_integer_field

    function cgns_field_name(snapshot, field_index) result(name)
        type(field_snapshot), intent(in) :: snapshot
        integer, intent(in) :: field_index
        character(len=32) :: name
        character(len=40) :: short_name

        ! Tecplot uses the NRG field short name directly.  Keep the CGNS
        ! FlowSolution_t DataArray_t names identical so post-processing and
        ! agent tooling see one format-independent field naming contract.
        short_name = trim(snapshot%fields(field_index)%short_name)

        if (len_trim(short_name) > len(name)) then
            write(*,'(A,A,A,I0,A)') &
                'CGNS output error: NRG short field name "', trim(short_name), &
                '" exceeds the CGNS 32-character node-name limit (', &
                len_trim(short_name), ' characters). '
            error stop 'CGNS field name is too long to preserve Tecplot naming exactly'
        end if

        name = trim(short_name)
    end function cgns_field_name


    function cgns_field_units(snapshot, field_index) result(units)
        type(field_snapshot), intent(in) :: snapshot
        integer, intent(in) :: field_index
        character(len=32) :: units
        character(len=80) :: long_name

        long_name = trim(snapshot%fields(field_index)%nrg_name)
        units = ''

        select case(trim(long_name))
            case('pressure', 'pressure_static', 'pressure_dynamic')
                units = 'Pa'
            case('temperature')
                units = 'K'
            case('density')
                units = 'kg/m^3'
            case default
                if (index(trim(long_name), 'velocity(') == 1) units = 'm/s'
                if (index(trim(long_name), 'specie_mass_fraction(') == 1) units = '1'
        end select
    end function cgns_field_units

    function to_c_string(text) result(c_text)
        character(len=*), intent(in) :: text
        character(kind=c_char), allocatable :: c_text(:)
        integer :: i, n

        n = len_trim(text)
        allocate(c_text(n+1))
        do i = 1, n
            c_text(i) = achar(iachar(text(i:i)), kind=c_char)
        end do
        c_text(n+1) = c_null_char
    end function to_c_string

    subroutine check_status(status, operation)
        integer(c_int), intent(in) :: status
        character(len=*), intent(in) :: operation

        if (status /= 0_c_int) then
            write(*,'(A)') 'CGNS writer failed while '//trim(operation)//'.'
            error stop 'CGNS writer failure; see CGNS diagnostic above'
        end if
    end subroutine check_status
#endif

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

end module cgns_writer_class
