module post_processor_operation_class

    use kind_parameters
    use computational_domain_class
    use data_manager_class
    use field_pointers
    use boundary_conditions_class

    implicit none

    private
    public :: post_processor_operation, post_processor_operation_c

    type :: post_processor_operation
        private
        type(field_scalar_cons_pointer), public :: operation_field_scal
        type(field_vector_cons_pointer), public :: operation_field_vect

        type(computational_domain) :: domain

        character(len=40) :: field_name
        character(len=20) :: operation_type
        character(len=16) :: location_mode = 'cell'
        integer :: grad_projection
        integer, dimension(3,2) :: operation_area
        integer, dimension(3) :: operation_area_center
        integer, dimension(3) :: operation_area_distance
    contains
        procedure, private :: read_properties
        procedure, private :: write_properties
        procedure, private :: set_properties

        procedure :: set_point
        procedure :: process_operation
        procedure :: get_values_number
        procedure :: get_output_variable_name
        procedure :: write_log
    end type post_processor_operation

    interface post_processor_operation_c
        module procedure constructor
        module procedure constructor_file
    end interface

contains

    type(post_processor_operation) function constructor( &
        manager, field_name, operation_type, operation_area, operation_area_distance, &
        grad_projection, post_processor_data_file_unit, location_mode)

        type(data_manager), intent(in) :: manager
        character(len=*), intent(in) :: field_name
        character(len=*), intent(in) :: operation_type
        integer, dimension(3,2), intent(in) :: operation_area
        integer, dimension(3), intent(in) :: operation_area_distance
        integer, intent(in) :: grad_projection
        integer, intent(in) :: post_processor_data_file_unit
        character(len=*), intent(in), optional :: location_mode

        call constructor%set_properties( &
            manager, field_name, operation_type, operation_area, &
            operation_area_distance, grad_projection, location_mode)
        call constructor%write_properties(post_processor_data_file_unit)
    end function constructor

    type(post_processor_operation) function constructor_file( &
        manager, post_processor_data_file_unit)

        type(data_manager), intent(in) :: manager
        integer, intent(in) :: post_processor_data_file_unit

        call constructor_file%read_properties(manager, post_processor_data_file_unit)
    end function constructor_file

    subroutine read_properties(this, manager, post_processor_data_file_unit)
        class(post_processor_operation), intent(inout) :: this
        type(data_manager), intent(in) :: manager
        integer, intent(in) :: post_processor_data_file_unit

        character(len=40) :: field_name
        character(len=20) :: operation_type
        character(len=16) :: location_mode
        integer, dimension(3,2) :: operation_area
        integer, dimension(3) :: operation_area_distance
        integer :: grad_projection

        namelist /pproc_operation/ field_name, operation_type, operation_area, &
            operation_area_distance, grad_projection, location_mode

        ! Backward-compatible default for old postprocessor files.
        location_mode = 'cell'
        read(unit=post_processor_data_file_unit, nml=pproc_operation)

        call this%set_properties( &
            manager, field_name, operation_type, operation_area, &
            operation_area_distance, grad_projection, location_mode)
    end subroutine read_properties

    subroutine write_properties(this, post_processor_data_file_unit)
        class(post_processor_operation), intent(in) :: this
        integer, intent(in) :: post_processor_data_file_unit

        character(len=40) :: field_name
        character(len=20) :: operation_type
        character(len=16) :: location_mode
        integer, dimension(3,2) :: operation_area
        integer, dimension(3) :: operation_area_distance
        integer :: grad_projection

        namelist /pproc_operation/ field_name, operation_type, operation_area, &
            operation_area_distance, grad_projection, location_mode

        field_name = trim(this%field_name)
        operation_type = this%operation_type
        location_mode = this%location_mode
        operation_area = this%operation_area
        operation_area_distance = this%operation_area_distance
        grad_projection = this%grad_projection

        write(unit=post_processor_data_file_unit, nml=pproc_operation)
    end subroutine write_properties

    subroutine set_properties( &
        this, manager, field_name, operation_type, operation_area, &
        operation_area_distance, grad_projection, location_mode)

        class(post_processor_operation), intent(inout) :: this
        type(data_manager), intent(in) :: manager
        character(len=*), intent(in) :: field_name
        character(len=*), intent(in) :: operation_type
        integer, dimension(3,2), intent(in) :: operation_area
        integer, dimension(3), intent(in) :: operation_area_distance
        integer, intent(in) :: grad_projection
        character(len=*), intent(in), optional :: location_mode

        integer :: dimensions
        integer, dimension(3,2) :: area
        integer, dimension(3) :: center
        integer, dimension(3) :: distance
        character(len=16) :: requested_location_mode

        type(field_scalar_cons_pointer) :: scal_ptr
        type(field_vector_cons_pointer) :: vect_ptr
        type(field_tensor_cons_pointer) :: tens_ptr

        dimensions = manager%domain%get_domain_dimensions()

        call manager%get_cons_field_pointer_by_name( &
            scal_ptr, vect_ptr, tens_ptr, field_name)

        if (associated(scal_ptr%s_ptr)) then
            this%operation_field_scal%s_ptr => scal_ptr%s_ptr
            this%operation_field_vect%v_ptr => null()
        end if
        if (associated(vect_ptr%v_ptr)) then
            this%operation_field_scal%s_ptr => null()
            this%operation_field_vect%v_ptr => vect_ptr%v_ptr
        end if

        area = operation_area
        distance = operation_area_distance
        center = 0.5*(operation_area(:,2) - abs(operation_area(:,1)))

        if (dimensions < 3) then
            area(dimensions+1:3,:) = 1
            center(dimensions+1:3) = 1
            distance(dimensions+1:3) = 0
        end if

        requested_location_mode = 'cell'
        if (present(location_mode)) requested_location_mode = adjustl(trim(location_mode))

        select case(trim(requested_location_mode))
        case('cell', 'quadratic')
        case default
            error stop 'Postprocessor: location_mode must be cell or quadratic'
        end select

        this%field_name = field_name
        this%operation_type = operation_type
        this%location_mode = requested_location_mode
        this%operation_area = area
        this%operation_area_distance = distance
        this%operation_area_center = center
        this%grad_projection = grad_projection
        this%domain = manager%domain
    end subroutine set_properties

    subroutine write_log(this, log_unit)
        class(post_processor_operation), intent(in) :: this
        integer, intent(in) :: log_unit

        write(log_unit,'(A,A)') 'Operation field name: ', this%field_name
        write(log_unit,'(A,A)') 'Operation type: ', this%operation_type
        write(log_unit,'(A,A)') 'Location mode: ', this%location_mode
        write(log_unit,'(A,I2)') &
            'Gradient projection (if type is gradient): ', this%grad_projection
        write(log_unit,'(A,6I8)') 'Operational area: ', this%operation_area
        write(log_unit,'(A,3I8)') &
            'Operational area center distance from the main point: ', &
            this%operation_area_distance
    end subroutine write_log

    subroutine process_operation(this, bc_ptr, point_indexes, value, point_indexes_real)

#ifdef mpi
        use MPI
#endif

        class(post_processor_operation), intent(in) :: this
        type(boundary_conditions_pointer), intent(in) :: bc_ptr
        integer, dimension(3), intent(out) :: point_indexes
        real(dp), intent(out) :: value
        real(dp), dimension(3), intent(out), optional :: point_indexes_real

        integer, dimension(:), allocatable, save :: point_indexes_array
        real(dp), dimension(:), allocatable, save :: point_indexes_real_array
        real(dp), dimension(:), allocatable, save :: values_array
        integer, dimension(:), allocatable, save :: valid_flags_array
        integer, dimension(:), allocatable, save :: counts_array

        real(dp) :: interm_value(1)
        integer :: interm_point_indexes(3)
        real(dp) :: interm_point_indexes_real(3)
        integer :: interm_valid(1)
        integer :: interm_count(1)

        integer :: selected_point_indexes(3)
        real(dp) :: selected_point_indexes_real(3)
        real(dp) :: selected_value

        integer :: transducer_cell(3), local_cell(3)
        integer :: local_operation_area(3,2)
        integer :: global_inner_cells_bounds(3,2)
        integer :: global_offset(3)

        integer :: dimensions
        integer :: processor_rank, processor_number, mpi_communicator
        integer :: proc, winning_proc, error
        integer :: total_count
        real(dp) :: weighted_sum
        logical :: has_local_area

        processor_rank = this%domain%get_processor_rank()
        mpi_communicator = this%domain%get_mpi_communicator()
        processor_number = this%domain%get_mpi_communicator_size()
        dimensions = this%domain%get_domain_dimensions()
        global_inner_cells_bounds = this%domain%get_global_inner_cells_bounds()
        global_offset = this%domain%get_global_offset()

        if (allocated(point_indexes_array)) then
            if (size(values_array) /= processor_number) then
                deallocate(point_indexes_array, point_indexes_real_array, values_array)
                deallocate(valid_flags_array, counts_array)
            end if
        end if
        if (.not. allocated(point_indexes_array)) then
            allocate(point_indexes_array(3*processor_number))
            allocate(point_indexes_real_array(3*processor_number))
            allocate(values_array(processor_number))
            allocate(valid_flags_array(processor_number))
            allocate(counts_array(processor_number))
        end if

        interm_value(1) = 0.0_dp
        interm_point_indexes = 0
        interm_point_indexes_real = 0.0_dp
        interm_valid(1) = 0
        interm_count(1) = 0

        call make_local_operation_area(this, local_operation_area, has_local_area)

        if (associated(this%operation_field_scal%s_ptr)) then
            select case(trim(this%operation_type))
            case('max')
                if (has_local_area) then
                    if (trim(this%location_mode) == 'quadratic') then
                        call this%operation_field_scal%s_ptr%get_field_max( &
                            bc_ptr, local_operation_area, interm_value(1), &
                            interm_point_indexes, interm_point_indexes_real)
                    else
                        call this%operation_field_scal%s_ptr%get_field_max( &
                            bc_ptr, local_operation_area, interm_value(1), &
                            interm_point_indexes)
                        interm_point_indexes_real = real(interm_point_indexes,dp)
                    end if
                    if (all(interm_point_indexes /= 0)) interm_valid(1) = 1
                    call local_point_to_global( &
                        interm_point_indexes, interm_point_indexes_real, &
                        global_offset, dimensions)
                end if

            case('min')
                if (has_local_area) then
                    if (trim(this%location_mode) == 'quadratic') then
                        call this%operation_field_scal%s_ptr%get_field_min( &
                            bc_ptr, local_operation_area, interm_value(1), &
                            interm_point_indexes, interm_point_indexes_real)
                    else
                        call this%operation_field_scal%s_ptr%get_field_min( &
                            bc_ptr, local_operation_area, interm_value(1), &
                            interm_point_indexes)
                        interm_point_indexes_real = real(interm_point_indexes,dp)
                    end if
                    if (all(interm_point_indexes /= 0)) interm_valid(1) = 1
                    call local_point_to_global( &
                        interm_point_indexes, interm_point_indexes_real, &
                        global_offset, dimensions)
                end if

            case('max_grad')
                if (has_local_area) then
                    if (trim(this%location_mode) == 'quadratic') then
                        call this%operation_field_scal%s_ptr%get_field_max_grad( &
                            bc_ptr, local_operation_area, interm_value(1), &
                            this%grad_projection, interm_point_indexes, &
                            interm_point_indexes_real)
                    else
                        call this%operation_field_scal%s_ptr%get_field_max_grad( &
                            bc_ptr, local_operation_area, interm_value(1), &
                            this%grad_projection, interm_point_indexes)
                        interm_point_indexes_real = real(interm_point_indexes,dp)
                    end if
                    if (all(interm_point_indexes /= 0)) interm_valid(1) = 1
                    call local_point_to_global( &
                        interm_point_indexes, interm_point_indexes_real, &
                        global_offset, dimensions)
                end if

            case('min_grad')
                if (has_local_area) then
                    if (trim(this%location_mode) == 'quadratic') then
                        call this%operation_field_scal%s_ptr%get_field_min_grad( &
                            bc_ptr, local_operation_area, interm_value(1), &
                            this%grad_projection, interm_point_indexes, &
                            interm_point_indexes_real)
                    else
                        call this%operation_field_scal%s_ptr%get_field_min_grad( &
                            bc_ptr, local_operation_area, interm_value(1), &
                            this%grad_projection, interm_point_indexes)
                        interm_point_indexes_real = real(interm_point_indexes,dp)
                    end if
                    if (all(interm_point_indexes /= 0)) interm_valid(1) = 1
                    call local_point_to_global( &
                        interm_point_indexes, interm_point_indexes_real, &
                        global_offset, dimensions)
                end if

            case('transducer')
                transducer_cell = (/ &
                    max(this%operation_area(1,1), global_inner_cells_bounds(1,1)), &
                    max(this%operation_area(2,1), global_inner_cells_bounds(2,1)), &
                    max(this%operation_area(3,1), global_inner_cells_bounds(3,1))/)
                transducer_cell = (/ &
                    min(transducer_cell(1), global_inner_cells_bounds(1,2)), &
                    min(transducer_cell(2), global_inner_cells_bounds(2,2)), &
                    min(transducer_cell(3), global_inner_cells_bounds(3,2))/)

                if (this%domain%cell_is_inside_local_domain(transducer_cell)) then
                    interm_point_indexes = transducer_cell
                    interm_point_indexes_real = real(transducer_cell,dp)
                    local_cell = transducer_cell
                    call this%domain%convert_global_cell_to_local(local_cell)
                    interm_value(1) = this%operation_field_scal%s_ptr%cells( &
                        local_cell(1), local_cell(2), local_cell(3))
                    interm_valid(1) = 1
                end if

            case('sum')
                if (has_local_area) then
                    call this%operation_field_scal%s_ptr%get_field_sum( &
                        bc_ptr, local_operation_area, interm_value(1))
                    interm_valid(1) = 1
                end if

            case('mean')
                if (has_local_area) then
                    call this%operation_field_scal%s_ptr%get_field_mean( &
                        bc_ptr, local_operation_area, interm_value(1), interm_count(1))
                    if (interm_count(1) > 0) interm_valid(1) = 1
                end if
            end select
        end if

        ! Vector operations are not currently implemented in this class. Chemical
        ! vector projections that are registered by long name are resolved by the
        ! data manager as scalar fields and therefore use the scalar path above.

        values_array(processor_rank+1) = interm_value(1)
        point_indexes_array(3*processor_rank+1:3*processor_rank+3) = &
            interm_point_indexes
        point_indexes_real_array(3*processor_rank+1:3*processor_rank+3) = &
            interm_point_indexes_real
        valid_flags_array(processor_rank+1) = interm_valid(1)
        counts_array(processor_rank+1) = interm_count(1)

#ifdef mpi
        call mpi_gather(interm_value, 1, MPI_DOUBLE_PRECISION, &
            values_array, 1, MPI_DOUBLE_PRECISION, 0, mpi_communicator, error)
        call mpi_gather(interm_point_indexes, 3, MPI_INTEGER, &
            point_indexes_array, 3, MPI_INTEGER, 0, mpi_communicator, error)
        call mpi_gather(interm_point_indexes_real, 3, MPI_DOUBLE_PRECISION, &
            point_indexes_real_array, 3, MPI_DOUBLE_PRECISION, &
            0, mpi_communicator, error)
        call mpi_gather(interm_valid, 1, MPI_INTEGER, &
            valid_flags_array, 1, MPI_INTEGER, 0, mpi_communicator, error)
        call mpi_gather(interm_count, 1, MPI_INTEGER, &
            counts_array, 1, MPI_INTEGER, 0, mpi_communicator, error)
#endif

        selected_value = 0.0_dp
        selected_point_indexes = 0
        selected_point_indexes_real = 0.0_dp

        if (processor_rank == 0) then
            select case(trim(this%operation_type))
            case('max', 'max_grad')
                winning_proc = first_valid_processor(valid_flags_array)
                if (winning_proc >= 0) then
                    do proc = winning_proc + 1, processor_number - 1
                        if (valid_flags_array(proc+1) == 0) cycle
                        if (values_array(proc+1) >= values_array(winning_proc+1)) then
                            winning_proc = proc
                        end if
                    end do
                    call copy_processor_result( &
                        winning_proc, values_array, point_indexes_array, &
                        point_indexes_real_array, selected_value, &
                        selected_point_indexes, selected_point_indexes_real)
                end if

            case('min', 'min_grad')
                winning_proc = first_valid_processor(valid_flags_array)
                if (winning_proc >= 0) then
                    do proc = winning_proc + 1, processor_number - 1
                        if (valid_flags_array(proc+1) == 0) cycle
                        if (values_array(proc+1) <= values_array(winning_proc+1)) then
                            winning_proc = proc
                        end if
                    end do
                    call copy_processor_result( &
                        winning_proc, values_array, point_indexes_array, &
                        point_indexes_real_array, selected_value, &
                        selected_point_indexes, selected_point_indexes_real)
                end if

            case('transducer')
                winning_proc = first_valid_processor(valid_flags_array)
                if (winning_proc >= 0) then
                    call copy_processor_result( &
                        winning_proc, values_array, point_indexes_array, &
                        point_indexes_real_array, selected_value, &
                        selected_point_indexes, selected_point_indexes_real)
                end if

            case('sum')
                selected_value = 0.0_dp
                do proc = 0, processor_number - 1
                    if (valid_flags_array(proc+1) /= 0) then
                        selected_value = selected_value + values_array(proc+1)
                    end if
                end do

            case('mean')
                weighted_sum = 0.0_dp
                total_count = 0
                do proc = 0, processor_number - 1
                    if (valid_flags_array(proc+1) == 0) cycle
                    weighted_sum = weighted_sum + &
                        values_array(proc+1)*real(counts_array(proc+1),dp)
                    total_count = total_count + counts_array(proc+1)
                end do
                if (total_count > 0) then
                    selected_value = weighted_sum/real(total_count,dp)
                end if
            end select
        end if

#ifdef mpi
        call mpi_bcast(selected_value, 1, MPI_DOUBLE_PRECISION, &
            0, mpi_communicator, error)
        call mpi_bcast(selected_point_indexes, 3, MPI_INTEGER, &
            0, mpi_communicator, error)
        call mpi_bcast(selected_point_indexes_real, 3, MPI_DOUBLE_PRECISION, &
            0, mpi_communicator, error)
#endif

        value = selected_value
        point_indexes = selected_point_indexes
        if (present(point_indexes_real)) then
            point_indexes_real = selected_point_indexes_real
        end if
    end subroutine process_operation

    integer function get_values_number(this)
        class(post_processor_operation), intent(in) :: this
        get_values_number = 1
    end function get_values_number

    function get_output_variable_name(this) result(variable_name)
        class(post_processor_operation), intent(in) :: this
        character(len=128) :: variable_name

        character(len=5), dimension(3) :: axis_names
        character(len=32) :: offset_text
        integer :: dimensions
        integer :: dim

        dimensions = this%domain%get_domain_dimensions()
        axis_names = this%domain%get_axis_names()

        variable_name = trim(this%field_name) // '_' // trim(this%operation_type)

        if ((trim(this%operation_type) == 'max_grad' .or. &
             trim(this%operation_type) == 'min_grad') .and. &
            this%grad_projection >= 1 .and. this%grad_projection <= dimensions) then
            variable_name = trim(variable_name) // '_' // &
                trim(axis_names(this%grad_projection))
        end if

        do dim = 1, dimensions
            if (this%operation_area_distance(dim) == 0) cycle
            write(offset_text,'(SP,I0)') this%operation_area_distance(dim)
            variable_name = trim(variable_name) // '_d' // trim(axis_names(dim)) // &
                trim(adjustl(offset_text))
        end do
    end function get_output_variable_name

    subroutine set_point(this, point_indexes)
        class(post_processor_operation), intent(inout) :: this
        integer, dimension(3), intent(in) :: point_indexes
        integer, dimension(3) :: new_center
        integer, dimension(3) :: shift

        new_center = point_indexes + this%operation_area_distance
        shift = new_center - this%operation_area_center
        this%operation_area_center = new_center
        this%operation_area(:,1) = this%operation_area(:,1) + shift
        this%operation_area(:,2) = this%operation_area(:,2) + shift
    end subroutine set_point

    subroutine make_local_operation_area(this, local_area, has_local_area)
        class(post_processor_operation), intent(in) :: this
        integer, dimension(3,2), intent(out) :: local_area
        logical, intent(out) :: has_local_area

        integer :: dimensions, dim
        integer :: global_offset(3)
        integer :: local_inner_bounds(3,2)
        integer :: global_local_low, global_local_high
        integer :: overlap_low, overlap_high

        dimensions = this%domain%get_domain_dimensions()
        global_offset = this%domain%get_global_offset()
        local_inner_bounds = this%domain%get_local_inner_cells_bounds()

        local_area = 1
        has_local_area = .true.

        do dim = 1, dimensions
            global_local_low = global_offset(dim) + local_inner_bounds(dim,1)
            global_local_high = global_offset(dim) + local_inner_bounds(dim,2)

            overlap_low = max(this%operation_area(dim,1), global_local_low)
            overlap_high = min(this%operation_area(dim,2), global_local_high)

            if (overlap_low > overlap_high) then
                has_local_area = .false.
                return
            end if

            local_area(dim,1) = overlap_low - global_offset(dim)
            local_area(dim,2) = overlap_high - global_offset(dim)
        end do
    end subroutine make_local_operation_area

    subroutine local_point_to_global( &
        point_indexes, point_indexes_real, global_offset, dimensions)

        integer, dimension(3), intent(inout) :: point_indexes
        real(dp), dimension(3), intent(inout) :: point_indexes_real
        integer, dimension(3), intent(in) :: global_offset
        integer, intent(in) :: dimensions
        integer :: dim

        if (any(point_indexes(1:dimensions) == 0)) return

        do dim = 1, dimensions
            point_indexes(dim) = point_indexes(dim) + global_offset(dim)
            point_indexes_real(dim) = &
                point_indexes_real(dim) + real(global_offset(dim),dp)
        end do
    end subroutine local_point_to_global

    integer function first_valid_processor(valid_flags)
        integer, dimension(:), intent(in) :: valid_flags
        integer :: proc

        first_valid_processor = -1
        do proc = 0, size(valid_flags) - 1
            if (valid_flags(proc+1) /= 0) then
                first_valid_processor = proc
                return
            end if
        end do
    end function first_valid_processor

    subroutine copy_processor_result( &
        proc, values, indexes, indexes_real, selected_value, &
        selected_indexes, selected_indexes_real)

        integer, intent(in) :: proc
        real(dp), dimension(:), intent(in) :: values
        integer, dimension(:), intent(in) :: indexes
        real(dp), dimension(:), intent(in) :: indexes_real
        real(dp), intent(out) :: selected_value
        integer, dimension(3), intent(out) :: selected_indexes
        real(dp), dimension(3), intent(out) :: selected_indexes_real

        selected_value = values(proc+1)
        selected_indexes = indexes(3*proc+1:3*proc+3)
        selected_indexes_real = indexes_real(3*proc+1:3*proc+3)
    end subroutine copy_processor_result

end module post_processor_operation_class
