module field_scalar_class

    use global_data
    use kind_parameters
    use computational_domain_class
    use boundary_conditions_class

    implicit none

    private
    public :: field_scalar, field_scalar_flow, field_scalar_cons
    public :: field_scalar_flow_c, field_scalar_cons_c

    type :: field_scalar
        character(len=80) :: name_long
        character(len=40) :: name_short
    end type field_scalar

    type, extends(field_scalar) :: field_scalar_cons
        real(dp), dimension(:,:,:), allocatable :: cells
        integer :: dimensions
    contains
        procedure :: get_field_max
        procedure :: get_field_min
        procedure :: get_field_mean
        procedure :: get_field_max_grad
        procedure :: get_field_min_grad
        procedure :: get_field_sum
    end type field_scalar_cons

    type, extends(field_scalar) :: field_scalar_flow
        real(dp), dimension(:,:,:,:), allocatable :: cells
    end type field_scalar_flow

    interface field_scalar_flow_c
        module procedure constructor_flow
    end interface

    interface field_scalar_cons_c
        module procedure constructor_cons
    end interface

    interface field_scalar_c
        module procedure constructor
    end interface

contains

    type(field_scalar) function constructor(field_name_long, field_name_short)
        character(len=*), intent(in) :: field_name_long
        character(len=*), intent(in) :: field_name_short

        constructor%name_short = field_name_short
        constructor%name_long = field_name_long
    end function constructor

    type(field_scalar_cons) function constructor_cons(field_name_long, field_name_short, domain)
        character(len=*), intent(in) :: field_name_long
        character(len=*), intent(in) :: field_name_short
        type(computational_domain), intent(in) :: domain

        integer :: dimensions
        integer, dimension(3,2) :: allocation_bounds

        constructor_cons%field_scalar = field_scalar_c(field_name_long, field_name_short)

        dimensions = domain%get_domain_dimensions()
        allocation_bounds = domain%get_local_utter_cells_bounds()

        allocate(constructor_cons%cells( &
            allocation_bounds(1,1):allocation_bounds(1,2), &
            allocation_bounds(2,1):allocation_bounds(2,2), &
            allocation_bounds(3,1):allocation_bounds(3,2)))

        constructor_cons%cells = 0.0_dp
        constructor_cons%dimensions = dimensions
    end function constructor_cons

    type(field_scalar_flow) function constructor_flow(field_name_long, field_name_short, domain)
        character(len=*), intent(in) :: field_name_long
        character(len=*), intent(in) :: field_name_short
        type(computational_domain), intent(in) :: domain

        integer :: dimensions
        integer, dimension(3,2) :: allocation_bounds

        constructor_flow%field_scalar = field_scalar_c(field_name_long, field_name_short)

        dimensions = domain%get_domain_dimensions()
        allocation_bounds = domain%get_local_utter_faces_bounds()

        allocate(constructor_flow%cells( &
            dimensions, &
            allocation_bounds(1,1):allocation_bounds(1,2), &
            allocation_bounds(2,1):allocation_bounds(2,2), &
            allocation_bounds(3,1):allocation_bounds(3,2)))

        constructor_flow%cells = 0.0_dp
    end function constructor_flow

    subroutine get_field_max(this, bc_ptr, bound, max_value, max_indexes, max_indexes_real)
        class(field_scalar_cons), intent(in) :: this
        type(boundary_conditions_pointer), intent(in) :: bc_ptr
        integer, dimension(3,2), intent(in), optional :: bound
        real(dp), intent(out) :: max_value
        integer, dimension(3), intent(out), optional :: max_indexes
        real(dp), dimension(3), intent(out), optional :: max_indexes_real

        integer, dimension(3,2) :: bound_sl
        integer, dimension(3) :: indexes
        integer :: i, j, k, dim
        logical :: found

        call scalar_bounds(this, bound_sl, bound)

        found = .false.
        max_value = -huge(1.0_dp)
        indexes = 0

        do k = bound_sl(3,1), bound_sl(3,2)
        do j = bound_sl(2,1), bound_sl(2,2)
        do i = bound_sl(1,1), bound_sl(1,2)
            if (bc_ptr%bc_ptr%bc_markers(i,j,k) /= 0) cycle
            if ((.not. found) .or. this%cells(i,j,k) >= max_value) then
                max_value = this%cells(i,j,k)
                indexes = (/i,j,k/)
                found = .true.
            end if
        end do
        end do
        end do

        if (.not. found) then
            max_value = 0.0_dp
            if (present(max_indexes)) max_indexes = 0
            if (present(max_indexes_real)) max_indexes_real = 0.0_dp
            return
        end if

        if (present(max_indexes)) max_indexes = indexes
        if (present(max_indexes_real)) then
            max_indexes_real = real(indexes, dp)
            call refine_scalar_extremum( &
                this, bc_ptr, bound_sl, indexes, 'max', max_indexes_real)
        end if
    end subroutine get_field_max

    subroutine get_field_min(this, bc_ptr, bound, min_value, min_indexes, min_indexes_real)
        class(field_scalar_cons), intent(in) :: this
        type(boundary_conditions_pointer), intent(in) :: bc_ptr
        integer, dimension(3,2), intent(in), optional :: bound
        real(dp), intent(out) :: min_value
        integer, dimension(3), intent(out), optional :: min_indexes
        real(dp), dimension(3), intent(out), optional :: min_indexes_real

        integer, dimension(3,2) :: bound_sl
        integer, dimension(3) :: indexes
        integer :: i, j, k
        logical :: found

        call scalar_bounds(this, bound_sl, bound)

        found = .false.
        min_value = huge(1.0_dp)
        indexes = 0

        do k = bound_sl(3,1), bound_sl(3,2)
        do j = bound_sl(2,1), bound_sl(2,2)
        do i = bound_sl(1,1), bound_sl(1,2)
            if (bc_ptr%bc_ptr%bc_markers(i,j,k) /= 0) cycle
            if ((.not. found) .or. this%cells(i,j,k) <= min_value) then
                min_value = this%cells(i,j,k)
                indexes = (/i,j,k/)
                found = .true.
            end if
        end do
        end do
        end do

        if (.not. found) then
            min_value = 0.0_dp
            if (present(min_indexes)) min_indexes = 0
            if (present(min_indexes_real)) min_indexes_real = 0.0_dp
            return
        end if

        if (present(min_indexes)) min_indexes = indexes
        if (present(min_indexes_real)) then
            min_indexes_real = real(indexes, dp)
            call refine_scalar_extremum( &
                this, bc_ptr, bound_sl, indexes, 'min', min_indexes_real)
        end if
    end subroutine get_field_min

    subroutine get_field_mean(this, bc_ptr, bound, field_mean_value, cells_count)
        class(field_scalar_cons), intent(in) :: this
        type(boundary_conditions_pointer), intent(in) :: bc_ptr
        integer, dimension(3,2), intent(in) :: bound
        real(dp), intent(out) :: field_mean_value
        integer, intent(out), optional :: cells_count

        real(dp) :: field_sum
        integer, dimension(3,2) :: bound_sl
        integer :: i, j, k, dim, count

        do dim = 1, 3
            bound_sl(dim,1) = max(bound(dim,1), lbound(this%cells,dim))
            bound_sl(dim,2) = min(bound(dim,2), ubound(this%cells,dim))
        end do

        field_sum = 0.0_dp
        count = 0

        do k = bound_sl(3,1), bound_sl(3,2)
        do j = bound_sl(2,1), bound_sl(2,2)
        do i = bound_sl(1,1), bound_sl(1,2)
            if (bc_ptr%bc_ptr%bc_markers(i,j,k) /= 0) cycle
            field_sum = field_sum + this%cells(i,j,k)
            count = count + 1
        end do
        end do
        end do

        if (count > 0) then
            field_mean_value = field_sum / real(count, dp)
        else
            field_mean_value = 0.0_dp
        end if

        if (present(cells_count)) cells_count = count
    end subroutine get_field_mean

    subroutine get_field_max_grad( &
        this, bc_ptr, bound, max_grad_value, grad_proj, max_grad_indexes, max_grad_indexes_real)

        class(field_scalar_cons), intent(in) :: this
        type(boundary_conditions_pointer), intent(in) :: bc_ptr
        integer, dimension(3,2), intent(in) :: bound
        real(dp), intent(out) :: max_grad_value
        integer, intent(in), optional :: grad_proj
        integer, dimension(3), intent(out) :: max_grad_indexes
        real(dp), dimension(3), intent(out), optional :: max_grad_indexes_real

        integer, dimension(3,2) :: bound_sl
        integer, dimension(3) :: indexes
        real(dp) :: candidate
        integer :: i, j, k
        logical :: found

        call gradient_bounds(this, bound, bound_sl)

        found = .false.
        max_grad_value = -huge(1.0_dp)
        indexes = 0

        do k = bound_sl(3,1), bound_sl(3,2)
        do j = bound_sl(2,1), bound_sl(2,2)
        do i = bound_sl(1,1), bound_sl(1,2)
            if (bc_ptr%bc_ptr%bc_markers(i,j,k) /= 0) cycle

            candidate = gradient_metric(this, bc_ptr, (/i,j,k/), grad_proj)
            if ((.not. found) .or. candidate >= max_grad_value) then
                max_grad_value = candidate
                indexes = (/i,j,k/)
                found = .true.
            end if
        end do
        end do
        end do

        if (.not. found) then
            max_grad_value = 0.0_dp
            max_grad_indexes = 0
            if (present(max_grad_indexes_real)) max_grad_indexes_real = 0.0_dp
            return
        end if

        max_grad_indexes = indexes
        if (present(max_grad_indexes_real)) then
            max_grad_indexes_real = real(indexes, dp)
            call refine_gradient_extremum( &
                this, bc_ptr, bound_sl, indexes, grad_proj, 'max', max_grad_indexes_real)
        end if
    end subroutine get_field_max_grad

    subroutine get_field_min_grad( &
        this, bc_ptr, bound, min_grad_value, grad_proj, min_grad_indexes, min_grad_indexes_real)

        class(field_scalar_cons), intent(in) :: this
        type(boundary_conditions_pointer), intent(in) :: bc_ptr
        integer, dimension(3,2), intent(in) :: bound
        real(dp), intent(out) :: min_grad_value
        integer, intent(in), optional :: grad_proj
        integer, dimension(3), intent(out) :: min_grad_indexes
        real(dp), dimension(3), intent(out), optional :: min_grad_indexes_real

        integer, dimension(3,2) :: bound_sl
        integer, dimension(3) :: indexes
        real(dp) :: candidate
        integer :: i, j, k
        logical :: found

        call gradient_bounds(this, bound, bound_sl)

        found = .false.
        min_grad_value = huge(1.0_dp)
        indexes = 0

        do k = bound_sl(3,1), bound_sl(3,2)
        do j = bound_sl(2,1), bound_sl(2,2)
        do i = bound_sl(1,1), bound_sl(1,2)
            if (bc_ptr%bc_ptr%bc_markers(i,j,k) /= 0) cycle

            candidate = gradient_metric(this, bc_ptr, (/i,j,k/), grad_proj)
            if ((.not. found) .or. candidate <= min_grad_value) then
                min_grad_value = candidate
                indexes = (/i,j,k/)
                found = .true.
            end if
        end do
        end do
        end do

        if (.not. found) then
            min_grad_value = 0.0_dp
            min_grad_indexes = 0
            if (present(min_grad_indexes_real)) min_grad_indexes_real = 0.0_dp
            return
        end if

        min_grad_indexes = indexes
        if (present(min_grad_indexes_real)) then
            min_grad_indexes_real = real(indexes, dp)
            call refine_gradient_extremum( &
                this, bc_ptr, bound_sl, indexes, grad_proj, 'min', min_grad_indexes_real)
        end if
    end subroutine get_field_min_grad

    subroutine get_field_sum(this, bc_ptr, bound, sum_value)
        class(field_scalar_cons), intent(in) :: this
        type(boundary_conditions_pointer), intent(in) :: bc_ptr
        integer, dimension(3,2), intent(in) :: bound
        real(dp), intent(out) :: sum_value

        integer, dimension(3,2) :: bound_sl
        integer :: i, j, k, dim

        sum_value = 0.0_dp
        bound_sl = 1

        do dim = 1, this%dimensions
            bound_sl(dim,1) = max(bound(dim,1), lbound(this%cells,dim)+1)
            bound_sl(dim,2) = min(bound(dim,2), ubound(this%cells,dim)-1)
        end do

        do k = bound_sl(3,1), bound_sl(3,2)
        do j = bound_sl(2,1), bound_sl(2,2)
        do i = bound_sl(1,1), bound_sl(1,2)
            if (bc_ptr%bc_ptr%bc_markers(i,j,k) == 0) then
                sum_value = sum_value + this%cells(i,j,k)
            end if
        end do
        end do
        end do
    end subroutine get_field_sum

    subroutine scalar_bounds(this, bound_sl, bound)
        class(field_scalar_cons), intent(in) :: this
        integer, dimension(3,2), intent(out) :: bound_sl
        integer, dimension(3,2), intent(in), optional :: bound
        integer :: dim

        if (present(bound)) then
            do dim = 1, 3
                bound_sl(dim,1) = max(bound(dim,1), lbound(this%cells,dim))
                bound_sl(dim,2) = min(bound(dim,2), ubound(this%cells,dim))
            end do
        else
            do dim = 1, 3
                bound_sl(dim,1) = lbound(this%cells,dim)
                bound_sl(dim,2) = ubound(this%cells,dim)
            end do
        end if
    end subroutine scalar_bounds

    subroutine gradient_bounds(this, bound, bound_sl)
        class(field_scalar_cons), intent(in) :: this
        integer, dimension(3,2), intent(in) :: bound
        integer, dimension(3,2), intent(out) :: bound_sl
        integer :: dim

        bound_sl = 1
        do dim = 1, this%dimensions
            bound_sl(dim,1) = max(bound(dim,1), lbound(this%cells,dim)+1)
            bound_sl(dim,2) = min(bound(dim,2), ubound(this%cells,dim)-1)
        end do
    end subroutine gradient_bounds

    pure subroutine quadratic_extremum_offset( &
        value_minus, value_center, value_plus, extremum_type, offset, valid)

        real(dp), intent(in) :: value_minus
        real(dp), intent(in) :: value_center
        real(dp), intent(in) :: value_plus
        character(len=*), intent(in) :: extremum_type
        real(dp), intent(out) :: offset
        logical, intent(out) :: valid

        real(dp) :: denominator, scale

        offset = 0.0_dp
        valid = .false.

        denominator = value_minus - 2.0_dp*value_center + value_plus
        scale = max(abs(value_minus), abs(value_center), abs(value_plus), tiny(1.0_dp))

        if (abs(denominator) <= 100.0_dp*epsilon(1.0_dp)*scale) return

        select case(trim(extremum_type))
        case('min')
            if (denominator <= 0.0_dp) return
            if (value_center > value_minus .or. value_center > value_plus) return
        case('max')
            if (denominator >= 0.0_dp) return
            if (value_center < value_minus .or. value_center < value_plus) return
        case default
            return
        end select

        offset = 0.5_dp*(value_minus - value_plus)/denominator
        if (abs(offset) > 0.5_dp) then
            offset = 0.0_dp
            return
        end if

        valid = .true.
    end subroutine quadratic_extremum_offset

    subroutine refine_scalar_extremum( &
        this, bc_ptr, bound, indexes, extremum_type, indexes_real)

        class(field_scalar_cons), intent(in) :: this
        type(boundary_conditions_pointer), intent(in) :: bc_ptr
        integer, dimension(3,2), intent(in) :: bound
        integer, dimension(3), intent(in) :: indexes
        character(len=*), intent(in) :: extremum_type
        real(dp), dimension(3), intent(inout) :: indexes_real

        integer, dimension(3) :: index_minus, index_plus
        real(dp) :: value_minus, value_center, value_plus, offset
        integer :: dim
        logical :: valid

        value_center = this%cells(indexes(1), indexes(2), indexes(3))

        do dim = 1, this%dimensions
            index_minus = indexes - I_m(dim,:)
            index_plus = indexes + I_m(dim,:)

            if (.not. point_inside_bound(index_minus, bound, this%dimensions)) cycle
            if (.not. point_inside_bound(index_plus, bound, this%dimensions)) cycle
            if (.not. point_inside_array(this, index_minus)) cycle
            if (.not. point_inside_array(this, index_plus)) cycle
            if (bc_ptr%bc_ptr%bc_markers( &
                index_minus(1), index_minus(2), index_minus(3)) /= 0) cycle
            if (bc_ptr%bc_ptr%bc_markers( &
                index_plus(1), index_plus(2), index_plus(3)) /= 0) cycle

            value_minus = this%cells(index_minus(1), index_minus(2), index_minus(3))
            value_plus = this%cells(index_plus(1), index_plus(2), index_plus(3))

            call quadratic_extremum_offset( &
                value_minus, value_center, value_plus, extremum_type, offset, valid)
            if (valid) indexes_real(dim) = real(indexes(dim),dp) + offset
        end do
    end subroutine refine_scalar_extremum

    subroutine refine_gradient_extremum( &
        this, bc_ptr, bound, indexes, grad_proj, extremum_type, indexes_real)

        class(field_scalar_cons), intent(in) :: this
        type(boundary_conditions_pointer), intent(in) :: bc_ptr
        integer, dimension(3,2), intent(in) :: bound
        integer, dimension(3), intent(in) :: indexes
        integer, intent(in), optional :: grad_proj
        character(len=*), intent(in) :: extremum_type
        real(dp), dimension(3), intent(inout) :: indexes_real

        integer, dimension(3) :: index_minus, index_plus
        real(dp) :: value_minus, value_center, value_plus, offset
        integer :: dim
        logical :: valid

        value_center = gradient_metric(this, bc_ptr, indexes, grad_proj)

        do dim = 1, this%dimensions
            index_minus = indexes - I_m(dim,:)
            index_plus = indexes + I_m(dim,:)

            if (.not. point_inside_bound(index_minus, bound, this%dimensions)) cycle
            if (.not. point_inside_bound(index_plus, bound, this%dimensions)) cycle
            if (.not. point_inside_array(this, index_minus)) cycle
            if (.not. point_inside_array(this, index_plus)) cycle
            if (bc_ptr%bc_ptr%bc_markers( &
                index_minus(1), index_minus(2), index_minus(3)) /= 0) cycle
            if (bc_ptr%bc_ptr%bc_markers( &
                index_plus(1), index_plus(2), index_plus(3)) /= 0) cycle

            value_minus = gradient_metric(this, bc_ptr, index_minus, grad_proj)
            value_plus = gradient_metric(this, bc_ptr, index_plus, grad_proj)

            call quadratic_extremum_offset( &
                value_minus, value_center, value_plus, extremum_type, offset, valid)
            if (valid) indexes_real(dim) = real(indexes(dim),dp) + offset
        end do
    end subroutine refine_gradient_extremum

    real(dp) function gradient_metric(this, bc_ptr, indexes, grad_proj)
        class(field_scalar_cons), intent(in) :: this
        type(boundary_conditions_pointer), intent(in) :: bc_ptr
        integer, dimension(3), intent(in) :: indexes
        integer, intent(in), optional :: grad_proj

        real(dp) :: grad
        integer :: dim

        if (present(grad_proj)) then
            if (grad_proj >= 1 .and. grad_proj <= this%dimensions) then
                gradient_metric = directional_gradient(this, bc_ptr, indexes, grad_proj)
                return
            end if
        end if

        gradient_metric = 0.0_dp
        do dim = 1, this%dimensions
            grad = directional_gradient(this, bc_ptr, indexes, dim)
            gradient_metric = gradient_metric + grad*grad
        end do
    end function gradient_metric

    real(dp) function directional_gradient(this, bc_ptr, indexes, dim)
        class(field_scalar_cons), intent(in) :: this
        type(boundary_conditions_pointer), intent(in) :: bc_ptr
        integer, dimension(3), intent(in) :: indexes
        integer, intent(in) :: dim

        integer, dimension(3) :: index_minus, index_plus
        real(dp) :: field_minus, field_plus, field_center

        index_minus = indexes - I_m(dim,:)
        index_plus = indexes + I_m(dim,:)
        field_center = this%cells(indexes(1), indexes(2), indexes(3))

        field_minus = field_center
        if (point_inside_array(this, index_minus)) then
            if (bc_ptr%bc_ptr%bc_markers( &
                index_minus(1), index_minus(2), index_minus(3)) == 0) then
                field_minus = this%cells(index_minus(1), index_minus(2), index_minus(3))
            end if
        end if

        field_plus = field_center
        if (point_inside_array(this, index_plus)) then
            if (bc_ptr%bc_ptr%bc_markers( &
                index_plus(1), index_plus(2), index_plus(3)) == 0) then
                field_plus = this%cells(index_plus(1), index_plus(2), index_plus(3))
            end if
        end if

        ! Keep the legacy postprocessor definition: no division by 2*dx.
        directional_gradient = field_plus - field_minus
    end function directional_gradient

    logical function point_inside_bound(indexes, bound, dimensions)
        integer, dimension(3), intent(in) :: indexes
        integer, dimension(3,2), intent(in) :: bound
        integer, intent(in) :: dimensions
        integer :: dim

        point_inside_bound = .true.
        do dim = 1, dimensions
            if (indexes(dim) < bound(dim,1) .or. indexes(dim) > bound(dim,2)) then
                point_inside_bound = .false.
                return
            end if
        end do
    end function point_inside_bound

    logical function point_inside_array(this, indexes)
        class(field_scalar_cons), intent(in) :: this
        integer, dimension(3), intent(in) :: indexes
        integer :: dim

        point_inside_array = .true.
        do dim = 1, 3
            if (indexes(dim) < lbound(this%cells,dim) .or. &
                indexes(dim) > ubound(this%cells,dim)) then
                point_inside_array = .false.
                return
            end if
        end do
    end function point_inside_array

end module field_scalar_class
