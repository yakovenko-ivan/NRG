!> Variable-coefficient finite-volume elliptic solver with geometric multigrid.
!!
!! The module solves the cell-centred problem
!!
!!     A(phi) = s,
!!
!! where A is assembled from non-negative face conductances.  For an active
!! cell P, the integrated discrete equation is
!!
!!     sum_f G_f (phi_P - phi_N) = s_P V_P + b_P.
!!
!! Dirichlet data are imposed at a face half a cell from the cell centre;
!! therefore the boundary coefficient is 2 G_f.  Neumann values are supplied
!! as already integrated outward fluxes and enter b_P directly.
!!
!! The production path uses red-black Gauss-Seidel smoothing, conservative
!! piecewise-constant residual restriction, Kwak prolongation, and a V-cycle.
!! A fully active one-dimensional problem is solved directly by a tridiagonal
!! algorithm.  Pure-Neumann systems are made compatible and fixed to a
!! volume-weighted zero-mean gauge.
module elliptic_multigrid_solver_class

    use kind_parameters, only: dp
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite

    implicit none
    private

    public :: elliptic_multigrid_solver
    public :: elliptic_solver_statistics
    public :: elliptic_boundary_data
    public :: elliptic_bc_internal
    public :: elliptic_bc_neumann
    public :: elliptic_bc_dirichlet

    integer, parameter :: elliptic_bc_internal  = 0
    integer, parameter :: elliptic_bc_neumann   = 1
    integer, parameter :: elliptic_bc_dirichlet = 2

    integer, parameter :: maximum_hierarchy_levels = 12
    integer, parameter :: minimum_coarsenable_extent = 4

    !> Boundary metadata stored on x-, y-, and z-normal faces.
    !!
    !! A Dirichlet value is the prescribed scalar value at the boundary face.
    !! A Neumann value is the integrated outward flux for that face.  Internal
    !! faces normally use elliptic_bc_internal and a zero value.
    type elliptic_boundary_data
        integer, allocatable :: type_x(:,:,:)
        integer, allocatable :: type_y(:,:,:)
        integer, allocatable :: type_z(:,:,:)
        real(dp), allocatable :: value_x(:,:,:)
        real(dp), allocatable :: value_y(:,:,:)
        real(dp), allocatable :: value_z(:,:,:)
    end type elliptic_boundary_data

    !> Convergence and null-space information returned by solve().
    type elliptic_solver_statistics
        integer :: cycles = 0
        integer :: smoothing_iterations = 0
        real(dp) :: initial_residual = 0.0_dp
        real(dp) :: final_residual = 0.0_dp
        real(dp) :: relative_residual = 0.0_dp
        real(dp) :: compatibility_correction = 0.0_dp
        logical :: converged = .false.
    end type elliptic_solver_statistics

    !> One level of the geometric multigrid hierarchy.
    type multigrid_level
        integer :: nx = 1
        integer :: ny = 1
        integer :: nz = 1
        integer :: dimensions = 1
        real(dp), allocatable :: x(:,:,:)
        real(dp), allocatable :: rhs(:,:,:)
        real(dp), allocatable :: residual(:,:,:)
        real(dp), allocatable :: volume(:,:,:)
        real(dp), allocatable :: conductance_x(:,:,:)
        real(dp), allocatable :: conductance_y(:,:,:)
        real(dp), allocatable :: conductance_z(:,:,:)
        logical, allocatable :: active(:,:,:)
        type(elliptic_boundary_data) :: boundary
    end type multigrid_level

    !> User-facing multigrid solver and tuning parameters.
    type elliptic_multigrid_solver
        integer :: pre_smoothing_steps = 2
        integer :: post_smoothing_steps = 2
        integer :: coarse_smoothing_steps = 40
        integer :: maximum_cycles = 20
        real(dp) :: relative_tolerance = 1.0e-4_dp
        real(dp) :: absolute_tolerance = 1.0e-10_dp
        type(multigrid_level), allocatable :: level(:)
    contains
        procedure :: solve => solve_elliptic_problem
        procedure, private :: build_hierarchy
        procedure, private :: load_finest_level
        procedure, private :: coarsen_hierarchy_data
        procedure, private :: v_cycle
        procedure, private :: smooth
        procedure, private :: calculate_residual
        procedure, private :: restrict_residual
        procedure, private :: prolongate_kwak_and_add
    end type elliptic_multigrid_solver

contains

    !==========================================================================
    ! Public solve driver
    !==========================================================================

    !> Solve a one-, two-, or three-dimensional finite-volume elliptic system.
    subroutine solve_elliptic_problem(this, solution, rhs, conductance_x, conductance_y, conductance_z, &
        cell_volume, active, boundary, dimensions, use_initial_guess, statistics)

        class(elliptic_multigrid_solver), intent(inout) :: this
        real(dp), intent(inout) :: solution(:,:,:)
        real(dp), intent(in) :: rhs(:,:,:)
        real(dp), intent(in) :: conductance_x(:,:,:)
        real(dp), intent(in) :: conductance_y(:,:,:)
        real(dp), intent(in) :: conductance_z(:,:,:)
        real(dp), intent(in) :: cell_volume(:,:,:)
        logical, intent(in) :: active(:,:,:)
        type(elliptic_boundary_data), intent(in) :: boundary
        integer, intent(in) :: dimensions
        logical, intent(in), optional :: use_initial_guess
        type(elliptic_solver_statistics), intent(out), optional :: statistics

        type(elliptic_solver_statistics) :: stat
        integer :: cycle
        logical :: warm_start, pure_neumann
        real(dp) :: forcing_norm, tolerance

        stat = elliptic_solver_statistics()

        call validate_problem(solution, rhs, conductance_x, conductance_y, conductance_z, &
            cell_volume, active, boundary, dimensions)
        call validate_solver_controls(this)

        warm_start = .false.
        if (present(use_initial_guess)) warm_start = use_initial_guess
        pure_neumann = .not. has_active_dirichlet_boundary(active, boundary, dimensions)

        ! A direct solve is cheaper and more accurate for a fully active 1-D row.
        if (dimensions == 1 .and. all(active)) then
            call solve_tridiagonal_1d(solution, rhs, conductance_x, cell_volume, boundary, &
                pure_neumann, this%relative_tolerance, this%absolute_tolerance, stat)
            if (present(statistics)) statistics = stat
            return
        end if

        call this%build_hierarchy(size(rhs,1), size(rhs,2), size(rhs,3), dimensions)
        call this%load_finest_level(rhs, conductance_x, conductance_y, conductance_z, &
            cell_volume, active, boundary)
        call this%coarsen_hierarchy_data()

        if (warm_start) then
            this%level(1)%x = solution
        else
            this%level(1)%x = 0.0_dp
        end if
        where (.not. this%level(1)%active) this%level(1)%x = 0.0_dp

        if (pure_neumann) then
            call enforce_neumann_compatibility(this%level(1), stat%compatibility_correction)
            ! Fix the gauge before the first convergence test as well as after
            ! every cycle.  This also gives deterministic early-exit results.
            call remove_volume_weighted_mean(this%level(1))
        end if

        call this%calculate_residual(1)
        stat%initial_residual = active_max_abs(this%level(1)%residual, this%level(1)%active)
        forcing_norm = calculate_forcing_norm(this%level(1))
        tolerance = this%absolute_tolerance + this%relative_tolerance*forcing_norm

        if (stat%initial_residual <= tolerance) then
            stat%converged = .true.
        else
            do cycle = 1, this%maximum_cycles
                call this%v_cycle(1, stat%smoothing_iterations)
                if (pure_neumann) call remove_volume_weighted_mean(this%level(1))

                call this%calculate_residual(1)
                stat%final_residual = active_max_abs(this%level(1)%residual, this%level(1)%active)
                stat%cycles = cycle

                if (stat%final_residual <= tolerance) then
                    stat%converged = .true.
                    exit
                end if
            end do
        end if

        if (stat%cycles == 0) stat%final_residual = stat%initial_residual
        stat%relative_residual = stat%final_residual/max(forcing_norm, tiny(1.0_dp))

        solution = this%level(1)%x
        where (.not. active) solution = 0.0_dp
        if (present(statistics)) statistics = stat
    end subroutine solve_elliptic_problem


    !==========================================================================
    ! Direct one-dimensional path
    !==========================================================================

    !> Solve the fully active one-dimensional operator by Thomas elimination.
    subroutine solve_tridiagonal_1d(solution, rhs, conductance_x, cell_volume, boundary, &
        pure_neumann, relative_tolerance, absolute_tolerance, statistics)

        real(dp), intent(inout) :: solution(:,:,:)
        real(dp), intent(in) :: rhs(:,:,:)
        real(dp), intent(in) :: conductance_x(:,:,:)
        real(dp), intent(in) :: cell_volume(:,:,:)
        type(elliptic_boundary_data), intent(in) :: boundary
        logical, intent(in) :: pure_neumann
        real(dp), intent(in) :: relative_tolerance, absolute_tolerance
        type(elliptic_solver_statistics), intent(out) :: statistics

        real(dp), allocatable :: lower(:), diagonal(:), upper(:)
        real(dp), allocatable :: source(:), compatible_rhs(:)
        real(dp) :: factor, forcing_norm, tolerance
        integer :: n, i

        statistics = elliptic_solver_statistics()
        n = size(rhs,1)
        allocate(lower(n), diagonal(n), upper(n), source(n), compatible_rhs(n))

        lower = 0.0_dp
        diagonal = 0.0_dp
        upper = 0.0_dp
        compatible_rhs = rhs(:,1,1)

        if (pure_neumann) then
            statistics%compatibility_correction = &
                (sum(compatible_rhs*cell_volume(:,1,1)) + boundary_neumann_flux_1d(boundary))/ &
                sum(cell_volume(:,1,1))
            compatible_rhs = compatible_rhs - statistics%compatibility_correction
        end if

        source = compatible_rhs*cell_volume(:,1,1)
        do i = 1, n
            if (i > 1) then
                lower(i) = -conductance_x(i,1,1)
                diagonal(i) = diagonal(i) + conductance_x(i,1,1)
            end if
            if (i < n) then
                upper(i) = -conductance_x(i+1,1,1)
                diagonal(i) = diagonal(i) + conductance_x(i+1,1,1)
            end if
        end do

        call add_boundary_equation_1d(conductance_x(1,1,1), boundary%type_x(1,1,1), &
            boundary%value_x(1,1,1), diagonal(1), source(1))
        call add_boundary_equation_1d(conductance_x(n+1,1,1), boundary%type_x(n+1,1,1), &
            boundary%value_x(n+1,1,1), diagonal(n), source(n))

        forcing_norm = maxval(abs(source/cell_volume(:,1,1)))
        statistics%initial_residual = forcing_norm
        tolerance = absolute_tolerance + relative_tolerance*forcing_norm

        ! A pure-Neumann matrix has a one-dimensional null space.  Pinning one
        ! algebraic row makes elimination nonsingular; the physical zero-mean
        ! gauge is restored after back substitution.
        if (pure_neumann) then
            diagonal(1) = 1.0_dp
            upper(1) = 0.0_dp
            source(1) = 0.0_dp
        end if

        do i = 2, n
            if (abs(diagonal(i-1)) <= tiny(1.0_dp)) &
                error stop 'Elliptic tridiagonal solver: zero pivot'
            factor = lower(i)/diagonal(i-1)
            diagonal(i) = diagonal(i) - factor*upper(i-1)
            source(i) = source(i) - factor*source(i-1)
        end do
        if (abs(diagonal(n)) <= tiny(1.0_dp)) &
            error stop 'Elliptic tridiagonal solver: zero final pivot'

        solution(n,1,1) = source(n)/diagonal(n)
        do i = n-1, 1, -1
            solution(i,1,1) = (source(i)-upper(i)*solution(i+1,1,1))/diagonal(i)
        end do

        if (pure_neumann) then
            solution(:,1,1) = solution(:,1,1) - &
                sum(solution(:,1,1)*cell_volume(:,1,1))/sum(cell_volume(:,1,1))
        end if

        statistics%cycles = 1
        statistics%final_residual = residual_norm_1d(solution(:,1,1), compatible_rhs, &
            conductance_x(:,1,1), cell_volume(:,1,1), boundary)
        statistics%relative_residual = statistics%final_residual/max(forcing_norm, tiny(1.0_dp))
        statistics%converged = statistics%final_residual <= tolerance
    end subroutine solve_tridiagonal_1d


    !> Add one physical boundary contribution to a 1-D matrix row.
    subroutine add_boundary_equation_1d(conductance, boundary_type, boundary_value, diagonal, source)
        real(dp), intent(in) :: conductance, boundary_value
        integer, intent(in) :: boundary_type
        real(dp), intent(inout) :: diagonal, source

        select case (boundary_type)
        case (elliptic_bc_dirichlet)
            diagonal = diagonal + 2.0_dp*conductance
            source = source + 2.0_dp*conductance*boundary_value
        case (elliptic_bc_neumann)
            source = source + boundary_value
        case (elliptic_bc_internal)
            continue
        case default
            error stop 'Elliptic tridiagonal solver: invalid boundary type'
        end select
    end subroutine add_boundary_equation_1d


    !> Infinity norm of the physical 1-D residual after a direct solve.
    real(dp) function residual_norm_1d(x, rhs, conductance, volume, boundary) result(norm)
        real(dp), intent(in) :: x(:), rhs(:), conductance(:), volume(:)
        type(elliptic_boundary_data), intent(in) :: boundary

        real(dp) :: diagonal, neighbour_sum, imposed_source, residual
        integer :: i, n

        n = size(x)
        norm = 0.0_dp
        do i = 1, n
            diagonal = 0.0_dp
            neighbour_sum = 0.0_dp
            imposed_source = 0.0_dp

            if (i == 1) then
                call add_boundary_terms(conductance(1), boundary%type_x(1,1,1), &
                    boundary%value_x(1,1,1), diagonal, imposed_source)
            else
                diagonal = diagonal + conductance(i)
                neighbour_sum = neighbour_sum + conductance(i)*x(max(1,i-1))
            end if

            if (i == n) then
                call add_boundary_terms(conductance(n+1), boundary%type_x(n+1,1,1), &
                    boundary%value_x(n+1,1,1), diagonal, imposed_source)
            else
                diagonal = diagonal + conductance(i+1)
                neighbour_sum = neighbour_sum + conductance(i+1)*x(min(n,i+1))
            end if

            residual = rhs(i) + (imposed_source + neighbour_sum - diagonal*x(i))/volume(i)
            norm = max(norm, abs(residual))
        end do
    end function residual_norm_1d


    !==========================================================================
    ! Hierarchy construction and coefficient coarsening
    !==========================================================================

    !> Allocate or reuse the largest regular 2:1 hierarchy supported by the grid.
    subroutine build_hierarchy(this, nx, ny, nz, dimensions)
        class(elliptic_multigrid_solver), intent(inout) :: this
        integer, intent(in) :: nx, ny, nz, dimensions

        integer :: number_of_levels, level_index
        integer :: level_nx, level_ny, level_nz

        number_of_levels = hierarchy_level_count(nx, ny, nz, dimensions)

        if (allocated(this%level)) then
            if (size(this%level) == number_of_levels .and. &
                this%level(1)%nx == nx .and. this%level(1)%ny == ny .and. &
                this%level(1)%nz == nz .and. this%level(1)%dimensions == dimensions) return
            deallocate(this%level)
        end if

        allocate(this%level(number_of_levels))
        level_nx = nx
        level_ny = ny
        level_nz = nz

        do level_index = 1, number_of_levels
            call allocate_level(this%level(level_index), level_nx, level_ny, level_nz, dimensions)
            if (level_index < number_of_levels) then
                level_nx = max(1, level_nx/2)
                if (dimensions >= 2) level_ny = max(1, level_ny/2)
                if (dimensions >= 3) level_nz = max(1, level_nz/2)
            end if
        end do
    end subroutine build_hierarchy


    !> Copy the current physical operator into the finest reusable level.
    subroutine load_finest_level(this, rhs, conductance_x, conductance_y, conductance_z, &
        cell_volume, active, boundary)

        class(elliptic_multigrid_solver), intent(inout) :: this
        real(dp), intent(in) :: rhs(:,:,:)
        real(dp), intent(in) :: conductance_x(:,:,:), conductance_y(:,:,:), conductance_z(:,:,:)
        real(dp), intent(in) :: cell_volume(:,:,:)
        logical, intent(in) :: active(:,:,:)
        type(elliptic_boundary_data), intent(in) :: boundary

        this%level(1)%rhs = rhs
        this%level(1)%volume = cell_volume
        this%level(1)%active = active
        this%level(1)%conductance_x = conductance_x
        this%level(1)%conductance_y = conductance_y
        this%level(1)%conductance_z = conductance_z
        call copy_boundary_data(this%level(1)%boundary, boundary)
    end subroutine load_finest_level


    !> Return the number of regular 2:1 levels available in every active dimension.
    integer function hierarchy_level_count(nx, ny, nz, dimensions) result(number_of_levels)
        integer, intent(in) :: nx, ny, nz, dimensions

        integer :: level_nx, level_ny, level_nz, minimum_extent

        number_of_levels = 1
        level_nx = nx
        level_ny = ny
        level_nz = nz

        do while (number_of_levels < maximum_hierarchy_levels)
            minimum_extent = level_nx
            if (dimensions >= 2) minimum_extent = min(minimum_extent, level_ny)
            if (dimensions >= 3) minimum_extent = min(minimum_extent, level_nz)
            if (minimum_extent < minimum_coarsenable_extent) exit

            if (mod(level_nx,2) /= 0) exit
            if (dimensions >= 2 .and. mod(level_ny,2) /= 0) exit
            if (dimensions >= 3 .and. mod(level_nz,2) /= 0) exit

            level_nx = max(1, level_nx/2)
            if (dimensions >= 2) level_ny = max(1, level_ny/2)
            if (dimensions >= 3) level_nz = max(1, level_nz/2)
            number_of_levels = number_of_levels + 1
        end do
    end function hierarchy_level_count


    !> Allocate all fields of one multigrid level.
    subroutine allocate_level(level, nx, ny, nz, dimensions)
        type(multigrid_level), intent(inout) :: level
        integer, intent(in) :: nx, ny, nz, dimensions

        level%nx = nx
        level%ny = ny
        level%nz = nz
        level%dimensions = dimensions

        allocate(level%x(nx,ny,nz), level%rhs(nx,ny,nz), level%residual(nx,ny,nz))
        allocate(level%volume(nx,ny,nz), level%active(nx,ny,nz))
        allocate(level%conductance_x(nx+1,ny,nz))
        allocate(level%conductance_y(nx,ny+1,nz))
        allocate(level%conductance_z(nx,ny,nz+1))
        call allocate_boundary_data(level%boundary, nx, ny, nz)

        level%x = 0.0_dp
        level%rhs = 0.0_dp
        level%residual = 0.0_dp
        level%volume = 1.0_dp
        level%active = .true.
        level%conductance_x = 0.0_dp
        level%conductance_y = 0.0_dp
        level%conductance_z = 0.0_dp
    end subroutine allocate_level


    !> Allocate face boundary arrays for one grid level.
    subroutine allocate_boundary_data(boundary, nx, ny, nz)
        type(elliptic_boundary_data), intent(inout) :: boundary
        integer, intent(in) :: nx, ny, nz

        allocate(boundary%type_x(nx+1,ny,nz), boundary%value_x(nx+1,ny,nz))
        allocate(boundary%type_y(nx,ny+1,nz), boundary%value_y(nx,ny+1,nz))
        allocate(boundary%type_z(nx,ny,nz+1), boundary%value_z(nx,ny,nz+1))

        boundary%type_x = elliptic_bc_internal
        boundary%type_y = elliptic_bc_internal
        boundary%type_z = elliptic_bc_internal
        boundary%value_x = 0.0_dp
        boundary%value_y = 0.0_dp
        boundary%value_z = 0.0_dp
    end subroutine allocate_boundary_data


    !> Copy caller-owned boundary metadata into an allocated hierarchy level.
    subroutine copy_boundary_data(destination, source)
        type(elliptic_boundary_data), intent(inout) :: destination
        type(elliptic_boundary_data), intent(in) :: source

        destination%type_x = source%type_x
        destination%type_y = source%type_y
        destination%type_z = source%type_z
        destination%value_x = source%value_x
        destination%value_y = source%value_y
        destination%value_z = source%value_z
    end subroutine copy_boundary_data


    !> Coarsen active masks, volumes, conductances, and boundary types.
    !!
    !! Coarse levels solve error equations, so physical non-zero boundary values
    !! are intentionally replaced by homogeneous values.  Conductances are
    !! summed over parallel fine faces and divided by two for the doubled normal
    !! distance of the rediscretized coarse operator.
    subroutine coarsen_hierarchy_data(this)
        class(elliptic_multigrid_solver), intent(inout) :: this

        integer :: level_index, i, j, k
        integer :: fine_i, fine_j, fine_k, ii, jj, kk

        do level_index = 2, size(this%level)
            this%level(level_index)%active = .false.
            this%level(level_index)%volume = 0.0_dp

            do k = 1, this%level(level_index)%nz
                do j = 1, this%level(level_index)%ny
                    do i = 1, this%level(level_index)%nx
                        ii = 2*i - 1
                        jj = 2*j - 1
                        kk = 2*k - 1
                        do fine_k = kk, min(kk+1, this%level(level_index-1)%nz)
                            do fine_j = jj, min(jj+1, this%level(level_index-1)%ny)
                                do fine_i = ii, min(ii+1, this%level(level_index-1)%nx)
                                    if (.not. this%level(level_index-1)%active(fine_i,fine_j,fine_k)) cycle
                                    this%level(level_index)%active(i,j,k) = .true.
                                    this%level(level_index)%volume(i,j,k) = &
                                        this%level(level_index)%volume(i,j,k) + &
                                        this%level(level_index-1)%volume(fine_i,fine_j,fine_k)
                                end do
                            end do
                        end do
                    end do
                end do
            end do

            call coarsen_x_faces(this%level(level_index-1), this%level(level_index))
            call coarsen_y_faces(this%level(level_index-1), this%level(level_index))
            call coarsen_z_faces(this%level(level_index-1), this%level(level_index))
        end do
    end subroutine coarsen_hierarchy_data


    !> Rediscretize x-normal face conductances and homogeneous error BCs.
    subroutine coarsen_x_faces(fine, coarse)
        type(multigrid_level), intent(in) :: fine
        type(multigrid_level), intent(inout) :: coarse

        integer :: i, j, k, fine_i, fine_j, fine_k, j_offset, k_offset, boundary_type
        real(dp) :: internal_conductance, boundary_conductance
        logical :: has_dirichlet

        coarse%conductance_x = 0.0_dp
        coarse%boundary%type_x = elliptic_bc_internal
        coarse%boundary%value_x = 0.0_dp

        do k = 1, coarse%nz
            do j = 1, coarse%ny
                do i = 1, coarse%nx+1
                    fine_i = min(2*i-1, fine%nx+1)
                    fine_j = 2*j-1
                    fine_k = 2*k-1
                    internal_conductance = 0.0_dp
                    boundary_conductance = 0.0_dp
                    has_dirichlet = .false.

                    do k_offset = fine_k, min(fine_k+1, fine%nz)
                        do j_offset = fine_j, min(fine_j+1, fine%ny)
                            boundary_type = fine%boundary%type_x(fine_i,j_offset,k_offset)
                            if (boundary_type == elliptic_bc_internal) then
                                internal_conductance = internal_conductance + &
                                    fine%conductance_x(fine_i,j_offset,k_offset)
                            else
                                boundary_conductance = boundary_conductance + &
                                    fine%conductance_x(fine_i,j_offset,k_offset)
                                has_dirichlet = has_dirichlet .or. &
                                    boundary_type == elliptic_bc_dirichlet
                            end if
                        end do
                    end do

                    call set_coarse_face(internal_conductance, boundary_conductance, has_dirichlet, &
                        coarse%conductance_x(i,j,k), coarse%boundary%type_x(i,j,k))
                end do
            end do
        end do
    end subroutine coarsen_x_faces


    !> Rediscretize y-normal face conductances and homogeneous error BCs.
    subroutine coarsen_y_faces(fine, coarse)
        type(multigrid_level), intent(in) :: fine
        type(multigrid_level), intent(inout) :: coarse

        integer :: i, j, k, fine_i, fine_j, fine_k, i_offset, k_offset, boundary_type
        real(dp) :: internal_conductance, boundary_conductance
        logical :: has_dirichlet

        coarse%conductance_y = 0.0_dp
        coarse%boundary%type_y = elliptic_bc_internal
        coarse%boundary%value_y = 0.0_dp
        if (coarse%dimensions < 2) return

        do k = 1, coarse%nz
            do j = 1, coarse%ny+1
                do i = 1, coarse%nx
                    fine_i = 2*i-1
                    fine_j = min(2*j-1, fine%ny+1)
                    fine_k = 2*k-1
                    internal_conductance = 0.0_dp
                    boundary_conductance = 0.0_dp
                    has_dirichlet = .false.

                    do k_offset = fine_k, min(fine_k+1, fine%nz)
                        do i_offset = fine_i, min(fine_i+1, fine%nx)
                            boundary_type = fine%boundary%type_y(i_offset,fine_j,k_offset)
                            if (boundary_type == elliptic_bc_internal) then
                                internal_conductance = internal_conductance + &
                                    fine%conductance_y(i_offset,fine_j,k_offset)
                            else
                                boundary_conductance = boundary_conductance + &
                                    fine%conductance_y(i_offset,fine_j,k_offset)
                                has_dirichlet = has_dirichlet .or. &
                                    boundary_type == elliptic_bc_dirichlet
                            end if
                        end do
                    end do

                    call set_coarse_face(internal_conductance, boundary_conductance, has_dirichlet, &
                        coarse%conductance_y(i,j,k), coarse%boundary%type_y(i,j,k))
                end do
            end do
        end do
    end subroutine coarsen_y_faces


    !> Rediscretize z-normal face conductances and homogeneous error BCs.
    subroutine coarsen_z_faces(fine, coarse)
        type(multigrid_level), intent(in) :: fine
        type(multigrid_level), intent(inout) :: coarse

        integer :: i, j, k, fine_i, fine_j, fine_k, i_offset, j_offset, boundary_type
        real(dp) :: internal_conductance, boundary_conductance
        logical :: has_dirichlet

        coarse%conductance_z = 0.0_dp
        coarse%boundary%type_z = elliptic_bc_internal
        coarse%boundary%value_z = 0.0_dp
        if (coarse%dimensions < 3) return

        do k = 1, coarse%nz+1
            do j = 1, coarse%ny
                do i = 1, coarse%nx
                    fine_i = 2*i-1
                    fine_j = 2*j-1
                    fine_k = min(2*k-1, fine%nz+1)
                    internal_conductance = 0.0_dp
                    boundary_conductance = 0.0_dp
                    has_dirichlet = .false.

                    do j_offset = fine_j, min(fine_j+1, fine%ny)
                        do i_offset = fine_i, min(fine_i+1, fine%nx)
                            boundary_type = fine%boundary%type_z(i_offset,j_offset,fine_k)
                            if (boundary_type == elliptic_bc_internal) then
                                internal_conductance = internal_conductance + &
                                    fine%conductance_z(i_offset,j_offset,fine_k)
                            else
                                boundary_conductance = boundary_conductance + &
                                    fine%conductance_z(i_offset,j_offset,fine_k)
                                has_dirichlet = has_dirichlet .or. &
                                    boundary_type == elliptic_bc_dirichlet
                            end if
                        end do
                    end do

                    call set_coarse_face(internal_conductance, boundary_conductance, has_dirichlet, &
                        coarse%conductance_z(i,j,k), coarse%boundary%type_z(i,j,k))
                end do
            end do
        end do
    end subroutine coarsen_z_faces


    !> Select the rediscretized conductance and boundary type of one coarse face.
    subroutine set_coarse_face(internal_conductance, boundary_conductance, has_dirichlet, &
        coarse_conductance, coarse_boundary_type)

        real(dp), intent(in) :: internal_conductance, boundary_conductance
        logical, intent(in) :: has_dirichlet
        real(dp), intent(out) :: coarse_conductance
        integer, intent(out) :: coarse_boundary_type

        if (internal_conductance > tiny(1.0_dp)) then
            coarse_conductance = 0.5_dp*internal_conductance
            coarse_boundary_type = elliptic_bc_internal
        else
            coarse_conductance = 0.5_dp*boundary_conductance
            if (boundary_conductance > tiny(1.0_dp)) then
                coarse_boundary_type = merge(elliptic_bc_dirichlet, elliptic_bc_neumann, has_dirichlet)
            else
                coarse_boundary_type = elliptic_bc_internal
            end if
        end if
    end subroutine set_coarse_face


    !==========================================================================
    ! V-cycle, smoother, residual, and transfer operators
    !==========================================================================

    !> Apply one recursive V-cycle to level_index.
    recursive subroutine v_cycle(this, level_index, smoothing_count)
        class(elliptic_multigrid_solver), intent(inout) :: this
        integer, intent(in) :: level_index
        integer, intent(inout) :: smoothing_count

        if (level_index == size(this%level)) then
            call this%smooth(level_index, this%coarse_smoothing_steps)
            smoothing_count = smoothing_count + this%coarse_smoothing_steps
            return
        end if

        call this%smooth(level_index, this%pre_smoothing_steps)
        smoothing_count = smoothing_count + this%pre_smoothing_steps

        call this%calculate_residual(level_index)
        call this%restrict_residual(level_index)
        this%level(level_index+1)%x = 0.0_dp

        call this%v_cycle(level_index+1, smoothing_count)
        call this%prolongate_kwak_and_add(level_index)

        call this%smooth(level_index, this%post_smoothing_steps)
        smoothing_count = smoothing_count + this%post_smoothing_steps
    end subroutine v_cycle


    !> Red-black Gauss-Seidel smoothing on one hierarchy level.
    subroutine smooth(this, level_index, iterations)
        class(elliptic_multigrid_solver), intent(inout) :: this
        integer, intent(in) :: level_index, iterations

        integer :: iteration, color, i, j, k
        real(dp) :: diagonal, neighbour_sum, imposed_source, algebraic_residual

        do iteration = 1, iterations
            do color = 0, 1
                !$omp parallel do collapse(3) schedule(static) &
                !$omp& private(i,j,k,diagonal,neighbour_sum,imposed_source,algebraic_residual)
                do k = 1, this%level(level_index)%nz
                    do j = 1, this%level(level_index)%ny
                        do i = 1, this%level(level_index)%nx
                            if (mod(i+j+k,2) /= color) cycle
                            if (.not. this%level(level_index)%active(i,j,k)) cycle

                            call stencil_terms(this%level(level_index), i, j, k, &
                                diagonal, neighbour_sum, imposed_source)
                            if (diagonal <= tiny(1.0_dp)) cycle

                            algebraic_residual = &
                                this%level(level_index)%rhs(i,j,k)*this%level(level_index)%volume(i,j,k) + &
                                imposed_source + neighbour_sum - &
                                diagonal*this%level(level_index)%x(i,j,k)
                            this%level(level_index)%x(i,j,k) = &
                                this%level(level_index)%x(i,j,k) + algebraic_residual/diagonal
                        end do
                    end do
                end do
                !$omp end parallel do
            end do
        end do
    end subroutine smooth


    !> Compute the cell-volume-normalized residual on one level.
    subroutine calculate_residual(this, level_index)
        class(elliptic_multigrid_solver), intent(inout) :: this
        integer, intent(in) :: level_index

        integer :: i, j, k
        real(dp) :: diagonal, neighbour_sum, imposed_source, integrated_residual

        this%level(level_index)%residual = 0.0_dp
        !$omp parallel do collapse(3) schedule(static) &
        !$omp& private(i,j,k,diagonal,neighbour_sum,imposed_source,integrated_residual)
        do k = 1, this%level(level_index)%nz
            do j = 1, this%level(level_index)%ny
                do i = 1, this%level(level_index)%nx
                    if (.not. this%level(level_index)%active(i,j,k)) cycle

                    call stencil_terms(this%level(level_index), i, j, k, &
                        diagonal, neighbour_sum, imposed_source)
                    integrated_residual = &
                        this%level(level_index)%rhs(i,j,k)*this%level(level_index)%volume(i,j,k) + &
                        imposed_source + neighbour_sum - &
                        diagonal*this%level(level_index)%x(i,j,k)
                    this%level(level_index)%residual(i,j,k) = &
                        integrated_residual/this%level(level_index)%volume(i,j,k)
                end do
            end do
        end do
        !$omp end parallel do
    end subroutine calculate_residual


    !> Assemble the diagonal, active-neighbour term, and imposed boundary source.
    subroutine stencil_terms(level, i, j, k, diagonal, neighbour_sum, imposed_source)
        type(multigrid_level), intent(in) :: level
        integer, intent(in) :: i, j, k
        real(dp), intent(out) :: diagonal, neighbour_sum, imposed_source

        diagonal = 0.0_dp
        neighbour_sum = 0.0_dp
        imposed_source = 0.0_dp

        call add_face(i-1, j, k, level%conductance_x(i,j,k), &
            level%boundary%type_x(i,j,k), level%boundary%value_x(i,j,k))
        call add_face(i+1, j, k, level%conductance_x(i+1,j,k), &
            level%boundary%type_x(i+1,j,k), level%boundary%value_x(i+1,j,k))

        if (level%dimensions >= 2) then
            call add_face(i, j-1, k, level%conductance_y(i,j,k), &
                level%boundary%type_y(i,j,k), level%boundary%value_y(i,j,k))
            call add_face(i, j+1, k, level%conductance_y(i,j+1,k), &
                level%boundary%type_y(i,j+1,k), level%boundary%value_y(i,j+1,k))
        end if

        if (level%dimensions >= 3) then
            call add_face(i, j, k-1, level%conductance_z(i,j,k), &
                level%boundary%type_z(i,j,k), level%boundary%value_z(i,j,k))
            call add_face(i, j, k+1, level%conductance_z(i,j,k+1), &
                level%boundary%type_z(i,j,k+1), level%boundary%value_z(i,j,k+1))
        end if

    contains

        subroutine add_face(neighbour_i, neighbour_j, neighbour_k, conductance, &
            boundary_type, boundary_value)

            integer, intent(in) :: neighbour_i, neighbour_j, neighbour_k, boundary_type
            real(dp), intent(in) :: conductance, boundary_value
            logical :: neighbour_is_active

            neighbour_is_active = neighbour_i >= 1 .and. neighbour_i <= level%nx .and. &
                neighbour_j >= 1 .and. neighbour_j <= level%ny .and. &
                neighbour_k >= 1 .and. neighbour_k <= level%nz
            if (neighbour_is_active) &
                neighbour_is_active = level%active(neighbour_i,neighbour_j,neighbour_k)

            if (neighbour_is_active) then
                diagonal = diagonal + conductance
                neighbour_sum = neighbour_sum + &
                    conductance*level%x(neighbour_i,neighbour_j,neighbour_k)
            else
                call add_boundary_terms(conductance, boundary_type, boundary_value, &
                    diagonal, imposed_source)
            end if
        end subroutine add_face

    end subroutine stencil_terms


    !> Add one Dirichlet or Neumann contribution to a cell equation.
    subroutine add_boundary_terms(conductance, boundary_type, boundary_value, diagonal, imposed_source)
        real(dp), intent(in) :: conductance, boundary_value
        integer, intent(in) :: boundary_type
        real(dp), intent(inout) :: diagonal, imposed_source

        select case (boundary_type)
        case (elliptic_bc_dirichlet)
            diagonal = diagonal + 2.0_dp*conductance
            imposed_source = imposed_source + 2.0_dp*conductance*boundary_value
        case (elliptic_bc_neumann)
            imposed_source = imposed_source + boundary_value
        case (elliptic_bc_internal)
            continue
        case default
            error stop 'Elliptic solver: invalid boundary type in stencil'
        end select
    end subroutine add_boundary_terms


    !> Restrict the fine residual as a conservative piecewise-constant source.
    !!
    !! The residual is stored per unit cell volume.  Therefore the coarse value
    !! is the sum of integrated child residuals divided by the active coarse
    !! volume.  On a uniform Cartesian grid this reduces to the arithmetic
    !! average used by the original implementation.
    subroutine restrict_residual(this, level_index)
        class(elliptic_multigrid_solver), intent(inout) :: this
        integer, intent(in) :: level_index

        integer :: i, j, k, fine_i, fine_j, fine_k, ii, jj, kk
        real(dp) :: integrated_residual, active_volume

        this%level(level_index+1)%rhs = 0.0_dp
        !$omp parallel do collapse(3) schedule(static) &
        !$omp& private(i,j,k,fine_i,fine_j,fine_k,ii,jj,kk,integrated_residual,active_volume)
        do k = 1, this%level(level_index+1)%nz
            do j = 1, this%level(level_index+1)%ny
                do i = 1, this%level(level_index+1)%nx
                    ii = 2*i - 1
                    jj = 2*j - 1
                    kk = 2*k - 1
                    integrated_residual = 0.0_dp
                    active_volume = 0.0_dp

                    do fine_k = kk, min(kk+1, this%level(level_index)%nz)
                        do fine_j = jj, min(jj+1, this%level(level_index)%ny)
                            do fine_i = ii, min(ii+1, this%level(level_index)%nx)
                                if (.not. this%level(level_index)%active(fine_i,fine_j,fine_k)) cycle
                                integrated_residual = integrated_residual + &
                                    this%level(level_index)%residual(fine_i,fine_j,fine_k)* &
                                    this%level(level_index)%volume(fine_i,fine_j,fine_k)
                                active_volume = active_volume + &
                                    this%level(level_index)%volume(fine_i,fine_j,fine_k)
                            end do
                        end do
                    end do

                    if (active_volume > 0.0_dp) then
                        this%level(level_index+1)%rhs(i,j,k) = integrated_residual/active_volume
                    end if
                end do
            end do
        end do
        !$omp end parallel do
    end subroutine restrict_residual


    !> Prolongate a coarse correction by the Kwak stencil and add it in place.
    subroutine prolongate_kwak_and_add(this, level_index)
        class(elliptic_multigrid_solver), intent(inout) :: this
        integer, intent(in) :: level_index

        integer :: i, j, k, coarse_i, coarse_j, coarse_k
        integer :: x_direction, y_direction, z_direction
        real(dp) :: center_value, correction

        !$omp parallel do collapse(3) schedule(static) &
        !$omp& private(i,j,k,coarse_i,coarse_j,coarse_k,x_direction,y_direction,z_direction, &
        !$omp& center_value,correction)
        do k = 1, this%level(level_index)%nz
            do j = 1, this%level(level_index)%ny
                do i = 1, this%level(level_index)%nx
                    if (.not. this%level(level_index)%active(i,j,k)) cycle

                    coarse_i = min((i+1)/2, this%level(level_index+1)%nx)
                    coarse_j = min((j+1)/2, this%level(level_index+1)%ny)
                    coarse_k = min((k+1)/2, this%level(level_index+1)%nz)
                    center_value = this%level(level_index+1)%x(coarse_i,coarse_j,coarse_k)

                    select case (this%level(level_index)%dimensions)
                    case (1)
                        correction = center_value
                    case (2)
                        x_direction = merge(-1, 1, mod(i,2) == 1)
                        y_direction = merge(-1, 1, mod(j,2) == 1)
                        correction = 0.5_dp*center_value + 0.25_dp*( &
                            coarse_neighbour_value(this%level(level_index+1), &
                                coarse_i,coarse_j,coarse_k,1,x_direction) + &
                            coarse_neighbour_value(this%level(level_index+1), &
                                coarse_i,coarse_j,coarse_k,2,y_direction))
                    case (3)
                        x_direction = merge(-1, 1, mod(i,2) == 1)
                        y_direction = merge(-1, 1, mod(j,2) == 1)
                        z_direction = merge(-1, 1, mod(k,2) == 1)
                        correction = 0.5_dp*center_value + ( &
                            coarse_neighbour_value(this%level(level_index+1), &
                                coarse_i,coarse_j,coarse_k,1,x_direction) + &
                            coarse_neighbour_value(this%level(level_index+1), &
                                coarse_i,coarse_j,coarse_k,2,y_direction) + &
                            coarse_neighbour_value(this%level(level_index+1), &
                                coarse_i,coarse_j,coarse_k,3,z_direction))/6.0_dp
                    end select

                    this%level(level_index)%x(i,j,k) = &
                        this%level(level_index)%x(i,j,k) + correction
                end do
            end do
        end do
        !$omp end parallel do
    end subroutine prolongate_kwak_and_add


    !> Return a neighbouring coarse correction or a homogeneous boundary ghost.
    real(dp) function coarse_neighbour_value(level, i, j, k, dimension, direction) result(value)
        type(multigrid_level), intent(in) :: level
        integer, intent(in) :: i, j, k, dimension, direction

        integer :: neighbour_i, neighbour_j, neighbour_k, boundary_type
        real(dp) :: center_value

        neighbour_i = i
        neighbour_j = j
        neighbour_k = k
        select case (dimension)
        case (1)
            neighbour_i = i + direction
        case (2)
            neighbour_j = j + direction
        case (3)
            neighbour_k = k + direction
        end select

        center_value = level%x(i,j,k)
        if (neighbour_i >= 1 .and. neighbour_i <= level%nx .and. &
            neighbour_j >= 1 .and. neighbour_j <= level%ny .and. &
            neighbour_k >= 1 .and. neighbour_k <= level%nz) then
            if (level%active(neighbour_i,neighbour_j,neighbour_k)) then
                value = level%x(neighbour_i,neighbour_j,neighbour_k)
                return
            end if
        end if

        boundary_type = face_boundary_type(level, i, j, k, dimension, direction)
        if (boundary_type == elliptic_bc_dirichlet) then
            ! Homogeneous Dirichlet correction at a halfway boundary.
            value = -center_value
        else
            ! Homogeneous Neumann correction.
            value = center_value
        end if
    end function coarse_neighbour_value


    !> Boundary type of the requested face; an unspecified face is treated as Neumann.
    integer function face_boundary_type(level, i, j, k, dimension, direction) result(boundary_type)
        type(multigrid_level), intent(in) :: level
        integer, intent(in) :: i, j, k, dimension, direction

        select case (dimension)
        case (1)
            boundary_type = level%boundary%type_x(i+merge(0,1,direction<0),j,k)
        case (2)
            boundary_type = level%boundary%type_y(i,j+merge(0,1,direction<0),k)
        case (3)
            boundary_type = level%boundary%type_z(i,j,k+merge(0,1,direction<0))
        case default
            boundary_type = elliptic_bc_neumann
        end select

        if (boundary_type == elliptic_bc_internal) boundary_type = elliptic_bc_neumann
    end function face_boundary_type


    !==========================================================================
    ! Compatibility, norms, and boundary topology
    !==========================================================================

    !> Subtract the source mean required by a pure-Neumann compatibility condition.
    subroutine enforce_neumann_compatibility(level, correction)
        type(multigrid_level), intent(inout) :: level
        real(dp), intent(out) :: correction

        real(dp) :: integrated_source, total_volume

        total_volume = sum(level%volume, mask=level%active)
        integrated_source = sum(level%rhs*level%volume, mask=level%active) + &
            total_active_neumann_flux(level)
        correction = integrated_source/total_volume
        where (level%active) level%rhs = level%rhs - correction
    end subroutine enforce_neumann_compatibility


    !> Remove the volume-weighted constant null-space mode from a Neumann solution.
    subroutine remove_volume_weighted_mean(level)
        type(multigrid_level), intent(inout) :: level

        real(dp) :: mean_value, total_volume

        total_volume = sum(level%volume, mask=level%active)
        mean_value = sum(level%x*level%volume, mask=level%active)/total_volume
        where (level%active) level%x = level%x - mean_value
    end subroutine remove_volume_weighted_mean


    !> Norm of the imposed differential source, including physical boundary data.
    real(dp) function calculate_forcing_norm(level) result(norm)
        type(multigrid_level), intent(in) :: level

        integer :: i, j, k
        real(dp) :: diagonal, neighbour_sum, imposed_source, local_forcing

        norm = 0.0_dp
        do k = 1, level%nz
            do j = 1, level%ny
                do i = 1, level%nx
                    if (.not. level%active(i,j,k)) cycle
                    call stencil_terms(level, i, j, k, diagonal, neighbour_sum, imposed_source)
                    ! Remove the current-solution neighbour term: only rhs and
                    ! prescribed boundary data define the convergence scale.
                    local_forcing = abs(level%rhs(i,j,k) + imposed_source/level%volume(i,j,k))
                    norm = max(norm, local_forcing)
                end do
            end do
        end do
        norm = max(norm, tiny(1.0_dp))
    end function calculate_forcing_norm


    !> Maximum absolute value over active cells.
    real(dp) function active_max_abs(field, active) result(norm)
        real(dp), intent(in) :: field(:,:,:)
        logical, intent(in) :: active(:,:,:)

        norm = maxval(abs(field), mask=active)
    end function active_max_abs


    !> True when at least one Dirichlet face borders an active cell.
    logical function has_active_dirichlet_boundary(active, boundary, dimensions) result(found)
        logical, intent(in) :: active(:,:,:)
        type(elliptic_boundary_data), intent(in) :: boundary
        integer, intent(in) :: dimensions

        integer :: i, j, k, nx, ny, nz

        nx = size(active,1)
        ny = size(active,2)
        nz = size(active,3)
        found = .false.

        do k = 1, nz
            do j = 1, ny
                do i = 1, nx+1
                    if (boundary%type_x(i,j,k) /= elliptic_bc_dirichlet) cycle
                    if (face_touches_active_x(active,i,j,k)) then
                        found = .true.
                        return
                    end if
                end do
            end do
        end do

        if (dimensions >= 2) then
            do k = 1, nz
                do j = 1, ny+1
                    do i = 1, nx
                        if (boundary%type_y(i,j,k) /= elliptic_bc_dirichlet) cycle
                        if (face_touches_active_y(active,i,j,k)) then
                            found = .true.
                            return
                        end if
                    end do
                end do
            end do
        end if

        if (dimensions >= 3) then
            do k = 1, nz+1
                do j = 1, ny
                    do i = 1, nx
                        if (boundary%type_z(i,j,k) /= elliptic_bc_dirichlet) cycle
                        if (face_touches_active_z(active,i,j,k)) then
                            found = .true.
                            return
                        end if
                    end do
                end do
            end do
        end if
    end function has_active_dirichlet_boundary


    !> Sum prescribed Neumann fluxes only on faces adjacent to active cells.
    real(dp) function total_active_neumann_flux(level) result(total_flux)
        type(multigrid_level), intent(in) :: level

        integer :: i, j, k

        total_flux = 0.0_dp
        do k = 1, level%nz
            do j = 1, level%ny
                do i = 1, level%nx+1
                    if (level%boundary%type_x(i,j,k) == elliptic_bc_neumann .and. &
                        face_touches_active_x(level%active,i,j,k)) then
                        total_flux = total_flux + level%boundary%value_x(i,j,k)
                    end if
                end do
            end do
        end do

        if (level%dimensions >= 2) then
            do k = 1, level%nz
                do j = 1, level%ny+1
                    do i = 1, level%nx
                        if (level%boundary%type_y(i,j,k) == elliptic_bc_neumann .and. &
                            face_touches_active_y(level%active,i,j,k)) then
                            total_flux = total_flux + level%boundary%value_y(i,j,k)
                        end if
                    end do
                end do
            end do
        end if

        if (level%dimensions >= 3) then
            do k = 1, level%nz+1
                do j = 1, level%ny
                    do i = 1, level%nx
                        if (level%boundary%type_z(i,j,k) == elliptic_bc_neumann .and. &
                            face_touches_active_z(level%active,i,j,k)) then
                            total_flux = total_flux + level%boundary%value_z(i,j,k)
                        end if
                    end do
                end do
            end do
        end if
    end function total_active_neumann_flux


    !> Sum the two endpoint Neumann fluxes in the direct 1-D path.
    real(dp) function boundary_neumann_flux_1d(boundary) result(total_flux)
        type(elliptic_boundary_data), intent(in) :: boundary

        integer :: upper_face

        upper_face = size(boundary%type_x,1)
        total_flux = 0.0_dp
        if (boundary%type_x(1,1,1) == elliptic_bc_neumann) &
            total_flux = total_flux + boundary%value_x(1,1,1)
        if (boundary%type_x(upper_face,1,1) == elliptic_bc_neumann) &
            total_flux = total_flux + boundary%value_x(upper_face,1,1)
    end function boundary_neumann_flux_1d


    logical function face_touches_active_x(active, i, j, k) result(touches)
        logical, intent(in) :: active(:,:,:)
        integer, intent(in) :: i, j, k

        touches = .false.
        if (i > 1) touches = touches .or. active(i-1,j,k)
        if (i <= size(active,1)) touches = touches .or. active(i,j,k)
    end function face_touches_active_x


    logical function face_touches_active_y(active, i, j, k) result(touches)
        logical, intent(in) :: active(:,:,:)
        integer, intent(in) :: i, j, k

        touches = .false.
        if (j > 1) touches = touches .or. active(i,j-1,k)
        if (j <= size(active,2)) touches = touches .or. active(i,j,k)
    end function face_touches_active_y


    logical function face_touches_active_z(active, i, j, k) result(touches)
        logical, intent(in) :: active(:,:,:)
        integer, intent(in) :: i, j, k

        touches = .false.
        if (k > 1) touches = touches .or. active(i,j,k-1)
        if (k <= size(active,3)) touches = touches .or. active(i,j,k)
    end function face_touches_active_z


    !==========================================================================
    ! Input validation
    !==========================================================================

    !> Validate solver tuning parameters before entering an iterative cycle.
    subroutine validate_solver_controls(this)
        class(elliptic_multigrid_solver), intent(in) :: this

        if (this%pre_smoothing_steps < 0 .or. this%post_smoothing_steps < 0) &
            error stop 'Elliptic solver: smoothing counts must be non-negative'
        if (this%coarse_smoothing_steps < 1) &
            error stop 'Elliptic solver: coarse smoothing count must be positive'
        if (this%maximum_cycles < 1) &
            error stop 'Elliptic solver: maximum cycle count must be positive'
        if (this%relative_tolerance < 0.0_dp .or. this%absolute_tolerance < 0.0_dp) &
            error stop 'Elliptic solver: tolerances must be non-negative'
    end subroutine validate_solver_controls


    !> Validate all cell, face, boundary, and dimensional array contracts.
    subroutine validate_problem(solution, rhs, conductance_x, conductance_y, conductance_z, &
        cell_volume, active, boundary, dimensions)

        real(dp), intent(in) :: solution(:,:,:), rhs(:,:,:)
        real(dp), intent(in) :: conductance_x(:,:,:), conductance_y(:,:,:), conductance_z(:,:,:)
        real(dp), intent(in) :: cell_volume(:,:,:)
        logical, intent(in) :: active(:,:,:)
        type(elliptic_boundary_data), intent(in) :: boundary
        integer, intent(in) :: dimensions

        integer :: nx, ny, nz

        nx = size(rhs,1)
        ny = size(rhs,2)
        nz = size(rhs,3)

        if (dimensions < 1 .or. dimensions > 3) &
            error stop 'Elliptic solver: dimensions must be 1, 2 or 3'
        if (dimensions == 1 .and. (ny /= 1 .or. nz /= 1)) &
            error stop 'Elliptic solver: a 1-D problem requires ny=nz=1'
        if (dimensions == 2 .and. nz /= 1) &
            error stop 'Elliptic solver: a 2-D problem requires nz=1'
        if (.not. any(active)) &
            error stop 'Elliptic solver: the domain contains no active cells'

        if (any(shape(solution) /= shape(rhs)) .or. &
            any(shape(cell_volume) /= shape(rhs)) .or. any(shape(active) /= shape(rhs))) &
            error stop 'Elliptic solver: incompatible cell-array shapes'
        if (any(shape(conductance_x) /= [nx+1,ny,nz])) &
            error stop 'Elliptic solver: invalid x-face conductance shape'
        if (any(shape(conductance_y) /= [nx,ny+1,nz])) &
            error stop 'Elliptic solver: invalid y-face conductance shape'
        if (any(shape(conductance_z) /= [nx,ny,nz+1])) &
            error stop 'Elliptic solver: invalid z-face conductance shape'

        call validate_boundary_shapes(boundary, nx, ny, nz)
        call validate_boundary_types(boundary, dimensions)

        if (any(cell_volume <= 0.0_dp .and. active)) &
            error stop 'Elliptic solver: non-positive active-cell volume'
        if (any(conductance_x < 0.0_dp) .or. any(conductance_y < 0.0_dp) .or. &
            any(conductance_z < 0.0_dp)) &
            error stop 'Elliptic solver: negative face conductance'

        if (.not. all(ieee_is_finite(solution)) .or. .not. all(ieee_is_finite(rhs)) .or. &
            .not. all(ieee_is_finite(cell_volume))) &
            error stop 'Elliptic solver: non-finite cell data'
        if (.not. all(ieee_is_finite(conductance_x)) .or. &
            .not. all(ieee_is_finite(conductance_y)) .or. &
            .not. all(ieee_is_finite(conductance_z))) &
            error stop 'Elliptic solver: non-finite conductance data'
    end subroutine validate_problem


    !> Check allocation and exact shapes of all boundary arrays.
    subroutine validate_boundary_shapes(boundary, nx, ny, nz)
        type(elliptic_boundary_data), intent(in) :: boundary
        integer, intent(in) :: nx, ny, nz

        if (.not. allocated(boundary%type_x) .or. .not. allocated(boundary%type_y) .or. &
            .not. allocated(boundary%type_z) .or. .not. allocated(boundary%value_x) .or. &
            .not. allocated(boundary%value_y) .or. .not. allocated(boundary%value_z)) &
            error stop 'Elliptic solver: boundary arrays are not allocated'

        if (any(shape(boundary%type_x) /= [nx+1,ny,nz]) .or. &
            any(shape(boundary%value_x) /= [nx+1,ny,nz])) &
            error stop 'Elliptic solver: invalid x-boundary array shape'
        if (any(shape(boundary%type_y) /= [nx,ny+1,nz]) .or. &
            any(shape(boundary%value_y) /= [nx,ny+1,nz])) &
            error stop 'Elliptic solver: invalid y-boundary array shape'
        if (any(shape(boundary%type_z) /= [nx,ny,nz+1]) .or. &
            any(shape(boundary%value_z) /= [nx,ny,nz+1])) &
            error stop 'Elliptic solver: invalid z-boundary array shape'

        if (.not. all(ieee_is_finite(boundary%value_x)) .or. &
            .not. all(ieee_is_finite(boundary%value_y)) .or. &
            .not. all(ieee_is_finite(boundary%value_z))) &
            error stop 'Elliptic solver: non-finite boundary data'
    end subroutine validate_boundary_shapes


    !> Check that every used boundary marker has a supported integer value.
    subroutine validate_boundary_types(boundary, dimensions)
        type(elliptic_boundary_data), intent(in) :: boundary
        integer, intent(in) :: dimensions

        if (any(boundary%type_x < elliptic_bc_internal) .or. &
            any(boundary%type_x > elliptic_bc_dirichlet)) &
            error stop 'Elliptic solver: invalid x-boundary type'
        if (dimensions >= 2) then
            if (any(boundary%type_y < elliptic_bc_internal) .or. &
                any(boundary%type_y > elliptic_bc_dirichlet)) &
                error stop 'Elliptic solver: invalid y-boundary type'
        end if
        if (dimensions >= 3) then
            if (any(boundary%type_z < elliptic_bc_internal) .or. &
                any(boundary%type_z > elliptic_bc_dirichlet)) &
                error stop 'Elliptic solver: invalid z-boundary type'
        end if
    end subroutine validate_boundary_types

end module elliptic_multigrid_solver_class
