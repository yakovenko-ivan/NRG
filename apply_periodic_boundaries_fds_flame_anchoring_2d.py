#!/usr/bin/env python3
"""Apply NRG periodic-boundary support to the shared elliptic solver and FDS low-Mach solver.

Prerequisite: run apply_periodic_boundaries_stage1.py first.
Target base: feature/flame-anchoring-2d @ 82f470052860e7d0985ce9640a23f4c81fad0fd4.

Scope:
  * shared elliptic operator accepts periodic(3) topology without introducing a
    fake physical boundary type;
  * periodic neighbours wrap on every multigrid level;
  * 1-D periodic all-active systems use a direct cyclic/null-space-safe solve;
  * red/black smoothing is forced serial on odd periodic levels, avoiding a
    same-colour seam race;
  * FDS passes domain periodicity to its prepared pressure operator;
  * FDS synchronizes periodic cell and staggered face states;
  * FDS CHARM/EOS-corrected species reconstruction wraps its full four-point
    stencil instead of dropping to the ordinary physical-boundary fallback.

The script is defensive: source edits must match the expected snapshot exactly.
"""
from __future__ import annotations

from pathlib import Path
import re

ROOT = Path.cwd()


def read(rel: str) -> str:
    p = ROOT / rel
    if not p.exists():
        raise SystemExit(f"missing expected file: {p}")
    return p.read_text(encoding="utf-8")


def write(rel: str, text: str) -> None:
    (ROOT / rel).write_text(text, encoding="utf-8")
    print(f"updated {rel}")


def replace_once(text: str, old: str, new: str, label: str) -> str:
    n = text.count(old)
    if n != 1:
        raise SystemExit(f"{label}: expected exactly one match, found {n}")
    return text.replace(old, new, 1)


def replace_exact_count(text: str, old: str, new: str, expected: int, label: str) -> str:
    n = text.count(old)
    if n != expected:
        raise SystemExit(f"{label}: expected {expected} matches, found {n}")
    return text.replace(old, new)


def replace_regex_once(text: str, pattern: str, repl: str, label: str, flags=0) -> str:
    out, n = re.subn(pattern, repl, text, count=1, flags=flags)
    if n != 1:
        raise SystemExit(f"{label}: expected exactly one regex match, found {n}")
    return out


def replace_procedure(text: str, name: str, new_body: str, kind: str = "subroutine") -> str:
    if kind == "subroutine":
        pattern = rf"(?ims)^\s*(?:recursive\s+)?subroutine\s+{re.escape(name)}\b.*?^\s*end\s+subroutine(?:\s+{re.escape(name)})?\s*$"
    elif kind == "real_function":
        pattern = rf"(?ims)^\s*real\(dp\)\s+function\s+{re.escape(name)}\b.*?^\s*end\s+function(?:\s+{re.escape(name)})?\s*$"
    elif kind == "logical_function":
        pattern = rf"(?ims)^\s*logical\s+function\s+{re.escape(name)}\b.*?^\s*end\s+function(?:\s+{re.escape(name)})?\s*$"
    else:
        raise ValueError(kind)
    out, n = re.subn(pattern, "\n" + new_body.strip("\n"), text, count=1)
    if n != 1:
        raise SystemExit(f"could not replace {kind} {name}; matches={n}")
    return out


def insert_before_end_module(text: str, module_name: str, insertion: str) -> str:
    marker = f"end module {module_name}"
    return replace_once(text, marker, insertion.rstrip() + "\n\n" + marker, f"insert before {module_name} end")


# =============================================================================
# Shared elliptic solver
# =============================================================================
rel = "computing_module/src/current_build/elliptic_multigrid_solver.f90"
s = read(rel)

s = replace_once(
    s,
    "        integer :: dimensions = 1\n        real(dp), allocatable :: x(:,:,:)\n",
    "        integer :: dimensions = 1\n"
    "        logical, dimension(3) :: periodic = .false.\n"
    "        real(dp), allocatable :: x(:,:,:)\n",
    "elliptic: level periodic member",
)

s = replace_once(
    s,
    "        logical :: prepared_pure_neumann = .false.\n        integer :: prepared_dimensions = 0\n",
    "        logical :: prepared_pure_neumann = .false.\n"
    "        logical, dimension(3) :: prepared_periodic = .false.\n"
    "        integer :: prepared_dimensions = 0\n",
    "elliptic: prepared periodic member",
)

solve_driver = r'''
    subroutine solve_elliptic_problem(this, solution, rhs, conductance_x, conductance_y, conductance_z, &
        cell_volume, active, boundary, dimensions, use_initial_guess, statistics, periodic)

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
        logical, dimension(3), intent(in), optional :: periodic

        type(elliptic_solver_statistics) :: stat
        logical :: warm_start, pure_neumann
        logical, dimension(3) :: periodic_axes

        periodic_axes = .false.
        if (present(periodic)) periodic_axes = periodic
        if (dimensions < 3) periodic_axes(dimensions+1:3) = .false.

        call validate_problem(solution, rhs, conductance_x, conductance_y, conductance_z, &
            cell_volume, active, boundary, dimensions)
        call validate_periodic_boundary_data(boundary, periodic_axes, dimensions)
        call validate_solver_controls(this)

        warm_start = .false.
        if (present(use_initial_guess)) warm_start = use_initial_guess
        pure_neumann = .not. has_active_dirichlet_boundary(active, boundary, dimensions)

        if (dimensions == 1 .and. all(active)) then
            call this%invalidate()
            if (periodic_axes(1)) then
                call solve_periodic_tridiagonal_1d(solution, rhs, conductance_x, cell_volume, &
                    this%relative_tolerance, this%absolute_tolerance, stat)
            else
                call solve_tridiagonal_1d(solution, rhs, conductance_x, cell_volume, boundary, &
                    pure_neumann, this%relative_tolerance, this%absolute_tolerance, stat)
            end if
            if (present(statistics)) statistics = stat
            return
        end if

        call this%prepare(conductance_x, conductance_y, conductance_z, cell_volume, active, &
            boundary, dimensions, periodic=periodic_axes)
        call this%solve_prepared(solution, rhs, warm_start, stat)
        if (present(statistics)) statistics = stat
    end subroutine solve_elliptic_problem
'''
s = replace_procedure(s, "solve_elliptic_problem", solve_driver)

prepare = r'''
    subroutine prepare_elliptic_operator(this, conductance_x, conductance_y, conductance_z, &
        cell_volume, active, boundary, dimensions, periodic)

        class(elliptic_multigrid_solver), intent(inout) :: this
        real(dp), intent(in) :: conductance_x(:,:,:), conductance_y(:,:,:), conductance_z(:,:,:)
        real(dp), intent(in) :: cell_volume(:,:,:)
        logical, intent(in) :: active(:,:,:)
        type(elliptic_boundary_data), intent(in) :: boundary
        integer, intent(in) :: dimensions
        logical, dimension(3), intent(in), optional :: periodic

        integer :: nx, ny, nz, level_index
        logical, dimension(3) :: periodic_axes

        nx = size(active,1)
        ny = size(active,2)
        nz = size(active,3)
        periodic_axes = .false.
        if (present(periodic)) periodic_axes = periodic
        if (dimensions < 3) periodic_axes(dimensions+1:3) = .false.

        call validate_operator_data(conductance_x, conductance_y, conductance_z, cell_volume, &
            active, boundary, dimensions)
        call validate_periodic_boundary_data(boundary, periodic_axes, dimensions)
        call validate_solver_controls(this)

        call this%build_hierarchy(nx, ny, nz, dimensions)
        do level_index = 1, size(this%level)
            this%level(level_index)%periodic = periodic_axes
        end do

        this%level(1)%rhs = 0.0_dp
        this%level(1)%volume = cell_volume
        this%level(1)%active = active
        this%level(1)%conductance_x = conductance_x
        this%level(1)%conductance_y = conductance_y
        this%level(1)%conductance_z = conductance_z
        call copy_boundary_data(this%level(1)%boundary, boundary)
        call this%coarsen_hierarchy_data()
        call assemble_hierarchy_operator(this)

        this%prepared_pure_neumann = &
            .not. has_active_dirichlet_boundary(active, boundary, dimensions)
        this%prepared_periodic = periodic_axes
        this%prepared_dimensions = dimensions
        this%operator_prepared = .true.
    end subroutine prepare_elliptic_operator
'''
s = replace_procedure(s, "prepare_elliptic_operator", prepare)

prepared_solve = r'''
    subroutine solve_prepared_elliptic_problem(this, solution, rhs, use_initial_guess, statistics)
        class(elliptic_multigrid_solver), intent(inout) :: this
        real(dp), intent(inout) :: solution(:,:,:)
        real(dp), intent(in) :: rhs(:,:,:)
        logical, intent(in), optional :: use_initial_guess
        type(elliptic_solver_statistics), intent(out), optional :: statistics

        type(elliptic_solver_statistics) :: stat
        integer :: cycle
        logical :: warm_start
        real(dp) :: forcing_norm, tolerance

        if (.not. this%operator_prepared) &
            error stop 'Elliptic solver: solve_prepared called without prepare'
        if (any(shape(solution) /= [this%level(1)%nx,this%level(1)%ny,this%level(1)%nz]) .or. &
            any(shape(rhs) /= [this%level(1)%nx,this%level(1)%ny,this%level(1)%nz])) &
            error stop 'Elliptic solver: prepared solve has incompatible cell-array shape'
        if (.not. all(ieee_is_finite(solution)) .or. .not. all(ieee_is_finite(rhs))) &
            error stop 'Elliptic solver: non-finite prepared solve data'
        call validate_solver_controls(this)

        stat = elliptic_solver_statistics()
        stat%hierarchy_levels = size(this%level)
        warm_start = .false.
        if (present(use_initial_guess)) warm_start = use_initial_guess

        if (this%prepared_dimensions == 1 .and. all(this%level(1)%active)) then
            if (this%prepared_periodic(1)) then
                call solve_periodic_tridiagonal_1d(solution, rhs, this%level(1)%conductance_x, &
                    this%level(1)%volume, this%relative_tolerance, this%absolute_tolerance, stat)
            else
                call solve_tridiagonal_1d(solution, rhs, this%level(1)%conductance_x, &
                    this%level(1)%volume, this%level(1)%boundary, this%prepared_pure_neumann, &
                    this%relative_tolerance, this%absolute_tolerance, stat)
            end if
            if (present(statistics)) statistics = stat
            return
        end if

        this%level(1)%rhs = rhs
        if (warm_start) then
            this%level(1)%x = solution
        else
            this%level(1)%x = 0.0_dp
        end if
        where (.not. this%level(1)%active) this%level(1)%x = 0.0_dp

        if (this%prepared_pure_neumann) then
            call enforce_neumann_compatibility(this%level(1), stat%compatibility_correction)
            call remove_volume_weighted_mean(this%level(1))
        end if

        call this%calculate_residual(1)
        stat%residual_evaluations = stat%residual_evaluations + 1
        stat%initial_residual = active_max_abs(this%level(1)%residual, this%level(1)%active)
        forcing_norm = calculate_forcing_norm(this%level(1))
        tolerance = this%absolute_tolerance + this%relative_tolerance*forcing_norm

        if (stat%initial_residual <= tolerance) then
            stat%converged = .true.
        else
            do cycle = 1, this%maximum_cycles
                call this%v_cycle(1, stat%smoothing_iterations, &
                    stat%parallel_smoothing_iterations, stat%serial_smoothing_iterations, &
                    stat%residual_evaluations)
                if (this%prepared_pure_neumann) call remove_volume_weighted_mean(this%level(1))

                call this%calculate_residual(1)
                stat%residual_evaluations = stat%residual_evaluations + 1
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
        where (.not. this%level(1)%active) solution = 0.0_dp
        if (present(statistics)) statistics = stat
    end subroutine solve_prepared_elliptic_problem
'''
s = replace_procedure(s, "solve_prepared_elliptic_problem", prepared_solve)

invalidate = r'''
    subroutine invalidate_prepared_operator(this)
        class(elliptic_multigrid_solver), intent(inout) :: this
        this%operator_prepared = .false.
        this%prepared_pure_neumann = .false.
        this%prepared_periodic = .false.
        this%prepared_dimensions = 0
    end subroutine invalidate_prepared_operator
'''
s = replace_procedure(s, "invalidate_prepared_operator", invalidate)

# Periodic direct 1-D path: pin one cell, solve the remaining tridiagonal block,
# then restore the volume-weighted zero-mean gauge.
periodic_direct = r'''

    !> Direct solve for an all-active one-dimensional periodic operator.
    !! The periodic Laplacian has the same constant null space as a pure-Neumann
    !! problem.  Compatibility is enforced, cell 1 is temporarily pinned to
    !! zero, the remaining (n-1)x(n-1) tridiagonal block is solved, and the
    !! volume-weighted mean is removed afterwards.
    subroutine solve_periodic_tridiagonal_1d(solution, rhs, conductance_x, cell_volume, &
        relative_tolerance, absolute_tolerance, statistics)

        real(dp), intent(inout) :: solution(:,:,:)
        real(dp), intent(in) :: rhs(:,:,:)
        real(dp), intent(in) :: conductance_x(:,:,:)
        real(dp), intent(in) :: cell_volume(:,:,:)
        real(dp), intent(in) :: relative_tolerance, absolute_tolerance
        type(elliptic_solver_statistics), intent(out) :: statistics

        real(dp), allocatable :: lower(:), diagonal(:), upper(:), source(:)
        real(dp), allocatable :: compatible_rhs(:), x_line(:)
        real(dp) :: factor, forcing_norm, tolerance, seam_conductance
        real(dp) :: left_g, right_g, residual
        integer :: n, m, i, left_i, right_i

        statistics = elliptic_solver_statistics()
        n = size(rhs,1)
        if (n < 1) error stop 'Elliptic periodic 1-D solver: empty system'

        allocate(compatible_rhs(n), x_line(n))
        compatible_rhs = rhs(:,1,1)
        statistics%compatibility_correction = &
            sum(compatible_rhs*cell_volume(:,1,1))/sum(cell_volume(:,1,1))
        compatible_rhs = compatible_rhs - statistics%compatibility_correction

        if (n == 1) then
            solution(1,1,1) = 0.0_dp
            statistics%cycles = 1
            statistics%initial_residual = abs(compatible_rhs(1))
            statistics%final_residual = 0.0_dp
            statistics%relative_residual = 0.0_dp
            statistics%converged = .true.
            return
        end if

        seam_conductance = 0.5_dp*(conductance_x(1,1,1) + conductance_x(n+1,1,1))
        if (seam_conductance <= tiny(1.0_dp)) &
            error stop 'Elliptic periodic 1-D solver: non-positive seam conductance'

        allocate(lower(n-1), diagonal(n-1), upper(n-1), source(n-1))
        lower = 0.0_dp
        diagonal = 0.0_dp
        upper = 0.0_dp
        source = 0.0_dp

        ! Unknowns are original cells 2:n; original cell 1 is pinned to zero.
        do i = 2, n
            m = i - 1
            left_g = conductance_x(i,1,1)
            if (i == n) then
                right_g = seam_conductance
            else
                right_g = conductance_x(i+1,1,1)
            end if
            diagonal(m) = left_g + right_g
            if (i > 2) lower(m) = -left_g
            if (i < n) upper(m) = -right_g
            source(m) = compatible_rhs(i)*cell_volume(i,1,1)
        end do

        do m = 2, n-1
            if (abs(diagonal(m-1)) <= tiny(1.0_dp)) &
                error stop 'Elliptic periodic 1-D solver: zero pivot'
            factor = lower(m)/diagonal(m-1)
            diagonal(m) = diagonal(m) - factor*upper(m-1)
            source(m) = source(m) - factor*source(m-1)
        end do
        if (abs(diagonal(n-1)) <= tiny(1.0_dp)) &
            error stop 'Elliptic periodic 1-D solver: zero final pivot'

        x_line = 0.0_dp
        x_line(n) = source(n-1)/diagonal(n-1)
        do m = n-2, 1, -1
            i = m + 1
            x_line(i) = (source(m)-upper(m)*x_line(i+1))/diagonal(m)
        end do

        x_line = x_line - sum(x_line*cell_volume(:,1,1))/sum(cell_volume(:,1,1))
        solution(:,1,1) = x_line

        forcing_norm = max(maxval(abs(compatible_rhs)), tiny(1.0_dp))
        statistics%initial_residual = forcing_norm
        tolerance = absolute_tolerance + relative_tolerance*forcing_norm
        statistics%final_residual = 0.0_dp
        do i = 1, n
            left_i = i-1
            if (left_i < 1) left_i = n
            right_i = i+1
            if (right_i > n) right_i = 1
            if (i == 1) then
                left_g = seam_conductance
            else
                left_g = conductance_x(i,1,1)
            end if
            if (i == n) then
                right_g = seam_conductance
            else
                right_g = conductance_x(i+1,1,1)
            end if
            residual = compatible_rhs(i) + &
                (left_g*x_line(left_i) + right_g*x_line(right_i) - &
                 (left_g+right_g)*x_line(i))/cell_volume(i,1,1)
            statistics%final_residual = max(statistics%final_residual, abs(residual))
        end do
        statistics%cycles = 1
        statistics%relative_residual = statistics%final_residual/forcing_norm
        statistics%converged = statistics%final_residual <= tolerance
    end subroutine solve_periodic_tridiagonal_1d
'''
s = replace_once(
    s,
    "\n\n    !> Add one physical boundary contribution to a 1-D matrix row.\n",
    periodic_direct + "\n\n    !> Add one physical boundary contribution to a 1-D matrix row.\n",
    "elliptic: periodic direct solver",
)

accumulate = r'''
    subroutine accumulate_operator_face(level, neighbour_i, neighbour_j, neighbour_k, conductance, &
        boundary_type, boundary_value, diagonal_value, source_value, regular_cell)

        type(multigrid_level), intent(in) :: level
        integer, intent(in) :: neighbour_i, neighbour_j, neighbour_k, boundary_type
        real(dp), intent(in) :: conductance, boundary_value
        real(dp), intent(inout) :: diagonal_value, source_value
        logical, intent(inout) :: regular_cell
        logical :: neighbour_is_active, wrapped
        integer :: mapped_i, mapped_j, mapped_k

        mapped_i = neighbour_i
        mapped_j = neighbour_j
        mapped_k = neighbour_k
        wrapped = .false.

        if ((mapped_i < 1 .or. mapped_i > level%nx) .and. level%periodic(1)) then
            mapped_i = 1 + modulo(mapped_i-1, level%nx)
            wrapped = .true.
        end if
        if ((mapped_j < 1 .or. mapped_j > level%ny) .and. level%periodic(2)) then
            mapped_j = 1 + modulo(mapped_j-1, level%ny)
            wrapped = .true.
        end if
        if ((mapped_k < 1 .or. mapped_k > level%nz) .and. level%periodic(3)) then
            mapped_k = 1 + modulo(mapped_k-1, level%nz)
            wrapped = .true.
        end if

        neighbour_is_active = mapped_i >= 1 .and. mapped_i <= level%nx .and. &
            mapped_j >= 1 .and. mapped_j <= level%ny .and. &
            mapped_k >= 1 .and. mapped_k <= level%nz
        if (neighbour_is_active) neighbour_is_active = level%active(mapped_i,mapped_j,mapped_k)

        if (neighbour_is_active) then
            diagonal_value = diagonal_value + conductance
            if (wrapped) regular_cell = .false.
        else
            regular_cell = .false.
            call add_boundary_terms(conductance, boundary_type, boundary_value, &
                diagonal_value, source_value)
        end if
    end subroutine accumulate_operator_face
'''
s = replace_procedure(s, "accumulate_operator_face", accumulate)

neighbour_sum = r'''
    real(dp) function preassembled_neighbour_sum(level, i, j, k) result(neighbour_sum)
        type(multigrid_level), intent(in) :: level
        integer, intent(in) :: i, j, k

        neighbour_sum = 0.0_dp

        if (level%regular_interior(i,j,k)) then
            neighbour_sum = level%conductance_x(i,j,k)*level%x(i-1,j,k) + &
                level%conductance_x(i+1,j,k)*level%x(i+1,j,k)
            if (level%dimensions >= 2) then
                neighbour_sum = neighbour_sum + &
                    level%conductance_y(i,j,k)*level%x(i,j-1,k) + &
                    level%conductance_y(i,j+1,k)*level%x(i,j+1,k)
            end if
            if (level%dimensions >= 3) then
                neighbour_sum = neighbour_sum + &
                    level%conductance_z(i,j,k)*level%x(i,j,k-1) + &
                    level%conductance_z(i,j,k+1)*level%x(i,j,k+1)
            end if
            return
        end if

        call add_periodic_or_internal_neighbour(level, i,j,k, 1,-1, level%conductance_x(i,j,k), neighbour_sum)
        call add_periodic_or_internal_neighbour(level, i,j,k, 1, 1, level%conductance_x(i+1,j,k), neighbour_sum)
        if (level%dimensions >= 2) then
            call add_periodic_or_internal_neighbour(level, i,j,k, 2,-1, level%conductance_y(i,j,k), neighbour_sum)
            call add_periodic_or_internal_neighbour(level, i,j,k, 2, 1, level%conductance_y(i,j+1,k), neighbour_sum)
        end if
        if (level%dimensions >= 3) then
            call add_periodic_or_internal_neighbour(level, i,j,k, 3,-1, level%conductance_z(i,j,k), neighbour_sum)
            call add_periodic_or_internal_neighbour(level, i,j,k, 3, 1, level%conductance_z(i,j,k+1), neighbour_sum)
        end if
    end function preassembled_neighbour_sum
'''
s = replace_procedure(s, "preassembled_neighbour_sum", neighbour_sum, kind="real_function")

add_neighbor_helper = r'''

    subroutine add_periodic_or_internal_neighbour(level, i, j, k, dimension, direction, conductance, neighbour_sum)
        type(multigrid_level), intent(in) :: level
        integer, intent(in) :: i, j, k, dimension, direction
        real(dp), intent(in) :: conductance
        real(dp), intent(inout) :: neighbour_sum
        integer :: ni, nj, nk

        ni = i
        nj = j
        nk = k
        select case (dimension)
        case (1)
            ni = i + direction
            if ((ni < 1 .or. ni > level%nx) .and. level%periodic(1)) ni = 1 + modulo(ni-1, level%nx)
        case (2)
            nj = j + direction
            if ((nj < 1 .or. nj > level%ny) .and. level%periodic(2)) nj = 1 + modulo(nj-1, level%ny)
        case (3)
            nk = k + direction
            if ((nk < 1 .or. nk > level%nz) .and. level%periodic(3)) nk = 1 + modulo(nk-1, level%nz)
        end select

        if (ni < 1 .or. ni > level%nx .or. nj < 1 .or. nj > level%ny .or. &
            nk < 1 .or. nk > level%nz) return
        if (level%active(ni,nj,nk)) neighbour_sum = neighbour_sum + conductance*level%x(ni,nj,nk)
    end subroutine add_periodic_or_internal_neighbour
'''
s = replace_once(
    s,
    "\n\n    !> Add one Dirichlet or Neumann contribution to a cell equation.\n",
    add_neighbor_helper + "\n\n    !> Add one Dirichlet or Neumann contribution to a cell equation.\n",
    "elliptic: neighbour helper",
)

coarse_neighbour = r'''
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
            if ((neighbour_i < 1 .or. neighbour_i > level%nx) .and. level%periodic(1)) &
                neighbour_i = 1 + modulo(neighbour_i-1, level%nx)
        case (2)
            neighbour_j = j + direction
            if ((neighbour_j < 1 .or. neighbour_j > level%ny) .and. level%periodic(2)) &
                neighbour_j = 1 + modulo(neighbour_j-1, level%ny)
        case (3)
            neighbour_k = k + direction
            if ((neighbour_k < 1 .or. neighbour_k > level%nz) .and. level%periodic(3)) &
                neighbour_k = 1 + modulo(neighbour_k-1, level%nz)
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
            value = -center_value
        else
            value = center_value
        end if
    end function coarse_neighbour_value
'''
s = replace_procedure(s, "coarse_neighbour_value", coarse_neighbour, kind="real_function")

# Red-black periodic seam is parallel-safe only when every periodic extent is even.
s = replace_once(
    s,
    "        use_parallel = this%level(level_index)%active_cells >= this%minimum_parallel_cells\n\n        do iteration = 1, iterations\n",
    "        use_parallel = this%level(level_index)%active_cells >= this%minimum_parallel_cells\n"
    "        use_parallel = use_parallel .and. periodic_red_black_safe(this%level(level_index))\n\n"
    "        do iteration = 1, iterations\n",
    "elliptic: safe periodic smoother",
)

periodic_safety = r'''

    logical function periodic_red_black_safe(level) result(safe)
        type(multigrid_level), intent(in) :: level

        safe = .true.
        if (level%periodic(1) .and. mod(level%nx,2) /= 0) safe = .false.
        if (level%dimensions >= 2 .and. level%periodic(2) .and. mod(level%ny,2) /= 0) safe = .false.
        if (level%dimensions >= 3 .and. level%periodic(3) .and. mod(level%nz,2) /= 0) safe = .false.
    end function periodic_red_black_safe
'''
s = replace_once(
    s,
    "\n\n    !> Compute the cell-volume-normalized residual on one level.\n",
    periodic_safety + "\n\n    !> Compute the cell-volume-normalized residual on one level.\n",
    "elliptic: periodic coloring helper",
)

# Validate that a periodic axis is connectivity, never a physical boundary face.
validation = r'''

    subroutine validate_periodic_boundary_data(boundary, periodic, dimensions)
        type(elliptic_boundary_data), intent(in) :: boundary
        logical, dimension(3), intent(in) :: periodic
        integer, intent(in) :: dimensions

        if (periodic(1)) then
            if (any(boundary%type_x(1,:,:) /= elliptic_bc_internal) .or. &
                any(boundary%type_x(size(boundary%type_x,1),:,:) /= elliptic_bc_internal)) &
                error stop 'Elliptic solver: periodic x faces must be marked internal'
        end if
        if (dimensions >= 2 .and. periodic(2)) then
            if (any(boundary%type_y(:,1,:) /= elliptic_bc_internal) .or. &
                any(boundary%type_y(:,size(boundary%type_y,2),:) /= elliptic_bc_internal)) &
                error stop 'Elliptic solver: periodic y faces must be marked internal'
        end if
        if (dimensions >= 3 .and. periodic(3)) then
            if (any(boundary%type_z(:,:,1) /= elliptic_bc_internal) .or. &
                any(boundary%type_z(:,:,size(boundary%type_z,3)) /= elliptic_bc_internal)) &
                error stop 'Elliptic solver: periodic z faces must be marked internal'
        end if
    end subroutine validate_periodic_boundary_data
'''
# Place near the end; no dependency on later routines.
s = insert_before_end_module(s, "elliptic_multigrid_solver_class", validation)
write(rel, s)


# =============================================================================
# FDS low-Mach solver
# =============================================================================
rel = "computing_module/src/current_build/fds_low_mach_solver.f90"
s = read(rel)

s = replace_once(
    s,
    "    use computational_mesh_class\n    use boundary_conditions_class\n",
    "    use computational_mesh_class\n"
    "    use mpi_communications_class\n"
    "    use boundary_conditions_class\n",
    "fds: use mpi communications",
)

s = replace_once(
    s,
    "        type(computational_domain) :: domain\n        type(thermophysical_properties_pointer) :: thermo\n",
    "        type(computational_domain) :: domain\n"
    "        type(mpi_communications) :: mpi_support\n"
    "        type(thermophysical_properties_pointer) :: thermo\n",
    "fds: mpi support member",
)

s = replace_once(
    s,
    "        procedure, private :: calculate_interm_Y_corrector\n        procedure, private :: apply_boundary_conditions\n",
    "        procedure, private :: calculate_interm_Y_corrector\n"
    "        procedure, private :: apply_boundary_conditions\n"
    "        procedure, private :: synchronize_periodic_cell_state\n"
    "        procedure, private :: synchronize_periodic_face_velocity\n"
    "        procedure, private :: synchronize_periodic_pressure_fields\n",
    "fds: periodic helper bindings",
)

s = replace_once(
    s,
    "        constructor%domain                = manager%domain\n        constructor%thermo%thermo_ptr    => manager%thermophysics%thermo_ptr\n",
    "        constructor%domain                = manager%domain\n"
    "        constructor%mpi_support           = manager%mpi_communications\n"
    "        constructor%thermo%thermo_ptr    => manager%thermophysics%thermo_ptr\n",
    "fds: copy mpi support",
)

# After restart/initial data input but before staggered velocity reconstruction,
# wrap the primary cell state so the initial periodic face values see the true neighbour.
s = replace_once(
    s,
    "        constructor%load_counter    = problem_data_io%get_load_counter()\n\n        dimensions        = manager%domain%get_domain_dimensions()\n",
    "        constructor%load_counter    = problem_data_io%get_load_counter()\n\n"
    "        call constructor%synchronize_periodic_cell_state(.false.)\n\n"
    "        dimensions        = manager%domain%get_domain_dimensions()\n",
    "fds: initial periodic primary sync",
)

# Re-sync after EOS initialization, then make the two representations of each
# periodic staggered interface identical.
s = replace_once(
    s,
    "            call constructor%state_eq%apply_boundary_conditions_for_initial_conditions()\n        end if\n\n        constructor%time                = calculation_time\n",
    "            call constructor%state_eq%apply_boundary_conditions_for_initial_conditions()\n"
    "        end if\n\n"
    "        call constructor%synchronize_periodic_cell_state(.true.)\n"
    "        call constructor%synchronize_periodic_face_velocity()\n\n"
    "        constructor%time                = calculation_time\n",
    "fds: post-EOS periodic initialization",
)

# Time-step synchronization points.
s = replace_once(
    s,
    "        call fds_timer%tic()\n\n        this%time = this%time + this%time_step\n",
    "        call fds_timer%tic()\n\n"
    "        call this%synchronize_periodic_cell_state(.true.)\n"
    "        call this%synchronize_periodic_face_velocity()\n\n"
    "        this%time = this%time + this%time_step\n",
    "fds: step-start periodic sync",
)

s = replace_once(
    s,
    "        call this%apply_boundary_conditions(this%time_step,predictor=.true.)\n        call this%calculate_divergence_v        (this%time_step,predictor=.true.)\n",
    "        call this%apply_boundary_conditions(this%time_step,predictor=.true.)\n"
    "        call this%synchronize_periodic_cell_state(.true.)\n"
    "        call this%synchronize_periodic_face_velocity()\n"
    "        call this%calculate_divergence_v        (this%time_step,predictor=.true.)\n",
    "fds: predictor pre-projection sync",
)

s = replace_once(
    s,
    "        call this%calculate_velocity            (this%time_step,predictor=.true.)\n        call fds_gas_dynamics_timer%toc(new_iter=.true.)\n",
    "        call this%calculate_velocity            (this%time_step,predictor=.true.)\n"
    "        call this%synchronize_periodic_face_velocity()\n"
    "        call fds_gas_dynamics_timer%toc(new_iter=.true.)\n",
    "fds: predictor post-velocity sync",
)

s = replace_once(
    s,
    "        call this%apply_boundary_conditions(this%time_step,predictor=.false.)\n        call this%calculate_divergence_v        (this%time_step,predictor=.false.)\n",
    "        call this%apply_boundary_conditions(this%time_step,predictor=.false.)\n"
    "        call this%synchronize_periodic_cell_state(.true.)\n"
    "        call this%synchronize_periodic_face_velocity()\n"
    "        call this%calculate_divergence_v        (this%time_step,predictor=.false.)\n",
    "fds: corrector pre-projection sync",
)

s = replace_once(
    s,
    "        call this%calculate_velocity            (this%time_step,predictor=.false.)\n\n        if (this%CFL_condition_flag) then\n",
    "        call this%calculate_velocity            (this%time_step,predictor=.false.)\n"
    "        call this%synchronize_periodic_face_velocity()\n\n"
    "        if (this%CFL_condition_flag) then\n",
    "fds: corrector post-velocity sync",
)

# Predictor/corrector know which directions wrap.  The current
# feature/flame-anchoring-2d branch has slightly different formatting in the
# corrector and an additional boundary-style condition elsewhere in FDS, so
# patch these two procedures independently instead of relying on global counts.
def patch_species_transport_procedure(text: str, name: str) -> str:
    pattern = rf"(?ims)^\s*subroutine\s+{re.escape(name)}\b.*?^\s*end\s+subroutine(?:\s+{re.escape(name)})?\s*$"
    match = re.search(pattern, text)
    if not match:
        raise SystemExit(f"fds: could not locate {name}")
    proc = match.group(0)

    proc = regex_once(
        proc,
        r"(?m)^(\s*integer\s*,dimension\(3,2\)\s*::\s*cons_inner_loop\s*)$",
        r"\1\n        logical    ,dimension(3)      :: periodic",
        f"fds: {name} periodic declaration",
    )
    proc = regex_once(
        proc,
        r"(?m)^(\s*cons_inner_loop\s*=\s*this%domain%get_local_inner_cells_bounds\(\)\s*)$",
        r"\1\n        periodic           = this%domain%get_periodic_directions()",
        f"fds: {name} periodic assignment",
    )

    proc = replace_once(
        proc,
        "if ((i*I_m(dim,1) + j*I_m(dim,2)  + k*I_m(dim,3)) < cons_inner_loop(dim,2)) then",
        "if (((i*I_m(dim,1) + j*I_m(dim,2) + k*I_m(dim,3)) < cons_inner_loop(dim,2)) .or. periodic(dim)) then",
        f"fds: {name} periodic right-face CHARM path",
    )
    proc = replace_once(
        proc,
        "if ((i*I_m(dim,1) + j*I_m(dim,2)  + k*I_m(dim,3)) > cons_inner_loop(dim,1)) then",
        "if (((i*I_m(dim,1) + j*I_m(dim,2) + k*I_m(dim,3)) > cons_inner_loop(dim,1)) .or. periodic(dim)) then",
        f"fds: {name} periodic left-face CHARM path",
    )
    proc = replace_once(
        proc,
        "v_f%pr(dim)%cells(dim,i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3)),flux_right_vec)",
        "v_f%pr(dim)%cells(dim,i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3)),periodic,cons_inner_loop,flux_right_vec)",
        f"fds: {name} right CHARM call metadata",
    )
    proc = replace_once(
        proc,
        "v_f%pr(dim)%cells(dim,i,j,k),flux_left_vec)",
        "v_f%pr(dim)%cells(dim,i,j,k),periodic,cons_inner_loop,flux_left_vec)",
        f"fds: {name} left CHARM call metadata",
    )

    return text[:match.start()] + proc + text[match.end():]


s = patch_species_transport_procedure(s, "calculate_interm_Y_predictor")
s = patch_species_transport_procedure(s, "calculate_interm_Y_corrector")

# Replace the helper with a periodic-aware four-point accessor.
pattern = r"(?ims)^\s*subroutine\s+eos_corrected_species_face_vector\b.*?^\s*end\s+subroutine(?:\s+eos_corrected_species_face_vector)?\s*$"
match = re.search(pattern, s)
if not match:
    raise SystemExit("fds: could not locate eos_corrected_species_face_vector")
old_func = match.group(0)
# Preserve the existing body algorithm while changing only its signature and its
# two stencil-index construction sites.  This avoids duplicating the validated
# EOS active-set correction logic in this patcher.
new_func = old_func
new_func = replace_once(
    new_func,
    "subroutine eos_corrected_species_face_vector(rho_field,Y_field,molar_masses,species_number,dim,i,j,k,face_side,velocity,phi)",
    "subroutine eos_corrected_species_face_vector(rho_field,Y_field,molar_masses,species_number,dim,i,j,k,face_side,velocity, &\n"
    "        periodic,inner_bounds,phi)",
    "fds: CHARM signature",
)
new_func = replace_once(
    new_func,
    "        real(dp), intent(in) :: velocity\n        real(dp), dimension(:), intent(out) :: phi\n",
    "        real(dp), intent(in) :: velocity\n"
    "        logical, dimension(3), intent(in) :: periodic\n"
    "        integer, dimension(3,2), intent(in) :: inner_bounds\n"
    "        real(dp), dimension(:), intent(out) :: phi\n",
    "fds: CHARM periodic args",
)
new_func = replace_once(
    new_func,
    "        ii = i + gamma_offset*I_m(dim,1)\n        jj = j + gamma_offset*I_m(dim,2)\n        kk = k + gamma_offset*I_m(dim,3)\n",
    "        call map_periodic_stencil_index(dim,i,j,k,gamma_offset,periodic,inner_bounds,ii,jj,kk)\n",
    "fds: CHARM gamma wrapped index",
)
new_func = replace_once(
    new_func,
    "                ii = i + offset*I_m(dim,1)\n                jj = j + offset*I_m(dim,2)\n                kk = k + offset*I_m(dim,3)\n",
    "                call map_periodic_stencil_index(dim,i,j,k,offset,periodic,inner_bounds,ii,jj,kk)\n",
    "fds: CHARM stencil wrapped index",
)
s = s[:match.start()] + new_func + s[match.end():]

# Pressure operator gets periodic connectivity.
s = replace_once(
    s,
    "        call this%pressure_solver%prepare(conductance_x, conductance_y, conductance_z, &\n            cell_volume, active, elliptic_boundary, dimensions)\n",
    "        call this%pressure_solver%prepare(conductance_x, conductance_y, conductance_z, &\n"
    "            cell_volume, active, elliptic_boundary, dimensions, &\n"
    "            periodic=this%domain%get_periodic_directions())\n",
    "fds: periodic pressure operator",
)

# Synchronize H immediately before dynamic-pressure recovery and p_dyn immediately
# afterwards, inside every nonlinear pressure iteration.
s = replace_regex_once(
    s,
    r"(?m)^(\s*)call this%calculate_dynamic_pressure\(time_step,predictor\)\s*$",
    r"\1call this%synchronize_periodic_pressure_fields()\n"
    r"\1call this%calculate_dynamic_pressure(time_step,predictor)\n"
    r"\1call this%synchronize_periodic_pressure_fields()",
    "fds: periodic pressure field sync",
)

# Periodic helper implementations.
helpers = r'''

    !--------------------------------------------------------------------------
    ! Periodic synchronization.  Stage 1 deliberately restricts periodic runs
    ! to one MPI rank, so the conservative exchange entry points reduce to local
    ! halo wrapping while preserving a future multi-rank API.
    !--------------------------------------------------------------------------
    subroutine synchronize_periodic_cell_state(this, include_intermediate)
        class(fds_solver), intent(inout) :: this
        logical, intent(in) :: include_intermediate
        logical, dimension(3) :: periodic

        periodic = this%domain%get_periodic_directions()
        if (.not. any(periodic)) return

        call this%mpi_support%exchange_conservative_scalar_field(this%rho%s_ptr)
        call this%mpi_support%exchange_conservative_scalar_field(this%T%s_ptr)
        call this%mpi_support%exchange_conservative_scalar_field(this%p%s_ptr)
        call this%mpi_support%exchange_conservative_vector_field(this%v%v_ptr)
        call this%mpi_support%exchange_conservative_vector_field(this%Y%v_ptr)

        if (include_intermediate) then
            call this%mpi_support%exchange_conservative_scalar_field(this%rho_int%s_ptr)
            call this%mpi_support%exchange_conservative_scalar_field(this%rho_old%s_ptr)
            call this%mpi_support%exchange_conservative_scalar_field(this%h_s%s_ptr)
            call this%mpi_support%exchange_conservative_scalar_field(this%mix_mol_mass%s_ptr)
            call this%mpi_support%exchange_conservative_vector_field(this%Y_int%v_ptr)
            call this%mpi_support%exchange_conservative_vector_field(this%Y_old%v_ptr)
        end if
    end subroutine synchronize_periodic_cell_state


    subroutine synchronize_periodic_pressure_fields(this)
        class(fds_solver), intent(inout) :: this
        logical, dimension(3) :: periodic

        periodic = this%domain%get_periodic_directions()
        if (.not. any(periodic)) return

        call this%mpi_support%exchange_conservative_scalar_field(this%H%s_ptr)
        call this%mpi_support%exchange_conservative_scalar_field(this%H_old%s_ptr)
        call this%mpi_support%exchange_conservative_scalar_field(this%p_dyn%s_ptr)
    end subroutine synchronize_periodic_pressure_fields


    !> Synchronize the staggered velocity representation across one-rank
    !! periodic seams.  Longitudinal storage has two representations of the
    !! same periodic interface; transverse storage follows cell-centred wrapping.
    subroutine synchronize_periodic_face_velocity(this)
        class(fds_solver), intent(inout) :: this
        logical, dimension(3) :: periodic
        integer, dimension(3,2) :: flow_inner
        integer :: dimensions, component, axis

        periodic = this%domain%get_periodic_directions()
        if (.not. any(periodic)) return

        dimensions = this%domain%get_domain_dimensions()
        flow_inner = this%domain%get_local_inner_faces_bounds()

        do component = 1, dimensions
            do axis = 1, dimensions
                if (.not. periodic(axis)) cycle
                call synchronize_face_component_axis(this%v_f%v_ptr%pr(component)%cells, &
                    component, axis, flow_inner)
                call synchronize_face_component_axis(this%v_f_old%v_ptr%pr(component)%cells, &
                    component, axis, flow_inner)
            end do
        end do
    end subroutine synchronize_periodic_face_velocity


    subroutine synchronize_face_component_axis(face_field, component, axis, flow_inner)
        real(dp), dimension(:,0:,0:,0:), intent(inout) :: face_field
        integer, intent(in) :: component, axis
        integer, dimension(3,2), intent(in) :: flow_inner
        integer :: lo, hi

        lo = flow_inner(axis,1)
        hi = flow_inner(axis,2)

        select case (axis)
        case (1)
            if (component == axis) then
                face_field(component,lo,:,:) = 0.5_dp*(face_field(component,lo,:,:) + face_field(component,hi,:,:))
                face_field(component,hi,:,:) = face_field(component,lo,:,:)
                face_field(component,lo-1,:,:) = face_field(component,hi-1,:,:)
                face_field(component,hi+1,:,:) = face_field(component,lo+1,:,:)
            else
                face_field(component,lo-1,:,:) = face_field(component,hi-1,:,:)
                face_field(component,hi,:,:) = face_field(component,lo,:,:)
            end if
        case (2)
            if (component == axis) then
                face_field(component,:,lo,:) = 0.5_dp*(face_field(component,:,lo,:) + face_field(component,:,hi,:))
                face_field(component,:,hi,:) = face_field(component,:,lo,:)
                face_field(component,:,lo-1,:) = face_field(component,:,hi-1,:)
                face_field(component,:,hi+1,:) = face_field(component,:,lo+1,:)
            else
                face_field(component,:,lo-1,:) = face_field(component,:,hi-1,:)
                face_field(component,:,hi,:) = face_field(component,:,lo,:)
            end if
        case (3)
            if (component == axis) then
                face_field(component,:,:,lo) = 0.5_dp*(face_field(component,:,:,lo) + face_field(component,:,:,hi))
                face_field(component,:,:,hi) = face_field(component,:,:,lo)
                face_field(component,:,:,lo-1) = face_field(component,:,:,hi-1)
                face_field(component,:,:,hi+1) = face_field(component,:,:,lo+1)
            else
                face_field(component,:,:,lo-1) = face_field(component,:,:,hi-1)
                face_field(component,:,:,hi) = face_field(component,:,:,lo)
            end if
        end select
    end subroutine synchronize_face_component_axis


    subroutine map_periodic_stencil_index(dim, i, j, k, offset, periodic, inner_bounds, ii, jj, kk)
        integer, intent(in) :: dim, i, j, k, offset
        logical, dimension(3), intent(in) :: periodic
        integer, dimension(3,2), intent(in) :: inner_bounds
        integer, intent(out) :: ii, jj, kk
        integer :: extent

        ii = i + offset*I_m(dim,1)
        jj = j + offset*I_m(dim,2)
        kk = k + offset*I_m(dim,3)
        if (.not. periodic(dim)) return

        extent = inner_bounds(dim,2) - inner_bounds(dim,1) + 1
        select case (dim)
        case (1)
            ii = inner_bounds(1,1) + modulo(ii-inner_bounds(1,1), extent)
        case (2)
            jj = inner_bounds(2,1) + modulo(jj-inner_bounds(2,1), extent)
        case (3)
            kk = inner_bounds(3,1) + modulo(kk-inner_bounds(3,1), extent)
        end select
    end subroutine map_periodic_stencil_index
'''
s = insert_before_end_module(s, "fds_low_mach_solver_class", helpers)
write(rel, s)

print("\nFDS + elliptic periodic-boundary implementation applied successfully.")
