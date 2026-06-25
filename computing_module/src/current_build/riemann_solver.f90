module riemann_solver_class

    use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
    use kind_parameters
    use global_data
    use data_manager_class
    use thermophysical_properties_class

    implicit none

#ifdef OMP
    include "omp_lib.h"
#endif

    private
    public  riemann_solver, riemann_solver_c

    real(dp), parameter :: riemann_temperature_floor = 1.0_dp
    real(dp), parameter :: riemann_temperature_ceiling = 20000.0_dp
    ! Shock Hugoniot solves must stay inside the thermodynamic model validity
    ! range.  Trial pressures that require a hotter shock state are treated as
    ! invalid by the wave curve, so the outer p* bracket can avoid them.
    real(dp), parameter :: riemann_shock_temperature_ceiling = tables_temperature_ceiling
    real(dp), parameter :: riemann_entropy_p_ref = 101325.0_dp
    integer,  parameter :: riemann_rarefaction_quadrature_points = 48
    integer,  parameter :: riemann_max_iterations = 100
    integer,  parameter :: riemann_newton_max_iterations = 30
    real(dp), parameter :: riemann_pressure_derivative_rel_step = 1.0e-5_dp
    real(dp), parameter :: riemann_newton_bracket_eps = 1.0e-12_dp

    type    :: riemann_solver
        ! Legacy/frozen-gamma input state.
        real(dp)    :: rho_l = 0.0_dp, rho_r = 0.0_dp
        real(dp)    :: p_l = 0.0_dp, p_r = 0.0_dp
        real(dp)    :: gamma_l = 0.0_dp, gamma_r = 0.0_dp
        real(dp)    :: u_l = 0.0_dp, u_r = 0.0_dp

        ! Sampled solution state at xi = x/t = 0 by default.
        real(dp)    :: rho = 0.0_dp, p = 0.0_dp, gamma = 0.0_dp, u = 0.0_dp
        real(dp)    :: T = 0.0_dp
        real(dp)    :: u_contact = 0.0_dp
        real(dp)    :: p_star = 0.0_dp

        ! Legacy exact-gamma coefficients.
        real(dp)    :: mu_l = 0.0_dp, mu_r = 0.0_dp
        real(dp)    :: h_l = 0.0_dp, h_r = 0.0_dp
        real(dp)    :: s_l = 0.0_dp, s_r = 0.0_dp
        real(dp)    :: G_l = 0.0_dp, G_r = 0.0_dp

        ! Thermally-perfect frozen-composition mode.
        logical     :: thermally_perfect_mode = .false.
        type(thermophysical_properties), pointer :: thermo_ptr => null()
        integer     :: species_number = 0
        real(dp), dimension(:), allocatable :: Y_l, Y_r, Y
        real(dp)    :: T_l = 0.0_dp, T_r = 0.0_dp
        real(dp)    :: M_l = 0.0_dp, M_r = 0.0_dp
        real(dp)    :: s_tp_l = 0.0_dp, s_tp_r = 0.0_dp
        real(dp)    :: h_tp_l = 0.0_dp, h_tp_r = 0.0_dp
        real(dp)    :: c_l = 0.0_dp, c_r = 0.0_dp

        real(dp)    :: p_floor = 1.0e-12_dp, rho_floor = 1.0e-12_dp
        logical     :: success = .false.
        integer     :: counter = 0
        integer     :: failure_counter = 0
    contains
        ! Backward compatible frozen-gamma API.
        procedure               :: set_parameters

        ! New thermally-perfect API.  The composition is frozen on each side of
        ! the contact; shocks and rarefactions use the JANAF/thermo object.
        procedure               :: set_thermally_perfect_parameters

        procedure               :: clear
        procedure               :: solve
        procedure, private      :: solve_gamma
        procedure, private      :: solve_thermally_perfect

        ! Legacy frozen-gamma wave curves.
        procedure, private      :: F_lr
        procedure, private      :: F_l
        procedure, private      :: F_r
        procedure, private      :: a_l
        procedure, private      :: a_r

        ! Thermally-perfect wave curves and state functions.
        procedure, private      :: F_lr_tp
        procedure, private      :: F_l_tp
        procedure, private      :: F_r_tp
        procedure, private      :: rarefaction_integral_tp
        procedure, private      :: shock_density_tp
        procedure, private      :: shock_temperature_tp
        procedure, private      :: isentropic_temperature_tp
        procedure, private      :: temperature_from_entropy_pressure_Y
        procedure, private      :: mixture_molar_mass_from_Y
        procedure, private      :: mole_fractions_from_Y
        procedure, private      :: normalize_Y
        procedure, private      :: specific_entropy_tp
        procedure, private      :: specific_enthalpy_tp
        procedure, private      :: specific_cp_tp
        procedure, private      :: gamma_from_TY
        procedure, private      :: sound_speed_tp
        procedure, private      :: temperature_from_p_rho_Y
        procedure, private      :: sample_thermally_perfect
        procedure, private      :: sample_gamma
        procedure, private      :: fan_pressure_left_tp
        procedure, private      :: fan_pressure_right_tp
        procedure, private      :: finite_positive

        ! Public sampling/reference API.
        procedure               :: sample_solution
        procedure               :: write_shock_tube_reference

        procedure               :: get_density
        procedure               :: get_pressure
        procedure               :: get_gamma
        procedure               :: get_velocity
        procedure               :: get_temperature
        procedure               :: get_contact_velocity
        procedure               :: get_star_pressure
        procedure               :: get_mass_fractions
        procedure               :: get_activation_count
        procedure               :: get_failure_count
        procedure               :: get_success
    end type

    interface riemann_solver_c
        module procedure constructor
    end interface

contains

    type(riemann_solver) function constructor (manager)
        type(data_manager), intent(inout) :: manager
        call constructor%clear(reset_counter=.true.)
    end function constructor


    subroutine set_parameters(this, rho_l, rho_r, p_l, p_r, gamma_l, gamma_r, u_l, u_r, p_floor, rho_floor)
        class(riemann_solver), intent(inout) :: this
        real(dp), intent(in) :: rho_l, rho_r, p_l, p_r, gamma_l, gamma_r, u_l, u_r
        real(dp), intent(in), optional :: p_floor, rho_floor

        this%thermally_perfect_mode = .false.
        nullify(this%thermo_ptr)
        this%species_number = 0
        if (allocated(this%Y_l)) deallocate(this%Y_l)
        if (allocated(this%Y_r)) deallocate(this%Y_r)
        if (allocated(this%Y))   deallocate(this%Y)

        if (present(p_floor)) then
            this%p_floor = max(p_floor, tiny(1.0_dp))
        else
            this%p_floor = max(1.0e-12_dp*max(1.0_dp, abs(p_l), abs(p_r)), tiny(1.0_dp))
        end if

        if (present(rho_floor)) then
            this%rho_floor = max(rho_floor, tiny(1.0_dp))
        else
            this%rho_floor = max(1.0e-12_dp*max(1.0_dp, abs(rho_l), abs(rho_r)), tiny(1.0_dp))
        end if

        this%rho_l = max(rho_l, this%rho_floor)
        this%rho_r = max(rho_r, this%rho_floor)
        this%p_l = max(p_l, this%p_floor)
        this%p_r = max(p_r, this%p_floor)
        this%gamma_l = max(gamma_l, 1.000001_dp)
        this%gamma_r = max(gamma_r, 1.000001_dp)
        this%u_l = u_l
        this%u_r = u_r

        this%mu_l = (this%gamma_l - 1.0_dp) / (2.0_dp * this%gamma_l)
        this%mu_r = (this%gamma_r - 1.0_dp) / (2.0_dp * this%gamma_r)
        this%h_l = (this%gamma_l - 1.0_dp) / (this%gamma_l + 1.0_dp)
        this%h_r = (this%gamma_r - 1.0_dp) / (this%gamma_r + 1.0_dp)
        this%s_l = this%p_l * this%rho_l**(-this%gamma_l)
        this%s_r = this%p_r * this%rho_r**(-this%gamma_r)
        this%G_l = 2.0_dp * sqrt(this%gamma_l) * this%s_l**(0.5_dp/this%gamma_l) / (this%gamma_l - 1.0_dp)
        this%G_r = 2.0_dp * sqrt(this%gamma_r) * this%s_r**(0.5_dp/this%gamma_r) / (this%gamma_r - 1.0_dp)

        this%rho = 0.0_dp
        this%p = 0.0_dp
        this%gamma = 0.0_dp
        this%u = 0.0_dp
        this%T = 0.0_dp
        this%u_contact = 0.0_dp
        this%p_star = 0.0_dp
        this%success = .false.
    end subroutine set_parameters


    subroutine set_thermally_perfect_parameters(this, thermo, rho_l, rho_r, p_l, p_r, u_l, u_r, Y_l, Y_r, p_floor, rho_floor)
        class(riemann_solver), intent(inout) :: this
        type(thermophysical_properties), target, intent(in) :: thermo
        real(dp), intent(in) :: rho_l, rho_r, p_l, p_r, u_l, u_r
        real(dp), dimension(:), intent(in) :: Y_l, Y_r
        real(dp), intent(in), optional :: p_floor, rho_floor

        this%thermally_perfect_mode = .true.
        this%thermo_ptr => thermo
        this%species_number = min(size(Y_l), size(Y_r))

        if (allocated(this%Y_l)) deallocate(this%Y_l)
        if (allocated(this%Y_r)) deallocate(this%Y_r)
        if (allocated(this%Y))   deallocate(this%Y)
        allocate(this%Y_l(this%species_number), this%Y_r(this%species_number), this%Y(this%species_number))
        this%Y_l = Y_l(1:this%species_number)
        this%Y_r = Y_r(1:this%species_number)
        call this%normalize_Y(this%Y_l)
        call this%normalize_Y(this%Y_r)
        this%Y = 0.0_dp

        if (present(p_floor)) then
            this%p_floor = max(p_floor, tiny(1.0_dp))
        else
            this%p_floor = max(1.0e-12_dp*max(1.0_dp, abs(p_l), abs(p_r)), tiny(1.0_dp))
        end if

        if (present(rho_floor)) then
            this%rho_floor = max(rho_floor, tiny(1.0_dp))
        else
            this%rho_floor = max(1.0e-12_dp*max(1.0_dp, abs(rho_l), abs(rho_r)), tiny(1.0_dp))
        end if

        this%rho_l = max(rho_l, this%rho_floor)
        this%rho_r = max(rho_r, this%rho_floor)
        this%p_l = max(p_l, this%p_floor)
        this%p_r = max(p_r, this%p_floor)
        this%u_l = u_l
        this%u_r = u_r

        this%M_l = this%mixture_molar_mass_from_Y(this%Y_l)
        this%M_r = this%mixture_molar_mass_from_Y(this%Y_r)
        this%T_l = this%temperature_from_p_rho_Y(this%p_l, this%rho_l, this%Y_l)
        this%T_r = this%temperature_from_p_rho_Y(this%p_r, this%rho_r, this%Y_r)
        this%gamma_l = this%gamma_from_TY(this%T_l, this%Y_l)
        this%gamma_r = this%gamma_from_TY(this%T_r, this%Y_r)
        this%c_l = this%sound_speed_tp(this%T_l, this%p_l, this%Y_l)
        this%c_r = this%sound_speed_tp(this%T_r, this%p_r, this%Y_r)
        this%s_tp_l = this%specific_entropy_tp(this%T_l, this%p_l, this%Y_l)
        this%s_tp_r = this%specific_entropy_tp(this%T_r, this%p_r, this%Y_r)
        this%h_tp_l = this%specific_enthalpy_tp(this%T_l, this%Y_l)
        this%h_tp_r = this%specific_enthalpy_tp(this%T_r, this%Y_r)

        ! Fill legacy fields as a fallback/initial estimate only.
        this%mu_l = (this%gamma_l - 1.0_dp) / (2.0_dp * this%gamma_l)
        this%mu_r = (this%gamma_r - 1.0_dp) / (2.0_dp * this%gamma_r)
        this%h_l = (this%gamma_l - 1.0_dp) / (this%gamma_l + 1.0_dp)
        this%h_r = (this%gamma_r - 1.0_dp) / (this%gamma_r + 1.0_dp)
        this%s_l = this%p_l * this%rho_l**(-this%gamma_l)
        this%s_r = this%p_r * this%rho_r**(-this%gamma_r)
        this%G_l = 2.0_dp * sqrt(this%gamma_l) * this%s_l**(0.5_dp/this%gamma_l) / max(this%gamma_l - 1.0_dp, tiny(1.0_dp))
        this%G_r = 2.0_dp * sqrt(this%gamma_r) * this%s_r**(0.5_dp/this%gamma_r) / max(this%gamma_r - 1.0_dp, tiny(1.0_dp))

        this%rho = 0.0_dp
        this%p = 0.0_dp
        this%gamma = 0.0_dp
        this%u = 0.0_dp
        this%T = 0.0_dp
        this%u_contact = 0.0_dp
        this%p_star = 0.0_dp
        this%success = .false.
    end subroutine set_thermally_perfect_parameters


    subroutine clear(this, reset_counter)
        class(riemann_solver), intent(inout) :: this
        logical, intent(in), optional :: reset_counter

        this%rho_l = 0.0_dp
        this%rho_r = 0.0_dp
        this%p_l = 0.0_dp
        this%p_r = 0.0_dp
        this%gamma_l = 0.0_dp
        this%gamma_r = 0.0_dp
        this%u_l = 0.0_dp
        this%u_r = 0.0_dp
        this%mu_l = 0.0_dp
        this%mu_r = 0.0_dp
        this%h_l = 0.0_dp
        this%h_r = 0.0_dp
        this%s_l = 0.0_dp
        this%s_r = 0.0_dp
        this%G_l = 0.0_dp
        this%G_r = 0.0_dp
        this%rho = 0.0_dp
        this%p = 0.0_dp
        this%gamma = 0.0_dp
        this%u = 0.0_dp
        this%T = 0.0_dp
        this%u_contact = 0.0_dp
        this%p_star = 0.0_dp
        this%thermally_perfect_mode = .false.
        nullify(this%thermo_ptr)
        this%species_number = 0
        this%T_l = 0.0_dp
        this%T_r = 0.0_dp
        this%M_l = 0.0_dp
        this%M_r = 0.0_dp
        this%s_tp_l = 0.0_dp
        this%s_tp_r = 0.0_dp
        this%h_tp_l = 0.0_dp
        this%h_tp_r = 0.0_dp
        this%c_l = 0.0_dp
        this%c_r = 0.0_dp
        if (allocated(this%Y_l)) deallocate(this%Y_l)
        if (allocated(this%Y_r)) deallocate(this%Y_r)
        if (allocated(this%Y))   deallocate(this%Y)
        this%success = .false.
        if (present(reset_counter)) then
            if (reset_counter) then
                this%counter = 0
                this%failure_counter = 0
            end if
        end if
    end subroutine clear


    subroutine solve(this)
        class(riemann_solver), intent(inout) :: this
        if (this%thermally_perfect_mode) then
            call this%solve_thermally_perfect()
        else
            call this%solve_gamma()
        end if
    end subroutine solve


    subroutine solve_gamma(this)
        class(riemann_solver), intent(inout) :: this

        real(dp) :: delta_p, p_0, p_c_old, p_c, u_c
        real(dp) :: func_value, func_deriv_value
        real(dp) :: p_min, p_scale, tol_abs, tol_rel
        integer  :: iter
        logical  :: converged

        this%counter = this%counter + 1
        this%success = .false.

        if (.not. this%finite_positive(this%rho_l) .or. .not. this%finite_positive(this%rho_r) .or. &
            .not. this%finite_positive(this%p_l)   .or. .not. this%finite_positive(this%p_r)) then
            this%failure_counter = this%failure_counter + 1
            return
        end if

        p_min = max(this%p_floor, tiny(1.0_dp))
        p_scale = max(this%p_l, this%p_r, p_min)

        this%c_l = sqrt(max(this%gamma_l*this%p_l/this%rho_l, 0.0_dp))
        this%c_r = sqrt(max(this%gamma_r*this%p_r/this%rho_r, 0.0_dp))
        p_0 = 0.5_dp*(this%p_l + this%p_r) - 0.125_dp*(this%u_r - this%u_l)*(this%rho_l + this%rho_r)*(this%c_l + this%c_r)
        p_0 = max(p_min, p_0)

        tol_rel = 1.0e-8_dp
        tol_abs = 1.0e-10_dp*p_scale
        converged = .false.
        p_c = p_0

        do iter = 1, 50
            delta_p = max(1.0e-6_dp*max(abs(p_0), p_min), p_min)
            func_value = this%F_lr(p_0)
            func_deriv_value = (this%F_lr(p_0 + delta_p) - this%F_lr(max(p_min, p_0 - delta_p))) / &
                               (p_0 + delta_p - max(p_min, p_0 - delta_p))

            if (func_deriv_value /= func_deriv_value .or. abs(func_deriv_value) <= tiny(1.0_dp)) exit

            p_c_old = p_0
            p_c = p_0 - func_value/func_deriv_value
            if (p_c /= p_c .or. p_c <= p_min) then
                p_c = 0.5_dp*(p_0 + p_min)
            end if
            p_c = max(p_min, p_c)

            if (abs(p_c - p_c_old) <= max(tol_abs, tol_rel*max(abs(p_c), p_min))) then
                converged = .true.
                exit
            end if
            p_0 = p_c
        end do

        if (.not. converged) then
            this%failure_counter = this%failure_counter + 1
            return
        end if

        this%p_star = p_c
        this%u_contact = this%F_l(p_c)
        if (this%u_contact /= this%u_contact) then
            this%failure_counter = this%failure_counter + 1
            return
        end if

        call this%sample_gamma(0.0_dp, this%rho, this%u, this%p, this%gamma)
        this%T = 0.0_dp
        this%success = this%finite_positive(this%rho) .and. this%finite_positive(this%p) .and. (this%u == this%u)
        if (.not. this%success) this%failure_counter = this%failure_counter + 1
    end subroutine solve_gamma


    subroutine solve_thermally_perfect(this)
        class(riemann_solver), intent(inout) :: this

        real(dp) :: p_low, p_high, p_trial, p_new
        real(dp) :: p_candidate, p_bad
        real(dp) :: f_low, f_high, f_trial, f_new, f_candidate
        real(dp) :: p_scale, tol_abs, tol_rel, vel_scale
        real(dp) :: p_pvrs, dp_fd, p_minus, p_plus, f_minus, f_plus, deriv
        integer :: iter, iter_bracket
        logical :: bracketed, use_newton, converged

        this%counter = this%counter + 1
        this%success = .false.

        if (.not. associated(this%thermo_ptr) .or. this%species_number <= 0) then
            this%failure_counter = this%failure_counter + 1
            return
        end if
        if (.not. this%finite_positive(this%rho_l) .or. .not. this%finite_positive(this%rho_r) .or. &
            .not. this%finite_positive(this%p_l)   .or. .not. this%finite_positive(this%p_r)) then
            this%failure_counter = this%failure_counter + 1
            return
        end if

        p_scale = max(this%p_l, this%p_r, this%p_floor)
        vel_scale = max(this%c_l, this%c_r, abs(this%u_l), abs(this%u_r), 1.0_dp)
        tol_rel = 1.0e-8_dp
        tol_abs = 1.0e-10_dp*p_scale

        ! PVRS estimate provides the first Newton point and a useful upper-bound scale.
        p_pvrs = 0.5_dp*(this%p_l + this%p_r) - &
                 0.125_dp*(this%u_r - this%u_l)*(this%rho_l + this%rho_r)*(this%c_l + this%c_r)
        p_pvrs = max(this%p_floor, p_pvrs)

        p_low = max(this%p_floor, 1.0e-12_dp*p_scale)
        f_low = this%F_lr_tp(p_low)
        if (f_low /= f_low .or. f_low < 0.0_dp) then
            this%failure_counter = this%failure_counter + 1
            return
        end if

        ! Build the upper pressure bracket cautiously.  Starting from
        ! 2*max(p_l,p_r) is too aggressive for light-gas/heavy-gas problems:
        ! for H2(250 atm)->air(1 atm), it asks the right shock Hugoniot for a
        ! 500-atm trial shock and drives the thermally-perfect temperature far
        ! beyond the valid JANAF range.  We therefore expand upward from the
        ! lower side and treat high-temperature/NaN wave-curve evaluations as
        ! invalid trial pressures, not as fatal Riemann failures.
        p_high = max(min(this%p_l, this%p_r), 10.0_dp*p_low, this%p_floor)
        f_high = this%F_lr_tp(p_high)
        bracketed = .false.
        if (f_high == f_high) then
            if (f_high <= 0.0_dp) then
                bracketed = .true.
            else
                p_low = p_high
                f_low = f_high
            end if
        else
            this%failure_counter = this%failure_counter + 1
            return
        end if

        do iter = 1, 80
            if (bracketed) exit

            p_candidate = min(2.0_dp*p_low, 1.0e12_dp*p_scale)
            if (p_candidate <= p_low*(1.0_dp + 100.0_dp*epsilon(1.0_dp))) exit

            f_candidate = this%F_lr_tp(p_candidate)
            if (f_candidate == f_candidate) then
                if (f_candidate <= 0.0_dp) then
                    p_high = p_candidate
                    f_high = f_candidate
                    bracketed = .true.
                else
                    p_low = p_candidate
                    f_low = f_candidate
                end if
            else
                ! The trial pressure is beyond the thermodynamic validity of
                ! at least one shock curve.  Search between the last finite
                ! positive point and the invalid point for a finite negative
                ! value.  If none exists, the physical p* is outside the valid
                ! thermodynamic range of this Riemann solve.
                p_bad = p_candidate
                do iter_bracket = 1, 80
                    p_candidate = sqrt(max(p_low, this%p_floor)*max(p_bad, this%p_floor))
                    if (p_candidate <= p_low*(1.0_dp + 100.0_dp*epsilon(1.0_dp)) .or. &
                        p_candidate >= p_bad*(1.0_dp - 100.0_dp*epsilon(1.0_dp))) then
                        p_candidate = 0.5_dp*(p_low + p_bad)
                    end if
                    if (p_candidate <= p_low .or. p_candidate >= p_bad) exit

                    f_candidate = this%F_lr_tp(p_candidate)
                    if (f_candidate == f_candidate) then
                        if (f_candidate <= 0.0_dp) then
                            p_high = p_candidate
                            f_high = f_candidate
                            bracketed = .true.
                            exit
                        else
                            p_low = p_candidate
                            f_low = f_candidate
                        end if
                    else
                        p_bad = p_candidate
                    end if
                end do
                exit
            end if
        end do

        if (.not. bracketed) then
            this%failure_counter = this%failure_counter + 1
            return
        end if

        p_trial = min(max(p_pvrs, p_low), p_high)
        f_trial = this%F_lr_tp(p_trial)
        if (f_trial /= f_trial) then
            p_trial = 0.5_dp*(p_low + p_high)
            f_trial = this%F_lr_tp(p_trial)
            if (f_trial /= f_trial) then
                this%failure_counter = this%failure_counter + 1
                return
            end if
        end if

        converged = .false.
        do iter = 1, riemann_max_iterations
            if (abs(f_trial) <= 1.0e-10_dp*vel_scale) then
                converged = .true.
                exit
            end if
            if (abs(p_high - p_low) <= max(tol_abs, tol_rel*max(abs(p_trial), this%p_floor))) then
                converged = .true.
                exit
            end if

            ! Safeguarded Newton step.  F_lr_tp is monotonically decreasing in p,
            ! but the derivative is evaluated numerically to avoid mistakes at
            ! shock/rarefaction transitions and for variable thermodynamics.
            use_newton = .false.
            dp_fd = max(riemann_pressure_derivative_rel_step*max(abs(p_trial), p_scale), &
                        1.0e-12_dp*p_scale)
            p_minus = max(p_low,  p_trial - dp_fd)
            p_plus  = min(p_high, p_trial + dp_fd)
            if (p_plus > p_minus) then
                f_minus = this%F_lr_tp(p_minus)
                f_plus  = this%F_lr_tp(p_plus)
                if (f_minus == f_minus .and. f_plus == f_plus) then
                    deriv = (f_plus - f_minus)/(p_plus - p_minus)
                    if (deriv == deriv .and. abs(deriv) > tiny(1.0_dp)) then
                        p_new = p_trial - f_trial/deriv
                        if (p_new == p_new .and. p_new > p_low  + riemann_newton_bracket_eps*(p_high-p_low) .and. &
                                            p_new < p_high - riemann_newton_bracket_eps*(p_high-p_low)) then
                            use_newton = .true.
                        end if
                    end if
                end if
            end if

            if (.not. use_newton) p_new = 0.5_dp*(p_low + p_high)
            f_new = this%F_lr_tp(p_new)
            if (f_new /= f_new) then
                p_new = 0.5_dp*(p_low + p_high)
                f_new = this%F_lr_tp(p_new)
                if (f_new /= f_new) then
                    this%failure_counter = this%failure_counter + 1
                    return
                end if
            end if

            if (f_new > 0.0_dp) then
                p_low = p_new
                f_low = f_new
            else
                p_high = p_new
                f_high = f_new
            end if

            p_trial = p_new
            f_trial = f_new
        end do

        if (.not. converged) then
            ! The bracket itself may already be tight enough even if the residual
            ! tolerance was not reached in a difficult near-vacuum/strong-shock case.
            if (abs(p_high - p_low) <= max(10.0_dp*tol_abs, 10.0_dp*tol_rel*max(abs(p_trial), this%p_floor))) then
                converged = .true.
            end if
        end if

        if (.not. converged) then
            this%failure_counter = this%failure_counter + 1
            return
        end if

        this%p_star = max(p_trial, this%p_floor)
        this%u_contact = this%F_l_tp(this%p_star)
        if (this%u_contact /= this%u_contact) then
            this%failure_counter = this%failure_counter + 1
            return
        end if

        call this%sample_thermally_perfect(0.0_dp, this%rho, this%u, this%p, this%T, this%gamma, this%Y)
        this%success = this%finite_positive(this%rho) .and. this%finite_positive(this%p) .and. &
                       this%finite_positive(this%T) .and. (this%u == this%u)
        if (.not. this%success) this%failure_counter = this%failure_counter + 1
    end subroutine solve_thermally_perfect


    logical function finite_positive(this, value)
        class(riemann_solver), intent(in) :: this
        real(dp), intent(in) :: value
        finite_positive = (value == value) .and. (value > 0.0_dp)
    end function finite_positive


    function F_lr(this, p) result(F_lr_value)
        class(riemann_solver), intent(in) :: this
        real(dp), intent(in) :: p
        real(dp) :: F_lr_value
        F_lr_value = this%F_l(p) - this%F_r(p)
    end function F_lr


    function F_l(this, p) result(F_l_value)
        class(riemann_solver), intent(in) :: this
        real(dp), intent(in) :: p
        real(dp) :: F_l_value
        real(dp) :: pressure
        pressure = max(p, this%p_floor)
        if (pressure < this%p_l) then
            F_l_value = this%u_l - this%G_l * (pressure**this%mu_l - this%p_l**this%mu_l)
        else
            F_l_value = this%u_l - sqrt(max((pressure - this%p_l) * (1.0_dp - this%a_l(pressure)) / this%rho_l, 0.0_dp))
        end if
    end function F_l


    function F_r(this, p) result(F_r_value)
        class(riemann_solver), intent(in) :: this
        real(dp), intent(in) :: p
        real(dp) :: F_r_value
        real(dp) :: pressure
        pressure = max(p, this%p_floor)
        if (pressure < this%p_r) then
            F_r_value = this%u_r + this%G_r * (pressure**this%mu_r - this%p_r**this%mu_r)
        else
            F_r_value = this%u_r + sqrt(max((pressure - this%p_r) * (1.0_dp - this%a_r(pressure)) / this%rho_r, 0.0_dp))
        end if
    end function F_r


    function a_l(this, p) result(a_l_value)
        class(riemann_solver), intent(in) :: this
        real(dp), intent(in) :: p
        real(dp) :: a_l_value
        real(dp) :: pressure
        pressure = max(p, this%p_floor)
        a_l_value = (this%p_l + this%h_l * pressure) / max(pressure + this%h_l * this%p_l, tiny(1.0_dp))
    end function a_l


    function a_r(this, p) result(a_r_value)
        class(riemann_solver), intent(in) :: this
        real(dp), intent(in) :: p
        real(dp) :: a_r_value
        real(dp) :: pressure
        pressure = max(p, this%p_floor)
        a_r_value = (this%p_r + this%h_r * pressure) / max(pressure + this%h_r * this%p_r, tiny(1.0_dp))
    end function a_r


    function F_lr_tp(this, p) result(value)
        class(riemann_solver), intent(in) :: this
        real(dp), intent(in) :: p
        real(dp) :: value
        value = this%F_l_tp(p) - this%F_r_tp(p)
    end function F_lr_tp


    function F_l_tp(this, p) result(value)
        class(riemann_solver), intent(in) :: this
        real(dp), intent(in) :: p
        real(dp) :: value, pressure, rho_star, dv
        pressure = max(p, this%p_floor)
        if (pressure < this%p_l) then
            value = this%u_l - this%rarefaction_integral_tp(.true., this%p_l, pressure)
        else
            rho_star = this%shock_density_tp(.true., pressure)
            if (rho_star /= rho_star) then
                value = ieee_value(value, ieee_quiet_nan)
                return
            end if
            dv = max(1.0_dp/this%rho_l - 1.0_dp/max(rho_star,this%rho_floor), 0.0_dp)
            value = this%u_l - sqrt(max((pressure - this%p_l)*dv, 0.0_dp))
        end if
    end function F_l_tp


    function F_r_tp(this, p) result(value)
        class(riemann_solver), intent(in) :: this
        real(dp), intent(in) :: p
        real(dp) :: value, pressure, rho_star, dv
        pressure = max(p, this%p_floor)
        if (pressure < this%p_r) then
            value = this%u_r + this%rarefaction_integral_tp(.false., this%p_r, pressure)
        else
            rho_star = this%shock_density_tp(.false., pressure)
            if (rho_star /= rho_star) then
                value = ieee_value(value, ieee_quiet_nan)
                return
            end if
            dv = max(1.0_dp/this%rho_r - 1.0_dp/max(rho_star,this%rho_floor), 0.0_dp)
            value = this%u_r + sqrt(max((pressure - this%p_r)*dv, 0.0_dp))
        end if
    end function F_r_tp


    function rarefaction_integral_tp(this, use_left, p_from, p_to) result(value)
        class(riemann_solver), intent(in) :: this
        logical, intent(in) :: use_left
        real(dp), intent(in) :: p_from, p_to
        real(dp) :: value
        real(dp) :: log_a, log_b, dlogp, pressure, dp_seg, temp, rho_loc, c_loc
        integer :: n

        value = 0.0_dp
        if (abs(p_to - p_from) <= tiny(1.0_dp)*max(abs(p_from),abs(p_to),1.0_dp)) return

        log_a = log(max(p_from, this%p_floor))
        log_b = log(max(p_to, this%p_floor))
        dlogp = (log_b - log_a) / real(riemann_rarefaction_quadrature_points, dp)
        do n = 1, riemann_rarefaction_quadrature_points
            pressure = exp(log_a + (real(n,dp) - 0.5_dp)*dlogp)
            if (use_left) then
                temp = this%isentropic_temperature_tp(.true., pressure)
                rho_loc = pressure*this%M_l/(r_gase_J*max(temp,riemann_temperature_floor))
                c_loc = this%sound_speed_tp(temp, pressure, this%Y_l)
            else
                temp = this%isentropic_temperature_tp(.false., pressure)
                rho_loc = pressure*this%M_r/(r_gase_J*max(temp,riemann_temperature_floor))
                c_loc = this%sound_speed_tp(temp, pressure, this%Y_r)
            end if
            dp_seg = pressure*dlogp
            value = value + dp_seg/(max(rho_loc,this%rho_floor)*max(c_loc,tiny(1.0_dp)))
        end do
    end function rarefaction_integral_tp


    function shock_density_tp(this, use_left, pressure) result(rho_star)
        class(riemann_solver), intent(in) :: this
        logical, intent(in) :: use_left
        real(dp), intent(in) :: pressure
        real(dp) :: rho_star, T_star, M_mix
        T_star = this%shock_temperature_tp(use_left, pressure)
        if (T_star /= T_star) then
            rho_star = ieee_value(rho_star, ieee_quiet_nan)
            return
        end if
        if (use_left) then
            M_mix = this%M_l
        else
            M_mix = this%M_r
        end if
        rho_star = max(pressure, this%p_floor)*M_mix/(r_gase_J*max(T_star, riemann_temperature_floor))
        rho_star = max(rho_star, this%rho_floor)
    end function shock_density_tp


    function shock_temperature_tp(this, use_left, pressure) result(T_star)
        class(riemann_solver), intent(in) :: this
        logical, intent(in) :: use_left
        real(dp), intent(in) :: pressure
        real(dp) :: T_star
        real(dp) :: T1, rho1, p1, h1, M_mix, T_low, T_high, T_trial, T_new
        real(dp) :: f_low, f_high, f_trial, f_new, deriv, heat_capacity
        real(dp), allocatable :: Y_side(:)
        integer :: iter
        logical :: use_newton

        T_star = ieee_value(T_star, ieee_quiet_nan)
        allocate(Y_side(this%species_number))
        if (use_left) then
            T1 = this%T_l; rho1 = this%rho_l; p1 = this%p_l; h1 = this%h_tp_l; M_mix = this%M_l
            Y_side = this%Y_l
        else
            T1 = this%T_r; rho1 = this%rho_r; p1 = this%p_r; h1 = this%h_tp_r; M_mix = this%M_r
            Y_side = this%Y_r
        end if

        if (pressure <= p1) then
            T_star = this%isentropic_temperature_tp(use_left, pressure)
            deallocate(Y_side)
            return
        end if

        T_low = riemann_temperature_floor
        T_high = min(max(1.2_dp*T1 + 1.0_dp, T1 + 10.0_dp), riemann_shock_temperature_ceiling)
        if (T_high <= T_low) then
            deallocate(Y_side)
            return
        end if

        f_low = this%specific_enthalpy_tp(T_low, Y_side) - h1 - &
                0.5_dp*(pressure - p1)*(1.0_dp/rho1 + r_gase_J*T_low/(pressure*M_mix))
        f_high = this%specific_enthalpy_tp(T_high, Y_side) - h1 - &
                 0.5_dp*(pressure - p1)*(1.0_dp/rho1 + r_gase_J*T_high/(pressure*M_mix))

        do iter = 1, 80
            if (f_high >= 0.0_dp .or. f_high /= f_high) exit
            if (T_high >= riemann_shock_temperature_ceiling) exit
            T_high = min(2.0_dp*T_high, riemann_shock_temperature_ceiling)
            f_high = this%specific_enthalpy_tp(T_high, Y_side) - h1 - &
                     0.5_dp*(pressure - p1)*(1.0_dp/rho1 + r_gase_J*T_high/(pressure*M_mix))
        end do

        ! If the Hugoniot temperature is above the validity ceiling, this trial
        ! pressure is not usable by the thermally-perfect wave curve.  Return
        ! NaN so the outer pressure solver can lower/avoid this trial pressure.
        if (f_low > 0.0_dp .or. f_high < 0.0_dp .or. f_low /= f_low .or. f_high /= f_high) then
            deallocate(Y_side)
            return
        end if

        T_trial = min(max(max(T1, 0.5_dp*(T_low + T_high)), T_low), T_high)
        f_trial = this%specific_enthalpy_tp(T_trial, Y_side) - h1 - &
                  0.5_dp*(pressure - p1)*(1.0_dp/rho1 + r_gase_J*T_trial/(pressure*M_mix))
        if (f_trial /= f_trial) then
            T_trial = 0.5_dp*(T_low + T_high)
            f_trial = this%specific_enthalpy_tp(T_trial, Y_side) - h1 - &
                      0.5_dp*(pressure - p1)*(1.0_dp/rho1 + r_gase_J*T_trial/(pressure*M_mix))
        end if
        if (f_trial /= f_trial) then
            deallocate(Y_side)
            return
        end if

        do iter = 1, riemann_newton_max_iterations
            if (abs(f_trial) <= 1.0e-8_dp*max(abs(h1), 1.0_dp)) exit
            if (abs(T_high - T_low) <= 1.0e-8_dp*max(T_trial, 1.0_dp)) exit

            heat_capacity = this%specific_cp_tp(T_trial, Y_side)
            deriv = heat_capacity - 0.5_dp*(pressure - p1)*r_gase_J/(pressure*M_mix)
            use_newton = .false.
            if (deriv == deriv .and. deriv > tiny(1.0_dp)) then
                T_new = T_trial - f_trial/deriv
                if (T_new == T_new .and. T_new > T_low  + riemann_newton_bracket_eps*(T_high-T_low) .and. &
                                      T_new < T_high - riemann_newton_bracket_eps*(T_high-T_low)) then
                    use_newton = .true.
                end if
            end if
            if (.not. use_newton) T_new = 0.5_dp*(T_low + T_high)

            f_new = this%specific_enthalpy_tp(T_new, Y_side) - h1 - &
                    0.5_dp*(pressure - p1)*(1.0_dp/rho1 + r_gase_J*T_new/(pressure*M_mix))
            if (f_new /= f_new) then
                T_new = 0.5_dp*(T_low + T_high)
                f_new = this%specific_enthalpy_tp(T_new, Y_side) - h1 - &
                        0.5_dp*(pressure - p1)*(1.0_dp/rho1 + r_gase_J*T_new/(pressure*M_mix))
                if (f_new /= f_new) then
                    deallocate(Y_side)
                    return
                end if
            end if

            if (f_new > 0.0_dp) then
                T_high = T_new
                f_high = f_new
            else
                T_low = T_new
                f_low = f_new
            end if
            T_trial = T_new
            f_trial = f_new
        end do

        T_star = min(max(T_trial, riemann_temperature_floor), riemann_shock_temperature_ceiling)
        deallocate(Y_side)
    end function shock_temperature_tp


    function isentropic_temperature_tp(this, use_left, pressure) result(temp)
        class(riemann_solver), intent(in) :: this
        logical, intent(in) :: use_left
        real(dp), intent(in) :: pressure
        real(dp) :: temp
        if (use_left) then
            temp = this%temperature_from_entropy_pressure_Y(this%s_tp_l, pressure, this%Y_l, this%T_l)
        else
            temp = this%temperature_from_entropy_pressure_Y(this%s_tp_r, pressure, this%Y_r, this%T_r)
        end if
    end function isentropic_temperature_tp


    function temperature_from_entropy_pressure_Y(this, s_target, pressure, Y_mass, T_guess) result(T_state)
        class(riemann_solver), intent(in) :: this
        real(dp), intent(in) :: s_target, pressure, T_guess
        real(dp), dimension(:), intent(in) :: Y_mass
        real(dp) :: T_state
        real(dp) :: T_low, T_high, T_trial, T_new
        real(dp) :: f_low, f_high, f_trial, f_new, deriv, heat_capacity
        integer :: iter
        logical :: use_newton

        T_low = riemann_temperature_floor
        T_high = riemann_temperature_ceiling
        f_low = this%specific_entropy_tp(T_low, pressure, Y_mass) - s_target
        f_high = this%specific_entropy_tp(T_high, pressure, Y_mass) - s_target

        ! For normal thermally-perfect gases s(T,p,Y) is monotone increasing in T.
        ! If the requested entropy is outside the supported temperature interval,
        ! return the nearest bound rather than iterating on an invalid bracket.
        if (f_low >= 0.0_dp) then
            T_state = T_low
            return
        end if
        if (f_high <= 0.0_dp) then
            T_state = T_high
            return
        end if

        T_trial = min(max(T_guess, T_low), T_high)
        f_trial = this%specific_entropy_tp(T_trial, pressure, Y_mass) - s_target
        if (f_trial /= f_trial) then
            T_trial = 0.5_dp*(T_low + T_high)
            f_trial = this%specific_entropy_tp(T_trial, pressure, Y_mass) - s_target
        end if

        do iter = 1, riemann_newton_max_iterations
            if (abs(f_trial) <= 1.0e-9_dp*max(abs(s_target), 1.0_dp)) exit
            if (abs(T_high - T_low) <= 1.0e-8_dp*max(T_trial, 1.0_dp)) exit

            heat_capacity = this%specific_cp_tp(T_trial, Y_mass)
            deriv = heat_capacity/max(T_trial, riemann_temperature_floor)
            use_newton = .false.
            if (deriv == deriv .and. deriv > tiny(1.0_dp)) then
                T_new = T_trial - f_trial/deriv
                if (T_new == T_new .and. T_new > T_low  + riemann_newton_bracket_eps*(T_high-T_low) .and. &
                                      T_new < T_high - riemann_newton_bracket_eps*(T_high-T_low)) then
                    use_newton = .true.
                end if
            end if
            if (.not. use_newton) T_new = 0.5_dp*(T_low + T_high)

            f_new = this%specific_entropy_tp(T_new, pressure, Y_mass) - s_target
            if (f_new /= f_new) then
                T_new = 0.5_dp*(T_low + T_high)
                f_new = this%specific_entropy_tp(T_new, pressure, Y_mass) - s_target
            end if

            if (f_new > 0.0_dp) then
                T_high = T_new
                f_high = f_new
            else
                T_low = T_new
                f_low = f_new
            end if
            T_trial = T_new
            f_trial = f_new
        end do
        T_state = min(max(T_trial, riemann_temperature_floor), riemann_temperature_ceiling)
    end function temperature_from_entropy_pressure_Y


    function mixture_molar_mass_from_Y(this, Y_mass) result(M_mix)
        class(riemann_solver), intent(in) :: this
        real(dp), dimension(:), intent(in) :: Y_mass
        real(dp) :: M_mix, denom
        integer :: spec, ns
        ns = min(size(Y_mass), size(this%thermo_ptr%molar_masses))
        denom = 0.0_dp
        do spec = 1, ns
            if (this%thermo_ptr%molar_masses(spec) > tiny(1.0_dp)) then
                denom = denom + max(Y_mass(spec), 0.0_dp)/this%thermo_ptr%molar_masses(spec)
            end if
        end do
        M_mix = 1.0_dp/max(denom, tiny(1.0_dp))
    end function mixture_molar_mass_from_Y


    subroutine mole_fractions_from_Y(this, Y_mass, X)
        class(riemann_solver), intent(in) :: this
        real(dp), dimension(:), intent(in) :: Y_mass
        real(dp), dimension(:), intent(out) :: X
        real(dp) :: M_mix
        integer :: spec, ns
        ns = min(size(Y_mass), size(X), size(this%thermo_ptr%molar_masses))
        X = 0.0_dp
        M_mix = this%mixture_molar_mass_from_Y(Y_mass)
        do spec = 1, ns
            if (this%thermo_ptr%molar_masses(spec) > tiny(1.0_dp)) then
                X(spec) = max(Y_mass(spec), 0.0_dp)*M_mix/this%thermo_ptr%molar_masses(spec)
            end if
        end do
        call this%normalize_Y(X)
    end subroutine mole_fractions_from_Y


    subroutine normalize_Y(this, Y_mass)
        class(riemann_solver), intent(in) :: this
        real(dp), dimension(:), intent(inout) :: Y_mass
        real(dp) :: sumY
        integer :: spec
        do spec = 1, size(Y_mass)
            if (Y_mass(spec) /= Y_mass(spec) .or. Y_mass(spec) < 0.0_dp) Y_mass(spec) = 0.0_dp
        end do
        sumY = sum(Y_mass)
        if (sumY > tiny(1.0_dp)) then
            Y_mass = Y_mass/sumY
        else
            Y_mass = 0.0_dp
            if (size(Y_mass) > 0) Y_mass(1) = 1.0_dp
        end if
    end subroutine normalize_Y


    function specific_entropy_tp(this, temperature, pressure, Y_mass) result(s_spec)
        class(riemann_solver), intent(in) :: this
        real(dp), intent(in) :: temperature, pressure
        real(dp), dimension(:), intent(in) :: Y_mass
        real(dp) :: s_spec, s_molar, M_mix, x_spec, T_loc, p_loc
        integer :: spec, ns

        ns = min(size(Y_mass), size(this%thermo_ptr%molar_masses))
        T_loc = min(max(temperature, riemann_temperature_floor), riemann_temperature_ceiling)
        p_loc = max(pressure, this%p_floor)
        M_mix = this%mixture_molar_mass_from_Y(Y_mass)
        s_molar = -r_gase_J*log(p_loc/riemann_entropy_p_ref)
        do spec = 1, ns
            if (this%thermo_ptr%molar_masses(spec) > tiny(1.0_dp)) then
                x_spec = max(Y_mass(spec), 0.0_dp)*M_mix/this%thermo_ptr%molar_masses(spec)
                if (x_spec > tiny(1.0_dp)) then
                    s_molar = s_molar + x_spec*this%thermo_ptr%specie_entropy_molar(T_loc, spec) &
                                      - r_gase_J*x_spec*log(x_spec)
                end if
            end if
        end do
        s_spec = s_molar/max(M_mix, tiny(1.0_dp))
    end function specific_entropy_tp


    function specific_enthalpy_tp(this, temperature, Y_mass) result(h_spec)
        class(riemann_solver), intent(in) :: this
        real(dp), intent(in) :: temperature
        real(dp), dimension(:), intent(in) :: Y_mass
        real(dp) :: h_spec, M_mix, T_loc
        real(dp), dimension(:), allocatable :: X
        T_loc = min(max(temperature, riemann_temperature_floor), riemann_temperature_ceiling)
        M_mix = this%mixture_molar_mass_from_Y(Y_mass)
        allocate(X(size(Y_mass)))
        call this%mole_fractions_from_Y(Y_mass, X)
        h_spec = this%thermo_ptr%mixture_enthalpy_molar(T_loc, X)/max(M_mix, tiny(1.0_dp))
        deallocate(X)
    end function specific_enthalpy_tp


    function specific_cp_tp(this, temperature, Y_mass) result(cp_spec)
        class(riemann_solver), intent(in) :: this
        real(dp), intent(in) :: temperature
        real(dp), dimension(:), intent(in) :: Y_mass
        real(dp) :: cp_spec, cp_molar, M_mix, T_loc
        real(dp), dimension(:), allocatable :: X
        T_loc = min(max(temperature, riemann_temperature_floor), riemann_temperature_ceiling)
        M_mix = this%mixture_molar_mass_from_Y(Y_mass)
        allocate(X(size(Y_mass)))
        call this%mole_fractions_from_Y(Y_mass, X)
        cp_molar = this%thermo_ptr%mixture_cp_molar(T_loc, X)
        cp_spec = cp_molar/max(M_mix, tiny(1.0_dp))
        deallocate(X)
    end function specific_cp_tp


    function gamma_from_TY(this, temperature, Y_mass) result(gamma_value)
        class(riemann_solver), intent(in) :: this
        real(dp), intent(in) :: temperature
        real(dp), dimension(:), intent(in) :: Y_mass
        real(dp) :: gamma_value, cp_molar, cv_molar, T_loc
        real(dp), dimension(:), allocatable :: X
        T_loc = min(max(temperature, riemann_temperature_floor), riemann_temperature_ceiling)
        allocate(X(size(Y_mass)))
        call this%mole_fractions_from_Y(Y_mass, X)
        cp_molar = this%thermo_ptr%mixture_cp_molar(T_loc, X)
        cv_molar = max(cp_molar - r_gase_J, tiny(1.0_dp))
        gamma_value = max(cp_molar/cv_molar, 1.000001_dp)
        deallocate(X)
    end function gamma_from_TY


    function sound_speed_tp(this, temperature, pressure, Y_mass) result(c)
        class(riemann_solver), intent(in) :: this
        real(dp), intent(in) :: temperature, pressure
        real(dp), dimension(:), intent(in) :: Y_mass
        real(dp) :: c, gamma_value, rho_value, M_mix, T_loc
        T_loc = min(max(temperature, riemann_temperature_floor), riemann_temperature_ceiling)
        M_mix = this%mixture_molar_mass_from_Y(Y_mass)
        rho_value = max(pressure, this%p_floor)*M_mix/(r_gase_J*T_loc)
        gamma_value = this%gamma_from_TY(T_loc, Y_mass)
        c = sqrt(max(gamma_value*max(pressure,this%p_floor)/max(rho_value,this%rho_floor), 0.0_dp))
    end function sound_speed_tp


    function temperature_from_p_rho_Y(this, pressure, density, Y_mass) result(temp)
        class(riemann_solver), intent(in) :: this
        real(dp), intent(in) :: pressure, density
        real(dp), dimension(:), intent(in) :: Y_mass
        real(dp) :: temp, M_mix
        M_mix = this%mixture_molar_mass_from_Y(Y_mass)
        temp = max(pressure, this%p_floor)*M_mix/(r_gase_J*max(density, this%rho_floor))
        temp = min(max(temp, riemann_temperature_floor), riemann_temperature_ceiling)
    end function temperature_from_p_rho_Y


    function fan_pressure_left_tp(this, xi) result(pressure)
        class(riemann_solver), intent(in) :: this
        real(dp), intent(in) :: xi
        real(dp) :: pressure, p_low, p_high, p_mid, g_mid
        integer :: iter
        p_low = min(this%p_star, this%p_l)
        p_high = max(this%p_star, this%p_l)
        p_mid = 0.5_dp*(p_low + p_high)
        do iter = 1, 80
            p_mid = 0.5_dp*(p_low + p_high)
            g_mid = this%F_l_tp(p_mid) - this%sound_speed_tp(this%isentropic_temperature_tp(.true., p_mid), p_mid, this%Y_l) - xi
            if (abs(g_mid) <= 1.0e-8_dp*max(this%c_l,1.0_dp)) exit
            if (g_mid > 0.0_dp) then
                p_low = p_mid
            else
                p_high = p_mid
            end if
        end do
        pressure = max(p_mid, this%p_floor)
    end function fan_pressure_left_tp


    function fan_pressure_right_tp(this, xi) result(pressure)
        class(riemann_solver), intent(in) :: this
        real(dp), intent(in) :: xi
        real(dp) :: pressure, p_low, p_high, p_mid, g_mid
        integer :: iter
        p_low = min(this%p_star, this%p_r)
        p_high = max(this%p_star, this%p_r)
        p_mid = 0.5_dp*(p_low + p_high)
        do iter = 1, 80
            p_mid = 0.5_dp*(p_low + p_high)
            g_mid = this%F_r_tp(p_mid) + this%sound_speed_tp(this%isentropic_temperature_tp(.false., p_mid), p_mid, this%Y_r) - xi
            if (abs(g_mid) <= 1.0e-8_dp*max(this%c_r,1.0_dp)) exit
            if (g_mid > 0.0_dp) then
                p_high = p_mid
            else
                p_low = p_mid
            end if
        end do
        pressure = max(p_mid, this%p_floor)
    end function fan_pressure_right_tp


    subroutine sample_thermally_perfect(this, xi, rho_out, u_out, p_out, T_out, gamma_out, Y_out)
        class(riemann_solver), intent(in) :: this
        real(dp), intent(in) :: xi
        real(dp), intent(out) :: rho_out, u_out, p_out, T_out, gamma_out
        real(dp), dimension(:), intent(out) :: Y_out
        real(dp) :: rho_star, T_star, c_star, mass_flux, wave_speed, pressure

        if (xi <= this%u_contact) then
            Y_out = this%Y_l
            if (this%p_star > this%p_l) then
                rho_star = this%shock_density_tp(.true., this%p_star)
                mass_flux = sqrt(max((this%p_star - this%p_l)/(1.0_dp/this%rho_l - 1.0_dp/rho_star), 0.0_dp))
                wave_speed = this%u_l - mass_flux/this%rho_l
                if (xi < wave_speed) then
                    rho_out = this%rho_l; u_out = this%u_l; p_out = this%p_l; T_out = this%T_l
                else
                    rho_out = rho_star; u_out = this%u_contact; p_out = this%p_star
                    T_out = this%temperature_from_p_rho_Y(p_out, rho_out, Y_out)
                end if
            else
                T_star = this%isentropic_temperature_tp(.true., this%p_star)
                c_star = this%sound_speed_tp(T_star, this%p_star, this%Y_l)
                if (xi < this%u_l - this%c_l) then
                    rho_out = this%rho_l; u_out = this%u_l; p_out = this%p_l; T_out = this%T_l
                elseif (xi > this%u_contact - c_star) then
                    p_out = this%p_star; u_out = this%u_contact; T_out = T_star
                    rho_out = p_out*this%M_l/(r_gase_J*T_out)
                else
                    pressure = this%fan_pressure_left_tp(xi)
                    p_out = pressure
                    T_out = this%isentropic_temperature_tp(.true., p_out)
                    u_out = this%F_l_tp(p_out)
                    rho_out = p_out*this%M_l/(r_gase_J*T_out)
                end if
            end if
        else
            Y_out = this%Y_r
            if (this%p_star > this%p_r) then
                rho_star = this%shock_density_tp(.false., this%p_star)
                mass_flux = sqrt(max((this%p_star - this%p_r)/(1.0_dp/this%rho_r - 1.0_dp/rho_star), 0.0_dp))
                wave_speed = this%u_r + mass_flux/this%rho_r
                if (xi > wave_speed) then
                    rho_out = this%rho_r; u_out = this%u_r; p_out = this%p_r; T_out = this%T_r
                else
                    rho_out = rho_star; u_out = this%u_contact; p_out = this%p_star
                    T_out = this%temperature_from_p_rho_Y(p_out, rho_out, Y_out)
                end if
            else
                T_star = this%isentropic_temperature_tp(.false., this%p_star)
                c_star = this%sound_speed_tp(T_star, this%p_star, this%Y_r)
                if (xi > this%u_r + this%c_r) then
                    rho_out = this%rho_r; u_out = this%u_r; p_out = this%p_r; T_out = this%T_r
                elseif (xi < this%u_contact + c_star) then
                    p_out = this%p_star; u_out = this%u_contact; T_out = T_star
                    rho_out = p_out*this%M_r/(r_gase_J*T_out)
                else
                    pressure = this%fan_pressure_right_tp(xi)
                    p_out = pressure
                    T_out = this%isentropic_temperature_tp(.false., p_out)
                    u_out = this%F_r_tp(p_out)
                    rho_out = p_out*this%M_r/(r_gase_J*T_out)
                end if
            end if
        end if
        rho_out = max(rho_out, this%rho_floor)
        p_out = max(p_out, this%p_floor)
        T_out = min(max(T_out, riemann_temperature_floor), riemann_temperature_ceiling)
        gamma_out = this%gamma_from_TY(T_out, Y_out)
    end subroutine sample_thermally_perfect


    subroutine sample_gamma(this, xi, rho_out, u_out, p_out, gamma_out)
        class(riemann_solver), intent(in) :: this
        real(dp), intent(in) :: xi
        real(dp), intent(out) :: rho_out, u_out, p_out, gamma_out
        real(dp) :: rho_star, c_side, c_star, mass_flux, wave_speed
        real(dp) :: gamma_side, p_side, rho_side, u_side, s_side, pfan, cfan

        if (xi <= this%u_contact) then
            gamma_side = this%gamma_l; p_side = this%p_l; rho_side = this%rho_l; u_side = this%u_l; s_side = this%s_l
            gamma_out = this%gamma_l
            if (this%p_star > this%p_l) then
                rho_star = this%rho_l/max(this%a_l(this%p_star), tiny(1.0_dp))
                mass_flux = sqrt(max((this%p_star - this%p_l)/(1.0_dp/this%rho_l - 1.0_dp/rho_star), 0.0_dp))
                wave_speed = this%u_l - mass_flux/this%rho_l
                if (xi < wave_speed) then
                    rho_out = this%rho_l; u_out = this%u_l; p_out = this%p_l
                else
                    rho_out = rho_star; u_out = this%u_contact; p_out = this%p_star
                end if
            else
                c_side = sqrt(this%gamma_l*this%p_l/this%rho_l)
                rho_star = (this%p_star/this%s_l)**(1.0_dp/this%gamma_l)
                c_star = sqrt(this%gamma_l*this%p_star/rho_star)
                if (xi < this%u_l - c_side) then
                    rho_out = this%rho_l; u_out = this%u_l; p_out = this%p_l
                elseif (xi > this%u_contact - c_star) then
                    rho_out = rho_star; u_out = this%u_contact; p_out = this%p_star
                else
                    u_out = 2.0_dp/(gamma_side + 1.0_dp)*(c_side + 0.5_dp*(gamma_side - 1.0_dp)*u_side + xi)
                    cfan = u_out - xi
                    pfan = p_side*(cfan/c_side)**(2.0_dp*gamma_side/(gamma_side - 1.0_dp))
                    p_out = max(pfan, this%p_floor)
                    rho_out = (p_out/s_side)**(1.0_dp/gamma_side)
                end if
            end if
        else
            gamma_side = this%gamma_r; p_side = this%p_r; rho_side = this%rho_r; u_side = this%u_r; s_side = this%s_r
            gamma_out = this%gamma_r
            if (this%p_star > this%p_r) then
                rho_star = this%rho_r/max(this%a_r(this%p_star), tiny(1.0_dp))
                mass_flux = sqrt(max((this%p_star - this%p_r)/(1.0_dp/this%rho_r - 1.0_dp/rho_star), 0.0_dp))
                wave_speed = this%u_r + mass_flux/this%rho_r
                if (xi > wave_speed) then
                    rho_out = this%rho_r; u_out = this%u_r; p_out = this%p_r
                else
                    rho_out = rho_star; u_out = this%u_contact; p_out = this%p_star
                end if
            else
                c_side = sqrt(this%gamma_r*this%p_r/this%rho_r)
                rho_star = (this%p_star/this%s_r)**(1.0_dp/this%gamma_r)
                c_star = sqrt(this%gamma_r*this%p_star/rho_star)
                if (xi > this%u_r + c_side) then
                    rho_out = this%rho_r; u_out = this%u_r; p_out = this%p_r
                elseif (xi < this%u_contact + c_star) then
                    rho_out = rho_star; u_out = this%u_contact; p_out = this%p_star
                else
                    u_out = 2.0_dp/(gamma_side + 1.0_dp)*(-c_side + 0.5_dp*(gamma_side - 1.0_dp)*u_side + xi)
                    cfan = xi - u_out
                    pfan = p_side*(cfan/c_side)**(2.0_dp*gamma_side/(gamma_side - 1.0_dp))
                    p_out = max(pfan, this%p_floor)
                    rho_out = (p_out/s_side)**(1.0_dp/gamma_side)
                end if
            end if
        end if
        rho_out = max(rho_out, this%rho_floor)
        p_out = max(p_out, this%p_floor)
    end subroutine sample_gamma


    subroutine sample_solution(this, xi, rho_out, u_out, p_out, T_out, gamma_out, Y_out, success)
        class(riemann_solver), intent(in) :: this
        real(dp), intent(in) :: xi
        real(dp), intent(out) :: rho_out, u_out, p_out, T_out, gamma_out
        real(dp), dimension(:), intent(out), optional :: Y_out
        logical, intent(out), optional :: success
        real(dp), dimension(:), allocatable :: Y_tmp
        logical :: ok

        ok = this%success
        rho_out = 0.0_dp; u_out = 0.0_dp; p_out = 0.0_dp; T_out = 0.0_dp; gamma_out = 0.0_dp
        if (.not. this%success) then
            if (present(success)) success = .false.
            if (present(Y_out)) Y_out = 0.0_dp
            return
        end if

        if (this%thermally_perfect_mode) then
            allocate(Y_tmp(this%species_number))
            call this%sample_thermally_perfect(xi, rho_out, u_out, p_out, T_out, gamma_out, Y_tmp)
            if (present(Y_out)) then
                Y_out = 0.0_dp
                Y_out(1:min(size(Y_out),size(Y_tmp))) = Y_tmp(1:min(size(Y_out),size(Y_tmp)))
            end if
            deallocate(Y_tmp)
        else
            call this%sample_gamma(xi, rho_out, u_out, p_out, gamma_out)
            T_out = 0.0_dp
            if (present(Y_out)) Y_out = 0.0_dp
        end if

        ok = this%finite_positive(rho_out) .and. this%finite_positive(p_out) .and. (u_out == u_out)
        if (this%thermally_perfect_mode) ok = ok .and. this%finite_positive(T_out)
        if (present(success)) success = ok
    end subroutine sample_solution


    subroutine write_shock_tube_reference(this, file_name, x0, time, x_min, x_max, n_points)
        class(riemann_solver), intent(in) :: this
        character(len=*), intent(in) :: file_name
        real(dp), intent(in) :: x0, time, x_min, x_max
        integer, intent(in) :: n_points
        integer :: unit_id, n, spec, ns
        real(dp) :: x, xi, dx, rho_s, u_s, p_s, T_s, gamma_s
        real(dp), dimension(:), allocatable :: Y_s
        logical :: ok

        ns = max(this%species_number, 0)
        allocate(Y_s(max(ns,1)))
        open(newunit=unit_id, file=trim(file_name), status='replace', form='formatted')
        write(unit_id,'(A)',advance='no') '# x xi rho u p T gamma'
        do spec = 1, ns
            write(unit_id,'(A,I0)',advance='no') ' Y', spec
        end do
        write(unit_id,*)

        if (n_points <= 1) then
            dx = 0.0_dp
        else
            dx = (x_max - x_min)/real(n_points - 1, dp)
        end if

        do n = 1, max(n_points,1)
            x = x_min + real(n - 1, dp)*dx
            if (time > 0.0_dp) then
                xi = (x - x0)/time
                call this%sample_solution(xi, rho_s, u_s, p_s, T_s, gamma_s, Y_s, ok)
            else
                xi = 0.0_dp
                if (x <= x0) then
                    rho_s = this%rho_l; u_s = this%u_l; p_s = this%p_l; gamma_s = this%gamma_l; T_s = this%T_l
                    if (ns > 0) Y_s(1:ns) = this%Y_l(1:ns)
                else
                    rho_s = this%rho_r; u_s = this%u_r; p_s = this%p_r; gamma_s = this%gamma_r; T_s = this%T_r
                    if (ns > 0) Y_s(1:ns) = this%Y_r(1:ns)
                end if
                ok = .true.
            end if
            if (ok) then
                write(unit_id,'(*(ES24.16,1X))') x, xi, rho_s, u_s, p_s, T_s, gamma_s, (Y_s(spec), spec=1,ns)
            end if
        end do
        close(unit_id)
        deallocate(Y_s)
    end subroutine write_shock_tube_reference


    function get_density(this) result(get_density_value)
        class(riemann_solver), intent(in) :: this
        real(dp) :: get_density_value
        get_density_value = this%rho
    end function get_density


    function get_pressure(this) result(get_pressure_value)
        class(riemann_solver), intent(in) :: this
        real(dp) :: get_pressure_value
        get_pressure_value = this%p
    end function get_pressure


    function get_gamma(this) result(get_gamma_value)
        class(riemann_solver), intent(in) :: this
        real(dp) :: get_gamma_value
        get_gamma_value = this%gamma
    end function get_gamma


    function get_velocity(this) result(get_velocity_value)
        class(riemann_solver), intent(in) :: this
        real(dp) :: get_velocity_value
        get_velocity_value = this%u
    end function get_velocity


    function get_temperature(this) result(get_temperature_value)
        class(riemann_solver), intent(in) :: this
        real(dp) :: get_temperature_value
        get_temperature_value = this%T
    end function get_temperature


    function get_contact_velocity(this) result(get_contact_velocity_value)
        class(riemann_solver), intent(in) :: this
        real(dp) :: get_contact_velocity_value
        get_contact_velocity_value = this%u_contact
    end function get_contact_velocity


    function get_star_pressure(this) result(get_star_pressure_value)
        class(riemann_solver), intent(in) :: this
        real(dp) :: get_star_pressure_value
        get_star_pressure_value = this%p_star
    end function get_star_pressure


    subroutine get_mass_fractions(this, Y_out)
        class(riemann_solver), intent(in) :: this
        real(dp), dimension(:), intent(out) :: Y_out
        integer :: ncopy
        Y_out = 0.0_dp
        if (allocated(this%Y)) then
            ncopy = min(size(Y_out), size(this%Y))
            if (ncopy > 0) Y_out(1:ncopy) = this%Y(1:ncopy)
        end if
    end subroutine get_mass_fractions


    function get_activation_count(this) result(get_activation_count_value)
        class(riemann_solver), intent(in) :: this
        integer :: get_activation_count_value
        get_activation_count_value = this%counter
    end function get_activation_count


    function get_failure_count(this) result(get_failure_count_value)
        class(riemann_solver), intent(in) :: this
        integer :: get_failure_count_value
        get_failure_count_value = this%failure_counter
    end function get_failure_count


    function get_success(this) result(get_success_value)
        class(riemann_solver), intent(in) :: this
        logical :: get_success_value
        get_success_value = this%success
    end function get_success

end module riemann_solver_class
