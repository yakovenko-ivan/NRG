module chemical_equilibrium_class

    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use kind_parameters, only: dp
    use global_data, only: r_gase_J, P_atm, tables_temperature_ceiling
    use chemical_properties_class, only: chemical_properties
    use thermophysical_properties_class, only: thermophysical_properties

    implicit none

    private
    public :: chemical_equilibrium, chemical_equilibrium_c
    public :: equilibrium_state, equilibrium_diagnostics

    real(dp), parameter :: minimum_equilibrium_temperature = 200.0_dp
    real(dp), parameter :: default_newton_tolerance = 1.0e-10_dp
    real(dp), parameter :: default_enthalpy_tolerance = 1.0e-9_dp
    integer, parameter :: default_newton_iterations = 100
    integer, parameter :: default_hp_iterations = 80

    type :: equilibrium_state
        real(dp) :: temperature = 0.0_dp
        real(dp) :: pressure = 0.0_dp
        real(dp) :: total_moles = 0.0_dp
        real(dp), allocatable :: mole_numbers(:)
        real(dp), allocatable :: mole_fractions(:)
        real(dp), allocatable :: mass_fractions(:)
    contains
        procedure :: clear => clear_equilibrium_state
    end type equilibrium_state

    type :: equilibrium_diagnostics
        logical :: converged = .false.
        integer :: iterations = 0
        real(dp) :: max_element_relative_residual = huge(1.0_dp)
        real(dp) :: mass_relative_residual = huge(1.0_dp)
        real(dp) :: mole_fraction_sum_residual = huge(1.0_dp)
        real(dp) :: stationarity_residual = huge(1.0_dp)
        real(dp) :: enthalpy_relative_residual = 0.0_dp
    end type equilibrium_diagnostics

    type :: chemical_equilibrium
        type(chemical_properties), pointer :: chemistry => null()
        type(thermophysical_properties), pointer :: thermophysics => null()
        integer :: elements_number = 0
        character(len=2), allocatable :: element_names(:)
        real(dp), allocatable :: species_element_counts(:,:)
    contains
        procedure :: equilibrate_tp
        procedure :: equilibrate_hp
        procedure :: get_element_index
        procedure :: get_species_element_count
        procedure :: element_totals
        procedure :: clear => clear_chemical_equilibrium
    end type chemical_equilibrium

    interface chemical_equilibrium_c
        module procedure constructor
    end interface chemical_equilibrium_c

contains

    type(chemical_equilibrium) function constructor(chemistry,thermophysics)
        type(chemical_properties), target, intent(in) :: chemistry
        type(thermophysical_properties), target, intent(in) :: thermophysics

        constructor%chemistry => chemistry
        constructor%thermophysics => thermophysics
        if (.not. allocated(thermophysics%element_names) .or. &
            .not. allocated(thermophysics%species_element_counts)) then
            error stop 'Chemical equilibrium: thermophysics has no elemental composition'
        end if
        constructor%elements_number = thermophysics%elements_number
        allocate(constructor%element_names(constructor%elements_number))
        allocate(constructor%species_element_counts(constructor%elements_number, &
            chemistry%species_number))
        constructor%element_names = thermophysics%element_names
        constructor%species_element_counts = thermophysics%species_element_counts
    end function constructor


    subroutine clear_chemical_equilibrium(this)
        class(chemical_equilibrium), intent(inout) :: this

        nullify(this%chemistry)
        nullify(this%thermophysics)
        this%elements_number = 0
        if (allocated(this%element_names)) deallocate(this%element_names)
        if (allocated(this%species_element_counts)) &
            deallocate(this%species_element_counts)
    end subroutine clear_chemical_equilibrium


    subroutine clear_equilibrium_state(this)
        class(equilibrium_state), intent(inout) :: this

        this%temperature = 0.0_dp
        this%pressure = 0.0_dp
        this%total_moles = 0.0_dp
        if (allocated(this%mole_numbers)) deallocate(this%mole_numbers)
        if (allocated(this%mole_fractions)) deallocate(this%mole_fractions)
        if (allocated(this%mass_fractions)) deallocate(this%mass_fractions)
    end subroutine clear_equilibrium_state


    integer function get_element_index(this,name) result(index_value)
        class(chemical_equilibrium), intent(in) :: this
        character(len=*), intent(in) :: name
        integer :: element
        character(len=2) :: requested

        index_value = 0
        requested = to_upper_ascii(adjustl(name))
        do element = 1, this%elements_number
            if (trim(requested) == trim(to_upper_ascii(this%element_names(element)))) then
                index_value = element
                return
            end if
        end do
    end function get_element_index


    real(dp) function get_species_element_count(this,element,specie) result(count)
        class(chemical_equilibrium), intent(in) :: this
        integer, intent(in) :: element,specie

        if (element < 1 .or. element > this%elements_number) then
            error stop 'Chemical equilibrium: invalid element index'
        end if
        if (specie < 1 .or. specie > this%chemistry%species_number) then
            error stop 'Chemical equilibrium: invalid species index'
        end if
        count = this%species_element_counts(element,specie)
    end function get_species_element_count


    subroutine element_totals(this,mole_numbers,totals)
        class(chemical_equilibrium), intent(in) :: this
        real(dp), dimension(:), intent(in) :: mole_numbers
        real(dp), dimension(:), intent(out) :: totals

        if (size(mole_numbers) /= this%chemistry%species_number) then
            error stop 'Chemical equilibrium: mole-number size mismatch'
        end if
        if (size(totals) /= this%elements_number) then
            error stop 'Chemical equilibrium: element-total size mismatch'
        end if
        totals = matmul(this%species_element_counts,mole_numbers)
    end subroutine element_totals


    subroutine equilibrate_tp(this,temperature,pressure,mole_numbers_initial, &
            state,diagnostics,mole_numbers_guess)
        class(chemical_equilibrium), intent(in) :: this
        real(dp), intent(in) :: temperature,pressure
        real(dp), dimension(:), intent(in) :: mole_numbers_initial
        type(equilibrium_state), intent(inout) :: state
        type(equilibrium_diagnostics), intent(out) :: diagnostics
        real(dp), dimension(:), intent(in), optional :: mole_numbers_guess

        integer :: species_number,active_count,element,specie,iter,line_iter
        integer :: variable_count
        integer, allocatable :: active_elements(:)
        logical, allocatable :: allowed_species(:)
        real(dp), allocatable :: element_inventory(:),lambda(:),variables(:)
        real(dp), allocatable :: residual(:),trial_residual(:),jacobian(:,:)
        real(dp), allocatable :: delta(:),trial_variables(:),log_q(:),q_scaled(:)
        real(dp), allocatable :: gibbs_rt(:),element_sum_scaled(:)
        real(dp) :: total_initial,log_scale,residual_norm,trial_norm,step
        real(dp) :: active_threshold
        real(dp) :: total_moles,mass_initial,mass_final,element_final
        real(dp) :: denominator,mu_residual,pressure_ratio
        logical :: solved

        diagnostics = equilibrium_diagnostics()
        species_number = this%chemistry%species_number

        if (.not. associated(this%chemistry) .or. &
            .not. associated(this%thermophysics)) then
            error stop 'Chemical equilibrium: uninitialized calculator'
        end if
        if (size(mole_numbers_initial) /= species_number) then
            error stop 'Chemical equilibrium: mole-number size mismatch'
        end if
        if (.not. ieee_is_finite(temperature) .or. temperature <= 0.0_dp) then
            error stop 'Chemical equilibrium: invalid temperature'
        end if
        if (.not. ieee_is_finite(pressure) .or. pressure <= 0.0_dp) then
            error stop 'Chemical equilibrium: invalid pressure'
        end if
        if (temperature > tables_temperature_ceiling) then
            error stop 'Chemical equilibrium: temperature exceeds thermo ceiling'
        end if
        if (any(mole_numbers_initial < 0.0_dp)) then
            error stop 'Chemical equilibrium: negative initial mole number'
        end if
        if (present(mole_numbers_guess)) then
            if (size(mole_numbers_guess) /= species_number) then
                error stop 'Chemical equilibrium: mole-number guess size mismatch'
            end if
            if (any(mole_numbers_guess < 0.0_dp)) then
                error stop 'Chemical equilibrium: negative mole-number guess'
            end if
            if (sum(mole_numbers_guess) <= tiny(1.0_dp)) then
                error stop 'Chemical equilibrium: empty mole-number guess'
            end if
        end if

        total_initial = sum(mole_numbers_initial)
        if (total_initial <= tiny(1.0_dp)) then
            error stop 'Chemical equilibrium: empty initial composition'
        end if

        allocate(element_inventory(this%elements_number))
        call this%element_totals(mole_numbers_initial,element_inventory)
        active_threshold = 100.0_dp*epsilon(1.0_dp)* &
            max(maxval(element_inventory),tiny(1.0_dp))
        active_count = count(element_inventory > active_threshold)
        if (active_count <= 0) then
            error stop 'Chemical equilibrium: no active elements'
        end if

        allocate(active_elements(active_count))
        active_count = 0
        do element = 1, this%elements_number
            if (element_inventory(element) > active_threshold) then
                active_count = active_count + 1
                active_elements(active_count) = element
            end if
        end do

        allocate(allowed_species(species_number))
        allowed_species = .true.
        do specie = 1, species_number
            do element = 1, this%elements_number
                if (element_inventory(element) > active_threshold) cycle
                if (this%species_element_counts(element,specie) > 0.0_dp) then
                    allowed_species(specie) = .false.
                    exit
                end if
            end do
        end do
        if (.not. any(allowed_species)) then
            error stop 'Chemical equilibrium: no species compatible with elemental inventory'
        end if

        variable_count = active_count+1
        allocate(lambda(active_count),variables(variable_count))
        allocate(residual(variable_count),trial_residual(variable_count))
        allocate(jacobian(variable_count,variable_count),delta(variable_count))
        allocate(trial_variables(variable_count))
        allocate(log_q(species_number),q_scaled(species_number))
        allocate(gibbs_rt(species_number),element_sum_scaled(active_count))

        pressure_ratio = pressure/P_atm
        do specie = 1, species_number
            gibbs_rt(specie) = ( &
                this%thermophysics%specie_enthalpy_molar(temperature,specie) - &
                temperature*this%thermophysics%specie_entropy_molar( &
                    temperature,specie))/(r_gase_J*temperature)
        end do

        if (present(mole_numbers_guess)) then
            call initialize_element_potentials(this,active_elements, &
                mole_numbers_guess,allowed_species,pressure,gibbs_rt,lambda)
        else
            call initialize_element_potentials(this,active_elements, &
                mole_numbers_initial,allowed_species,pressure,gibbs_rt,lambda)
        end if
        variables(1:active_count) = lambda
        variables(variable_count) = log(total_initial)

        solved = .false.
        do iter = 1, default_newton_iterations
            call equilibrium_residual_and_jacobian(this,active_elements, &
                allowed_species,element_inventory,gibbs_rt,pressure_ratio,variables, &
                residual,jacobian,log_q,q_scaled,element_sum_scaled)
            residual_norm = maxval(abs(residual))
            diagnostics%iterations = iter
            if (residual_norm <= default_newton_tolerance) then
                solved = .true.
                exit
            end if

            delta = -residual
            call solve_dense_linear_system(jacobian,delta,solved)
            if (.not. solved) exit

            step = 1.0_dp
            solved = .false.
            do line_iter = 1, 24
                trial_variables = variables + step*delta
                call equilibrium_residual_only(this,active_elements, &
                    allowed_species,element_inventory,gibbs_rt,pressure_ratio,trial_variables, &
                    trial_residual)
                trial_norm = maxval(abs(trial_residual))
                if (ieee_is_finite(trial_norm) .and. &
                    trial_norm < residual_norm) then
                    variables = trial_variables
                    solved = .true.
                    exit
                end if
                step = 0.5_dp*step
            end do
            if (.not. solved) exit
            solved = .false.
        end do

        ! Re-evaluate at the accepted state and construct mole numbers from the
        ! same log-abundance representation used by Newton.
        call equilibrium_residual_and_jacobian(this,active_elements, &
            allowed_species,element_inventory,gibbs_rt,pressure_ratio,variables, &
            residual,jacobian,log_q,q_scaled,element_sum_scaled)
        residual_norm = maxval(abs(residual))
        diagnostics%converged = residual_norm <= default_newton_tolerance

        log_scale = maxval(log_q)
        q_scaled = 0.0_dp
        where (allowed_species) q_scaled = exp(log_q-log_scale)
        denominator = sum(q_scaled)
        total_moles = exp(variables(variable_count))

        call state%clear()
        allocate(state%mole_numbers(species_number))
        allocate(state%mole_fractions(species_number))
        allocate(state%mass_fractions(species_number))
        state%temperature = temperature
        state%pressure = pressure
        state%total_moles = total_moles
        state%mole_fractions = q_scaled/denominator
        state%mole_numbers = total_moles*state%mole_fractions

        mass_final = sum(state%mole_numbers*this%thermophysics%molar_masses)
        if (mass_final <= tiny(1.0_dp)) then
            error stop 'Chemical equilibrium: non-positive equilibrium mass'
        end if
        state%mass_fractions = state%mole_numbers* &
            this%thermophysics%molar_masses/mass_final

        diagnostics%mole_fraction_sum_residual = abs(sum(state%mole_fractions)-1.0_dp)
        diagnostics%max_element_relative_residual = 0.0_dp
        do element = 1, active_count
            specie = active_elements(element)
            element_final = dot_product( &
                this%species_element_counts(specie,:),state%mole_numbers)
            diagnostics%max_element_relative_residual = max( &
                diagnostics%max_element_relative_residual, &
                abs(element_final-element_inventory(specie))/ &
                max(abs(element_inventory(specie)),tiny(1.0_dp)))
        end do

        mass_initial = sum(mole_numbers_initial*this%thermophysics%molar_masses)
        diagnostics%mass_relative_residual = abs(mass_final-mass_initial)/ &
            max(abs(mass_initial),tiny(1.0_dp))

        diagnostics%stationarity_residual = 0.0_dp
        do specie = 1, species_number
            if (.not. allowed_species(specie)) cycle
            if (state%mole_fractions(specie) <= 1.0e-30_dp) cycle
            mu_residual = gibbs_rt(specie) + &
                log(state%mole_fractions(specie)*pressure_ratio)
            do element = 1, active_count
                mu_residual = mu_residual + variables(element)* &
                    this%species_element_counts(active_elements(element),specie)
            end do
            diagnostics%stationarity_residual = max( &
                diagnostics%stationarity_residual,abs(mu_residual))
        end do

        diagnostics%converged = diagnostics%converged .and. &
            diagnostics%max_element_relative_residual <= 1.0e-8_dp .and. &
            diagnostics%mass_relative_residual <= 1.0e-8_dp

        deallocate(element_inventory,active_elements,allowed_species,lambda,variables,residual, &
            trial_residual,jacobian,delta,trial_variables,log_q,q_scaled, &
            gibbs_rt,element_sum_scaled)
    end subroutine equilibrate_tp


    subroutine equilibrate_hp(this,initial_temperature,pressure, &
            mole_numbers_initial,state,diagnostics)
        class(chemical_equilibrium), intent(in) :: this
        real(dp), intent(in) :: initial_temperature,pressure
        real(dp), dimension(:), intent(in) :: mole_numbers_initial
        type(equilibrium_state), intent(inout) :: state
        type(equilibrium_diagnostics), intent(out) :: diagnostics

        type(equilibrium_state) :: trial_state,best_state
        type(equilibrium_state) :: previous_state,left_state,right_state
        type(equilibrium_diagnostics) :: trial_diag,best_diag
        integer, parameter :: scan_points = 80
        integer :: specie,point,iter
        real(dp) :: initial_enthalpy,t_low,t_high,t_left,t_right,t_trial
        real(dp) :: f_left,f_right,f_trial,best_abs_residual,enthalpy_scale
        real(dp) :: scan_t,scan_f,previous_t,previous_f
        logical :: bracket_found,best_found,previous_converged
        logical :: have_continuation_seed

        if (initial_temperature <= 0.0_dp .or. &
            initial_temperature > tables_temperature_ceiling) then
            error stop 'Chemical equilibrium HP: invalid initial temperature'
        end if
        if (size(mole_numbers_initial) /= this%chemistry%species_number) then
            error stop 'Chemical equilibrium HP: mole-number size mismatch'
        end if

        initial_enthalpy = 0.0_dp
        do specie = 1, this%chemistry%species_number
            initial_enthalpy = initial_enthalpy + mole_numbers_initial(specie)* &
                this%thermophysics%specie_enthalpy_molar( &
                    initial_temperature,specie)
        end do
        enthalpy_scale = max(abs(initial_enthalpy),r_gase_J* &
            max(initial_temperature,300.0_dp)*max(sum(mole_numbers_initial),1.0_dp))

        t_low = minimum_equilibrium_temperature
        t_high = tables_temperature_ceiling
        bracket_found = .false.
        best_abs_residual = huge(1.0_dp)
        best_found = .false.
        previous_converged = .false.
        have_continuation_seed = .false.

        ! Search from high to low temperature. Once a TP state converges, use
        ! that equilibrium composition only as the initializer for the next
        ! lower-temperature TP solve. The elemental inventory remains the one
        ! defined by mole_numbers_initial. This continuation is especially
        ! important for cold, very lean mixtures where direct initialization
        ! from the unburned reactants can be far from the equilibrium solution.
        do point = 1, scan_points
            scan_t = t_high - (t_high-t_low)*real(point-1,dp)/ &
                real(scan_points-1,dp)

            if (have_continuation_seed) then
                call this%equilibrate_tp(scan_t,pressure,mole_numbers_initial, &
                    trial_state,trial_diag,previous_state%mole_numbers)
                if (.not. trial_diag%converged) then
                    ! Retain the original unseeded path as a fallback.
                    call this%equilibrate_tp(scan_t,pressure,mole_numbers_initial, &
                        trial_state,trial_diag)
                end if
            else
                call this%equilibrate_tp(scan_t,pressure,mole_numbers_initial, &
                    trial_state,trial_diag)
            end if

            if (.not. trial_diag%converged) cycle

            scan_f = total_enthalpy(this,trial_state)-initial_enthalpy
            if (abs(scan_f) < best_abs_residual) then
                best_abs_residual = abs(scan_f)
                best_state = trial_state
                best_diag = trial_diag
                best_found = .true.
            end if

            if (previous_converged .and. &
                opposite_or_zero_sign(previous_f,scan_f)) then
                ! Scanning downward means scan_t is the lower-temperature
                ! endpoint and previous_t is the higher-temperature endpoint.
                t_left = scan_t
                f_left = scan_f
                left_state = trial_state
                t_right = previous_t
                f_right = previous_f
                right_state = previous_state
                bracket_found = .true.
                exit
            end if

            previous_t = scan_t
            previous_f = scan_f
            previous_state = trial_state
            previous_converged = .true.
            have_continuation_seed = .true.
        end do

        if (.not. bracket_found) then
            if (.not. best_found) then
                call state%clear()
                state%temperature = initial_temperature
                state%pressure = pressure
                diagnostics = equilibrium_diagnostics()
                diagnostics%converged = .false.
                return
            end if
            state = best_state
            diagnostics = best_diag
            diagnostics%converged = .false.
            diagnostics%enthalpy_relative_residual = &
                best_abs_residual/enthalpy_scale
            return
        end if

        do iter = 1, default_hp_iterations
            if (abs(f_right-f_left) > tiny(1.0_dp)) then
                t_trial = t_right-f_right*(t_right-t_left)/(f_right-f_left)
            else
                t_trial = 0.5_dp*(t_left+t_right)
            end if
            if (t_trial <= t_left .or. t_trial >= t_right .or. &
                .not. ieee_is_finite(t_trial)) then
                t_trial = 0.5_dp*(t_left+t_right)
            end if

            ! Warm-start from the closer bracket endpoint. If that seeded solve
            ! fails, try the opposite endpoint before falling back to the
            ! original unseeded initialization.
            if (abs(t_trial-t_left) <= abs(t_right-t_trial)) then
                call this%equilibrate_tp(t_trial,pressure,mole_numbers_initial, &
                    trial_state,trial_diag,left_state%mole_numbers)
                if (.not. trial_diag%converged) then
                    call this%equilibrate_tp(t_trial,pressure,mole_numbers_initial, &
                        trial_state,trial_diag,right_state%mole_numbers)
                end if
            else
                call this%equilibrate_tp(t_trial,pressure,mole_numbers_initial, &
                    trial_state,trial_diag,right_state%mole_numbers)
                if (.not. trial_diag%converged) then
                    call this%equilibrate_tp(t_trial,pressure,mole_numbers_initial, &
                        trial_state,trial_diag,left_state%mole_numbers)
                end if
            end if
            if (.not. trial_diag%converged) then
                call this%equilibrate_tp(t_trial,pressure,mole_numbers_initial, &
                    trial_state,trial_diag)
            end if
            if (.not. trial_diag%converged) then
                ! A failed TP evaluation cannot safely update the bracket.
                ! Use a strict bisection point and retry both continuation seeds.
                t_trial = 0.5_dp*(t_left+t_right)
                call this%equilibrate_tp(t_trial,pressure,mole_numbers_initial, &
                    trial_state,trial_diag,left_state%mole_numbers)
                if (.not. trial_diag%converged) then
                    call this%equilibrate_tp(t_trial,pressure,mole_numbers_initial, &
                        trial_state,trial_diag,right_state%mole_numbers)
                end if
            end if
            if (.not. trial_diag%converged) exit

            f_trial = total_enthalpy(this,trial_state)-initial_enthalpy
            if (abs(f_trial)/enthalpy_scale < default_enthalpy_tolerance) exit

            if (opposite_or_zero_sign(f_left,f_trial)) then
                t_right = t_trial
                f_right = f_trial
                right_state = trial_state
            else
                t_left = t_trial
                f_left = f_trial
                left_state = trial_state
            end if
        end do

        if (.not. trial_diag%converged) then
            ! Preserve the best physically valid TP state for diagnostics.
            state = best_state
            diagnostics = best_diag
            diagnostics%converged = .false.
            diagnostics%enthalpy_relative_residual = &
                best_abs_residual/enthalpy_scale
            return
        end if

        state = trial_state
        diagnostics = trial_diag
        diagnostics%iterations = diagnostics%iterations + iter
        diagnostics%enthalpy_relative_residual = abs(f_trial)/enthalpy_scale
        diagnostics%converged = diagnostics%converged .and. &
            diagnostics%enthalpy_relative_residual < default_enthalpy_tolerance
    end subroutine equilibrate_hp


    real(dp) function total_enthalpy(this,state) result(enthalpy)
        class(chemical_equilibrium), intent(in) :: this
        type(equilibrium_state), intent(in) :: state
        integer :: specie

        enthalpy = 0.0_dp
        do specie = 1, this%chemistry%species_number
            enthalpy = enthalpy + state%mole_numbers(specie)* &
                this%thermophysics%specie_enthalpy_molar( &
                    state%temperature,specie)
        end do
    end function total_enthalpy


    subroutine equilibrium_residual_only(this,active_elements,allowed_species, &
            element_inventory,gibbs_rt,pressure_ratio,variables,residual)
        class(chemical_equilibrium), intent(in) :: this
        integer, dimension(:), intent(in) :: active_elements
        logical, dimension(:), intent(in) :: allowed_species
        real(dp), dimension(:), intent(in) :: element_inventory,gibbs_rt,variables
        real(dp), intent(in) :: pressure_ratio
        real(dp), dimension(:), intent(out) :: residual

        real(dp), allocatable :: dummy_jacobian(:,:),log_q(:),q_scaled(:)
        real(dp), allocatable :: element_sum_scaled(:)
        integer :: nvar,nspecies

        nvar = size(variables)
        nspecies = size(gibbs_rt)
        allocate(dummy_jacobian(nvar,nvar),log_q(nspecies),q_scaled(nspecies))
        allocate(element_sum_scaled(size(active_elements)))
        call equilibrium_residual_and_jacobian(this,active_elements,allowed_species, &
            element_inventory,gibbs_rt,pressure_ratio,variables,residual, &
            dummy_jacobian,log_q,q_scaled,element_sum_scaled)
        deallocate(dummy_jacobian,log_q,q_scaled,element_sum_scaled)
    end subroutine equilibrium_residual_only


    subroutine equilibrium_residual_and_jacobian(this,active_elements,allowed_species, &
            element_inventory,gibbs_rt,pressure_ratio,variables,residual, &
            jacobian,log_q,q_scaled,element_sum_scaled)
        class(chemical_equilibrium), intent(in) :: this
        integer, dimension(:), intent(in) :: active_elements
        logical, dimension(:), intent(in) :: allowed_species
        real(dp), dimension(:), intent(in) :: element_inventory,gibbs_rt,variables
        real(dp), intent(in) :: pressure_ratio
        real(dp), dimension(:), intent(out) :: residual
        real(dp), dimension(:,:), intent(out) :: jacobian
        real(dp), dimension(:), intent(out) :: log_q,q_scaled
        real(dp), dimension(:), intent(out) :: element_sum_scaled

        integer :: active_count,species_number,e,j,specie
        real(dp) :: log_scale,q_sum_scaled,cross_sum

        active_count = size(active_elements)
        species_number = size(gibbs_rt)

        do specie = 1, species_number
            if (.not. allowed_species(specie)) then
                log_q(specie) = -huge(1.0_dp)
                cycle
            end if
            log_q(specie) = -gibbs_rt(specie)-log(pressure_ratio)
            do e = 1, active_count
                log_q(specie) = log_q(specie)-variables(e)* &
                    this%species_element_counts(active_elements(e),specie)
            end do
        end do

        log_scale = maxval(log_q,mask=allowed_species)
        q_scaled = 0.0_dp
        where (allowed_species) q_scaled = exp(log_q-log_scale)
        q_sum_scaled = sum(q_scaled)
        residual(1) = log_scale+log(q_sum_scaled)

        do e = 1, active_count
            element_sum_scaled(e) = dot_product( &
                this%species_element_counts(active_elements(e),:),q_scaled)
            if (element_sum_scaled(e) <= tiny(1.0_dp)) then
                residual(e+1) = huge(1.0_dp)
            else
                residual(e+1) = variables(active_count+1)+log_scale+ &
                    log(element_sum_scaled(e))- &
                    log(element_inventory(active_elements(e)))
            end if
        end do

        jacobian = 0.0_dp
        do j = 1, active_count
            jacobian(1,j) = -dot_product(q_scaled, &
                this%species_element_counts(active_elements(j),:))/ &
                q_sum_scaled
        end do
        do e = 1, active_count
            jacobian(e+1,active_count+1) = 1.0_dp
            if (element_sum_scaled(e) <= tiny(1.0_dp)) cycle
            do j = 1, active_count
                cross_sum = 0.0_dp
                do specie = 1, species_number
                    cross_sum = cross_sum + &
                        this%species_element_counts(active_elements(e),specie)* &
                        q_scaled(specie)* &
                        this%species_element_counts(active_elements(j),specie)
                end do
                jacobian(e+1,j) = -cross_sum/element_sum_scaled(e)
            end do
        end do
    end subroutine equilibrium_residual_and_jacobian


    subroutine initialize_element_potentials(this,active_elements, &
            mole_numbers_initial,allowed_species,pressure,gibbs_rt,lambda)
        class(chemical_equilibrium), intent(in) :: this
        integer, dimension(:), intent(in) :: active_elements
        logical, dimension(:), intent(in) :: allowed_species
        real(dp), dimension(:), intent(in) :: mole_numbers_initial,gibbs_rt
        real(dp), intent(in) :: pressure
        real(dp), dimension(:), intent(out) :: lambda

        integer :: active_count,species_number,e,j,specie
        real(dp), allocatable :: normal_matrix(:,:),rhs(:),x_initial(:)
        real(dp) :: total_moles,target,weight
        logical :: solved

        active_count = size(active_elements)
        species_number = size(mole_numbers_initial)
        allocate(normal_matrix(active_count,active_count),rhs(active_count))
        allocate(x_initial(species_number))
        normal_matrix = 0.0_dp
        rhs = 0.0_dp
        total_moles = sum(mole_numbers_initial)
        x_initial = mole_numbers_initial/total_moles

        do specie = 1, species_number
            if (.not. allowed_species(specie)) cycle
            if (x_initial(specie) <= 1.0e-14_dp) cycle
            target = -gibbs_rt(specie)-log(x_initial(specie)*pressure/P_atm)
            weight = max(x_initial(specie),1.0e-6_dp)
            do e = 1, active_count
                rhs(e) = rhs(e)+weight* &
                    this%species_element_counts(active_elements(e),specie)*target
                do j = 1, active_count
                    normal_matrix(e,j) = normal_matrix(e,j)+weight* &
                        this%species_element_counts(active_elements(e),specie)* &
                        this%species_element_counts(active_elements(j),specie)
                end do
            end do
        end do

        do e = 1, active_count
            normal_matrix(e,e) = normal_matrix(e,e) + 1.0e-12_dp
        end do
        lambda = rhs
        call solve_dense_linear_system(normal_matrix,lambda,solved)
        if (.not. solved .or. any(.not. ieee_is_finite(lambda))) lambda = 0.0_dp

        deallocate(normal_matrix,rhs,x_initial)
    end subroutine initialize_element_potentials


    subroutine solve_dense_linear_system(matrix,rhs,solved)
        real(dp), dimension(:,:), intent(in) :: matrix
        real(dp), dimension(:), intent(inout) :: rhs
        logical, intent(out) :: solved

        real(dp), allocatable :: a(:,:)
        real(dp) :: pivot_value,factor,temp
        integer :: n,i,j,k,pivot

        n = size(rhs)
        if (size(matrix,1) /= n .or. size(matrix,2) /= n) then
            error stop 'Chemical equilibrium: dense solve size mismatch'
        end if
        allocate(a(n,n))
        a = matrix
        solved = .false.

        do k = 1, n-1
            pivot = k-1+maxloc(abs(a(k:n,k)),dim=1)
            pivot_value = abs(a(pivot,k))
            if (pivot_value <= 100.0_dp*epsilon(1.0_dp)* &
                max(1.0_dp,maxval(abs(a)))) then
                deallocate(a)
                return
            end if
            if (pivot /= k) then
                do j = k, n
                    temp = a(k,j)
                    a(k,j) = a(pivot,j)
                    a(pivot,j) = temp
                end do
                temp = rhs(k)
                rhs(k) = rhs(pivot)
                rhs(pivot) = temp
            end if
            do i = k+1, n
                factor = a(i,k)/a(k,k)
                a(i,k) = 0.0_dp
                do j = k+1, n
                    a(i,j) = a(i,j)-factor*a(k,j)
                end do
                rhs(i) = rhs(i)-factor*rhs(k)
            end do
        end do

        if (abs(a(n,n)) <= 100.0_dp*epsilon(1.0_dp)* &
            max(1.0_dp,maxval(abs(a)))) then
            deallocate(a)
            return
        end if

        do i = n, 1, -1
            if (i < n) rhs(i) = rhs(i)-dot_product(a(i,i+1:n),rhs(i+1:n))
            rhs(i) = rhs(i)/a(i,i)
        end do
        solved = all(ieee_is_finite(rhs))
        deallocate(a)
    end subroutine solve_dense_linear_system


    pure function to_upper_ascii(value) result(upper_value)
        character(len=*), intent(in) :: value
        character(len=len(value)) :: upper_value
        integer :: i,code

        upper_value = value
        do i = 1, len(value)
            code = iachar(upper_value(i:i))
            if (code >= iachar('a') .and. code <= iachar('z')) then
                upper_value(i:i) = achar(code-iachar('a')+iachar('A'))
            end if
        end do
    end function to_upper_ascii


    pure logical function opposite_or_zero_sign(a,b) result(brackets)
        real(dp), intent(in) :: a,b

        brackets = (a <= 0.0_dp .and. b >= 0.0_dp) .or. &
            (a >= 0.0_dp .and. b <= 0.0_dp)
    end function opposite_or_zero_sign

end module chemical_equilibrium_class
