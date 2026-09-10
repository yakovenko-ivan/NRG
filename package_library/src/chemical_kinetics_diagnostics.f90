module chemical_kinetics_diagnostics

    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite

    use kind_parameters, only: dp
    use global_data, only: r_gase_J
    use chemical_kinetics_core_class, only: chemical_kinetics_core, &
        chemical_rate_state, reaction_rate_diagnostics

    implicit none

    private
    public :: hydrogen_crossover_diagnostics
    public :: hydrogen_crossover_temperature

    type :: hydrogen_crossover_diagnostics
        logical :: converged = .false.
        integer :: iterations = 0
        integer, allocatable :: branching_reactions(:)
        integer, allocatable :: stabilization_reactions(:)
        real(dp) :: branching_rate_coefficient = 0.0_dp
        real(dp) :: stabilization_rate_coefficient = 0.0_dp
        real(dp) :: relative_balance_residual = huge(1.0_dp)
        type(reaction_rate_diagnostics), allocatable :: branching_details(:)
        type(reaction_rate_diagnostics), allocatable :: stabilization_details(:)
    end type hydrogen_crossover_diagnostics

contains

    subroutine hydrogen_crossover_temperature(core,pressure,mole_fractions, &
            crossover_temperature,diagnostics,temperature_min,temperature_max)
        type(chemical_kinetics_core), intent(in) :: core
        real(dp), intent(in) :: pressure
        real(dp), dimension(:), intent(in) :: mole_fractions
        real(dp), intent(out) :: crossover_temperature
        type(hydrogen_crossover_diagnostics), intent(out) :: diagnostics
        real(dp), intent(in), optional :: temperature_min,temperature_max

        integer, parameter :: scan_points = 80
        integer, parameter :: max_root_iterations = 80
        real(dp), parameter :: root_tolerance = 1.0e-10_dp

        integer :: h_index,o2_index,o_index,oh_index,ho2_index
        integer :: reaction,branch_count,stabilization_count,point,iter
        integer, allocatable :: branch_temp(:),stabilization_temp(:)
        integer, dimension(2) :: reactants,reactant_nu,branch_products
        integer, dimension(2) :: branch_product_nu
        integer, dimension(1) :: stabilization_products,stabilization_product_nu
        real(dp), allocatable :: x(:)
        real(dp) :: t_min,t_max,t_previous,t_current,t_left,t_right,t_trial
        real(dp) :: f_previous,f_current,f_left,f_right,f_trial
        real(dp) :: k_branch,k_stabilization,best_t,best_relative,current_relative
        logical :: bracket_found
        type(chemical_rate_state) :: rate_state

        diagnostics = hydrogen_crossover_diagnostics()
        crossover_temperature = 0.0_dp

        if (.not. associated(core%chemistry)) then
            error stop 'Hydrogen crossover: uninitialized kinetics core'
        end if
        if (size(mole_fractions) /= core%species_number) then
            error stop 'Hydrogen crossover: mole-fraction size mismatch'
        end if
        if (.not. ieee_is_finite(pressure) .or. pressure <= 0.0_dp) then
            error stop 'Hydrogen crossover: invalid pressure'
        end if
        if (any(mole_fractions < 0.0_dp) .or. &
            sum(mole_fractions) <= tiny(1.0_dp)) then
            error stop 'Hydrogen crossover: invalid bath composition'
        end if

        allocate(x(core%species_number))
        x = mole_fractions/sum(mole_fractions)

        h_index = species_index(core,'H')
        o2_index = species_index(core,'O2')
        o_index = species_index(core,'O')
        oh_index = species_index(core,'OH')
        ho2_index = species_index(core,'HO2')
        if (min(h_index,o2_index,o_index,oh_index,ho2_index) <= 0) then
            error stop 'Hydrogen crossover: H/O2/O/OH/HO2 species are required'
        end if

        reactants = (/h_index,o2_index/)
        reactant_nu = (/1,1/)
        branch_products = (/o_index,oh_index/)
        branch_product_nu = (/1,1/)
        stabilization_products = (/ho2_index/)
        stabilization_product_nu = (/1/)

        allocate(branch_temp(core%reactions_number))
        allocate(stabilization_temp(core%reactions_number))
        branch_count = 0
        stabilization_count = 0
        do reaction = 1, core%reactions_number
            if (core%reaction_matches(reaction,reactants,reactant_nu, &
                branch_products,branch_product_nu)) then
                branch_count = branch_count+1
                branch_temp(branch_count) = reaction
            end if
            if (core%reaction_matches(reaction,reactants,reactant_nu, &
                stabilization_products,stabilization_product_nu)) then
                stabilization_count = stabilization_count+1
                stabilization_temp(stabilization_count) = reaction
            end if
        end do
        if (branch_count <= 0) then
            error stop 'Hydrogen crossover: H+O2 -> O+OH reaction not found'
        end if
        if (stabilization_count <= 0) then
            error stop 'Hydrogen crossover: H+O2 -> HO2 reaction not found'
        end if

        allocate(diagnostics%branching_reactions(branch_count))
        allocate(diagnostics%stabilization_reactions(stabilization_count))
        diagnostics%branching_reactions = branch_temp(1:branch_count)
        diagnostics%stabilization_reactions = &
            stabilization_temp(1:stabilization_count)
        deallocate(branch_temp,stabilization_temp)

        t_min = 400.0_dp
        t_max = 3000.0_dp
        if (present(temperature_min)) t_min = temperature_min
        if (present(temperature_max)) t_max = temperature_max
        if (t_min <= 0.0_dp .or. t_max <= t_min) then
            error stop 'Hydrogen crossover: invalid temperature bracket'
        end if

        best_relative = huge(1.0_dp)
        best_t = t_min
        bracket_found = .false.
        t_previous = t_min
        call evaluate_crossover_balance(core,rate_state,pressure,x, &
            diagnostics%branching_reactions, &
            diagnostics%stabilization_reactions,t_previous, &
            f_previous,k_branch,k_stabilization)
        best_relative = relative_balance(f_previous,k_branch,k_stabilization)

        do point = 2, scan_points
            t_current = t_min + (t_max-t_min)*real(point-1,dp)/ &
                real(scan_points-1,dp)
            call evaluate_crossover_balance(core,rate_state,pressure,x, &
                diagnostics%branching_reactions, &
                diagnostics%stabilization_reactions,t_current, &
                f_current,k_branch,k_stabilization)
            current_relative = &
                relative_balance(f_current,k_branch,k_stabilization)
            if (current_relative < best_relative) then
                best_relative = current_relative
                best_t = t_current
            end if
            if (opposite_or_zero_sign(f_previous,f_current)) then
                t_left = t_previous
                f_left = f_previous
                t_right = t_current
                f_right = f_current
                bracket_found = .true.
                exit
            end if
            t_previous = t_current
            f_previous = f_current
        end do

        if (.not. bracket_found) then
            crossover_temperature = best_t
            call populate_final_diagnostics(core,rate_state,pressure,x,best_t, &
                diagnostics)
            diagnostics%converged = .false.
            deallocate(x)
            return
        end if

        do iter = 1, max_root_iterations
            if (abs(f_right-f_left) > tiny(1.0_dp)) then
                t_trial = t_right-f_right*(t_right-t_left)/(f_right-f_left)
            else
                t_trial = 0.5_dp*(t_left+t_right)
            end if
            if (t_trial <= t_left .or. t_trial >= t_right .or. &
                .not. ieee_is_finite(t_trial)) then
                t_trial = 0.5_dp*(t_left+t_right)
            end if

            call evaluate_crossover_balance(core,rate_state,pressure,x, &
                diagnostics%branching_reactions, &
                diagnostics%stabilization_reactions,t_trial, &
                f_trial,k_branch,k_stabilization)
            current_relative = relative_balance(f_trial,k_branch,k_stabilization)
            if (current_relative <= root_tolerance .or. &
                abs(t_right-t_left) <= 1.0e-8_dp*max(1.0_dp,t_trial)) exit

            if (opposite_or_zero_sign(f_left,f_trial)) then
                t_right = t_trial
                f_right = f_trial
            else
                t_left = t_trial
                f_left = f_trial
            end if
        end do

        crossover_temperature = t_trial
        diagnostics%iterations = iter
        call populate_final_diagnostics(core,rate_state,pressure,x,t_trial, &
            diagnostics)
        diagnostics%converged = &
            diagnostics%relative_balance_residual <= root_tolerance
        deallocate(x)
    end subroutine hydrogen_crossover_temperature


    subroutine evaluate_crossover_balance(core,state,pressure,x, &
            branch_reactions,stabilization_reactions,temperature,balance, &
            k_branch,k_stabilization)
        type(chemical_kinetics_core), intent(in) :: core
        type(chemical_rate_state), intent(inout) :: state
        real(dp), intent(in) :: pressure,temperature
        real(dp), dimension(:), intent(in) :: x
        integer, dimension(:), intent(in) :: branch_reactions
        integer, dimension(:), intent(in) :: stabilization_reactions
        real(dp), intent(out) :: balance,k_branch,k_stabilization

        real(dp), allocatable :: concentration(:)
        integer :: i

        allocate(concentration(core%species_number))
        concentration = x*pressure/(r_gase_J*temperature)
        call core%prepare_rate_state(temperature,state)

        k_branch = 0.0_dp
        do i = 1, size(branch_reactions)
            k_branch = k_branch + core%effective_forward_rate_coefficient( &
                state,branch_reactions(i),concentration)
        end do
        k_stabilization = 0.0_dp
        do i = 1, size(stabilization_reactions)
            k_stabilization = k_stabilization + &
                core%effective_forward_rate_coefficient( &
                    state,stabilization_reactions(i),concentration)
        end do
        balance = k_branch-k_stabilization
        deallocate(concentration)
    end subroutine evaluate_crossover_balance


    subroutine populate_final_diagnostics(core,state,pressure,x,temperature, &
            diagnostics)
        type(chemical_kinetics_core), intent(in) :: core
        type(chemical_rate_state), intent(inout) :: state
        real(dp), intent(in) :: pressure,temperature
        real(dp), dimension(:), intent(in) :: x
        type(hydrogen_crossover_diagnostics), intent(inout) :: diagnostics

        real(dp), allocatable :: concentration(:)
        real(dp) :: balance
        integer :: i

        allocate(concentration(core%species_number))
        concentration = x*pressure/(r_gase_J*temperature)
        call core%prepare_rate_state(temperature,state)

        if (allocated(diagnostics%branching_details)) &
            deallocate(diagnostics%branching_details)
        if (allocated(diagnostics%stabilization_details)) &
            deallocate(diagnostics%stabilization_details)
        allocate(diagnostics%branching_details( &
            size(diagnostics%branching_reactions)))
        allocate(diagnostics%stabilization_details( &
            size(diagnostics%stabilization_reactions)))

        diagnostics%branching_rate_coefficient = 0.0_dp
        do i = 1, size(diagnostics%branching_reactions)
            diagnostics%branching_rate_coefficient = &
                diagnostics%branching_rate_coefficient + &
                core%effective_forward_rate_coefficient( &
                    state,diagnostics%branching_reactions(i),concentration, &
                    diagnostics%branching_details(i))
        end do
        diagnostics%stabilization_rate_coefficient = 0.0_dp
        do i = 1, size(diagnostics%stabilization_reactions)
            diagnostics%stabilization_rate_coefficient = &
                diagnostics%stabilization_rate_coefficient + &
                core%effective_forward_rate_coefficient( &
                    state,diagnostics%stabilization_reactions(i),concentration, &
                    diagnostics%stabilization_details(i))
        end do

        balance = diagnostics%branching_rate_coefficient - &
            diagnostics%stabilization_rate_coefficient
        diagnostics%relative_balance_residual = relative_balance(balance, &
            diagnostics%branching_rate_coefficient, &
            diagnostics%stabilization_rate_coefficient)
        deallocate(concentration)
    end subroutine populate_final_diagnostics


    real(dp) function relative_balance(balance,k_branch,k_stabilization) &
            result(relative_value)
        real(dp), intent(in) :: balance,k_branch,k_stabilization
        relative_value = abs(balance)/max(abs(k_branch),abs(k_stabilization), &
            tiny(1.0_dp))
    end function relative_balance


    integer function species_index(core,name) result(index_value)
        type(chemical_kinetics_core), intent(in) :: core
        character(len=*), intent(in) :: name
        integer :: specie

        index_value = 0
        do specie = 1, core%species_number
            if (trim(core%chemistry%species_names(specie)) == trim(name)) then
                index_value = specie
                return
            end if
        end do
    end function species_index


    pure logical function opposite_or_zero_sign(a,b) result(brackets)
        real(dp), intent(in) :: a,b

        brackets = (a <= 0.0_dp .and. b >= 0.0_dp) .or. &
            (a >= 0.0_dp .and. b <= 0.0_dp)
    end function opposite_or_zero_sign

end module chemical_kinetics_diagnostics
