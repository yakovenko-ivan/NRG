module chemical_kinetics_core_class

    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite

    use kind_parameters, only: dp
    use global_data, only: r_gase_J, P_atm
    use chemical_properties_class, only: chemical_properties
    use thermophysical_properties_class, only: thermophysical_properties

    implicit none

    private
    public :: chemical_kinetics_core, chemical_kinetics_core_c
    public :: chemical_rate_state, reaction_rate_diagnostics

    real(dp), parameter :: maximum_rate_temperature = 10000.0_dp
    integer, parameter :: maximum_reaction_components = 3
    integer, parameter :: maximum_net_reaction_species = &
        2*maximum_reaction_components

    type :: chemical_rate_state
        real(dp) :: temperature = 0.0_dp
        real(dp), allocatable :: high_pressure_rate(:)
        real(dp), allocatable :: low_pressure_rate(:)
        real(dp), allocatable :: reverse_factor(:)
        real(dp), allocatable :: troe_f_center(:)
        real(dp), allocatable :: troe_c(:)
        real(dp), allocatable :: troe_n(:)
        real(dp), allocatable :: entropy(:)
        real(dp), allocatable :: enthalpy(:)
    contains
        procedure :: clear => clear_rate_state
    end type chemical_rate_state

    type :: reaction_rate_diagnostics
        real(dp) :: high_pressure_rate = 0.0_dp
        real(dp) :: low_pressure_rate = 0.0_dp
        real(dp) :: third_body_concentration = 0.0_dp
        real(dp) :: reduced_pressure = 0.0_dp
        real(dp) :: troe_factor = 1.0_dp
        real(dp) :: forward_rate_coefficient = 0.0_dp
        real(dp) :: reverse_rate_coefficient = 0.0_dp
        real(dp) :: reverse_factor = 0.0_dp
    end type reaction_rate_diagnostics

    type :: chemical_kinetics_core
        type(chemical_properties), pointer :: chemistry => null()
        type(thermophysical_properties), pointer :: thermophysics => null()

        integer :: species_number = 0
        integer :: reactions_number = 0
        real(dp) :: default_third_body_efficiency = 0.0_dp
        logical :: has_any_third_body_reaction = .false.

        integer, allocatable :: reaction_type(:)
        logical, allocatable :: reaction_uses_third_body(:)
        logical, allocatable :: reaction_uses_falloff(:)
        logical, allocatable :: reaction_uses_troe(:)
        logical, allocatable :: reaction_reversible(:)
        integer, allocatable :: forward_third_body_power(:)
        integer, allocatable :: reverse_third_body_power(:)
        integer, allocatable :: reactant_count(:)
        integer, allocatable :: product_count(:)
        integer, allocatable :: net_count(:)
        integer, allocatable :: reactant_species(:,:)
        integer, allocatable :: reactant_multiplicity(:,:)
        integer, allocatable :: product_species(:,:)
        integer, allocatable :: product_multiplicity(:,:)
        integer, allocatable :: net_species(:,:)
        integer, allocatable :: net_stoichiometry(:,:)
        integer, allocatable :: third_body_offset(:)
        integer, allocatable :: third_body_species(:)
        real(dp), allocatable :: third_body_efficiency_delta(:)
    contains
        procedure :: prepare_rate_state
        procedure :: calculate_species_rates
        procedure :: third_body_concentration
        procedure :: effective_rate_constants
        procedure :: effective_forward_rate_coefficient
        procedure :: reaction_matches
        procedure :: reaction_equation_string
        procedure :: clear => clear_kinetics_core
    end type chemical_kinetics_core

    interface chemical_kinetics_core_c
        module procedure constructor
    end interface chemical_kinetics_core_c

contains

    type(chemical_kinetics_core) function constructor(chemistry,thermophysics)
        type(chemical_properties), target, intent(in) :: chemistry
        type(thermophysical_properties), target, intent(in) :: thermophysics

        constructor%chemistry => chemistry
        constructor%thermophysics => thermophysics
        constructor%species_number = chemistry%species_number
        constructor%reactions_number = chemistry%reactions_number

        if (constructor%species_number <= 0) then
            error stop 'Chemical kinetics core: mechanism contains no species'
        end if
        if (constructor%reactions_number <= 0) then
            error stop 'Chemical kinetics core: mechanism contains no reactions'
        end if

        call preprocess_mechanism(constructor)
    end function constructor


    subroutine clear_kinetics_core(this)
        class(chemical_kinetics_core), intent(inout) :: this

        nullify(this%chemistry)
        nullify(this%thermophysics)
        this%species_number = 0
        this%reactions_number = 0
        this%default_third_body_efficiency = 0.0_dp
        this%has_any_third_body_reaction = .false.

        if (allocated(this%reaction_type)) deallocate(this%reaction_type)
        if (allocated(this%reaction_uses_third_body)) &
            deallocate(this%reaction_uses_third_body)
        if (allocated(this%reaction_uses_falloff)) &
            deallocate(this%reaction_uses_falloff)
        if (allocated(this%reaction_uses_troe)) &
            deallocate(this%reaction_uses_troe)
        if (allocated(this%reaction_reversible)) &
            deallocate(this%reaction_reversible)
        if (allocated(this%forward_third_body_power)) &
            deallocate(this%forward_third_body_power)
        if (allocated(this%reverse_third_body_power)) &
            deallocate(this%reverse_third_body_power)
        if (allocated(this%reactant_count)) deallocate(this%reactant_count)
        if (allocated(this%product_count)) deallocate(this%product_count)
        if (allocated(this%net_count)) deallocate(this%net_count)
        if (allocated(this%reactant_species)) deallocate(this%reactant_species)
        if (allocated(this%reactant_multiplicity)) &
            deallocate(this%reactant_multiplicity)
        if (allocated(this%product_species)) deallocate(this%product_species)
        if (allocated(this%product_multiplicity)) &
            deallocate(this%product_multiplicity)
        if (allocated(this%net_species)) deallocate(this%net_species)
        if (allocated(this%net_stoichiometry)) &
            deallocate(this%net_stoichiometry)
        if (allocated(this%third_body_offset)) deallocate(this%third_body_offset)
        if (allocated(this%third_body_species)) deallocate(this%third_body_species)
        if (allocated(this%third_body_efficiency_delta)) &
            deallocate(this%third_body_efficiency_delta)
    end subroutine clear_kinetics_core


    subroutine clear_rate_state(this)
        class(chemical_rate_state), intent(inout) :: this

        this%temperature = 0.0_dp
        if (allocated(this%high_pressure_rate)) deallocate(this%high_pressure_rate)
        if (allocated(this%low_pressure_rate)) deallocate(this%low_pressure_rate)
        if (allocated(this%reverse_factor)) deallocate(this%reverse_factor)
        if (allocated(this%troe_f_center)) deallocate(this%troe_f_center)
        if (allocated(this%troe_c)) deallocate(this%troe_c)
        if (allocated(this%troe_n)) deallocate(this%troe_n)
        if (allocated(this%entropy)) deallocate(this%entropy)
        if (allocated(this%enthalpy)) deallocate(this%enthalpy)
    end subroutine clear_rate_state


    subroutine preprocess_mechanism(this)
        type(chemical_kinetics_core), intent(inout) :: this

        integer :: reaction, component, specie, index, position
        integer :: left_raw_count, right_raw_count, sparse_count
        integer :: net_position, net_value, reaction_kind
        real(dp) :: efficiency

        allocate(this%reaction_type(this%reactions_number))
        allocate(this%reaction_uses_third_body(this%reactions_number))
        allocate(this%reaction_uses_falloff(this%reactions_number))
        allocate(this%reaction_uses_troe(this%reactions_number))
        allocate(this%reaction_reversible(this%reactions_number))
        allocate(this%forward_third_body_power(this%reactions_number))
        allocate(this%reverse_third_body_power(this%reactions_number))
        allocate(this%reactant_count(this%reactions_number))
        allocate(this%product_count(this%reactions_number))
        allocate(this%net_count(this%reactions_number))
        allocate(this%reactant_species(maximum_reaction_components, &
            this%reactions_number))
        allocate(this%reactant_multiplicity(maximum_reaction_components, &
            this%reactions_number))
        allocate(this%product_species(maximum_reaction_components, &
            this%reactions_number))
        allocate(this%product_multiplicity(maximum_reaction_components, &
            this%reactions_number))
        allocate(this%net_species(maximum_net_reaction_species, &
            this%reactions_number))
        allocate(this%net_stoichiometry(maximum_net_reaction_species, &
            this%reactions_number))
        allocate(this%third_body_offset(0:this%reactions_number))

        this%default_third_body_efficiency = &
            this%chemistry%default_enhanced_efficiencies
        this%reaction_type = this%chemistry%reactions_type
        this%reaction_uses_third_body = .false.
        this%reaction_uses_falloff = .false.
        this%reaction_uses_troe = .false.
        this%reaction_reversible = .false.
        this%forward_third_body_power = 0
        this%reverse_third_body_power = 0
        this%reactant_count = 0
        this%product_count = 0
        this%net_count = 0
        this%reactant_species = 0
        this%reactant_multiplicity = 0
        this%product_species = 0
        this%product_multiplicity = 0
        this%net_species = 0
        this%net_stoichiometry = 0

        do reaction = 1, this%reactions_number
            reaction_kind = this%reaction_type(reaction)
            this%reaction_uses_falloff(reaction) = &
                reaction_kind == 2 .or. reaction_kind == 3 .or. &
                reaction_kind == 7
            this%reaction_uses_troe(reaction) = &
                reaction_kind == 3 .or. &
                (reaction_kind == 7 .and. &
                any(this%chemistry%Troe_coeffs(reaction,:) /= 0.0_dp))
            this%reaction_reversible(reaction) = &
                reaction_kind >= 0 .and. reaction_kind <= 4

            left_raw_count = this%chemistry%chemical_coeffs(1,reaction,1)
            right_raw_count = this%chemistry%chemical_coeffs(1,reaction,2)

            do component = 2, left_raw_count + 1
                specie = this%chemistry%chemical_coeffs(component,reaction,1)
                if (specie == this%species_number + 1) then
                    this%forward_third_body_power(reaction) = &
                        this%forward_third_body_power(reaction) + 1
                    cycle
                end if
                if (specie < 1 .or. specie > this%species_number) cycle
                position = 0
                do index = 1, this%reactant_count(reaction)
                    if (this%reactant_species(index,reaction) == specie) then
                        position = index
                        exit
                    end if
                end do
                if (position == 0) then
                    this%reactant_count(reaction) = &
                        this%reactant_count(reaction) + 1
                    position = this%reactant_count(reaction)
                    if (position > maximum_reaction_components) then
                        error stop 'Chemical kinetics core: too many reactants'
                    end if
                    this%reactant_species(position,reaction) = specie
                end if
                this%reactant_multiplicity(position,reaction) = &
                    this%reactant_multiplicity(position,reaction) + 1
            end do

            do component = 2, right_raw_count + 1
                specie = this%chemistry%chemical_coeffs(component,reaction,2)
                if (specie == this%species_number + 1) then
                    this%reverse_third_body_power(reaction) = &
                        this%reverse_third_body_power(reaction) + 1
                    cycle
                end if
                if (specie < 1 .or. specie > this%species_number) cycle
                position = 0
                do index = 1, this%product_count(reaction)
                    if (this%product_species(index,reaction) == specie) then
                        position = index
                        exit
                    end if
                end do
                if (position == 0) then
                    this%product_count(reaction) = &
                        this%product_count(reaction) + 1
                    position = this%product_count(reaction)
                    if (position > maximum_reaction_components) then
                        error stop 'Chemical kinetics core: too many products'
                    end if
                    this%product_species(position,reaction) = specie
                end if
                this%product_multiplicity(position,reaction) = &
                    this%product_multiplicity(position,reaction) + 1
            end do

            do index = 1, this%reactant_count(reaction)
                specie = this%reactant_species(index,reaction)
                net_value = -this%reactant_multiplicity(index,reaction)
                net_position = 0
                do position = 1, this%net_count(reaction)
                    if (this%net_species(position,reaction) == specie) then
                        net_position = position
                        exit
                    end if
                end do
                if (net_position == 0) then
                    this%net_count(reaction) = this%net_count(reaction) + 1
                    net_position = this%net_count(reaction)
                    if (net_position > maximum_net_reaction_species) then
                        error stop 'Chemical kinetics core: too many net species'
                    end if
                    this%net_species(net_position,reaction) = specie
                end if
                this%net_stoichiometry(net_position,reaction) = &
                    this%net_stoichiometry(net_position,reaction) + net_value
            end do

            do index = 1, this%product_count(reaction)
                specie = this%product_species(index,reaction)
                net_value = this%product_multiplicity(index,reaction)
                net_position = 0
                do position = 1, this%net_count(reaction)
                    if (this%net_species(position,reaction) == specie) then
                        net_position = position
                        exit
                    end if
                end do
                if (net_position == 0) then
                    this%net_count(reaction) = this%net_count(reaction) + 1
                    net_position = this%net_count(reaction)
                    if (net_position > maximum_net_reaction_species) then
                        error stop 'Chemical kinetics core: too many net species'
                    end if
                    this%net_species(net_position,reaction) = specie
                end if
                this%net_stoichiometry(net_position,reaction) = &
                    this%net_stoichiometry(net_position,reaction) + net_value
            end do

            this%reaction_uses_third_body(reaction) = &
                this%reaction_uses_falloff(reaction) .or. &
                this%forward_third_body_power(reaction) > 0 .or. &
                this%reverse_third_body_power(reaction) > 0 .or. &
                reaction_kind == 1 .or. reaction_kind == 6
        end do

        this%has_any_third_body_reaction = &
            any(this%reaction_uses_third_body)

        sparse_count = 0
        this%third_body_offset(0) = 1
        do reaction = 1, this%reactions_number
            if (this%reaction_uses_third_body(reaction)) then
                do specie = 1, this%species_number
                    efficiency = this%chemistry%enhanced_efficiencies( &
                        reaction,specie)
                    if (efficiency /= this%default_third_body_efficiency) then
                        sparse_count = sparse_count + 1
                    end if
                end do
            end if
            this%third_body_offset(reaction) = sparse_count + 1
        end do

        allocate(this%third_body_species(sparse_count))
        allocate(this%third_body_efficiency_delta(sparse_count))
        sparse_count = 0
        do reaction = 1, this%reactions_number
            if (.not. this%reaction_uses_third_body(reaction)) cycle
            do specie = 1, this%species_number
                efficiency = this%chemistry%enhanced_efficiencies(reaction,specie)
                if (efficiency == this%default_third_body_efficiency) cycle
                sparse_count = sparse_count + 1
                this%third_body_species(sparse_count) = specie
                this%third_body_efficiency_delta(sparse_count) = &
                    efficiency-this%default_third_body_efficiency
            end do
        end do
    end subroutine preprocess_mechanism


    subroutine prepare_rate_state(this,input_temperature,state)
        class(chemical_kinetics_core), intent(in) :: this
        real(dp), intent(in) :: input_temperature
        type(chemical_rate_state), intent(inout) :: state

        integer :: reaction, component, specie
        real(dp) :: temperature, entropy_change, enthalpy_change
        real(dp) :: mole_change, equilibrium_pressure, equilibrium_concentration
        real(dp) :: exponent_argument, alpha, temperature_1
        real(dp) :: temperature_2, temperature_3, f_center, log_f_center

        if (.not. associated(this%chemistry) .or. &
            .not. associated(this%thermophysics)) then
            error stop 'Chemical kinetics core: uninitialized core'
        end if

        temperature = min(input_temperature,maximum_rate_temperature)
        if (.not. ieee_is_finite(temperature) .or. temperature <= 0.0_dp) then
            error stop 'Chemical kinetics core: invalid rate temperature'
        end if

        call ensure_rate_state_size(state,this%species_number, &
            this%reactions_number)

        do specie = 1, this%species_number
            state%entropy(specie) = this%thermophysics%specie_entropy_molar( &
                temperature,specie)
            state%enthalpy(specie) = this%thermophysics%specie_enthalpy_molar( &
                temperature,specie)
        end do

        state%temperature = temperature
        do reaction = 1, this%reactions_number
            state%high_pressure_rate(reaction) = &
                this%chemistry%A(reaction)*temperature** &
                this%chemistry%beta(reaction)* &
                exp(-this%chemistry%E_act(reaction)/(r_gase_J*temperature))

            state%low_pressure_rate(reaction) = 0.0_dp
            if (this%chemistry%A_low(reaction) > 0.0_dp) then
                state%low_pressure_rate(reaction) = &
                    this%chemistry%A_low(reaction)*temperature** &
                    this%chemistry%beta_low(reaction)* &
                    exp(-this%chemistry%E_act_low(reaction)/ &
                    (r_gase_J*temperature))
            end if

            entropy_change = 0.0_dp
            enthalpy_change = 0.0_dp
            mole_change = 0.0_dp
            do component = 1, this%reactant_count(reaction)
                specie = this%reactant_species(component,reaction)
                entropy_change = entropy_change - real( &
                    this%reactant_multiplicity(component,reaction),dp)* &
                    state%entropy(specie)
                enthalpy_change = enthalpy_change - real( &
                    this%reactant_multiplicity(component,reaction),dp)* &
                    state%enthalpy(specie)
                mole_change = mole_change - real( &
                    this%reactant_multiplicity(component,reaction),dp)
            end do
            do component = 1, this%product_count(reaction)
                specie = this%product_species(component,reaction)
                entropy_change = entropy_change + real( &
                    this%product_multiplicity(component,reaction),dp)* &
                    state%entropy(specie)
                enthalpy_change = enthalpy_change + real( &
                    this%product_multiplicity(component,reaction),dp)* &
                    state%enthalpy(specie)
                mole_change = mole_change + real( &
                    this%product_multiplicity(component,reaction),dp)
            end do

            state%reverse_factor(reaction) = 0.0_dp
            if (this%reaction_reversible(reaction)) then
                exponent_argument = entropy_change/r_gase_J - &
                    enthalpy_change/(r_gase_J*temperature)
                equilibrium_pressure = exp(exponent_argument)
                equilibrium_concentration = equilibrium_pressure* &
                    (P_atm/(r_gase_J*temperature))**mole_change
                if (.not. ieee_is_finite(equilibrium_concentration) .or. &
                    equilibrium_concentration <= 0.0_dp) then
                    error stop 'Chemical kinetics core: invalid equilibrium constant'
                end if
                state%reverse_factor(reaction) = &
                    1.0_dp/equilibrium_concentration
            end if

            state%troe_f_center(reaction) = 1.0_dp
            state%troe_c(reaction) = 0.0_dp
            state%troe_n(reaction) = 1.0_dp
            if (this%reaction_uses_troe(reaction)) then
                alpha = this%chemistry%Troe_coeffs(reaction,1)
                temperature_1 = this%chemistry%Troe_coeffs(reaction,2)
                temperature_2 = this%chemistry%Troe_coeffs(reaction,3)
                temperature_3 = this%chemistry%Troe_coeffs(reaction,4)
                f_center = 0.0_dp
                if (temperature_1 > 0.0_dp) then
                    f_center = f_center + (1.0_dp-alpha)* &
                        exp(-temperature/temperature_1)
                end if
                if (temperature_2 > 0.0_dp) then
                    f_center = f_center + alpha* &
                        exp(-temperature/temperature_2)
                end if
                if (temperature_3 > 0.0_dp) then
                    f_center = f_center + exp(-temperature_3/temperature)
                end if
                f_center = min(max(f_center,tiny(1.0_dp)),1.0_dp)
                log_f_center = log10(f_center)
                state%troe_f_center(reaction) = f_center
                state%troe_c(reaction) = -0.4_dp-0.67_dp*log_f_center
                state%troe_n(reaction) = 0.75_dp-1.27_dp*log_f_center
            end if

            if (.not. ieee_is_finite(state%high_pressure_rate(reaction)) .or. &
                state%high_pressure_rate(reaction) < 0.0_dp) then
                error stop 'Chemical kinetics core: invalid high-pressure rate'
            end if
            if (.not. ieee_is_finite(state%low_pressure_rate(reaction)) .or. &
                state%low_pressure_rate(reaction) < 0.0_dp) then
                error stop 'Chemical kinetics core: invalid low-pressure rate'
            end if
        end do

    end subroutine prepare_rate_state


    subroutine ensure_rate_state_size(state,species_number,reactions_number)
        type(chemical_rate_state), intent(inout) :: state
        integer, intent(in) :: species_number,reactions_number

        logical :: rebuild

        rebuild = .not. allocated(state%high_pressure_rate)
        if (.not. rebuild) rebuild = &
            size(state%high_pressure_rate) /= reactions_number
        if (.not. rebuild) rebuild = .not. allocated(state%entropy)
        if (.not. rebuild) rebuild = size(state%entropy) /= species_number
        if (.not. rebuild) return

        call state%clear()
        allocate(state%high_pressure_rate(reactions_number))
        allocate(state%low_pressure_rate(reactions_number))
        allocate(state%reverse_factor(reactions_number))
        allocate(state%troe_f_center(reactions_number))
        allocate(state%troe_c(reactions_number))
        allocate(state%troe_n(reactions_number))
        allocate(state%entropy(species_number))
        allocate(state%enthalpy(species_number))
        state%high_pressure_rate = 0.0_dp
        state%low_pressure_rate = 0.0_dp
        state%reverse_factor = 0.0_dp
        state%troe_f_center = 1.0_dp
        state%troe_c = 0.0_dp
        state%troe_n = 1.0_dp
        state%entropy = 0.0_dp
        state%enthalpy = 0.0_dp
    end subroutine ensure_rate_state_size


    real(dp) function third_body_concentration(this,reaction,concentration, &
            total_concentration) result(third_body)
        class(chemical_kinetics_core), intent(in) :: this
        integer, intent(in) :: reaction
        real(dp), dimension(:), intent(in) :: concentration
        real(dp), intent(in), optional :: total_concentration

        integer :: index, specie_index
        real(dp) :: concentration_sum
        real(dp) :: cancellation_scale, correction, direct_sum

        if (reaction < 1 .or. reaction > this%reactions_number) then
            error stop 'Chemical kinetics core: invalid reaction index'
        end if
        if (size(concentration) /= this%species_number) then
            error stop 'Chemical kinetics core: concentration size mismatch'
        end if

        if (.not. this%reaction_uses_third_body(reaction)) then
            third_body = 0.0_dp
            return
        end if

        if (present(total_concentration)) then
            concentration_sum = total_concentration
        else
            concentration_sum = sum(max(concentration,0.0_dp))
        end if

        third_body = this%default_third_body_efficiency*concentration_sum
        cancellation_scale = abs(third_body)

        do index = this%third_body_offset(reaction-1), &
                this%third_body_offset(reaction)-1
            specie_index = this%third_body_species(index)
            correction = this%third_body_efficiency_delta(index)* &
                max(concentration(specie_index),0.0_dp)
            third_body = third_body + correction
            cancellation_scale = cancellation_scale + abs(correction)
        end do

        ! The sparse delta form above is fast for the usual case where most
        ! efficiencies equal the default value.  For collider-specific
        ! channels, however, it can evaluate an exact zero as
        !
        !   C_total - sum(C_k)
        !
        ! and leave a roundoff-sized positive residue.  This is particularly
        ! visible for the Ar-only / He-only duplicate HO2 channels when Ar and
        ! He are absent.  If the result is at the scale of the accumulated
        ! floating-point cancellation error, recompute the mathematically
        ! equivalent direct weighted sum.  This path is rare and avoids any
        ! chemistry-dependent cutoff.
        if (cancellation_scale > 0.0_dp) then
            if (abs(third_body) <= &
                    epsilon(1.0_dp)* &
                    real(max(1,1 + this%third_body_offset(reaction) - &
                        this%third_body_offset(reaction-1)),dp)* &
                    cancellation_scale) then

                direct_sum = 0.0_dp
                do specie_index = 1, this%species_number
                    direct_sum = direct_sum + &
                        this%chemistry%enhanced_efficiencies( &
                            reaction,specie_index)* &
                        max(concentration(specie_index),0.0_dp)
                end do
                third_body = direct_sum
            end if
        end if

        third_body = max(third_body,0.0_dp)
    end function third_body_concentration


    subroutine effective_rate_constants(this,state,reaction,third_body, &
            forward_constant,reverse_constant,diagnostics)
        class(chemical_kinetics_core), intent(in) :: this
        type(chemical_rate_state), intent(in) :: state
        integer, intent(in) :: reaction
        real(dp), intent(in) :: third_body
        real(dp), intent(out) :: forward_constant,reverse_constant
        type(reaction_rate_diagnostics), intent(out), optional :: diagnostics

        real(dp) :: high_rate, low_rate, reduced_pressure, troe_factor

        if (.not. allocated(state%high_pressure_rate)) then
            error stop 'Chemical kinetics core: rate state is not prepared'
        end if
        if (reaction < 1 .or. reaction > this%reactions_number) then
            error stop 'Chemical kinetics core: invalid reaction index'
        end if

        high_rate = state%high_pressure_rate(reaction)
        low_rate = state%low_pressure_rate(reaction)
        reduced_pressure = 0.0_dp
        troe_factor = 1.0_dp

        if (this%reaction_uses_falloff(reaction)) then
            if (high_rate <= 0.0_dp .or. low_rate <= 0.0_dp .or. &
                third_body <= 0.0_dp) then
                forward_constant = 0.0_dp
            else
                reduced_pressure = low_rate*third_body/high_rate
                forward_constant = high_rate*reduced_pressure/ &
                    (1.0_dp+reduced_pressure)
                if (this%reaction_uses_troe(reaction)) then
                    troe_factor = troe_falloff_factor( &
                        state%troe_f_center(reaction), &
                        state%troe_c(reaction),state%troe_n(reaction), &
                        reduced_pressure)
                    forward_constant = forward_constant*troe_factor
                end if
            end if
        else
            forward_constant = high_rate
        end if

        if (this%reaction_reversible(reaction)) then
            reverse_constant = forward_constant*state%reverse_factor(reaction)
        else
            reverse_constant = 0.0_dp
        end if

        if (present(diagnostics)) then
            diagnostics%high_pressure_rate = high_rate
            diagnostics%low_pressure_rate = low_rate
            diagnostics%third_body_concentration = third_body
            diagnostics%reduced_pressure = reduced_pressure
            diagnostics%troe_factor = troe_factor
            diagnostics%forward_rate_coefficient = forward_constant
            diagnostics%reverse_rate_coefficient = reverse_constant
            diagnostics%reverse_factor = state%reverse_factor(reaction)
        end if
    end subroutine effective_rate_constants


    real(dp) function effective_forward_rate_coefficient(this,state,reaction, &
            concentration,diagnostics) result(forward_constant)
        class(chemical_kinetics_core), intent(in) :: this
        type(chemical_rate_state), intent(in) :: state
        integer, intent(in) :: reaction
        real(dp), dimension(:), intent(in) :: concentration
        type(reaction_rate_diagnostics), intent(out), optional :: diagnostics

        real(dp) :: reverse_constant, third_body, total_concentration

        total_concentration = sum(max(concentration,0.0_dp))
        third_body = this%third_body_concentration( &
            reaction,concentration,total_concentration)
        call this%effective_rate_constants(state,reaction,third_body, &
            forward_constant,reverse_constant,diagnostics)
    end function effective_forward_rate_coefficient


    subroutine calculate_species_rates(this,state,concentration,concentration_rate)
        class(chemical_kinetics_core), intent(in) :: this
        type(chemical_rate_state), intent(in) :: state
        real(dp), dimension(:), intent(in) :: concentration
        real(dp), dimension(:), intent(out) :: concentration_rate

        integer :: reaction, component, multiplicity, specie_index
        integer :: stoichiometry
        real(dp) :: forward_rate, reverse_rate, net_rate
        real(dp) :: forward_constant, reverse_constant
        real(dp) :: third_body, total_concentration, positive_concentration

        if (size(concentration) /= this%species_number .or. &
            size(concentration_rate) /= this%species_number) then
            error stop 'Chemical kinetics core: species-rate array size mismatch'
        end if
        if (.not. allocated(state%high_pressure_rate)) then
            error stop 'Chemical kinetics core: rate state is not prepared'
        end if

        concentration_rate = 0.0_dp
        total_concentration = 0.0_dp
        if (this%has_any_third_body_reaction) then
            total_concentration = sum(max(concentration,0.0_dp))
        end if

        do reaction = 1, this%reactions_number
            third_body = 0.0_dp
            if (this%reaction_uses_third_body(reaction)) then
                third_body = this%third_body_concentration( &
                    reaction,concentration,total_concentration)
            end if

            call this%effective_rate_constants(state,reaction,third_body, &
                forward_constant,reverse_constant)

            forward_rate = forward_constant
            do component = 1, this%reactant_count(reaction)
                specie_index = this%reactant_species(component,reaction)
                positive_concentration = max(concentration(specie_index),0.0_dp)
                do multiplicity = 1, &
                        this%reactant_multiplicity(component,reaction)
                    forward_rate = forward_rate*positive_concentration
                end do
            end do
            do multiplicity = 1, this%forward_third_body_power(reaction)
                forward_rate = forward_rate*third_body
            end do

            reverse_rate = reverse_constant
            do component = 1, this%product_count(reaction)
                specie_index = this%product_species(component,reaction)
                positive_concentration = max(concentration(specie_index),0.0_dp)
                do multiplicity = 1, &
                        this%product_multiplicity(component,reaction)
                    reverse_rate = reverse_rate*positive_concentration
                end do
            end do
            do multiplicity = 1, this%reverse_third_body_power(reaction)
                reverse_rate = reverse_rate*third_body
            end do

            net_rate = forward_rate-reverse_rate

            do component = 1, this%net_count(reaction)
                specie_index = this%net_species(component,reaction)
                stoichiometry = this%net_stoichiometry(component,reaction)
                concentration_rate(specie_index) = &
                    concentration_rate(specie_index) + &
                    real(stoichiometry,dp)*net_rate
            end do
        end do
    end subroutine calculate_species_rates


    logical function reaction_matches(this,reaction,reactants,reactant_nu, &
            products,product_nu) result(matches)
        class(chemical_kinetics_core), intent(in) :: this
        integer, intent(in) :: reaction
        integer, dimension(:), intent(in) :: reactants,reactant_nu
        integer, dimension(:), intent(in) :: products,product_nu

        integer :: i,j
        logical :: found

        matches = .false.
        if (size(reactants) /= size(reactant_nu) .or. &
            size(products) /= size(product_nu)) return
        if (reaction < 1 .or. reaction > this%reactions_number) return
        if (this%reactant_count(reaction) /= size(reactants)) return
        if (this%product_count(reaction) /= size(products)) return

        do i = 1, size(reactants)
            found = .false.
            do j = 1, this%reactant_count(reaction)
                if (this%reactant_species(j,reaction) == reactants(i) .and. &
                    this%reactant_multiplicity(j,reaction) == reactant_nu(i)) then
                    found = .true.
                    exit
                end if
            end do
            if (.not. found) return
        end do

        do i = 1, size(products)
            found = .false.
            do j = 1, this%product_count(reaction)
                if (this%product_species(j,reaction) == products(i) .and. &
                    this%product_multiplicity(j,reaction) == product_nu(i)) then
                    found = .true.
                    exit
                end if
            end do
            if (.not. found) return
        end do

        matches = .true.
    end function reaction_matches


    character(len=256) function reaction_equation_string(this,reaction) &
            result(equation)
        class(chemical_kinetics_core), intent(in) :: this
        integer, intent(in) :: reaction

        integer :: component, multiplicity
        character(len=20) :: specie_name
        character(len=16) :: coefficient

        equation = ''
        if (reaction < 1 .or. reaction > this%reactions_number) return

        do component = 1, this%reactant_count(reaction)
            if (component > 1) equation = trim(equation)//' + '
            multiplicity = this%reactant_multiplicity(component,reaction)
            if (multiplicity > 1) then
                write(coefficient,'(I0)') multiplicity
                equation = trim(equation)//trim(coefficient)
            end if
            specie_name = this%chemistry%species_names( &
                this%reactant_species(component,reaction))
            equation = trim(equation)//trim(specie_name)
        end do

        if (this%reaction_reversible(reaction)) then
            equation = trim(equation)//' = '
        else
            equation = trim(equation)//' => '
        end if

        do component = 1, this%product_count(reaction)
            if (component > 1) equation = trim(equation)//' + '
            multiplicity = this%product_multiplicity(component,reaction)
            if (multiplicity > 1) then
                write(coefficient,'(I0)') multiplicity
                equation = trim(equation)//trim(coefficient)
            end if
            specie_name = this%chemistry%species_names( &
                this%product_species(component,reaction))
            equation = trim(equation)//trim(specie_name)
        end do
    end function reaction_equation_string


    pure real(dp) function troe_falloff_factor(f_center,c_troe,n_troe, &
            reduced_pressure) result(factor)
        real(dp), intent(in) :: f_center,c_troe,n_troe,reduced_pressure

        real(dp), parameter :: d_troe = 0.14_dp
        real(dp) :: log_pressure, denominator, exponent

        if (reduced_pressure <= 0.0_dp) then
            factor = 1.0_dp
            return
        end if

        log_pressure = log10(reduced_pressure)
        denominator = n_troe-d_troe*(log_pressure+c_troe)
        if (abs(denominator) <= tiny(1.0_dp)) then
            factor = f_center
            return
        end if
        exponent = 1.0_dp/(1.0_dp+ &
            ((log_pressure+c_troe)/denominator)**2)
        factor = f_center**exponent
    end function troe_falloff_factor

end module chemical_kinetics_core_class
