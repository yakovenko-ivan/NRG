module chemical_kinetics_solver_class

    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use, intrinsic :: iso_fortran_env, only: error_unit, output_unit, int64

    use kind_parameters, only: dp
    use global_data, only: r_gase_J, P_atm, T_ref, task_setup_folder, &
        fold_sep, chemical_mechanisms_folder
    use field_pointers
    use boundary_conditions_class
    use data_manager_class
    use computational_domain_class
    use thermophysical_properties_class
    use chemical_properties_class

    implicit none

    private
    public :: chemical_kinetics_solver, chemical_kinetics_solver_c

    real(dp), parameter :: default_activation_temperature = 305.0_dp
    real(dp), parameter :: default_slatec_accuracy = 1.0e-9_dp
    real(dp), parameter :: default_slatec_error_weight = 1.0e-4_dp
    real(dp), parameter :: default_slatec_max_internal_step = 1.0e-7_dp
    integer, parameter :: default_slatec_max_steps = 10000
    real(dp), parameter :: default_negative_concentration_tolerance = 1.0e-10_dp
    real(dp), parameter :: default_mass_balance_tolerance = 1.0e-6_dp
    real(dp), parameter :: default_table_start_temperature = 300.5_dp
    real(dp), parameter :: default_table_start_temperature_max = 310.0_dp
    real(dp), parameter :: maximum_rate_temperature = 10000.0_dp
    integer, parameter :: maximum_reaction_components = 3
    integer, parameter :: maximum_net_reaction_species = &
        2*maximum_reaction_components


    ! DDRIV3 accepts an external right-hand-side callback and provides no user
    ! context argument.  The only module state retained by the refactored solver
    ! is therefore one private workspace per OpenMP thread.  Registered fields,
    ! configuration and table data remain owned by each solver instance.
    type :: kinetics_thread_workspace
        integer :: species_number = 0
        integer :: reactions_number = 0
        type(chemical_properties), pointer :: chemistry => null()
        real(dp) :: temperature = 0.0_dp
        real(dp) :: default_third_body_efficiency = 0.0_dp
        logical :: has_any_third_body_reaction = .false.

        real(dp), allocatable :: concentration_initial(:)
        real(dp), allocatable :: concentration_final(:)
        real(dp), allocatable :: high_pressure_rate(:)
        real(dp), allocatable :: low_pressure_rate(:)
        real(dp), allocatable :: reverse_factor(:)
        real(dp), allocatable :: troe_f_center(:)
        real(dp), allocatable :: troe_c(:)
        real(dp), allocatable :: troe_n(:)
        real(dp), allocatable :: entropy(:)
        real(dp), allocatable :: enthalpy(:)

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

        real(dp), allocatable :: mass_fraction_cell(:)
        real(dp), allocatable :: species_source_cell(:)
        real(dp), allocatable :: concentration_increment_cell(:)
        real(dp), allocatable :: work(:)
        integer, allocatable :: iwork(:)
    end type kinetics_thread_workspace

    type(kinetics_thread_workspace), save :: thread_workspace
!$omp threadprivate(thread_workspace)

    interface
        subroutine ddriv3(n, t, y, f, nstate, tout, ntask, nroot, eps, &
            ewt, ierror, mint, miter, impl, ml, mu, mxord, hmax, work, &
            lenw, iwork, leniw, jacobn, fa, nde, mxstep, g, users, ierflg)
            implicit none
            external :: f, jacobn, fa, users
            double precision, external :: g
            double precision :: eps, ewt(*), hmax, t, tout, work(*), y(*)
            integer :: ierror, ierflg, impl, iwork(*), leniw, lenw
            integer :: mint, miter, ml, mu, mxord, mxstep, n, nde
            integer :: nroot, nstate, ntask
        end subroutine ddriv3
    end interface

    type :: chemical_kinetics_solver
        private

        type(field_scalar_cons_pointer) :: temperature
        type(field_scalar_cons_pointer) :: density
        type(field_scalar_cons_pointer) :: energy_source
        type(field_vector_cons_pointer) :: mass_fraction
        type(field_vector_cons_pointer) :: species_source

        type(field_scalar_cons), pointer :: energy_source_store => null()
        type(field_vector_cons), pointer :: species_source_store => null()

        type(computational_domain) :: domain
        type(thermophysical_properties_pointer) :: thermophysics
        type(chemical_properties_pointer) :: chemistry
        type(boundary_conditions_pointer) :: boundary

        integer :: species_number = 0
        integer :: reactions_number = 0

        real(dp), allocatable :: inverse_molar_mass(:)
        real(dp), allocatable :: reference_enthalpy_mass(:)
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

        character(len=20) :: ode_solver = 'slatec'
        real(dp) :: activation_temperature = default_activation_temperature
        real(dp) :: slatec_accuracy = default_slatec_accuracy
        real(dp) :: slatec_error_weight = default_slatec_error_weight
        real(dp) :: slatec_max_internal_step = &
            default_slatec_max_internal_step
        integer :: slatec_max_steps = default_slatec_max_steps
        real(dp) :: negative_concentration_tolerance = &
            default_negative_concentration_tolerance
        real(dp) :: mass_balance_tolerance = default_mass_balance_tolerance
        real(dp) :: table_start_temperature = &
            default_table_start_temperature
        real(dp) :: table_start_temperature_max = &
            default_table_start_temperature_max

        logical :: record_concentration_increment = .false.
        real(dp), allocatable :: concentration_increment(:,:,:,:)
        real(dp), allocatable :: table_temperature(:)
        real(dp), allocatable :: table_concentration_increment(:,:)

        integer(int64) :: total_solve_calls = 0_int64
        integer(int64) :: total_active_cells = 0_int64
        integer(int64) :: total_integrator_calls = 0_int64
        integer(int64) :: total_internal_steps = 0_int64
        integer(int64) :: total_rhs_evaluations = 0_int64
        integer(int64) :: total_jacobian_evaluations = 0_int64
        integer(int64) :: last_active_cells = 0_int64
        integer(int64) :: last_integrator_calls = 0_int64
        integer(int64) :: last_internal_steps = 0_int64
        integer(int64) :: last_rhs_evaluations = 0_int64
        integer(int64) :: last_jacobian_evaluations = 0_int64
        real(dp) :: total_rate_preparation_time = 0.0_dp
        real(dp) :: total_integration_time = 0.0_dp
        real(dp) :: total_source_assembly_time = 0.0_dp
        real(dp) :: last_rate_preparation_time = 0.0_dp
        real(dp) :: last_integration_time = 0.0_dp
        real(dp) :: last_source_assembly_time = 0.0_dp
    contains
        procedure :: solve_chemical_kinetics
        procedure :: write_chemical_kinetics_table
        procedure :: set_activation_temperature
        procedure :: set_slatec_controls
        procedure :: configure_table_approximated
        procedure :: use_detailed_kinetics
        procedure :: set_concentration_increment_recording
        procedure :: reset_performance_statistics
        procedure :: write_performance_statistics
        procedure, private :: preprocess_mechanism
        procedure, private :: allocate_concentration_increment_storage
        procedure, private :: ensure_thread_workspace
        procedure, private :: prepare_cell_rate_coefficients
        procedure, private :: solve_cell_detailed_kinetics
        procedure, private :: assemble_cell_sources
        procedure, private :: read_chemical_kinetics_table
        procedure, private :: interpolate_table_increment
        procedure, private :: validate_configuration
    end type chemical_kinetics_solver

    interface chemical_kinetics_solver_c
        module procedure constructor
    end interface chemical_kinetics_solver_c

contains

    type(chemical_kinetics_solver) function constructor(manager, &
            activation_temperature, ode_solver, table_file, &
            record_concentration_increment)
        type(data_manager), intent(inout) :: manager
        real(dp), intent(in), optional :: activation_temperature
        character(len=*), intent(in), optional :: ode_solver
        character(len=*), intent(in), optional :: table_file
        logical, intent(in), optional :: record_concentration_increment

        type(field_scalar_cons_pointer) :: scalar_pointer
        type(field_vector_cons_pointer) :: vector_pointer
        type(field_tensor_cons_pointer) :: tensor_pointer
        integer :: specie

        constructor%domain = manager%domain
        constructor%boundary%bc_ptr => &
            manager%boundary_conditions_pointer%bc_ptr
        constructor%thermophysics%thermo_ptr => &
            manager%thermophysics%thermo_ptr
        constructor%chemistry%chem_ptr => manager%chemistry%chem_ptr

        constructor%species_number = &
            constructor%chemistry%chem_ptr%species_number
        constructor%reactions_number = &
            constructor%chemistry%chem_ptr%reactions_number

        call manager%get_cons_field_pointer_by_name( &
            scalar_pointer, vector_pointer, tensor_pointer, 'temperature')
        constructor%temperature%s_ptr => scalar_pointer%s_ptr
        call manager%get_cons_field_pointer_by_name( &
            scalar_pointer, vector_pointer, tensor_pointer, 'density')
        constructor%density%s_ptr => scalar_pointer%s_ptr
        call manager%get_cons_field_pointer_by_name( &
            scalar_pointer, vector_pointer, tensor_pointer, &
            'specie_mass_fraction')
        constructor%mass_fraction%v_ptr => vector_pointer%v_ptr

        allocate(constructor%energy_source_store)
        allocate(constructor%species_source_store)

        call manager%create_scalar_field( &
            constructor%energy_source_store, &
            'energy_production_chemistry', 'E_f_prod_chem')
        constructor%energy_source%s_ptr => constructor%energy_source_store

        call manager%create_vector_field( &
            constructor%species_source_store, &
            'specie_production_chemistry', 'Y_prod_chem', 'chemical')
        constructor%species_source%v_ptr => constructor%species_source_store

        constructor%energy_source%s_ptr%cells = 0.0_dp
        call zero_vector_field(constructor%species_source%v_ptr)

        if (present(record_concentration_increment)) then
            call constructor%set_concentration_increment_recording( &
                record_concentration_increment)
        end if

        allocate(constructor%inverse_molar_mass(constructor%species_number))
        allocate(constructor%reference_enthalpy_mass( &
            constructor%species_number))
        do specie = 1, constructor%species_number
            if (constructor%thermophysics%thermo_ptr%molar_masses(specie) <= &
                0.0_dp) then
                error stop 'Chemical kinetics: non-positive species molar mass'
            end if
            constructor%inverse_molar_mass(specie) = 1.0_dp/ &
                constructor%thermophysics%thermo_ptr%molar_masses(specie)
            constructor%reference_enthalpy_mass(specie) = &
                constructor%thermophysics%thermo_ptr%specie_enthalpy_molar( &
                    T_ref,specie)*constructor%inverse_molar_mass(specie)
        end do

        call constructor%preprocess_mechanism()

        if (present(activation_temperature)) then
            call constructor%set_activation_temperature(activation_temperature)
        end if

        if (present(ode_solver)) then
            select case (trim(adjustl(ode_solver)))
            case ('slatec')
                call constructor%use_detailed_kinetics()
            case ('table_approximated')
                if (.not. present(table_file)) then
                    error stop 'Chemical kinetics: table file was not supplied'
                end if
                call constructor%configure_table_approximated(table_file)
            case default
                error stop 'Chemical kinetics: unsupported ODE solver'
            end select
        else if (present(table_file)) then
            call constructor%configure_table_approximated(table_file)
        end if

        call constructor%validate_configuration()
    end function constructor


    subroutine validate_configuration(this)
        class(chemical_kinetics_solver), intent(in) :: this

        if (this%species_number <= 0) then
            error stop 'Chemical kinetics: mechanism contains no species'
        end if
        if (this%reactions_number <= 0) then
            error stop 'Chemical kinetics: mechanism contains no reactions'
        end if
        if (.not. ieee_is_finite(this%activation_temperature) .or. &
            this%activation_temperature <= 0.0_dp) then
            error stop 'Chemical kinetics: invalid activation temperature'
        end if
        if (.not. ieee_is_finite(this%slatec_accuracy) .or. &
            this%slatec_accuracy <= 0.0_dp) then
            error stop 'Chemical kinetics: invalid SLATEC accuracy'
        end if
        if (.not. ieee_is_finite(this%slatec_error_weight) .or. &
            this%slatec_error_weight <= 0.0_dp) then
            error stop 'Chemical kinetics: invalid SLATEC error weight'
        end if
        if (.not. ieee_is_finite(this%slatec_max_internal_step) .or. &
            this%slatec_max_internal_step <= 0.0_dp) then
            error stop 'Chemical kinetics: invalid maximum internal step'
        end if
        if (this%slatec_max_steps <= 0) then
            error stop 'Chemical kinetics: invalid maximum step count'
        end if
    end subroutine validate_configuration


    subroutine set_activation_temperature(this, temperature)
        class(chemical_kinetics_solver), intent(inout) :: this
        real(dp), intent(in) :: temperature

        if (.not. ieee_is_finite(temperature) .or. temperature <= 0.0_dp) then
            error stop 'Chemical kinetics: activation temperature must be positive'
        end if
        this%activation_temperature = temperature
    end subroutine set_activation_temperature


    subroutine set_slatec_controls(this, accuracy, error_weight, &
            maximum_internal_step, maximum_steps)
        class(chemical_kinetics_solver), intent(inout) :: this
        real(dp), intent(in), optional :: accuracy
        real(dp), intent(in), optional :: error_weight
        real(dp), intent(in), optional :: maximum_internal_step
        integer, intent(in), optional :: maximum_steps

        if (present(accuracy)) this%slatec_accuracy = accuracy
        if (present(error_weight)) this%slatec_error_weight = error_weight
        if (present(maximum_internal_step)) then
            this%slatec_max_internal_step = maximum_internal_step
        end if
        if (present(maximum_steps)) this%slatec_max_steps = maximum_steps
        call this%validate_configuration()
    end subroutine set_slatec_controls


    subroutine configure_table_approximated(this, table_file)
        class(chemical_kinetics_solver), intent(inout) :: this
        character(len=*), intent(in) :: table_file

        call this%read_chemical_kinetics_table(table_file)
        this%ode_solver = 'table_approximated'
    end subroutine configure_table_approximated


    subroutine use_detailed_kinetics(this)
        class(chemical_kinetics_solver), intent(inout) :: this

        this%ode_solver = 'slatec'
    end subroutine use_detailed_kinetics

    subroutine set_concentration_increment_recording(this, enabled)
        class(chemical_kinetics_solver), intent(inout) :: this
        logical, intent(in) :: enabled

        this%record_concentration_increment = enabled
        if (enabled) then
            call this%allocate_concentration_increment_storage()
        else if (allocated(this%concentration_increment)) then
            deallocate(this%concentration_increment)
        end if
    end subroutine set_concentration_increment_recording


    subroutine allocate_concentration_increment_storage(this)
        class(chemical_kinetics_solver), intent(inout) :: this
        integer, dimension(3,2) :: allocation_bounds

        if (allocated(this%concentration_increment)) return
        allocation_bounds = this%domain%get_local_utter_cells_bounds()
        allocate(this%concentration_increment(this%species_number, &
            allocation_bounds(1,1):allocation_bounds(1,2), &
            allocation_bounds(2,1):allocation_bounds(2,2), &
            allocation_bounds(3,1):allocation_bounds(3,2)))
        this%concentration_increment = 0.0_dp
    end subroutine allocate_concentration_increment_storage


    subroutine reset_performance_statistics(this)
        class(chemical_kinetics_solver), intent(inout) :: this

        this%total_solve_calls = 0_int64
        this%total_active_cells = 0_int64
        this%total_integrator_calls = 0_int64
        this%total_internal_steps = 0_int64
        this%total_rhs_evaluations = 0_int64
        this%total_jacobian_evaluations = 0_int64
        this%last_active_cells = 0_int64
        this%last_integrator_calls = 0_int64
        this%last_internal_steps = 0_int64
        this%last_rhs_evaluations = 0_int64
        this%last_jacobian_evaluations = 0_int64
        this%total_rate_preparation_time = 0.0_dp
        this%total_integration_time = 0.0_dp
        this%total_source_assembly_time = 0.0_dp
        this%last_rate_preparation_time = 0.0_dp
        this%last_integration_time = 0.0_dp
        this%last_source_assembly_time = 0.0_dp
    end subroutine reset_performance_statistics


    subroutine write_performance_statistics(this, unit)
        class(chemical_kinetics_solver), intent(in) :: this
        integer, intent(in), optional :: unit
        integer :: output

        output = output_unit
        if (present(unit)) output = unit
        write(output,'(A)') 'Chemical-kinetics performance statistics'
        write(output,'(A,I0)') '  solve calls: ',this%total_solve_calls
        write(output,'(A,I0)') '  active cells: ',this%total_active_cells
        write(output,'(A,I0)') '  integrator calls: ', &
            this%total_integrator_calls
        write(output,'(A,I0)') '  internal steps: ',this%total_internal_steps
        write(output,'(A,I0)') '  RHS evaluations: ', &
            this%total_rhs_evaluations
        write(output,'(A,I0)') '  Jacobian evaluations: ', &
            this%total_jacobian_evaluations
#ifdef CHEMISTRY_PROFILE
        write(output,'(A,ES14.6)') '  accumulated rate preparation time [s]: ', &
            this%total_rate_preparation_time
        write(output,'(A,ES14.6)') '  accumulated integration time [s]: ', &
            this%total_integration_time
        write(output,'(A,ES14.6)') '  accumulated source assembly time [s]: ', &
            this%total_source_assembly_time
#else
        write(output,'(A)') '  timing disabled; compile with CHEMISTRY_PROFILE'
#endif
    end subroutine write_performance_statistics


    subroutine preprocess_mechanism(this)
        class(chemical_kinetics_solver), intent(inout) :: this

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
            this%chemistry%chem_ptr%default_enhanced_efficiencies
        this%reaction_type = this%chemistry%chem_ptr%reactions_type
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

        associate(chemistry => this%chemistry%chem_ptr)
        do reaction = 1, this%reactions_number
            reaction_kind = this%reaction_type(reaction)
            this%reaction_uses_falloff(reaction) = &
                reaction_kind == 2 .or. reaction_kind == 3 .or. &
                reaction_kind == 7
            this%reaction_uses_troe(reaction) = &
                reaction_kind == 3 .or. &
                (reaction_kind == 7 .and. &
                any(chemistry%Troe_coeffs(reaction,:) /= 0.0_dp))
            this%reaction_reversible(reaction) = &
                reaction_kind >= 0 .and. reaction_kind <= 4

            left_raw_count = chemistry%chemical_coeffs(1,reaction,1)
            right_raw_count = chemistry%chemical_coeffs(1,reaction,2)

            do component = 2, left_raw_count + 1
                specie = chemistry%chemical_coeffs(component,reaction,1)
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
                        error stop 'Chemical kinetics: too many reactants'
                    end if
                    this%reactant_species(position,reaction) = specie
                end if
                this%reactant_multiplicity(position,reaction) = &
                    this%reactant_multiplicity(position,reaction) + 1
            end do

            do component = 2, right_raw_count + 1
                specie = chemistry%chemical_coeffs(component,reaction,2)
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
                        error stop 'Chemical kinetics: too many products'
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
                    efficiency = chemistry%enhanced_efficiencies( &
                        reaction,specie)
                    if (efficiency /= &
                        this%default_third_body_efficiency) then
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
                efficiency = chemistry%enhanced_efficiencies(reaction,specie)
                if (efficiency == this%default_third_body_efficiency) cycle
                sparse_count = sparse_count + 1
                this%third_body_species(sparse_count) = specie
                this%third_body_efficiency_delta(sparse_count) = &
                    efficiency-this%default_third_body_efficiency
            end do
        end do
        end associate
    end subroutine preprocess_mechanism



    subroutine solve_chemical_kinetics(this, time_step)
        class(chemical_kinetics_solver), intent(inout) :: this
        real(dp), intent(in) :: time_step

        integer, dimension(3,2) :: cell_loop
        integer :: i, j, k, specie
        real(dp) :: density_cell, temperature_cell, energy_source_cell
        real(dp) :: mass_fraction_sum, negative_tolerance
        integer(int64) :: active_cells_step, integrator_calls_step
        integer(int64) :: internal_steps_step, rhs_evaluations_step
        integer(int64) :: jacobian_evaluations_step
        integer(int64) :: cell_internal_steps, cell_rhs_evaluations
        integer(int64) :: cell_jacobian_evaluations
        real(dp) :: rate_preparation_time_step, integration_time_step
        real(dp) :: source_assembly_time_step
        real(dp) :: cell_rate_preparation_time, cell_integration_time
#ifdef CHEMISTRY_PROFILE
        real(dp) :: source_time_start
#endif

        if (.not. ieee_is_finite(time_step) .or. time_step <= 0.0_dp) then
            error stop 'Chemical kinetics: time step must be finite and positive'
        end if

        call this%validate_configuration()
        cell_loop = this%domain%get_local_inner_cells_bounds()

        this%energy_source%s_ptr%cells = 0.0_dp
        call zero_vector_field(this%species_source%v_ptr)
        if (this%record_concentration_increment) then
            call this%allocate_concentration_increment_storage()
            this%concentration_increment = 0.0_dp
        end if

        active_cells_step = 0_int64
        integrator_calls_step = 0_int64
        internal_steps_step = 0_int64
        rhs_evaluations_step = 0_int64
        jacobian_evaluations_step = 0_int64
        rate_preparation_time_step = 0.0_dp
        integration_time_step = 0.0_dp
        source_assembly_time_step = 0.0_dp

        associate( &
            temperature_field => this%temperature%s_ptr, &
            density_field => this%density%s_ptr, &
            mass_fraction_field => this%mass_fraction%v_ptr, &
            energy_source_field => this%energy_source%s_ptr, &
            species_source_field => this%species_source%v_ptr, &
            markers => this%boundary%bc_ptr%bc_markers)

!$omp parallel default(shared) &
!$omp private(i,j,k,specie,density_cell,temperature_cell,energy_source_cell) &
!$omp private(mass_fraction_sum,negative_tolerance) &
!$omp private(cell_internal_steps,cell_rhs_evaluations) &
!$omp private(cell_jacobian_evaluations,cell_rate_preparation_time) &
!$omp private(cell_integration_time) &
#ifdef CHEMISTRY_PROFILE
!$omp private(source_time_start) &
#endif
!$omp reduction(+:active_cells_step,integrator_calls_step) &
!$omp reduction(+:internal_steps_step,rhs_evaluations_step) &
!$omp reduction(+:jacobian_evaluations_step) &
!$omp reduction(+:rate_preparation_time_step,integration_time_step) &
!$omp reduction(+:source_assembly_time_step)
        call this%ensure_thread_workspace()

!$omp do collapse(3) schedule(guided,1)
        do k = cell_loop(3,1), cell_loop(3,2)
            do j = cell_loop(2,1), cell_loop(2,2)
                do i = cell_loop(1,1), cell_loop(1,2)
                    if (markers(i,j,k) /= 0) cycle

                    density_cell = density_field%cells(i,j,k)
                    temperature_cell = temperature_field%cells(i,j,k)
                    if (.not. ieee_is_finite(density_cell) .or. &
                        density_cell <= 0.0_dp) then
                        call report_invalid_cell('non-positive density',i,j,k)
                    end if
                    if (.not. ieee_is_finite(temperature_cell) .or. &
                        temperature_cell <= 0.0_dp) then
                        call report_invalid_cell('non-positive temperature',i,j,k)
                    end if
                    if (temperature_cell < this%activation_temperature) cycle
                    active_cells_step = active_cells_step + 1_int64

                    mass_fraction_sum = 0.0_dp
                    negative_tolerance = &
                        this%negative_concentration_tolerance
                    do specie = 1, this%species_number
                        thread_workspace%mass_fraction_cell(specie) = &
                            mass_fraction_field%pr(specie)%cells(i,j,k)
                        if (.not. ieee_is_finite( &
                            thread_workspace%mass_fraction_cell(specie))) then
                            call report_invalid_cell( &
                                'non-finite mass fraction',i,j,k)
                        end if
                        if (thread_workspace%mass_fraction_cell(specie) < &
                            -negative_tolerance) then
                            call report_invalid_cell('negative mass fraction',i,j,k)
                        end if
                        thread_workspace%mass_fraction_cell(specie) = max( &
                            thread_workspace%mass_fraction_cell(specie),0.0_dp)
                        mass_fraction_sum = mass_fraction_sum + &
                            thread_workspace%mass_fraction_cell(specie)
                    end do
                    if (mass_fraction_sum <= tiny(1.0_dp)) then
                        call report_invalid_cell('empty composition',i,j,k)
                    end if
                    thread_workspace%mass_fraction_cell = &
                        thread_workspace%mass_fraction_cell/mass_fraction_sum

                    cell_internal_steps = 0_int64
                    cell_rhs_evaluations = 0_int64
                    cell_jacobian_evaluations = 0_int64
                    cell_rate_preparation_time = 0.0_dp
                    cell_integration_time = 0.0_dp
                    select case (this%ode_solver)
                    case ('slatec')
                        call this%solve_cell_detailed_kinetics( &
                            density_cell,temperature_cell, &
                            thread_workspace%mass_fraction_cell,time_step, &
                            i,j,k,thread_workspace%concentration_increment_cell, &
                            cell_internal_steps,cell_rhs_evaluations, &
                            cell_jacobian_evaluations, &
                            cell_rate_preparation_time,cell_integration_time)
                        integrator_calls_step = integrator_calls_step + 1_int64
                        internal_steps_step = internal_steps_step + &
                            cell_internal_steps
                        rhs_evaluations_step = rhs_evaluations_step + &
                            cell_rhs_evaluations
                        jacobian_evaluations_step = &
                            jacobian_evaluations_step + &
                            cell_jacobian_evaluations
                        rate_preparation_time_step = &
                            rate_preparation_time_step + &
                            cell_rate_preparation_time
                        integration_time_step = integration_time_step + &
                            cell_integration_time
                    case ('table_approximated')
                        call this%interpolate_table_increment( &
                            temperature_cell, &
                            thread_workspace%concentration_increment_cell)
                    case default
                        call report_invalid_cell('unsupported ODE solver',i,j,k)
                    end select

#ifdef CHEMISTRY_PROFILE
                    source_time_start = chemistry_wall_time()
#endif
                    call this%assemble_cell_sources( &
                        density_cell,thread_workspace%mass_fraction_cell, &
                        time_step, &
                        thread_workspace%concentration_increment_cell, &
                        thread_workspace%species_source_cell, &
                        energy_source_cell,i,j,k)
#ifdef CHEMISTRY_PROFILE
                    source_assembly_time_step = source_assembly_time_step + &
                        chemistry_wall_time()-source_time_start
#endif

                    energy_source_field%cells(i,j,k) = energy_source_cell
                    do specie = 1, this%species_number
                        species_source_field%pr(specie)%cells(i,j,k) = &
                            thread_workspace%species_source_cell(specie)
                        if (this%record_concentration_increment) then
                            this%concentration_increment(specie,i,j,k) = &
                                thread_workspace%concentration_increment_cell(specie)
                        end if
                    end do
                end do
            end do
        end do
!$omp end do
!$omp end parallel

        end associate

        this%last_active_cells = active_cells_step
        this%last_integrator_calls = integrator_calls_step
        this%last_internal_steps = internal_steps_step
        this%last_rhs_evaluations = rhs_evaluations_step
        this%last_jacobian_evaluations = jacobian_evaluations_step
        this%last_rate_preparation_time = rate_preparation_time_step
        this%last_integration_time = integration_time_step
        this%last_source_assembly_time = source_assembly_time_step
        this%total_solve_calls = this%total_solve_calls + 1_int64
        this%total_active_cells = this%total_active_cells + active_cells_step
        this%total_integrator_calls = this%total_integrator_calls + &
            integrator_calls_step
        this%total_internal_steps = this%total_internal_steps + &
            internal_steps_step
        this%total_rhs_evaluations = this%total_rhs_evaluations + &
            rhs_evaluations_step
        this%total_jacobian_evaluations = &
            this%total_jacobian_evaluations + jacobian_evaluations_step
        this%total_rate_preparation_time = &
            this%total_rate_preparation_time + rate_preparation_time_step
        this%total_integration_time = this%total_integration_time + &
            integration_time_step
        this%total_source_assembly_time = &
            this%total_source_assembly_time + source_assembly_time_step
    end subroutine solve_chemical_kinetics


    subroutine ensure_thread_workspace(this)
        class(chemical_kinetics_solver), intent(in) :: this

        logical :: rebuild
        integer :: work_size

        rebuild = thread_workspace%species_number /= this%species_number .or. &
            thread_workspace%reactions_number /= this%reactions_number
        if (.not. associated(thread_workspace%chemistry, &
            this%chemistry%chem_ptr)) rebuild = .true.

        if (.not. rebuild) return

        call clear_thread_workspace()
        thread_workspace%species_number = this%species_number
        thread_workspace%reactions_number = this%reactions_number
        thread_workspace%chemistry => this%chemistry%chem_ptr
        thread_workspace%default_third_body_efficiency = &
            this%default_third_body_efficiency
        thread_workspace%has_any_third_body_reaction = &
            this%has_any_third_body_reaction

        allocate(thread_workspace%concentration_initial(this%species_number))
        allocate(thread_workspace%concentration_final(this%species_number))
        allocate(thread_workspace%high_pressure_rate(this%reactions_number))
        allocate(thread_workspace%low_pressure_rate(this%reactions_number))
        allocate(thread_workspace%reverse_factor(this%reactions_number))
        allocate(thread_workspace%troe_f_center(this%reactions_number))
        allocate(thread_workspace%troe_c(this%reactions_number))
        allocate(thread_workspace%troe_n(this%reactions_number))
        allocate(thread_workspace%entropy(this%species_number))
        allocate(thread_workspace%enthalpy(this%species_number))

        allocate(thread_workspace%reaction_type(this%reactions_number), &
            source=this%reaction_type)
        allocate(thread_workspace%reaction_uses_third_body( &
            this%reactions_number),source=this%reaction_uses_third_body)
        allocate(thread_workspace%reaction_uses_falloff( &
            this%reactions_number),source=this%reaction_uses_falloff)
        allocate(thread_workspace%reaction_uses_troe( &
            this%reactions_number),source=this%reaction_uses_troe)
        allocate(thread_workspace%reaction_reversible( &
            this%reactions_number),source=this%reaction_reversible)
        allocate(thread_workspace%forward_third_body_power( &
            this%reactions_number),source=this%forward_third_body_power)
        allocate(thread_workspace%reverse_third_body_power( &
            this%reactions_number),source=this%reverse_third_body_power)
        allocate(thread_workspace%reactant_count(this%reactions_number), &
            source=this%reactant_count)
        allocate(thread_workspace%product_count(this%reactions_number), &
            source=this%product_count)
        allocate(thread_workspace%net_count(this%reactions_number), &
            source=this%net_count)
        allocate(thread_workspace%reactant_species( &
            maximum_reaction_components,this%reactions_number), &
            source=this%reactant_species)
        allocate(thread_workspace%reactant_multiplicity( &
            maximum_reaction_components,this%reactions_number), &
            source=this%reactant_multiplicity)
        allocate(thread_workspace%product_species( &
            maximum_reaction_components,this%reactions_number), &
            source=this%product_species)
        allocate(thread_workspace%product_multiplicity( &
            maximum_reaction_components,this%reactions_number), &
            source=this%product_multiplicity)
        allocate(thread_workspace%net_species( &
            maximum_net_reaction_species,this%reactions_number), &
            source=this%net_species)
        allocate(thread_workspace%net_stoichiometry( &
            maximum_net_reaction_species,this%reactions_number), &
            source=this%net_stoichiometry)
        allocate(thread_workspace%third_body_offset( &
            0:this%reactions_number),source=this%third_body_offset)
        allocate(thread_workspace%third_body_species( &
            size(this%third_body_species)),source=this%third_body_species)
        allocate(thread_workspace%third_body_efficiency_delta( &
            size(this%third_body_efficiency_delta)), &
            source=this%third_body_efficiency_delta)

        allocate(thread_workspace%mass_fraction_cell(this%species_number))
        allocate(thread_workspace%species_source_cell(this%species_number))
        allocate(thread_workspace%concentration_increment_cell( &
            this%species_number))

        work_size = this%species_number**2 + 10*this%species_number + 250
        allocate(thread_workspace%work(work_size))
        allocate(thread_workspace%iwork(50 + this%species_number))

        thread_workspace%concentration_initial = 0.0_dp
        thread_workspace%concentration_final = 0.0_dp
        thread_workspace%high_pressure_rate = 0.0_dp
        thread_workspace%low_pressure_rate = 0.0_dp
        thread_workspace%reverse_factor = 0.0_dp
        thread_workspace%troe_f_center = 1.0_dp
        thread_workspace%troe_c = 0.0_dp
        thread_workspace%troe_n = 1.0_dp
        thread_workspace%entropy = 0.0_dp
        thread_workspace%enthalpy = 0.0_dp
        thread_workspace%mass_fraction_cell = 0.0_dp
        thread_workspace%species_source_cell = 0.0_dp
        thread_workspace%concentration_increment_cell = 0.0_dp
        thread_workspace%work = 0.0_dp
        thread_workspace%iwork = 0
    end subroutine ensure_thread_workspace


    subroutine clear_thread_workspace()
        if (allocated(thread_workspace%concentration_initial)) &
            deallocate(thread_workspace%concentration_initial)
        if (allocated(thread_workspace%concentration_final)) &
            deallocate(thread_workspace%concentration_final)
        if (allocated(thread_workspace%high_pressure_rate)) &
            deallocate(thread_workspace%high_pressure_rate)
        if (allocated(thread_workspace%low_pressure_rate)) &
            deallocate(thread_workspace%low_pressure_rate)
        if (allocated(thread_workspace%reverse_factor)) &
            deallocate(thread_workspace%reverse_factor)
        if (allocated(thread_workspace%troe_f_center)) &
            deallocate(thread_workspace%troe_f_center)
        if (allocated(thread_workspace%troe_c)) &
            deallocate(thread_workspace%troe_c)
        if (allocated(thread_workspace%troe_n)) &
            deallocate(thread_workspace%troe_n)
        if (allocated(thread_workspace%entropy)) &
            deallocate(thread_workspace%entropy)
        if (allocated(thread_workspace%enthalpy)) &
            deallocate(thread_workspace%enthalpy)
        if (allocated(thread_workspace%reaction_type)) &
            deallocate(thread_workspace%reaction_type)
        if (allocated(thread_workspace%reaction_uses_third_body)) &
            deallocate(thread_workspace%reaction_uses_third_body)
        if (allocated(thread_workspace%reaction_uses_falloff)) &
            deallocate(thread_workspace%reaction_uses_falloff)
        if (allocated(thread_workspace%reaction_uses_troe)) &
            deallocate(thread_workspace%reaction_uses_troe)
        if (allocated(thread_workspace%reaction_reversible)) &
            deallocate(thread_workspace%reaction_reversible)
        if (allocated(thread_workspace%forward_third_body_power)) &
            deallocate(thread_workspace%forward_third_body_power)
        if (allocated(thread_workspace%reverse_third_body_power)) &
            deallocate(thread_workspace%reverse_third_body_power)
        if (allocated(thread_workspace%reactant_count)) &
            deallocate(thread_workspace%reactant_count)
        if (allocated(thread_workspace%product_count)) &
            deallocate(thread_workspace%product_count)
        if (allocated(thread_workspace%net_count)) &
            deallocate(thread_workspace%net_count)
        if (allocated(thread_workspace%reactant_species)) &
            deallocate(thread_workspace%reactant_species)
        if (allocated(thread_workspace%reactant_multiplicity)) &
            deallocate(thread_workspace%reactant_multiplicity)
        if (allocated(thread_workspace%product_species)) &
            deallocate(thread_workspace%product_species)
        if (allocated(thread_workspace%product_multiplicity)) &
            deallocate(thread_workspace%product_multiplicity)
        if (allocated(thread_workspace%net_species)) &
            deallocate(thread_workspace%net_species)
        if (allocated(thread_workspace%net_stoichiometry)) &
            deallocate(thread_workspace%net_stoichiometry)
        if (allocated(thread_workspace%third_body_offset)) &
            deallocate(thread_workspace%third_body_offset)
        if (allocated(thread_workspace%third_body_species)) &
            deallocate(thread_workspace%third_body_species)
        if (allocated(thread_workspace%third_body_efficiency_delta)) &
            deallocate(thread_workspace%third_body_efficiency_delta)
        if (allocated(thread_workspace%mass_fraction_cell)) &
            deallocate(thread_workspace%mass_fraction_cell)
        if (allocated(thread_workspace%species_source_cell)) &
            deallocate(thread_workspace%species_source_cell)
        if (allocated(thread_workspace%concentration_increment_cell)) &
            deallocate(thread_workspace%concentration_increment_cell)
        if (allocated(thread_workspace%work)) &
            deallocate(thread_workspace%work)
        if (allocated(thread_workspace%iwork)) &
            deallocate(thread_workspace%iwork)
        nullify(thread_workspace%chemistry)
        thread_workspace%species_number = 0
        thread_workspace%reactions_number = 0
        thread_workspace%temperature = 0.0_dp
        thread_workspace%default_third_body_efficiency = 0.0_dp
        thread_workspace%has_any_third_body_reaction = .false.
    end subroutine clear_thread_workspace


    subroutine solve_cell_detailed_kinetics(this, density, temperature, &
            mass_fraction, time_step, i_cell, j_cell, k_cell, &
            concentration_increment, internal_steps, rhs_evaluations, &
            jacobian_evaluations, rate_preparation_time, integration_time)
        class(chemical_kinetics_solver), intent(in) :: this
        real(dp), intent(in) :: density, temperature, time_step
        real(dp), dimension(:), intent(in) :: mass_fraction
        integer, intent(in) :: i_cell, j_cell, k_cell
        real(dp), dimension(:), intent(out) :: concentration_increment
        integer(int64), intent(out) :: internal_steps, rhs_evaluations
        integer(int64), intent(out) :: jacobian_evaluations
        real(dp), intent(out) :: rate_preparation_time, integration_time

        integer :: n, nstate, ntask, nroot, ierror, mint, miter, impl
        integer :: ml, mu, mxord, lenw, leniw, nde, ierflg
        real(dp) :: eps, hmax, time_in, time_out
        real(dp), dimension(1) :: error_weight
        integer :: specie
        real(dp) :: concentration_scale, negative_limit
#ifdef CHEMISTRY_PROFILE
        real(dp) :: timer_start
#endif

        call this%ensure_thread_workspace()

        do specie = 1, this%species_number
            thread_workspace%concentration_initial(specie) = &
                density*mass_fraction(specie)* &
                this%inverse_molar_mass(specie)
        end do
        thread_workspace%concentration_final = &
            thread_workspace%concentration_initial
        thread_workspace%temperature = min(temperature,maximum_rate_temperature)

        rate_preparation_time = 0.0_dp
        integration_time = 0.0_dp
#ifdef CHEMISTRY_PROFILE
        timer_start = chemistry_wall_time()
#endif
        call this%prepare_cell_rate_coefficients(temperature)
#ifdef CHEMISTRY_PROFILE
        rate_preparation_time = chemistry_wall_time()-timer_start
        timer_start = chemistry_wall_time()
#endif

        n = this%species_number
        nstate = 1
        ntask = 1
        nroot = 0
        eps = this%slatec_accuracy
        error_weight(1) = this%slatec_error_weight
        ierror = 3
        mint = 2
        miter = 2
        impl = 0
        ml = 0
        mu = 0
        mxord = 5
        time_in = 0.0_dp
        time_out = time_step
        hmax = min(time_step,this%slatec_max_internal_step)
        lenw = size(thread_workspace%work)
        leniw = size(thread_workspace%iwork)
        nde = n
        ierflg = 0

        ! NSTATE=1 instructs DDRIV3 to initialize a new problem.  The complete
        ! WORK and IWORK arrays are initialized only when the thread workspace is
        ! allocated; clearing O(N_s^2) storage for every reactive cell is
        ! unnecessary and was not done by the original NRG implementation.
        call ddriv3(n,time_in,thread_workspace%concentration_final, &
            kinetics_rhs,nstate,time_out,ntask,nroot,eps,error_weight, &
            ierror,mint,miter,impl,ml,mu,mxord,hmax, &
            thread_workspace%work,lenw,thread_workspace%iwork,leniw, &
            kinetics_rhs,kinetics_rhs,nde,this%slatec_max_steps, &
            dummy_root,kinetics_rhs,ierflg)
#ifdef CHEMISTRY_PROFILE
        integration_time = chemistry_wall_time()-timer_start
#endif

        internal_steps = int(max(thread_workspace%iwork(3),0),int64)
        rhs_evaluations = int(max(thread_workspace%iwork(4),0),int64)
        jacobian_evaluations = int(max(thread_workspace%iwork(5),0),int64)

        if (ierflg /= 0 .or. nstate <= 0) then
            call report_slatec_failure(i_cell,j_cell,k_cell,nstate,ierflg)
        end if

        concentration_scale = max(1.0_dp, &
            maxval(thread_workspace%concentration_initial))
        negative_limit = this%negative_concentration_tolerance* &
            concentration_scale
        do specie = 1, this%species_number
            if (.not. ieee_is_finite( &
                thread_workspace%concentration_final(specie))) then
                call report_invalid_cell( &
                    'non-finite final concentration',i_cell,j_cell,k_cell)
            end if
            if (thread_workspace%concentration_final(specie) < &
                -negative_limit) then
                call report_invalid_cell( &
                    'negative final concentration',i_cell,j_cell,k_cell)
            end if
            thread_workspace%concentration_final(specie) = max( &
                thread_workspace%concentration_final(specie),0.0_dp)
            concentration_increment(specie) = &
                thread_workspace%concentration_final(specie) - &
                thread_workspace%concentration_initial(specie)
        end do
    end subroutine solve_cell_detailed_kinetics


    subroutine prepare_cell_rate_coefficients(this, input_temperature)
        class(chemical_kinetics_solver), intent(in) :: this
        real(dp), intent(in) :: input_temperature

        integer :: reaction, component, specie
        real(dp) :: temperature, entropy_change, enthalpy_change
        real(dp) :: mole_change, equilibrium_pressure, equilibrium_concentration
        real(dp) :: exponent_argument, alpha, temperature_1
        real(dp) :: temperature_2, temperature_3, f_center, log_f_center

        temperature = min(input_temperature,maximum_rate_temperature)
        if (temperature <= 0.0_dp) then
            error stop 'Chemical kinetics: invalid rate temperature'
        end if

        do specie = 1, this%species_number
            thread_workspace%entropy(specie) = &
                this%thermophysics%thermo_ptr%specie_entropy_molar( &
                    temperature,specie)
            thread_workspace%enthalpy(specie) = &
                this%thermophysics%thermo_ptr%specie_enthalpy_molar( &
                    temperature,specie)
        end do

        associate(chemistry => this%chemistry%chem_ptr)
        do reaction = 1, this%reactions_number
            thread_workspace%high_pressure_rate(reaction) = &
                chemistry%A(reaction)*temperature**chemistry%beta(reaction)* &
                exp(-chemistry%E_act(reaction)/(r_gase_J*temperature))
            thread_workspace%low_pressure_rate(reaction) = 0.0_dp
            if (chemistry%A_low(reaction) > 0.0_dp) then
                thread_workspace%low_pressure_rate(reaction) = &
                    chemistry%A_low(reaction)* &
                    temperature**chemistry%beta_low(reaction)* &
                    exp(-chemistry%E_act_low(reaction)/ &
                    (r_gase_J*temperature))
            end if

            entropy_change = 0.0_dp
            enthalpy_change = 0.0_dp
            mole_change = 0.0_dp
            do component = 1, this%reactant_count(reaction)
                specie = this%reactant_species(component,reaction)
                entropy_change = entropy_change - real( &
                    this%reactant_multiplicity(component,reaction),dp)* &
                    thread_workspace%entropy(specie)
                enthalpy_change = enthalpy_change - real( &
                    this%reactant_multiplicity(component,reaction),dp)* &
                    thread_workspace%enthalpy(specie)
                mole_change = mole_change - real( &
                    this%reactant_multiplicity(component,reaction),dp)
            end do
            do component = 1, this%product_count(reaction)
                specie = this%product_species(component,reaction)
                entropy_change = entropy_change + real( &
                    this%product_multiplicity(component,reaction),dp)* &
                    thread_workspace%entropy(specie)
                enthalpy_change = enthalpy_change + real( &
                    this%product_multiplicity(component,reaction),dp)* &
                    thread_workspace%enthalpy(specie)
                mole_change = mole_change + real( &
                    this%product_multiplicity(component,reaction),dp)
            end do

            thread_workspace%reverse_factor(reaction) = 0.0_dp
            if (this%reaction_reversible(reaction)) then
                exponent_argument = entropy_change/r_gase_J - &
                    enthalpy_change/(r_gase_J*temperature)
                equilibrium_pressure = exp(exponent_argument)
                equilibrium_concentration = equilibrium_pressure* &
                    (P_atm/(r_gase_J*temperature))**mole_change
                if (.not. ieee_is_finite(equilibrium_concentration) .or. &
                    equilibrium_concentration <= 0.0_dp) then
                    error stop 'Chemical kinetics: invalid equilibrium constant'
                end if
                thread_workspace%reverse_factor(reaction) = &
                    1.0_dp/equilibrium_concentration
            end if

            thread_workspace%troe_f_center(reaction) = 1.0_dp
            thread_workspace%troe_c(reaction) = 0.0_dp
            thread_workspace%troe_n(reaction) = 1.0_dp
            if (this%reaction_uses_troe(reaction)) then
                alpha = chemistry%Troe_coeffs(reaction,1)
                temperature_1 = chemistry%Troe_coeffs(reaction,2)
                temperature_2 = chemistry%Troe_coeffs(reaction,3)
                temperature_3 = chemistry%Troe_coeffs(reaction,4)
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
                thread_workspace%troe_f_center(reaction) = f_center
                thread_workspace%troe_c(reaction) = &
                    -0.4_dp-0.67_dp*log_f_center
                thread_workspace%troe_n(reaction) = &
                    0.75_dp-1.27_dp*log_f_center
            end if

            if (.not. ieee_is_finite( &
                thread_workspace%high_pressure_rate(reaction)) .or. &
                thread_workspace%high_pressure_rate(reaction) < 0.0_dp) then
                error stop 'Chemical kinetics: invalid high-pressure rate'
            end if
            if (.not. ieee_is_finite( &
                thread_workspace%low_pressure_rate(reaction)) .or. &
                thread_workspace%low_pressure_rate(reaction) < 0.0_dp) then
                error stop 'Chemical kinetics: invalid low-pressure rate'
            end if
        end do
        end associate
    end subroutine prepare_cell_rate_coefficients


    subroutine assemble_cell_sources(this, density, mass_fraction, time_step, &
            concentration_increment, species_source, energy_source, &
            i_cell, j_cell, k_cell)
        class(chemical_kinetics_solver), intent(in) :: this
        real(dp), intent(in) :: density, time_step
        real(dp), dimension(:), intent(in) :: mass_fraction
        real(dp), dimension(:), intent(inout) :: concentration_increment
        real(dp), dimension(:), intent(out) :: species_source
        real(dp), intent(out) :: energy_source
        integer, intent(in) :: i_cell, j_cell, k_cell

        integer :: specie, correction_specie
        real(dp) :: source_sum, source_l1
        real(dp) :: integrated_mass_defect, integrated_mass_activity
        real(dp) :: integrated_mass_tolerance, composition_sum
        real(dp) :: integrator_mass_floor

        do specie = 1, this%species_number
            species_source(specie) = concentration_increment(specie)/ &
                (this%inverse_molar_mass(specie)*time_step)
        end do

        ! Assess conservation in integrated mass units rather than in source-rate
        ! units.  The latter magnifies roundoff from
        !
        !   concentration_final - concentration_initial
        !
        ! by 1/time_step and can falsely reject otherwise conservative DDRIV3
        ! solutions at small CFD time steps.
        source_sum = sum(species_source)
        source_l1 = sum(abs(species_source))
        integrated_mass_defect = source_sum*time_step
        integrated_mass_activity = source_l1*time_step

        ! The first term detects a defect relative to the chemistry-induced mass
        ! redistribution.  The second term allows for the requested ODE accuracy
        ! and subtraction roundoff relative to the cell mass.
        integrator_mass_floor = density*max( &
            10.0_dp*this%slatec_accuracy, &
            1000.0_dp*epsilon(1.0_dp))
        integrated_mass_tolerance = max( &
            this%mass_balance_tolerance*integrated_mass_activity, &
            integrator_mass_floor)

        if (abs(integrated_mass_defect) > integrated_mass_tolerance) then
            call report_mass_imbalance( &
                i_cell,j_cell,k_cell,source_sum,source_l1, &
                integrated_mass_defect,integrated_mass_tolerance)
        end if

        ! Remove the integrator-level residual.  Normalize the correction weights
        ! defensively even though the caller already normalizes mass fractions.
        composition_sum = sum(mass_fraction)
        if (composition_sum <= tiny(1.0_dp)) then
            call report_invalid_cell('empty correction composition', &
                i_cell,j_cell,k_cell)
        end if
        species_source = species_source - &
            mass_fraction*(source_sum/composition_sum)

        ! Eliminate the final summation residual in the dominant component.
        correction_specie = maxloc(mass_fraction,dim=1)
        species_source(correction_specie) = &
            species_source(correction_specie) - sum(species_source)

        energy_source = 0.0_dp
        do specie = 1, this%species_number
            energy_source = energy_source - &
                this%reference_enthalpy_mass(specie)*species_source(specie)
            concentration_increment(specie) = species_source(specie)* &
                time_step*this%inverse_molar_mass(specie)
        end do

        if (.not. ieee_is_finite(energy_source)) then
            call report_invalid_cell( &
                'non-finite chemistry energy source',i_cell,j_cell,k_cell)
        end if
    end subroutine assemble_cell_sources


    subroutine write_chemical_kinetics_table(this, table_file)
        class(chemical_kinetics_solver), intent(in) :: this
        character(len=*), intent(in) :: table_file

        integer, dimension(3,2) :: cell_loop
        integer :: io_unit, i, j, k, specie, start_index, peak_index
        integer :: candidate, step, point_count
        real(dp) :: peak_temperature, previous_temperature
        real(dp) :: candidate_distance
        real(dp), allocatable :: row(:)
        character(len=1024) :: header
        character(len=20) :: specie_name

        if (.not. this%record_concentration_increment .or. &
            .not. allocated(this%concentration_increment)) then
            error stop 'Chemical kinetics table: increment recording is disabled'
        end if

        cell_loop = this%domain%get_local_inner_cells_bounds()
        j = (cell_loop(2,1) + cell_loop(2,2))/2
        k = (cell_loop(3,1) + cell_loop(3,2))/2

        peak_index = cell_loop(1,1)
        peak_temperature = -huge(1.0_dp)
        do i = cell_loop(1,1), cell_loop(1,2)
            if (this%temperature%s_ptr%cells(i,j,k) > peak_temperature) then
                peak_temperature = this%temperature%s_ptr%cells(i,j,k)
                peak_index = i
            end if
        end do

        start_index = 0
        candidate_distance = huge(1.0_dp)
        do i = cell_loop(1,1), cell_loop(1,2)
            if (this%temperature%s_ptr%cells(i,j,k) >= &
                    this%table_start_temperature .and. &
                this%temperature%s_ptr%cells(i,j,k) <= &
                    this%table_start_temperature_max) then
                if (real(abs(i-peak_index),dp) < candidate_distance) then
                    candidate_distance = real(abs(i-peak_index),dp)
                    start_index = i
                end if
            end if
        end do
        if (start_index == 0) then
            error stop 'Chemical kinetics table: no preheat-side start cell'
        end if

        step = merge(1,-1,peak_index >= start_index)
        allocate(row(this%species_number + 1))
        open(newunit=io_unit, file=trim(task_setup_folder)//trim(fold_sep)// &
            trim(chemical_mechanisms_folder)//trim(fold_sep)// &
            trim(table_file), status='replace', form='formatted')

        header = 'VARIABLES="T"'
        do specie = 1, this%species_number
            specie_name = this%chemistry%chem_ptr%get_chemical_specie_name( &
                specie)
            header = trim(header)//' "dC_'//trim(specie_name)//'"'
        end do
        write(io_unit,'(A)') trim(header)

        previous_temperature = -huge(1.0_dp)
        point_count = 0
        candidate = start_index
        do
            row(1) = this%temperature%s_ptr%cells(candidate,j,k)
            if (row(1) > previous_temperature + &
                100.0_dp*epsilon(max(1.0_dp,abs(row(1))))) then
                do specie = 1, this%species_number
                    row(specie+1) = &
                        this%concentration_increment(specie,candidate,j,k)
                end do
                write(io_unit,'(*(ES24.16E3,1X))') row
                previous_temperature = row(1)
                point_count = point_count + 1
            end if
            if (candidate == peak_index) exit
            candidate = candidate + step
        end do
        close(io_unit)
        deallocate(row)

        if (point_count < 2) then
            error stop 'Chemical kinetics table: fewer than two monotone points'
        end if
    end subroutine write_chemical_kinetics_table


    subroutine read_chemical_kinetics_table(this, table_file)
        class(chemical_kinetics_solver), intent(inout) :: this
        character(len=*), intent(in) :: table_file

        integer :: io_unit, io_status, table_size, point
        real(dp), allocatable :: row(:)
        character(len=2048) :: header
        character(len=:), allocatable :: file_path

        file_path = trim(task_setup_folder)//trim(fold_sep)// &
            trim(chemical_mechanisms_folder)//trim(fold_sep)// &
            trim(table_file)
        allocate(row(this%species_number + 1))

        open(newunit=io_unit,file=file_path,status='old',form='formatted', &
            action='read',iostat=io_status)
        if (io_status /= 0) then
            error stop 'Chemical kinetics: unable to open kinetics table'
        end if
        read(io_unit,'(A)',iostat=io_status) header
        if (io_status /= 0) then
            error stop 'Chemical kinetics: unable to read kinetics table header'
        end if

        table_size = 0
        do
            read(io_unit,*,iostat=io_status) row
            if (io_status /= 0) exit
            table_size = table_size + 1
        end do
        if (table_size < 2) then
            error stop 'Chemical kinetics: kinetics table is too short'
        end if

        if (allocated(this%table_temperature)) &
            deallocate(this%table_temperature)
        if (allocated(this%table_concentration_increment)) &
            deallocate(this%table_concentration_increment)
        allocate(this%table_temperature(table_size))
        allocate(this%table_concentration_increment( &
            this%species_number,table_size))

        rewind(io_unit)
        read(io_unit,'(A)') header
        do point = 1, table_size
            read(io_unit,*,iostat=io_status) row
            if (io_status /= 0) then
                error stop 'Chemical kinetics: malformed kinetics table row'
            end if
            this%table_temperature(point) = row(1)
            this%table_concentration_increment(:,point) = &
                row(2:this%species_number+1)
            if (point > 1) then
                if (this%table_temperature(point) <= &
                    this%table_temperature(point-1)) then
                    error stop 'Chemical kinetics: table temperatures not increasing'
                end if
            end if
        end do
        close(io_unit)
        deallocate(row)
    end subroutine read_chemical_kinetics_table


    subroutine interpolate_table_increment(this, temperature, increment)
        class(chemical_kinetics_solver), intent(in) :: this
        real(dp), intent(in) :: temperature
        real(dp), dimension(:), intent(out) :: increment

        integer :: lower, upper, middle
        real(dp) :: weight

        if (.not. allocated(this%table_temperature) .or. &
            .not. allocated(this%table_concentration_increment)) then
            error stop 'Chemical kinetics: table solver is not configured'
        end if
        if (size(increment) /= this%species_number) then
            error stop 'Chemical kinetics: table interpolation size mismatch'
        end if

        if (temperature < this%table_temperature(1)) then
            increment = 0.0_dp
            return
        end if
        if (temperature >= this%table_temperature( &
            size(this%table_temperature))) then
            increment = this%table_concentration_increment(:, &
                size(this%table_temperature))
            return
        end if

        lower = 1
        upper = size(this%table_temperature)
        do while (upper-lower > 1)
            middle = (lower+upper)/2
            if (temperature >= this%table_temperature(middle)) then
                lower = middle
            else
                upper = middle
            end if
        end do

        weight = (temperature-this%table_temperature(lower))/ &
            (this%table_temperature(upper)- &
            this%table_temperature(lower))
        increment = (1.0_dp-weight)* &
            this%table_concentration_increment(:,lower) + &
            weight*this%table_concentration_increment(:,upper)
    end subroutine interpolate_table_increment


    subroutine kinetics_rhs(n, time, concentration, concentration_rate)
        integer, intent(in) :: n
        real(dp), intent(in) :: time
        real(dp), dimension(*), intent(in) :: concentration
        real(dp), dimension(*), intent(out) :: concentration_rate

        integer :: reaction, specie, component, multiplicity, index
        integer :: specie_index, stoichiometry
        real(dp) :: forward_rate, reverse_rate, net_rate
        real(dp) :: forward_constant, reverse_constant
        real(dp) :: third_body, total_concentration, positive_concentration

        if (.not. associated(thread_workspace%chemistry)) then
            error stop 'Chemical kinetics RHS: thread workspace is not initialized'
        end if
        if (n /= thread_workspace%species_number) then
            error stop 'Chemical kinetics RHS: species-count mismatch'
        end if

        total_concentration = 0.0_dp
        do specie = 1, n
            concentration_rate(specie) = 0.0_dp
            if (thread_workspace%has_any_third_body_reaction) then
                total_concentration = total_concentration + &
                    max(concentration(specie),0.0_dp)
            end if
        end do

        do reaction = 1, thread_workspace%reactions_number
            third_body = 0.0_dp
            if (thread_workspace%reaction_uses_third_body(reaction)) then
                third_body = thread_workspace%default_third_body_efficiency* &
                    total_concentration
                do index = thread_workspace%third_body_offset(reaction-1), &
                        thread_workspace%third_body_offset(reaction)-1
                    specie_index = &
                        thread_workspace%third_body_species(index)
                    third_body = third_body + &
                        thread_workspace%third_body_efficiency_delta(index)* &
                        max(concentration(specie_index),0.0_dp)
                end do
            end if

            call effective_rate_constants(reaction,third_body, &
                forward_constant,reverse_constant)

            forward_rate = forward_constant
            do component = 1, thread_workspace%reactant_count(reaction)
                specie_index = &
                    thread_workspace%reactant_species(component,reaction)
                positive_concentration = max(concentration(specie_index),0.0_dp)
                do multiplicity = 1, &
                        thread_workspace%reactant_multiplicity( &
                        component,reaction)
                    forward_rate = forward_rate*positive_concentration
                end do
            end do
            do multiplicity = 1, &
                    thread_workspace%forward_third_body_power(reaction)
                forward_rate = forward_rate*third_body
            end do

            reverse_rate = reverse_constant
            do component = 1, thread_workspace%product_count(reaction)
                specie_index = &
                    thread_workspace%product_species(component,reaction)
                positive_concentration = max(concentration(specie_index),0.0_dp)
                do multiplicity = 1, &
                        thread_workspace%product_multiplicity( &
                        component,reaction)
                    reverse_rate = reverse_rate*positive_concentration
                end do
            end do
            do multiplicity = 1, &
                    thread_workspace%reverse_third_body_power(reaction)
                reverse_rate = reverse_rate*third_body
            end do

            net_rate = forward_rate-reverse_rate
#ifdef CHEMISTRY_STRICT_DIAGNOSTICS
            if (.not. ieee_is_finite(net_rate)) then
                error stop 'Chemical kinetics RHS: non-finite reaction rate'
            end if
#endif

            do component = 1, thread_workspace%net_count(reaction)
                specie_index = thread_workspace%net_species(component,reaction)
                stoichiometry = &
                    thread_workspace%net_stoichiometry(component,reaction)
                concentration_rate(specie_index) = &
                    concentration_rate(specie_index) + &
                    real(stoichiometry,dp)*net_rate
            end do
        end do

        if (time < -huge(1.0_dp)) concentration_rate(1) = &
            concentration_rate(1)
    end subroutine kinetics_rhs


    subroutine effective_rate_constants(reaction, third_body, &
            forward_constant, reverse_constant)
        integer, intent(in) :: reaction
        real(dp), intent(in) :: third_body
        real(dp), intent(out) :: forward_constant, reverse_constant

        real(dp) :: high_rate, low_rate, reduced_pressure

        high_rate = thread_workspace%high_pressure_rate(reaction)
        low_rate = thread_workspace%low_pressure_rate(reaction)

        if (thread_workspace%reaction_uses_falloff(reaction)) then
            if (high_rate <= 0.0_dp .or. low_rate <= 0.0_dp .or. &
                third_body <= 0.0_dp) then
                forward_constant = 0.0_dp
            else
                reduced_pressure = low_rate*third_body/high_rate
                forward_constant = high_rate*reduced_pressure/ &
                    (1.0_dp+reduced_pressure)
                if (thread_workspace%reaction_uses_troe(reaction)) then
                    forward_constant = forward_constant* &
                        troe_falloff_factor(reaction,reduced_pressure)
                end if
            end if
        else
            forward_constant = high_rate
        end if

        if (thread_workspace%reaction_reversible(reaction)) then
            reverse_constant = forward_constant* &
                thread_workspace%reverse_factor(reaction)
        else
            reverse_constant = 0.0_dp
        end if
    end subroutine effective_rate_constants


    real(dp) function troe_falloff_factor(reaction, reduced_pressure) &
            result(factor)
        integer, intent(in) :: reaction
        real(dp), intent(in) :: reduced_pressure

        real(dp) :: f_center, c_troe, n_troe, d_troe, log_pressure
        real(dp) :: denominator, exponent

        if (reduced_pressure <= 0.0_dp) then
            factor = 1.0_dp
            return
        end if

        f_center = thread_workspace%troe_f_center(reaction)
        c_troe = thread_workspace%troe_c(reaction)
        n_troe = thread_workspace%troe_n(reaction)
        d_troe = 0.14_dp
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


#ifdef CHEMISTRY_PROFILE
    real(dp) function chemistry_wall_time() result(time_value)
#ifdef OMP
        use omp_lib, only: omp_get_wtime
        time_value = omp_get_wtime()
#else
        integer(int64) :: count, rate
        call system_clock(count,rate)
        if (rate > 0_int64) then
            time_value = real(count,dp)/real(rate,dp)
        else
            time_value = 0.0_dp
        end if
#endif
    end function chemistry_wall_time
#endif


    double precision function dummy_root(n, time, state, root_index)
        integer, intent(in) :: n, root_index
        double precision, intent(in) :: time
        double precision, dimension(*), intent(in) :: state

        dummy_root = time
        if (n < 0 .or. root_index < 0) dummy_root = state(1)
    end function dummy_root


    subroutine report_invalid_cell(message, i, j, k)
        character(len=*), intent(in) :: message
        integer, intent(in) :: i, j, k

!$omp critical(chemical_kinetics_error_output)
        write(error_unit,'(A,1X,A,1X,A,3(I0,1X))') &
            'Chemical kinetics error:',trim(message),'cell',i,j,k
!$omp end critical(chemical_kinetics_error_output)
        error stop 'Chemical kinetics solver failed'
    end subroutine report_invalid_cell


    subroutine report_mass_imbalance(i, j, k, source_sum, source_l1, &
            integrated_defect, integrated_tolerance)
        integer, intent(in) :: i, j, k
        real(dp), intent(in) :: source_sum, source_l1
        real(dp), intent(in) :: integrated_defect, integrated_tolerance

!$omp critical(chemical_kinetics_error_output)
        write(error_unit,'(A,3(I0,1X))') &
            'Chemical kinetics mass imbalance at cell ',i,j,k
        write(error_unit,'(A,ES24.16)') &
            '  sum(species_source) [kg m-3 s-1] = ',source_sum
        write(error_unit,'(A,ES24.16)') &
            '  sum(abs(species_source))          = ',source_l1
        write(error_unit,'(A,ES24.16)') &
            '  integrated mass defect [kg m-3]   = ',integrated_defect
        write(error_unit,'(A,ES24.16)') &
            '  allowed integrated defect         = ',integrated_tolerance
!$omp end critical(chemical_kinetics_error_output)
        error stop 'Chemical kinetics solver failed'
    end subroutine report_mass_imbalance


    subroutine report_slatec_failure(i, j, k, nstate, ierflg)
        integer, intent(in) :: i, j, k, nstate, ierflg

!$omp critical(chemical_kinetics_error_output)
        write(error_unit,'(A,3(I0,1X),A,I0,A,I0)') &
            'SLATEC chemistry failure at cell ',i,j,k, &
            ' NSTATE=',nstate,' IERFLG=',ierflg
!$omp end critical(chemical_kinetics_error_output)
        error stop 'Chemical kinetics SLATEC integration failed'
    end subroutine report_slatec_failure


    subroutine zero_vector_field(field)
        type(field_vector_cons), pointer, intent(inout) :: field

        integer :: component

        do component = 1, size(field%pr)
            field%pr(component)%cells = 0.0_dp
        end do
    end subroutine zero_vector_field

end module chemical_kinetics_solver_class
