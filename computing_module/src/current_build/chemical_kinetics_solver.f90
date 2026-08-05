module chemical_kinetics_solver_class

    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use, intrinsic :: iso_fortran_env, only: error_unit

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

    ! DDRIV3 accepts an external right-hand-side callback and provides no user
    ! context argument.  The only module state retained by the refactored solver
    ! is therefore one private workspace per OpenMP thread.  Registered fields,
    ! configuration and table data remain owned by each solver instance.
    type :: kinetics_thread_workspace
        integer :: species_number = 0
        integer :: reactions_number = 0
        type(chemical_properties), pointer :: chemistry => null()
        real(dp) :: temperature = 0.0_dp
        real(dp), allocatable :: concentration_initial(:)
        real(dp), allocatable :: concentration_final(:)
        real(dp), allocatable :: third_body_concentration(:)
        real(dp), allocatable :: high_pressure_rate(:)
        real(dp), allocatable :: low_pressure_rate(:)
        real(dp), allocatable :: reverse_factor(:)
        real(dp), allocatable :: entropy(:)
        real(dp), allocatable :: enthalpy(:)
        real(dp), allocatable :: work(:)
        integer, allocatable :: iwork(:)
    end type kinetics_thread_workspace

    type(kinetics_thread_workspace), save :: thread_workspace
!$omp threadprivate(thread_workspace)

    interface
        subroutine ddriv3(n, t, y, f, nstate, tout, ntask, nroot, eps, &
            ewt, ierror, mint, miter, impl, ml, mu, mxord, hmax, work, &
            lenw, iwork, leniw, jacobn, fa, nde, mxstep, g, users, ierflg)
            external :: f, jacobn, fa, g, users
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

        real(dp), allocatable :: concentration_increment(:,:,:,:)
        real(dp), allocatable :: table_temperature(:)
        real(dp), allocatable :: table_concentration_increment(:,:)
    contains
        procedure :: solve_chemical_kinetics
        procedure :: write_chemical_kinetics_table
        procedure :: set_activation_temperature
        procedure :: set_slatec_controls
        procedure :: configure_table_approximated
        procedure :: use_detailed_kinetics
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
            activation_temperature, ode_solver, table_file)
        type(data_manager), intent(inout) :: manager
        real(dp), intent(in), optional :: activation_temperature
        character(len=*), intent(in), optional :: ode_solver
        character(len=*), intent(in), optional :: table_file

        type(field_scalar_cons_pointer) :: scalar_pointer
        type(field_vector_cons_pointer) :: vector_pointer
        type(field_tensor_cons_pointer) :: tensor_pointer
        integer, dimension(3,2) :: allocation_bounds
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

        allocation_bounds = manager%domain%get_local_utter_cells_bounds()
        allocate(constructor%concentration_increment( &
            constructor%species_number, &
            allocation_bounds(1,1):allocation_bounds(1,2), &
            allocation_bounds(2,1):allocation_bounds(2,2), &
            allocation_bounds(3,1):allocation_bounds(3,2)))
        constructor%concentration_increment = 0.0_dp

        do specie = 1, constructor%species_number
            if (constructor%thermophysics%thermo_ptr%molar_masses(specie) <= &
                0.0_dp) then
                error stop 'Chemical kinetics: non-positive species molar mass'
            end if
        end do

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


    subroutine solve_chemical_kinetics(this, time_step)
        class(chemical_kinetics_solver), intent(inout) :: this
        real(dp), intent(in) :: time_step

        integer, dimension(3,2) :: cell_loop
        integer :: i, j, k, specie
        real(dp) :: density_cell, temperature_cell, energy_source_cell
        real(dp), allocatable :: mass_fraction_cell(:)
        real(dp), allocatable :: species_source_cell(:)
        real(dp), allocatable :: concentration_increment_cell(:)
        real(dp) :: mass_fraction_sum, negative_tolerance

        if (.not. ieee_is_finite(time_step) .or. time_step <= 0.0_dp) then
            error stop 'Chemical kinetics: time step must be finite and positive'
        end if

        call this%validate_configuration()
        cell_loop = this%domain%get_local_inner_cells_bounds()

        this%energy_source%s_ptr%cells = 0.0_dp
        call zero_vector_field(this%species_source%v_ptr)
        this%concentration_increment = 0.0_dp

        associate( &
            temperature_field => this%temperature%s_ptr, &
            density_field => this%density%s_ptr, &
            mass_fraction_field => this%mass_fraction%v_ptr, &
            energy_source_field => this%energy_source%s_ptr, &
            species_source_field => this%species_source%v_ptr, &
            markers => this%boundary%bc_ptr%bc_markers)

!$omp parallel default(shared) &
!$omp private(i,j,k,specie,density_cell,temperature_cell,energy_source_cell) &
!$omp private(mass_fraction_cell,species_source_cell) &
!$omp private(concentration_increment_cell,mass_fraction_sum) &
!$omp private(negative_tolerance)
        allocate(mass_fraction_cell(this%species_number))
        allocate(species_source_cell(this%species_number))
        allocate(concentration_increment_cell(this%species_number))
        if (this%ode_solver == 'slatec') call this%ensure_thread_workspace()

!$omp do collapse(3) schedule(dynamic)
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

                    mass_fraction_sum = 0.0_dp
                    negative_tolerance = &
                        this%negative_concentration_tolerance
                    do specie = 1, this%species_number
                        mass_fraction_cell(specie) = &
                            mass_fraction_field%pr(specie)%cells(i,j,k)
                        if (.not. ieee_is_finite(mass_fraction_cell(specie))) then
                            call report_invalid_cell( &
                                'non-finite mass fraction',i,j,k)
                        end if
                        if (mass_fraction_cell(specie) < -negative_tolerance) then
                            call report_invalid_cell('negative mass fraction',i,j,k)
                        end if
                        mass_fraction_cell(specie) = &
                            max(mass_fraction_cell(specie),0.0_dp)
                        mass_fraction_sum = mass_fraction_sum + &
                            mass_fraction_cell(specie)
                    end do
                    if (mass_fraction_sum <= tiny(1.0_dp)) then
                        call report_invalid_cell('empty composition',i,j,k)
                    end if
                    mass_fraction_cell = &
                        mass_fraction_cell/mass_fraction_sum

                    select case (this%ode_solver)
                    case ('slatec')
                        call this%solve_cell_detailed_kinetics( &
                            density_cell,temperature_cell,mass_fraction_cell, &
                            time_step,i,j,k,concentration_increment_cell)
                    case ('table_approximated')
                        call this%interpolate_table_increment( &
                            temperature_cell,concentration_increment_cell)
                    case default
                        call report_invalid_cell('unsupported ODE solver',i,j,k)
                    end select

                    call this%assemble_cell_sources( &
                        density_cell,mass_fraction_cell,time_step, &
                        concentration_increment_cell,species_source_cell, &
                        energy_source_cell,i,j,k)

                    energy_source_field%cells(i,j,k) = energy_source_cell
                    do specie = 1, this%species_number
                        species_source_field%pr(specie)%cells(i,j,k) = &
                            species_source_cell(specie)
                        this%concentration_increment(specie,i,j,k) = &
                            concentration_increment_cell(specie)
                    end do
                end do
            end do
        end do
!$omp end do

        deallocate(mass_fraction_cell)
        deallocate(species_source_cell)
        deallocate(concentration_increment_cell)
!$omp end parallel

        end associate
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

        allocate(thread_workspace%concentration_initial(this%species_number))
        allocate(thread_workspace%concentration_final(this%species_number))
        allocate(thread_workspace%third_body_concentration( &
            this%reactions_number))
        allocate(thread_workspace%high_pressure_rate(this%reactions_number))
        allocate(thread_workspace%low_pressure_rate(this%reactions_number))
        allocate(thread_workspace%reverse_factor(this%reactions_number))
        allocate(thread_workspace%entropy(this%species_number))
        allocate(thread_workspace%enthalpy(this%species_number))

        work_size = this%species_number**2 + 10*this%species_number + 250
        allocate(thread_workspace%work(work_size))
        allocate(thread_workspace%iwork(50 + this%species_number))

        thread_workspace%concentration_initial = 0.0_dp
        thread_workspace%concentration_final = 0.0_dp
        thread_workspace%third_body_concentration = 0.0_dp
        thread_workspace%high_pressure_rate = 0.0_dp
        thread_workspace%low_pressure_rate = 0.0_dp
        thread_workspace%reverse_factor = 0.0_dp
        thread_workspace%entropy = 0.0_dp
        thread_workspace%enthalpy = 0.0_dp
        thread_workspace%work = 0.0_dp
        thread_workspace%iwork = 0
    end subroutine ensure_thread_workspace


    subroutine clear_thread_workspace()
        if (allocated(thread_workspace%concentration_initial)) &
            deallocate(thread_workspace%concentration_initial)
        if (allocated(thread_workspace%concentration_final)) &
            deallocate(thread_workspace%concentration_final)
        if (allocated(thread_workspace%third_body_concentration)) &
            deallocate(thread_workspace%third_body_concentration)
        if (allocated(thread_workspace%high_pressure_rate)) &
            deallocate(thread_workspace%high_pressure_rate)
        if (allocated(thread_workspace%low_pressure_rate)) &
            deallocate(thread_workspace%low_pressure_rate)
        if (allocated(thread_workspace%reverse_factor)) &
            deallocate(thread_workspace%reverse_factor)
        if (allocated(thread_workspace%entropy)) &
            deallocate(thread_workspace%entropy)
        if (allocated(thread_workspace%enthalpy)) &
            deallocate(thread_workspace%enthalpy)
        if (allocated(thread_workspace%work)) &
            deallocate(thread_workspace%work)
        if (allocated(thread_workspace%iwork)) &
            deallocate(thread_workspace%iwork)
        nullify(thread_workspace%chemistry)
        thread_workspace%species_number = 0
        thread_workspace%reactions_number = 0
        thread_workspace%temperature = 0.0_dp
    end subroutine clear_thread_workspace


    subroutine solve_cell_detailed_kinetics(this, density, temperature, &
            mass_fraction, time_step, i_cell, j_cell, k_cell, &
            concentration_increment)
        class(chemical_kinetics_solver), intent(in) :: this
        real(dp), intent(in) :: density, temperature, time_step
        real(dp), dimension(:), intent(in) :: mass_fraction
        integer, intent(in) :: i_cell, j_cell, k_cell
        real(dp), dimension(:), intent(out) :: concentration_increment

        integer :: n, nstate, ntask, nroot, ierror, mint, miter, impl
        integer :: ml, mu, mxord, lenw, leniw, nde, ierflg
        real(dp) :: eps, hmax, time_in, time_out
        real(dp), dimension(1) :: error_weight
        integer :: specie
        real(dp) :: concentration_scale, negative_limit

        call this%ensure_thread_workspace()

        do specie = 1, this%species_number
            thread_workspace%concentration_initial(specie) = &
                density*mass_fraction(specie)/ &
                this%thermophysics%thermo_ptr%molar_masses(specie)
        end do
        thread_workspace%concentration_final = &
            thread_workspace%concentration_initial
        thread_workspace%temperature = min(temperature,maximum_rate_temperature)

        call this%prepare_cell_rate_coefficients(temperature)

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
        thread_workspace%work = 0.0_dp
        thread_workspace%iwork = 0

        call ddriv3(n,time_in,thread_workspace%concentration_final, &
            kinetics_rhs,nstate,time_out,ntask,nroot,eps,error_weight, &
            ierror,mint,miter,impl,ml,mu,mxord,hmax, &
            thread_workspace%work,lenw,thread_workspace%iwork,leniw, &
            kinetics_rhs,kinetics_rhs,nde,this%slatec_max_steps, &
            dummy_root,kinetics_rhs,ierflg)

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

        integer :: reaction, component, specie, reaction_type
        integer :: left_count, right_count, specie_index
        real(dp) :: temperature, entropy_change, enthalpy_change
        real(dp) :: mole_change, equilibrium_pressure, equilibrium_concentration
        real(dp) :: exponent_argument

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
            left_count = chemistry%chemical_coeffs(1,reaction,1)
            right_count = chemistry%chemical_coeffs(1,reaction,2)

            do component = 2, left_count + 1
                specie_index = chemistry%chemical_coeffs(component,reaction,1)
                if (specie_index < 1 .or. &
                    specie_index > this%species_number) cycle
                entropy_change = entropy_change - &
                    thread_workspace%entropy(specie_index)
                enthalpy_change = enthalpy_change - &
                    thread_workspace%enthalpy(specie_index)
                mole_change = mole_change - 1.0_dp
            end do
            do component = 2, right_count + 1
                specie_index = chemistry%chemical_coeffs(component,reaction,2)
                if (specie_index < 1 .or. &
                    specie_index > this%species_number) cycle
                entropy_change = entropy_change + &
                    thread_workspace%entropy(specie_index)
                enthalpy_change = enthalpy_change + &
                    thread_workspace%enthalpy(specie_index)
                mole_change = mole_change + 1.0_dp
            end do

            thread_workspace%reverse_factor(reaction) = 0.0_dp
            reaction_type = chemistry%reactions_type(reaction)
            if (reaction_type == 0 .or. reaction_type == 1 .or. &
                reaction_type == 2 .or. reaction_type == 3 .or. &
                reaction_type == 4) then
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

        integer :: specie
        real(dp) :: source_sum, source_scale, source_tolerance
        real(dp) :: reference_enthalpy_mass

        do specie = 1, this%species_number
            species_source(specie) = concentration_increment(specie)* &
                this%thermophysics%thermo_ptr%molar_masses(specie)/ &
                time_step
        end do

        source_sum = sum(species_source)
        source_scale = max(sum(abs(species_source)), &
            density/time_step*100.0_dp*epsilon(1.0_dp))
        source_tolerance = this%mass_balance_tolerance*source_scale
        if (abs(source_sum) > source_tolerance) then
            call report_invalid_cell( &
                'chemical mass imbalance exceeds tolerance', &
                i_cell,j_cell,k_cell)
        end if

        ! Remove only the integrator-level roundoff residual.  The correction is
        ! distributed with the normalized beginning-of-step composition so the
        ! conservative chemistry source satisfies sum_k S_{rho Y_k}=0 exactly.
        species_source = species_source - mass_fraction*source_sum

        energy_source = 0.0_dp
        do specie = 1, this%species_number
            reference_enthalpy_mass = &
                this%thermophysics%thermo_ptr%specie_enthalpy_molar( &
                    T_ref,specie)/ &
                this%thermophysics%thermo_ptr%molar_masses(specie)
            energy_source = energy_source - &
                reference_enthalpy_mass*species_source(specie)
            concentration_increment(specie) = species_source(specie)* &
                time_step/ &
                this%thermophysics%thermo_ptr%molar_masses(specie)
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

        integer :: reaction, specie, component, specie_index
        integer :: left_count, right_count
        real(dp) :: forward_rate, reverse_rate, net_rate
        real(dp) :: forward_constant, reverse_constant
        real(dp) :: third_body

        if (.not. associated(thread_workspace%chemistry)) then
            error stop 'Chemical kinetics RHS: thread workspace is not initialized'
        end if
        if (n /= thread_workspace%species_number) then
            error stop 'Chemical kinetics RHS: species-count mismatch'
        end if

        do specie = 1, n
            concentration_rate(specie) = 0.0_dp
        end do

        associate(chemistry => thread_workspace%chemistry)
        do reaction = 1, thread_workspace%reactions_number
            third_body = 0.0_dp
            do specie = 1, n
                third_body = third_body + max(concentration(specie),0.0_dp)* &
                    chemistry%enhanced_efficiencies(reaction,specie)
            end do
            thread_workspace%third_body_concentration(reaction) = third_body

            call effective_rate_constants(reaction,third_body, &
                forward_constant,reverse_constant)

            forward_rate = forward_constant
            left_count = chemistry%chemical_coeffs(1,reaction,1)
            do component = 2, left_count + 1
                specie_index = chemistry%chemical_coeffs(component,reaction,1)
                if (specie_index == n+1) then
                    forward_rate = forward_rate*third_body
                else if (specie_index >= 1 .and. specie_index <= n) then
                    forward_rate = forward_rate* &
                        max(concentration(specie_index),0.0_dp)
                end if
            end do

            reverse_rate = reverse_constant
            right_count = chemistry%chemical_coeffs(1,reaction,2)
            do component = 2, right_count + 1
                specie_index = chemistry%chemical_coeffs(component,reaction,2)
                if (specie_index == n+1) then
                    reverse_rate = reverse_rate*third_body
                else if (specie_index >= 1 .and. specie_index <= n) then
                    reverse_rate = reverse_rate* &
                        max(concentration(specie_index),0.0_dp)
                end if
            end do

            net_rate = forward_rate-reverse_rate
            if (.not. ieee_is_finite(net_rate)) then
                error stop 'Chemical kinetics RHS: non-finite reaction rate'
            end if

            do component = 2, left_count + 1
                specie_index = chemistry%chemical_coeffs(component,reaction,1)
                if (specie_index >= 1 .and. specie_index <= n) then
                    concentration_rate(specie_index) = &
                        concentration_rate(specie_index)-net_rate
                end if
            end do
            do component = 2, right_count + 1
                specie_index = chemistry%chemical_coeffs(component,reaction,2)
                if (specie_index >= 1 .and. specie_index <= n) then
                    concentration_rate(specie_index) = &
                        concentration_rate(specie_index)+net_rate
                end if
            end do
        end do
        end associate

        if (time < -huge(1.0_dp)) concentration_rate(1) = &
            concentration_rate(1)
    end subroutine kinetics_rhs


    subroutine effective_rate_constants(reaction, third_body, &
            forward_constant, reverse_constant)
        integer, intent(in) :: reaction
        real(dp), intent(in) :: third_body
        real(dp), intent(out) :: forward_constant, reverse_constant

        integer :: reaction_type
        real(dp) :: high_rate, low_rate, reduced_pressure
        real(dp) :: falloff_factor

        reaction_type = thread_workspace%chemistry%reactions_type(reaction)
        high_rate = thread_workspace%high_pressure_rate(reaction)
        low_rate = thread_workspace%low_pressure_rate(reaction)

        select case (reaction_type)
        case (2,3,7)
            if (high_rate <= 0.0_dp .or. low_rate <= 0.0_dp .or. &
                third_body <= 0.0_dp) then
                forward_constant = 0.0_dp
            else
                reduced_pressure = low_rate*third_body/high_rate
                forward_constant = high_rate*reduced_pressure/ &
                    (1.0_dp+reduced_pressure)
                if (reaction_type == 3 .or. &
                    (reaction_type == 7 .and. has_troe_parameters(reaction))) then
                    falloff_factor = troe_falloff_factor( &
                        reaction,reduced_pressure)
                    forward_constant = forward_constant*falloff_factor
                end if
            end if
        case default
            forward_constant = high_rate
        end select

        select case (reaction_type)
        case (0,1,2,3,4)
            reverse_constant = forward_constant* &
                thread_workspace%reverse_factor(reaction)
        case default
            reverse_constant = 0.0_dp
        end select
    end subroutine effective_rate_constants


    logical function has_troe_parameters(reaction) result(has_troe)
        integer, intent(in) :: reaction

        has_troe = maxval(abs( &
            thread_workspace%chemistry%Troe_coeffs(reaction,:))) > 0.0_dp
    end function has_troe_parameters


    real(dp) function troe_falloff_factor(reaction, reduced_pressure) &
            result(factor)
        integer, intent(in) :: reaction
        real(dp), intent(in) :: reduced_pressure

        real(dp) :: alpha, temperature_1, temperature_2, temperature_3
        real(dp) :: f_center, c_troe, n_troe, d_troe, log_pressure
        real(dp) :: denominator, exponent

        alpha = thread_workspace%chemistry%Troe_coeffs(reaction,1)
        temperature_1 = thread_workspace%chemistry%Troe_coeffs(reaction,2)
        temperature_2 = thread_workspace%chemistry%Troe_coeffs(reaction,3)
        temperature_3 = thread_workspace%chemistry%Troe_coeffs(reaction,4)

        f_center = 0.0_dp
        if (temperature_1 > 0.0_dp) then
            f_center = f_center + (1.0_dp-alpha)* &
                exp(-thread_workspace%temperature/temperature_1)
        end if
        if (temperature_2 > 0.0_dp) then
            f_center = f_center + alpha* &
                exp(-thread_workspace%temperature/temperature_2)
        end if
        if (temperature_3 > 0.0_dp) then
            f_center = f_center + &
                exp(-temperature_3/thread_workspace%temperature)
        end if
        f_center = min(max(f_center,tiny(1.0_dp)),1.0_dp)

        if (reduced_pressure <= 0.0_dp) then
            factor = 1.0_dp
            return
        end if

        c_troe = -0.4_dp-0.67_dp*log10(f_center)
        n_troe = 0.75_dp-1.27_dp*log10(f_center)
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
