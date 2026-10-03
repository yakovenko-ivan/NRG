module chemical_kinetics_solver_class

    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use, intrinsic :: iso_fortran_env, only: error_unit, output_unit, int64
#ifdef NRG_ENABLE_CVODE
    use, intrinsic :: iso_c_binding, only: c_int, c_long, c_int64_t, &
        c_double, c_ptr, c_funptr, c_funloc, c_associated, c_null_ptr, &
        c_loc, c_f_pointer
    use fsundials_core_mod
    use fcvode_mod
    use fnvector_serial_mod
    use fsunmatrix_dense_mod
    use fsunlinsol_dense_mod
#endif

    use kind_parameters, only: dp
    use global_data, only: r_gase_J, P_atm, T_ref, task_setup_folder, &
        fold_sep, chemical_mechanisms_folder
    use field_pointers
    use boundary_conditions_class
    use data_manager_class
    use computational_domain_class
    use thermophysical_properties_class
    use chemical_properties_class
    use chemical_kinetics_core_class, only: chemical_kinetics_core, &
        chemical_kinetics_core_c, chemical_rate_state
#ifdef OMP
    use omp_lib, only: omp_get_thread_num, omp_get_max_threads
#endif

    implicit none

    private
    public :: chemical_kinetics_solver, chemical_kinetics_solver_c

    real(dp), parameter :: default_activation_temperature = 305.0_dp
    real(dp), parameter :: default_slatec_accuracy = 1.0e-9_dp
    real(dp), parameter :: default_slatec_error_weight = 1.0e-4_dp
    real(dp), parameter :: default_slatec_max_internal_step = 1.0e-7_dp
    integer, parameter :: default_slatec_max_steps = 10000
    real(dp), parameter :: default_cvode_relative_tolerance = 1.0e-8_dp
    real(dp), parameter :: default_cvode_absolute_tolerance = 1.0e-12_dp
    integer, parameter :: default_cvode_max_steps = 10000
    ! KIN 4.04.02 / K.INP-inspired QSS1 defaults.
    real(dp), parameter :: default_qss1_relative_change_limit = 5.0e-2_dp
    real(dp), parameter :: default_qss1_minimum_internal_step = 1.0e-9_dp
    real(dp), parameter :: default_qss1_active_concentration_fraction = &
        1.0e-7_dp
    real(dp), parameter :: default_qss1_step_growth_factor = 1.04_dp
    integer, parameter :: default_qss1_max_steps = 100000
    real(dp), parameter :: qss1_concentration_floor = 1.0e-30_dp
    real(dp), parameter :: qss1_small_loss_argument = 1.0e-4_dp
    real(dp), parameter :: qss1_conservation_rank_tolerance = 1.0e-12_dp
    real(dp), parameter :: qss1_conservation_check_tolerance = 1.0e-10_dp
    real(dp), parameter :: qss1_projection_negative_factor = 4096.0_dp
    real(dp), parameter :: default_negative_concentration_tolerance = 1.0e-10_dp
    real(dp), parameter :: default_mass_balance_tolerance = 1.0e-6_dp
#ifdef CHEMISTRY_PROFILE
    ! Periodic profiling output cadence. Change to 500 for quieter long runs.
    integer, parameter :: chemistry_profile_output_interval = 100
#endif
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
        type(chemical_kinetics_core) :: kinetics_core
        type(chemical_rate_state) :: rate_state
        real(dp) :: temperature = 0.0_dp
        real(dp) :: default_third_body_efficiency = 0.0_dp
        logical :: has_any_third_body_reaction = .false.

        real(dp), allocatable :: concentration_initial(:)
        real(dp), allocatable :: concentration_final(:)
        real(dp), allocatable :: qss_production(:)
        real(dp), allocatable :: qss_destruction(:)
        real(dp), allocatable :: qss_trial(:)
        real(dp), allocatable :: qss_projected(:)
        real(dp), allocatable :: qss_projection_weights(:)
        real(dp), allocatable :: qss_projection_target(:)
        real(dp), allocatable :: qss_projection_residual(:)
        real(dp), allocatable :: qss_projection_lambda(:)
        real(dp), allocatable :: qss_projection_factor(:,:)
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

#ifdef NRG_ENABLE_CVODE
    ! CVODE state is kept explicitly per OpenMP worker rather than in
    ! THREADPRIVATE storage.  This follows the SUNDIALS multi-threading model:
    ! one independent SUNContext and one complete solver object graph per thread.
    type, bind(C) :: cvode_user_context
        integer(c_int) :: worker_index = 0_c_int
    end type cvode_user_context

    type :: cvode_worker_workspace
        integer :: species_number = 0
        integer :: reactions_number = 0
        type(chemical_kinetics_core) :: kinetics_core
        type(chemical_rate_state) :: rate_state
        real(dp), allocatable :: concentration_initial(:)
        real(dp), allocatable :: concentration_final(:)
        type(c_ptr) :: cvode_mem = c_null_ptr
        type(c_ptr) :: sunctx = c_null_ptr
        type(N_Vector), pointer :: cvode_y => null()
        real(c_double), pointer :: cvode_y_data(:) => null()
        type(SUNMatrix), pointer :: cvode_matrix => null()
        type(SUNLinearSolver), pointer :: cvode_linear_solver => null()
        logical :: initialized = .false.
    end type cvode_worker_workspace

    type(cvode_worker_workspace), allocatable, save :: cvode_workers(:)
    type(cvode_user_context), allocatable, target, save :: cvode_user_contexts(:)
#endif

#ifdef NRG_ENABLE_CVODE
    ! Direct interfaces to the native SUNDIALS C API used in the concurrent
    ! chemistry hot path.  The generated F2003/SWIG wrappers remain useful for
    ! one-time object construction, but are intentionally bypassed for calls
    ! made concurrently by OpenMP workers.
    interface
        integer(c_int) function nrg_cvode_init_native( &
                cvode_mem, rhs_function, t0, y0) bind(C,name='CVodeInit')
            import :: c_int, c_double, c_ptr, c_funptr
            type(c_ptr), value :: cvode_mem
            type(c_funptr), value :: rhs_function
            real(c_double), value :: t0
            type(c_ptr), value :: y0
        end function nrg_cvode_init_native

        integer(c_int) function nrg_cvode_set_user_data_native( &
                cvode_mem, user_data) bind(C,name='CVodeSetUserData')
            import :: c_int, c_ptr
            type(c_ptr), value :: cvode_mem
            type(c_ptr), value :: user_data
        end function nrg_cvode_set_user_data_native

        integer(c_int) function nrg_cvode_reinit_native( &
                cvode_mem, t0, y0) bind(C,name='CVodeReInit')
            import :: c_int, c_double, c_ptr
            type(c_ptr), value :: cvode_mem
            real(c_double), value :: t0
            type(c_ptr), value :: y0
        end function nrg_cvode_reinit_native

        integer(c_int) function nrg_cvode_step_native( &
                cvode_mem, tout, yout, tret, itask) bind(C,name='CVode')
            import :: c_int, c_double, c_ptr
            type(c_ptr), value :: cvode_mem
            real(c_double), value :: tout
            type(c_ptr), value :: yout
            real(c_double), intent(out) :: tret
            integer(c_int), value :: itask
        end function nrg_cvode_step_native

        integer(c_int) function nrg_cvode_get_num_steps_native( &
                cvode_mem, nsteps) bind(C,name='CVodeGetNumSteps')
            import :: c_int, c_long, c_ptr
            type(c_ptr), value :: cvode_mem
            integer(c_long), intent(out) :: nsteps
        end function nrg_cvode_get_num_steps_native

        integer(c_int) function nrg_cvode_get_num_rhs_evals_native( &
                cvode_mem, nfe) bind(C,name='CVodeGetNumRhsEvals')
            import :: c_int, c_long, c_ptr
            type(c_ptr), value :: cvode_mem
            integer(c_long), intent(out) :: nfe
        end function nrg_cvode_get_num_rhs_evals_native

        integer(c_int) function nrg_cvode_get_num_lin_rhs_evals_native( &
                cvode_mem, nfe_ls) bind(C,name='CVodeGetNumLinRhsEvals')
            import :: c_int, c_long, c_ptr
            type(c_ptr), value :: cvode_mem
            integer(c_long), intent(out) :: nfe_ls
        end function nrg_cvode_get_num_lin_rhs_evals_native

        integer(c_int) function nrg_cvode_get_num_jac_evals_native( &
                cvode_mem, nje) bind(C,name='CVodeGetNumJacEvals')
            import :: c_int, c_long, c_ptr
            type(c_ptr), value :: cvode_mem
            integer(c_long), intent(out) :: nje
        end function nrg_cvode_get_num_jac_evals_native

        function nrg_nvector_get_array_pointer_native(vector) &
                result(data_pointer) bind(C,name='N_VGetArrayPointer')
            import :: c_ptr
            type(c_ptr), value :: vector
            type(c_ptr) :: data_pointer
        end function nrg_nvector_get_array_pointer_native
    end interface
#endif

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
        type(chemical_kinetics_core) :: kinetics_core
        type(boundary_conditions_pointer) :: boundary

        integer :: species_number = 0
        integer :: reactions_number = 0

        real(dp), allocatable :: inverse_molar_mass(:)
        real(dp), allocatable :: reference_enthalpy_mass(:)
        integer :: qss_stoichiometric_rank = 0
        integer :: qss_conservation_count = 0
        integer :: qss_projection_constraint_count = 0
        logical :: qss_has_independent_mass_constraint = .false.
        real(dp), allocatable :: qss_conservation_basis(:,:)
        real(dp), allocatable :: qss_projection_basis(:,:)
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

        character(len=20) :: ode_solver = 'qss1'
        real(dp) :: activation_temperature = default_activation_temperature
        real(dp) :: slatec_accuracy = default_slatec_accuracy
        real(dp) :: slatec_error_weight = default_slatec_error_weight
        real(dp) :: slatec_max_internal_step = &
            default_slatec_max_internal_step
        integer :: slatec_max_steps = default_slatec_max_steps
        real(dp) :: cvode_relative_tolerance = default_cvode_relative_tolerance
        real(dp) :: cvode_absolute_tolerance = default_cvode_absolute_tolerance
        integer :: cvode_max_steps = default_cvode_max_steps
        real(dp) :: qss1_relative_change_limit = &
            default_qss1_relative_change_limit
        real(dp) :: qss1_minimum_internal_step = &
            default_qss1_minimum_internal_step
        real(dp) :: qss1_active_concentration_fraction = &
            default_qss1_active_concentration_fraction
        real(dp) :: qss1_step_growth_factor = &
            default_qss1_step_growth_factor
        integer :: qss1_max_steps = default_qss1_max_steps
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
        integer(int64) :: total_qss_projection_clips = 0_int64
        integer(int64) :: total_qss_projection_clipped_components = 0_int64
        real(dp) :: maximum_qss_projection_clip_magnitude = 0.0_dp
        integer :: maximum_qss_projection_clip_species = 0
        integer(int64) :: last_active_cells = 0_int64
        integer(int64) :: last_integrator_calls = 0_int64
        integer(int64) :: last_internal_steps = 0_int64
        integer(int64) :: last_rhs_evaluations = 0_int64
        integer(int64) :: last_jacobian_evaluations = 0_int64
        integer(int64) :: last_qss_projection_clips = 0_int64
        real(dp) :: total_cvode_packing_time = 0.0_dp
        real(dp) :: total_rate_preparation_time = 0.0_dp
        real(dp) :: total_cvode_reinitialization_time = 0.0_dp
        real(dp) :: total_integration_time = 0.0_dp
        real(dp) :: total_qss_projection_time = 0.0_dp
        real(dp) :: maximum_qss_pre_projection_conservation_residual = 0.0_dp
        real(dp) :: maximum_qss_post_projection_conservation_residual = 0.0_dp
        real(dp) :: maximum_qss_post_source_conservation_residual = 0.0_dp
        real(dp) :: total_cvode_statistics_time = 0.0_dp
        real(dp) :: total_cvode_result_processing_time = 0.0_dp
        real(dp) :: total_source_assembly_time = 0.0_dp
        real(dp) :: total_solver_wall_time = 0.0_dp
        real(dp) :: last_cvode_packing_time = 0.0_dp
        real(dp) :: last_rate_preparation_time = 0.0_dp
        real(dp) :: last_cvode_reinitialization_time = 0.0_dp
        real(dp) :: last_integration_time = 0.0_dp
        real(dp) :: last_qss_projection_time = 0.0_dp
        real(dp) :: last_cvode_statistics_time = 0.0_dp
        real(dp) :: last_cvode_result_processing_time = 0.0_dp
        real(dp) :: last_source_assembly_time = 0.0_dp
        real(dp) :: last_solver_wall_time = 0.0_dp
    contains
        procedure :: solve_chemical_kinetics
        procedure :: write_chemical_kinetics_table
        procedure :: can_write_chemical_kinetics_table
        procedure :: set_activation_temperature
        procedure :: set_slatec_controls
        procedure :: set_cvode_controls
        procedure :: set_qss1_controls
        procedure :: use_cvode_kinetics
        procedure :: configure_table_approximated
        procedure :: use_detailed_kinetics
        procedure :: use_qss1_kinetics
        procedure :: set_concentration_increment_recording
        procedure :: reset_performance_statistics
        procedure :: write_performance_statistics
        procedure, private :: preprocess_mechanism
        procedure, private :: build_qss_conservation_basis
        procedure, private :: project_qss_final_state
        procedure, private :: qss_increment_conservation_residual
        procedure, private :: allocate_concentration_increment_storage
        procedure, private :: ensure_thread_workspace
        procedure, private :: prepare_cell_rate_coefficients
        procedure, private :: solve_cell_detailed_kinetics
        procedure, private :: solve_cell_qss1
#ifdef NRG_ENABLE_CVODE
        procedure, private :: ensure_cvode_workers
        procedure, private :: solve_cell_cvode_kinetics
#endif
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
        call constructor%build_qss_conservation_basis()

        if (present(activation_temperature)) then
            call constructor%set_activation_temperature(activation_temperature)
        end if

        ! Backend controls belong to solver_options so package interfaces can
        ! persist them into solver_data and every CFD solver sees the same
        ! chemistry configuration without backend-specific constructor edits.
        call constructor%set_slatec_controls( &
            accuracy=manager%solver_options%get_chemistry_slatec_accuracy(), &
            error_weight= &
                manager%solver_options%get_chemistry_slatec_error_weight(), &
            maximum_internal_step= &
                manager%solver_options%get_chemistry_slatec_max_internal_step(), &
            maximum_steps= &
                manager%solver_options%get_chemistry_slatec_max_steps())
        call constructor%set_cvode_controls( &
            relative_tolerance= &
                manager%solver_options%get_chemistry_cvode_relative_tolerance(), &
            absolute_tolerance= &
                manager%solver_options%get_chemistry_cvode_absolute_tolerance(), &
            maximum_steps= &
                manager%solver_options%get_chemistry_cvode_max_steps())
        call constructor%set_qss1_controls( &
            relative_change_limit= &
                manager%solver_options%get_chemistry_qss1_relative_change_limit(), &
            minimum_internal_step= &
                manager%solver_options%get_chemistry_qss1_minimum_internal_step(), &
            active_concentration_fraction=manager%solver_options% &
                get_chemistry_qss1_active_concentration_fraction(), &
            step_growth_factor= &
                manager%solver_options%get_chemistry_qss1_step_growth_factor(), &
            maximum_steps= &
                manager%solver_options%get_chemistry_qss1_max_steps())

        if (present(ode_solver)) then
            select case (trim(adjustl(ode_solver)))
            case ('slatec')
                call constructor%use_detailed_kinetics()
            case ('cvode')
                call constructor%use_cvode_kinetics()
            case ('qss1')
                call constructor%use_qss1_kinetics()
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
        else
            select case(trim(adjustl( &
                manager%solver_options%get_chemistry_backend())))
            case ('slatec')
                call constructor%use_detailed_kinetics()
            case ('cvode')
                call constructor%use_cvode_kinetics()
            case ('qss1')
                call constructor%use_qss1_kinetics()
            case default
                error stop 'Chemical kinetics: unsupported configured backend'
            end select
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
        if (.not. ieee_is_finite(this%cvode_relative_tolerance) .or. &
            this%cvode_relative_tolerance <= 0.0_dp) then
            error stop 'Chemical kinetics: invalid CVODE relative tolerance'
        end if
        if (.not. ieee_is_finite(this%cvode_absolute_tolerance) .or. &
            this%cvode_absolute_tolerance <= 0.0_dp) then
            error stop 'Chemical kinetics: invalid CVODE absolute tolerance'
        end if
        if (this%cvode_max_steps <= 0) then
            error stop 'Chemical kinetics: invalid CVODE maximum step count'
        end if
        if (.not. ieee_is_finite(this%qss1_relative_change_limit) .or. &
            this%qss1_relative_change_limit <= 0.0_dp) then
            error stop 'Chemical kinetics: invalid QSS1 change limit'
        end if
        if (.not. ieee_is_finite(this%qss1_minimum_internal_step) .or. &
            this%qss1_minimum_internal_step <= 0.0_dp) then
            error stop 'Chemical kinetics: invalid QSS1 minimum step'
        end if
        if (.not. ieee_is_finite( &
            this%qss1_active_concentration_fraction) .or. &
            this%qss1_active_concentration_fraction <= 0.0_dp) then
            error stop 'Chemical kinetics: invalid QSS1 active fraction'
        end if
        if (.not. ieee_is_finite(this%qss1_step_growth_factor) .or. &
            this%qss1_step_growth_factor < 1.0_dp) then
            error stop 'Chemical kinetics: invalid QSS1 growth factor'
        end if
        if (this%qss1_max_steps <= 0) then
            error stop 'Chemical kinetics: invalid QSS1 maximum step count'
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

    subroutine set_cvode_controls(this, relative_tolerance, &
            absolute_tolerance, maximum_steps)
        class(chemical_kinetics_solver), intent(inout) :: this
        real(dp), intent(in), optional :: relative_tolerance
        real(dp), intent(in), optional :: absolute_tolerance
        integer, intent(in), optional :: maximum_steps

        if (present(relative_tolerance)) then
            this%cvode_relative_tolerance = relative_tolerance
        end if
        if (present(absolute_tolerance)) then
            this%cvode_absolute_tolerance = absolute_tolerance
        end if
        if (present(maximum_steps)) this%cvode_max_steps = maximum_steps
        call this%validate_configuration()
    end subroutine set_cvode_controls


    subroutine set_qss1_controls(this, relative_change_limit, &
            minimum_internal_step, active_concentration_fraction, &
            step_growth_factor, maximum_steps)
        class(chemical_kinetics_solver), intent(inout) :: this
        real(dp), intent(in), optional :: relative_change_limit
        real(dp), intent(in), optional :: minimum_internal_step
        real(dp), intent(in), optional :: active_concentration_fraction
        real(dp), intent(in), optional :: step_growth_factor
        integer, intent(in), optional :: maximum_steps

        if (present(relative_change_limit)) then
            this%qss1_relative_change_limit = relative_change_limit
        end if
        if (present(minimum_internal_step)) then
            this%qss1_minimum_internal_step = minimum_internal_step
        end if
        if (present(active_concentration_fraction)) then
            this%qss1_active_concentration_fraction = &
                active_concentration_fraction
        end if
        if (present(step_growth_factor)) then
            this%qss1_step_growth_factor = step_growth_factor
        end if
        if (present(maximum_steps)) this%qss1_max_steps = maximum_steps
        call this%validate_configuration()
    end subroutine set_qss1_controls


    subroutine use_cvode_kinetics(this)
        class(chemical_kinetics_solver), intent(inout) :: this

#ifdef NRG_ENABLE_CVODE
        this%ode_solver = 'cvode'
#else
        error stop 'Chemical kinetics: CVODE backend requested, but NRG was '// &
            'built with NRG_ENABLE_CVODE=OFF'
#endif
    end subroutine use_cvode_kinetics


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


    subroutine use_qss1_kinetics(this)
        class(chemical_kinetics_solver), intent(inout) :: this

        this%ode_solver = 'qss1'
    end subroutine use_qss1_kinetics


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
        this%total_qss_projection_clips = 0_int64
        this%total_qss_projection_clipped_components = 0_int64
        this%maximum_qss_projection_clip_magnitude = 0.0_dp
        this%maximum_qss_projection_clip_species = 0
        this%last_active_cells = 0_int64
        this%last_integrator_calls = 0_int64
        this%last_internal_steps = 0_int64
        this%last_rhs_evaluations = 0_int64
        this%last_jacobian_evaluations = 0_int64
        this%last_qss_projection_clips = 0_int64
        this%total_cvode_packing_time = 0.0_dp
        this%total_rate_preparation_time = 0.0_dp
        this%total_cvode_reinitialization_time = 0.0_dp
        this%total_integration_time = 0.0_dp
        this%total_qss_projection_time = 0.0_dp
        this%maximum_qss_pre_projection_conservation_residual = 0.0_dp
        this%maximum_qss_post_projection_conservation_residual = 0.0_dp
        this%maximum_qss_post_source_conservation_residual = 0.0_dp
        this%total_cvode_statistics_time = 0.0_dp
        this%total_cvode_result_processing_time = 0.0_dp
        this%total_source_assembly_time = 0.0_dp
        this%total_solver_wall_time = 0.0_dp
        this%last_cvode_packing_time = 0.0_dp
        this%last_rate_preparation_time = 0.0_dp
        this%last_cvode_reinitialization_time = 0.0_dp
        this%last_integration_time = 0.0_dp
        this%last_qss_projection_time = 0.0_dp
        this%last_cvode_statistics_time = 0.0_dp
        this%last_cvode_result_processing_time = 0.0_dp
        this%last_source_assembly_time = 0.0_dp
        this%last_solver_wall_time = 0.0_dp
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
        if (trim(this%ode_solver) == 'qss1') then
            write(output,'(A,I0)') '  QSS1 stoichiometric rank: ', &
                this%qss_stoichiometric_rank
            write(output,'(A,I0)') '  QSS1 invariant count: ', &
                this%qss_conservation_count
            write(output,'(A,I0)') '  QSS1 projection constraints: ', &
                this%qss_projection_constraint_count
            write(output,'(A,L1)') '  QSS1 independent mass row: ', &
                this%qss_has_independent_mass_constraint
        end if
#ifdef CHEMISTRY_PROFILE
        if (trim(this%ode_solver) == 'cvode') then
            write(output,'(A,ES14.6)') &
                '  summed CVODE concentration packing time [s]: ', &
                this%total_cvode_packing_time
        end if
        write(output,'(A,ES14.6)') &
            '  summed rate preparation time [s]: ', &
            this%total_rate_preparation_time
        if (trim(this%ode_solver) == 'cvode') then
            write(output,'(A,ES14.6)') &
                '  summed CVODE reinitialization time [s]: ', &
                this%total_cvode_reinitialization_time
        end if
        write(output,'(A,ES14.6)') '  summed integration time [s]: ', &
            this%total_integration_time
        if (trim(this%ode_solver) == 'qss1') then
            write(output,'(A,ES14.6)') &
                '  summed QSS1 final projection time [s]: ', &
                this%total_qss_projection_time
            write(output,'(A,I0)') &
                '  QSS1 projection positivity clips: ', &
                this%total_qss_projection_clips
            write(output,'(A,I0)') &
                '  QSS1 projection clipped components: ', &
                this%total_qss_projection_clipped_components
            write(output,'(A,ES14.6)') &
                '  QSS1 maximum clipped negative concentration [mol/m3]: ', &
                this%maximum_qss_projection_clip_magnitude
            write(output,'(A,I0)') &
                '  QSS1 maximum clip species index: ', &
                this%maximum_qss_projection_clip_species
            if (this%maximum_qss_projection_clip_species > 0) then
                write(output,'(A,A)') &
                    '  QSS1 maximum clip species name: ', &
                    trim(this%chemistry%chem_ptr%species_names( &
                        this%maximum_qss_projection_clip_species))
            end if
            write(output,'(A,ES14.6)') &
                '  QSS1 max relative invariant residual before projection: ', &
                this%maximum_qss_pre_projection_conservation_residual
            write(output,'(A,ES14.6)') &
                '  QSS1 max relative invariant residual after projection: ', &
                this%maximum_qss_post_projection_conservation_residual
            write(output,'(A,ES14.6)') &
                '  QSS1 max relative invariant residual after source assembly: ', &
                this%maximum_qss_post_source_conservation_residual
        end if
        if (trim(this%ode_solver) == 'cvode') then
            write(output,'(A,ES14.6)') &
                '  summed CVODE statistics-query time [s]: ', &
                this%total_cvode_statistics_time
            write(output,'(A,ES14.6)') &
                '  summed CVODE result-processing time [s]: ', &
                this%total_cvode_result_processing_time
        end if
        write(output,'(A,ES14.6)') &
            '  summed source assembly time [s]: ', &
            this%total_source_assembly_time
        write(output,'(A,ES14.6)') &
            '  chemistry solver wall time [s]: ', &
            this%total_solver_wall_time
#else
        write(output,'(A)') &
            '  detailed counters/timing disabled; compile with CHEMISTRY_PROFILE'
#endif
    end subroutine write_performance_statistics


    subroutine preprocess_mechanism(this)
        class(chemical_kinetics_solver), intent(inout) :: this

        this%kinetics_core = chemical_kinetics_core_c( &
            this%chemistry%chem_ptr,this%thermophysics%thermo_ptr)
    end subroutine preprocess_mechanism



    subroutine solve_chemical_kinetics(this, time_step)
        class(chemical_kinetics_solver), intent(inout) :: this
        real(dp), intent(in) :: time_step

        integer, dimension(3,2) :: cell_loop
        integer :: i, j, k, specie, cvode_worker_index
        real(dp) :: density_cell, temperature_cell, energy_source_cell
        real(dp) :: mass_fraction_sum, negative_tolerance
        integer(int64) :: active_cells_step, integrator_calls_step
        integer(int64) :: internal_steps_step, rhs_evaluations_step
        integer(int64) :: jacobian_evaluations_step
        integer(int64) :: qss_projection_clips_step
        integer(int64) :: qss_projection_clipped_components_step
        integer(int64) :: cell_internal_steps, cell_rhs_evaluations
        integer(int64) :: cell_jacobian_evaluations
        integer(int64) :: cell_qss_projection_clips
        integer(int64) :: cell_qss_projection_clipped_components
        integer :: qss_max_clip_species_step
        integer :: cell_qss_projection_clip_species
        real(dp) :: cvode_packing_time_step, rate_preparation_time_step
        real(dp) :: cvode_reinitialization_time_step, integration_time_step
        real(dp) :: cvode_statistics_time_step
        real(dp) :: cvode_result_processing_time_step
        real(dp) :: qss_projection_time_step
        real(dp) :: qss_max_clip_magnitude_step
        real(dp) :: qss_pre_projection_residual_step
        real(dp) :: qss_post_projection_residual_step
        real(dp) :: qss_post_source_residual_step
        real(dp) :: source_assembly_time_step
        real(dp) :: cell_cvode_packing_time, cell_rate_preparation_time
        real(dp) :: cell_cvode_reinitialization_time, cell_integration_time
        real(dp) :: cell_cvode_statistics_time
        real(dp) :: cell_cvode_result_processing_time
        real(dp) :: cell_qss_projection_time
        real(dp) :: cell_qss_projection_clip_magnitude
        real(dp) :: cell_qss_pre_projection_residual
        real(dp) :: cell_qss_post_projection_residual
        real(dp) :: cell_qss_post_source_residual
        real(dp) :: source_time_start
#ifdef CHEMISTRY_PROFILE
        real(dp) :: solver_wall_time_start, solver_wall_time_step
#endif

        if (.not. ieee_is_finite(time_step) .or. time_step <= 0.0_dp) then
            error stop 'Chemical kinetics: time step must be finite and positive'
        end if

#ifdef CHEMISTRY_PROFILE
        solver_wall_time_start = chemistry_wall_time()
#endif
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
        qss_projection_clips_step = 0_int64
        qss_projection_clipped_components_step = 0_int64
        qss_max_clip_magnitude_step = 0.0_dp
        qss_max_clip_species_step = 0
        cvode_packing_time_step = 0.0_dp
        rate_preparation_time_step = 0.0_dp
        cvode_reinitialization_time_step = 0.0_dp
        integration_time_step = 0.0_dp
        cvode_statistics_time_step = 0.0_dp
        cvode_result_processing_time_step = 0.0_dp
        qss_projection_time_step = 0.0_dp
        qss_pre_projection_residual_step = 0.0_dp
        qss_post_projection_residual_step = 0.0_dp
        qss_post_source_residual_step = 0.0_dp
        source_assembly_time_step = 0.0_dp
#ifdef CHEMISTRY_PROFILE
        solver_wall_time_step = 0.0_dp
#endif

#ifdef NRG_ENABLE_CVODE
        if (trim(this%ode_solver) == 'cvode') then
            call this%ensure_cvode_workers()
        end if
#endif

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
!$omp private(cell_jacobian_evaluations,cell_cvode_packing_time) &
!$omp private(cell_rate_preparation_time,cell_cvode_reinitialization_time) &
!$omp private(cell_integration_time,cell_cvode_statistics_time) &
!$omp private(cell_cvode_result_processing_time,cell_qss_projection_time) &
!$omp private(cell_qss_projection_clips,cell_qss_projection_clipped_components) &
!$omp private(cell_qss_projection_clip_magnitude,cell_qss_projection_clip_species) &
!$omp private(cell_qss_pre_projection_residual,cell_qss_post_projection_residual) &
!$omp private(cell_qss_post_source_residual) &
!$omp private(source_time_start) &
!$omp private(cvode_worker_index) &
!$omp reduction(+:active_cells_step,integrator_calls_step) &
!$omp reduction(+:internal_steps_step,rhs_evaluations_step) &
!$omp reduction(+:jacobian_evaluations_step) &
!$omp reduction(+:qss_projection_clips_step,qss_projection_clipped_components_step) &
!$omp reduction(max:qss_pre_projection_residual_step) &
!$omp reduction(max:qss_post_projection_residual_step,qss_post_source_residual_step) &
!$omp reduction(+:cvode_packing_time_step,rate_preparation_time_step) &
!$omp reduction(+:cvode_reinitialization_time_step,integration_time_step) &
!$omp reduction(+:cvode_statistics_time_step,cvode_result_processing_time_step) &
!$omp reduction(+:qss_projection_time_step,source_assembly_time_step)
        call this%ensure_thread_workspace()
        cvode_worker_index = 1
#ifdef OMP
        if (trim(this%ode_solver) == 'cvode') then
            cvode_worker_index = omp_get_thread_num() + 1
        end if
#endif

!$omp do collapse(3) schedule(dynamic,2)
        do k = cell_loop(3,1), cell_loop(3,2)
            do j = cell_loop(2,1), cell_loop(2,2)
                do i = cell_loop(1,1), cell_loop(1,2)
                    if (markers(i,j,k) /= 0) cycle

                    density_cell = density_field%cells(i,j,k)
                    temperature_cell = temperature_field%cells(i,j,k)
                    if (.not. ieee_is_finite(density_cell) .or. &
                        density_cell <= 0.0_dp) then
                        call report_invalid_cell( &
                            this,'non-positive density',i,j,k, &
                            offending_value=density_cell, &
                            threshold=0.0_dp,time_step=time_step)
                    end if
                    if (.not. ieee_is_finite(temperature_cell) .or. &
                        temperature_cell <= 0.0_dp) then
                        call report_invalid_cell( &
                            this,'non-positive temperature',i,j,k, &
                            offending_value=temperature_cell, &
                            threshold=0.0_dp,time_step=time_step)
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
                                this,'non-finite mass fraction',i,j,k, &
                                specie_index=specie, &
                                offending_value= &
                                    thread_workspace%mass_fraction_cell(specie), &
                                time_step=time_step)
                        end if
                        if (thread_workspace%mass_fraction_cell(specie) < &
                            -negative_tolerance) then
                            call report_invalid_cell( &
                                this,'negative mass fraction',i,j,k, &
                                specie_index=specie, &
                                offending_value= &
                                    thread_workspace%mass_fraction_cell(specie), &
                                threshold=-negative_tolerance, &
                                time_step=time_step)
                        end if
                        thread_workspace%mass_fraction_cell(specie) = max( &
                            thread_workspace%mass_fraction_cell(specie),0.0_dp)
                        mass_fraction_sum = mass_fraction_sum + &
                            thread_workspace%mass_fraction_cell(specie)
                    end do
                    if (mass_fraction_sum <= tiny(1.0_dp)) then
                        call report_invalid_cell( &
                            this,'empty composition',i,j,k, &
                            offending_value=mass_fraction_sum, &
                            threshold=tiny(1.0_dp), &
                            time_step=time_step, &
                            composition_sum=mass_fraction_sum)
                    end if
                    thread_workspace%mass_fraction_cell = &
                        thread_workspace%mass_fraction_cell/mass_fraction_sum

                    cell_internal_steps = 0_int64
                    cell_rhs_evaluations = 0_int64
                    cell_jacobian_evaluations = 0_int64
                    cell_cvode_packing_time = 0.0_dp
                    cell_rate_preparation_time = 0.0_dp
                    cell_cvode_reinitialization_time = 0.0_dp
                    cell_integration_time = 0.0_dp
                    cell_cvode_statistics_time = 0.0_dp
                    cell_cvode_result_processing_time = 0.0_dp
                    cell_qss_projection_time = 0.0_dp
                    cell_qss_projection_clips = 0_int64
                    cell_qss_projection_clipped_components = 0_int64
                    cell_qss_projection_clip_magnitude = 0.0_dp
                    cell_qss_projection_clip_species = 0
                    cell_qss_pre_projection_residual = 0.0_dp
                    cell_qss_post_projection_residual = 0.0_dp
                    cell_qss_post_source_residual = 0.0_dp
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
#ifdef NRG_ENABLE_CVODE
                    case ('cvode')
                        call this%solve_cell_cvode_kinetics( &
                            cvode_worker_index,density_cell,temperature_cell, &
                            thread_workspace%mass_fraction_cell,time_step, &
                            i,j,k,thread_workspace%concentration_increment_cell, &
                            cell_internal_steps,cell_rhs_evaluations, &
                            cell_jacobian_evaluations, &
                            cell_cvode_packing_time,cell_rate_preparation_time, &
                            cell_cvode_reinitialization_time, &
                            cell_integration_time,cell_cvode_statistics_time, &
                            cell_cvode_result_processing_time)
                        integrator_calls_step = integrator_calls_step + 1_int64
                        internal_steps_step = internal_steps_step + cell_internal_steps
                        rhs_evaluations_step = rhs_evaluations_step + cell_rhs_evaluations
                        jacobian_evaluations_step = jacobian_evaluations_step + &
                            cell_jacobian_evaluations
                        cvode_packing_time_step = cvode_packing_time_step + &
                            cell_cvode_packing_time
                        rate_preparation_time_step = rate_preparation_time_step + &
                            cell_rate_preparation_time
                        cvode_reinitialization_time_step = &
                            cvode_reinitialization_time_step + &
                            cell_cvode_reinitialization_time
                        integration_time_step = integration_time_step + &
                            cell_integration_time
                        cvode_statistics_time_step = cvode_statistics_time_step + &
                            cell_cvode_statistics_time
                        cvode_result_processing_time_step = &
                            cvode_result_processing_time_step + &
                            cell_cvode_result_processing_time
#endif
                    case ('qss1')
                        call this%solve_cell_qss1( &
                            density_cell,temperature_cell, &
                            thread_workspace%mass_fraction_cell,time_step, &
                            i,j,k,thread_workspace%concentration_increment_cell, &
                            cell_internal_steps,cell_rhs_evaluations, &
                            cell_jacobian_evaluations, &
                            cell_rate_preparation_time,cell_integration_time, &
                            cell_qss_projection_time, &
                            cell_qss_projection_clips, &
                            cell_qss_projection_clipped_components, &
                            cell_qss_projection_clip_magnitude, &
                            cell_qss_projection_clip_species, &
                            cell_qss_pre_projection_residual, &
                            cell_qss_post_projection_residual)
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
                        qss_projection_time_step = qss_projection_time_step + &
                            cell_qss_projection_time
                        qss_projection_clips_step = &
                            qss_projection_clips_step + &
                            cell_qss_projection_clips
                        qss_projection_clipped_components_step = &
                            qss_projection_clipped_components_step + &
                            cell_qss_projection_clipped_components
                        if (cell_qss_projection_clip_magnitude > 0.0_dp) then
!$omp critical(qss1_clip_maximum)
                            if (cell_qss_projection_clip_magnitude > &
                                    qss_max_clip_magnitude_step) then
                                qss_max_clip_magnitude_step = &
                                    cell_qss_projection_clip_magnitude
                                qss_max_clip_species_step = &
                                    cell_qss_projection_clip_species
                            end if
!$omp end critical(qss1_clip_maximum)
                        end if
                        qss_pre_projection_residual_step = max( &
                            qss_pre_projection_residual_step, &
                            cell_qss_pre_projection_residual)
                        qss_post_projection_residual_step = max( &
                            qss_post_projection_residual_step, &
                            cell_qss_post_projection_residual)
                    case ('table_approximated')
                        call this%interpolate_table_increment( &
                            temperature_cell, &
                            thread_workspace%concentration_increment_cell)
                    case default
                        call report_invalid_cell( &
                            this,'unsupported ODE solver',i,j,k, &
                            time_step=time_step)
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
                    if (trim(this%ode_solver) == 'qss1') then
                        cell_qss_post_source_residual = &
                            this%qss_increment_conservation_residual( &
                                density_cell, &
                                thread_workspace%mass_fraction_cell, &
                                thread_workspace%concentration_increment_cell)
                        qss_post_source_residual_step = max( &
                            qss_post_source_residual_step, &
                            cell_qss_post_source_residual)
                    end if
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
        this%last_cvode_packing_time = cvode_packing_time_step
        this%last_rate_preparation_time = rate_preparation_time_step
        this%last_cvode_reinitialization_time = &
            cvode_reinitialization_time_step
        this%last_integration_time = integration_time_step
        this%last_qss_projection_time = qss_projection_time_step
        this%last_qss_projection_clips = qss_projection_clips_step
        this%last_cvode_statistics_time = cvode_statistics_time_step
        this%last_cvode_result_processing_time = &
            cvode_result_processing_time_step
        this%last_source_assembly_time = source_assembly_time_step
#ifdef CHEMISTRY_PROFILE
        solver_wall_time_step = chemistry_wall_time()-solver_wall_time_start
        this%last_solver_wall_time = solver_wall_time_step
#else
        this%last_solver_wall_time = 0.0_dp
#endif
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
        this%total_cvode_packing_time = this%total_cvode_packing_time + &
            cvode_packing_time_step
        this%total_rate_preparation_time = &
            this%total_rate_preparation_time + rate_preparation_time_step
        this%total_cvode_reinitialization_time = &
            this%total_cvode_reinitialization_time + &
            cvode_reinitialization_time_step
        this%total_integration_time = this%total_integration_time + &
            integration_time_step
        this%total_qss_projection_time = this%total_qss_projection_time + &
            qss_projection_time_step
        this%total_qss_projection_clips = &
            this%total_qss_projection_clips + &
            qss_projection_clips_step
        this%total_qss_projection_clipped_components = &
            this%total_qss_projection_clipped_components + &
            qss_projection_clipped_components_step
        if (qss_max_clip_magnitude_step > &
                this%maximum_qss_projection_clip_magnitude) then
            this%maximum_qss_projection_clip_magnitude = &
                qss_max_clip_magnitude_step
            this%maximum_qss_projection_clip_species = &
                qss_max_clip_species_step
        end if
        this%maximum_qss_pre_projection_conservation_residual = max( &
            this%maximum_qss_pre_projection_conservation_residual, &
            qss_pre_projection_residual_step)
        this%maximum_qss_post_projection_conservation_residual = max( &
            this%maximum_qss_post_projection_conservation_residual, &
            qss_post_projection_residual_step)
        this%maximum_qss_post_source_conservation_residual = max( &
            this%maximum_qss_post_source_conservation_residual, &
            qss_post_source_residual_step)
        this%total_cvode_statistics_time = &
            this%total_cvode_statistics_time + cvode_statistics_time_step
        this%total_cvode_result_processing_time = &
            this%total_cvode_result_processing_time + &
            cvode_result_processing_time_step
        this%total_source_assembly_time = &
            this%total_source_assembly_time + source_assembly_time_step
#ifdef CHEMISTRY_PROFILE
        this%total_solver_wall_time = this%total_solver_wall_time + &
            solver_wall_time_step
        if (chemistry_profile_output_interval > 0) then
            if (mod(this%total_solve_calls, &
                    int(chemistry_profile_output_interval,int64)) == 0_int64) then
                call this%write_performance_statistics()
            end if
        end if
#endif
    end subroutine solve_chemical_kinetics


    subroutine ensure_thread_workspace(this)
        class(chemical_kinetics_solver), intent(in) :: this

        logical :: rebuild
        integer :: work_size

        rebuild = thread_workspace%species_number /= this%species_number .or. &
            thread_workspace%reactions_number /= this%reactions_number
        if (.not. associated(thread_workspace%chemistry, &
            this%chemistry%chem_ptr)) rebuild = .true.
        if (.not. associated(thread_workspace%kinetics_core%thermophysics, &
            this%thermophysics%thermo_ptr)) rebuild = .true.

        if (.not. rebuild) return

        call clear_thread_workspace()
        thread_workspace%species_number = this%species_number
        thread_workspace%reactions_number = this%reactions_number
        thread_workspace%chemistry => this%chemistry%chem_ptr
        thread_workspace%kinetics_core = this%kinetics_core
        allocate(thread_workspace%concentration_initial(this%species_number))
        allocate(thread_workspace%concentration_final(this%species_number))
        allocate(thread_workspace%qss_production(this%species_number))
        allocate(thread_workspace%qss_destruction(this%species_number))
        allocate(thread_workspace%qss_trial(this%species_number))
        allocate(thread_workspace%qss_projected(this%species_number))
        allocate(thread_workspace%qss_projection_weights(this%species_number))
        if (this%qss_projection_constraint_count > 0) then
            allocate(thread_workspace%qss_projection_target( &
                this%qss_projection_constraint_count))
            allocate(thread_workspace%qss_projection_residual( &
                this%qss_projection_constraint_count))
            allocate(thread_workspace%qss_projection_lambda( &
                this%qss_projection_constraint_count))
            allocate(thread_workspace%qss_projection_factor( &
                this%qss_projection_constraint_count, &
                this%qss_projection_constraint_count))
        end if

        allocate(thread_workspace%mass_fraction_cell(this%species_number))
        allocate(thread_workspace%species_source_cell(this%species_number))
        allocate(thread_workspace%concentration_increment_cell( &
            this%species_number))

        work_size = this%species_number**2 + 10*this%species_number + 250
        allocate(thread_workspace%work(work_size))
        allocate(thread_workspace%iwork(50 + this%species_number))

        thread_workspace%concentration_initial = 0.0_dp
        thread_workspace%concentration_final = 0.0_dp
        thread_workspace%qss_production = 0.0_dp
        thread_workspace%qss_destruction = 0.0_dp
        thread_workspace%qss_trial = 0.0_dp
        thread_workspace%qss_projected = 0.0_dp
        thread_workspace%qss_projection_weights = 0.0_dp
        if (allocated(thread_workspace%qss_projection_target)) &
            thread_workspace%qss_projection_target = 0.0_dp
        if (allocated(thread_workspace%qss_projection_residual)) &
            thread_workspace%qss_projection_residual = 0.0_dp
        if (allocated(thread_workspace%qss_projection_lambda)) &
            thread_workspace%qss_projection_lambda = 0.0_dp
        if (allocated(thread_workspace%qss_projection_factor)) &
            thread_workspace%qss_projection_factor = 0.0_dp
        thread_workspace%mass_fraction_cell = 0.0_dp
        thread_workspace%species_source_cell = 0.0_dp
        thread_workspace%concentration_increment_cell = 0.0_dp
        thread_workspace%work = 0.0_dp
        thread_workspace%iwork = 0
    end subroutine ensure_thread_workspace


    subroutine clear_thread_workspace()
        call thread_workspace%rate_state%clear()
        call thread_workspace%kinetics_core%clear()
        if (allocated(thread_workspace%concentration_initial)) &
            deallocate(thread_workspace%concentration_initial)
        if (allocated(thread_workspace%concentration_final)) &
            deallocate(thread_workspace%concentration_final)
        if (allocated(thread_workspace%qss_production)) &
            deallocate(thread_workspace%qss_production)
        if (allocated(thread_workspace%qss_destruction)) &
            deallocate(thread_workspace%qss_destruction)
        if (allocated(thread_workspace%qss_trial)) &
            deallocate(thread_workspace%qss_trial)
        if (allocated(thread_workspace%qss_projected)) &
            deallocate(thread_workspace%qss_projected)
        if (allocated(thread_workspace%qss_projection_weights)) &
            deallocate(thread_workspace%qss_projection_weights)
        if (allocated(thread_workspace%qss_projection_target)) &
            deallocate(thread_workspace%qss_projection_target)
        if (allocated(thread_workspace%qss_projection_residual)) &
            deallocate(thread_workspace%qss_projection_residual)
        if (allocated(thread_workspace%qss_projection_lambda)) &
            deallocate(thread_workspace%qss_projection_lambda)
        if (allocated(thread_workspace%qss_projection_factor)) &
            deallocate(thread_workspace%qss_projection_factor)
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
        hmax = time_step !min(time_step,this%slatec_max_internal_step)
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
            call report_slatec_failure( &
                this,i_cell,j_cell,k_cell,nstate,ierflg,time_step, &
                thread_workspace%concentration_initial, &
                thread_workspace%concentration_final)
        end if

        concentration_scale = max(1.0_dp, &
            maxval(thread_workspace%concentration_initial))
        negative_limit = this%negative_concentration_tolerance* &
            concentration_scale
        do specie = 1, this%species_number
            if (.not. ieee_is_finite( &
                thread_workspace%concentration_final(specie))) then
                call report_invalid_cell( &
                    this,'non-finite final concentration', &
                    i_cell,j_cell,k_cell,specie_index=specie, &
                    offending_value= &
                        thread_workspace%concentration_final(specie), &
                    time_step=time_step, &
                    initial_value= &
                        thread_workspace%concentration_initial(specie), &
                    final_value= &
                        thread_workspace%concentration_final(specie), &
                    concentration_initial= &
                        thread_workspace%concentration_initial, &
                    concentration_final= &
                        thread_workspace%concentration_final)
            end if
            if (thread_workspace%concentration_final(specie) < &
                -negative_limit) then
                call report_invalid_cell( &
                    this,'negative final concentration', &
                    i_cell,j_cell,k_cell,specie_index=specie, &
                    offending_value= &
                        thread_workspace%concentration_final(specie), &
                    threshold=-negative_limit,time_step=time_step, &
                    initial_value= &
                        thread_workspace%concentration_initial(specie), &
                    final_value= &
                        thread_workspace%concentration_final(specie), &
                    concentration_initial= &
                        thread_workspace%concentration_initial, &
                    concentration_final= &
                        thread_workspace%concentration_final)
            end if
            thread_workspace%concentration_final(specie) = max( &
                thread_workspace%concentration_final(specie),0.0_dp)
            concentration_increment(specie) = &
                thread_workspace%concentration_final(specie) - &
                thread_workspace%concentration_initial(specie)
        end do
    end subroutine solve_cell_detailed_kinetics


    subroutine build_qss_conservation_basis(this)
        class(chemical_kinetics_solver), intent(inout) :: this

        real(dp), allocatable :: reaction_basis(:,:), candidate(:)
        real(dp), allocatable :: molar_mass(:),mass_component(:)
        integer :: reaction,component,basis_index,pass,specie
        integer :: rank,invariant,conservation_count,constraint_count
        real(dp) :: coefficient,vector_norm,reference_norm,mass_norm
        real(dp) :: mass_component_norm,maximum_null_residual

        allocate(reaction_basis(this%species_number,this%species_number))
        allocate(candidate(this%species_number))
        allocate(molar_mass(this%species_number))
        allocate(mass_component(this%species_number))
        reaction_basis = 0.0_dp
        candidate = 0.0_dp
        molar_mass = 1.0_dp/this%inverse_molar_mass

        ! Build an orthonormal basis of col(S) by two-pass modified
        ! Gram-Schmidt over the sparse net-stoichiometry vectors.
        rank = 0
        do reaction = 1,this%reactions_number
            candidate = 0.0_dp
            do component = 1,this%kinetics_core%net_count(reaction)
                specie = this%kinetics_core%net_species(component,reaction)
                candidate(specie) = candidate(specie) + real( &
                    this%kinetics_core%net_stoichiometry(component,reaction),dp)
            end do
            reference_norm = sqrt(max(dot_product(candidate,candidate),0.0_dp))
            if (reference_norm <= tiny(1.0_dp)) cycle

            do pass = 1,2
                do basis_index = 1,rank
                    coefficient = dot_product( &
                        reaction_basis(:,basis_index),candidate)
                    candidate = candidate - &
                        coefficient*reaction_basis(:,basis_index)
                end do
            end do
            vector_norm = sqrt(max(dot_product(candidate,candidate),0.0_dp))
            if (vector_norm > qss1_conservation_rank_tolerance* &
                    max(1.0_dp,reference_norm)) then
                rank = rank + 1
                reaction_basis(:,rank) = candidate/vector_norm
            end if
        end do

        this%qss_stoichiometric_rank = rank
        conservation_count = this%species_number-rank
        if (conservation_count <= 0) then
            error stop 'QSS1: mechanism has no stoichiometric invariants'
        end if
        this%qss_conservation_count = conservation_count
        if (allocated(this%qss_conservation_basis)) &
            deallocate(this%qss_conservation_basis)
        allocate(this%qss_conservation_basis( &
            conservation_count,this%species_number))
        this%qss_conservation_basis = 0.0_dp

        ! Complete the reaction-space basis to an orthonormal basis of R^Ns.
        ! The added vectors span null(S^T), hence C*S = 0.  This automatically
        ! includes elemental inventories and independently inert species.
        invariant = 0
        do specie = 1,this%species_number
            candidate = 0.0_dp
            candidate(specie) = 1.0_dp
            do pass = 1,2
                do basis_index = 1,rank
                    coefficient = dot_product( &
                        reaction_basis(:,basis_index),candidate)
                    candidate = candidate - &
                        coefficient*reaction_basis(:,basis_index)
                end do
                do basis_index = 1,invariant
                    coefficient = dot_product( &
                        this%qss_conservation_basis(basis_index,:),candidate)
                    candidate = candidate - coefficient* &
                        this%qss_conservation_basis(basis_index,:)
                end do
            end do
            vector_norm = sqrt(max(dot_product(candidate,candidate),0.0_dp))
            if (vector_norm > qss1_conservation_rank_tolerance) then
                invariant = invariant + 1
                this%qss_conservation_basis(invariant,:) = &
                    candidate/vector_norm
                if (invariant == conservation_count) exit
            end if
        end do
        if (invariant /= conservation_count) then
            error stop 'QSS1: failed to construct conservation basis'
        end if

        ! Verify C*S=0 using every mechanism reaction vector.
        maximum_null_residual = 0.0_dp
        do reaction = 1,this%reactions_number
            candidate = 0.0_dp
            do component = 1,this%kinetics_core%net_count(reaction)
                specie = this%kinetics_core%net_species(component,reaction)
                candidate(specie) = candidate(specie) + real( &
                    this%kinetics_core%net_stoichiometry(component,reaction),dp)
            end do
            do invariant = 1,conservation_count
                maximum_null_residual = max(maximum_null_residual,abs( &
                    dot_product(this%qss_conservation_basis(invariant,:), &
                    candidate)))
            end do
        end do
        if (.not. ieee_is_finite(maximum_null_residual) .or. &
            maximum_null_residual > &
                100.0_dp*qss1_conservation_rank_tolerance) then
            error stop 'QSS1: inaccurate stoichiometric conservation basis'
        end if

        ! NRG's molar-mass table can be intentionally rounded.  Usually the
        ! mass vector already belongs to null(S^T), but if the rounded masses
        ! contain a small component outside that nullspace we append exactly
        ! one independent mass constraint.  This lets the QSS projection obey
        ! both exact stoichiometric invariants and the CFD cell-mass convention.
        mass_component = molar_mass
        do pass = 1,2
            do invariant = 1,conservation_count
                coefficient = dot_product( &
                    this%qss_conservation_basis(invariant,:),mass_component)
                mass_component = mass_component - coefficient* &
                    this%qss_conservation_basis(invariant,:)
            end do
        end do
        mass_norm = sqrt(max(dot_product(molar_mass,molar_mass),0.0_dp))
        mass_component_norm = sqrt(max( &
            dot_product(mass_component,mass_component),0.0_dp))
        this%qss_has_independent_mass_constraint = &
            mass_component_norm > qss1_conservation_rank_tolerance* &
                max(mass_norm,tiny(1.0_dp))

        constraint_count = conservation_count
        if (this%qss_has_independent_mass_constraint) &
            constraint_count = constraint_count + 1
        this%qss_projection_constraint_count = constraint_count
        if (allocated(this%qss_projection_basis)) &
            deallocate(this%qss_projection_basis)
        allocate(this%qss_projection_basis( &
            constraint_count,this%species_number))
        this%qss_projection_basis(1:conservation_count,:) = &
            this%qss_conservation_basis
        if (this%qss_has_independent_mass_constraint) then
            this%qss_projection_basis(constraint_count,:) = &
                mass_component/mass_component_norm
        end if

        deallocate(reaction_basis,candidate,molar_mass,mass_component)
    end subroutine build_qss_conservation_basis


    subroutine project_qss_final_state(this,active_threshold,i_cell,j_cell, &
            k_cell,time_step,clipped_negative,clipped_components, &
            maximum_clip_magnitude,maximum_clip_species, &
            pre_projection_relative_residual,post_projection_relative_residual)
        class(chemical_kinetics_solver), intent(in) :: this
        real(dp), intent(in) :: active_threshold,time_step
        integer, intent(in) :: i_cell,j_cell,k_cell
        logical, intent(out) :: clipped_negative
        integer(int64), intent(out) :: clipped_components
        real(dp), intent(out) :: maximum_clip_magnitude
        integer, intent(out) :: maximum_clip_species
        real(dp), intent(out) :: pre_projection_relative_residual
        real(dp), intent(out) :: post_projection_relative_residual

        integer :: specie,row,column,k
        real(dp) :: correction,residual_value,post_residual
        real(dp) :: gram_value,diagonal_value,sum_value
        real(dp) :: total_concentration,negative_tolerance,conservation_scale

        clipped_negative = .false.
        clipped_components = 0_int64
        maximum_clip_magnitude = 0.0_dp
        maximum_clip_species = 0
        pre_projection_relative_residual = 0.0_dp
        post_projection_relative_residual = 0.0_dp
        total_concentration = sum(thread_workspace%concentration_initial)
        conservation_scale = max(total_concentration,1.0_dp)
        negative_tolerance = qss1_projection_negative_factor* &
            epsilon(1.0_dp)*conservation_scale

        ! Weighted final projection is the normal path.  It is performed once
        ! per complete chemistry cell-call, not once per internal QSS1 step.
        ! Abundant species carry most of the conservation correction while
        ! trace radicals/intermediates are protected from large relative shifts.
        do specie = 1,this%species_number
            thread_workspace%qss_projection_weights(specie) = max( &
                thread_workspace%concentration_final(specie), &
                thread_workspace%concentration_initial(specie), &
                active_threshold,qss1_concentration_floor)
        end do

        ! Enforce B c_projected = B c_initial.
        thread_workspace%qss_projection_residual = -matmul( &
            this%qss_projection_basis, &
            thread_workspace%concentration_final - &
                thread_workspace%concentration_initial)
        pre_projection_relative_residual = maxval(abs( &
            thread_workspace%qss_projection_residual))/conservation_scale

        ! Form and factor G = B W B^T.  The number of constraints is normally
        ! only 3--6, so this remains a tiny cell-local solve and is done once
        ! per chemistry call.
        thread_workspace%qss_projection_factor = 0.0_dp
        do row = 1,this%qss_projection_constraint_count
            do column = 1,row
                gram_value = 0.0_dp
                do specie = 1,this%species_number
                    gram_value = gram_value + &
                        this%qss_projection_basis(row,specie)* &
                        thread_workspace%qss_projection_weights(specie)* &
                        this%qss_projection_basis(column,specie)
                end do
                do k = 1,column-1
                    gram_value = gram_value - &
                        thread_workspace%qss_projection_factor(row,k)* &
                        thread_workspace%qss_projection_factor(column,k)
                end do
                if (row == column) then
                    diagonal_value = gram_value
                    if (.not. ieee_is_finite(diagonal_value) .or. &
                        diagonal_value <= tiny(1.0_dp)) then
                        call report_invalid_cell( &
                            this,'singular weighted QSS1 final projection', &
                            i_cell,j_cell,k_cell, &
                            offending_value=diagonal_value, &
                            threshold=tiny(1.0_dp),time_step=time_step, &
                            concentration_initial= &
                                thread_workspace%concentration_initial, &
                            concentration_final= &
                                thread_workspace%concentration_final)
                    end if
                    thread_workspace%qss_projection_factor(row,column) = &
                        sqrt(diagonal_value)
                else
                    thread_workspace%qss_projection_factor(row,column) = &
                        gram_value/ &
                        thread_workspace%qss_projection_factor(column,column)
                end if
            end do
        end do

        ! Forward/back substitutions for G lambda = residual.
        thread_workspace%qss_projection_lambda = &
            thread_workspace%qss_projection_residual
        do row = 1,this%qss_projection_constraint_count
            sum_value = thread_workspace%qss_projection_lambda(row)
            do k = 1,row-1
                sum_value = sum_value - &
                    thread_workspace%qss_projection_factor(row,k)* &
                    thread_workspace%qss_projection_lambda(k)
            end do
            thread_workspace%qss_projection_lambda(row) = sum_value/ &
                thread_workspace%qss_projection_factor(row,row)
        end do
        do row = this%qss_projection_constraint_count,1,-1
            sum_value = thread_workspace%qss_projection_lambda(row)
            do k = row+1,this%qss_projection_constraint_count
                sum_value = sum_value - &
                    thread_workspace%qss_projection_factor(k,row)* &
                    thread_workspace%qss_projection_lambda(k)
            end do
            thread_workspace%qss_projection_lambda(row) = sum_value/ &
                thread_workspace%qss_projection_factor(row,row)
        end do

        do specie = 1,this%species_number
            correction = dot_product( &
                this%qss_projection_basis(:,specie), &
                thread_workspace%qss_projection_lambda)
            thread_workspace%qss_projected(specie) = &
                thread_workspace%concentration_final(specie) + &
                thread_workspace%qss_projection_weights(specie)*correction
        end do

        ! The weighted projection should normally remain positive because the
        ! raw QSS1 state is positive and trace species have small weights.  A
        ! materially negative projected state is treated as a real numerical
        ! failure rather than silently switching to the Euclidean projector.
        do specie = 1,this%species_number
            if (.not. ieee_is_finite(thread_workspace%qss_projected(specie))) then
                call report_invalid_cell( &
                    this,'non-finite weighted QSS1 final projection', &
                    i_cell,j_cell,k_cell,specie_index=specie, &
                    offending_value=thread_workspace%qss_projected(specie), &
                    time_step=time_step, &
                    concentration_initial= &
                        thread_workspace%concentration_initial, &
                    concentration_final=thread_workspace%qss_projected)
            end if
            if (thread_workspace%qss_projected(specie) < &
                -negative_tolerance) then
                call report_invalid_cell( &
                    this,'weighted QSS1 final projection violates positivity', &
                    i_cell,j_cell,k_cell,specie_index=specie, &
                    offending_value=thread_workspace%qss_projected(specie), &
                    threshold=-negative_tolerance,time_step=time_step, &
                    concentration_initial= &
                        thread_workspace%concentration_initial, &
                    concentration_final=thread_workspace%qss_projected)
            end if
            if (thread_workspace%qss_projected(specie) < 0.0_dp) then
                clipped_components = clipped_components + 1_int64
                if (-thread_workspace%qss_projected(specie) > &
                        maximum_clip_magnitude) then
                    maximum_clip_magnitude = &
                        -thread_workspace%qss_projected(specie)
                    maximum_clip_species = specie
                end if
                thread_workspace%qss_projected(specie) = 0.0_dp
                clipped_negative = .true.
            end if
        end do

        ! Verify every stoichiometric invariant (and the independent mass row,
        ! if one was required when the basis was constructed).
        thread_workspace%qss_projection_residual = matmul( &
            this%qss_projection_basis, &
            thread_workspace%qss_projected - &
                thread_workspace%concentration_initial)
        post_residual = maxval(abs(thread_workspace%qss_projection_residual))
        post_projection_relative_residual = post_residual/conservation_scale
        if (.not. ieee_is_finite(post_residual) .or. &
            post_residual > qss1_conservation_check_tolerance* &
                conservation_scale) then
            call report_invalid_cell( &
                this,'QSS1 final concentration increment is not conservative', &
                i_cell,j_cell,k_cell,offending_value=post_residual, &
                threshold=qss1_conservation_check_tolerance* &
                    conservation_scale,time_step=time_step, &
                concentration_initial=thread_workspace%concentration_initial, &
                concentration_final=thread_workspace%qss_projected)
        end if
    end subroutine project_qss_final_state




    subroutine solve_cell_qss1(this, density, temperature, mass_fraction, &
            time_step, i_cell, j_cell, k_cell, concentration_increment, &
            internal_steps, rhs_evaluations, jacobian_evaluations, &
            rate_preparation_time, integration_time,projection_time, &
            projection_clips,projection_clipped_components, &
            projection_clip_magnitude,projection_clip_species, &
            pre_projection_relative_residual,post_projection_relative_residual)
        class(chemical_kinetics_solver), intent(in) :: this
        real(dp), intent(in) :: density, temperature, time_step
        real(dp), dimension(:), intent(in) :: mass_fraction
        integer, intent(in) :: i_cell, j_cell, k_cell
        real(dp), dimension(:), intent(out) :: concentration_increment
        integer(int64), intent(out) :: internal_steps, rhs_evaluations
        integer(int64), intent(out) :: jacobian_evaluations
        integer(int64), intent(out) :: projection_clips
        integer(int64), intent(out) :: projection_clipped_components
        real(dp), intent(out) :: rate_preparation_time, integration_time
        real(dp), intent(out) :: projection_time
        real(dp), intent(out) :: projection_clip_magnitude
        integer, intent(out) :: projection_clip_species
        real(dp), intent(out) :: pre_projection_relative_residual
        real(dp), intent(out) :: post_projection_relative_residual

        integer :: specie
        real(dp) :: elapsed_time,remaining_time,trial_step,accepted_step
        real(dp) :: total_concentration,active_threshold,relative_change
        real(dp) :: old_concentration,loss_coefficient,loss_argument,phi
        real(dp) :: decay_factor,trial_mass_density,mass_scale
        real(dp) :: projected_mass_density,mass_relative_error
        logical :: accept_step,projection_clipped
#ifdef CHEMISTRY_PROFILE
        real(dp) :: timer_start,projection_timer_start
#endif

        call this%ensure_thread_workspace()

        do specie = 1,this%species_number
            thread_workspace%concentration_initial(specie) = &
                density*mass_fraction(specie)*this%inverse_molar_mass(specie)
        end do
        thread_workspace%concentration_final = &
            thread_workspace%concentration_initial
        thread_workspace%temperature = min(temperature,maximum_rate_temperature)

        total_concentration = sum(thread_workspace%concentration_initial)
        active_threshold = this%qss1_active_concentration_fraction* &
            max(total_concentration,qss1_concentration_floor)

        rate_preparation_time = 0.0_dp
        integration_time = 0.0_dp
        projection_time = 0.0_dp
        projection_clips = 0_int64
        projection_clipped_components = 0_int64
        projection_clip_magnitude = 0.0_dp
        projection_clip_species = 0
        pre_projection_relative_residual = 0.0_dp
        post_projection_relative_residual = 0.0_dp
        internal_steps = 0_int64
        rhs_evaluations = 0_int64
        jacobian_evaluations = 0_int64
#ifdef CHEMISTRY_PROFILE
        timer_start = chemistry_wall_time()
#endif
        call this%prepare_cell_rate_coefficients(temperature)
#ifdef CHEMISTRY_PROFILE
        rate_preparation_time = chemistry_wall_time()-timer_start
        timer_start = chemistry_wall_time()
#endif

        elapsed_time = 0.0_dp
        trial_step = time_step

        do while (elapsed_time < time_step)
            if (internal_steps >= int(this%qss1_max_steps,int64)) then
                call report_invalid_cell( &
                    this,'QSS1 maximum internal step count exceeded', &
                    i_cell,j_cell,k_cell,time_step=time_step, &
                    concentration_initial= &
                        thread_workspace%concentration_initial, &
                    concentration_final= &
                        thread_workspace%concentration_final)
            end if

            remaining_time = time_step-elapsed_time
            trial_step = min(trial_step,remaining_time)

            call thread_workspace%kinetics_core%calculate_production_loss( &
                thread_workspace%rate_state, &
                thread_workspace%concentration_final, &
                thread_workspace%qss_production, &
                thread_workspace%qss_destruction)
            rhs_evaluations = rhs_evaluations + 1_int64

            do
                accept_step = .true.
                do specie = 1,this%species_number
                    old_concentration = max( &
                        thread_workspace%concentration_final(specie), &
                        qss1_concentration_floor)
                    loss_coefficient = &
                        thread_workspace%qss_destruction(specie)/ &
                        old_concentration
                    loss_argument = loss_coefficient*trial_step

                    if (loss_argument <= qss1_small_loss_argument) then
                        phi = 1.0_dp - 0.5_dp*loss_argument + &
                            loss_argument*loss_argument/6.0_dp - &
                            loss_argument*loss_argument*loss_argument/24.0_dp
                        thread_workspace%qss_trial(specie) = &
                            thread_workspace%concentration_final(specie) + &
                            trial_step*( &
                            thread_workspace%qss_production(specie) - &
                            loss_coefficient* &
                            thread_workspace%concentration_final(specie))*phi
                    else
                        decay_factor = exp(-loss_argument)
                        thread_workspace%qss_trial(specie) = &
                            thread_workspace%concentration_final(specie)* &
                            decay_factor + &
                            thread_workspace%qss_production(specie)/ &
                            loss_coefficient*(1.0_dp-decay_factor)
                    end if

                    if (.not. ieee_is_finite(thread_workspace%qss_trial(specie))) then
                        call report_invalid_cell( &
                            this,'non-finite QSS1 trial concentration', &
                            i_cell,j_cell,k_cell,specie_index=specie, &
                            offending_value=thread_workspace%qss_trial(specie), &
                            time_step=time_step, &
                            initial_value= &
                                thread_workspace%concentration_final(specie), &
                            final_value=thread_workspace%qss_trial(specie), &
                            concentration_initial= &
                                thread_workspace%concentration_initial, &
                            concentration_final=thread_workspace%qss_trial)
                    end if
                    thread_workspace%qss_trial(specie) = max( &
                        thread_workspace%qss_trial(specie),0.0_dp)

                    if (thread_workspace%qss_trial(specie) > active_threshold) then
                        relative_change = abs( &
                            thread_workspace%qss_trial(specie)/ &
                            old_concentration - 1.0_dp)
                        if (relative_change > &
                            this%qss1_relative_change_limit .and. &
                            trial_step > this%qss1_minimum_internal_step) then
                            accept_step = .false.
                            exit
                        end if
                    end if
                end do

                if (accept_step) exit
                trial_step = 0.5_dp*trial_step
            end do

            ! Preserve the baseline QSS1 internal-stage behavior.  Conservation
            ! of the full CFD chemistry increment is enforced once after all
            ! internal QSS1 substeps have completed.
            trial_mass_density = 0.0_dp
            do specie = 1,this%species_number
                trial_mass_density = trial_mass_density + &
                    thread_workspace%qss_trial(specie)/ &
                    this%inverse_molar_mass(specie)
            end do
            if (.not. ieee_is_finite(trial_mass_density) .or. &
                trial_mass_density <= tiny(1.0_dp)) then
                call report_invalid_cell( &
                    this,'invalid QSS1 trial mass density', &
                    i_cell,j_cell,k_cell, &
                    offending_value=trial_mass_density, &
                    threshold=tiny(1.0_dp),time_step=time_step, &
                    concentration_initial= &
                        thread_workspace%concentration_initial, &
                    concentration_final=thread_workspace%qss_trial)
            end if
            mass_scale = density/trial_mass_density
            thread_workspace%concentration_final = &
                thread_workspace%qss_trial*mass_scale

            accepted_step = trial_step
            elapsed_time = elapsed_time + accepted_step
            internal_steps = internal_steps + 1_int64

            remaining_time = time_step-elapsed_time
            if (remaining_time <= &
                max(epsilon(time_step)*time_step,tiny(1.0_dp))) exit
            trial_step = min(remaining_time, &
                accepted_step*this%qss1_step_growth_factor)
        end do

#ifdef CHEMISTRY_PROFILE
        integration_time = chemistry_wall_time()-timer_start
        projection_timer_start = chemistry_wall_time()
#endif
        call this%project_qss_final_state( &
            active_threshold,i_cell,j_cell,k_cell,time_step, &
            projection_clipped,projection_clipped_components, &
            projection_clip_magnitude,projection_clip_species, &
            pre_projection_relative_residual, &
            post_projection_relative_residual)
#ifdef CHEMISTRY_PROFILE
        projection_time = chemistry_wall_time()-projection_timer_start
#endif
        if (projection_clipped) projection_clips = 1_int64

        ! The returned chemistry state, not the internal QSS stages, must obey
        ! both the stoichiometric invariants and the NRG cell-mass convention.
        projected_mass_density = 0.0_dp
        do specie = 1,this%species_number
            projected_mass_density = projected_mass_density + &
                thread_workspace%qss_projected(specie)/ &
                this%inverse_molar_mass(specie)
        end do
        mass_relative_error = abs(projected_mass_density-density)/ &
            max(density,tiny(1.0_dp))
        if (.not. ieee_is_finite(mass_relative_error) .or. &
            mass_relative_error > qss1_conservation_check_tolerance) then
            call report_invalid_cell( &
                this,'QSS1 final projected state violates cell mass', &
                i_cell,j_cell,k_cell, &
                offending_value=mass_relative_error, &
                threshold=qss1_conservation_check_tolerance, &
                time_step=time_step, &
                concentration_initial=thread_workspace%concentration_initial, &
                concentration_final=thread_workspace%qss_projected)
        end if

        thread_workspace%concentration_final = thread_workspace%qss_projected
        do specie = 1,this%species_number
            if (.not. ieee_is_finite(thread_workspace%concentration_final(specie))) then
                call report_invalid_cell( &
                    this,'non-finite QSS1 final concentration', &
                    i_cell,j_cell,k_cell,specie_index=specie, &
                    offending_value= &
                        thread_workspace%concentration_final(specie), &
                    time_step=time_step, &
                    initial_value=thread_workspace%concentration_initial(specie), &
                    final_value=thread_workspace%concentration_final(specie), &
                    concentration_initial=thread_workspace%concentration_initial, &
                    concentration_final=thread_workspace%concentration_final)
            end if
            concentration_increment(specie) = &
                thread_workspace%concentration_final(specie) - &
                thread_workspace%concentration_initial(specie)
        end do
    end subroutine solve_cell_qss1


#ifdef NRG_ENABLE_CVODE
    subroutine ensure_cvode_workers(this)
        class(chemical_kinetics_solver), intent(in) :: this

        integer :: worker_count, worker
        logical :: rebuild

        worker_count = 1
#ifdef OMP
        worker_count = max(1,omp_get_max_threads())
#endif

        rebuild = .not. allocated(cvode_workers)
        if (.not. rebuild) then
            rebuild = size(cvode_workers) /= worker_count
        end if
        if (.not. rebuild .and. worker_count > 0) then
            rebuild = cvode_workers(1)%species_number /= this%species_number .or. &
                cvode_workers(1)%reactions_number /= this%reactions_number
            if (.not. associated(cvode_workers(1)%kinetics_core%chemistry, &
                    this%chemistry%chem_ptr)) rebuild = .true.
            if (.not. associated(cvode_workers(1)%kinetics_core%thermophysics, &
                    this%thermophysics%thermo_ptr)) rebuild = .true.
        end if

        if (.not. rebuild) return

        call clear_cvode_workers()
        allocate(cvode_workers(worker_count))
        allocate(cvode_user_contexts(worker_count))

        ! Initialize sequentially before entering the OpenMP cell loop.  Besides
        ! avoiding lazy-init races, this makes the ownership of every SUNDIALS
        ! object deterministic and mirrors the SUNDIALS documented pattern for
        ! one integrator per worker.
        do worker = 1, worker_count
            call initialize_cvode_worker(this,worker)
        end do
    end subroutine ensure_cvode_workers


    subroutine initialize_cvode_worker(this, worker_index)
        class(chemical_kinetics_solver), intent(in) :: this
        integer, intent(in) :: worker_index

        integer(c_int) :: status
        integer(c_int64_t) :: neq

        cvode_workers(worker_index)%species_number = this%species_number
        cvode_workers(worker_index)%reactions_number = this%reactions_number
        cvode_workers(worker_index)%kinetics_core = this%kinetics_core
        allocate(cvode_workers(worker_index)%concentration_initial( &
            this%species_number))
        allocate(cvode_workers(worker_index)%concentration_final( &
            this%species_number))
        cvode_workers(worker_index)%concentration_initial = 0.0_dp
        cvode_workers(worker_index)%concentration_final = 0.0_dp

        neq = int(this%species_number,c_int64_t)
        status = FSUNContext_Create(SUN_COMM_NULL, &
            cvode_workers(worker_index)%sunctx)
        if (status /= 0) then
            error stop 'Chemical kinetics: FSUNContext_Create failed'
        end if

        allocate(cvode_workers(worker_index)%cvode_y_data( &
            this%species_number))
        cvode_workers(worker_index)%cvode_y_data = 0.0_c_double
        cvode_workers(worker_index)%cvode_y => FN_VMake_Serial(neq, &
            cvode_workers(worker_index)%cvode_y_data, &
            cvode_workers(worker_index)%sunctx)
        if (.not. associated(cvode_workers(worker_index)%cvode_y)) then
            error stop 'Chemical kinetics: FN_VMake_Serial failed'
        end if

        cvode_workers(worker_index)%cvode_matrix => FSUNDenseMatrix( &
            neq,neq,cvode_workers(worker_index)%sunctx)
        if (.not. associated(cvode_workers(worker_index)%cvode_matrix)) then
            error stop 'Chemical kinetics: FSUNDenseMatrix failed'
        end if

        cvode_workers(worker_index)%cvode_linear_solver => FSUNLinSol_Dense( &
            cvode_workers(worker_index)%cvode_y, &
            cvode_workers(worker_index)%cvode_matrix, &
            cvode_workers(worker_index)%sunctx)
        if (.not. associated( &
                cvode_workers(worker_index)%cvode_linear_solver)) then
            error stop 'Chemical kinetics: FSUNLinSol_Dense failed'
        end if

        cvode_workers(worker_index)%cvode_mem = FCVodeCreate( &
            CV_BDF,cvode_workers(worker_index)%sunctx)
        if (.not. c_associated(cvode_workers(worker_index)%cvode_mem)) then
            error stop 'Chemical kinetics: FCVodeCreate failed'
        end if

        status = nrg_cvode_init_native( &
            cvode_workers(worker_index)%cvode_mem, &
            c_funloc(cvode_rhs_native),0.0_c_double, &
            c_loc(cvode_workers(worker_index)%cvode_y))
        if (status /= 0) then
            error stop 'Chemical kinetics: native CVodeInit failed'
        end if

        cvode_user_contexts(worker_index)%worker_index = &
            int(worker_index,c_int)
        status = nrg_cvode_set_user_data_native( &
            cvode_workers(worker_index)%cvode_mem, &
            c_loc(cvode_user_contexts(worker_index)))
        if (status /= 0) then
            error stop 'Chemical kinetics: native CVodeSetUserData failed'
        end if

        status = FCVodeSStolerances(cvode_workers(worker_index)%cvode_mem, &
            real(this%cvode_relative_tolerance,c_double), &
            real(this%cvode_absolute_tolerance,c_double))
        if (status /= 0) then
            error stop 'Chemical kinetics: FCVodeSStolerances failed'
        end if

        status = FCVodeSetLinearSolver( &
            cvode_workers(worker_index)%cvode_mem, &
            cvode_workers(worker_index)%cvode_linear_solver, &
            cvode_workers(worker_index)%cvode_matrix)
        if (status /= 0) then
            error stop 'Chemical kinetics: FCVodeSetLinearSolver failed'
        end if

        status = FCVodeSetMaxNumSteps( &
            cvode_workers(worker_index)%cvode_mem, &
            int(this%cvode_max_steps,c_long))
        if (status /= 0) then
            error stop 'Chemical kinetics: FCVodeSetMaxNumSteps failed'
        end if

        cvode_workers(worker_index)%initialized = .true.
    end subroutine initialize_cvode_worker


    subroutine clear_cvode_workers()
        integer :: worker
        integer(c_int) :: cvode_status

        if (.not. allocated(cvode_workers)) then
            if (allocated(cvode_user_contexts)) deallocate(cvode_user_contexts)
            return
        end if

        do worker = 1, size(cvode_workers)
            if (c_associated(cvode_workers(worker)%cvode_mem)) then
                call FCVodeFree(cvode_workers(worker)%cvode_mem)
                cvode_workers(worker)%cvode_mem = c_null_ptr
            end if
            if (associated(cvode_workers(worker)%cvode_linear_solver)) then
                cvode_status = FSUNLinSolFree( &
                    cvode_workers(worker)%cvode_linear_solver)
                nullify(cvode_workers(worker)%cvode_linear_solver)
            end if
            if (associated(cvode_workers(worker)%cvode_matrix)) then
                call FSUNMatDestroy(cvode_workers(worker)%cvode_matrix)
                nullify(cvode_workers(worker)%cvode_matrix)
            end if
            if (associated(cvode_workers(worker)%cvode_y)) then
                call FN_VDestroy(cvode_workers(worker)%cvode_y)
                nullify(cvode_workers(worker)%cvode_y)
            end if
            if (associated(cvode_workers(worker)%cvode_y_data)) then
                deallocate(cvode_workers(worker)%cvode_y_data)
                nullify(cvode_workers(worker)%cvode_y_data)
            end if
            if (c_associated(cvode_workers(worker)%sunctx)) then
                cvode_status = FSUNContext_Free(cvode_workers(worker)%sunctx)
                cvode_workers(worker)%sunctx = c_null_ptr
            end if
            call cvode_workers(worker)%rate_state%clear()
            call cvode_workers(worker)%kinetics_core%clear()
            if (allocated(cvode_workers(worker)%concentration_initial)) &
                deallocate(cvode_workers(worker)%concentration_initial)
            if (allocated(cvode_workers(worker)%concentration_final)) &
                deallocate(cvode_workers(worker)%concentration_final)
            cvode_workers(worker)%species_number = 0
            cvode_workers(worker)%reactions_number = 0
            cvode_workers(worker)%initialized = .false.
        end do

        deallocate(cvode_workers)
        if (allocated(cvode_user_contexts)) deallocate(cvode_user_contexts)
    end subroutine clear_cvode_workers


    subroutine solve_cell_cvode_kinetics(this, worker_index, density, &
            temperature, mass_fraction, time_step, i_cell, j_cell, k_cell, &
            concentration_increment, internal_steps, rhs_evaluations, &
            jacobian_evaluations, packing_time, rate_preparation_time, &
            reinitialization_time, integration_time, statistics_time, &
            result_processing_time)
        class(chemical_kinetics_solver), intent(in) :: this
        integer, intent(in) :: worker_index
        real(dp), intent(in) :: density, temperature, time_step
        real(dp), dimension(:), intent(in) :: mass_fraction
        integer, intent(in) :: i_cell, j_cell, k_cell
        real(dp), dimension(:), intent(out) :: concentration_increment
        integer(int64), intent(out) :: internal_steps, rhs_evaluations
        integer(int64), intent(out) :: jacobian_evaluations
        real(dp), intent(out) :: packing_time, rate_preparation_time
        real(dp), intent(out) :: reinitialization_time, integration_time
        real(dp), intent(out) :: statistics_time, result_processing_time

        integer :: specie
        integer(c_int) :: status
#ifdef CHEMISTRY_PROFILE
        integer(c_long) :: nsteps, nfe, nfe_ls, nje
#endif
        real(c_double) :: tret
        real(dp) :: concentration_scale, negative_limit
#ifdef CHEMISTRY_PROFILE
        real(dp) :: timer_start
#endif

        if (.not. allocated(cvode_workers)) then
            error stop 'Chemical kinetics: CVODE workers are not initialized'
        end if
        if (worker_index < 1 .or. worker_index > size(cvode_workers)) then
            error stop 'Chemical kinetics: invalid CVODE worker index'
        end if
        if (.not. cvode_workers(worker_index)%initialized) then
            error stop 'Chemical kinetics: CVODE worker is not initialized'
        end if

        packing_time = 0.0_dp
        rate_preparation_time = 0.0_dp
        reinitialization_time = 0.0_dp
        integration_time = 0.0_dp
        statistics_time = 0.0_dp
        result_processing_time = 0.0_dp
        internal_steps = 0_int64
        rhs_evaluations = 0_int64
        jacobian_evaluations = 0_int64

#ifdef CHEMISTRY_PROFILE
        timer_start = chemistry_wall_time()
#endif
        do specie = 1, this%species_number
            cvode_workers(worker_index)%concentration_initial(specie) = &
                density*mass_fraction(specie)*this%inverse_molar_mass(specie)
        end do
        cvode_workers(worker_index)%concentration_final = &
            cvode_workers(worker_index)%concentration_initial
        cvode_workers(worker_index)%cvode_y_data = &
            real(cvode_workers(worker_index)%concentration_initial,c_double)
#ifdef CHEMISTRY_PROFILE
        packing_time = chemistry_wall_time()-timer_start
        timer_start = chemistry_wall_time()
#endif
        call cvode_workers(worker_index)%kinetics_core%prepare_rate_state( &
            temperature,cvode_workers(worker_index)%rate_state)
#ifdef CHEMISTRY_PROFILE
        rate_preparation_time = chemistry_wall_time()-timer_start
#endif

        ! The concurrent hot path calls the native SUNDIALS C API directly.
        ! Each worker owns an independent CVODE object graph and an independent
        ! user-owned state array wrapped by its N_Vector.  This avoids the
        ! generated F2003/SWIG pointer-return wrappers that were found to be
        ! unsafe under concurrent Windows/ifx execution.
#ifdef CHEMISTRY_PROFILE
        timer_start = chemistry_wall_time()
#endif
        status = nrg_cvode_reinit_native( &
            cvode_workers(worker_index)%cvode_mem,0.0_c_double, &
            c_loc(cvode_workers(worker_index)%cvode_y))
#ifdef CHEMISTRY_PROFILE
        reinitialization_time = chemistry_wall_time()-timer_start
#endif
        if (status /= 0) then
            call report_invalid_cell(this,'CVODE reinitialization failure', &
                i_cell,j_cell,k_cell,offending_value=real(status,dp), &
                time_step=time_step)
        end if

#ifdef CHEMISTRY_PROFILE
        timer_start = chemistry_wall_time()
#endif
        status = nrg_cvode_step_native( &
            cvode_workers(worker_index)%cvode_mem, &
            real(time_step,c_double), &
            c_loc(cvode_workers(worker_index)%cvode_y),tret,CV_NORMAL)
#ifdef CHEMISTRY_PROFILE
        integration_time = chemistry_wall_time()-timer_start
#endif
        if (status < 0) then
            call report_invalid_cell(this,'CVODE integration failure', &
                i_cell,j_cell,k_cell,offending_value=real(status,dp), &
                time_step=time_step, &
                concentration_initial= &
                    cvode_workers(worker_index)%concentration_initial)
        end if

#ifdef CHEMISTRY_PROFILE
        timer_start = chemistry_wall_time()
        status = nrg_cvode_get_num_steps_native( &
            cvode_workers(worker_index)%cvode_mem,nsteps)
        if (status /= 0) nsteps = 0_c_long
        status = nrg_cvode_get_num_rhs_evals_native( &
            cvode_workers(worker_index)%cvode_mem,nfe)
        if (status /= 0) nfe = 0_c_long
        status = nrg_cvode_get_num_lin_rhs_evals_native( &
            cvode_workers(worker_index)%cvode_mem,nfe_ls)
        if (status /= 0) nfe_ls = 0_c_long
        status = nrg_cvode_get_num_jac_evals_native( &
            cvode_workers(worker_index)%cvode_mem,nje)
        if (status /= 0) nje = 0_c_long

        internal_steps = int(max(nsteps,0_c_long),int64)
        rhs_evaluations = int(max(nfe,0_c_long),int64) + &
            int(max(nfe_ls,0_c_long),int64)
        jacobian_evaluations = int(max(nje,0_c_long),int64)
        statistics_time = chemistry_wall_time()-timer_start
        timer_start = chemistry_wall_time()
#endif

        cvode_workers(worker_index)%concentration_final = &
            real(cvode_workers(worker_index)%cvode_y_data,dp)

        concentration_scale = max(1.0_dp, &
            maxval(cvode_workers(worker_index)%concentration_initial))
        negative_limit = this%negative_concentration_tolerance* &
            concentration_scale
        do specie = 1, this%species_number
            if (.not. ieee_is_finite( &
                    cvode_workers(worker_index)%concentration_final(specie))) then
                call report_invalid_cell(this, &
                    'non-finite CVODE final concentration', &
                    i_cell,j_cell,k_cell,specie_index=specie, &
                    offending_value= &
                        cvode_workers(worker_index)%concentration_final(specie), &
                    time_step=time_step, &
                    concentration_initial= &
                        cvode_workers(worker_index)%concentration_initial, &
                    concentration_final= &
                        cvode_workers(worker_index)%concentration_final)
            end if
            if (cvode_workers(worker_index)%concentration_final(specie) < &
                    -negative_limit) then
                call report_invalid_cell(this, &
                    'negative CVODE final concentration', &
                    i_cell,j_cell,k_cell,specie_index=specie, &
                    offending_value= &
                        cvode_workers(worker_index)%concentration_final(specie), &
                    threshold=-negative_limit,time_step=time_step, &
                    concentration_initial= &
                        cvode_workers(worker_index)%concentration_initial, &
                    concentration_final= &
                        cvode_workers(worker_index)%concentration_final)
            end if
            cvode_workers(worker_index)%concentration_final(specie) = max( &
                cvode_workers(worker_index)%concentration_final(specie),0.0_dp)
            concentration_increment(specie) = &
                cvode_workers(worker_index)%concentration_final(specie) - &
                cvode_workers(worker_index)%concentration_initial(specie)
        end do
#ifdef CHEMISTRY_PROFILE
        result_processing_time = chemistry_wall_time()-timer_start
#endif
    end subroutine solve_cell_cvode_kinetics


    integer(c_int) function cvode_rhs_native( &
            t,sunvec_y,sunvec_f,user_data) result(status) &
            bind(C,name='nrg_cvode_rhs_native')
        real(c_double), value :: t
        type(c_ptr), value :: sunvec_y
        type(c_ptr), value :: sunvec_f
        type(c_ptr), value :: user_data

        type(cvode_user_context), pointer :: context
        type(c_ptr) :: y_data, ydot_data
        real(c_double), pointer :: y(:), ydot(:)
        integer :: worker_index

        if (.not. c_associated(user_data)) then
            status = -1_c_int
            return
        end if
        call c_f_pointer(user_data,context)
        if (.not. associated(context)) then
            status = -1_c_int
            return
        end if

        worker_index = int(context%worker_index)
        if (.not. allocated(cvode_workers)) then
            status = -1_c_int
            return
        end if
        if (worker_index < 1 .or. worker_index > size(cvode_workers)) then
            status = -1_c_int
            return
        end if
        if (.not. cvode_workers(worker_index)%initialized) then
            status = -1_c_int
            return
        end if

        y_data = nrg_nvector_get_array_pointer_native(sunvec_y)
        ydot_data = nrg_nvector_get_array_pointer_native(sunvec_f)
        if (.not. c_associated(y_data) .or. &
                .not. c_associated(ydot_data)) then
            status = -1_c_int
            return
        end if

        call c_f_pointer(y_data,y, &
            [cvode_workers(worker_index)%species_number])
        call c_f_pointer(ydot_data,ydot, &
            [cvode_workers(worker_index)%species_number])
        if (.not. associated(y) .or. .not. associated(ydot)) then
            status = -1_c_int
            return
        end if

        call cvode_workers(worker_index)%kinetics_core%calculate_species_rates( &
            cvode_workers(worker_index)%rate_state,y,ydot)
        status = 0_c_int

        if (t < -huge(1.0_c_double)) status = status
    end function cvode_rhs_native
#endif


    subroutine prepare_cell_rate_coefficients(this, input_temperature)
        class(chemical_kinetics_solver), intent(in) :: this
        real(dp), intent(in) :: input_temperature

        if (this%species_number /= thread_workspace%kinetics_core%species_number) then
            error stop 'Chemical kinetics: shared-core species-count mismatch'
        end if
        call thread_workspace%kinetics_core%prepare_rate_state( &
            input_temperature,thread_workspace%rate_state)
    end subroutine prepare_cell_rate_coefficients



    real(dp) function qss_increment_conservation_residual( &
            this,density,mass_fraction,concentration_increment) result(value)
        class(chemical_kinetics_solver), intent(in) :: this
        real(dp), intent(in) :: density
        real(dp), dimension(:), intent(in) :: mass_fraction
        real(dp), dimension(:), intent(in) :: concentration_increment

        integer :: specie
        real(dp) :: total_concentration,conservation_scale

        if (this%qss_projection_constraint_count <= 0) then
            value = 0.0_dp
            return
        end if

        total_concentration = 0.0_dp
        do specie = 1,this%species_number
            total_concentration = total_concentration + &
                density*max(mass_fraction(specie),0.0_dp)* &
                this%inverse_molar_mass(specie)
        end do
        conservation_scale = max(total_concentration,1.0_dp)

        thread_workspace%qss_projection_residual = matmul( &
            this%qss_projection_basis,concentration_increment)
        value = maxval(abs(thread_workspace%qss_projection_residual))/ &
            conservation_scale
    end function qss_increment_conservation_residual


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
        ! DDRIV3 controls the local weighted state error, not the global
        ! conservation invariant.  Accumulation over internal BDF steps and
        ! subtraction of nearly equal initial/final concentrations can make the
        ! invariant defect several tens of EPS times the cell mass.  A factor of
        ! 100 remains stringent (1e-7 of cell mass for the default EPS=1e-9)
        ! while avoiding false failures in strongly reacting cells.
        select case (this%ode_solver)
        case ('cvode')
            integrator_mass_floor = density*max( &
                100.0_dp*this%cvode_relative_tolerance, &
                1000.0_dp*epsilon(1.0_dp))
        case default
            integrator_mass_floor = density*max( &
                100.0_dp*this%slatec_accuracy, &
                1000.0_dp*epsilon(1.0_dp))
        end select
        integrated_mass_tolerance = max( &
            this%mass_balance_tolerance*integrated_mass_activity, &
            integrator_mass_floor)

        if (abs(integrated_mass_defect) > integrated_mass_tolerance) then
            call report_mass_imbalance( &
                this,i_cell,j_cell,k_cell,density,source_sum,source_l1, &
                integrated_mass_defect,integrated_mass_activity, &
                integrated_mass_tolerance,time_step)
        end if

        ! Remove the integrator-level residual.  Normalize the correction weights
        ! defensively even though the caller already normalizes mass fractions.
        composition_sum = sum(mass_fraction)
        if (composition_sum <= tiny(1.0_dp)) then
            call report_invalid_cell( &
                this,'empty correction composition', &
                i_cell,j_cell,k_cell,offending_value=composition_sum, &
                threshold=tiny(1.0_dp),time_step=time_step, &
                composition_sum=composition_sum)
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
                this,'non-finite chemistry energy source', &
                i_cell,j_cell,k_cell,offending_value=energy_source, &
                time_step=time_step)
        end if
    end subroutine assemble_cell_sources


    logical function can_write_chemical_kinetics_table(this)
        class(chemical_kinetics_solver), intent(in) :: this

        can_write_chemical_kinetics_table = this%record_concentration_increment .and. &
            allocated(this%concentration_increment)
    end function can_write_chemical_kinetics_table


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

        if (.not. associated(thread_workspace%chemistry)) then
            error stop 'Chemical kinetics RHS: thread workspace is not initialized'
        end if
        if (n /= thread_workspace%species_number) then
            error stop 'Chemical kinetics RHS: species-count mismatch'
        end if

        call thread_workspace%kinetics_core%calculate_species_rates( &
            thread_workspace%rate_state,concentration(1:n), &
            concentration_rate(1:n))

        if (time < -huge(1.0_dp)) concentration_rate(1) = &
            concentration_rate(1)
    end subroutine kinetics_rhs




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


    subroutine report_invalid_cell(this, message, i, j, k, &
            specie_index, offending_value, threshold, time_step, &
            initial_value, final_value, composition_sum, &
            concentration_initial, concentration_final)
        class(chemical_kinetics_solver), intent(in) :: this
        character(len=*), intent(in) :: message
        integer, intent(in) :: i, j, k
        integer, intent(in), optional :: specie_index
        real(dp), intent(in), optional :: offending_value, threshold
        real(dp), intent(in), optional :: time_step
        real(dp), intent(in), optional :: initial_value, final_value
        real(dp), intent(in), optional :: composition_sum
        real(dp), dimension(:), intent(in), optional :: concentration_initial
        real(dp), dimension(:), intent(in), optional :: concentration_final

        integer :: specie
        character(len=20) :: specie_name

!$omp critical(chemical_kinetics_error_output)
        write(error_unit,'(A)') ''
        write(error_unit,'(A)') &
            '============================================================'
        write(error_unit,'(A)') 'CHEMICAL KINETICS ERROR'
        write(error_unit,'(A)') &
            '============================================================'
        write(error_unit,'(A,A)') 'Reason                    : ',trim(message)

        if (present(time_step)) then
            call write_cell_diagnostic_context(this,i,j,k,time_step)
        else
            call write_cell_diagnostic_context(this,i,j,k)
        end if

        if (present(specie_index)) then
            if (specie_index >= 1 .and. specie_index <= this%species_number) then
                specie_name = this%chemistry%chem_ptr%species_names(specie_index)
                write(error_unit,'(A,I0)') &
                    'Species index             : ',specie_index
                write(error_unit,'(A,A)') &
                    'Species                   : ',trim(specie_name)
                write(error_unit,'(A,ES24.16)') &
                    'Species mass fraction     : ', &
                    this%mass_fraction%v_ptr%pr(specie_index)%cells(i,j,k)
                write(error_unit,'(A,ES24.16)') &
                    'Species molar mass [kg/mol]: ', &
                    this%thermophysics%thermo_ptr%molar_masses(specie_index)
            else
                write(error_unit,'(A,I0)') &
                    'Invalid species index     : ',specie_index
            end if
        end if

        if (present(offending_value)) then
            write(error_unit,'(A,ES24.16)') &
                'Offending value           : ',offending_value
        end if
        if (present(threshold)) then
            write(error_unit,'(A,ES24.16)') &
                'Threshold / tolerance     : ',threshold
        end if
        if (present(initial_value)) then
            write(error_unit,'(A,ES24.16)') &
                'Initial value             : ',initial_value
        end if
        if (present(final_value)) then
            write(error_unit,'(A,ES24.16)') &
                'Final value               : ',final_value
        end if
        if (present(composition_sum)) then
            write(error_unit,'(A,ES24.16)') &
                'Composition sum           : ',composition_sum
        end if

        if (present(concentration_initial) .and. &
                present(concentration_final)) then
            call write_cell_species_state(this,i,j,k, &
                concentration_initial,concentration_final)
        else
            call write_cell_species_state(this,i,j,k)
        end if

        write(error_unit,'(A)') &
            '============================================================'
!$omp end critical(chemical_kinetics_error_output)

        error stop 'Chemical kinetics solver failed'
    end subroutine report_invalid_cell


    subroutine report_mass_imbalance(this, i, j, k, density, source_sum, &
            source_l1, integrated_defect, integrated_activity, &
            integrated_tolerance, time_step)
        class(chemical_kinetics_solver), intent(in) :: this
        integer, intent(in) :: i, j, k
        real(dp), intent(in) :: density, source_sum, source_l1
        real(dp), intent(in) :: integrated_defect, integrated_activity
        real(dp), intent(in) :: integrated_tolerance, time_step
        real(dp) :: relative_to_density, relative_to_activity

        relative_to_density = abs(integrated_defect)/max(density,tiny(1.0_dp))
        relative_to_activity = abs(integrated_defect)/ &
            max(integrated_activity,tiny(1.0_dp))

!$omp critical(chemical_kinetics_error_output)
        write(error_unit,'(A)') ''
        write(error_unit,'(A)') &
            '============================================================'
        write(error_unit,'(A)') 'CHEMICAL KINETICS MASS IMBALANCE'
        write(error_unit,'(A)') &
            '============================================================'

        call write_cell_diagnostic_context(this,i,j,k,time_step)

        write(error_unit,'(A,ES24.16)') &
            'sum(species_source) [kg m-3 s-1] = ',source_sum
        write(error_unit,'(A,ES24.16)') &
            'sum(abs(species_source))          = ',source_l1
        write(error_unit,'(A,ES24.16)') &
            'integrated mass defect [kg m-3]   = ',integrated_defect
        write(error_unit,'(A,ES24.16)') &
            'defect / cell density             = ',relative_to_density
        write(error_unit,'(A,ES24.16)') &
            'defect / chemistry activity       = ',relative_to_activity
        write(error_unit,'(A,ES24.16)') &
            'allowed integrated defect         = ',integrated_tolerance

        call write_cell_species_state(this,i,j,k)

        write(error_unit,'(A)') &
            '============================================================'
!$omp end critical(chemical_kinetics_error_output)

        error stop 'Chemical kinetics solver failed'
    end subroutine report_mass_imbalance


    subroutine report_slatec_failure(this, i, j, k, nstate, ierflg, &
            time_step, concentration_initial, concentration_final)
        class(chemical_kinetics_solver), intent(in) :: this
        integer, intent(in) :: i, j, k, nstate, ierflg
        real(dp), intent(in) :: time_step
        real(dp), dimension(:), intent(in) :: concentration_initial
        real(dp), dimension(:), intent(in) :: concentration_final

!$omp critical(chemical_kinetics_error_output)
        write(error_unit,'(A)') ''
        write(error_unit,'(A)') &
            '============================================================'
        write(error_unit,'(A)') 'SLATEC CHEMISTRY INTEGRATION FAILURE'
        write(error_unit,'(A)') &
            '============================================================'

        call write_cell_diagnostic_context(this,i,j,k,time_step)

        write(error_unit,'(A,I0)') 'NSTATE                    : ',nstate
        write(error_unit,'(A,I0)') 'IERFLG                    : ',ierflg
        write(error_unit,'(A,I0)') &
            'Internal steps            : ',max(thread_workspace%iwork(3),0)
        write(error_unit,'(A,I0)') &
            'RHS evaluations           : ',max(thread_workspace%iwork(4),0)
        write(error_unit,'(A,I0)') &
            'Jacobian evaluations      : ',max(thread_workspace%iwork(5),0)

        call write_cell_species_state(this,i,j,k, &
            concentration_initial,concentration_final)

        write(error_unit,'(A)') &
            '============================================================'
!$omp end critical(chemical_kinetics_error_output)

        error stop 'Chemical kinetics SLATEC integration failed'
    end subroutine report_slatec_failure


    subroutine write_cell_diagnostic_context(this, i, j, k, time_step)
        class(chemical_kinetics_solver), intent(in) :: this
        integer, intent(in) :: i, j, k
        real(dp), intent(in), optional :: time_step

        write(error_unit,'(A,3(I0,1X))') &
            'Cell (i,j,k)              : ',i,j,k
        write(error_unit,'(A,I0)') &
            'Boundary marker           : ', &
            this%boundary%bc_ptr%bc_markers(i,j,k)
        write(error_unit,'(A,ES24.16)') &
            'Density [kg/m3]           : ',this%density%s_ptr%cells(i,j,k)
        write(error_unit,'(A,ES24.16)') &
            'Temperature [K]           : ', &
            this%temperature%s_ptr%cells(i,j,k)
        write(error_unit,'(A,ES24.16)') &
            'Activation temperature [K]: ',this%activation_temperature
        if (present(time_step)) then
            write(error_unit,'(A,ES24.16)') &
                'CFD time step [s]         : ',time_step
        end if
        write(error_unit,'(A,A)') &
            'ODE solver                : ',trim(this%ode_solver)
        write(error_unit,'(A,A)') &
            'Mechanism                 : ', &
            trim(this%chemistry%chem_ptr%chemical_mechanism_file_name)
        write(error_unit,'(A,I0)') &
            'Species number            : ',this%species_number
        write(error_unit,'(A,I0)') &
            'Reactions number          : ',this%reactions_number
        write(error_unit,'(A,ES24.16)') &
            'Negative concentration tol: ', &
            this%negative_concentration_tolerance
        write(error_unit,'(A,ES24.16)') &
            'SLATEC accuracy           : ',this%slatec_accuracy
        write(error_unit,'(A,ES24.16)') &
            'SLATEC error weight       : ',this%slatec_error_weight
        write(error_unit,'(A,ES24.16)') &
            'SLATEC max internal step  : ', &
            this%slatec_max_internal_step
        write(error_unit,'(A,I0)') &
            'SLATEC max steps          : ',this%slatec_max_steps
        write(error_unit,'(A,ES24.16)') &
            'QSS1 relative change lim : ',this%qss1_relative_change_limit
        write(error_unit,'(A,ES24.16)') &
            'QSS1 minimum step [s]    : ',this%qss1_minimum_internal_step
        write(error_unit,'(A,ES24.16)') &
            'QSS1 active fraction     : ', &
            this%qss1_active_concentration_fraction
        write(error_unit,'(A,ES24.16)') &
            'QSS1 step growth factor  : ',this%qss1_step_growth_factor
        write(error_unit,'(A,I0)') &
            'QSS1 max steps           : ',this%qss1_max_steps
        write(error_unit,'(A,I0)') &
            'QSS1 stoichiometric rank : ',this%qss_stoichiometric_rank
        write(error_unit,'(A,I0)') &
            'QSS1 invariant count     : ',this%qss_conservation_count
        write(error_unit,'(A,I0)') &
            'QSS1 projection constraints: ', &
            this%qss_projection_constraint_count
        write(error_unit,'(A,L1)') &
            'QSS1 independent mass row: ', &
            this%qss_has_independent_mass_constraint
    end subroutine write_cell_diagnostic_context


    subroutine write_cell_species_state(this, i, j, k, &
            concentration_initial, concentration_final)
        class(chemical_kinetics_solver), intent(in) :: this
        integer, intent(in) :: i, j, k
        real(dp), dimension(:), intent(in), optional :: concentration_initial
        real(dp), dimension(:), intent(in), optional :: concentration_final

        integer :: specie
        logical :: have_initial, have_final

        have_initial = present(concentration_initial)
        have_final = present(concentration_final)

        if (have_initial) then
            have_initial = size(concentration_initial) >= this%species_number
        end if
        if (have_final) then
            have_final = size(concentration_final) >= this%species_number
        end if

        write(error_unit,'(A)') 'Species state:'

        if (have_initial .and. have_final) then
            write(error_unit,'(A)') &
                '  idx  species                    Y_k'// &
                '       C_initial [mol/m3]       C_final [mol/m3]'
            do specie = 1, this%species_number
                write(error_unit,'(2X,I4,2X,A20,3(2X,ES24.16))') &
                    specie,trim(this%chemistry%chem_ptr%species_names(specie)), &
                    this%mass_fraction%v_ptr%pr(specie)%cells(i,j,k), &
                    concentration_initial(specie),concentration_final(specie)
            end do
        else
            write(error_unit,'(A)') &
                '  idx  species                    Y_k'
            do specie = 1, this%species_number
                write(error_unit,'(2X,I4,2X,A20,2X,ES24.16)') &
                    specie,trim(this%chemistry%chem_ptr%species_names(specie)), &
                    this%mass_fraction%v_ptr%pr(specie)%cells(i,j,k)
            end do
        end if
    end subroutine write_cell_species_state


    subroutine zero_vector_field(field)
        type(field_vector_cons), pointer, intent(inout) :: field

        integer :: component

        do component = 1, size(field%pr)
            field%pr(component)%cells = 0.0_dp
        end do
    end subroutine zero_vector_field

end module chemical_kinetics_solver_class
