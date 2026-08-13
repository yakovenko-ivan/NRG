module fds_low_mach_solver_class

!===============================================================================
! Variable-density low-Mach solver using the FDS pressure-potential formulation.
!
! The thermodynamic pressure p0 is spatially uniform.  Species density is
! advanced conservatively with CHARM-limited face fluxes, while temperature and
! density are closed by the thermally perfect mixture EOS.  Heat release,
! diffusion, conduction, radiation, particles, and p0 evolution define the
! target velocity divergence.  A predictor/corrector pressure projection then
! advances face velocity through
!
!   du/dt = -grad(H) - F_a - F_b,
!   H = p'/rho + |u|^2/2,
!   F_b = -p' grad(1/rho),
!
! with F_a collecting vorticity transport, reduced-gravity buoyancy, viscous
! momentum diffusion, and particle momentum exchange.  The nonlinear F_b term
! is converged by fixed-point iteration; each linear defect equation is solved
! by elliptic_multigrid_solver_class.
!===============================================================================

    use kind_parameters
    use global_data
    use data_manager_class
    use data_io_class
    use computational_domain_class
    use computational_mesh_class
    use boundary_conditions_class
    use field_pointers
    use table_approximated_real_gas_class
    use thermophysical_properties_class
    use chemical_kinetics_solver_class
    use chemical_properties_class
    use viscosity_solver_class
    use fourier_heat_transfer_solver_class
    use fickean_diffusion_solver_class
    use dispersed_phase_solver_class, only: dispersed_phase_solver, &
        dispersed_phase_solver_c
    use thermal_radiation_solver_class
    use solver_options_class
    use benchmarking
    use supplementary_routines
    use elliptic_multigrid_solver_class

    implicit none

    private
    public :: fds_solver, fds_solver_c

    ! Cell-centred work fields used by the low-Mach predictor/corrector.
    type(field_scalar_cons), target :: p_dyn, p_stat, p_stat_old, dp_stat_dt
    type(field_scalar_cons), target :: rho_int, rho_old
    type(field_scalar_cons), target :: div_v_int, ddiv_v_dt
    type(field_scalar_cons), target :: H, H_old, R
    type(field_vector_cons), target :: Y_int, Y_old

    ! Face-centred velocity and projection acceleration work fields.
    type(field_scalar_flow), target :: F_a, F_b
    type(field_vector_flow), target :: v_f, v_f_old

    type(timer) :: fds_timer
    type(timer) :: fds_gas_dynamics_timer
    type(timer) :: fds_multigrid_timer
    type(timer) :: fds_eos_timer
    type(timer) :: fds_chemistry_timer
    type(timer) :: fds_diffusion_timer
    type(timer) :: fds_heattransfer_timer
    type(timer) :: fds_radiation_timer
    type(timer) :: fds_viscosity_timer
    type :: fds_solver
        ! Enabled physical submodels.
        logical :: diffusion_flag = .false.
        logical :: viscosity_flag = .false.
        logical :: heat_trans_flag = .false.
        logical :: radiation_flag = .false.
        logical :: reactive_flag = .false.
        logical :: CFL_condition_flag = .false.
        logical :: all_Neumann_flag = .true.
        logical :: stabilizing_inlet_flag = .false.

        real(dp) :: courant_fraction = 1.0_dp
        real(dp) :: time = 0.0_dp
        real(dp) :: time_step = 0.0_dp
        real(dp) :: initial_time_step = 0.0_dp
        real(dp) :: inlet_velocity = 0.0_dp
        real(dp) :: buoyancy_reference_density = 0.0_dp
        real(dp), dimension(3) :: g = 0.0_dp
        integer :: additional_particles_phases_number = 0
        integer :: load_counter = 0

        ! Lightweight pressure-solver profiling.  One CSV row is emitted for
        ! each predictor/corrector projection call when enabled.
        logical :: pressure_profiling_enabled = .true.
        integer :: pressure_profile_unit = 0

        ! FDS-only defect-correction smoothing around the shared multigrid solve.
        ! The proper MG-only experiment sets both counts to zero.  Residual
        ! construction/scaling and MG correction accumulation remain unconditional.
        ! settings do not affect the generic elliptic solver used by CABARET.
        integer :: pressure_legacy_pre_smoothing_sweeps = 0
        integer :: pressure_legacy_post_smoothing_sweeps = 0

        ! Coupled physical models and the shared geometric multigrid backend.
        type(viscosity_solver) :: visc_solver
        type(heat_transfer_solver) :: heat_trans_solver
        type(diffusion_solver) :: diff_solver
        type(chemical_kinetics_solver) :: chem_kin_solver
        type(table_approximated_real_gas) :: state_eq
        type(thermal_radiation_solver) :: radiation_solver
        type(elliptic_multigrid_solver) :: pressure_solver
        type(dispersed_phase_solver), dimension(:), allocatable :: particles_solver

        type(computational_domain) :: domain
        type(thermophysical_properties_pointer) :: thermo
        type(chemical_properties_pointer) :: chem
        type(computational_mesh_pointer) :: mesh
        type(boundary_conditions_pointer) :: boundary

        ! Primary and thermodynamic fields owned by the data manager.
        type(field_scalar_cons_pointer) :: rho, rho_int, rho_old
        type(field_scalar_cons_pointer) :: T, p, v_s, mix_mol_mass
        type(field_scalar_cons_pointer) :: E_f, h_s, gamma, mixture_cp
        type(field_scalar_cons_pointer) :: p_stat, p_stat_old, dp_stat_dt, p_dyn
        type(field_scalar_cons_pointer) :: div_v_int, ddiv_v_dt, H, H_old, R
        type(field_scalar_cons_pointer) :: E_f_prod_chem, E_f_prod_heat
        type(field_scalar_cons_pointer) :: E_f_prod_diff, E_f_prod_rad
        type(field_scalar_cons_pointer) :: nu, kappa
        type(field_scalar_flow_pointer) :: F_a, F_b

        type(field_vector_cons_pointer) :: v, v_prod_visc
        type(field_vector_cons_pointer) :: Y, Y_prod_diff, Y_prod_chem
        type(field_vector_cons_pointer) :: Y_int, Y_old, D
        type(field_vector_flow_pointer) :: v_f, v_f_old

        type(field_scalar_cons_pointer), dimension(:), allocatable :: rho_prod_particles
        type(field_scalar_cons_pointer), dimension(:), allocatable :: E_f_prod_particles
        type(field_vector_cons_pointer), dimension(:), allocatable :: Y_prod_particles
        type(field_vector_cons_pointer), dimension(:), allocatable :: v_prod_particles

        ! Curl and acceleration-divergence work arrays for the FDS projection.
        real(dp), dimension(:,:,:,:), allocatable :: vorticity
        real(dp), dimension(:,:,:,:), allocatable :: grad_F_a, grad_F_b
        ! Finest-grid defect and homogeneous correction used by the shared
        ! multigrid adapter.  The obsolete in-class hierarchy has been removed.
        real(dp), dimension(:,:,:), allocatable :: pressure_residual
        real(dp), dimension(:,:,:), allocatable :: pressure_correction
        ! Compact work arrays for repeated prepared elliptic solves.  The
        ! operator hierarchy itself is owned by pressure_solver.
        real(dp), dimension(:,:,:), allocatable :: pressure_adapter_rhs
        real(dp), dimension(:,:,:), allocatable :: pressure_adapter_correction
        logical :: pressure_operator_prepared = .false.
        integer :: pressure_operator_preparations = 0
        integer :: pressure_operator_dimensions = 0
        integer :: pressure_operator_geometry_power = -1
        integer, dimension(3,2) :: pressure_operator_loop = 0
        real(dp), dimension(3) :: pressure_operator_dx = 0.0_dp
        real(dp) :: pressure_operator_base_volume = 0.0_dp

    contains
        procedure, private :: calculate_interm_Y_predictor
        procedure, private :: calculate_divergence_v
        procedure, private :: calculate_pressure_poisson
        procedure, private :: solve_shared_elliptic_correction
        procedure, private :: prepare_shared_elliptic_operator
        procedure, private :: invalidate_pressure_operator_cache
        procedure, private :: write_pressure_profile
        procedure, private :: calculate_dynamic_pressure
        procedure, private :: calculate_velocity
        procedure, private :: calculate_interm_Y_corrector
        procedure, private :: apply_boundary_conditions
        procedure, private :: stabilizing_inlet_1D
        procedure, private :: write_data_table
        procedure :: solve_problem
        procedure :: calculate_time_step
        procedure :: get_time_step
        procedure :: get_time
    end type fds_solver
    interface    fds_solver_c
        module procedure    constructor
    end interface

contains
    !--------------------------------------------------------------------------
    ! Bind data-manager fields, construct enabled physical subsolvers, initialize
    ! restart state and face velocities, and classify the pressure problem as
    ! pure Neumann (closed) or mixed/Dirichlet (contains an outlet).
    !--------------------------------------------------------------------------
    type(fds_solver) function constructor(manager,problem_data_io)

        type(data_manager), intent(inout) :: manager
        type(data_io), intent(inout) :: problem_data_io

        type(field_scalar_cons_pointer) :: scal_ptr
        type(field_vector_cons_pointer) :: vect_ptr
        type(field_tensor_cons_pointer) :: tens_ptr
        type(field_vector_flow_pointer) :: vect_f_ptr
        type(particles_phase) :: particles_params

        integer, dimension(3,2) :: cons_allocation_bounds
        integer, dimension(3,2) :: cons_inner_loop, flow_inner_loop, loop
        integer :: particles_phase_counter
        integer :: bound, dimensions, dim, dim1, i, j, k
        real(dp) :: calculation_time, farfield_velocity
        character(len=20) :: boundary_type_name
        character(len=40) :: var_name


        constructor%diffusion_flag        = manager%solver_options%get_molecular_diffusion_flag()
        constructor%viscosity_flag        = manager%solver_options%get_viscosity_flag()
        constructor%heat_trans_flag        = manager%solver_options%get_heat_transfer_flag()
        constructor%reactive_flag        = manager%solver_options%get_chemical_reaction_flag()
        constructor%radiation_flag      = manager%solver_options%get_thermal_radiation_flag()
        constructor%courant_fraction    = manager%solver_options%get_CFL_condition_coefficient()
        constructor%CFL_condition_flag    = manager%solver_options%get_CFL_condition_flag()
        constructor%stabilizing_inlet_flag = .false.

        constructor%g                   = manager%solver_options%get_grav_acc()

        constructor%additional_particles_phases_number    = manager%solver_options%get_additional_particles_phases_number()

        constructor%domain                = manager%domain
        constructor%thermo%thermo_ptr    => manager%thermophysics%thermo_ptr
        constructor%chem%chem_ptr        => manager%chemistry%chem_ptr
        constructor%boundary%bc_ptr        => manager%boundary_conditions_pointer%bc_ptr
        constructor%mesh%mesh_ptr        => manager%computational_mesh_pointer%mesh_ptr

        call manager%get_cons_field_pointer_by_name(scal_ptr,vect_ptr,tens_ptr,'density')
        constructor%rho%s_ptr                    => scal_ptr%s_ptr
        call manager%get_cons_field_pointer_by_name(scal_ptr,vect_ptr,tens_ptr,'temperature')
        constructor%T%s_ptr                    => scal_ptr%s_ptr
        call manager%get_cons_field_pointer_by_name(scal_ptr,vect_ptr,tens_ptr,'pressure')
        constructor%p%s_ptr                    => scal_ptr%s_ptr
        call manager%get_cons_field_pointer_by_name(scal_ptr,vect_ptr,tens_ptr,'full_energy')
        constructor%E_f%s_ptr                    => scal_ptr%s_ptr
        call manager%get_cons_field_pointer_by_name(scal_ptr,vect_ptr,tens_ptr,'velocity')
        constructor%v%v_ptr                        => vect_ptr%v_ptr
        call manager%get_cons_field_pointer_by_name(scal_ptr,vect_ptr,tens_ptr,'specie_mass_fraction')
        constructor%Y%v_ptr                        => vect_ptr%v_ptr
        call manager%get_cons_field_pointer_by_name(scal_ptr,vect_ptr,tens_ptr,'mixture_molar_mass')
        constructor%mix_mol_mass%s_ptr        => scal_ptr%s_ptr

        call manager%create_scalar_field(p_dyn        ,'pressure_dynamic'                        ,'p_dyn')
        constructor%p_dyn%s_ptr                        => p_dyn
        call manager%create_scalar_field(p_stat        ,'pressure_static'                        ,'p_stat')
        constructor%p_stat%s_ptr                    => p_stat
        call manager%create_scalar_field(p_stat_old    ,'pressure_static_old'                    ,'p_stat_old')
        constructor%p_stat_old%s_ptr                => p_stat_old
        call manager%create_scalar_field(H            ,'stagnation_energy'                    ,'H')
        constructor%H%s_ptr                            => H
        call manager%create_scalar_field(dp_stat_dt    ,'pressure_static_change'                ,'dp_stat_dt')
        constructor%dp_stat_dt%s_ptr                => dp_stat_dt

        call manager%create_scalar_field(rho_int    ,'density_interm'                        ,'rho_int')
        constructor%rho_int%s_ptr                => rho_int
        call manager%create_scalar_field(rho_old    ,'density_old'                            ,'rho_old')
        constructor%rho_old%s_ptr                => rho_old
        call manager%create_scalar_field(div_v_int    ,'velocity_divergence_interm'            ,'div_v_int')
        constructor%div_v_int%s_ptr => div_v_int
        call manager%create_scalar_field(ddiv_v_dt    ,'velocity_divergence_change'            ,'ddiv_v_dt')
        constructor%ddiv_v_dt%s_ptr                => ddiv_v_dt
        call manager%create_scalar_field(H_old        ,'stagnation_energy_old'                ,'H_old')
        constructor%H_old%s_ptr                    => H_old
        call manager%create_scalar_field(R            ,'residual'                                ,'R')
        constructor%R%s_ptr                            => R

        call manager%create_vector_field(Y_int        ,'specie_mass_fraction_interm'    ,'Y_int'    ,'chemical')
        constructor%Y_int%v_ptr                    => Y_int
        call manager%create_vector_field(Y_old        ,'specie_mass_fraction_old'        ,'Y_old'    ,'chemical')
        constructor%Y_old%v_ptr                    => Y_old

        call manager%create_scalar_field(F_a    ,'F_a'                ,'F_a')
        constructor%F_a%s_ptr     => F_a
        call manager%create_scalar_field(F_b    ,'F_b'                ,'F_b')
        constructor%F_b%s_ptr     => F_b

        call manager%create_vector_field(v_f    ,'velocity_flow'                        ,'v_f'        ,'spatial')
        constructor%v_f%v_ptr => v_f
        call manager%create_vector_field(v_f_old,'velocity_flow_old'                    ,'v_f_old'    ,'spatial')
        constructor%v_f_old%v_ptr => v_f_old

        constructor%state_eq    =    table_approximated_real_gas_c(manager)
        call manager%get_cons_field_pointer_by_name(scal_ptr,vect_ptr,tens_ptr,'velocity_of_sound')
        constructor%v_s%s_ptr            => scal_ptr%s_ptr
        call manager%get_cons_field_pointer_by_name(scal_ptr,vect_ptr,tens_ptr,'sensible_enthalpy')
        constructor%h_s%s_ptr            => scal_ptr%s_ptr
        call manager%get_cons_field_pointer_by_name(scal_ptr,vect_ptr,tens_ptr,'adiabatic_index')
        constructor%gamma%s_ptr            => scal_ptr%s_ptr
        call manager%get_cons_field_pointer_by_name(scal_ptr, vect_ptr, tens_ptr, 'mixture_cp')
        constructor%mixture_cp%s_ptr => scal_ptr%s_ptr

        if (constructor%viscosity_flag .or. &
            constructor%additional_particles_phases_number > 0) then
            constructor%visc_solver            = viscosity_solver_c(manager)
            call manager%get_cons_field_pointer_by_name(scal_ptr,vect_ptr,tens_ptr,'velocity_production_viscosity')
            constructor%v_prod_visc%v_ptr            => vect_ptr%v_ptr
            call manager%get_cons_field_pointer_by_name(scal_ptr,vect_ptr,tens_ptr,'viscosity')
            constructor%nu%s_ptr                    => scal_ptr%s_ptr
        end if

        if (constructor%heat_trans_flag .or. &
            constructor%additional_particles_phases_number > 0) then
            constructor%heat_trans_solver    = heat_transfer_solver_c(manager)
            call manager%get_cons_field_pointer_by_name(scal_ptr,vect_ptr,tens_ptr,'energy_production_heat_transfer')
            constructor%E_f_prod_heat%s_ptr        => scal_ptr%s_ptr
            call manager%get_cons_field_pointer_by_name(scal_ptr,vect_ptr,tens_ptr,'thermal_conductivity')
            constructor%kappa%s_ptr                => scal_ptr%s_ptr
        end if

        if (constructor%diffusion_flag) then
            constructor%diff_solver            = diffusion_solver_c(manager)
            call manager%get_cons_field_pointer_by_name(scal_ptr,vect_ptr,tens_ptr,'energy_production_diffusion')
            constructor%E_f_prod_diff%s_ptr            => scal_ptr%s_ptr
            call manager%get_cons_field_pointer_by_name(scal_ptr,vect_ptr,tens_ptr,'specie_production_diffusion')
            constructor%Y_prod_diff%v_ptr            => vect_ptr%v_ptr
            call manager%get_cons_field_pointer_by_name(scal_ptr,vect_ptr,tens_ptr,'diffusivity')
            constructor%D%v_ptr                        => vect_ptr%v_ptr
        end if

        if (constructor%radiation_flag) then
            constructor%radiation_solver    = thermal_radiation_solver_c(manager)
            call manager%get_cons_field_pointer_by_name(scal_ptr,vect_ptr,tens_ptr,'energy_production_radiation')
            constructor%E_f_prod_rad%s_ptr        => scal_ptr%s_ptr
        end if
        if(constructor%additional_particles_phases_number /= 0) then
            allocate(constructor%particles_solver(constructor%additional_particles_phases_number))
            allocate(constructor%rho_prod_particles(constructor%additional_particles_phases_number))
            allocate(constructor%E_f_prod_particles(constructor%additional_particles_phases_number))
            allocate(constructor%v_prod_particles(constructor%additional_particles_phases_number))
            allocate(constructor%Y_prod_particles(constructor%additional_particles_phases_number))
            do particles_phase_counter = 1, constructor%additional_particles_phases_number
                particles_params = manager%solver_options%get_particles_params(particles_phase_counter)
                constructor%particles_solver(particles_phase_counter) = dispersed_phase_solver_c( &
                    manager, particles_params, particles_phase_counter)
                write(var_name,'(A,I2.2)') 'energy_production_particles', particles_phase_counter
                call manager%get_cons_field_pointer_by_name(scal_ptr,vect_ptr,tens_ptr,var_name)
                constructor%E_f_prod_particles(particles_phase_counter)%s_ptr    => scal_ptr%s_ptr
                write(var_name,'(A,I2.2)') 'density_production_particles', particles_phase_counter
                call manager%get_cons_field_pointer_by_name(scal_ptr,vect_ptr,tens_ptr,var_name)
                constructor%rho_prod_particles(particles_phase_counter)%s_ptr    => scal_ptr%s_ptr
                write(var_name,'(A,I2.2)') 'velocity_production_particles', particles_phase_counter
                call manager%get_cons_field_pointer_by_name(scal_ptr,vect_ptr,tens_ptr,var_name)
                constructor%v_prod_particles(particles_phase_counter)%v_ptr        => vect_ptr%v_ptr
                write(var_name,'(A,I2.2)') 'concentration_production_particles', particles_phase_counter
                call manager%get_cons_field_pointer_by_name(scal_ptr,vect_ptr,tens_ptr,var_name)
                constructor%Y_prod_particles(particles_phase_counter)%v_ptr        => vect_ptr%v_ptr
            end do
        end if

        if (constructor%reactive_flag) then
            constructor%chem_kin_solver        = chemical_kinetics_solver_c(manager)
            call manager%get_cons_field_pointer_by_name(scal_ptr,vect_ptr,tens_ptr,'energy_production_chemistry')
            constructor%E_f_prod_chem%s_ptr            => scal_ptr%s_ptr
            call manager%get_cons_field_pointer_by_name(scal_ptr,vect_ptr,tens_ptr,'specie_production_chemistry')
            constructor%Y_prod_chem%v_ptr            => vect_ptr%v_ptr
        end if

        cons_allocation_bounds        = manager%domain%get_local_utter_cells_bounds()

        allocate(constructor%vorticity(        3                        , &
                                            cons_allocation_bounds(1,1):cons_allocation_bounds(1,2), &
                                            cons_allocation_bounds(2,1):cons_allocation_bounds(2,2), &
                                            cons_allocation_bounds(3,1):cons_allocation_bounds(3,2)))

        allocate(constructor%grad_F_a(        3                        , &
                                            cons_allocation_bounds(1,1):cons_allocation_bounds(1,2), &
                                            cons_allocation_bounds(2,1):cons_allocation_bounds(2,2), &
                                            cons_allocation_bounds(3,1):cons_allocation_bounds(3,2)))

        allocate(constructor%grad_F_b(        3                        , &
                                            cons_allocation_bounds(1,1):cons_allocation_bounds(1,2), &
                                            cons_allocation_bounds(2,1):cons_allocation_bounds(2,2), &
                                            cons_allocation_bounds(3,1):cons_allocation_bounds(3,2)))

        constructor%vorticity= 0.0_dp
        constructor%grad_F_a = 0.0_dp
        constructor%grad_F_b = 0.0_dp
        allocate(constructor%pressure_residual( &
            cons_allocation_bounds(1,1):cons_allocation_bounds(1,2), &
            cons_allocation_bounds(2,1):cons_allocation_bounds(2,2), &
            cons_allocation_bounds(3,1):cons_allocation_bounds(3,2)))
        allocate(constructor%pressure_correction, mold=constructor%pressure_residual)
        constructor%pressure_residual = 0.0_dp
        constructor%pressure_correction = 0.0_dp

        cons_inner_loop = manager%domain%get_local_inner_cells_bounds()
        flow_inner_loop    = manager%domain%get_local_inner_faces_bounds()

        problem_data_io                = data_io_c(manager,calculation_time)

        if(problem_data_io%get_load_counter() /= 0) then
            call problem_data_io%add_io_scalar_cons_field(constructor%p_dyn)
        end if

        call problem_data_io%input_all_data()

        constructor%load_counter    = problem_data_io%get_load_counter()

        dimensions        = manager%domain%get_domain_dimensions()

        do dim = 1, dimensions
            loop = flow_inner_loop
            do dim1 = 1,dimensions
                loop(dim1,2) = flow_inner_loop(dim1,2) - (1 - I_m(dim1,dim))
            end do
            do k = loop(3,1),loop(3,2)
            do j = loop(2,1),loop(2,2)
            do i = loop(1,1),loop(1,2)
                do dim1 = 1, dimensions
                    constructor%v_f%v_ptr%pr(dim1)%cells(dim,i,j,k) = 0.5_dp * (constructor%v%v_ptr%pr(dim1)%cells(i,j,k) + &
                        constructor%v%v_ptr%pr(dim1)%cells(i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3)) )
                end do
            end do
            end do
            end do
        end do

        constructor%p_stat%s_ptr%cells    = constructor%p%s_ptr%cells

        if(constructor%load_counter == 1) then
            call problem_data_io%add_io_scalar_cons_field(constructor%p_dyn)

            constructor%p_dyn%s_ptr%cells        = 0.0_dp

            if(constructor%additional_particles_phases_number /= 0) then
                do particles_phase_counter = 1, constructor%additional_particles_phases_number
                    call constructor%particles_solver(particles_phase_counter)%set_initial_distributions()
                end do
            end if
        end if

         if(constructor%load_counter == 1) then
            call constructor%state_eq%apply_state_equation_for_initial_conditions()
        else
            call &
                constructor%state_eq%apply_state_equation_low_mach_fds(manager%solver_options%get_initial_time_step(), &
                predictor=.true.)
            call constructor%state_eq%apply_boundary_conditions_for_initial_conditions()
        end if

        constructor%time                = calculation_time
        constructor%initial_time_step    = manager%solver_options%get_initial_time_step()
        constructor%time_step            = constructor%initial_time_step

        farfield_velocity = 0.0_dp
        constructor%all_Neumann_flag = .true.
        do bound = 1, size(constructor%boundary%bc_ptr%boundary_types)
            boundary_type_name = constructor%boundary%bc_ptr%boundary_types(bound)%get_type_name()
            if (boundary_type_name == 'inlet') then
                farfield_velocity = constructor%boundary%bc_ptr%boundary_types(bound)%get_farfield_velocity()
            else if (boundary_type_name == 'outlet') then
                constructor%all_Neumann_flag = .false.
            end if
        end do

        if (constructor%load_counter == 1) then
            constructor%inlet_velocity = farfield_velocity
        else
            constructor%inlet_velocity = constructor%v%v_ptr%pr(1)%cells( &
                cons_inner_loop(1,1), &
                (cons_inner_loop(2,1)+cons_inner_loop(2,2))/2, &
                cons_inner_loop(3,1))
        end if

        ! The reduced-gravity term uses a reference density representative of
        ! the unburned/far-field state.  This preserves the validated branch.
        constructor%buoyancy_reference_density = constructor%rho%s_ptr%cells( &
            cons_inner_loop(1,2), cons_inner_loop(2,2), 1)

        call manager%create_timer(fds_timer                 ,'FDS solver time'              , 'sol_t')
        call manager%create_timer(fds_gas_dynamics_timer    ,'FDS gas dynamics time'        , 'gd_t')
        call manager%create_timer(fds_eos_timer             ,'FDS eos solver time'          , 'eos_t')
        call manager%create_timer(fds_chemistry_timer       ,'FDS chemistry solver time'    , 'chem_t')
        call manager%create_timer(fds_diffusion_timer       ,'FDS diffusion solver time'    , 'diff_t')
        call manager%create_timer(fds_heattransfer_timer    ,'FDS heattransfer solver time' , 'ht_t')
        call manager%create_timer(fds_radiation_timer       ,'FDS radiation solver time'    , 'rad_t')
        call manager%create_timer(fds_viscosity_timer       ,'FDS viscosity solver time'    , 'visc_t')
        call manager%create_timer(fds_multigrid_timer       ,'Multigrid solver time'        , 'mg_t')

    end function
    !--------------------------------------------------------------------------
    ! Advance one complete predictor/corrector step.  Continuum source rates are
    ! evaluated first, followed by conservative species transport, EOS closure,
    ! divergence construction, and two pressure-potential projections.
    !--------------------------------------------------------------------------
    subroutine solve_problem(this,iteration,stop_flag)
        class(fds_solver)    ,intent(inout)    :: this
        integer                ,intent(in)        :: iteration
        logical                ,intent(inout)    :: stop_flag

        integer    :: particles_phase_counter

        logical    :: stabilized_flag

        call fds_timer%tic()

        this%time = this%time + this%time_step

        call fds_chemistry_timer%tic()
        if (this%reactive_flag)                call this%chem_kin_solver%solve_chemical_kinetics(this%time_step)
        call fds_chemistry_timer%toc(new_iter=.true.)

        call fds_diffusion_timer%tic()
        if (this%diffusion_flag)            call this%diff_solver%solve_diffusion(this%time_step)
        call fds_diffusion_timer%toc(new_iter=.true.)

        call fds_viscosity_timer%tic()
        if (this%viscosity_flag)            call this%visc_solver%solve_viscosity(this%time_step)
        call fds_viscosity_timer%toc(new_iter=.true.)

        call fds_heattransfer_timer%tic()
        if (this%heat_trans_flag)            call this%heat_trans_solver%solve_heat_transfer(this%time_step)
        call fds_heattransfer_timer%toc(new_iter=.true.)

        call fds_radiation_timer%tic()
        if (this%radiation_flag)            call this%radiation_solver%solve_radiation(this%time_step)
        call fds_radiation_timer%toc(new_iter=.true.)

        if (this%additional_particles_phases_number > 0) then
            ! Particle drag and heat transfer require current gas transport
            ! coefficients even when the corresponding gas terms are disabled.
            if (.not. this%heat_trans_flag) &
                call this%heat_trans_solver%update_thermal_conductivity()
            if (.not. this%viscosity_flag) &
                call this%visc_solver%update_dynamic_viscosity()
            do particles_phase_counter = 1, this%additional_particles_phases_number
                call this%particles_solver(particles_phase_counter)%advance( &
                    this%time_step, this%time - this%time_step)

            end do
        end if

        call fds_gas_dynamics_timer%tic()
        call this%calculate_interm_Y_predictor(this%time_step)
        call fds_gas_dynamics_timer%toc()

        call fds_eos_timer%tic()
        call this%state_eq%apply_state_equation_low_mach_fds(this%time_step,predictor=.true.)
        call fds_eos_timer%toc(new_iter=.true.)

        call fds_gas_dynamics_timer%tic()
        call this%apply_boundary_conditions(this%time_step,predictor=.true.)
        call this%calculate_divergence_v        (this%time_step,predictor=.true.)
        call this%calculate_pressure_poisson    (this%time_step,predictor=.true.)
        call this%calculate_velocity            (this%time_step,predictor=.true.)
        call fds_gas_dynamics_timer%toc(new_iter=.true.)

        call fds_viscosity_timer%tic()
        if (this%viscosity_flag)        call this%visc_solver%solve_viscosity(this%time_step)
        call fds_viscosity_timer%toc()

        call fds_gas_dynamics_timer%tic()
        call this%calculate_interm_Y_corrector(this%time_step)
        call fds_gas_dynamics_timer%toc()

        call fds_eos_timer%tic()
        call this%state_eq%apply_state_equation_low_mach_fds(this%time_step,predictor=.false.)
        call fds_eos_timer%toc(new_iter=.true.)

        call fds_gas_dynamics_timer%tic()
        call this%apply_boundary_conditions(this%time_step,predictor=.false.)
        call this%calculate_divergence_v        (this%time_step,predictor=.false.)
        call this%calculate_pressure_poisson    (this%time_step,predictor=.false.)
        call this%calculate_velocity            (this%time_step,predictor=.false.)

        if (this%CFL_condition_flag) then
            call this%calculate_time_step()
        end if

        call fds_gas_dynamics_timer%toc(new_iter=.true.)

        if (this%stabilizing_inlet_flag) then
            call this%stabilizing_inlet_1D(this%time, stabilized_flag)
            if (stabilized_flag) stop_flag = .true.
        end if

        call fds_timer%toc(new_iter=.true.)

    end subroutine
    !--------------------------------------------------------------------------
    ! Conservative Euler predictor for rho*Y_k.  CHARM face reconstruction is
    ! corrected so that the interpolated species densities preserve rho/M at
    ! each face, avoiding pressure errors at material contacts.
    !--------------------------------------------------------------------------
    subroutine calculate_interm_Y_predictor(this,time_step)
        class(fds_solver)    ,intent(inout)    :: this
        real(dp)            ,intent(in)        :: time_step

                real(dp) :: spec_summ

        real(dp)    ,dimension(3)    :: cell_size
        real(dp), dimension (3,3)    :: lame_coeffs
        real(dp), allocatable        :: flux_left_vec(:), flux_right_vec(:), rhs_vec(:), rhoY_new(:)

        integer    :: dimensions, species_number, coord_id
        integer    ,dimension(3,2)    :: cons_inner_loop
        character(len=20)            :: coordinate_system

        integer :: i,j,k,dim,spec,particles_phase_counter

        dimensions        = this%domain%get_domain_dimensions()
        species_number    = this%chem%chem_ptr%species_number

        cons_inner_loop    = this%domain%get_local_inner_cells_bounds()

        cell_size        = this%mesh%mesh_ptr%get_cell_edges_length()
        if (time_step <= 0.0_dp) error stop 'FDS species update: non-positive time step'

        coordinate_system    = this%domain%get_coordinate_system_name()
        select case(coordinate_system)
            case ('cartesian')
                coord_id = 0
            case ('cylindrical')
                coord_id = 1
            case ('spherical')
                coord_id = 2
            case default
                coord_id = 0
        end select

        associate (    rho                    => this%rho%s_ptr            , &
                    rho_int                => this%rho_int%s_ptr        , &
                    rho_old             => this%rho_old%s_ptr        , &
                    v_f                    => this%v_f%v_ptr            , &
                    Y                    => this%Y%v_ptr                , &
                    Y_int                => this%Y_int%v_ptr            , &
                    Y_old               => this%Y_old%v_ptr         , &
                    mesh                => this%mesh%mesh_ptr        , &
                    bc                    => this%boundary%bc_ptr        , &
                    thermo                => this%thermo%thermo_ptr)

            ! Optimized scalar predictor:
            !   * species face values are computed as vectors once per face/direction;
            !   * EOS correction is applied once per face, not once per species;
            !   * boundedness and reconstruction are done in the same cell loop;
            !   * coordinate-system string checks are moved out of the inner loop.
            !$omp parallel default(shared) &
            !$omp& private(i,j,k,dim,spec,spec_summ,lame_coeffs,flux_left_vec,flux_right_vec,rhs_vec,rhoY_new, &
            !$omp& particles_phase_counter)
            allocate(flux_left_vec(species_number), flux_right_vec(species_number), rhs_vec(species_number), &
                rhoY_new(species_number))

            !$omp do collapse(3) schedule(static)
            do k = cons_inner_loop(3,1),cons_inner_loop(3,2)
            do j = cons_inner_loop(2,1),cons_inner_loop(2,2)
            do i = cons_inner_loop(1,1),cons_inner_loop(1,2)
                if(bc%bc_markers(i,j,k) == 0) then

                    lame_coeffs        = 1.0_dp
                    if (coord_id == 1) then
                        ! cylindrical: x -> r, y -> z
                        lame_coeffs(1,1)    =  mesh%mesh(1,i,j,k) - 0.5_dp*cell_size(1)
                        lame_coeffs(1,2)    =  mesh%mesh(1,i,j,k)
                        lame_coeffs(1,3)    =  mesh%mesh(1,i,j,k) + 0.5_dp*cell_size(1)
                    else if (coord_id == 2) then
                        ! spherical: x -> r
                        lame_coeffs(1,1)    =  (mesh%mesh(1,i,j,k) - 0.5_dp*cell_size(1))**2
                        lame_coeffs(1,2)    =  (mesh%mesh(1,i,j,k))**2
                        lame_coeffs(1,3)    =  (mesh%mesh(1,i,j,k) + 0.5_dp*cell_size(1))**2
                    end if

                    rhs_vec = 0.0_dp
                    do dim = 1, dimensions
                        if ((i*I_m(dim,1) + j*I_m(dim,2)  + k*I_m(dim,3)) < cons_inner_loop(dim,2)) then
                            call eos_corrected_species_face_vector(rho,Y,thermo%molar_masses,species_number,dim,i,j,k,1, &
                                 v_f%pr(dim)%cells(dim,i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3)),flux_right_vec)
                        else
                            if (v_f%pr(dim)%cells(dim,i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3)) > 0.0_dp) then
                                do spec = 1,species_number
                                    flux_right_vec(spec) = rho%cells(i,j,k) * Y%pr(spec)%cells(i,j,k)
                                end do
                            else
                                do spec = 1,species_number
                                    flux_right_vec(spec) = rho%cells(i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3)) * &
                                        Y%pr(spec)%cells(i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3))
                                end do
                            end if
                        end if

                        if ((i*I_m(dim,1) + j*I_m(dim,2)  + k*I_m(dim,3)) > cons_inner_loop(dim,1)) then
                            call eos_corrected_species_face_vector(rho,Y,thermo%molar_masses,species_number,dim,i,j,k,-1, &
                                 v_f%pr(dim)%cells(dim,i,j,k),flux_left_vec)
                        else
                            if (v_f%pr(dim)%cells(dim,i,j,k) > 0.0_dp) then
                                do spec = 1,species_number
                                    flux_left_vec(spec) = rho%cells(i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3)) * &
                                        Y%pr(spec)%cells(i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3))
                                end do
                            else
                                do spec = 1,species_number
                                    flux_left_vec(spec) = rho%cells(i,j,k) * Y%pr(spec)%cells(i,j,k)
                                end do
                            end if
                        end if

                        do spec = 1,species_number
                            rhs_vec(spec) = rhs_vec(spec) -(     flux_right_vec(spec) * &
                                v_f%pr(dim)%cells(dim,i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3)) * lame_coeffs(dim,3) &
                                                        - flux_left_vec(spec)  * v_f%pr(dim)%cells(dim,i,j,k) * &
                                                            lame_coeffs(dim,1)) /cell_size(dim)/lame_coeffs(dim,2)
                        end do
                    end do

                    do spec = 1, species_number
                        rhoY_new(spec) = rho%cells(i,j,k)*Y%pr(spec)%cells(i,j,k) + rhs_vec(spec)*time_step
                        if (this%diffusion_flag) rhoY_new(spec) = rhoY_new(spec) + &
                            this%Y_prod_diff%v_ptr%pr(spec)%cells(i,j,k) * time_step
                        if (this%reactive_flag)  rhoY_new(spec) = rhoY_new(spec) + &
                            this%Y_prod_chem%v_ptr%pr(spec)%cells(i,j,k) * time_step
                        if (this%additional_particles_phases_number /= 0) then
                            do particles_phase_counter = 1, this%additional_particles_phases_number
                                rhoY_new(spec) = rhoY_new(spec) + &
                                    this%Y_prod_particles(particles_phase_counter)%v_ptr%pr(spec)%cells(i,j,k) * time_step
                            end do
                        end if
                    end do

                    ! Conservative boundedness correction and reconstruction.
                    spec_summ = 0.0_dp
                    do spec = 1,species_number
                        rhoY_new(spec) = max(rhoY_new(spec), 0.0_dp)
                        spec_summ = spec_summ + rhoY_new(spec)
                    end do
                    if (spec_summ <= tiny(1.0_dp)) then
                        error stop 'FDS species update: non-positive total density'
                    else
                        rho_int%cells(i,j,k)    = spec_summ
                        do spec = 1, species_number
                            Y_int%pr(spec)%cells(i,j,k) = rhoY_new(spec) / rho_int%cells(i,j,k)
                        end do
                    end if
                end if
            end do
            end do
            end do
            !$omp end do

            ! The final copy cannot be merged with the flux-computation loop, because
            ! predictor fluxes must see the old rho,Y values in neighbouring cells.
            !$omp do collapse(3) schedule(static)
            do k = cons_inner_loop(3,1),cons_inner_loop(3,2)
            do j = cons_inner_loop(2,1),cons_inner_loop(2,2)
            do i = cons_inner_loop(1,1),cons_inner_loop(1,2)
                if(bc%bc_markers(i,j,k) == 0) then
                    rho_old%cells(i,j,k)    = rho%cells(i,j,k)
                    rho%cells(i,j,k)        = rho_int%cells(i,j,k)
                    do spec = 1, species_number
                        Y_old%pr(spec)%cells(i,j,k)    = Y%pr(spec)%cells(i,j,k)
                        Y%pr(spec)%cells(i,j,k)        = Y_int%pr(spec)%cells(i,j,k)
                    end do
                end if
            end do
            end do
            end do
            !$omp end do

            deallocate(flux_left_vec, flux_right_vec, rhs_vec, rhoY_new)
            !$omp end parallel
            end associate
    end subroutine
    !--------------------------------------------------------------------------
    ! Assemble the low-Mach divergence target from sensible-enthalpy and species
    ! source terms.  In a closed domain, solve the global compatibility relation
    ! for dp0/dt; an outlet instead fixes p0 and therefore dp0/dt=0.
    !--------------------------------------------------------------------------
    subroutine calculate_divergence_v(this,time_step,predictor)
        class(fds_solver)    ,intent(inout)    :: this
        real(dp)            ,intent(in)        :: time_step
        logical                ,intent(in)        :: predictor

        real(dp) :: flux_left, flux_right, specie_enthalpy, mixture_cp
        real(dp) :: D_sum, P_sum, U_sum
        real(dp), dimension(this%chem%chem_ptr%species_number) :: species_flux_left
        real(dp), dimension(this%chem%chem_ptr%species_number) :: species_flux_right
        real(dp), dimension(this%chem%chem_ptr%species_number) :: species_coefficient
        real(dp), dimension(this%chem%chem_ptr%species_number) :: species_reference_enthalpy

        real(dp)                    :: cell_volume, base_cell_volume
        real(dp)    ,dimension(3)    :: cell_size, cell_surface_area, base_cell_surface_area

        real(dp), dimension (3,3)    :: lame_coeffs
        character(len=20)            :: coordinate_system

        integer :: dimensions, species_number
        integer    ,dimension(3,2)    :: cons_inner_loop

        integer    :: bound_number, plus, sign
        integer :: i,j,k,dim,spec,particles_phase_counter

        dimensions        = this%domain%get_domain_dimensions()
        species_number    = this%chem%chem_ptr%species_number

        coordinate_system    = this%domain%get_coordinate_system_name()

        cons_inner_loop    = this%domain%get_local_inner_cells_bounds()

        cell_size             = this%mesh%mesh_ptr%get_cell_edges_length()
        base_cell_volume      = this%mesh%mesh_ptr%get_cell_volume()
        base_cell_surface_area = this%mesh%mesh_ptr%get_cell_surface_area()
        if (time_step <= 0.0_dp) error stop 'FDS divergence: non-positive time step'

        ! Species reference enthalpies depend only on the thermodynamic database
        ! and T_ref.  Evaluate them once per divergence assembly rather than once
        ! per cell and per spatial direction.
        do spec = 1, species_number
            species_reference_enthalpy(spec) = &
                this%thermo%thermo_ptr%specie_enthalpy_molar(T_ref,spec)
        end do

        D_sum = 0.0_dp
        P_sum = 0.0_dp
        U_sum = 0.0_dp


        associate (    div_v_int            => this%div_v_int%s_ptr        , &
                    T                    => this%T%s_ptr                , &
                    rho                    => this%rho%s_ptr            , &
                    p_stat                => this%p_stat%s_ptr        , &
                    h_s                    => this%h_s%s_ptr            , &
                    dp_stat_dt            => this%dp_stat_dt%s_ptr    , &
                    v_f                    => this%v_f%v_ptr            , &
                    Y                    => this%Y%v_ptr                , &
                    mix_mol_mass        => this%mix_mol_mass%s_ptr  , &
                    mixture_cp_field    => this%mixture_cp%s_ptr    , &
                    thermo                => this%thermo%thermo_ptr    , &
                    mesh                => this%mesh%mesh_ptr        , &
                    bc                    => this%boundary%bc_ptr)

            !$omp parallel default(shared) &
            !$omp& private(flux_right,flux_left,i,j,k,dim,spec,mixture_cp,specie_enthalpy,plus,sign,bound_number,cell_volume, &
            !$omp& cell_surface_area,lame_coeffs,species_flux_left,species_flux_right,species_coefficient)

            !$omp do collapse(3) schedule(static)    reduction(+:D_sum,P_sum)
            do k = cons_inner_loop(3,1),cons_inner_loop(3,2)
            do j = cons_inner_loop(2,1),cons_inner_loop(2,2)
            do i = cons_inner_loop(1,1),cons_inner_loop(1,2)
                if(bc%bc_markers(i,j,k) == 0) then

                    cell_volume    = base_cell_volume
                    lame_coeffs    = 1.0_dp

                    select case(coordinate_system)
                        case ('cartesian')

                        case ('cylindrical')
                            ! x -> r, y -> z
                            lame_coeffs(1,1)    =  mesh%mesh(1,i,j,k) - 0.5_dp*cell_size(1)
                            lame_coeffs(1,2)    =  mesh%mesh(1,i,j,k)
                            lame_coeffs(1,3)    =  mesh%mesh(1,i,j,k) + 0.5_dp*cell_size(1)
                            cell_volume            = cell_volume * mesh%mesh(1,i,j,k)
                        case ('spherical')
                            ! x -> r
                            lame_coeffs(1,1)    =  (mesh%mesh(1,i,j,k) - 0.5_dp*cell_size(1))**2
                            lame_coeffs(1,2)    =  (mesh%mesh(1,i,j,k))**2
                            lame_coeffs(1,3)    =  (mesh%mesh(1,i,j,k) + 0.5_dp*cell_size(1))**2
                            cell_volume            = cell_volume * mesh%mesh(1,i,j,k)**2
                    end select

                    div_v_int%cells(i,j,k) = 0.0_dp

                    if (this%heat_trans_flag)    div_v_int%cells(i,j,k) = div_v_int%cells(i,j,k) + &
                        this%E_f_prod_heat%s_ptr%cells(i,j,k)
                    ! [J/m^3/s]
                    if (this%diffusion_flag)    div_v_int%cells(i,j,k) = div_v_int%cells(i,j,k) + &
                        this%E_f_prod_diff%s_ptr%cells(i,j,k)
                    ! [J/m^3/s]
                    if (this%reactive_flag)        div_v_int%cells(i,j,k) = div_v_int%cells(i,j,k) + &
                        this%E_f_prod_chem%s_ptr%cells(i,j,k)
                    ! [J/m^3/s]
                    if (this%radiation_flag)    div_v_int%cells(i,j,k) = div_v_int%cells(i,j,k) + &
                        this%E_f_prod_rad%s_ptr%cells(i,j,k)

                    if (this%additional_particles_phases_number /= 0) then
                        do particles_phase_counter = 1, this%additional_particles_phases_number
                            div_v_int%cells(i,j,k) = div_v_int%cells(i,j,k) + &
                                this%E_f_prod_particles(particles_phase_counter)%s_ptr%cells(i,j,k)
                            ! [J/m^3/s]
                        end do
                    end if

                    do dim = 1,dimensions

                        if ((i*I_m(dim,1) + j*I_m(dim,2)  + k*I_m(dim,3)) < cons_inner_loop(dim,2)) then
                            flux_right = charm_face_value( &
                                rho%cells(i-I_m(dim,1):i+2*I_m(dim,1),j-I_m(dim,2):j+2*I_m(dim,2),k-I_m(dim,3):k+2*I_m(dim,3)) &
                                * &
                                                                    h_s%cells(i-I_m(dim,1):i+2*I_m(dim,1),j-I_m(dim,2):j+2* &
                                                                        I_m(dim,2),k-I_m(dim,3):k+2*I_m(dim,3)), &
                                                                    v_f%pr(dim)%cells(dim,i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim, &
                                                                        3)))
                        else
                            if (v_f%pr(dim)%cells(dim,i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3)) > 0.0_dp) then
                                flux_right = rho%cells(i,j,k) * h_s%cells(i,j,k)
                            else
                                flux_right = rho%cells(i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3)) * &
                                    h_s%cells(i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3))
                            end if
                        end if

                        if ((i*I_m(dim,1) + j*I_m(dim,2)  + k*I_m(dim,3)) > cons_inner_loop(dim,1)) then
                            flux_left = charm_face_value( &
                                rho%cells(i-2*I_m(dim,1):i+I_m(dim,1),j-2*I_m(dim,2):j+I_m(dim,2),k-2*I_m(dim,3):k+I_m(dim,3)) &
                                * &
                                                                    h_s%cells(i-2*I_m(dim,1):i+I_m(dim,1),j-2*I_m(dim,2):j+ &
                                                                        I_m(dim,2),k-2*I_m(dim,3):k+I_m(dim,3)), &
                                                                    v_f%pr(dim)%cells(dim,i,j,k))
                        else
                            if (v_f%pr(dim)%cells(dim,i,j,k) > 0.0_dp) then
                                flux_left = rho%cells(i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3)) * &
                                    h_s%cells(i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3))
                            else
                                flux_left = rho%cells(i,j,k) * h_s%cells(i,j,k)
                            end if
                        end if

                        div_v_int%cells(i,j,k) = div_v_int%cells(i,j,k)  -  ( &
                            v_f%pr(dim)%cells(dim,i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3))    * lame_coeffs(dim,3) * &
                            (flux_right    -  rho%cells(i,j,k) * h_s%cells(i,j,k)) / cell_size(dim) &
                                                                              -    v_f%pr(dim)%cells(dim,i,j,k) &
                                                                                  * lame_coeffs(dim,1) * (flux_left    - &
                                                                                  rho%cells(i,j,k) * h_s%cells(i,j,k)) / &
                                                                                  cell_size(dim)) / lame_coeffs(dim,2)
                    end do

                    ! The EOS is called immediately before this routine and has
                    ! already evaluated the molar mixture heat capacity for the
                    ! current T,Y state.  Reuse that field instead of rebuilding
                    ! mole fractions and evaluating mixture_cp_molar again.
                    mixture_cp = mixture_cp_field%cells(i,j,k)

                    div_v_int%cells(i,j,k) = div_v_int%cells(i,j,k) / &
                        (rho%cells(i,j,k) * mixture_cp * T%cells(i,j,k) / mix_mol_mass%cells(i,j,k))

                    ! Composition changes alter density through both the mixture
                    ! molar mass and the sensible-enthalpy equation.  The
                    ! thermodynamic coefficient is a cell/species quantity, not a
                    ! face-direction quantity, so form it once and reuse it for
                    ! convective, diffusive, reactive and particle source terms.
                    do spec = 1, species_number
                        specie_enthalpy = (thermo%specie_enthalpy_molar(T%cells(i,j,k),spec) - &
                            species_reference_enthalpy(spec)) / thermo%molar_masses(spec)
                        species_coefficient(spec) = mix_mol_mass%cells(i,j,k) / &
                            (thermo%molar_masses(spec)*rho%cells(i,j,k)) - &
                            specie_enthalpy / (rho%cells(i,j,k)*mixture_cp*T%cells(i,j,k)/ &
                            mix_mol_mass%cells(i,j,k))
                    end do

                    ! Reconstruct all conservative species densities once per face
                    ! and apply the already-assembled cell/species coefficients.
                    do dim = 1, dimensions
                        if ((i*I_m(dim,1) + j*I_m(dim,2) + k*I_m(dim,3)) < &
                            cons_inner_loop(dim,2)) then
                            call eos_corrected_species_face_vector(rho, Y, thermo%molar_masses, &
                                species_number, dim, i, j, k, 1, &
                                v_f%pr(dim)%cells(dim,i+I_m(dim,1),j+I_m(dim,2), &
                                k+I_m(dim,3)), species_flux_right)
                        else
                            do spec = 1, species_number
                                if (v_f%pr(dim)%cells(dim,i+I_m(dim,1),j+I_m(dim,2), &
                                    k+I_m(dim,3)) > 0.0_dp) then
                                    species_flux_right(spec) = rho%cells(i,j,k) * &
                                        Y%pr(spec)%cells(i,j,k)
                                else
                                    species_flux_right(spec) = rho%cells(i+I_m(dim,1), &
                                        j+I_m(dim,2),k+I_m(dim,3)) * &
                                        Y%pr(spec)%cells(i+I_m(dim,1),j+I_m(dim,2), &
                                        k+I_m(dim,3))
                                end if
                            end do
                        end if

                        if ((i*I_m(dim,1) + j*I_m(dim,2) + k*I_m(dim,3)) > &
                            cons_inner_loop(dim,1)) then
                            call eos_corrected_species_face_vector(rho, Y, thermo%molar_masses, &
                                species_number, dim, i, j, k, -1, &
                                v_f%pr(dim)%cells(dim,i,j,k), species_flux_left)
                        else
                            do spec = 1, species_number
                                if (v_f%pr(dim)%cells(dim,i,j,k) > 0.0_dp) then
                                    species_flux_left(spec) = rho%cells(i-I_m(dim,1), &
                                        j-I_m(dim,2),k-I_m(dim,3)) * &
                                        Y%pr(spec)%cells(i-I_m(dim,1),j-I_m(dim,2), &
                                        k-I_m(dim,3))
                                else
                                    species_flux_left(spec) = rho%cells(i,j,k) * &
                                        Y%pr(spec)%cells(i,j,k)
                                end if
                            end do
                        end if

                        do spec = 1, species_number
                            div_v_int%cells(i,j,k) = div_v_int%cells(i,j,k) + &
                                species_coefficient(spec) * (-( &
                                v_f%pr(dim)%cells(dim,i+I_m(dim,1),j+I_m(dim,2), &
                                k+I_m(dim,3))*lame_coeffs(dim,3) * &
                                (species_flux_right(spec)-rho%cells(i,j,k)* &
                                Y%pr(spec)%cells(i,j,k))/cell_size(dim) - &
                                v_f%pr(dim)%cells(dim,i,j,k)*lame_coeffs(dim,1) * &
                                (species_flux_left(spec)-rho%cells(i,j,k)* &
                                Y%pr(spec)%cells(i,j,k))/cell_size(dim)) / &
                                lame_coeffs(dim,2))
                        end do
                    end do

                    do spec = 1, species_number
                        if (this%diffusion_flag) then
                            div_v_int%cells(i,j,k) = div_v_int%cells(i,j,k) + &
                                species_coefficient(spec)*this%Y_prod_diff%v_ptr%pr(spec)%cells(i,j,k)
                        end if
                        if (this%reactive_flag) then
                            div_v_int%cells(i,j,k) = div_v_int%cells(i,j,k) + &
                                species_coefficient(spec)*this%Y_prod_chem%v_ptr%pr(spec)%cells(i,j,k)
                        end if
                        do particles_phase_counter = 1, this%additional_particles_phases_number
                            div_v_int%cells(i,j,k) = div_v_int%cells(i,j,k) + &
                                species_coefficient(spec) * &
                                this%Y_prod_particles(particles_phase_counter)%v_ptr%pr(spec)%cells(i,j,k)
                        end do
                    end do

                    D_sum = D_sum + div_v_int%cells(i,j,k) * cell_volume
                    P_sum = P_sum + (1.0_dp / p_stat%cells(i,j,k) - 1.0_dp / (rho%cells(i,j,k) * mixture_cp * T%cells(i,j,k) / &
                        mix_mol_mass%cells(i,j,k))) * cell_volume
                end if
            end do
            end do
            end do
            !$omp end do

            !$omp do collapse(3) schedule(static)    reduction(+:U_sum)
            do k = cons_inner_loop(3,1),cons_inner_loop(3,2)
            do j = cons_inner_loop(2,1),cons_inner_loop(2,2)
            do i = cons_inner_loop(1,1),cons_inner_loop(1,2)
                if(bc%bc_markers(i,j,k) == 0) then
                    do dim = 1,dimensions
                        do plus = 1,2
                            sign            = (-1)**plus
                            bound_number    = bc%bc_markers(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3))
                            if( bound_number /= 0 ) then
                                cell_surface_area    = base_cell_surface_area
                                select case(coordinate_system)
                                    case ('cartesian')
                                        cell_surface_area    = cell_surface_area
                                    case ('cylindrical')
                                        ! x -> r, y -> z
                                        if(dim==1) cell_surface_area(dim) = cell_surface_area(dim) * mesh%mesh(1,i,j,k)
                                        if(dim==2) cell_surface_area(dim) = cell_surface_area(dim) * mesh%mesh(1,i,j,k)
                                        ! - 0.5_dp*cell_size(1)
                                    case ('spherical')
                                        ! x -> r
                                        if(dim==1) cell_surface_area(dim) = cell_surface_area(dim) * (mesh%mesh(1,i,j,k))**2
                                        ! - 0.5_dp*cell_size(1)
                                end select
                                U_sum = U_sum + sign * &
                                    v_f%pr(dim)%cells(dim,i+max(sign,0)*I_m(dim,1),j+max(sign,0)*I_m(dim,2),k+max(sign,0)* &
                                    I_m(dim,3)) * cell_surface_area(dim)
                            end if
                        end do
                    end do
                end if
            end do
            end do
            end do
            !$omp end do

            !$omp do collapse(3) schedule(static)
            do k = cons_inner_loop(3,1),cons_inner_loop(3,2)
            do j = cons_inner_loop(2,1),cons_inner_loop(2,2)
            do i = cons_inner_loop(1,1),cons_inner_loop(1,2)
                if(bc%bc_markers(i,j,k) == 0) then

                    mixture_cp = mixture_cp_field%cells(i,j,k)
                    if (this%all_Neumann_flag) then
                        if (abs(P_sum) <= tiny(1.0_dp)) then
                            error stop 'FDS divergence constraint: singular pressure compatibility denominator'
                        end if
                        dp_stat_dt%cells(i,j,k) = (D_sum-U_sum)/P_sum
                        div_v_int%cells(i,j,k) = div_v_int%cells(i,j,k) - &
                            (1.0_dp/p_stat%cells(i,j,k) - 1.0_dp/(rho%cells(i,j,k) * &
                            mixture_cp*T%cells(i,j,k)/mix_mol_mass%cells(i,j,k))) * &
                            dp_stat_dt%cells(i,j,k)
                    else
                        ! An outlet fixes the background thermodynamic pressure.
                        dp_stat_dt%cells(i,j,k) = 0.0_dp
                    end if

                end if
            end do
            end do
            end do
            !$omp end do

            !$omp end parallel

            end associate

    end subroutine
    !--------------------------------------------------------------------------
    ! Solve the FDS pressure-potential equation.  The explicit acceleration F_a
    ! and baroclinic term F_b are assembled on faces.  Red-black/Jacobi-like
    ! smoothing is followed by a homogeneous multigrid defect correction from
    ! elliptic_multigrid_solver_class.  The outer fixed-point loop updates
    ! F_b=-p' grad(1/rho) until its face acceleration changes sufficiently little.
    !--------------------------------------------------------------------------
    subroutine calculate_pressure_poisson(this, time_step, predictor)
		class(fds_solver)	,intent(inout)	:: this
		real(dp)			,intent(in)		:: time_step
		logical				,intent(in)		:: predictor
		real(dp) :: sub_F_summ, r_summ1, r_summ2
		real(dp) :: a_norm_init, a_norm, sum_ddiv_v_dt
		real(dp) :: grad_F_a_summ, grad_F_b_summ
		real(dp) :: farfield_velocity
		real(dp)	:: rho_face
		real(dp)	:: F_b_candidate
		real(dp)	:: F_b_change_max_local, F_b_change_max
		real(dp)	:: F_b_norm_max_local, F_b_norm_max
		real(dp)	:: F_a_norm_max_local, F_a_norm_max
		real(dp) :: F_b_convergence_scale
		real(dp), parameter :: F_b_relative_tolerance = 1.0e-3_dp
		real(dp), parameter :: F_b_absolute_tolerance = 1.0e-10_dp
		integer, parameter :: maximum_pressure_iterations = 200

		real(dp), dimension (3,3)	:: lame_coeffs
		character(len=20)				:: coordinate_system
		real(dp)	,dimension(3)	:: cell_size
		real(dp) :: beta

		integer :: dimensions
		integer	,dimension(3,2)	:: cons_inner_loop, cons_utter_loop, flow_inner_loop
		integer	,dimension(3,2)	:: loop
		character(len=20)		:: boundary_type_name,boundary_type_name1,boundary_type_name2,boundary_type_name3
		integer	:: particles_phase_counter
		integer	:: sign, bound_number, bound_number1, bound_number2, bound_number3
		integer :: i, j, k, dim, dim1, dim2, plus

		integer :: legacy_smoothing_iteration, pressure_iteration, v_cycle_iteration
		integer			:: nu_0
		real(dp)		:: tolerance

        type(elliptic_solver_statistics) :: multigrid_statistics
        integer :: profile_legacy_smoothing_sweeps, profile_multigrid_calls
        integer :: profile_multigrid_cycles, profile_multigrid_smoothing_iterations
        integer :: profile_multigrid_parallel_smoothing, profile_multigrid_serial_smoothing
        integer :: profile_multigrid_residual_evaluations, profile_max_hierarchy_levels
        real(dp) :: profile_multigrid_initial_residual_max
        real(dp) :: profile_multigrid_final_residual_max
        real(dp) :: profile_multigrid_relative_residual_max

		logical :: v_cycle_converged = .false.
		logical			:: pressure_converged = .false.

		dimensions		= this%domain%get_domain_dimensions()

        if (this%pressure_legacy_pre_smoothing_sweeps < 0 .or. &
            this%pressure_legacy_post_smoothing_sweeps < 0) then
            error stop 'FDS pressure legacy smoothing sweep counts must be non-negative.'
        end if

        profile_legacy_smoothing_sweeps = 0
        profile_multigrid_calls = 0
        profile_multigrid_cycles = 0
        profile_multigrid_smoothing_iterations = 0
        profile_multigrid_parallel_smoothing = 0
        profile_multigrid_serial_smoothing = 0
        profile_multigrid_residual_evaluations = 0
        profile_max_hierarchy_levels = 1
        profile_multigrid_initial_residual_max = 0.0_dp
        profile_multigrid_final_residual_max = 0.0_dp
        profile_multigrid_relative_residual_max = 0.0_dp

		coordinate_system	= this%domain%get_coordinate_system_name()

		cons_inner_loop	= this%domain%get_local_inner_cells_bounds()
		cons_utter_loop = this%domain%get_local_utter_cells_bounds()


		flow_inner_loop	= this%domain%get_local_inner_faces_bounds()

		cell_size		= this%mesh%mesh_ptr%get_cell_edges_length()
            associate (	vorticity		=> this%vorticity			, &
                        ddiv_v_dt		=> this%ddiv_v_dt%s_ptr		, &
					    F_a				=> this%F_a%s_ptr			, &
					    F_b				=> this%F_b%s_ptr)

			ddiv_v_dt%cells	= 0.0_dp
            F_a%cells		= 0.0_dp
			F_b%cells		= 0.0_dp
            vorticity		= 0.0_dp

            end associate

			sum_ddiv_v_dt	= 0.0_dp
			grad_F_a_summ	= 0.0_dp
			r_summ1			= 0.0_dp
            r_summ2			= 0.0_dp

            associate ( rho_old			    => this%rho_old%s_ptr		    , &
                        rho_int			    => this%rho_int%s_ptr		    , &
                        div_v_int		    => this%div_v_int%s_ptr		    , &
                        ddiv_v_dt		    => this%ddiv_v_dt%s_ptr		    , &
                        F_a				    => this%F_a%s_ptr			    , &
                        vorticity		    => this%vorticity			    , &
                        grad_F_a		    => this%grad_F_a			    , &
                        v_f				    => this%v_f%v_ptr			    , &
                        v_f_old			    => this%v_f_old%v_ptr		    , &
                        mesh			    => this%mesh%mesh_ptr		    , &
                        bc				    => this%boundary%bc_ptr		)

            !$omp parallel default(shared) &
            !$omp& private(i,j,k,dim,dim1,dim2,loop,lame_coeffs,farfield_velocity,rho_face,sign,bound_number, &
            !$omp& bound_number1,bound_number2,bound_number3,plus,boundary_type_name,boundary_type_name1, &
            !$omp& boundary_type_name2,boundary_type_name3)

            !$omp do collapse(3) schedule(static)	reduction(+:sum_ddiv_v_dt,r_summ1)
			do k = cons_inner_loop(3,1),cons_inner_loop(3,2)
			do j = cons_inner_loop(2,1),cons_inner_loop(2,2)
			do i = cons_inner_loop(1,1),cons_inner_loop(1,2)
				if(bc%bc_markers(i,j,k) == 0) then

					lame_coeffs		= 1.0_dp

					select case(coordinate_system)
						case ('cartesian')
							lame_coeffs			= 1.0_dp

                            r_summ1				= r_summ1 + 1.0_dp
						case ('cylindrical')
							! x -> r, y -> z
							lame_coeffs(1,1)	=  mesh%mesh(1,i,j,k) - 0.5_dp*cell_size(1)
							lame_coeffs(1,2)	=  mesh%mesh(1,i,j,k)
							lame_coeffs(1,3)	=  mesh%mesh(1,i,j,k) + 0.5_dp*cell_size(1)

                            r_summ1				= r_summ1 + mesh%mesh(1,i,j,k)
						case ('spherical')
							! x -> r
							lame_coeffs(1,1)	=  (mesh%mesh(1,i,j,k) - 0.5_dp*cell_size(1))**2
							lame_coeffs(1,2)	=  (mesh%mesh(1,i,j,k))**2
							lame_coeffs(1,3)	=  (mesh%mesh(1,i,j,k) + 0.5_dp*cell_size(1))**2

                            r_summ1				= r_summ1 + mesh%mesh(1,i,j,k)**2
					end select

					if (predictor) then
						ddiv_v_dt%cells(i,j,k)	= div_v_int%cells(i,j,k) / time_step

						do dim = 1, dimensions
							ddiv_v_dt%cells(i,j,k)	= ddiv_v_dt%cells(i,j,k) - (v_f%pr(dim)%cells(dim,i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3)) * &
							    lame_coeffs(dim,3)  - v_f%pr(dim)%cells(dim,i,j,k) * lame_coeffs(dim,1)) / cell_size(1) / lame_coeffs(dim,2) / &
							    time_step
						end do
					else
						ddiv_v_dt%cells(i,j,k)	= div_v_int%cells(i,j,k) / (0.5_dp * time_step)

						do dim = 1, dimensions
							ddiv_v_dt%cells(i,j,k)	= ddiv_v_dt%cells(i,j,k)  - 0.5_dp * ( &
							    (v_f%pr(dim)%cells(dim,i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3)) * lame_coeffs(dim,3) - v_f%pr(dim)%cells(dim,i,j,k) &
							    * lame_coeffs(dim,1)) / lame_coeffs(dim,2) / cell_size(1) &
																							+ (v_f_old%pr(dim)%cells(dim,i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3)) * lame_coeffs(dim,3) - &
																							    v_f_old%pr(dim)%cells(dim,i,j,k) * lame_coeffs(dim,1)) / lame_coeffs(dim,2) / cell_size(1)) / &
																							    (0.5_dp * time_step)
						end do
					end if

					sum_ddiv_v_dt = sum_ddiv_v_dt + ddiv_v_dt%cells(i,j,k) * lame_coeffs(1,2)
				end if
			end do
			end do
			end do
			!$omp end do

			if (this%all_Neumann_flag) then
			!$omp do collapse(3) schedule(static)
			do k = cons_inner_loop(3,1),cons_inner_loop(3,2)
			do j = cons_inner_loop(2,1),cons_inner_loop(2,2)
			do i = cons_inner_loop(1,1),cons_inner_loop(1,2)
				if(bc%bc_markers(i,j,k) == 0) then
					ddiv_v_dt%cells(i,j,k)	=	ddiv_v_dt%cells(i,j,k) - sum_ddiv_v_dt/r_summ1
				end if
			end do
			end do
			end do
			!$omp end do
			end if

			!$omp do collapse(3) schedule(static)
			do k = flow_inner_loop(3,1),flow_inner_loop(3,2)
			do j = flow_inner_loop(2,1),flow_inner_loop(2,2)
			do i = flow_inner_loop(1,1),flow_inner_loop(1,2)
				do dim =  1, 3
				do dim1 = 1, dimensions
				do dim2 = 1, dimensions
						if ((dim1 /= dim).and.(dim2 /= dim1).and.(dim2 /= dim)) then
							if(((dim1-dim) == 1).or.((dim1-dim) == -2)) then
								vorticity(dim,i,j,k) = vorticity(dim,i,j,k) - v_f%pr(dim1)%cells(dim1,i,j,k) + &
								    v_f%pr(dim1)%cells(dim1,i-I_m(dim2,1),j-I_m(dim2,2),k-I_m(dim2,3))
							else
								vorticity(dim,i,j,k) = vorticity(dim,i,j,k) + v_f%pr(dim1)%cells(dim1,i,j,k) - &
								    v_f%pr(dim1)%cells(dim1,i-I_m(dim2,1),j-I_m(dim2,2),k-I_m(dim2,3))
							end if
						end if
				end do
				end do

					vorticity(dim,i,j,k) = vorticity(dim,i,j,k) / cell_size(1)
					if(abs(vorticity(dim,i,j,k)) < 1e-03) vorticity(dim,i,j,k) = 0.0_dp
				end do


			end do
			end do
			end do
			!$omp end do

			!$omp do collapse(3) schedule(static)
			do k = cons_inner_loop(3,1),cons_inner_loop(3,2)
			do j = cons_inner_loop(2,1),cons_inner_loop(2,2)
			do i = cons_inner_loop(1,1),cons_inner_loop(1,2)
				if(bc%bc_markers(i,j,k) == 0) then
					do dim = 1, dimensions
						do plus = 1,2
							sign			= (-1)**plus
							bound_number1	= bc%bc_markers(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3))
							do dim2 = 1, dimensions
								if(dim2 /= dim) then
									bound_number2	= &
									    bc%bc_markers(i+sign*I_m(dim,1)-I_m(dim2,1),j+sign*I_m(dim,2)-I_m(dim2,2),k+sign*I_m(dim,3)-I_m(dim2,3))
									bound_number3	= &
									    bc%bc_markers(i+sign*I_m(dim,1)+I_m(dim2,1),j+sign*I_m(dim,2)+I_m(dim2,2),k+sign*I_m(dim,3)+I_m(dim2,3))
									exit
								else
									bound_number2	= 0
									bound_number3	= 0
								end if
							end do
							if ((( bound_number1 /= 0 ).and.(bound_number2 /= 0)).or.(( bound_number1 /= 0).and.(bound_number3/= 0))) then

								boundary_type_name1 = bc%boundary_types(bound_number1)%get_type_name()
								if (bound_number2 /= 0) boundary_type_name2 = bc%boundary_types(bound_number2)%get_type_name()
								if (bound_number3 /= 0) boundary_type_name3 = bc%boundary_types(bound_number3)%get_type_name()

								boundary_type_name = 'outlet'

								if (( bound_number1 /= 0 ).and.(bound_number2 /= 0)) then
									if ((boundary_type_name1=='wall').or.(boundary_type_name2=='wall')) boundary_type_name = 'wall'
									if ((boundary_type_name1=='inlet').or.(boundary_type_name2=='inlet')) boundary_type_name = 'inlet'


									select case(boundary_type_name)
										case('wall')
										case('inlet')
											continue
										case('outlet')
											if (sign < 0) then
												vorticity(3,i,j,k) = 0.0_dp
												! 2.0*vorticity(3,i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3)) -
												! vorticity(3,i+2*I_m(dim,1),j+2*I_m(dim,2),k+2*I_m(dim,3))
												! !0.5_dp*(vorticity(3,i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3)) +
												! vorticity(3,i+2*I_m(dim,1),j+2*I_m(dim,2),k+2*I_m(dim,3)))
												! !2.0*vorticity(3,i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3)) -
												! vorticity(3,i+2*I_m(dim,1),j+2*I_m(dim,2),k+2*I_m(dim,3))
												if (dim2 == 1) then
													v_f%pr(dim2)%cells(dim2,i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3)) =  cell_size(1)*(vorticity(3,i,j,k) - &
													    (v_f%pr(dim)%cells(dim,i,j,k) - &
													    v_f%pr(dim)%cells(dim,i-I_m(dim2,1),j-I_m(dim2,2),k-I_m(dim2,3)))/cell_size(1) + &
													    v_f%pr(dim2)%cells(dim2,i,j,k)/cell_size(1))
												else
													v_f%pr(dim2)%cells(dim2,i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3)) = -cell_size(1)*(vorticity(3,i,j,k) + &
													    (v_f%pr(dim)%cells(dim,i,j,k) - &
													    v_f%pr(dim)%cells(dim,i-I_m(dim2,1),j-I_m(dim2,2),k-I_m(dim2,3)))/cell_size(1) - &
													    v_f%pr(dim2)%cells(dim2,i,j,k)/cell_size(1))
												end if
											else
												vorticity(3,i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3)) = 0.0_dp
												! 2.0*vorticity(3,i,j,k) - vorticity(3,i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3)) !0.5_dp*(vorticity(3,i,j,k) +
												! vorticity(3,i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3))) !2.0*vorticity(3,i,j,k) -
												! vorticity(3,i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3))
												if (dim2 == 1) then
													v_f%pr(dim2)%cells(dim2,i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3)) = &
													    -cell_size(1)*(vorticity(3,i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3)) - &
													    (v_f%pr(dim)%cells(dim,i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3)) - &
													    v_f%pr(dim)%cells(dim,i-I_m(dim2,1)+I_m(dim,1),j-I_m(dim2,2)+I_m(dim,2),k-I_m(dim2,3)+I_m(dim,3)))/ &
													    cell_size(1) - v_f%pr(dim2)%cells(dim2,i,j,k)/cell_size(1))
												else
													v_f%pr(dim2)%cells(dim2,i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3)) = &
													    cell_size(1)*(vorticity(3,i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3)) + &
													    (v_f%pr(dim)%cells(dim,i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3)) - &
													    v_f%pr(dim)%cells(dim,i-I_m(dim2,1)+I_m(dim,1),j-I_m(dim2,2)+I_m(dim,2),k-I_m(dim2,3)+I_m(dim,3)))/ &
													    cell_size(1) + v_f%pr(dim2)%cells(dim2,i,j,k)/cell_size(1))
												end if
											end if
									end select
								end if

								if (( bound_number1 /= 0 ).and.(bound_number3 /= 0)) then
									if ((boundary_type_name1=='wall').or.(boundary_type_name3=='wall')) boundary_type_name = 'wall'
									if ((boundary_type_name1=='inlet').or.(boundary_type_name3=='inlet')) boundary_type_name = 'inlet'

									select case(boundary_type_name)
										case('wall')
										case('inlet')
											continue
										case('outlet')
											if (sign < 0) then
												vorticity(3,i+I_m(dim2,1),j+I_m(dim2,2),k+I_m(dim2,3)) = 0.0_dp
												! 2.0*vorticity(3,i+I_m(dim,1)+I_m(dim2,1),j+I_m(dim,2)+I_m(dim2,2),k+I_m(dim,3)+I_m(dim2,3)) -
												! vorticity(3,i+2*I_m(dim,1)+I_m(dim2,1),j+2*I_m(dim,2)+I_m(dim2,2),k+2*I_m(dim,3)+I_m(dim2,3))
												! !vorticity(3,i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3))!2.0*vorticity(3,i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3)) -
												! vorticity(3,i+2*I_m(dim,1),j+2*I_m(dim,2),k+2*I_m(dim,3))
												! !0.5_dp*(vorticity(3,i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3)) +
												! vorticity(3,i+2*I_m(dim,1),j+2*I_m(dim,2),k+2*I_m(dim,3)))
												! !2.0*vorticity(3,i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3)) -
												! vorticity(3,i+2*I_m(dim,1),j+2*I_m(dim,2),k+2*I_m(dim,3))
												if (dim2 == 1) then
													v_f%pr(dim2)%cells(dim2,i-I_m(dim,1)+I_m(dim2,1),j-I_m(dim,2)+I_m(dim2,2),k-I_m(dim,3)+I_m(dim2,3)) = &
													    cell_size(1)*(vorticity(3,i+I_m(dim2,1),j+I_m(dim2,2),k+I_m(dim2,3)) - &
													    (v_f%pr(dim)%cells(dim,i+I_m(dim2,1),j+I_m(dim2,2),k+I_m(dim2,3)) - &
													    v_f%pr(dim)%cells(dim,i,j,k))/cell_size(1) + &
													    v_f%pr(dim2)%cells(dim2,i+I_m(dim2,1),j+I_m(dim2,2),k+I_m(dim2,3))/cell_size(1))
												else
													v_f%pr(dim2)%cells(dim2,i-I_m(dim,1)+I_m(dim2,1),j-I_m(dim,2)+I_m(dim2,2),k-I_m(dim,3)+I_m(dim2,3)) = &
													    -cell_size(1)*(vorticity(3,i+I_m(dim2,1),j+I_m(dim2,2),k+I_m(dim2,3)) + &
													    (v_f%pr(dim)%cells(dim,i+I_m(dim2,1),j+I_m(dim2,2),k+I_m(dim2,3)) - &
													    v_f%pr(dim)%cells(dim,i,j,k))/cell_size(1) - &
													    v_f%pr(dim2)%cells(dim2,i+I_m(dim2,1),j+I_m(dim2,2),k+I_m(dim2,3))/cell_size(1))
												end if
											else
												vorticity(3,i+I_m(dim,1)+I_m(dim2,1),j+I_m(dim,2)+I_m(dim2,2),k+I_m(dim,3)+I_m(dim2,3)) = 0.0_dp
												! 2.0*vorticity(3,i+I_m(dim2,1),j+I_m(dim2,2),k+I_m(dim2,3)) -
												! vorticity(3,i-I_m(dim,1)+I_m(dim2,1),j-I_m(dim,2)+I_m(dim2,2),k-I_m(dim,3)+I_m(dim2,3))
												! !vorticity(3,i,j,k)!2.0*vorticity(3,i,j,k) - vorticity(3,i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3))
												! !0.5_dp*(vorticity(3,i,j,k) + vorticity(3,i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3))) !2.0*vorticity(3,i,j,k) -
												! vorticity(3,i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3))
												if (dim2 == 1) then
													v_f%pr(dim2)%cells(dim2,i+I_m(dim,1)+I_m(dim2,1),j+I_m(dim,2)+I_m(dim2,2),k+I_m(dim,3)+I_m(dim2,3)) = &
													    -cell_size(1)*(vorticity(3,i+I_m(dim,1)+I_m(dim2,1),j+I_m(dim,2)+I_m(dim2,2),k+I_m(dim,3)+I_m(dim2,3)) - &
													    (v_f%pr(dim)%cells(dim,i+I_m(dim2,1)+I_m(dim,1),j+I_m(dim2,2)+I_m(dim,2),k+I_m(dim2,3)+I_m(dim,3)) - &
													    v_f%pr(dim)%cells(dim,i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3)))/cell_size(1) - &
													    v_f%pr(dim2)%cells(dim2,i+I_m(dim2,1),j+I_m(dim2,2),k+I_m(dim2,3))/cell_size(1))
												else
													v_f%pr(dim2)%cells(dim2,i+I_m(dim,1)+I_m(dim2,1),j+I_m(dim,2)+I_m(dim2,2),k+I_m(dim,3)+I_m(dim2,3)) = &
													    cell_size(1)*(vorticity(3,i+I_m(dim,1)+I_m(dim2,1),j+I_m(dim,2)+I_m(dim2,2),k+I_m(dim,3)+I_m(dim2,3)) + &
													    (v_f%pr(dim)%cells(dim,i+I_m(dim2,1)+I_m(dim,1),j+I_m(dim2,2)+I_m(dim,2),k+I_m(dim2,3)+I_m(dim,3)) - &
													    v_f%pr(dim)%cells(dim,i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3)))/cell_size(1) + &
													    v_f%pr(dim2)%cells(dim2,i+I_m(dim2,1),j+I_m(dim2,2),k+I_m(dim2,3))/cell_size(1))
												end if
											end if
									end select
								end if
							end if
						end do
					end do
				end if
			end do
			end do
			end do
			!$omp end do

			do dim = 1, dimensions
				loop(3,1) = cons_inner_loop(3,1)
				loop(3,2) = cons_utter_loop(3,2)*I_m(dim,3) + cons_inner_loop(3,2)*(1 - I_m(dim,3))

				loop(2,1) = cons_inner_loop(2,1)
				loop(2,2) = cons_utter_loop(2,2)*I_m(dim,2) + cons_inner_loop(2,2)*(1 - I_m(dim,2))

				loop(1,1) = cons_inner_loop(1,1)
				loop(1,2) = cons_utter_loop(1,2)*I_m(dim,1) + cons_inner_loop(1,2)*(1 - I_m(dim,1))

				!$omp do collapse(3) schedule(static)
				do k = loop(3,1),loop(3,2)
				do j = loop(2,1),loop(2,2)
				do i = loop(1,1),loop(1,2)
					if((bc%bc_markers(i,j,k) == 0).or.(bc%bc_markers(i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3)) == 0)) then
						if (predictor) then
							rho_face = 0.5_dp*(rho_old%cells(i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3)) + rho_old%cells(i,j,k))
							if (abs(this%buoyancy_reference_density - rho_face) > 1e-010_dp) then
								! Reduced gravity: ((rho-rho_ref)/rho) g on the face.
								F_a%cells(dim,i,j,k) = F_a%cells(dim,i,j,k) - ((this%buoyancy_reference_density - rho_face) * this%g(dim)) / rho_face
							end if

							if (this%viscosity_flag)		F_a%cells(dim,i,j,k)	=  F_a%cells(dim,i,j,k) - &
							    (1.0_dp/(0.5_dp*(rho_old%cells(i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3)) + rho_old%cells(i,j,k))) * &
							    (0.5_dp*(this%v_prod_visc%v_ptr%pr(dim)%cells(i,j,k) + &
							    this%v_prod_visc%v_ptr%pr(dim)%cells(i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3)))))

                            if (this%additional_particles_phases_number /= 0) then
							    do particles_phase_counter = 1, this%additional_particles_phases_number
								    ! Common contract: particle momentum coupling is gas acceleration.
								    F_a%cells(dim,i,j,k) = F_a%cells(dim,i,j,k) - 0.5_dp * ( &
								        this%v_prod_particles(particles_phase_counter)%v_ptr%pr(dim)%cells(i,j,k) + &
								        this%v_prod_particles(particles_phase_counter)%v_ptr%pr(dim)%cells( &
								        i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3)))
							    end do
						    end if
                        else
							rho_face = 0.5_dp*(rho_int%cells(i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3)) + rho_int%cells(i,j,k))
							if (abs(this%buoyancy_reference_density - rho_face) > 1e-010_dp) then
								! Reduced gravity: ((rho-rho_ref)/rho) g on the face.
								F_a%cells(dim,i,j,k) = F_a%cells(dim,i,j,k) - ((this%buoyancy_reference_density - rho_face) * this%g(dim)) / rho_face
							end if

							if (this%viscosity_flag)		F_a%cells(dim,i,j,k)	=  F_a%cells(dim,i,j,k) - &
							    (1.0_dp/(0.5_dp*(rho_int%cells(i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3)) + rho_int%cells(i,j,k))) * &
							    (0.5_dp*(this%v_prod_visc%v_ptr%pr(dim)%cells(i,j,k) + &
							    this%v_prod_visc%v_ptr%pr(dim)%cells(i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3)))))
                            if (this%additional_particles_phases_number /= 0) then
							    do particles_phase_counter = 1, this%additional_particles_phases_number
								    ! Corrector uses the same physical acceleration source.
								    F_a%cells(dim,i,j,k) = F_a%cells(dim,i,j,k) - 0.5_dp * ( &
								        this%v_prod_particles(particles_phase_counter)%v_ptr%pr(dim)%cells(i,j,k) + &
								        this%v_prod_particles(particles_phase_counter)%v_ptr%pr(dim)%cells( &
								        i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3)))
							    end do
						    end if
						end if
					end if
				end do
				end do
				end do
				!$omp end do
			end do

			do dim = 1, dimensions
				loop(3,1) = cons_inner_loop(3,1)
				loop(3,2) = cons_utter_loop(3,2)*I_m(dim,3) + cons_inner_loop(3,2)*(1 - I_m(dim,3))

				loop(2,1) = cons_inner_loop(2,1)
				loop(2,2) = cons_utter_loop(2,2)*I_m(dim,2) + cons_inner_loop(2,2)*(1 - I_m(dim,2))

				loop(1,1) = cons_inner_loop(1,1)
				loop(1,2) = cons_utter_loop(1,2)*I_m(dim,1) + cons_inner_loop(1,2)*(1 - I_m(dim,1))

				!$omp do collapse(3) schedule(static)
				do k = loop(3,1),loop(3,2)
				do j = loop(2,1),loop(2,2)
				do i = loop(1,1),loop(1,2)

					if((bc%bc_markers(i,j,k) == 0).or.(bc%bc_markers(i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3)) == 0)) then

						do dim1 = 1, 3
						do dim2 = 1, dimensions
								if ((dim1 /= dim).and.(dim2 /= dim1).and.(dim2 /= dim)) then
								if(((dim-dim1) == 1).or.((dim-dim1) == -2)) then
									F_a%cells(dim,i,j,k)	=  F_a%cells(dim,i,j,k) - 0.25_dp * (	vorticity(dim1,i,j,k)										* &
									    (v_f%pr(dim2)%cells(dim2,i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3)) + v_f%pr(dim2)%cells(dim2,i,j,k)) + &
																									vorticity(dim1,i+I_m(dim2,1),j+I_m(dim2,2),k+I_m(dim2,3))	* &
																									    (v_f%pr(dim2)%cells(dim2,i-I_m(dim,1)+I_m(dim2,1),j-I_m(dim,2)+I_m(dim2,2),k-I_m(dim,3)+I_m(dim2, &
																									    3)) + v_f%pr(dim2)%cells(dim2,i+I_m(dim2,1),j+I_m(dim2,2),k+I_m(dim2,3))))
								else
									F_a%cells(dim,i,j,k)	=  F_a%cells(dim,i,j,k) + 0.25_dp * (	vorticity(dim1,i,j,k)										* &
									    (v_f%pr(dim2)%cells(dim2,i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3)) + v_f%pr(dim2)%cells(dim2,i,j,k)) + &
																									vorticity(dim1,i+I_m(dim2,1),j+I_m(dim2,2),k+I_m(dim2,3))	* &
																									    (v_f%pr(dim2)%cells(dim2,i-I_m(dim,1)+I_m(dim2,1),j-I_m(dim,2)+I_m(dim2,2),k-I_m(dim,3)+I_m(dim2, &
																									    3)) + v_f%pr(dim2)%cells(dim2,i+I_m(dim2,1),j+I_m(dim2,2),k+I_m(dim2,3))))
								end if
								end if
						end do
						end do
					end if
				end do
				end do
				end do
				!$omp end do
			end do

			!				sign			= (-1)**plus
			!				bound_number	= bc%bc_markers(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3))
			!					boundary_type_name = bc%boundary_types(bound_number)%get_type_name()
			!					select case(boundary_type_name)
			!						case('outlet')
			!							!if (sign < 0) then
			!							!	F_a%cells(dim,i,j,k) = F_a%cells(dim,i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3))
			!							!else
			!							!	F_a%cells(dim,i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3)) = F_a%cells(dim,i,j,k)
			!							!end if
			!							!do dim2 = 1, dimensions
			!							!	if (dim /= dim2) then
			!							!		F_a%cells(dim2,i,j,k) = F_a%cells(dim,i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3))
			!							!	end if
			!							!end do
			!					end select

!$omp do collapse(3) schedule(static) reduction(+:grad_F_a_summ, r_summ2)
			do k = cons_inner_loop(3,1),cons_inner_loop(3,2)
			do j = cons_inner_loop(2,1),cons_inner_loop(2,2)
			do i = cons_inner_loop(1,1),cons_inner_loop(1,2)
				if(bc%bc_markers(i,j,k) == 0) then

					lame_coeffs		= 1.0_dp

					select case(coordinate_system)
						case ('cartesian')
							lame_coeffs			= 1.0_dp

                            r_summ2				= r_summ2 + 1.0_dp
						case ('cylindrical')
							! x -> r, y -> z
							lame_coeffs(1,1)	=  mesh%mesh(1,i,j,k) - 0.5_dp*cell_size(1)
							lame_coeffs(1,2)	=  mesh%mesh(1,i,j,k)
							lame_coeffs(1,3)	=  mesh%mesh(1,i,j,k) + 0.5_dp*cell_size(1)

                            r_summ2				= r_summ2 + mesh%mesh(1,i,j,k)
						case ('spherical')
							! x -> r
							lame_coeffs(1,1)	=  (mesh%mesh(1,i,j,k) - 0.5_dp*cell_size(1))**2
							lame_coeffs(1,2)	=  (mesh%mesh(1,i,j,k))**2
							lame_coeffs(1,3)	=  (mesh%mesh(1,i,j,k) + 0.5_dp*cell_size(1))**2

                            r_summ2				= r_summ2 + mesh%mesh(1,i,j,k)**2
					end select

					do dim = 1, dimensions
						grad_F_a(dim,i,j,k)	= (F_a%cells(dim,i+i_m(dim,1),j+i_m(dim,2),k+i_m(dim,3)) * lame_coeffs(dim,3) - F_a%cells(dim,i,j,k) &
						    * lame_coeffs(dim,1)) /  lame_coeffs(dim,2)
						grad_F_a_summ		= grad_F_a_summ + grad_F_a(dim,i,j,k) * lame_coeffs(1,2)
					end do

				end if
			end do
			end do
			end do
			!$omp end do


			if (this%all_Neumann_flag) then
			!$omp do collapse(3) schedule(static)
			do k = cons_inner_loop(3,1),cons_inner_loop(3,2)
			do j = cons_inner_loop(2,1),cons_inner_loop(2,2)
			do i = cons_inner_loop(1,1),cons_inner_loop(1,2)
				if(bc%bc_markers(i,j,k) == 0) then
					do dim = 1, dimensions
						grad_F_a(dim,i,j,k)  = grad_F_a(dim,i,j,k)  - grad_F_a_summ/r_summ2
					end do
				end if
			end do
			end do
			end do
			!$omp end do
			end if

			!$omp end parallel

            end associate
			pressure_iteration	= 0
			pressure_converged	= .false.

			do while ((.not.pressure_converged).and.(pressure_iteration < maximum_pressure_iterations))
				a_norm_init = 0.0_dp

				grad_F_b_summ = 0.0_dp

                r_summ1 = 0.0_dp

				associate ( p_dyn			=> this%p_dyn%s_ptr			, &
                            rho_old			=> this%rho_old%s_ptr		, &
                            rho_int			=> this%rho_int%s_ptr		, &
                            ddiv_v_dt		=> this%ddiv_v_dt%s_ptr		, &
                            F_a				=> this%F_a%s_ptr			, &
                            F_b				=> this%F_b%s_ptr			, &
                            H_old			=> this%H_old%s_ptr			, &
                            R				=> this%R%s_ptr				, &

                            grad_F_a		=> this%grad_F_a			, &
                            grad_F_b		=> this%grad_F_b			, &
                            v				=> this%v%v_ptr				, &
                            v_f				=> this%v_f%v_ptr			, &
                            v_f_old			=> this%v_f_old%v_ptr		, &

                            mesh			=> this%mesh%mesh_ptr		, &
                            bc				=> this%boundary%bc_ptr		)

				!$omp parallel default(shared) &
				!$omp& private(i,j,k,dim,dim2,loop,plus,sign,bound_number,boundary_type_name,lame_coeffs,farfield_velocity,rho_face)

				do dim = 1, dimensions
					loop(3,1) = cons_inner_loop(3,1)
					loop(3,2) = cons_utter_loop(3,2)*I_m(dim,3) + cons_inner_loop(3,2)*(1 - I_m(dim,3))

					loop(2,1) = cons_inner_loop(2,1)
					loop(2,2) = cons_utter_loop(2,2)*I_m(dim,2) + cons_inner_loop(2,2)*(1 - I_m(dim,2))

					loop(1,1) = cons_inner_loop(1,1)
					loop(1,2) = cons_utter_loop(1,2)*I_m(dim,1) + cons_inner_loop(1,2)*(1 - I_m(dim,1))

					!$omp do collapse(3) schedule(static)
					do k = loop(3,1),loop(3,2)
					do j = loop(2,1),loop(2,2)
					do i = loop(1,1),loop(1,2)

						if((bc%bc_markers(i,j,k) == 0).or.(bc%bc_markers(i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3)) == 0)) then
							if (predictor) then
								F_b%cells(dim,i,j,k)=	-	(p_dyn%cells(i,j,k)	*rho_old%cells(i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3))		&
														+	p_dyn%cells(i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3))	*rho_old%cells(i,j,k))		&
														/	(rho_old%cells(i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3))	+ rho_old%cells(i,j,k))	&
														*	(1.0_dp/rho_old%cells(i,j,k)	- 1.0_dp/rho_old%cells(i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3)))	/ cell_size(1)
														! F_b = -p_hat * grad(1/rho), [m/s^2]
							else
								F_b%cells(dim,i,j,k)=	-	(p_dyn%cells(i,j,k)	*rho_int%cells(i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3))		&
														+	p_dyn%cells(i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3))	*rho_int%cells(i,j,k))		&
														/	(rho_int%cells(i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3))	+ rho_int%cells(i,j,k))	&
														*	(1.0_dp/rho_int%cells(i,j,k)	- 1.0_dp/rho_int%cells(i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3)))	/ cell_size(1)
														! F_b = -p_hat * grad(1/rho), [m/s^2]
							end if

						end if
					end do
					end do
					end do
					!$omp end do

				end do

				!$omp do collapse(3) schedule(static) reduction(+:grad_F_b_summ,r_summ1)
				do k = cons_inner_loop(3,1),cons_inner_loop(3,2)
				do j = cons_inner_loop(2,1),cons_inner_loop(2,2)
				do i = cons_inner_loop(1,1),cons_inner_loop(1,2)
					if(bc%bc_markers(i,j,k) == 0) then

						lame_coeffs		= 1.0_dp

						select case(coordinate_system)
							case ('cartesian')
								lame_coeffs			= 1.0_dp

                                r_summ1				= r_summ1 + 1.0_dp
							case ('cylindrical')
								! x -> r, y -> z
								lame_coeffs(1,1)	=  mesh%mesh(1,i,j,k) - 0.5_dp*cell_size(1)
								lame_coeffs(1,2)	=  mesh%mesh(1,i,j,k)
								lame_coeffs(1,3)	=  mesh%mesh(1,i,j,k) + 0.5_dp*cell_size(1)

                                r_summ1				= r_summ1 + mesh%mesh(1,i,j,k)
							case ('spherical')
								! x -> r
								lame_coeffs(1,1)	=  (mesh%mesh(1,i,j,k) - 0.5_dp*cell_size(1))**2
								lame_coeffs(1,2)	=  (mesh%mesh(1,i,j,k))**2
								lame_coeffs(1,3)	=  (mesh%mesh(1,i,j,k) + 0.5_dp*cell_size(1))**2

                                r_summ1				= r_summ1 + mesh%mesh(1,i,j,k)**2
						end select

						do dim = 1, dimensions
							grad_F_b(dim,i,j,k)	= (F_b%cells(dim,i+i_m(dim,1),j+i_m(dim,2),k+i_m(dim,3)) * lame_coeffs(dim,3) - &
							    F_b%cells(dim,i,j,k) * lame_coeffs(dim,1)) /  lame_coeffs(dim,2)
							grad_F_b_summ	= grad_F_b_summ + grad_F_b(dim,i,j,k) * lame_coeffs(1,2)
						end do

					end if
				end do
				end do
				end do
				!$omp end do

				if (this%all_Neumann_flag) then
				!$omp do collapse(3) schedule(static)
				do k = cons_inner_loop(3,1),cons_inner_loop(3,2)
				do j = cons_inner_loop(2,1),cons_inner_loop(2,2)
				do i = cons_inner_loop(1,1),cons_inner_loop(1,2)
					if(bc%bc_markers(i,j,k) == 0) then
						do dim = 1, dimensions
							grad_F_b(dim,i,j,k)  = grad_F_b(dim,i,j,k)  - grad_F_b_summ/r_summ1
						end do
					end if
				end do
				end do
				end do
				!$omp end do
				end if

				!$omp do collapse(3) schedule(static) reduction(+:a_norm_init)
				do k = cons_inner_loop(3,1),cons_inner_loop(3,2)
				do j = cons_inner_loop(2,1),cons_inner_loop(2,2)
				do i = cons_inner_loop(1,1),cons_inner_loop(1,2)
					if(bc%bc_markers(i,j,k) == 0) then

						lame_coeffs		= 1.0_dp

						select case(coordinate_system)
							case ('cartesian')
								lame_coeffs			= 1.0_dp
							case ('cylindrical')
								! x -> r, y -> z
								lame_coeffs(1,1)	=  mesh%mesh(1,i,j,k) - 0.5_dp*cell_size(1)
								lame_coeffs(1,2)	=  mesh%mesh(1,i,j,k)
								lame_coeffs(1,3)	=  mesh%mesh(1,i,j,k) + 0.5_dp*cell_size(1)
							case ('spherical')
								! x -> r
								lame_coeffs(1,1)	=  (mesh%mesh(1,i,j,k) - 0.5_dp*cell_size(1))**2
								lame_coeffs(1,2)	=  (mesh%mesh(1,i,j,k))**2
								lame_coeffs(1,3)	=  (mesh%mesh(1,i,j,k) + 0.5_dp*cell_size(1))**2
						end select

						R%cells(i,j,k)	= 0.0_dp

						R%cells(i,j,k)	= R%cells(i,j,k) + cell_size(1)*cell_size(1)*ddiv_v_dt%cells(i,j,k)

						do dim = 1, dimensions

							R%cells(i,j,k)	= R%cells(i,j,k) +	(H_old%cells(i+i_m(dim,1),j+i_m(dim,2),k+i_m(dim,3)) - H_old%cells(i,j,k))* &
							    lame_coeffs(dim,3) / lame_coeffs(dim,2)

							R%cells(i,j,k)	= R%cells(i,j,k) -	(H_old%cells(i,j,k) - H_old%cells(i-i_m(dim,1),j-i_m(dim,2),k-i_m(dim,3)))* &
							    lame_coeffs(dim,1) / lame_coeffs(dim,2)

							R%cells(i,j,k)	= R%cells(i,j,k) +	cell_size(1)*grad_F_a(dim,i,j,k)

							R%cells(i,j,k)	= R%cells(i,j,k) +	cell_size(1)*grad_F_b(dim,i,j,k)

							do plus = 1,2
								sign			= (-1)**plus
								bound_number	= bc%bc_markers(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3))
								if( bound_number /= 0 ) then
									boundary_type_name = bc%boundary_types(bound_number)%get_type_name()
									select case(boundary_type_name)
										case('wall')
                                        if (predictor) then
                                            R%cells(i,j,k) = R%cells(i,j,k) + (H_old%cells(i,j,k)                             &
                                                - H_old%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3))          &
                                                - &
                                                    sign*(F_a%cells(dim,i+max(sign,0)*I_m(dim,1),j+max(sign,0)*I_m(dim,2),k+ &
                                                    max(sign,0)*I_m(dim,3)) &
                                                + &
                                                    F_b%cells(dim,i+max(sign,0)*I_m(dim,1),j+max(sign,0)*I_m(dim,2),k+ &
                                                    max(sign,0)*I_m(dim,3)))*cell_size(1) &
                                                - &
                                                    sign*(0.0_dp-v_f%pr(dim)%cells(dim,i+max(sign,0)*I_m(dim,1),j+max(sign,0)* &
                                                    I_m(dim,2),k+max(sign,0)*I_m(dim,3)))*cell_size(1)/time_step) &
                                                *lame_coeffs(dim,2+sign)/lame_coeffs(dim,2)
                                        else
                                            R%cells(i,j,k) = R%cells(i,j,k) + (H_old%cells(i,j,k)                             &
                                                - H_old%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3))          &
                                                - &
                                                    sign*(F_a%cells(dim,i+max(sign,0)*I_m(dim,1),j+max(sign,0)*I_m(dim,2),k+ &
                                                    max(sign,0)*I_m(dim,3)) &
                                                + &
                                                    F_b%cells(dim,i+max(sign,0)*I_m(dim,1),j+max(sign,0)*I_m(dim,2),k+ &
                                                    max(sign,0)*I_m(dim,3)))*cell_size(1) &
                                                - &
                                                    sign*(0.0_dp-0.5_dp*(v_f%pr(dim)%cells(dim,i+max(sign,0)*I_m(dim,1),j+ &
                                                    max(sign,0)*I_m(dim,2),k+max(sign,0)*I_m(dim,3)) &
                                                + &
                                                    v_f_old%pr(dim)%cells(dim,i+max(sign,0)*I_m(dim,1),j+max(sign,0)*I_m(dim, &
                                                    2),k+max(sign,0)*I_m(dim,3))))*cell_size(1)/(0.5_dp*time_step)) &
                                                *lame_coeffs(dim,2+sign)/lame_coeffs(dim,2)
                                        end if
                                        case('inlet')
											farfield_velocity = this%inlet_velocity


											if(predictor) then
												R%cells(i,j,k) = R%cells(i,j,k) + (H_old%cells(i,j,k) - &
												    H_old%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3)) &
																				- sign*(	F_a%cells(dim,i+max(sign,0)*I_m(dim,1),j+max(sign,0)*I_m(dim,2),k+max(sign,0)*I_m(dim,3)) &
																						+	F_b%cells(dim,i+max(sign,0)*I_m(dim,1),j+max(sign,0)*I_m(dim,2),k+max(sign,0)*I_m(dim,3))) * 1.0_dp * &
																						    cell_size(1) &
																				- sign*(	farfield_velocity - &
																				    v_f%pr(dim)%cells(dim,i+max(sign,0)*I_m(dim,1),j+max(sign,0)*I_m(dim,2),k+max(sign,0)*I_m(dim,3))) * &
																				    1.0_dp* cell_size(1)/time_step) * lame_coeffs(dim,2+sign) / lame_coeffs(dim,2)

											else
												R%cells(i,j,k) = R%cells(i,j,k) + (H_old%cells(i,j,k) - &
												    H_old%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3)) &
																				- sign*(	F_a%cells(dim,i+max(sign,0)*I_m(dim,1),j+max(sign,0)*I_m(dim,2),k+max(sign,0)*I_m(dim,3)) &
																						+	F_b%cells(dim,i+max(sign,0)*I_m(dim,1),j+max(sign,0)*I_m(dim,2),k+max(sign,0)*I_m(dim,3))) * 1.0_dp * &
																						    cell_size(1) &
																				- sign*(	farfield_velocity - &
																				    0.5_dp*(v_f%pr(dim)%cells(dim,i+max(sign,0)*I_m(dim,1),j+max(sign,0)*I_m(dim,2),k+max(sign,0)*I_m(dim, &
																				    3)) &
																						+ &
																						    v_f_old%pr(dim)%cells(dim,i+max(sign,0)*I_m(dim,1),j+max(sign,0)*I_m(dim,2),k+max(sign,0)*I_m(dim,3) &
																						    ))) * 1.0_dp * cell_size(1)/(0.5_dp*time_step)) * lame_coeffs(dim,2+sign) / lame_coeffs(dim,2)
											end if
										case('outlet')
                                        R%cells(i,j,k) = R%cells(i,j,k) - H_old%cells(i,j,k)                                 &
                                            - H_old%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3))
                                        do dim2 = 1, dimensions
                                            R%cells(i,j,k) = R%cells(i,j,k)                                                   &
                                                + v%pr(dim2)%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3))**2
                                        end do
									end select
								end if
							end do
						end do

						a_norm_init = a_norm_init + abs(R%cells(i,j,k)*lame_coeffs(1,2))

					end if
				end do
				end do
				end do
				!$omp end do

				!$omp end parallel

                end associate

				beta				= 2.0_dp/3.0_dp
				v_cycle_converged	= .false.

				nu_0 = 5

				tolerance = 1e-03_dp

                call fds_multigrid_timer%tic()

				v_cycle_iteration = 0
				do while (((v_cycle_iteration <= 0).or.(.not.v_cycle_converged)).and.(v_cycle_iteration <= nu_0))

					associate (     ddiv_v_dt		=> this%ddiv_v_dt%s_ptr		, &
                                    F_a				=> this%F_a%s_ptr			, &
                                    F_b				=> this%F_b%s_ptr			, &
                                    H				=> this%H%s_ptr				, &
                                    H_old			=> this%H_old%s_ptr			, &
                                    R				=> this%R%s_ptr				, &
                                    grad_F_a		=> this%grad_F_a			, &
                                    grad_F_b		=> this%grad_F_b			, &
                                    v				=> this%v%v_ptr				, &
                                    v_f				=> this%v_f%v_ptr			, &
                                    v_f_old			=> this%v_f_old%v_ptr		, &
                                    mesh			=> this%mesh%mesh_ptr		, &
                                    bc				=> this%boundary%bc_ptr)
					do legacy_smoothing_iteration = 1, this%pressure_legacy_pre_smoothing_sweeps

						a_norm	= 0.0_dp

                        !$omp parallel default(shared) &
                        !$omp& private(i,j,k,dim,dim2,plus,sign,bound_number,boundary_type_name,lame_coeffs,farfield_velocity)

						!$omp do collapse(3) schedule(static) reduction(+:a_norm)
						do k = cons_inner_loop(3,1),cons_inner_loop(3,2)
						do j = cons_inner_loop(2,1),cons_inner_loop(2,2)
						do i = cons_inner_loop(1,1),cons_inner_loop(1,2)
							if(bc%bc_markers(i,j,k) == 0) then

								lame_coeffs		= 1.0_dp

								select case(coordinate_system)
									case ('cartesian')
										lame_coeffs			= 1.0_dp
									case ('cylindrical')
										! x -> r, y -> z
										lame_coeffs(1,1)	=  mesh%mesh(1,i,j,k) - 0.5_dp*cell_size(1)
										lame_coeffs(1,2)	=  mesh%mesh(1,i,j,k)
										lame_coeffs(1,3)	=  mesh%mesh(1,i,j,k) + 0.5_dp*cell_size(1)
									case ('spherical')
										! x -> r
										lame_coeffs(1,1)	=  (mesh%mesh(1,i,j,k) - 0.5_dp*cell_size(1))**2
										lame_coeffs(1,2)	=  (mesh%mesh(1,i,j,k))**2
										lame_coeffs(1,3)	=  (mesh%mesh(1,i,j,k) + 0.5_dp*cell_size(1))**2
								end select

								R%cells(i,j,k) = 0.0_dp

								R%cells(i,j,k) = R%cells(i,j,k) + cell_size(1)*cell_size(1)*ddiv_v_dt%cells(i,j,k)

								do dim = 1, dimensions

									R%cells(i,j,k)	= R%cells(i,j,k) +	(H_old%cells(i+i_m(dim,1),j+i_m(dim,2),k+i_m(dim,3)) - H_old%cells(i,j,k))* &
									    lame_coeffs(dim,3) / lame_coeffs(dim,2)

									R%cells(i,j,k)	= R%cells(i,j,k) -	(H_old%cells(i,j,k) - H_old%cells(i-i_m(dim,1),j-i_m(dim,2),k-i_m(dim,3)))* &
									    lame_coeffs(dim,1) / lame_coeffs(dim,2)

									R%cells(i,j,k)	= R%cells(i,j,k) +	cell_size(1)*grad_F_a(dim,i,j,k)

									R%cells(i,j,k)	= R%cells(i,j,k) +	cell_size(1)*grad_F_b(dim,i,j,k)

									do plus = 1,2
										sign			= (-1)**plus
										bound_number	= bc%bc_markers(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3))
										if( bound_number /= 0 ) then
											boundary_type_name = bc%boundary_types(bound_number)%get_type_name()
											select case(boundary_type_name)
												case('wall')
                                                    if (predictor) then
                                                        R%cells(i,j,k) = R%cells(i,j,k) + (H_old%cells(i,j,k) &
                                                            - &
                                                                H_old%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign* &
                                                                I_m(dim,3)) &
                                                            - &
                                                                sign*(F_a%cells(dim,i+max(sign,0)*I_m(dim,1),j+max(sign,0)* &
                                                                I_m(dim,2),k+max(sign,0)*I_m(dim,3)) &
                                                            + &
                                                                F_b%cells(dim,i+max(sign,0)*I_m(dim,1),j+max(sign,0)*I_m(dim, &
                                                                2),k+max(sign,0)*I_m(dim,3)))*cell_size(1) &
                                                            - &
                                                                sign*(0.0_dp-v_f%pr(dim)%cells(dim,i+max(sign,0)*I_m(dim,1),j+ &
                                                                max(sign,0)*I_m(dim,2),k+max(sign,0)*I_m(dim,3)))*cell_size(1) &
                                                                /time_step) &
                                                            *lame_coeffs(dim,2+sign)/lame_coeffs(dim,2)
                                                    else
                                                        R%cells(i,j,k) = R%cells(i,j,k) + (H_old%cells(i,j,k) &
                                                            - &
                                                                H_old%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign* &
                                                                I_m(dim,3)) &
                                                            - &
                                                                sign*(F_a%cells(dim,i+max(sign,0)*I_m(dim,1),j+max(sign,0)* &
                                                                I_m(dim,2),k+max(sign,0)*I_m(dim,3)) &
                                                            + &
                                                                F_b%cells(dim,i+max(sign,0)*I_m(dim,1),j+max(sign,0)*I_m(dim, &
                                                                2),k+max(sign,0)*I_m(dim,3)))*cell_size(1) &
                                                            - &
                                                                sign*(0.0_dp-0.5_dp*(v_f%pr(dim)%cells(dim,i+max(sign,0)* &
                                                                I_m(dim,1),j+max(sign,0)*I_m(dim,2),k+max(sign,0)*I_m(dim,3)) &
                                                            + &
                                                                v_f_old%pr(dim)%cells(dim,i+max(sign,0)*I_m(dim,1),j+max(sign, &
                                                                0)*I_m(dim,2),k+max(sign,0)*I_m(dim,3))))*cell_size(1)/ &
                                                                (0.5_dp*time_step)) &
                                                            *lame_coeffs(dim,2+sign)/lame_coeffs(dim,2)
                                                    end if
                                                case('inlet')
													farfield_velocity = this%inlet_velocity


													if(predictor) then
														R%cells(i,j,k) = R%cells(i,j,k) + (H_old%cells(i,j,k) - &
														    H_old%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3)) &
																						- sign*(	F_a%cells(dim,i+max(sign,0)*I_m(dim,1),j+max(sign,0)*I_m(dim,2),k+max(sign,0)*I_m(dim,3)) &
																								+	F_b%cells(dim,i+max(sign,0)*I_m(dim,1),j+max(sign,0)*I_m(dim,2),k+max(sign,0)*I_m(dim,3))) * 1.0_dp &
																								    * cell_size(1) &
																						- sign*(	farfield_velocity - &
																						    v_f%pr(dim)%cells(dim,i+max(sign,0)*I_m(dim,1),j+max(sign,0)*I_m(dim,2),k+max(sign,0)*I_m(dim,3))) * &
																						    1.0_dp* cell_size(1)/time_step) * lame_coeffs(dim,2+sign) / lame_coeffs(dim,2)

													else
														R%cells(i,j,k) = R%cells(i,j,k) + (H_old%cells(i,j,k) - &
														    H_old%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3)) &
																						- sign*(	F_a%cells(dim,i+max(sign,0)*I_m(dim,1),j+max(sign,0)*I_m(dim,2),k+max(sign,0)*I_m(dim,3)) &
																								+	F_b%cells(dim,i+max(sign,0)*I_m(dim,1),j+max(sign,0)*I_m(dim,2),k+max(sign,0)*I_m(dim,3))) * 1.0_dp &
																								    * cell_size(1) &
																						- sign*(	farfield_velocity - &
																						    0.5_dp*(v_f%pr(dim)%cells(dim,i+max(sign,0)*I_m(dim,1),j+max(sign,0)*I_m(dim,2),k+max(sign,0)* &
																						    I_m(dim,3)) &
																								+ &
																								    v_f_old%pr(dim)%cells(dim,i+max(sign,0)*I_m(dim,1),j+max(sign,0)*I_m(dim,2),k+max(sign,0)*I_m(dim, &
																								    3)))) * 1.0_dp * cell_size(1)/(0.5_dp*time_step)) * lame_coeffs(dim,2+sign) / lame_coeffs(dim,2)
													end if
												case('outlet')
                                        R%cells(i,j,k) = R%cells(i,j,k) - H_old%cells(i,j,k)                                 &
                                            - H_old%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3))
                                        do dim2 = 1, dimensions
                                            R%cells(i,j,k) = R%cells(i,j,k)                                                   &
                                                + v%pr(dim2)%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3))**2
                                        end do
											end select
										end if
									end do
								end do

								H%cells(i,j,k)	= H_old%cells(i,j,k) + 1.0_dp/(2.0_dp*dimensions)*beta*R%cells(i,j,k)

								a_norm = a_norm + abs(R%cells(i,j,k)*lame_coeffs(1,2))
							end if
						end do
						end do
						end do
						!$omp end do

						!$omp do collapse(3) schedule(static)
						do k = cons_utter_loop(3,1),cons_utter_loop(3,2)
						do j = cons_utter_loop(2,1),cons_utter_loop(2,2)
						do i = cons_utter_loop(1,1),cons_utter_loop(1,2)
							! Keep R in the dimensionless FDS stencil form.  The conversion
							! to the shared elliptic RHS is performed unconditionally below.
							H_old%cells(i,j,k) = H%cells(i,j,k)
						end do
						end do
						end do
						!$omp end do
                        !$omp end parallel


	                    profile_legacy_smoothing_sweeps = profile_legacy_smoothing_sweeps + 1
					end do

					a_norm		= 0.0_dp
                    sub_F_summ	= 0.0_dp
                    r_summ1		= 0.0_dp

                    !$omp parallel default(shared) &
                    !$omp& private(i,j,k,dim,dim2,plus,sign,bound_number,boundary_type_name,lame_coeffs,farfield_velocity)

					!$omp do collapse(3) schedule(static) reduction(+:sub_F_summ,r_summ1)
					do k = cons_inner_loop(3,1),cons_inner_loop(3,2)
					do j = cons_inner_loop(2,1),cons_inner_loop(2,2)
					do i = cons_inner_loop(1,1),cons_inner_loop(1,2)
                        if(bc%bc_markers(i,j,k) == 0) then
							! R is dx**2 times the dimensional elliptic residual.  Scale it here
							! regardless of whether legacy pre-smoothing is enabled.
							this%pressure_residual(i,j,k) = R%cells(i,j,k) / cell_size(1) / cell_size(1)

							select case(coordinate_system)
								case ('cartesian')
									sub_F_summ	= sub_F_summ + this%pressure_residual(i,j,k)
									r_summ1 = r_summ1 + 1.0_dp
								case ('cylindrical')
									! x -> r, y -> z
									sub_F_summ	= sub_F_summ + this%pressure_residual(i,j,k) * mesh%mesh(1,i,j,k)
									r_summ1		= r_summ1 + mesh%mesh(1,i,j,k)
								case ('spherical')
									! x -> r
									sub_F_summ	= sub_F_summ + this%pressure_residual(i,j,k) * mesh%mesh(1,i,j,k)**2
									r_summ1		= r_summ1 + mesh%mesh(1,i,j,k)**2
							end select
                        end if

					end do
					end do
                    end do
					!$omp end do

                    if (this%all_Neumann_flag) then
                        !$omp do collapse(3) schedule(static)
					do k = cons_inner_loop(3,1),cons_inner_loop(3,2)
					do j = cons_inner_loop(2,1),cons_inner_loop(2,2)
					do i = cons_inner_loop(1,1),cons_inner_loop(1,2)
                        if(bc%bc_markers(i,j,k) == 0) then
							this%pressure_residual(i,j,k) = this%pressure_residual(i,j,k) - sub_F_summ/r_summ1
                        end if
					end do
					end do
                    end do
					!$omp end do
                    end if

					!$omp end parallel
                    end associate

					this%pressure_correction = 0.0_dp
                    call this%solve_shared_elliptic_correction(multigrid_statistics)
                    profile_multigrid_calls = profile_multigrid_calls + 1
                    profile_multigrid_cycles = profile_multigrid_cycles + multigrid_statistics%cycles
                    profile_multigrid_smoothing_iterations = profile_multigrid_smoothing_iterations + &
                        multigrid_statistics%smoothing_iterations
                    profile_multigrid_parallel_smoothing = profile_multigrid_parallel_smoothing + &
                        multigrid_statistics%parallel_smoothing_iterations
                    profile_multigrid_serial_smoothing = profile_multigrid_serial_smoothing + &
                        multigrid_statistics%serial_smoothing_iterations
                    profile_multigrid_residual_evaluations = profile_multigrid_residual_evaluations + &
                        multigrid_statistics%residual_evaluations
                    profile_max_hierarchy_levels = max(profile_max_hierarchy_levels, &
                        multigrid_statistics%hierarchy_levels)
                    profile_multigrid_initial_residual_max = max(profile_multigrid_initial_residual_max, &
                        multigrid_statistics%initial_residual)
                    profile_multigrid_final_residual_max = max(profile_multigrid_final_residual_max, &
                        multigrid_statistics%final_residual)
                    profile_multigrid_relative_residual_max = max(profile_multigrid_relative_residual_max, &
                        multigrid_statistics%relative_residual)

					associate (     ddiv_v_dt		=> this%ddiv_v_dt%s_ptr		, &
                                    F_a				=> this%F_a%s_ptr			, &
                                    F_b				=> this%F_b%s_ptr			, &
                                    H				=> this%H%s_ptr				, &
                                    H_old			=> this%H_old%s_ptr			, &
                                    R				=> this%R%s_ptr				, &
                                    grad_F_a		=> this%grad_F_a			, &
                                    grad_F_b		=> this%grad_F_b			, &
                                    v				=> this%v%v_ptr				, &
                                    v_f				=> this%v_f%v_ptr			, &
                                    v_f_old			=> this%v_f_old%v_ptr		, &
                                    mesh			=> this%mesh%mesh_ptr		, &
                                    bc				=> this%boundary%bc_ptr)
					!$omp parallel default(shared) &
					!$omp& private(i,j,k,dim,dim2,plus,sign,bound_number,boundary_type_name,lame_coeffs,farfield_velocity)

					!$omp do collapse(3) schedule(static)
					do k = cons_inner_loop(3,1),cons_inner_loop(3,2)
					do j = cons_inner_loop(2,1),cons_inner_loop(2,2)
					do i = cons_inner_loop(1,1),cons_inner_loop(1,2)
						! Apply the MG correction cumulatively and keep H/H_old synchronized.
						! With legacy pre-smoothing enabled H==H_old here, so this is
						! algebraically equivalent to the previous H_old = H + correction.
						H%cells(i,j,k) = H_old%cells(i,j,k) + this%pressure_correction(i,j,k)
						H_old%cells(i,j,k) = H%cells(i,j,k)
					end do
					end do
					end do
					!$omp end do

					!$omp do collapse(3) schedule(static) reduction(+:a_norm)
					do k = cons_inner_loop(3,1),cons_inner_loop(3,2)
					do j = cons_inner_loop(2,1),cons_inner_loop(2,2)
					do i = cons_inner_loop(1,1),cons_inner_loop(1,2)
						if(bc%bc_markers(i,j,k) == 0) then

							lame_coeffs		= 1.0_dp

							select case(coordinate_system)
								case ('cartesian')
									lame_coeffs			= 1.0_dp
								case ('cylindrical')
									! x -> r, y -> z
									lame_coeffs(1,1)	=  mesh%mesh(1,i,j,k) - 0.5_dp*cell_size(1)
									lame_coeffs(1,2)	=  mesh%mesh(1,i,j,k)
									lame_coeffs(1,3)	=  mesh%mesh(1,i,j,k) + 0.5_dp*cell_size(1)
								case ('spherical')
									! x -> r
									lame_coeffs(1,1)	=  (mesh%mesh(1,i,j,k) - 0.5_dp*cell_size(1))**2
									lame_coeffs(1,2)	=  (mesh%mesh(1,i,j,k))**2
									lame_coeffs(1,3)	=  (mesh%mesh(1,i,j,k) + 0.5_dp*cell_size(1))**2
							end select

							R%cells(i,j,k) = 0.0_dp

							R%cells(i,j,k) = R%cells(i,j,k) +  ddiv_v_dt%cells(i,j,k) * cell_size(1) * cell_size(1)

							do dim = 1, dimensions

								R%cells(i,j,k)	= R%cells(i,j,k) +	(H_old%cells(i+i_m(dim,1),j+i_m(dim,2),k+i_m(dim,3)) - H_old%cells(i,j,k))* &
								    lame_coeffs(dim,3) / lame_coeffs(dim,2)

								R%cells(i,j,k)	= R%cells(i,j,k) -	(H_old%cells(i,j,k) - H_old%cells(i-i_m(dim,1),j-i_m(dim,2),k-i_m(dim,3)))* &
								    lame_coeffs(dim,1) / lame_coeffs(dim,2)

								R%cells(i,j,k)	= R%cells(i,j,k) +	grad_F_a(dim,i,j,k) * cell_size(1)

								R%cells(i,j,k)	= R%cells(i,j,k) +	grad_F_b(dim,i,j,k) * cell_size(1)

								do plus = 1,2
									sign			= (-1)**plus
									bound_number	= bc%bc_markers(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3))
									if( bound_number /= 0 ) then
										boundary_type_name = bc%boundary_types(bound_number)%get_type_name()
										select case(boundary_type_name)
											case('wall')
                                        if (predictor) then
                                            R%cells(i,j,k) = R%cells(i,j,k) + (H_old%cells(i,j,k)                             &
                                                - H_old%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3))          &
                                                - &
                                                    sign*(F_a%cells(dim,i+max(sign,0)*I_m(dim,1),j+max(sign,0)*I_m(dim,2),k+ &
                                                    max(sign,0)*I_m(dim,3)) &
                                                + &
                                                    F_b%cells(dim,i+max(sign,0)*I_m(dim,1),j+max(sign,0)*I_m(dim,2),k+ &
                                                    max(sign,0)*I_m(dim,3)))*cell_size(1) &
                                                - &
                                                    sign*(0.0_dp-v_f%pr(dim)%cells(dim,i+max(sign,0)*I_m(dim,1),j+max(sign,0)* &
                                                    I_m(dim,2),k+max(sign,0)*I_m(dim,3)))*cell_size(1)/time_step) &
                                                *lame_coeffs(dim,2+sign)/lame_coeffs(dim,2)
                                        else
                                            R%cells(i,j,k) = R%cells(i,j,k) + (H_old%cells(i,j,k)                             &
                                                - H_old%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3))          &
                                                - &
                                                    sign*(F_a%cells(dim,i+max(sign,0)*I_m(dim,1),j+max(sign,0)*I_m(dim,2),k+ &
                                                    max(sign,0)*I_m(dim,3)) &
                                                + &
                                                    F_b%cells(dim,i+max(sign,0)*I_m(dim,1),j+max(sign,0)*I_m(dim,2),k+ &
                                                    max(sign,0)*I_m(dim,3)))*cell_size(1) &
                                                - &
                                                    sign*(0.0_dp-0.5_dp*(v_f%pr(dim)%cells(dim,i+max(sign,0)*I_m(dim,1),j+ &
                                                    max(sign,0)*I_m(dim,2),k+max(sign,0)*I_m(dim,3)) &
                                                + &
                                                    v_f_old%pr(dim)%cells(dim,i+max(sign,0)*I_m(dim,1),j+max(sign,0)*I_m(dim, &
                                                    2),k+max(sign,0)*I_m(dim,3))))*cell_size(1)/(0.5_dp*time_step)) &
                                                *lame_coeffs(dim,2+sign)/lame_coeffs(dim,2)
                                        end if
                                    case('inlet')
												farfield_velocity = this%inlet_velocity


												if(predictor) then
													R%cells(i,j,k) = R%cells(i,j,k) + (H_old%cells(i,j,k) - &
													    H_old%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3)) &
																					- sign*(	F_a%cells(dim,i+max(sign,0)*I_m(dim,1),j+max(sign,0)*I_m(dim,2),k+max(sign,0)*I_m(dim,3)) &
																							+	F_b%cells(dim,i+max(sign,0)*I_m(dim,1),j+max(sign,0)*I_m(dim,2),k+max(sign,0)*I_m(dim,3))) * 1.0_dp * &
																							    cell_size(1) &
																					- sign*(	farfield_velocity - &
																					    v_f%pr(dim)%cells(dim,i+max(sign,0)*I_m(dim,1),j+max(sign,0)*I_m(dim,2),k+max(sign,0)*I_m(dim,3))) * &
																					    1.0_dp* cell_size(1)/time_step) * lame_coeffs(dim,2+sign) / lame_coeffs(dim,2)

												else
													R%cells(i,j,k) = R%cells(i,j,k) + (H_old%cells(i,j,k) - &
													    H_old%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3)) &
																					- sign*(	F_a%cells(dim,i+max(sign,0)*I_m(dim,1),j+max(sign,0)*I_m(dim,2),k+max(sign,0)*I_m(dim,3)) &
																							+	F_b%cells(dim,i+max(sign,0)*I_m(dim,1),j+max(sign,0)*I_m(dim,2),k+max(sign,0)*I_m(dim,3))) * 1.0_dp * &
																							    cell_size(1) &
																					- sign*(	farfield_velocity - &
																					    0.5_dp*(v_f%pr(dim)%cells(dim,i+max(sign,0)*I_m(dim,1),j+max(sign,0)*I_m(dim,2),k+max(sign,0)* &
																					    I_m(dim,3)) &
																							+ &
																							    v_f_old%pr(dim)%cells(dim,i+max(sign,0)*I_m(dim,1),j+max(sign,0)*I_m(dim,2),k+max(sign,0)*I_m(dim, &
																							    3)))) * 1.0_dp * cell_size(1)/(0.5_dp*time_step)) * lame_coeffs(dim,2+sign) / lame_coeffs(dim,2)
												end if
											case('outlet')
                                        R%cells(i,j,k) = R%cells(i,j,k) - H_old%cells(i,j,k)                                 &
                                            - H_old%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3))
                                        do dim2 = 1, dimensions
                                            R%cells(i,j,k) = R%cells(i,j,k)                                                   &
                                                + v%pr(dim2)%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3))**2
                                        end do
										end select
									end if
								end do
							end do

							a_norm = a_norm + abs(R%cells(i,j,k)*lame_coeffs(1,2))

						end if
					end do
					end do
					end do
					!$omp end do
                    !$omp end parallel

					do legacy_smoothing_iteration = 1, this%pressure_legacy_post_smoothing_sweeps

						a_norm	= 0.0_dp

						!$omp parallel default(shared) &
						!$omp& private(i,j,k,dim,dim2,plus,sign,bound_number,boundary_type_name,lame_coeffs,farfield_velocity)

						!$omp do collapse(3) schedule(static) reduction(+:a_norm)
						do k = cons_inner_loop(3,1),cons_inner_loop(3,2)
						do j = cons_inner_loop(2,1),cons_inner_loop(2,2)
						do i = cons_inner_loop(1,1),cons_inner_loop(1,2)
							if(bc%bc_markers(i,j,k) == 0) then

								lame_coeffs		= 1.0_dp

								select case(coordinate_system)
									case ('cartesian')
										lame_coeffs			= 1.0_dp
									case ('cylindrical')
										! x -> r, y -> z
										lame_coeffs(1,1)	=  mesh%mesh(1,i,j,k) - 0.5_dp*cell_size(1)
										lame_coeffs(1,2)	=  mesh%mesh(1,i,j,k)
										lame_coeffs(1,3)	=  mesh%mesh(1,i,j,k) + 0.5_dp*cell_size(1)
									case ('spherical')
										! x -> r
										lame_coeffs(1,1)	=  (mesh%mesh(1,i,j,k) - 0.5_dp*cell_size(1))**2
										lame_coeffs(1,2)	=  (mesh%mesh(1,i,j,k))**2
										lame_coeffs(1,3)	=  (mesh%mesh(1,i,j,k) + 0.5_dp*cell_size(1))**2
								end select

								R%cells(i,j,k) = 0.0_dp

								R%cells(i,j,k) = R%cells(i,j,k) + cell_size(1)*cell_size(1)*ddiv_v_dt%cells(i,j,k)

								do dim = 1, dimensions

									R%cells(i,j,k)	= R%cells(i,j,k) +	(H_old%cells(i+i_m(dim,1),j+i_m(dim,2),k+i_m(dim,3)) - H_old%cells(i,j,k))* &
									    lame_coeffs(dim,3) / lame_coeffs(dim,2)

									R%cells(i,j,k)	= R%cells(i,j,k) -	(H_old%cells(i,j,k) - H_old%cells(i-i_m(dim,1),j-i_m(dim,2),k-i_m(dim,3)))* &
									    lame_coeffs(dim,1) / lame_coeffs(dim,2)

									R%cells(i,j,k)	= R%cells(i,j,k) +	cell_size(1)*grad_F_a(dim,i,j,k)

									R%cells(i,j,k)	= R%cells(i,j,k) +	cell_size(1)*grad_F_b(dim,i,j,k)

									do plus = 1,2
										sign			= (-1)**plus
										bound_number	= bc%bc_markers(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3))
										if( bound_number /= 0 ) then
											boundary_type_name = bc%boundary_types(bound_number)%get_type_name()
											select case(boundary_type_name)
												case('wall')
                                        if (predictor) then
                                            R%cells(i,j,k) = R%cells(i,j,k) + (H_old%cells(i,j,k)                             &
                                                - H_old%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3))          &
                                                - &
                                                    sign*(F_a%cells(dim,i+max(sign,0)*I_m(dim,1),j+max(sign,0)*I_m(dim,2),k+ &
                                                    max(sign,0)*I_m(dim,3)) &
                                                + &
                                                    F_b%cells(dim,i+max(sign,0)*I_m(dim,1),j+max(sign,0)*I_m(dim,2),k+ &
                                                    max(sign,0)*I_m(dim,3)))*cell_size(1) &
                                                - &
                                                    sign*(0.0_dp-v_f%pr(dim)%cells(dim,i+max(sign,0)*I_m(dim,1),j+max(sign,0)* &
                                                    I_m(dim,2),k+max(sign,0)*I_m(dim,3)))*cell_size(1)/time_step) &
                                                *lame_coeffs(dim,2+sign)/lame_coeffs(dim,2)
                                        else
                                            R%cells(i,j,k) = R%cells(i,j,k) + (H_old%cells(i,j,k)                             &
                                                - H_old%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3))          &
                                                - &
                                                    sign*(F_a%cells(dim,i+max(sign,0)*I_m(dim,1),j+max(sign,0)*I_m(dim,2),k+ &
                                                    max(sign,0)*I_m(dim,3)) &
                                                + &
                                                    F_b%cells(dim,i+max(sign,0)*I_m(dim,1),j+max(sign,0)*I_m(dim,2),k+ &
                                                    max(sign,0)*I_m(dim,3)))*cell_size(1) &
                                                - &
                                                    sign*(0.0_dp-0.5_dp*(v_f%pr(dim)%cells(dim,i+max(sign,0)*I_m(dim,1),j+ &
                                                    max(sign,0)*I_m(dim,2),k+max(sign,0)*I_m(dim,3)) &
                                                + &
                                                    v_f_old%pr(dim)%cells(dim,i+max(sign,0)*I_m(dim,1),j+max(sign,0)*I_m(dim, &
                                                    2),k+max(sign,0)*I_m(dim,3))))*cell_size(1)/(0.5_dp*time_step)) &
                                                *lame_coeffs(dim,2+sign)/lame_coeffs(dim,2)
                                        end if
                                    case('inlet')
													farfield_velocity = this%inlet_velocity


													if(predictor) then
														R%cells(i,j,k) = R%cells(i,j,k) + (H_old%cells(i,j,k) - &
														    H_old%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3)) &
																						- sign*(	F_a%cells(dim,i+max(sign,0)*I_m(dim,1),j+max(sign,0)*I_m(dim,2),k+max(sign,0)*I_m(dim,3)) &
																								+	F_b%cells(dim,i+max(sign,0)*I_m(dim,1),j+max(sign,0)*I_m(dim,2),k+max(sign,0)*I_m(dim,3))) * 1.0_dp &
																								    * cell_size(1) &
																						- sign*(	farfield_velocity - &
																						    v_f%pr(dim)%cells(dim,i+max(sign,0)*I_m(dim,1),j+max(sign,0)*I_m(dim,2),k+max(sign,0)*I_m(dim,3))) * &
																						    1.0_dp* cell_size(1)/time_step) * lame_coeffs(dim,2+sign) / lame_coeffs(dim,2)

													else
														R%cells(i,j,k) = R%cells(i,j,k) + (H_old%cells(i,j,k) - &
														    H_old%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3)) &
																						- sign*(	F_a%cells(dim,i+max(sign,0)*I_m(dim,1),j+max(sign,0)*I_m(dim,2),k+max(sign,0)*I_m(dim,3)) &
																								+	F_b%cells(dim,i+max(sign,0)*I_m(dim,1),j+max(sign,0)*I_m(dim,2),k+max(sign,0)*I_m(dim,3))) * 1.0_dp &
																								    * cell_size(1) &
																						- sign*(	farfield_velocity - &
																						    0.5_dp*(v_f%pr(dim)%cells(dim,i+max(sign,0)*I_m(dim,1),j+max(sign,0)*I_m(dim,2),k+max(sign,0)* &
																						    I_m(dim,3)) &
																								+ &
																								    v_f_old%pr(dim)%cells(dim,i+max(sign,0)*I_m(dim,1),j+max(sign,0)*I_m(dim,2),k+max(sign,0)*I_m(dim, &
																								    3)))) * 1.0_dp * cell_size(1)/(0.5_dp*time_step)) * lame_coeffs(dim,2+sign) / lame_coeffs(dim,2)
													end if
												case('outlet')
                                        R%cells(i,j,k) = R%cells(i,j,k) - H_old%cells(i,j,k)                                 &
                                            - H_old%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3))
														do dim2 = 1, dimensions
                                            R%cells(i,j,k) = R%cells(i,j,k)                                                   &
                                                + v%pr(dim2)%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3))**2
                                                														end do
											end select
										end if
									end do
								end do

								H%cells(i,j,k)	= H_old%cells(i,j,k) + 1.0_dp/(2.0_dp*dimensions)*beta*R%cells(i,j,k)

								a_norm = a_norm + abs(R%cells(i,j,k)*lame_coeffs(1,2))


							end if
						end do
						end do
						end do
						!$omp end do

						!$omp do collapse(3) schedule(static)
						do k = cons_utter_loop(3,1),cons_utter_loop(3,2)
						do j = cons_utter_loop(2,1),cons_utter_loop(2,2)
						do i = cons_utter_loop(1,1),cons_utter_loop(1,2)
							H_old%cells(i,j,k) = H%cells(i,j,k)
						end do
						end do
						end do
						!$omp end do


					!$omp end parallel

                        profile_legacy_smoothing_sweeps = profile_legacy_smoothing_sweeps + 1

                    end do
                    end associate

                    v_cycle_iteration = v_cycle_iteration + 1
                    this%pressure_correction = 0.0_dp
					if(abs(a_norm_init) > 1e-10_dp) then
						if((a_norm/a_norm_init < tolerance).or.(a_norm < 100e-0)) then
							v_cycle_converged = .true.
						end if
					end if
				end do

                call fds_multigrid_timer%toc(new_iter=.true.)

				associate (	    p_dyn			=> this%p_dyn%s_ptr			, &
					            F_a				=> this%F_a%s_ptr			, &
					            F_b				=> this%F_b%s_ptr			, &
					            H				=> this%H%s_ptr				, &
					            v				=> this%v%v_ptr				, &
					            v_f				=> this%v_f%v_ptr			, &
					            v_f_old			=> this%v_f_old%v_ptr		, &
					            mesh			=> this%mesh%mesh_ptr		, &
					            bc				=> this%boundary%bc_ptr		)

				do k = cons_inner_loop(3,1),cons_inner_loop(3,2)
				do j = cons_inner_loop(2,1),cons_inner_loop(2,2)
				do i = cons_inner_loop(1,1),cons_inner_loop(1,2)

					if(bc%bc_markers(i,j,k) == 0) then
						do dim = 1,dimensions
							do plus = 1,2
								sign			= (-1)**plus
								bound_number	= bc%bc_markers(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3))
								if( bound_number /= 0 ) then
									boundary_type_name = bc%boundary_types(bound_number)%get_type_name()
									select case(boundary_type_name)
										case('wall')
                                        ! Dirichlet normal velocity u_n=0. The former condition only cancelled
                                        ! acceleration and therefore preserved any spurious wall-face velocity.
                                        if (predictor) then
                                            H%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3)) = H%cells(i,j,k) &
                                                - &
                                                    sign*(F_a%cells(dim,i+max(sign,0)*I_m(dim,1),j+max(sign,0)*I_m(dim,2),k+ &
                                                    max(sign,0)*I_m(dim,3)) &
                                                + &
                                                    F_b%cells(dim,i+max(sign,0)*I_m(dim,1),j+max(sign,0)*I_m(dim,2),k+ &
                                                    max(sign,0)*I_m(dim,3)))*cell_size(1) &
                                                - &
                                                    sign*(0.0_dp-v_f%pr(dim)%cells(dim,i+max(sign,0)*I_m(dim,1),j+max(sign,0)* &
                                                    I_m(dim,2),k+max(sign,0)*I_m(dim,3)))*cell_size(1)/time_step
                                        else
                                            H%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3)) = H%cells(i,j,k) &
                                                - &
                                                    sign*(F_a%cells(dim,i+max(sign,0)*I_m(dim,1),j+max(sign,0)*I_m(dim,2),k+ &
                                                    max(sign,0)*I_m(dim,3)) &
                                                + &
                                                    F_b%cells(dim,i+max(sign,0)*I_m(dim,1),j+max(sign,0)*I_m(dim,2),k+ &
                                                    max(sign,0)*I_m(dim,3)))*cell_size(1) &
                                                - &
                                                    sign*(0.0_dp-0.5_dp*(v_f%pr(dim)%cells(dim,i+max(sign,0)*I_m(dim,1),j+ &
                                                    max(sign,0)*I_m(dim,2),k+max(sign,0)*I_m(dim,3)) &
                                                + &
                                                    v_f_old%pr(dim)%cells(dim,i+max(sign,0)*I_m(dim,1),j+max(sign,0)*I_m(dim, &
                                                    2),k+max(sign,0)*I_m(dim,3))))*cell_size(1)/(0.5_dp*time_step)
                                        end if
                                        case('outlet')
                                        ! Fixed perturbational pressure p'_b=0. Since H=p'/rho+|u|^2/2,
                                        ! use H_g=2H_b-H_i=-H_i+|u_b|^2 for both outflow and backflow.
                                            H%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3)) = -H%cells(i,j,k)
											do dim2 = 1, dimensions
                                                H%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3))  = &
                                                    H%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3)) &
                                                                                                                + &
                                                                                                                    v%pr(dim2) &
                                                                                                                    %cells(i+ &
                                                                                                                    sign* &
                                                                                                                    I_m(dim,1) &
                                                                                                                    ,j+sign* &
                                                                                                                    I_m(dim,2) &
                                                                                                                    ,k+sign* &
                                                                                                                    I_m(dim,3) &
                                                                                                                    )**2
                                            end do

                                        case('inlet')
											farfield_velocity = this%inlet_velocity


											if (predictor) then
												H%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3))	= H%cells(i,j,k)														&
																												- sign*(	F_a%cells(dim,i+max(sign,0)*I_m(dim,1),j+max(sign,0)*I_m(dim,2),k+max(sign,0)*I_m(dim,3)) &
																														+	F_b%cells(dim,i+max(sign,0)*I_m(dim,1),j+max(sign,0)*I_m(dim,2),k+max(sign,0)*I_m(dim,3))) * &
																														    1.0_dp * cell_size(1) &
																												- sign*(	farfield_velocity - &
																												    v_f%pr(dim)%cells(dim,i+max(sign,0)*I_m(dim,1),j+max(sign,0)*I_m(dim,2),k+max(sign,0)*I_m(dim, &
																												    3))) * 1.0_dp * cell_size(1)/time_step
											else
												H%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3))	= H%cells(i,j,k)														&
																												- sign*(	F_a%cells(dim,i+max(sign,0)*I_m(dim,1),j+max(sign,0)*I_m(dim,2),k+max(sign,0)*I_m(dim,3)) &
																														+	F_b%cells(dim,i+max(sign,0)*I_m(dim,1),j+max(sign,0)*I_m(dim,2),k+max(sign,0)*I_m(dim,3))) * &
																														    1.0_dp * cell_size(1) &
																												- sign*(	farfield_velocity - &
																												    0.5_dp*(v_f%pr(dim)%cells(dim,i+max(sign,0)*I_m(dim,1),j+max(sign,0)*I_m(dim,2),k+max(sign,0)* &
																												    I_m(dim,3)) &
																														+ &
																														    v_f_old%pr(dim)%cells(dim,i+max(sign,0)*I_m(dim,1),j+max(sign,0)*I_m(dim,2),k+max(sign,0)* &
																														    I_m(dim,3)))) * 1.0_dp * cell_size(1)/(0.5_dp*time_step)
											end if
									end select
								end if
							end do
						end do
					end if
				end do
				end do
				end do

                end associate

				call this%calculate_dynamic_pressure(time_step,predictor)

				! Direct fixed-point convergence test for the baroclinic pressure term.
				! F_b currently contains F_b[p^m], which was used in the Poisson solve.
				! Recompute F_b[p^(m+1)] from the updated dynamic pressure, measure
				! its maximum change, and retain the updated field for the next outer
				! iteration (or for the final velocity update after convergence).
				F_b_change_max_local = 0.0_dp
				F_b_norm_max_local   = 0.0_dp
				F_a_norm_max_local   = 0.0_dp

                associate ( p_dyn            => this%p_dyn%s_ptr          , &
                            rho_old          => this%rho_old%s_ptr        , &
                            rho_int          => this%rho_int%s_ptr        , &
                            F_a              => this%F_a%s_ptr            , &
                            F_b              => this%F_b%s_ptr            , &
                            bc               => this%boundary%bc_ptr      )

                !$omp parallel default(shared) private(i,j,k,dim,loop,F_b_candidate) &
                !$omp& reduction(max:F_b_change_max_local,F_b_norm_max_local,F_a_norm_max_local)

                do dim = 1, dimensions
                    loop(3,1) = cons_inner_loop(3,1)
                    loop(3,2) = cons_utter_loop(3,2)*I_m(dim,3) + cons_inner_loop(3,2)*(1 - I_m(dim,3))

                    loop(2,1) = cons_inner_loop(2,1)
                    loop(2,2) = cons_utter_loop(2,2)*I_m(dim,2) + cons_inner_loop(2,2)*(1 - I_m(dim,2))

                    loop(1,1) = cons_inner_loop(1,1)
                    loop(1,2) = cons_utter_loop(1,2)*I_m(dim,1) + cons_inner_loop(1,2)*(1 - I_m(dim,1))

                    !$omp do collapse(3) schedule(static)
                    do k = loop(3,1),loop(3,2)
                    do j = loop(2,1),loop(2,2)
                    do i = loop(1,1),loop(1,2)
                        if ((bc%bc_markers(i,j,k) == 0).or. &
                            (bc%bc_markers(i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3)) == 0)) then

                            if (predictor) then
                                F_b_candidate = - (p_dyn%cells(i,j,k) * &
                                    rho_old%cells(i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3)) + &
                                    p_dyn%cells(i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3)) * &
                                    rho_old%cells(i,j,k)) / &
                                    (rho_old%cells(i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3)) + &
                                     rho_old%cells(i,j,k)) * &
                                    (1.0_dp/rho_old%cells(i,j,k) - &
                                     1.0_dp/rho_old%cells(i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3))) / &
                                    cell_size(1)
                            else
                                F_b_candidate = - (p_dyn%cells(i,j,k) * &
                                    rho_int%cells(i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3)) + &
                                    p_dyn%cells(i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3)) * &
                                    rho_int%cells(i,j,k)) / &
                                    (rho_int%cells(i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3)) + &
                                     rho_int%cells(i,j,k)) * &
                                    (1.0_dp/rho_int%cells(i,j,k) - &
                                     1.0_dp/rho_int%cells(i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3))) / &
                                    cell_size(1)
                            end if

                            F_b_change_max_local = max(F_b_change_max_local, &
                                                       abs(F_b_candidate - F_b%cells(dim,i,j,k)))
                            F_b_norm_max_local   = max(F_b_norm_max_local, abs(F_b_candidate))
                            F_a_norm_max_local   = max(F_a_norm_max_local, abs(F_a%cells(dim,i,j,k)))
                            F_b%cells(dim,i,j,k) = F_b_candidate
                        end if
                    end do
                    end do
                    end do
                    !$omp end do
                end do

                !$omp end parallel

                end associate

				F_b_change_max = F_b_change_max_local
				F_b_norm_max   = F_b_norm_max_local
				F_a_norm_max   = F_a_norm_max_local

				F_b_convergence_scale = max(maxval(abs(this%g)), F_b_norm_max, F_a_norm_max, tiny(1.0_dp))
				pressure_converged = F_b_change_max <= F_b_absolute_tolerance + &
				                     F_b_relative_tolerance * F_b_convergence_scale

				pressure_iteration = pressure_iteration + 1

			end do

            call this%write_pressure_profile(predictor, pressure_converged, pressure_iteration, &
                profile_legacy_smoothing_sweeps, profile_multigrid_calls, profile_multigrid_cycles, &
                profile_multigrid_smoothing_iterations, profile_multigrid_parallel_smoothing, &
                profile_multigrid_serial_smoothing, profile_multigrid_residual_evaluations, &
                profile_max_hierarchy_levels, profile_multigrid_initial_residual_max, &
                profile_multigrid_final_residual_max, profile_multigrid_relative_residual_max)

			if (.not.pressure_converged) then
				print *, 'WARNING: F_b pressure iteration did not converge within the configured iteration limit.'
			end if
    end subroutine calculate_pressure_poisson

    !--------------------------------------------------------------------------
    ! Translate the current finite-difference defect into the cell-volume and
    ! face-conductance contract of elliptic_multigrid_solver_class, solve the
    ! homogeneous correction problem, and return the finest-grid correction.
    !--------------------------------------------------------------------------
    subroutine solve_shared_elliptic_correction(this, statistics)
        class(fds_solver), intent(inout) :: this
        type(elliptic_solver_statistics), intent(out) :: statistics

        integer, dimension(3,2) :: loop
        real(dp), dimension(3) :: dx
        integer :: dimensions, geometry_power
        integer :: i, j, k, ii, jj, kk
        real(dp) :: base_volume
        character(len=20) :: coordinate_system

        dimensions = this%domain%get_domain_dimensions()
        loop = this%domain%get_local_inner_cells_bounds()
        dx = this%mesh%mesh_ptr%get_cell_edges_length()
        base_volume = this%mesh%mesh_ptr%get_cell_volume()
        coordinate_system = this%domain%get_coordinate_system_name()
        select case (coordinate_system)
        case ('cartesian')
            geometry_power = 0
        case ('cylindrical')
            geometry_power = 1
        case ('spherical')
            geometry_power = 2
        case default
            error stop 'FDS multigrid adapter: unsupported coordinate system'
        end select

        ! Cheap scalar/shape guards catch mesh or geometry changes.  Boundary
        ! topology in FDS is assumed static after setup; any future run-time BC
        ! topology mutation must call invalidate_pressure_operator_cache().
        if (this%pressure_operator_prepared) then
            if (dimensions /= this%pressure_operator_dimensions .or. &
                geometry_power /= this%pressure_operator_geometry_power .or. &
                any(loop /= this%pressure_operator_loop) .or. &
                any(dx /= this%pressure_operator_dx) .or. &
                base_volume /= this%pressure_operator_base_volume) then
                call this%invalidate_pressure_operator_cache()
            end if
        end if
        if (.not. this%pressure_operator_prepared) call this%prepare_shared_elliptic_operator()

        ! Only the finest RHS changes between homogeneous pressure corrections.
        !$omp parallel do collapse(3) schedule(static) private(i,j,k,ii,jj,kk)
        do k = loop(3,1), loop(3,2)
            do j = loop(2,1), loop(2,2)
                do i = loop(1,1), loop(1,2)
                    ii = i-loop(1,1)+1
                    jj = j-loop(2,1)+1
                    kk = k-loop(3,1)+1
                    this%pressure_adapter_rhs(ii,jj,kk) = this%pressure_residual(i,j,k)
                end do
            end do
        end do
        !$omp end parallel do

        this%pressure_solver%relative_tolerance = 1.0e-3_dp
        this%pressure_solver%absolute_tolerance = 1.0e-10_dp
        call this%pressure_solver%solve_prepared(this%pressure_adapter_correction, &
            this%pressure_adapter_rhs, use_initial_guess=.false., statistics=statistics)
        if (.not. statistics%converged) then
            error stop 'FDS shared prepared multigrid correction did not converge'
        end if

        !$omp parallel do collapse(3) schedule(static) private(i,j,k,ii,jj,kk)
        do k = loop(3,1), loop(3,2)
            do j = loop(2,1), loop(2,2)
                do i = loop(1,1), loop(1,2)
                    ii = i-loop(1,1)+1
                    jj = j-loop(2,1)+1
                    kk = k-loop(3,1)+1
                    this%pressure_correction(i,j,k) = this%pressure_adapter_correction(ii,jj,kk)
                end do
            end do
        end do
        !$omp end parallel do
    end subroutine solve_shared_elliptic_correction


    !--------------------------------------------------------------------------
    ! Build the FDS pressure-correction operator once and ask the shared solver
    ! to retain its fully coarsened/preassembled hierarchy.  Conductances,
    ! volumes, active masks and boundary topology are static for an ordinary FDS
    ! run, while the correction RHS changes hundreds of times.
    !--------------------------------------------------------------------------
    subroutine prepare_shared_elliptic_operator(this)
        class(fds_solver), intent(inout) :: this

        type(elliptic_boundary_data) :: elliptic_boundary
        integer, dimension(3,2) :: loop
        real(dp), dimension(3) :: dx
        real(dp), allocatable :: cell_volume(:,:,:)
        real(dp), allocatable :: conductance_x(:,:,:), conductance_y(:,:,:)
        real(dp), allocatable :: conductance_z(:,:,:)
        logical, allocatable :: active(:,:,:)
        integer :: dimensions, geometry_power
        integer :: i, j, k, ii, jj, kk
        integer :: left_marker, right_marker, boundary_marker
        real(dp) :: base_volume, radius, radius_face
        real(dp) :: metric_center, metric_face
        character(len=20) :: boundary_name, coordinate_system

        dimensions = this%domain%get_domain_dimensions()
        loop = this%domain%get_local_inner_cells_bounds()
        dx = this%mesh%mesh_ptr%get_cell_edges_length()
        base_volume = this%mesh%mesh_ptr%get_cell_volume()
        coordinate_system = this%domain%get_coordinate_system_name()
        select case (coordinate_system)
        case ('cartesian')
            geometry_power = 0
        case ('cylindrical')
            geometry_power = 1
        case ('spherical')
            geometry_power = 2
        case default
            error stop 'FDS multigrid adapter: unsupported coordinate system'
        end select

        if (allocated(this%pressure_adapter_rhs)) deallocate(this%pressure_adapter_rhs)
        if (allocated(this%pressure_adapter_correction)) deallocate(this%pressure_adapter_correction)
        allocate(this%pressure_adapter_rhs(loop(1,2)-loop(1,1)+1, &
            loop(2,2)-loop(2,1)+1, loop(3,2)-loop(3,1)+1))
        allocate(this%pressure_adapter_correction, mold=this%pressure_adapter_rhs)
        this%pressure_adapter_rhs = 0.0_dp
        this%pressure_adapter_correction = 0.0_dp

        allocate(cell_volume, mold=this%pressure_adapter_rhs)
        allocate(active(size(this%pressure_adapter_rhs,1), size(this%pressure_adapter_rhs,2), &
            size(this%pressure_adapter_rhs,3)))
        allocate(conductance_x(size(active,1)+1,size(active,2),size(active,3)))
        allocate(conductance_y(size(active,1),size(active,2)+1,size(active,3)))
        allocate(conductance_z(size(active,1),size(active,2),size(active,3)+1))
        allocate(elliptic_boundary%type_x(size(conductance_x,1),size(conductance_x,2),size(conductance_x,3)))
        allocate(elliptic_boundary%value_x(size(conductance_x,1),size(conductance_x,2),size(conductance_x,3)))
        allocate(elliptic_boundary%type_y(size(conductance_y,1),size(conductance_y,2),size(conductance_y,3)))
        allocate(elliptic_boundary%value_y(size(conductance_y,1),size(conductance_y,2),size(conductance_y,3)))
        allocate(elliptic_boundary%type_z(size(conductance_z,1),size(conductance_z,2),size(conductance_z,3)))
        allocate(elliptic_boundary%value_z(size(conductance_z,1),size(conductance_z,2),size(conductance_z,3)))

        cell_volume = base_volume
        active = .false.
        conductance_x = 0.0_dp
        conductance_y = 0.0_dp
        conductance_z = 0.0_dp
        elliptic_boundary%type_x = elliptic_bc_internal
        elliptic_boundary%type_y = elliptic_bc_internal
        elliptic_boundary%type_z = elliptic_bc_internal
        elliptic_boundary%value_x = 0.0_dp
        elliptic_boundary%value_y = 0.0_dp
        elliptic_boundary%value_z = 0.0_dp

        associate (bc => this%boundary%bc_ptr, mesh => this%mesh%mesh_ptr)
        do k = loop(3,1), loop(3,2)
            do j = loop(2,1), loop(2,2)
                do i = loop(1,1), loop(1,2)
                    ii = i-loop(1,1)+1
                    jj = j-loop(2,1)+1
                    kk = k-loop(3,1)+1
                    if (bc%bc_markers(i,j,k) /= 0) cycle
                    active(ii,jj,kk) = .true.
                    radius = mesh%mesh(1,i,j,k)
                    metric_center = 1.0_dp
                    if (geometry_power > 0) metric_center = max(radius,0.0_dp)**geometry_power
                    cell_volume(ii,jj,kk) = base_volume*metric_center
                end do
            end do
        end do

        do kk = 1, size(active,3)
            do jj = 1, size(active,2)
                do ii = 1, size(active,1)+1
                    i = loop(1,1)+ii-1
                    j = loop(2,1)+jj-1
                    k = loop(3,1)+kk-1
                    radius_face = mesh%mesh(1,loop(1,1),j,k) - 0.5_dp*dx(1) + real(ii-1,dp)*dx(1)
                    metric_face = 1.0_dp
                    if (geometry_power > 0) metric_face = max(radius_face,0.0_dp)**geometry_power
                    conductance_x(ii,jj,kk) = base_volume*metric_face/dx(1)**2
                    left_marker = bc%bc_markers(i-1,j,k)
                    right_marker = bc%bc_markers(i,j,k)
                    call classify_correction_boundary(left_marker, right_marker, &
                        elliptic_boundary%type_x(ii,jj,kk))
                end do
            end do
        end do

        if (dimensions >= 2) then
            do kk = 1, size(active,3)
                do jj = 1, size(active,2)+1
                    do ii = 1, size(active,1)
                        i = loop(1,1)+ii-1
                        j = loop(2,1)+jj-1
                        k = loop(3,1)+kk-1
                        radius = mesh%mesh(1,i,loop(2,1),k)
                        metric_center = 1.0_dp
                        if (geometry_power > 0) metric_center = max(radius,0.0_dp)**geometry_power
                        conductance_y(ii,jj,kk) = base_volume*metric_center/dx(2)**2
                        left_marker = bc%bc_markers(i,j-1,k)
                        right_marker = bc%bc_markers(i,j,k)
                        call classify_correction_boundary(left_marker, right_marker, &
                            elliptic_boundary%type_y(ii,jj,kk))
                    end do
                end do
            end do
        end if

        if (dimensions >= 3) then
            do kk = 1, size(active,3)+1
                do jj = 1, size(active,2)
                    do ii = 1, size(active,1)
                        i = loop(1,1)+ii-1
                        j = loop(2,1)+jj-1
                        k = loop(3,1)+kk-1
                        radius = mesh%mesh(1,i,j,loop(3,1))
                        metric_center = 1.0_dp
                        if (geometry_power > 0) metric_center = max(radius,0.0_dp)**geometry_power
                        conductance_z(ii,jj,kk) = base_volume*metric_center/dx(3)**2
                        left_marker = bc%bc_markers(i,j,k-1)
                        right_marker = bc%bc_markers(i,j,k)
                        call classify_correction_boundary(left_marker, right_marker, &
                            elliptic_boundary%type_z(ii,jj,kk))
                    end do
                end do
            end do
        end if
        end associate

        call this%pressure_solver%prepare(conductance_x, conductance_y, conductance_z, &
            cell_volume, active, elliptic_boundary, dimensions)
        this%pressure_operator_prepared = .true.
        this%pressure_operator_preparations = this%pressure_operator_preparations + 1
        this%pressure_operator_dimensions = dimensions
        this%pressure_operator_geometry_power = geometry_power
        this%pressure_operator_loop = loop
        this%pressure_operator_dx = dx
        this%pressure_operator_base_volume = base_volume

    contains
        subroutine classify_correction_boundary(marker_a, marker_b, boundary_type)
            integer, intent(in) :: marker_a, marker_b
            integer, intent(out) :: boundary_type

            boundary_type = elliptic_bc_internal
            if (marker_a == 0 .and. marker_b == 0) return
            boundary_marker = merge(marker_a, marker_b, marker_a /= 0)
            boundary_name = this%boundary%bc_ptr%boundary_types(boundary_marker)%get_type_name()
            if (boundary_name == 'outlet') then
                boundary_type = elliptic_bc_dirichlet
            else
                boundary_type = elliptic_bc_neumann
            end if
        end subroutine classify_correction_boundary
    end subroutine prepare_shared_elliptic_operator


    !> Explicitly invalidate the cached pressure operator after any topology or
    !! metric change not detectable by the inexpensive scalar guards above.
    subroutine invalidate_pressure_operator_cache(this)
        class(fds_solver), intent(inout) :: this
        this%pressure_operator_prepared = .false.
        call this%pressure_solver%invalidate()
    end subroutine invalidate_pressure_operator_cache


    !--------------------------------------------------------------------------
    ! Emit one compact pressure-performance record for each predictor/corrector
    ! projection.  The file is intentionally independent of the generic
    ! elliptic solver so CABARET can reuse that solver without FDS diagnostics.
    !--------------------------------------------------------------------------
    subroutine write_pressure_profile(this, predictor, pressure_converged, pressure_iterations, &
        legacy_smoothing_sweeps, multigrid_calls, multigrid_cycles, &
        multigrid_smoothing_iterations, multigrid_parallel_smoothing, &
        multigrid_serial_smoothing, multigrid_residual_evaluations, hierarchy_levels, &
        multigrid_initial_residual_max, multigrid_final_residual_max, &
        multigrid_relative_residual_max)

        class(fds_solver), intent(inout) :: this
        logical, intent(in) :: predictor, pressure_converged
        integer, intent(in) :: pressure_iterations, legacy_smoothing_sweeps
        integer, intent(in) :: multigrid_calls, multigrid_cycles, multigrid_smoothing_iterations
        integer, intent(in) :: multigrid_parallel_smoothing, multigrid_serial_smoothing
        integer, intent(in) :: multigrid_residual_evaluations, hierarchy_levels
        real(dp), intent(in) :: multigrid_initial_residual_max
        real(dp), intent(in) :: multigrid_final_residual_max
        real(dp), intent(in) :: multigrid_relative_residual_max

        character(len=9) :: stage

        if (.not. this%pressure_profiling_enabled) return

        if (this%pressure_profile_unit == 0) then
            open(newunit=this%pressure_profile_unit, file='fds_pressure_profile.csv', &
                status='replace', action='write')
            write(this%pressure_profile_unit,'(A)') &
                'time,stage,pressure_converged,nonlinear_pressure_iterations,'// &
                'legacy_pre_sweeps_per_correction,legacy_post_sweeps_per_correction,'// &
                'legacy_smoothing_sweeps,multigrid_calls,multigrid_cycles,'// &
                'multigrid_smoothing_iterations,multigrid_parallel_smoothing_iterations,'// &
                'multigrid_serial_smoothing_iterations,multigrid_residual_evaluations,'// &
                'multigrid_hierarchy_levels,prepared_operator_rebuilds,multigrid_initial_residual_max,'// &
                'multigrid_final_residual_max,multigrid_relative_residual_max'
        end if

        if (predictor) then
            stage = 'predictor'
        else
            stage = 'corrector'
        end if

        write(this%pressure_profile_unit, &
            '(ES24.16,",",A,",",L1,12(",",I0),3(",",ES24.16))') &
            this%time, trim(stage), pressure_converged, pressure_iterations, &
            this%pressure_legacy_pre_smoothing_sweeps, &
            this%pressure_legacy_post_smoothing_sweeps, &
            legacy_smoothing_sweeps, multigrid_calls, multigrid_cycles, &
            multigrid_smoothing_iterations, multigrid_parallel_smoothing, &
            multigrid_serial_smoothing, multigrid_residual_evaluations, hierarchy_levels, &
            this%pressure_operator_preparations, multigrid_initial_residual_max, &
            multigrid_final_residual_max, &
            multigrid_relative_residual_max
        flush(this%pressure_profile_unit)
    end subroutine write_pressure_profile

    !--------------------------------------------------------------------------
    ! Recover perturbational pressure from H = p'/rho + |u|^2/2 and populate
    ! ghost values consistent with wall/inlet velocity conditions and p'=0 at
    ! an outlet.
    !--------------------------------------------------------------------------
    subroutine calculate_dynamic_pressure(this,time_step,predictor)
        class(fds_solver)    ,intent(inout)    :: this
        real(dp)            ,intent(in)        :: time_step
        logical                ,intent(in)        :: predictor

        real(dp)    :: farfield_density, farfield_pressure, farfield_velocity
        real(dp)    :: vel_abs

        real(dp)    ,dimension(3)    :: cell_size

        integer    :: dimensions
        integer    ,dimension(3,2)    :: cons_inner_loop

        character(len=20)        :: boundary_type_name

        integer    :: sign, bound_number
        integer :: i,j,k,dim,dim1,dim2,spec,plus

        dimensions        = this%domain%get_domain_dimensions()

        cons_inner_loop    = this%domain%get_local_inner_cells_bounds()

        cell_size        = this%mesh%mesh_ptr%get_cell_edges_length()

        associate (    rho_old            => this%rho_old%s_ptr        , &
                    rho_int            => this%rho_int%s_ptr        , &
                    p_dyn            => this%p_dyn%s_ptr            , &
                    F_a                => this%F_a%s_ptr            , &
                    F_b                => this%F_b%s_ptr            , &
                    H                => this%H%s_ptr                , &
                    v                => this%v%v_ptr                , &
                    v_f                => this%v_f%v_ptr            , &
                    v_f_old            => this%v_f_old%v_ptr        , &
                    bc                => this%boundary%bc_ptr)

            !$omp parallel default(shared) &
            !$omp& private(i,j,k,dim,dim2,vel_abs,plus,sign,bound_number,boundary_type_name,farfield_density, &
            !$omp& farfield_pressure)

            !$omp do collapse(3) schedule(static)
            do k = cons_inner_loop(3,1),cons_inner_loop(3,2)
            do j = cons_inner_loop(2,1),cons_inner_loop(2,2)
            do i = cons_inner_loop(1,1),cons_inner_loop(1,2)
                if(bc%bc_markers(i,j,k) == 0) then

                    vel_abs                = 0.0_dp
                    do dim = 1,dimensions
                        vel_abs                = vel_abs + (0.5_dp*(v_f%pr(dim)%cells(dim,i,j,k) + &
                            v_f%pr(dim)%cells(dim,i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3))))**2
                    end do

                    if (predictor) then
                        p_dyn%cells(i,j,k)    =  (H%cells(i,j,k) - 0.5_dp * vel_abs)*rho_old%cells(i,j,k)
                    else
                        p_dyn%cells(i,j,k)    =  (H%cells(i,j,k) - 0.5_dp * vel_abs)*rho_int%cells(i,j,k)
                    end if

                    do dim = 1,dimensions
                        do plus = 1,2
                            sign            = (-1)**plus
                            bound_number    = bc%bc_markers(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3))
                            if( bound_number /= 0 ) then
                                boundary_type_name = bc%boundary_types(bound_number)%get_type_name()
                                select case(boundary_type_name)
                                    case('wall')
                                        p_dyn%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3))    = &
                                            H%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3))

                                        vel_abs = 0.0_dp
                                        do dim2 = 1,dimensions
                                            vel_abs = vel_abs + v%pr(dim2)%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2), &
                                                k+sign*I_m(dim,3))**2
                                        end do
                                        p_dyn%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3))    = &
                                            p_dyn%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3)) - 0.5_dp*vel_abs
                                        if (predictor) then
                                            p_dyn%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3))    = &
                                                p_dyn%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3))* &
                                                rho_old%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3))
                                        else
                                            p_dyn%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3))    = &
                                                p_dyn%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3))* &
                                                rho_int%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3))
                                        end if
                                    case('outlet')
                                        p_dyn%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3)) = -p_dyn%cells(i,j,k)

                                    case('inlet')
                                        p_dyn%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3))    = &
                                            H%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3))

                                        vel_abs = 0.0_dp
                                        do dim2 = 1,dimensions
                                            vel_abs = vel_abs + v%pr(dim2)%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2), &
                                                k+sign*I_m(dim,3))**2
                                        end do
                                        p_dyn%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3))    = &
                                            p_dyn%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3)) - 0.5_dp*vel_abs
                                        if (predictor) then
                                            p_dyn%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3))    = &
                                                p_dyn%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3))* &
                                                rho_old%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3))
                                        else
                                            p_dyn%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3))    = &
                                                p_dyn%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3))* &
                                                rho_int%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3))
                                        end if
                                end select
                            end if
                        end do
                    end do
                end if
            end do
            end do
            end do

            !$omp end do
            !$omp end parallel

        end associate

    end subroutine
    !--------------------------------------------------------------------------
    ! Apply the pressure-potential and explicit-acceleration update to face-normal
    ! velocity, then reconstruct cell velocity and its boundary ghost values.
    !--------------------------------------------------------------------------
    subroutine calculate_velocity(this,time_step,predictor)
        class(fds_solver)    ,intent(inout)    :: this
        real(dp)            ,intent(in)        :: time_step
        logical                ,intent(in)        :: predictor

        real(dp)    ,dimension(3)    :: cell_size

        integer    :: dimensions, iterations
        integer    ,dimension(3,2)    :: cons_inner_loop, cons_utter_loop, flow_inner_loop
        integer    ,dimension(3,2)    :: loop

        integer                    :: plus, sign, bound_number
        character(len=20)        :: boundary_type_name

        real(dp)                :: farfield_velocity

        integer :: i,j,k,dim,dim1,dim2,spec

        dimensions        = this%domain%get_domain_dimensions()

        cons_inner_loop    = this%domain%get_local_inner_cells_bounds()
        cons_utter_loop = this%domain%get_local_utter_cells_bounds()

        flow_inner_loop    = this%domain%get_local_inner_faces_bounds()

        cell_size        = this%mesh%mesh_ptr%get_cell_edges_length()

        associate (    v                => this%v%v_ptr                , &
                    v_f                => this%v_f%v_ptr            , &
                    v_f_old            => this%v_f_old%v_ptr        , &
                    F_a                => this%F_a%s_ptr            , &
                    F_b                => this%F_b%s_ptr            , &
                    H                => this%H%s_ptr                , &
                    bc                => this%boundary%bc_ptr)

        !$omp parallel default(shared) &
        !$omp& private(i,j,k,dim,dim1,plus,sign,bound_number,boundary_type_name,farfield_velocity,loop)

            do dim = 1, dimensions
                loop(3,1) = cons_inner_loop(3,1)
                loop(3,2) = cons_utter_loop(3,2)*I_m(dim,3) + cons_inner_loop(3,2)*(1 - I_m(dim,3))

                loop(2,1) = cons_inner_loop(2,1)
                loop(2,2) = cons_utter_loop(2,2)*I_m(dim,2) + cons_inner_loop(2,2)*(1 - I_m(dim,2))

                loop(1,1) = cons_inner_loop(1,1)
                loop(1,2) = cons_utter_loop(1,2)*I_m(dim,1) + cons_inner_loop(1,2)*(1 - I_m(dim,1))

                !$omp do collapse(3) schedule(static)
                do k = loop(3,1), loop(3,2)
                do j = loop(2,1), loop(2,2)
                do i = loop(1,1), loop(1,2)

                    if((bc%bc_markers(i,j,k) == 0).or.(bc%bc_markers(i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3)) == 0)) then

                        if (predictor)    then
                            v_f_old%pr(dim)%cells(dim,i,j,k) = v_f%pr(dim)%cells(dim,i,j,k)
                            v_f%pr(dim)%cells(dim,i,j,k) = v_f%pr(dim)%cells(dim,i,j,k) - time_step * (F_a%cells(dim,i,j,k) + &
                                F_b%cells(dim,i,j,k)  +  (H%cells(i,j,k) - &
                                H%cells(i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3)))/cell_size(dim))
                        else
                            v_f%pr(dim)%cells(dim,i,j,k) = 0.5_dp * (v_f%pr(dim)%cells(dim,i,j,k) + &
                                v_f_old%pr(dim)%cells(dim,i,j,k)) - (0.5_dp * time_step) * (F_a%cells(dim,i,j,k) + &
                                F_b%cells(dim,i,j,k)  +  (H%cells(i,j,k) - &
                                H%cells(i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3)))/cell_size(dim) )
                        end if

                    end if
                end do
                end do
                end do
                !$omp end do

            end do

            !$omp do collapse(3) schedule(static)
            do k = cons_inner_loop(3,1),cons_inner_loop(3,2)
            do j = cons_inner_loop(2,1),cons_inner_loop(2,2)
            do i = cons_inner_loop(1,1),cons_inner_loop(1,2)
                if(bc%bc_markers(i,j,k) == 0) then
                    do dim = 1, dimensions
                        v%pr(dim)%cells(i,j,k) = 0.5_dp * ( v_f%pr(dim)%cells(dim,i,j,k) + &
                            v_f%pr(dim)%cells(dim,i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3)))
                    end do
                end if
            end do
            end do
            end do
            !$omp end do

            !$omp do collapse(3) schedule(static)
            do k = cons_inner_loop(3,1),cons_inner_loop(3,2)
            do j = cons_inner_loop(2,1),cons_inner_loop(2,2)
            do i = cons_inner_loop(1,1),cons_inner_loop(1,2)
                if(bc%bc_markers(i,j,k) == 0) then
                    do dim = 1, dimensions
                        do plus = 1,2
                            sign            = (-1)**plus
                            bound_number    = bc%bc_markers(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3))
                            if( bound_number /= 0 ) then
                                boundary_type_name = bc%boundary_types(bound_number)%get_type_name()
                                select case(boundary_type_name)
                                    case('wall')
                                        do dim1 = 1, dimensions
                                            if(dim1 == dim) then
                                                v%pr(dim1)%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3)) &
                                                    = -v%pr(dim1)%cells(i,j,k)
                                            else
                                                v%pr(dim1)%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3)) &
                                                    =  v%pr(dim1)%cells(i,j,k)
                                            end if
                                        end do
                                    case('outlet')
                                        do dim1 = 1, dimensions
                                            if(dim1 == dim) then
                                                v%pr(dim1)%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3)) &
                                                    =  v%pr(dim1)%cells(i,j,k)
                                            else
                                                v%pr(dim1)%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3)) &
                                                    =  v%pr(dim1)%cells(i,j,k)
                                            end if
                                        end do
                                    case('inlet')
                                        farfield_velocity        =  this%inlet_velocity

                                        do dim1 = 1, dimensions
                                            if(dim1 == dim) then
                                                v%pr(dim1)%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3)) &
                                                    =  farfield_velocity
                                            else
                                                v%pr(dim1)%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3)) &
                                                    =  0.0_dp
                                            end if
                                        end do
                                end select
                            end if
                        end do
                    end do
                end if

            end do
            end do
            end do

            !$omp end do
        !$omp end parallel

        end associate
    end subroutine
    !--------------------------------------------------------------------------
    ! Trapezoidal corrector for conservative species densities using the
    ! projected face velocity and the predictor state stored in rho_old,Y_old.
    !--------------------------------------------------------------------------
    subroutine calculate_interm_Y_corrector(this,time_step)
        class(fds_solver)    ,intent(inout)    :: this
        real(dp)            ,intent(in)        :: time_step

                real(dp) :: spec_summ

        real(dp), dimension (3,3):: lame_coeffs
        real(dp)    ,dimension(3)    :: cell_size
        real(dp), allocatable        :: flux_left_vec(:), flux_right_vec(:), rhs_vec(:), rhoY_new(:)

        integer    :: dimensions, species_number, coord_id
        integer    ,dimension(3,2)    :: cons_inner_loop
        character(len=20)    :: coordinate_system

        integer    :: i,j,k,dim,spec,particles_phase_counter

        dimensions        = this%domain%get_domain_dimensions()
        species_number    = this%chem%chem_ptr%species_number

        cons_inner_loop    = this%domain%get_local_inner_cells_bounds()

        cell_size        = this%mesh%mesh_ptr%get_cell_edges_length()
        if (time_step <= 0.0_dp) error stop 'FDS species update: non-positive time step'

        coordinate_system    = this%domain%get_coordinate_system_name()
        select case(coordinate_system)
            case ('cartesian')
                coord_id = 0
            case ('cylindrical')
                coord_id = 1
            case ('spherical')
                coord_id = 2
            case default
                coord_id = 0
        end select

        associate (    rho                    => this%rho%s_ptr            , &
                    rho_int                => this%rho_int%s_ptr        , &
                    rho_old                => this%rho_old%s_ptr        , &
                    v_f                    => this%v_f%v_ptr            , &
                    Y                    => this%Y%v_ptr                , &
                    Y_int                => this%Y_int%v_ptr            , &
                    Y_old                => this%Y_old%v_ptr            , &
                    mesh                => this%mesh%mesh_ptr        , &
                    bc                    => this%boundary%bc_ptr        , &
                    thermo                => this%thermo%thermo_ptr)

            ! Optimized scalar corrector: face-vector EOS correction, integer
            ! coordinate selector, and boundedness/reconstruction in the same loop.
            !$omp parallel default(shared) &
            !$omp& private(i,j,k,dim,spec,spec_summ,lame_coeffs,flux_left_vec,flux_right_vec,rhs_vec,rhoY_new, &
            !$omp& particles_phase_counter)
            allocate(flux_left_vec(species_number), flux_right_vec(species_number), rhs_vec(species_number), &
                rhoY_new(species_number))

            !$omp do collapse(3) schedule(static)
            do k = cons_inner_loop(3,1),cons_inner_loop(3,2)
            do j = cons_inner_loop(2,1),cons_inner_loop(2,2)
            do i = cons_inner_loop(1,1),cons_inner_loop(1,2)
                if(bc%bc_markers(i,j,k) == 0) then
                    lame_coeffs        = 1.0_dp
                    if (coord_id == 1) then
                        ! cylindrical: x -> r, y -> z
                        lame_coeffs(1,1)    =  mesh%mesh(1,i,j,k) - 0.5_dp*cell_size(1)
                        lame_coeffs(1,2)    =  mesh%mesh(1,i,j,k)
                        lame_coeffs(1,3)    =  mesh%mesh(1,i,j,k) + 0.5_dp*cell_size(1)
                    else if (coord_id == 2) then
                        ! spherical: x -> r
                        lame_coeffs(1,1)    =  (mesh%mesh(1,i,j,k) - 0.5_dp*cell_size(1))**2
                        lame_coeffs(1,2)    =  (mesh%mesh(1,i,j,k))**2
                        lame_coeffs(1,3)    =  (mesh%mesh(1,i,j,k) + 0.5_dp*cell_size(1))**2
                    end if

                    rhs_vec = 0.0_dp
                    do dim = 1, dimensions
                        if ((i*I_m(dim,1) + j*I_m(dim,2)  + k*I_m(dim,3)) < cons_inner_loop(dim,2)) then
                            call &
                                eos_corrected_species_face_vector(rho_int,Y_int,thermo%molar_masses,species_number,dim,i,j,k, &
                                1, &
                                 v_f%pr(dim)%cells(dim,i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3)),flux_right_vec)
                        else
                            if (v_f%pr(dim)%cells(dim,i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3)) > 0.0_dp) then
                                do spec = 1,species_number
                                    flux_right_vec(spec) = rho_int%cells(i,j,k) * Y_int%pr(spec)%cells(i,j,k)
                                end do
                            else
                                do spec = 1,species_number
                                    flux_right_vec(spec) = rho_int%cells(i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3)) * &
                                        Y_int%pr(spec)%cells(i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3))
                                end do
                            end if
                        end if

                        if ((i*I_m(dim,1) + j*I_m(dim,2)  + k*I_m(dim,3)) > cons_inner_loop(dim,1)) then
                            call &
                                eos_corrected_species_face_vector(rho_int,Y_int,thermo%molar_masses,species_number,dim,i,j,k,- &
                                1, &
                                 v_f%pr(dim)%cells(dim,i,j,k),flux_left_vec)
                        else
                            if (v_f%pr(dim)%cells(dim,i,j,k) > 0.0_dp) then
                                do spec = 1,species_number
                                    flux_left_vec(spec) = rho_int%cells(i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3)) * &
                                        Y_int%pr(spec)%cells(i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3))
                                end do
                            else
                                do spec = 1,species_number
                                    flux_left_vec(spec) = rho_int%cells(i,j,k) * Y_int%pr(spec)%cells(i,j,k)
                                end do
                            end if
                        end if

                        do spec = 1,species_number
                            rhs_vec(spec) = rhs_vec(spec) - (     flux_right_vec(spec) * &
                                v_f%pr(dim)%cells(dim,i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3)) * lame_coeffs(dim,3) &
                                                        - flux_left_vec(spec) * v_f%pr(dim)%cells(dim,i,j,k)* &
                                                            lame_coeffs(dim,1)) /cell_size(dim)/lame_coeffs(dim,2)
                        end do
                    end do

                    do spec = 1, species_number
                        rhoY_new(spec) = 0.5_dp * ( rho_old%cells(i,j,k)*Y_old%pr(spec)%cells(i,j,k) + &
                            rho_int%cells(i,j,k)*Y_int%pr(spec)%cells(i,j,k)) &
                                         + rhs_vec(spec) * (0.5_dp * time_step)
                        if (this%diffusion_flag) rhoY_new(spec) = rhoY_new(spec) + 0.5_dp * &
                            this%Y_prod_diff%v_ptr%pr(spec)%cells(i,j,k) * time_step
                        if (this%reactive_flag)  rhoY_new(spec) = rhoY_new(spec) + 0.5_dp * &
                            this%Y_prod_chem%v_ptr%pr(spec)%cells(i,j,k) * time_step
                        if (this%additional_particles_phases_number /= 0) then
                            do particles_phase_counter = 1, this%additional_particles_phases_number
                                rhoY_new(spec) = rhoY_new(spec) + &
                                    0.5_dp*this%Y_prod_particles(particles_phase_counter)%v_ptr%pr(spec)%cells(i,j,k)*time_step
                            end do
                        end if
                    end do

                    ! Conservative boundedness correction and reconstruction.
                    spec_summ = 0.0_dp
                    do spec = 1,species_number
                        rhoY_new(spec) = max(rhoY_new(spec), 0.0_dp)
                        spec_summ = spec_summ + rhoY_new(spec)
                    end do
                    if (spec_summ <= tiny(1.0_dp)) then
                        error stop 'FDS species update: non-positive total density'
                    else
                        rho%cells(i,j,k) = spec_summ
                        do spec = 1,species_number
                            Y%pr(spec)%cells(i,j,k) = rhoY_new(spec) / rho%cells(i,j,k)
                        end do
                    end if
                end if
            end do
            end do
            end do
            !$omp end do

            deallocate(flux_left_vec, flux_right_vec, rhs_vec, rhoY_new)
            !$omp end parallel

            end associate
    end subroutine
    !--------------------------------------------------------------------------
    ! Optional flame anchoring and laminar-burning-velocity workflow.  The routine
    ! tracks a temperature-gradient front, adapts inlet velocity, can export a
    ! stabilized flamelet, and optionally measures open-loop front drift.
    !--------------------------------------------------------------------------
    subroutine stabilizing_inlet_1D(this,time,stabilized)
!        ! Generalized stabilizing inlet routine with two explicit workflows:
!        !   2) stabilize/anchor flame for long-term observation without flamelet or drift measurement.
!        ! Anchored-v3 revision: keeps the heat-release gate and spatial anchor of v2,
!        ! but adds a fast capture mode for cases where the initial inlet velocity
!        ! convects the flame toward the outlet before the fine controller can react.
!        ! The routine keeps the old public signature, but the
!        ! controlled observable is no longer a single leading point.  It is the
!        ! heat-release-weighted flame centroid projected on the inlet axis.
!        !
!        ! Main stabilizing changes relative to the first generalized version:
!        !   1. sample-and-hold control: after every inlet correction the velocity
!        !      history is cleared and the controller waits for inlet ramp/settling;
!        !   2. tracking and control intervals are separated;
!        !   3. the inlet velocity is changed through a target/applied ramp, not by
!        !      an instantaneous boundary jump;
!        !   4. proportional gain is small and adaptively damped after sign changes
!        !      or error growth;
!        !   5. when two measurements bracket V_f=0, the next target is chosen by a
!        !      safeguarded bisection step;
!        !   6. a damped secant step is used only when the measured response slope is
!        !      reliable; otherwise the code falls back to a small proportional step;
!        !   7. the controller now uses a domain-safe-window anchor.  Inside the
!        !      usable domain the control target is V_f only; positional feedback is
!        !      activated only when the flame approaches a guarded boundary.  The
!        !      initial flame coordinate is diagnostic and is not used as an exact
!        !      positional set point;
!        !   8. heat-release validity is checked explicitly.  H-radical and T-gradient
!        !      centroids are diagnostics/fallback locations only and are not allowed
!        !      to drive the feedback after the burning front is lost.
!        !
!        ! inlet normal along another direction, replace this by a boundary-derived
!        ! inlet-normal coordinate.
        class(fds_solver), intent(inout) :: this
        real(dp),          intent(in)    :: time
        logical,           intent(out)   :: stabilized

        integer, parameter :: max_hist_size = 120
        integer, parameter :: max_diag_hist_size = 200
        integer, parameter :: min_hist_for_control = 12
        integer, parameter :: min_hist_for_capture = 6
        integer, parameter :: stable_required_count = 50
        integer, parameter :: flamelet_required_count = 50
        integer, parameter :: post_flamelet_hold_required_count = 100

        real(dp), parameter :: time_delay_default         = 1.0e-03_dp
        real(dp), parameter :: time_track_default         = 1.0e-04_dp
        real(dp), parameter :: time_control_default       = 5.0e-04_dp
        real(dp), parameter :: time_control_capture       = 2.0e-04_dp
        real(dp), parameter :: response_settle_default    = 5.0e-04_dp
        real(dp), parameter :: response_settle_capture    = 2.0e-04_dp
        real(dp), parameter :: inlet_ramp_time_fine       = 2.0e-03_dp
        real(dp), parameter :: inlet_ramp_time_capture    = 1.0e-03_dp
        real(dp), parameter :: controller_gain_initial    = 2.5e-01_dp
        real(dp), parameter :: controller_gain_min        = 2.0e-02_dp
        real(dp), parameter :: controller_gain_max        = 5.0e-01_dp
        real(dp), parameter :: controller_gain_capture    = 5.0e-01_dp
        real(dp), parameter :: controller_max_fraction    = 5.0e-02_dp
        real(dp), parameter :: controller_max_fraction_capture = 1.0e-01_dp
        real(dp), parameter :: controller_max_fraction_emergency = 2.5e-01_dp
        real(dp), parameter :: min_abs_velocity_step      = 2.0e-05_dp
        real(dp), parameter :: min_abs_velocity_step_capture = 1.0e-03_dp
        real(dp), parameter :: emergency_min_fraction     = 1.0e-01_dp
        real(dp), parameter :: velocity_tolerance_on      = 2.0e-05_dp
        real(dp), parameter :: velocity_tolerance_off     = 5.0e-06_dp
        real(dp), parameter :: filter_alpha               = 1.0e-01_dp
        real(dp), parameter :: secant_relaxation          = 2.0e-01_dp
        real(dp), parameter :: feedback_sign_default      = -1.0_dp
        real(dp), parameter :: heat_release_cut_fraction  = 1.0e-08_dp
        real(dp), parameter :: heat_release_valid_relative = 1.0e-07_dp
        real(dp), parameter :: heat_release_valid_absolute = 1.0e-20_dp
        real(dp), parameter :: position_relaxation_time   = 5.0e-02_dp
        real(dp), parameter :: position_tolerance_cells   = 4.0_dp
        real(dp), parameter :: capture_position_tolerance_cells = 12.0_dp
        real(dp), parameter :: capture_velocity_threshold = 2.0e-02_dp
        real(dp), parameter :: outlet_guard_cells         = 20.0_dp
        real(dp), parameter :: outlet_guard_fraction      = 1.5e-01_dp
        real(dp), parameter :: min_secant_du              = 5.0e-05_dp
        real(dp), parameter :: max_response_slope_abs     = 1.0e+03_dp
        real(dp), parameter :: persistent_error_factor    = 5.0_dp
        real(dp), parameter :: ramp_settle_fraction       = 1.0e-03_dp
        real(dp), parameter :: tiny_weight                = tiny(1.0_dp)
        logical,  parameter :: use_secant_control         = .false.

!        ! Workflow selector.
!        ! workflow_flamelet_sl:
!        !   Stabilize the flame, write 1D flamelet/chemistry data, perturb the inlet,
!        !   measure the linear drift velocity, write the laminar-burning-velocity result,
!        ! workflow_anchor_observation:
!        !   Keep the flame anchored for long-term observation.  No flamelet output, no
!        !   open-loop drift measurement, and no final stop signal is generated.
        integer,  parameter :: workflow_flamelet_sl = 1
        integer,  parameter :: workflow_anchor_observation = 2
        integer,  parameter :: active_workflow = workflow_flamelet_sl
        logical,  parameter :: enable_flamelet_output = (active_workflow == workflow_flamelet_sl)
        logical,  parameter :: enable_drift_measurement = (active_workflow == workflow_flamelet_sl)
        logical,  parameter :: pause_after_sl_measurement = (active_workflow == workflow_flamelet_sl)
        integer,  parameter :: stage_anchor_control = 0
        integer,  parameter :: stage_flamelet_ready = 1
        integer,  parameter :: stage_measurement_ramp = 2
        integer,  parameter :: stage_drift_measurement = 3
        integer,  parameter :: stage_measurement_done = 4
        integer,  parameter :: stage_measurement_failed = 5
        integer,  parameter :: min_hist_for_sl_measurement = 40
        real(dp), parameter :: measurement_delta_fraction = 3.0e-02_dp
        real(dp), parameter :: measurement_delta_sign = 1.0_dp
        real(dp), parameter :: measurement_delta_min = 5.0e-03_dp
        real(dp), parameter :: measurement_delta_max_fraction = 3.0e-01_dp
        real(dp), parameter :: measurement_delta_growth = 2.0_dp
        real(dp), parameter :: measurement_max_duration = 2.0e-02_dp
        real(dp), parameter :: measurement_no_motion_displacement_cells = 0.5_dp
        integer,  parameter :: measurement_max_attempts = 5
        real(dp), parameter :: measurement_ramp_time = 1.0e-03_dp
        real(dp), parameter :: measurement_settle_time = 1.0e-03_dp
        real(dp), parameter :: measurement_min_displacement_cells = 6.0_dp
        real(dp), parameter :: measurement_min_duration = 5.0e-03_dp
        real(dp), parameter :: measurement_r2_min = 9.95e-01_dp
        real(dp), parameter :: measurement_split_slope_rel_tol = 2.5e-01_dp
        real(dp), parameter :: measurement_split_slope_abs_tol = 2.0e-04_dp
        real(dp), parameter :: measurement_residual_cells = 1.0_dp

        integer :: dimensions, species_number, boundary_types
        integer :: front_axis
        integer :: i, j, k, dim, bound_number, specie_number, specie_index
        integer :: H2_index, H_index
        integer :: cons_inner_loop(3,2)
        integer :: active_track_number

        real(dp) :: cell_size(3), cell_volume
        real(dp) :: current_flame_location(3)
        real(dp) :: heat_release_centroid(3)
        real(dp) :: H_centroid(3)
        real(dp) :: Tgrad_centroid(3)
        real(dp) :: current_front_coord
        real(dp) :: flame_velocity_lsq, flame_velocity_filtered
        real(dp) :: diag_flame_velocity_lsq, diag_flame_velocity_filtered
        real(dp) :: control_velocity, position_error, position_velocity, position_control_error
        real(dp) :: front_spread, heat_release_integral
        real(dp) :: heat_release_max, heat_release_valid_limit, H_max, Tgrad_max
        real(dp) :: measured_inlet_velocity, proposed_inlet_velocity
        real(dp) :: du_raw, du_limited, max_velocity_step, X_H2
        real(dp) :: time_delay, time_track, time_control, response_settle_time, inlet_ramp_time
        real(dp) :: time_control_effective, response_settle_effective
        real(dp) :: coord(3), weight, qdot_cut, qdot, hval, tgrad
        real(dp) :: sum_weight, sum_s, sum_s2
        real(dp) :: sum_H_weight, sum_Tgrad_weight
        real(dp) :: s_coord
        real(dp) :: ramp_elapsed, ramp_fraction, ramp_residual, ramp_tolerance
        real(dp) :: target_step_for_log
        real(dp) :: sl_displacement, linear_r2, linear_rms, split_slope_diff
        real(dp) :: measurement_elapsed, measurement_displacement, measurement_delta
        logical :: measurement_linear_ok
        real(dp) :: front_position_tolerance, capture_position_tolerance
        real(dp) :: domain_front_min, domain_front_max, domain_front_length
        real(dp) :: outlet_guard_distance, outlet_distance
        real(dp) :: front_safe_min, front_safe_max
        real(dp), allocatable, save :: time_hist(:), front_coord_hist(:)
        real(dp), allocatable, save :: diag_time_hist(:), diag_front_coord_hist(:)
        real(dp), allocatable, save :: farfield_concentrations(:), concs(:)
        character(len=10), allocatable, save :: farfield_species_names(:)
        character(len=5) :: axis_names(3)
        character(len=20) :: flame_data_file
        character(len=500) :: av_header
        character(len=200), save :: data_table_filename, chem_table_filename
        character(len=100) :: chemical_mechanism
        logical :: found_inlet_farfield, trace_success, flame_detected, control_performed
        logical :: measurement_enabled, ramp_settled, reset_history_after_log
        logical :: capture_mode, emergency_mode

        integer, save :: track_counter = 0
        integer, save :: correction_counter = 0
        integer, save :: stabilization_counter = 0
        integer, save :: hist_count = 0
        integer, save :: diag_hist_count = 0
        integer, save :: same_sign_error_counter = 0
        integer, save :: post_flamelet_hold_counter = 0
        integer, save :: measurement_attempt = 0
        integer, save :: flame_loc_unit = -1
        real(dp), save :: previous_correction_time = -huge(1.0_dp)
        real(dp), save :: filtered_velocity_save = 0.0_dp
        real(dp), save :: diag_filtered_velocity_save = 0.0_dp
        real(dp), save :: previous_front_coord = 0.0_dp
        real(dp), save :: adaptive_gain = 5.0e-02_dp
        real(dp), save :: inlet_velocity_target = 0.0_dp
        real(dp), save :: inlet_velocity_applied = 0.0_dp
        real(dp), save :: ramp_start_time = 0.0_dp
        real(dp), save :: ramp_start_velocity = 0.0_dp
        real(dp), save :: active_inlet_ramp_time = inlet_ramp_time_fine
        real(dp), save :: previous_control_velocity = 0.0_dp
        real(dp), save :: previous_control_inlet_velocity = 0.0_dp
        real(dp), save :: front_reference_coord = 0.0_dp
        real(dp), save :: heat_release_peak_save = 0.0_dp
        real(dp), save :: bracket_u_a = 0.0_dp, bracket_v_a = 0.0_dp
        real(dp), save :: bracket_u_b = 0.0_dp, bracket_v_b = 0.0_dp
        logical, save :: initialized = .false.
        logical, save :: inlet_velocity_initialized = .false.
        logical, save :: output_initialized = .false.
        logical, save :: final_output_written = .false.
        logical, save :: flamelet_output_written = .false.
        logical, save :: sl_output_written = .false.
        logical, save :: have_previous_control_point = .false.
        logical, save :: has_bracket = .false.
        logical, save :: front_reference_initialized = .false.
        integer, save :: control_stage = stage_anchor_control
        real(dp), save :: stabilized_inlet_velocity = 0.0_dp
        real(dp), save :: measurement_inlet_velocity = 0.0_dp
        real(dp), save :: measurement_start_time = 0.0_dp
        real(dp), save :: measurement_start_coord = 0.0_dp
        real(dp), save :: sl_displacement_save = 0.0_dp
        real(dp), save :: measurement_velocity_save = 0.0_dp
        real(dp), save :: measurement_r2_save = 0.0_dp
        real(dp), save :: measurement_rms_save = 0.0_dp
        real(dp), save :: measurement_split_slope_diff_save = 0.0_dp
        real(dp), save :: current_measurement_delta = 0.0_dp

        stabilized = .false.

        dimensions      = this%domain%get_domain_dimensions()
        axis_names      = this%domain%get_axis_names()
        species_number  = this%chem%chem_ptr%species_number
        boundary_types  = this%boundary%bc_ptr%get_boundary_types()
        cons_inner_loop = this%domain%get_local_inner_cells_bounds()
        cell_size       = this%mesh%mesh_ptr%get_cell_edges_length()

        front_axis = 1
        if (front_axis > dimensions) front_axis = 1

        cell_volume = 1.0_dp
        do dim = 1, dimensions
            cell_volume = cell_volume * cell_size(dim)
        end do
        front_position_tolerance = max(position_tolerance_cells * cell_size(front_axis), 1.0e-8_dp)
        capture_position_tolerance = max(capture_position_tolerance_cells * cell_size(front_axis), &
            front_position_tolerance)
        domain_front_min = (real(cons_inner_loop(front_axis,1),dp) - 0.5_dp) * cell_size(front_axis)
        domain_front_max = (real(cons_inner_loop(front_axis,2),dp) - 0.5_dp) * cell_size(front_axis)
        domain_front_length = max(domain_front_max - domain_front_min, cell_size(front_axis))
        outlet_guard_distance = max(outlet_guard_cells * cell_size(front_axis), &
            outlet_guard_fraction * domain_front_length)

        H2_index = this%chem%chem_ptr%get_chemical_specie_index('H2')
        H_index  = this%chem%chem_ptr%get_chemical_specie_index('H')

        if (.not. allocated(time_hist)) then
            allocate(time_hist(max_hist_size), front_coord_hist(max_hist_size))
            time_hist        = 0.0_dp
            front_coord_hist = 0.0_dp
        end if

        if (.not. allocated(diag_time_hist)) then
            allocate(diag_time_hist(max_diag_hist_size), diag_front_coord_hist(max_diag_hist_size))
            diag_time_hist        = 0.0_dp
            diag_front_coord_hist = 0.0_dp
        end if

        if (.not. initialized) then
            allocate(concs(species_number))
            concs = 0.0_dp

            found_inlet_farfield = .false.
            do bound_number = 1, boundary_types
                if (this%boundary%bc_ptr%boundary_types(bound_number)%get_type_name() == 'inlet') then
                    call this%boundary%bc_ptr%boundary_types(bound_number)%get_farfield_concentrations(farfield_concentrations)
                    call this%boundary%bc_ptr%boundary_types(bound_number)%get_farfield_species_names(farfield_species_names)
                    found_inlet_farfield = .true.
                    exit
                end if
            end do

            if (found_inlet_farfield .and. allocated(farfield_species_names)) then
                do specie_number = 1, size(farfield_species_names)
                    specie_index = this%chem%chem_ptr%get_chemical_specie_index(farfield_species_names(specie_number))
                    if (specie_index >= 1 .and. specie_index <= species_number) then
                        concs(specie_index) = farfield_concentrations(specie_number)
                    end if
                end do
            end if

            if (sum(concs) > 0.0_dp .and. H2_index >= 1 .and. H2_index <= species_number) then
                X_H2 = concs(H2_index) / sum(concs) * 100.0_dp
            else
                X_H2 = 0.0_dp
            end if

            chemical_mechanism = trim(this%chem%chem_ptr%get_chemical_mechanism())
            data_table_filename = 'H2-Air_flamelet_' // trim(chemical_mechanism) // '_' // &
                trim(str_r(X_H2)) // '_pcnt_' // trim(str_e(cell_size(1))) // '_dx.dat'
            chem_table_filename = 'H2-Air_chem_table_' // trim(chemical_mechanism) // '_' // &
                trim(str_r(X_H2)) // '_pcnt_' // trim(str_e(cell_size(1))) // '_dx.dat'

            previous_correction_time = time
            adaptive_gain = controller_gain_initial
            initialized = .true.
        end if

        if (.not. inlet_velocity_initialized) then
            inlet_velocity_target  = this%inlet_velocity
            inlet_velocity_applied = this%inlet_velocity
            ramp_start_velocity = inlet_velocity_applied
            ramp_start_time = time
            inlet_velocity_initialized = .true.
        end if

        time_delay           = time_delay_default
        time_track           = time_track_default
        time_control         = time_control_default
        response_settle_time = response_settle_default
        inlet_ramp_time      = active_inlet_ramp_time

!        ! Apply a finite-duration linear inlet ramp on every call, not only on
!        ! tracking calls.  This avoids an acoustic jump but also avoids the very
!        ! long tail of exponential relaxation.
        if (inlet_ramp_time > 0.0_dp) then
            ramp_elapsed = max(time - ramp_start_time, 0.0_dp)
            ramp_fraction = min(ramp_elapsed / inlet_ramp_time, 1.0_dp)
        else
            ramp_fraction = 1.0_dp
        end if
        inlet_velocity_applied = ramp_start_velocity + &
            ramp_fraction * (inlet_velocity_target - ramp_start_velocity)
        this%inlet_velocity = inlet_velocity_applied

        ramp_residual = abs(inlet_velocity_target - inlet_velocity_applied)
        ramp_tolerance = max(ramp_settle_fraction * max(abs(inlet_velocity_target), min_abs_velocity_step), &
            1.0e-12_dp)
        ramp_settled = (ramp_residual <= ramp_tolerance .or. ramp_fraction >= 1.0_dp)

        if ((time - time_delay) / time_track <= real(track_counter + 1, dp)) return

        associate (T => this%T%s_ptr, &
                   Y => this%Y%v_ptr, &
                   E_f_prod_chem => this%E_f_prod_chem%s_ptr, &
                   bc => this%boundary%bc_ptr)

            if (.not. output_initialized) call initialize_output_file()

!            ! --- First pass: maxima for robust thresholds and diagnostics.
            heat_release_max = 0.0_dp
            H_max            = 0.0_dp
            Tgrad_max        = 0.0_dp
            do k = cons_inner_loop(3,1), cons_inner_loop(3,2)
            do j = cons_inner_loop(2,1), cons_inner_loop(2,2)
            do i = cons_inner_loop(1,1), cons_inner_loop(1,2)
                if (bc%bc_markers(i,j,k) /= 0) cycle

                qdot = max(this%E_f_prod_chem%s_ptr%cells(i,j,k), 0.0_dp)
                heat_release_max = max(heat_release_max, qdot)

                if (H_index >= 1 .and. H_index <= species_number) then
                    hval = abs(Y%pr(H_index)%cells(i,j,k))
                    H_max = max(H_max, hval)
                end if

                tgrad = temperature_gradient_norm(i,j,k)
                Tgrad_max = max(Tgrad_max, tgrad)
            end do
            end do
            end do

            heat_release_peak_save = max(heat_release_peak_save, heat_release_max)
            heat_release_valid_limit = max(heat_release_valid_absolute, &
                heat_release_valid_relative * heat_release_peak_save)

!            ! --- Second pass: heat-release centroid and fallback centroids.
            heat_release_centroid = 0.0_dp
            H_centroid            = 0.0_dp
            Tgrad_centroid        = 0.0_dp
            sum_weight            = 0.0_dp
            sum_H_weight          = 0.0_dp
            sum_Tgrad_weight      = 0.0_dp
            sum_s                 = 0.0_dp
            sum_s2                = 0.0_dp
            heat_release_integral = 0.0_dp

            qdot_cut = heat_release_cut_fraction * heat_release_max

            do k = cons_inner_loop(3,1), cons_inner_loop(3,2)
            do j = cons_inner_loop(2,1), cons_inner_loop(2,2)
            do i = cons_inner_loop(1,1), cons_inner_loop(1,2)
                if (bc%bc_markers(i,j,k) /= 0) cycle

                coord = cell_center_coordinates(i,j,k)
                s_coord = coord(front_axis)

                qdot = max(this%E_f_prod_chem%s_ptr%cells(i,j,k), 0.0_dp)
                if (qdot > qdot_cut) then
                    weight = qdot * cell_volume
                    heat_release_centroid = heat_release_centroid + weight * coord
                    sum_weight = sum_weight + weight
                    sum_s  = sum_s  + weight * s_coord
                    sum_s2 = sum_s2 + weight * s_coord * s_coord
                    heat_release_integral = heat_release_integral + weight
                end if

                if (H_index >= 1 .and. H_index <= species_number) then
                    hval = max(Y%pr(H_index)%cells(i,j,k), 0.0_dp)
                    if (hval > 0.0_dp) then
                        weight = hval * cell_volume
                        H_centroid = H_centroid + weight * coord
                        sum_H_weight = sum_H_weight + weight
                    end if
                end if

                tgrad = temperature_gradient_norm(i,j,k)
                if (tgrad > 0.0_dp) then
                    weight = tgrad * cell_volume
                    Tgrad_centroid = Tgrad_centroid + weight * coord
                    sum_Tgrad_weight = sum_Tgrad_weight + weight
                end if
            end do
            end do
            end do

            trace_success = .false.
            flame_detected = .false.
            current_flame_location = 0.0_dp
            front_spread = 0.0_dp

            if (sum_weight > tiny_weight) then
                current_flame_location = heat_release_centroid / sum_weight
                front_spread = max(sum_s2 / sum_weight - (sum_s / sum_weight)**2, 0.0_dp)
                front_spread = sqrt(front_spread)
                trace_success = .true.
                flame_detected = (heat_release_max >= heat_release_valid_limit)
            else if (sum_H_weight > tiny_weight) then
!                ! Diagnostic fallback only.  Do not use this point for feedback control.
                current_flame_location = H_centroid / sum_H_weight
                trace_success = .true.
            else if (sum_Tgrad_weight > tiny_weight) then
!                ! Diagnostic fallback only.  Do not use this point for feedback control.
                current_flame_location = Tgrad_centroid / sum_Tgrad_weight
                trace_success = .true.
            end if

            if (.not. trace_success) then
                track_counter = track_counter + 1
                return
            end if

            current_front_coord = current_flame_location(front_axis)
            if (flame_detected .and. .not. front_reference_initialized) then
                front_reference_coord = current_front_coord
                front_reference_initialized = .true.
            end if

            front_safe_min = domain_front_min + outlet_guard_distance
            front_safe_max = domain_front_max - outlet_guard_distance
!            ! Do not intersect the safe window with front_reference_coord +/-
!            ! capture_position_tolerance.  In coarse-grid runs the flame can settle at a
!            ! nearby discrete equilibrium that is slightly shifted from the first detected
!            ! centroid.  Treating this harmless offset as a positional error produces
!            ! two-level inlet-velocity dithering.  The reference coordinate is therefore
!            ! kept for diagnostics, while positional feedback is used only as boundary
!            ! protection.
            if (front_safe_min >= front_safe_max) then
                front_safe_min = domain_front_min
                front_safe_max = domain_front_max
            end if

            diag_flame_velocity_lsq = 0.0_dp
            diag_flame_velocity_filtered = 0.0_dp
            if (flame_detected) then
                call append_diagnostic_history(time, current_front_coord)
                diag_flame_velocity_lsq = diagnostic_least_squares_velocity()
                if (diag_hist_count <= 2) then
                    diag_flame_velocity_filtered = diag_flame_velocity_lsq
                    diag_filtered_velocity_save = diag_flame_velocity_filtered
                else
                    diag_flame_velocity_filtered = (1.0_dp - filter_alpha) * diag_filtered_velocity_save + &
                        filter_alpha * diag_flame_velocity_lsq
                    diag_filtered_velocity_save = diag_flame_velocity_filtered
                end if
            else
                diag_hist_count = 0
                diag_filtered_velocity_save = 0.0_dp
            end if

            control_performed = .false.
            reset_history_after_log = .false.
            du_raw = 0.0_dp
            du_limited = 0.0_dp
            target_step_for_log = 0.0_dp
            flame_velocity_lsq = diag_flame_velocity_lsq
            flame_velocity_filtered = diag_flame_velocity_filtered
            position_error = 0.0_dp
            position_control_error = 0.0_dp
            position_velocity = 0.0_dp
            control_velocity = 0.0_dp
            sl_displacement = sl_displacement_save
            linear_r2 = measurement_r2_save
            linear_rms = measurement_rms_save
            split_slope_diff = measurement_split_slope_diff_save
            measurement_elapsed = 0.0_dp
            measurement_displacement = 0.0_dp
            measurement_linear_ok = .false.
            capture_mode = .false.
            emergency_mode = .false.

            if (front_reference_initialized) then
                position_error = current_front_coord - front_reference_coord
                position_control_error = safe_window_error(current_front_coord, front_safe_min, front_safe_max)
                position_velocity = position_control_error / max(position_relaxation_time, time_track)
                outlet_distance = domain_front_max - current_front_coord
                capture_mode = (position_control_error /= 0.0_dp) .or. &
                    (outlet_distance < outlet_guard_distance)
                emergency_mode = (outlet_distance < 0.5_dp * outlet_guard_distance) .or. &
                    (abs(position_control_error) > capture_position_tolerance)
            end if

            if (capture_mode) then
                response_settle_effective = response_settle_capture
                time_control_effective = time_control_capture
            else
                response_settle_effective = response_settle_time
                time_control_effective = time_control
            end if

!            ! Do not collect velocity statistics while the inlet is still ramping,
!            ! while the flow is inside the post-correction settling interval, or when
!            ! the actual heat-release front is lost.  H/T-gradient fallback points are
!            ! logged but never allowed to drive the inlet.
            measurement_enabled = flame_detected .and. front_reference_initialized .and. ramp_settled .and. &
                ((control_stage == stage_anchor_control) .or. &
                 (control_stage == stage_flamelet_ready)) .and. &
                ((time - previous_correction_time) >= response_settle_effective)

            if (measurement_enabled) then
                call append_front_history(time, current_front_coord)
                flame_velocity_lsq = least_squares_velocity()

                if (hist_count <= 2) then
                    flame_velocity_filtered = flame_velocity_lsq
                    filtered_velocity_save  = flame_velocity_filtered
                else
                    flame_velocity_filtered = (1.0_dp - filter_alpha) * filtered_velocity_save + &
                        filter_alpha * flame_velocity_lsq
                    filtered_velocity_save = flame_velocity_filtered
                end if

                position_error = current_front_coord - front_reference_coord
                position_control_error = safe_window_error(current_front_coord, front_safe_min, front_safe_max)
                position_velocity = position_control_error / max(position_relaxation_time, time_track)
                control_velocity = flame_velocity_filtered + position_velocity
                outlet_distance = domain_front_max - current_front_coord
                capture_mode = (position_control_error /= 0.0_dp) .or. &
                    (abs(flame_velocity_filtered) > capture_velocity_threshold) .or. &
                    (outlet_distance < outlet_guard_distance)
                emergency_mode = (outlet_distance < 0.5_dp * outlet_guard_distance) .or. &
                    (abs(position_control_error) > capture_position_tolerance)
                if (capture_mode) then
                    time_control_effective = time_control_capture
                else
                    time_control_effective = time_control
                end if
            else if (control_stage == stage_anchor_control) then
                call clear_control_history()
                if (.not. flame_detected) has_bracket = .false.
            else if (.not. flame_detected) then
                call clear_control_history()
                has_bracket = .false.
            end if

!            ! --- Open-loop drift measurement after the flame has first been stabilized.
!            !     In this stage the inlet target is fixed.  The flame-front coordinate
!            !     must be well approximated by a straight line in time.  The displacement
!            !     flame-speed estimate is then S_L = U_in - V_f for the current sign
!            !     convention.  For other coordinate conventions only this diagnostic sign
!            !     may need to be changed.
            if (control_stage == stage_measurement_ramp) then
                if (flame_detected .and. ramp_settled .and. &
                    (time - ramp_start_time) >= measurement_settle_time) then
                    call clear_control_history()
                    measurement_start_time = time
                    measurement_start_coord = current_front_coord
                    control_stage = stage_drift_measurement
                end if
            else if (control_stage == stage_drift_measurement) then
                if (flame_detected .and. ramp_settled) then
                    call append_front_history(time, current_front_coord)
                    flame_velocity_lsq = least_squares_velocity()
                    flame_velocity_filtered = flame_velocity_lsq
                    measurement_velocity_save = flame_velocity_lsq
                    control_velocity = flame_velocity_lsq
                    measurement_elapsed = time - measurement_start_time
                    measurement_displacement = current_front_coord - measurement_start_coord
                    call drift_linearity_diagnostics(linear_r2, linear_rms, split_slope_diff)
                    measurement_r2_save = linear_r2
                    measurement_rms_save = linear_rms
                    measurement_split_slope_diff_save = split_slope_diff
                    sl_displacement = inlet_velocity_target - flame_velocity_lsq
                    sl_displacement_save = sl_displacement
                    measurement_linear_ok = (hist_count >= min_hist_for_sl_measurement) .and. &
                        (measurement_elapsed >= measurement_min_duration) .and. &
                        (abs(measurement_displacement) >= measurement_min_displacement_cells * cell_size(front_axis)) .and. &
                        (linear_r2 >= measurement_r2_min) .and. &
                        (linear_rms <= measurement_residual_cells * cell_size(front_axis)) .and. &
                        (split_slope_diff <= max(measurement_split_slope_abs_tol, &
                            measurement_split_slope_rel_tol * max(abs(flame_velocity_lsq), velocity_tolerance_on)))
                    if (measurement_linear_ok) then
                        measurement_velocity_save = flame_velocity_lsq
                        control_stage = stage_measurement_done
                        stabilization_counter = stable_required_count
                        call write_laminar_velocity_once()
                    else if (measurement_elapsed >= measurement_max_duration .and. &
                        abs(measurement_displacement) < measurement_no_motion_displacement_cells * cell_size(front_axis)) then
!                        ! The perturbation did not move the heat-release centroid by even a
!                        ! fraction of a cell.  This usually means numerical pinning or that the
!                        ! imposed velocity perturbation is too small compared with the discrete
!                        ! front-location uncertainty.  Retry with a larger open-loop perturbation.
                        if (measurement_attempt < measurement_max_attempts) then
                            measurement_attempt = measurement_attempt + 1
                            current_measurement_delta = min(measurement_delta_growth * max(current_measurement_delta, &
                                measurement_delta_min), measurement_delta_max_fraction * &
                                max(abs(stabilized_inlet_velocity), min_abs_velocity_step))
                            measurement_inlet_velocity = max(stabilized_inlet_velocity + &
                                measurement_delta_sign * current_measurement_delta, 0.0_dp)
                            inlet_velocity_target = measurement_inlet_velocity
                            ramp_start_velocity = inlet_velocity_applied
                            ramp_start_time = time
                            active_inlet_ramp_time = measurement_ramp_time
                            control_stage = stage_measurement_ramp
                            sl_displacement_save = 0.0_dp
                            measurement_velocity_save = 0.0_dp
                            measurement_r2_save = 0.0_dp
                            measurement_rms_save = 0.0_dp
                            measurement_split_slope_diff_save = 0.0_dp
                            call clear_control_history()
                        else
                            control_stage = stage_measurement_failed
                        end if
                    end if
                else
                    call clear_control_history()
                end if
            end if

!            ! --- Feedback control.  The correction is rare and based on the response
!            !     measured after the previous target velocity has settled.
            if (measurement_enabled .and. &
                hist_count >= merge(min_hist_for_capture, min_hist_for_control, capture_mode)) then
                if ((time - previous_correction_time) >= time_control_effective) then
                    if (control_action_needed()) then
                        measured_inlet_velocity = inlet_velocity_target

                        call update_adaptive_gain(control_velocity)
                        call update_bracket(measured_inlet_velocity, control_velocity)
                        call choose_new_inlet_target(measured_inlet_velocity, control_velocity, &
                            proposed_inlet_velocity, du_raw, du_limited)

                        if (abs(proposed_inlet_velocity - measured_inlet_velocity) > 0.0_dp) then
                            inlet_velocity_target = proposed_inlet_velocity
                            target_step_for_log = inlet_velocity_target - measured_inlet_velocity
                            if (capture_mode) then
                                active_inlet_ramp_time = inlet_ramp_time_capture
                            else
                                active_inlet_ramp_time = inlet_ramp_time_fine
                            end if
                            ramp_start_velocity = inlet_velocity_applied
                            ramp_start_time = time
                            previous_correction_time = time
                            correction_counter = correction_counter + 1
                            control_performed = .true.
                            reset_history_after_log = .true.
                        end if

                        previous_control_inlet_velocity = measured_inlet_velocity
                        previous_control_velocity = control_velocity
                        have_previous_control_point = .true.
                    end if
                end if
            end if

            if (control_stage == stage_anchor_control) then
                if (measurement_enabled .and. ramp_settled .and. hist_count >= min_hist_for_control .and. &
                    abs(diag_flame_velocity_filtered) < velocity_tolerance_off .and. &
                    position_control_error == 0.0_dp) then
                    stabilization_counter = stabilization_counter + 1
                else if (control_performed .or. .not. flame_detected .or. &
                    abs(diag_flame_velocity_filtered) > velocity_tolerance_on .or. &
                    position_control_error /= 0.0_dp) then
                    stabilization_counter = 0
                end if
            else if (control_stage == stage_flamelet_ready) then
                if (measurement_enabled .and. ramp_settled .and. hist_count >= min_hist_for_control .and. &
                    abs(diag_flame_velocity_filtered) < velocity_tolerance_off .and. &
                    position_control_error == 0.0_dp) then
                    post_flamelet_hold_counter = post_flamelet_hold_counter + 1
                    stabilization_counter = post_flamelet_hold_counter
                else if (control_performed .or. .not. flame_detected .or. &
                    abs(diag_flame_velocity_filtered) > velocity_tolerance_on .or. &
                    position_control_error /= 0.0_dp) then
                    post_flamelet_hold_counter = 0
                    stabilization_counter = 0
                end if
            else if (control_stage == stage_measurement_done) then
                stabilization_counter = stable_required_count
            end if

            if (control_stage == stage_anchor_control .or. control_stage == stage_flamelet_ready) then
                flame_velocity_lsq = diag_flame_velocity_lsq
                flame_velocity_filtered = diag_flame_velocity_filtered
            else if (control_stage == stage_measurement_done) then
                flame_velocity_lsq = measurement_velocity_save
                flame_velocity_filtered = measurement_velocity_save
                control_velocity = measurement_velocity_save
                sl_displacement = sl_displacement_save
                linear_r2 = measurement_r2_save
                linear_rms = measurement_rms_save
                split_slope_diff = measurement_split_slope_diff_save
            end if

            active_track_number = track_counter
            call write_tracking_line()

            if (reset_history_after_log) call clear_control_history()

            previous_front_coord = current_front_coord
            track_counter = track_counter + 1

            if (control_stage == stage_anchor_control .and. &
                stabilization_counter >= flamelet_required_count) then
                if (active_workflow == workflow_flamelet_sl) then
!                    ! Scenario 1: the first stabilization stage is complete.  Write the
!                    ! flamelet before any open-loop velocity perturbation is applied.
                    if (enable_flamelet_output) call write_flamelet_tables_once()
                    control_stage = stage_flamelet_ready
                    stabilization_counter = 0
                    post_flamelet_hold_counter = 0
                    call clear_control_history()
                else
!                    ! Scenario 2: keep anchoring indefinitely.  Do not write flamelet
!                    ! data, do not perturb the inlet, and do not return a stop signal.
                    stabilized = .false.
                end if

            else if (control_stage == stage_flamelet_ready .and. &
                post_flamelet_hold_counter >= post_flamelet_hold_required_count) then
                if (enable_drift_measurement .and. dimensions == 1) then
                    call start_drift_measurement_ramp(time)
                else if (pause_after_sl_measurement) then
                    stabilized = .true.
                end if

            else if (control_stage == stage_measurement_done) then
                if (pause_after_sl_measurement) stabilized = .true.
            end if

        end associate

    contains

        subroutine initialize_output_file()
            integer :: local_dim

            write(flame_data_file,'(A,I0,A)') 'av_flame_data_', this%load_counter, '.dat'
            open(newunit = flame_loc_unit, file = flame_data_file, status = 'replace', form = 'formatted')

            av_header = 'VARIABLES="time" '
            do local_dim = 1, dimensions
                av_header = trim(av_header) // '"xf_' // trim(axis_names(local_dim)) // '" '
            end do
            av_header = trim(av_header) // &
                '"Vfl_lsq" "Vfl_filtered" "Vfl_diag_lsq" "Vfl_diag_filtered" "Vfl_measurement_lsq" "Vcontrol" "pos_error" "x_ref" '
            av_header = trim(av_header) // '"U_in_applied" "U_in_target" '
            av_header = trim(av_header) // '"dU_target" "adaptive_gain" "hist_count" ' // &
                '"measurement_on" "bracket_on" "capture_on" "emergency_on" "flame_detected" '
            av_header = trim(av_header) // '"front_spread" "Qint" "Qmax" "Qvalid" "Hmax" "Tgradmax" '
            av_header = trim(av_header) // '"stage" "SL_disp" "lin_R2" "lin_RMS" "split_dV" '
            av_header = trim(av_header) // '"corr_count" "stab_count" "track_count"'

            write(flame_loc_unit,'(A)') trim(av_header)
            output_initialized = .true.
        end subroutine initialize_output_file

        function cell_center_coordinates(ii,jj,kk) result(xc)
            integer, intent(in) :: ii, jj, kk
            real(dp) :: xc(3)

            xc = 0.0_dp
            xc(1) = (real(ii,dp) - 0.5_dp) * cell_size(1)
            if (dimensions >= 2) xc(2) = (real(jj,dp) - 0.5_dp) * cell_size(2)
            if (dimensions >= 3) xc(3) = (real(kk,dp) - 0.5_dp) * cell_size(3)
        end function cell_center_coordinates

        function temperature_gradient_norm(ii,jj,kk) result(grad_norm)
            integer, intent(in) :: ii, jj, kk
            real(dp) :: grad_norm
            real(dp) :: g2, gd

            g2 = 0.0_dp
            if (dimensions >= 1) then
                if (ii > cons_inner_loop(1,1) .and. ii < cons_inner_loop(1,2)) then
                    gd = (this%T%s_ptr%cells(ii+1,jj,kk) - this%T%s_ptr%cells(ii-1,jj,kk)) / (2.0_dp * cell_size(1))
                    g2 = g2 + gd * gd
                end if
            end if
            if (dimensions >= 2) then
                if (jj > cons_inner_loop(2,1) .and. jj < cons_inner_loop(2,2)) then
                    gd = (this%T%s_ptr%cells(ii,jj+1,kk) - this%T%s_ptr%cells(ii,jj-1,kk)) / (2.0_dp * cell_size(2))
                    g2 = g2 + gd * gd
                end if
            end if
            if (dimensions >= 3) then
                if (kk > cons_inner_loop(3,1) .and. kk < cons_inner_loop(3,2)) then
                    gd = (this%T%s_ptr%cells(ii,jj,kk+1) - this%T%s_ptr%cells(ii,jj,kk-1)) / (2.0_dp * cell_size(3))
                    g2 = g2 + gd * gd
                end if
            end if
            grad_norm = sqrt(g2)
        end function temperature_gradient_norm

        subroutine append_front_history(t_new, s_new)
            real(dp), intent(in) :: t_new, s_new

            if (hist_count < max_hist_size) then
                hist_count = hist_count + 1
                time_hist(hist_count) = t_new
                front_coord_hist(hist_count) = s_new
            else
                time_hist(1:max_hist_size-1) = time_hist(2:max_hist_size)
                front_coord_hist(1:max_hist_size-1) = front_coord_hist(2:max_hist_size)
                time_hist(max_hist_size) = t_new
                front_coord_hist(max_hist_size) = s_new
            end if
        end subroutine append_front_history

        subroutine append_diagnostic_history(t_new, s_new)
            real(dp), intent(in) :: t_new, s_new

            if (diag_hist_count < max_diag_hist_size) then
                diag_hist_count = diag_hist_count + 1
                diag_time_hist(diag_hist_count) = t_new
                diag_front_coord_hist(diag_hist_count) = s_new
            else
                diag_time_hist(1:max_diag_hist_size-1) = diag_time_hist(2:max_diag_hist_size)
                diag_front_coord_hist(1:max_diag_hist_size-1) = diag_front_coord_hist(2:max_diag_hist_size)
                diag_time_hist(max_diag_hist_size) = t_new
                diag_front_coord_hist(max_diag_hist_size) = s_new
            end if
        end subroutine append_diagnostic_history

        subroutine clear_control_history()
            hist_count = 0
            filtered_velocity_save = 0.0_dp
            time_hist = 0.0_dp
            front_coord_hist = 0.0_dp
        end subroutine clear_control_history

        function least_squares_velocity() result(vfit)
            real(dp) :: vfit
            integer :: n
            real(dp) :: t_av, s_av, numerator, denominator

            if (hist_count < 2) then
                vfit = 0.0_dp
                return
            end if

            t_av = sum(time_hist(1:hist_count)) / real(hist_count, dp)
            s_av = sum(front_coord_hist(1:hist_count)) / real(hist_count, dp)
            numerator   = 0.0_dp
            denominator = 0.0_dp
            do n = 1, hist_count
                numerator = numerator + (time_hist(n) - t_av) * (front_coord_hist(n) - s_av)
                denominator = denominator + (time_hist(n) - t_av)**2
            end do

            if (denominator > tiny(denominator)) then
                vfit = numerator / denominator
            else
                vfit = 0.0_dp
            end if
        end function least_squares_velocity

        function diagnostic_least_squares_velocity() result(vfit)
            real(dp) :: vfit
            integer :: n
            real(dp) :: t_av, s_av, numerator, denominator

            if (diag_hist_count < 2) then
                vfit = 0.0_dp
                return
            end if

            t_av = sum(diag_time_hist(1:diag_hist_count)) / real(diag_hist_count, dp)
            s_av = sum(diag_front_coord_hist(1:diag_hist_count)) / real(diag_hist_count, dp)
            numerator   = 0.0_dp
            denominator = 0.0_dp
            do n = 1, diag_hist_count
                numerator = numerator + (diag_time_hist(n) - t_av) * &
                    (diag_front_coord_hist(n) - s_av)
                denominator = denominator + (diag_time_hist(n) - t_av)**2
            end do

            if (denominator > tiny(denominator)) then
                vfit = numerator / denominator
            else
                vfit = 0.0_dp
            end if
        end function diagnostic_least_squares_velocity

        subroutine drift_linearity_diagnostics(r2_out, rms_out, split_slope_diff_out)
            real(dp), intent(out) :: r2_out, rms_out, split_slope_diff_out
            integer :: n, mid
            real(dp) :: v_all, t_av, s_av, intercept, ss_tot, ss_res, residual
            real(dp) :: v_first, v_second

            r2_out = 0.0_dp
            rms_out = huge(1.0_dp)
            split_slope_diff_out = huge(1.0_dp)
            if (hist_count < 4) return

            v_all = least_squares_velocity()
            t_av = sum(time_hist(1:hist_count)) / real(hist_count, dp)
            s_av = sum(front_coord_hist(1:hist_count)) / real(hist_count, dp)
            intercept = s_av - v_all * t_av
            ss_tot = 0.0_dp
            ss_res = 0.0_dp
            do n = 1, hist_count
                residual = front_coord_hist(n) - (intercept + v_all * time_hist(n))
                ss_res = ss_res + residual * residual
                ss_tot = ss_tot + (front_coord_hist(n) - s_av)**2
            end do
            rms_out = sqrt(ss_res / real(hist_count, dp))
            if (ss_tot > tiny(ss_tot)) then
                r2_out = max(0.0_dp, 1.0_dp - ss_res / ss_tot)
            else
                r2_out = 0.0_dp
            end if

            mid = hist_count / 2
            v_first = least_squares_velocity_range(1, mid)
            v_second = least_squares_velocity_range(mid + 1, hist_count)
            split_slope_diff_out = abs(v_second - v_first)
        end subroutine drift_linearity_diagnostics

        function least_squares_velocity_range(n_first, n_last) result(vfit)
            integer, intent(in) :: n_first, n_last
            real(dp) :: vfit
            integer :: n, n_local
            real(dp) :: t_av, s_av, numerator, denominator

            n_local = n_last - n_first + 1
            if (n_local < 2) then
                vfit = 0.0_dp
                return
            end if

            t_av = 0.0_dp
            s_av = 0.0_dp
            do n = n_first, n_last
                t_av = t_av + time_hist(n)
                s_av = s_av + front_coord_hist(n)
            end do
            t_av = t_av / real(n_local, dp)
            s_av = s_av / real(n_local, dp)

            numerator = 0.0_dp
            denominator = 0.0_dp
            do n = n_first, n_last
                numerator = numerator + (time_hist(n) - t_av) * (front_coord_hist(n) - s_av)
                denominator = denominator + (time_hist(n) - t_av)**2
            end do
            if (denominator > tiny(denominator)) then
                vfit = numerator / denominator
            else
                vfit = 0.0_dp
            end if
        end function least_squares_velocity_range

        pure function deadband_error(value, tolerance) result(error_out)
            real(dp), intent(in) :: value, tolerance
            real(dp) :: error_out

            if (abs(value) <= tolerance) then
                error_out = 0.0_dp
            else
                error_out = sign(abs(value) - tolerance, value)
            end if
        end function deadband_error

        pure function safe_window_error(value, lower_bound, upper_bound) result(error_out)
            real(dp), intent(in) :: value, lower_bound, upper_bound
            real(dp) :: error_out

            if (value < lower_bound) then
                error_out = value - lower_bound
            else if (value > upper_bound) then
                error_out = value - upper_bound
            else
                error_out = 0.0_dp
            end if
        end function safe_window_error

        logical function control_action_needed()
            control_action_needed = (abs(flame_velocity_filtered) > velocity_tolerance_on) .or. &
                (abs(position_control_error) > 0.0_dp .and. abs(control_velocity) > velocity_tolerance_on)
        end function control_action_needed

        subroutine update_adaptive_gain(v_current)
            real(dp), intent(in) :: v_current

!            The previous version damped the gain whenever the measured velocity
!            magnitude increased.  In the supplied logs this drove the gain to its
!            minimum while the flame was still drifting monotonically.  Here the gain
!            is reduced only after a true sign reversal of the anchored control error.
!            Same-sign errors are treated as persistent drift and the gain is allowed
!            to recover.  The finite step limiter and inlet ramp remain responsible for
!            preventing acoustic forcing.
            if (have_previous_control_point) then
                if (v_current * previous_control_velocity < 0.0_dp) then
                    same_sign_error_counter = 0
                    adaptive_gain = max(0.7_dp * adaptive_gain, controller_gain_min)
                else
                    same_sign_error_counter = same_sign_error_counter + 1
                    if (same_sign_error_counter >= 2) then
                        adaptive_gain = min(1.10_dp * adaptive_gain, controller_gain_max)
                    else
                        adaptive_gain = min(1.03_dp * adaptive_gain, controller_gain_max)
                    end if
                end if
            else
                same_sign_error_counter = 0
            end if
        end subroutine update_adaptive_gain

        subroutine update_bracket(u_current, v_current)
            real(dp), intent(in) :: u_current, v_current

            if (abs(v_current) <= velocity_tolerance_off) return

            if (have_previous_control_point) then
                if (abs(u_current - previous_control_inlet_velocity) >= min_secant_du .and. &
                    v_current * previous_control_velocity < 0.0_dp) then
                    bracket_u_a = previous_control_inlet_velocity
                    bracket_v_a = previous_control_velocity
                    bracket_u_b = u_current
                    bracket_v_b = v_current
                    has_bracket = .true.
                end if
            end if

            if (has_bracket) then
                if (v_current * bracket_v_a > 0.0_dp) then
                    bracket_u_a = u_current
                    bracket_v_a = v_current
                else if (v_current * bracket_v_b > 0.0_dp) then
                    bracket_u_b = u_current
                    bracket_v_b = v_current
                end if
                if (bracket_v_a * bracket_v_b > 0.0_dp) has_bracket = .false.
                if (abs(bracket_u_a - bracket_u_b) < min_secant_du) has_bracket = .false.
            end if
        end subroutine update_bracket

        subroutine choose_new_inlet_target(u_current, v_current, u_new, du_unlimited, du_final)
            real(dp), intent(in)  :: u_current, v_current
            real(dp), intent(out) :: u_new, du_unlimited, du_final
            real(dp) :: dU, dV, response_slope, secant_step, bracket_target
            real(dp) :: gain_effective, max_fraction_effective, min_step_effective

            if (emergency_mode) then
                gain_effective = controller_gain_capture
                max_fraction_effective = controller_max_fraction_emergency
                min_step_effective = max(min_abs_velocity_step_capture, &
                    emergency_min_fraction * max(abs(u_current), min_abs_velocity_step))
            else if (capture_mode) then
                gain_effective = controller_gain_capture
                max_fraction_effective = controller_max_fraction_capture
                min_step_effective = min_abs_velocity_step_capture
            else
                gain_effective = adaptive_gain
                max_fraction_effective = controller_max_fraction
                min_step_effective = min_abs_velocity_step
            end if

            max_velocity_step = max(max_fraction_effective * &
                max(abs(u_current), min_abs_velocity_step), min_step_effective)

            du_unlimited = feedback_sign_default * gain_effective * v_current
            if (emergency_mode) then
!                ! When the flame is close to leaving the useful domain, the controller
!                ! must reduce the inflow even if the filtered velocity estimate is
!                ! temporarily small.  The direction still follows feedback_sign_default.
                if (abs(du_unlimited) < min_step_effective) then
                    du_unlimited = sign(min_step_effective, feedback_sign_default * max(abs(v_current), velocity_tolerance_on))
                end if
            else if (has_bracket .and. .not. capture_mode) then
                bracket_target = 0.5_dp * (bracket_u_a + bracket_u_b)
                du_unlimited = bracket_target - u_current
            else if (use_secant_control .and. have_previous_control_point .and. .not. capture_mode) then
                dU = u_current - previous_control_inlet_velocity
                dV = v_current - previous_control_velocity
                if (abs(dU) >= min_secant_du .and. abs(dV) > velocity_tolerance_off) then
                    response_slope = dV / dU
!                    ! Accept the secant only if the observed response has the sign required
!                    ! by the stabilizing feedback.  Otherwise it is most likely contaminated
!                    ! by delayed flame/acoustic transients and would turn the controller into
!                    ! positive feedback.
                    if (response_slope * feedback_sign_default < 0.0_dp .and. &
                        abs(response_slope) > 1.0e-12_dp .and. &
                        abs(response_slope) < max_response_slope_abs) then
                        secant_step = -secant_relaxation * v_current / response_slope
                        if (abs(secant_step) <= 5.0_dp * max_velocity_step) then
                            du_unlimited = secant_step
                        end if
                    end if
                end if
            end if

!            ! Even a capture/bisection/secant target is passed through a step limiter.
!            ! The applied boundary value is additionally ramped in physical time.
            du_final = min(max(du_unlimited, -max_velocity_step), max_velocity_step)

!            ! Persistent nonzero drift must produce a finite control action.  This lower
!            ! bound is active only outside the velocity/position deadband.
            if (abs(control_velocity) > persistent_error_factor * velocity_tolerance_on .or. &
                abs(position_control_error) > 0.0_dp) then
                if (abs(du_final) < min_step_effective) then
                    if (du_unlimited /= 0.0_dp) then
                        du_final = sign(min_step_effective, du_unlimited)
                    else
                        du_final = feedback_sign_default * sign(min_step_effective, v_current)
                    end if
                end if
            end if

            u_new = max(u_current + du_final, 0.0_dp)
        end subroutine choose_new_inlet_target

        subroutine write_flamelet_tables_once()
            if (.not. flamelet_output_written) then
                call this%chem_kin_solver%write_chemical_kinetics_table(chem_table_filename)
                call this%write_data_table(data_table_filename)
                flamelet_output_written = .true.
                final_output_written = .true.
            end if
        end subroutine write_flamelet_tables_once

        subroutine write_laminar_velocity_once()
            integer :: sl_unit
            character(len=200) :: sl_file_name

            if (.not. sl_output_written) then
                write(sl_file_name,'(A,I0,A)') 'laminar_flame_velocity_', this%load_counter, '.dat'
                open(newunit = sl_unit, file = sl_file_name, status = 'replace', form = 'formatted')
                write(sl_unit,'(A)') &
                    'VARIABLES="time" "U_in" "Vfl_measurement" "SL_disp" "lin_R2" "lin_RMS" "split_dV" "attempt"'
                write(sl_unit,'(100E20.12)') time, measurement_inlet_velocity, measurement_velocity_save, &
                    sl_displacement_save, measurement_r2_save, measurement_rms_save, &
                    measurement_split_slope_diff_save, real(measurement_attempt,dp)
                close(sl_unit)
                sl_output_written = .true.
            end if
        end subroutine write_laminar_velocity_once

        subroutine start_drift_measurement_ramp(t_now)
            real(dp), intent(in) :: t_now

            stabilized_inlet_velocity = inlet_velocity_target
            measurement_attempt = 1
            current_measurement_delta = max(measurement_delta_fraction * &
                max(abs(stabilized_inlet_velocity), min_abs_velocity_step), measurement_delta_min)
            current_measurement_delta = min(current_measurement_delta, measurement_delta_max_fraction * &
                max(abs(stabilized_inlet_velocity), min_abs_velocity_step))
            measurement_inlet_velocity = max(stabilized_inlet_velocity + &
                measurement_delta_sign * current_measurement_delta, 0.0_dp)
            inlet_velocity_target = measurement_inlet_velocity
            ramp_start_velocity = inlet_velocity_applied
            ramp_start_time = t_now
            active_inlet_ramp_time = measurement_ramp_time
            sl_displacement_save = 0.0_dp
            measurement_velocity_save = 0.0_dp
            measurement_r2_save = 0.0_dp
            measurement_rms_save = 0.0_dp
            measurement_split_slope_diff_save = 0.0_dp
            control_stage = stage_measurement_ramp
            stabilization_counter = 0
            post_flamelet_hold_counter = 0
            call clear_control_history()
        end subroutine start_drift_measurement_ramp

        subroutine write_tracking_line()
            real(dp) :: measurement_flag, bracket_flag, capture_flag, emergency_flag, flame_detected_flag

            if (measurement_enabled) then
                measurement_flag = 1.0_dp
            else
                measurement_flag = 0.0_dp
            end if

            if (has_bracket) then
                bracket_flag = 1.0_dp
            else
                bracket_flag = 0.0_dp
            end if

            if (capture_mode) then
                capture_flag = 1.0_dp
            else
                capture_flag = 0.0_dp
            end if

            if (emergency_mode) then
                emergency_flag = 1.0_dp
            else
                emergency_flag = 0.0_dp
            end if

            if (flame_detected) then
                flame_detected_flag = 1.0_dp
            else
                flame_detected_flag = 0.0_dp
            end if

            write(flame_loc_unit,'(100E20.12)') &
                time, current_flame_location(1:dimensions), &
                flame_velocity_lsq, flame_velocity_filtered, diag_flame_velocity_lsq, &
                diag_flame_velocity_filtered, measurement_velocity_save, &
                control_velocity, position_error, front_reference_coord, &
                inlet_velocity_applied, inlet_velocity_target, &
                target_step_for_log, adaptive_gain, real(hist_count,dp), measurement_flag, bracket_flag, &
                capture_flag, emergency_flag, flame_detected_flag, &
                front_spread, heat_release_integral, heat_release_max, heat_release_valid_limit, H_max, Tgrad_max, &
                real(control_stage,dp), sl_displacement, linear_r2, linear_rms, split_slope_diff, &
                real(correction_counter,dp), real(stabilization_counter,dp), real(active_track_number,dp)
        end subroutine write_tracking_line

    end subroutine stabilizing_inlet_1D
    !--------------------------------------------------------------------------
    ! Populate thermochemical and velocity ghost cells for walls, inlets, and
    ! outlets before divergence and pressure assembly.  Conductive walls impose
    ! temperature while preserving pressure through the density ghost state.
    !--------------------------------------------------------------------------
    subroutine apply_boundary_conditions(this,time_step,predictor)
        class(fds_solver)    ,intent(inout)    :: this
        real(dp)            ,intent(in)        :: time_step
        logical, intent(in) :: predictor

        real(dp), dimension(this%chem%chem_ptr%species_number) :: mole_fractions

        real(dp)                    :: wall_temperature, farfield_temperature, farfield_pressure, farfield_velocity
        real(dp)                    :: farfield_density, average_molar_mass, normal_face_velocity
        real(dp), dimension(:), allocatable :: farfield_concentrations
        character(len=10), dimension(:), allocatable :: farfield_species_names
        logical                     :: outflow_flag

        integer :: dimensions, species_number
        integer, dimension(3,2) :: cons_inner_loop

        integer :: bound_number, sign
        integer :: i,j,k,plus,dim,dim1,spec,specie_index

        character(len=20) :: boundary_type_name

        if (time_step <= 0.0_dp) error stop 'FDS boundary update: non-positive time step'
        dimensions        = this%domain%get_domain_dimensions()
        species_number    = this%chem%chem_ptr%species_number

        cons_inner_loop    = this%domain%get_local_inner_cells_bounds()


        associate (    rho                => this%rho%s_ptr            , &
                    rho_int            => this%rho_int%s_ptr        , &
                    rho_old            => this%rho_old%s_ptr        , &
                    Y                => this%Y%v_ptr                , &
                    Y_int            => this%Y_int%v_ptr            , &
                    Y_old            => this%Y_old%v_ptr            , &
                    v                => this%v%v_ptr                , &
                    v_f                => this%v_f%v_ptr            , &
                    h_s                => this%h_s%s_ptr            , &
                    T                => this%T%s_ptr                , &
                    bc                => this%boundary%bc_ptr)

            do k = cons_inner_loop(3,1),cons_inner_loop(3,2)
            do j = cons_inner_loop(2,1),cons_inner_loop(2,2)
            do i = cons_inner_loop(1,1),cons_inner_loop(1,2)
                if(bc%bc_markers(i,j,k) == 0) then

                    do dim = 1,dimensions
                        do plus = 1,2
                            sign            = (-1)**plus
                            bound_number    = bc%bc_markers(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3))
                            if( bound_number /= 0 ) then
                                boundary_type_name = bc%boundary_types(bound_number)%get_type_name()

                                select case(boundary_type_name)
                                    case('wall')
                                        if(bc%boundary_types(bound_number)%is_conductive()) then
                                            wall_temperature = bc%boundary_types(bound_number)%get_wall_temperature()
                                            T%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3)) = wall_temperature
                                        else
                                            T%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3)) = T%cells(i,j,k)
                                        end if

                                        if (predictor) then
                                            rho_int%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3))    = &
                                                rho_int%cells(i,j,k) * T%cells(i,j,k) / &
                                                T%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3))
                                            rho_old%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3))    = &
                                                rho_old%cells(i,j,k) * T%cells(i,j,k) / &
                                                T%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3))
                                            rho%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3))        = &
                                                rho_int%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3))
                                            do spec = 1, species_number
                                                Y_int%pr(spec)%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3)) &
                                                    = Y_int%pr(spec)%cells(i,j,k)
                                                Y_old%pr(spec)%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3)) = &
                                                Y_old%pr(spec)%cells(i,j,k)
                                                Y%pr(spec)%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3)) &
                                                    = &
                                                    Y_int%pr(spec)%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3))
                                            end do
                                        else
                                            rho%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3)) &
                                                = rho%cells(i,j,k) * T%cells(i,j,k) / &
                                                T%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3))
                                            do spec = 1, species_number
                                                Y%pr(spec)%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3)) &
                                                    = Y%pr(spec)%cells(i,j,k)
                                            end do
                                        end if

                                        average_molar_mass = 0.0_dp
                                        do spec = 1, species_number
                                            average_molar_mass = average_molar_mass + &
                                                Y%pr(spec)%cells(i,j,k)/this%thermo%thermo_ptr%molar_masses(spec)
                                        end do
                                        average_molar_mass = 1.0_dp/average_molar_mass

                                        do spec = 1, species_number
                                            mole_fractions(spec) = &
                                                Y%pr(spec)%cells(i,j,k)*average_molar_mass/ &
                                                this%thermo%thermo_ptr%molar_masses(spec)
                                        end do
                                        h_s%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3)) = &
                                            (this%thermo%thermo_ptr%mixture_enthalpy_molar(T%cells(i+sign*I_m(dim,1),j+sign* &
                                            I_m(dim,2),k+sign*I_m(dim,3)),mole_fractions) &
                                                                                                            - &
                                                this%thermo%thermo_ptr%mixture_enthalpy_molar( &
                                                T_ref,mole_fractions))/average_molar_mass

                                        do dim1 = 1, dimensions
                                            if(dim1 == dim) then
                                            v%pr(dim1)%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3)) &
                                                = -v%pr(dim1)%cells(i,j,k)
                                            else
                                                v%pr(dim1)%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3)) &
                                                    =  v%pr(dim1)%cells(i,j,k)
                                            end if
                                        end do

                                    case('outlet')
                                        normal_face_velocity = v_f%pr(dim)%cells(dim, i+max(sign,0)*I_m(dim,1), &
                                            j+max(sign,0)*I_m(dim,2), k+max(sign,0)*I_m(dim,3))
                                        outflow_flag = (sign*normal_face_velocity >= 0.0_dp)

                                        if (outflow_flag) then
                                            ! Convective outflow: extrapolate thermodynamic and transport states.
                                            T%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3)) = T%cells(i,j,k)
                                            h_s%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3)) = h_s%cells(i,j,k)

                                            rho%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3)) = rho%cells(i,j,k)
                                            rho_int%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3)) = &
                                                rho_int%cells(i,j,k)
                                            rho_old%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3)) = &
                                                rho_old%cells(i,j,k)

                                            do spec = 1, species_number
                                                Y%pr(spec)%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3)) = &
                                                    Y%pr(spec)%cells(i,j,k)
                                                Y_int%pr(spec)%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3)) = &
                                                    Y_int%pr(spec)%cells(i,j,k)
                                                Y_old%pr(spec)%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3)) = &
                                                    Y_old%pr(spec)%cells(i,j,k)
                                            end do

                                            do dim1 = 1, dimensions
                                                v%pr(dim1)%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3)) = &
                                                    v%pr(dim1)%cells(i,j,k)
                                            end do
                                        else
                                            ! Backflow through a pressure outlet: restore the prescribed ambient state.
                                            farfield_temperature = bc%boundary_types(bound_number)%get_farfield_temperature()
                                            farfield_density     = bc%boundary_types(bound_number)%get_farfield_density()
                                            farfield_velocity    = bc%boundary_types(bound_number)%get_farfield_velocity()

                                            if (allocated(farfield_species_names)) deallocate(farfield_species_names)
                                            if (allocated(farfield_concentrations)) deallocate(farfield_concentrations)
                                            call &
                                                bc%boundary_types(bound_number) &
                                                %get_farfield_species_names(farfield_species_names)
                                            call &
                                                bc%boundary_types(bound_number) &
                                                %get_farfield_concentrations(farfield_concentrations)

                                            mole_fractions = 0.0_dp
                                            do spec = 1, size(farfield_species_names)
                                                specie_index = &
                                                    this%chem%chem_ptr%get_chemical_specie_index(farfield_species_names(spec))
                                                if ((specie_index >= 1).and.(specie_index <= species_number)) then
                                                    mole_fractions(specie_index) = farfield_concentrations(spec)
                                                end if
                                            end do
                                            call this%thermo%thermo_ptr%change_cell_units_mole_to_dimless(mole_fractions)

                                            if (farfield_density <= tiny(1.0_dp)) farfield_density = rho%cells(i,j,k)

                                            T%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3)) = &
                                                farfield_temperature
                                            rho%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3)) = farfield_density
                                            rho_int%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3)) = &
                                                farfield_density
                                            rho_old%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3)) = &
                                                farfield_density

                                            do spec = 1, species_number
                                                Y%pr(spec)%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3)) = &
                                                    mole_fractions(spec)
                                                Y_int%pr(spec)%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3)) = &
                                                    mole_fractions(spec)
                                                Y_old%pr(spec)%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3)) = &
                                                    mole_fractions(spec)
                                            end do

                                            average_molar_mass = 0.0_dp
                                            do spec = 1, species_number
                                                average_molar_mass = average_molar_mass + &
                                                    mole_fractions(spec)/this%thermo%thermo_ptr%molar_masses(spec)
                                            end do
                                            average_molar_mass = 1.0_dp/max(average_molar_mass,tiny(1.0_dp))

                                            ! Convert mass fractions to mole fractions for JANAF mixture enthalpy.
                                            do spec = 1, species_number
                                                mole_fractions(spec) = &
                                                    mole_fractions(spec)*average_molar_mass/ &
                                                    this%thermo%thermo_ptr%molar_masses(spec)
                                            end do
                                            h_s%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3)) = &
                                                (this%thermo%thermo_ptr%mixture_enthalpy_molar(farfield_temperature, &
                                                    mole_fractions) &
                                                -this%thermo%thermo_ptr%mixture_enthalpy_molar(T_ref,mole_fractions))/ &
                                                    average_molar_mass

                                            do dim1 = 1, dimensions
                                                if (dim1 == dim) then
                                                    v%pr(dim1)%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3)) = &
                                                        farfield_velocity
                                                else
                                                    v%pr(dim1)%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3)) = &
                                                        0.0_dp
                                                end if
                                            end do
                                        end if
                                    case('inlet')
                                        if(sign == 1) then
                                            if(v_f%pr(dim)%cells(dim,i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3)) > 0.0_dp) then
                                                T%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3))    = &
                                                    T%cells(i,j,k)
                                            else
                                                farfield_temperature                                            = &
                                                    this%boundary%bc_ptr%boundary_types(bound_number)%get_farfield_temperature()
                                                T%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3))    = &
                                                    farfield_temperature
                                            end if
                                        else
                                            if(v_f%pr(dim)%cells(dim,i,j,k) < 0.0_dp) then
                                                T%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3))    = &
                                                    T%cells(i,j,k)
                                            else
                                                farfield_temperature                                            = &
                                                    this%boundary%bc_ptr%boundary_types(bound_number)%get_farfield_temperature()
                                                T%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3))    = &
                                                    farfield_temperature
                                            end if
                                        end if

                                         farfield_velocity        =  this%inlet_velocity

                                        do dim1 = 1, dimensions
                                            if(dim1 == dim) then
                                                v%pr(dim1)%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3)) &
                                                    =  farfield_velocity
                                            else
                                                v%pr(dim1)%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3)) &
                                                    =  0.0_dp
                                            end if
                                        end do

                                        if (predictor) then
                                            rho_int%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3))    = &
                                                rho%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3))
                                            rho_old%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3))    = &
                                                rho%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3))
                                            do spec = 1, species_number
                                                Y_int%pr(spec)%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3)) &
                                                    = Y%pr(spec)%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3))
                                                Y_old%pr(spec)%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3)) &
                                                    = Y%pr(spec)%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3))
                                            end do
                                        end if
                                end select
                            end if
                        end do
                    end do
                end if
            end do
            end do
            end do

        end associate

    end subroutine
    pure function charm_face_value(scalar_array,velocity) result(phi)
        real(dp), dimension(4), intent(in) :: scalar_array
        real(dp), intent(in) :: velocity
        real(dp) :: phi
        real(dp) :: B_r, phi_loc, phi_up, r, s

        phi_up = 0.0_dp
        B_r = 0.0_dp
        phi_loc = scalar_array(3) - scalar_array(2)
        if (velocity > 0.0_dp) then
            phi_up = scalar_array(2) - scalar_array(1)
        else
            phi_up = scalar_array(4) - scalar_array(3)
        end if
        if (phi_loc /= 0.0_dp) then
            r = phi_up / phi_loc
        else
            r = 0.0_dp
        end if
        if (r > 1.0e-010_dp) then
            s = 1.0_dp / r
            B_r = s * (3.0_dp*s + 1.0_dp) / (s + 1.0_dp) / (s + 1.0_dp)
        else
            B_r = 0.0_dp
        end if
        if (velocity > 0.0_dp) then
            phi = scalar_array(2) + B_r * 0.5_dp * phi_up
        else
            phi = scalar_array(3) - B_r * 0.5_dp * phi_up
        end if
    end function charm_face_value

    subroutine eos_corrected_species_face_vector(rho_field,Y_field,molar_masses,species_number,dim,i,j,k,face_side,velocity,phi)
        type(field_scalar_cons), intent(in) :: rho_field
        type(field_vector_cons), intent(in) :: Y_field
        real(dp), dimension(:), intent(in) :: molar_masses
        integer, intent(in) :: species_number, dim, i, j, k, face_side
        real(dp), intent(in) :: velocity
        real(dp), dimension(:), intent(out) :: phi

        real(dp), dimension(4) :: scalar_array, rho_over_w_array
        real(dp) :: rhoY_max, rhoY_loc, rho_over_w_face, sum_not_gamma
        integer :: s, n, gamma, offset0, offset, gamma_offset
        integer :: ii, jj, kk

        if (face_side > 0) then
            offset0 = -1
        else
            offset0 = -2
        end if

        if (velocity > 0.0_dp) then
            gamma_offset = offset0 + 1
        else
            gamma_offset = offset0 + 2
        end if

        gamma = 1
        rhoY_max = -huge(1.0_dp)
        ii = i + gamma_offset*I_m(dim,1)
        jj = j + gamma_offset*I_m(dim,2)
        kk = k + gamma_offset*I_m(dim,3)
        do s = 1,species_number
            rhoY_loc = rho_field%cells(ii,jj,kk) * Y_field%pr(s)%cells(ii,jj,kk)
            if (rhoY_loc > rhoY_max) then
                rhoY_max = rhoY_loc
                gamma = s
            end if
        end do

        rho_over_w_array = 0.0_dp
        do s = 1,species_number
            do n = 1,4
                offset = offset0 + n - 1
                ii = i + offset*I_m(dim,1)
                jj = j + offset*I_m(dim,2)
                kk = k + offset*I_m(dim,3)
                scalar_array(n) = rho_field%cells(ii,jj,kk) * Y_field%pr(s)%cells(ii,jj,kk)
                rho_over_w_array(n) = rho_over_w_array(n) + scalar_array(n) / molar_masses(s)
            end do
            phi(s) = charm_face_value(scalar_array,velocity)
        end do

        ! FDS-style EOS preservation with a positivity-safe active-set repair.
        ! The face composition must satisfy sum_s(phi_s/M_s)=(rho/M)_face.
        ! Correct the dominant upwind carrier when possible.  If the other
        ! reconstructed species already exceed the target molar density, scale
        ! those positive contributions proportionally and set the carrier to zero.
        rho_over_w_face = max(charm_face_value(rho_over_w_array,velocity), 0.0_dp)
        phi = max(phi, 0.0_dp)
        sum_not_gamma = 0.0_dp
        do s = 1,species_number
            if (s /= gamma) sum_not_gamma = sum_not_gamma + phi(s) / molar_masses(s)
        end do
        if (sum_not_gamma <= rho_over_w_face) then
            phi(gamma) = molar_masses(gamma) * (rho_over_w_face - sum_not_gamma)
        else if (sum_not_gamma > tiny(1.0_dp)) then
            do s = 1,species_number
                if (s /= gamma) phi(s) = phi(s) * rho_over_w_face / sum_not_gamma
            end do
            phi(gamma) = 0.0_dp
        else
            phi = 0.0_dp
        end if
    end subroutine eos_corrected_species_face_vector
    !--------------------------------------------------------------------------
    ! Compute explicit advection/divergence and molecular-transport restrictions.
    ! The low-Mach method has no acoustic CFL condition.
    !--------------------------------------------------------------------------

    !--------------------------------------------------------------------------
    ! This production helper replaces the former heuristic window selection,
    ! which could use uninitialized front indices for nearly uniform fields.
    !--------------------------------------------------------------------------
    subroutine write_data_table(this, table_file)
        class(fds_solver), intent(in) :: this
        character(len=*), intent(in) :: table_file

        integer, dimension(3,2) :: inner
        integer :: unit_id, dimensions, species_number
        integer :: i, j0, k0, dim, spec

        dimensions = this%domain%get_domain_dimensions()
        species_number = this%chem%chem_ptr%species_number
        inner = this%domain%get_local_inner_cells_bounds()
        j0 = inner(2,1)
        k0 = inner(3,1)

        open(newunit=unit_id, file=trim(table_file), status='replace', &
            action='write', form='formatted')
        write(unit_id,'(A)',advance='no') '# x T rho'
        do dim = 1, dimensions
            write(unit_id,'(A,I0)',advance='no') ' u', dim
        end do
        do spec = 1, species_number
            write(unit_id,'(A,I0)',advance='no') ' Y', spec
        end do
        write(unit_id,*)

        do i = inner(1,1), inner(1,2)
            write(unit_id,'(*(ES24.16,1X))') &
                this%mesh%mesh_ptr%mesh(1,i,j0,k0), &
                this%T%s_ptr%cells(i,j0,k0), &
                this%rho%s_ptr%cells(i,j0,k0), &
                (this%v%v_ptr%pr(dim)%cells(i,j0,k0), dim=1,dimensions), &
                (this%Y%v_ptr%pr(spec)%cells(i,j0,k0), spec=1,species_number)
        end do
        close(unit_id)
    end subroutine write_data_table
    subroutine calculate_time_step(this)
        class(fds_solver), intent(inout) :: this

        real(dp), parameter :: maximum_time_step = 5.0e-3_dp
        real(dp), parameter :: diffusion_safety = 0.25_dp
        real(dp), dimension(3) :: cell_size
        real(dp) :: advective_limit, transport_limit
        real(dp) :: advective_rate, inverse_length_squared
        real(dp) :: cp_mass, diffusivity_value, thermal_diffusivity
        real(dp) :: candidate
        integer, dimension(3,2) :: inner
        integer :: dimensions, species_number
        integer :: i, j, k, dim, spec

        dimensions = this%domain%get_domain_dimensions()
        species_number = this%chem%chem_ptr%species_number
        inner = this%domain%get_local_inner_cells_bounds()
        cell_size = this%mesh%mesh_ptr%get_cell_edges_length()

        inverse_length_squared = 0.0_dp
        do dim = 1, dimensions
            if (cell_size(dim) <= tiny(1.0_dp)) then
                error stop 'FDS time step: non-positive cell size'
            end if
            inverse_length_squared = inverse_length_squared + 1.0_dp/cell_size(dim)**2
        end do

        ! Explicit low-Mach transport limit.  The target divergence is included
        ! because local expansion/compression changes a cell volume even when the
        ! translational velocity is small.
        advective_limit = huge(1.0_dp)
        associate (velocity => this%v%v_ptr, target_divergence => this%div_v_int%s_ptr, &
            bc => this%boundary%bc_ptr)
        !$omp parallel do collapse(3) schedule(static) &
        !$omp& private(i,j,k,dim,advective_rate) reduction(min:advective_limit)
        do k = inner(3,1), inner(3,2)
            do j = inner(2,1), inner(2,2)
                do i = inner(1,1), inner(1,2)
                    if (bc%bc_markers(i,j,k) /= 0) cycle
                    advective_rate = abs(target_divergence%cells(i,j,k))
                    do dim = 1, dimensions
                        advective_rate = advective_rate + &
                            abs(velocity%pr(dim)%cells(i,j,k))/cell_size(dim)
                    end do
                    if (advective_rate > tiny(1.0_dp)) then
                        advective_limit = min(advective_limit,1.0_dp/advective_rate)
                    end if
                end do
            end do
        end do
        !$omp end parallel do
        end associate

        transport_limit = huge(1.0_dp)

        if (this%diffusion_flag) then
            associate (diffusivity => this%D%v_ptr, bc => this%boundary%bc_ptr)
            !$omp parallel do collapse(3) schedule(static) &
            !$omp& private(i,j,k,spec,diffusivity_value) reduction(min:transport_limit)
            do k = inner(3,1), inner(3,2)
                do j = inner(2,1), inner(2,2)
                    do i = inner(1,1), inner(1,2)
                        if (bc%bc_markers(i,j,k) /= 0) cycle
                        do spec = 1, species_number
                            diffusivity_value = diffusivity%pr(spec)%cells(i,j,k)
                            if (diffusivity_value > tiny(1.0_dp)) then
                                transport_limit = min(transport_limit, &
                                    diffusion_safety/(diffusivity_value*inverse_length_squared))
                            end if
                        end do
                    end do
                end do
            end do
            !$omp end parallel do
            end associate
        end if

        if (this%viscosity_flag) then
            associate (density => this%rho%s_ptr, viscosity => this%nu%s_ptr, &
                bc => this%boundary%bc_ptr)
            !$omp parallel do collapse(3) schedule(static) &
            !$omp& private(i,j,k,candidate) reduction(min:transport_limit)
            do k = inner(3,1), inner(3,2)
                do j = inner(2,1), inner(2,2)
                    do i = inner(1,1), inner(1,2)
                        if (bc%bc_markers(i,j,k) /= 0) cycle
                        if (density%cells(i,j,k) <= tiny(1.0_dp)) cycle
                        candidate = viscosity%cells(i,j,k)/density%cells(i,j,k)
                        if (candidate > tiny(1.0_dp)) then
                            transport_limit = min(transport_limit, &
                                diffusion_safety/(candidate*inverse_length_squared))
                        end if
                    end do
                end do
            end do
            !$omp end parallel do
            end associate
        end if

        if (this%heat_trans_flag) then
            associate (density => this%rho%s_ptr, conductivity => this%kappa%s_ptr, &
                cp_molar => this%mixture_cp%s_ptr, molar_mass => this%mix_mol_mass%s_ptr, &
                bc => this%boundary%bc_ptr)
            !$omp parallel do collapse(3) schedule(static) &
            !$omp& private(i,j,k,cp_mass,thermal_diffusivity) reduction(min:transport_limit)
            do k = inner(3,1), inner(3,2)
                do j = inner(2,1), inner(2,2)
                    do i = inner(1,1), inner(1,2)
                        if (bc%bc_markers(i,j,k) /= 0) cycle
                        if (density%cells(i,j,k) <= tiny(1.0_dp)) cycle
                        if (molar_mass%cells(i,j,k) <= tiny(1.0_dp)) cycle
                        cp_mass = cp_molar%cells(i,j,k)/molar_mass%cells(i,j,k)
                        if (cp_mass <= tiny(1.0_dp)) cycle
                        thermal_diffusivity = conductivity%cells(i,j,k)/ &
                            (density%cells(i,j,k)*cp_mass)
                        if (thermal_diffusivity > tiny(1.0_dp)) then
                            transport_limit = min(transport_limit, &
                                diffusion_safety/(thermal_diffusivity*inverse_length_squared))
                        end if
                    end do
                end do
            end do
            !$omp end parallel do
            end associate
        end if

        candidate = min(maximum_time_step, transport_limit, &
            this%courant_fraction*advective_limit)
        if (candidate <= 0.0_dp .or. candidate >= 0.5_dp*huge(1.0_dp)) then
            candidate = min(maximum_time_step,this%initial_time_step)
        end if
        this%time_step = candidate
    end subroutine calculate_time_step
    pure function get_time_step(this)
        real(dp)                        :: get_time_step
        class(fds_solver)    ,intent(in)        :: this

        get_time_step = this%time_step
    end function

    pure function get_time(this)
        real(dp)                        :: get_time
        class(fds_solver)    ,intent(in)        :: this

        get_time = this%time
    end function
end module

