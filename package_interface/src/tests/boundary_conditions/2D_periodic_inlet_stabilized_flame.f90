!===============================================================================
! 2D inlet-stabilized hydrogen-air flame interface
!
! Production-oriented replacement for the legacy
! 2D_counterflow_channel_flame.f90 interface.
!
! Main changes:
!   * single namelist-driven case instead of six nested task loops;
!   * portable case-directory handling (no IFPORT/xcopy dependency);
!   * explicit chemistry backend selection, including QSS2;
!   * mechanism-specific KEROMNES thermo/transport files;
!   * no runtime dependency on molar_masses.dat;
!   * spatially uniform thermodynamic pressure for the low-Mach solver;
!   * element-conserving complete-combustion initialization for the hot zone;
!   * species are addressed by name, never by hard-coded mechanism indices;
!   * deterministic/reproducible front perturbation;
!   * active inlet anchoring through problem_controls/flame_stabilization_solver;
!   * Soret diffusion enabled by default for the production case.
!
! NOTE:
! The flame-stabilization implementation already computes multidimensional
! heat-release/temperature-gradient centroids and updates all inlet boundary
! types.  The current problem-control validation guard still blocks dimensions
! /= 1, so apply the accompanying validator patch before building this 2D case.
!===============================================================================

program package_interface

    use, intrinsic :: iso_c_binding, only: c_int, c_size_t, c_char, c_ptr, &
        c_null_char, c_associated
    use kind_parameters
    use global_data
    use nrg_build_info, only: write_nrg_source_revision
    use computational_domain_class
    use chemical_properties_class
    use thermophysical_properties_class
    use solver_options_class
    use problem_control_class
    use computational_mesh_class
    use mpi_communications_class
    use data_manager_class
    use boundary_conditions_class
    use field_scalar_class
    use field_vector_class
    use data_save_class
    use data_io_class
    use post_processor_manager_class

    implicit none

    interface
#ifdef WIN
        function c_getcwd(buffer, maxlen) bind(C, name="_getcwd") result(ptr)
            import :: c_ptr, c_char, c_int
            character(kind=c_char) :: buffer(*)
            integer(c_int), value :: maxlen
            type(c_ptr) :: ptr
        end function c_getcwd

        function c_chdir(path) bind(C, name="_chdir") result(status)
            import :: c_int, c_char
            character(kind=c_char), intent(in) :: path(*)
            integer(c_int) :: status
        end function c_chdir
#else
        function c_getcwd(buffer, maxlen) bind(C, name="getcwd") result(ptr)
            import :: c_ptr, c_char, c_size_t
            character(kind=c_char) :: buffer(*)
            integer(c_size_t), value :: maxlen
            type(c_ptr) :: ptr
        end function c_getcwd

        function c_chdir(path) bind(C, name="chdir") result(status)
            import :: c_int, c_char
            character(kind=c_char), intent(in) :: path(*)
            integer(c_int) :: status
        end function c_chdir
#endif
    end interface

    ! NRG objects.
    type(computational_domain)                :: problem_domain
    type(data_manager)                        :: problem_data_manager
    type(mpi_communications)                  :: problem_mpi_support
    type(chemical_properties), target         :: problem_chemistry
    type(thermophysical_properties), target   :: problem_thermophysics
    type(solver_options)                      :: problem_solver_options
    type(flame_stabilization_control)         :: problem_flame_control
    type(problem_controls)                    :: problem_controls_setup
    type(computational_mesh), target          :: problem_mesh
    type(boundary_conditions), target         :: problem_boundaries
    type(field_scalar_cons), target           :: p, T, rho
    type(field_vector_cons), target           :: v, Y
    type(post_processor_manager)              :: problem_post_proc_manager
    type(data_io)                             :: problem_data_io
    type(data_save)                           :: problem_data_save

    ! Case identity.
    character(len=64)  :: case_id
    character(len=512) :: results_root
    character(len=1024):: work_dir
    character(len=4096):: initial_work_dir
    character(len=4096):: config_file
    logical            :: overwrite_existing_case

    namelist /case_config/ case_id, results_root, overwrite_existing_case

    ! Physical mixture/state.
    real(dp) :: hydrogen_mole_percent
    real(dp) :: n2_o2_molar_ratio
    real(dp) :: operating_pressure
    real(dp) :: fresh_temperature
    real(dp) :: burned_temperature
    real(dp) :: inlet_velocity

    namelist /mixture_config/ hydrogen_mole_percent, n2_o2_molar_ratio, &
        operating_pressure, fresh_temperature, burned_temperature, inlet_velocity

    ! Geometry / perturbation.
    real(dp) :: domain_length
    real(dp) :: domain_width
    real(dp) :: cell_size_xy
    real(dp) :: front_position_fraction
    real(dp) :: perturbation_amplitude
    real(dp) :: perturbation_phase
    real(dp) :: random_jitter_amplitude
    integer  :: perturbation_mode
    integer  :: random_seed_value
    logical  :: transverse_periodic

    namelist /geometry_config/ domain_length, domain_width, cell_size_xy, &
        front_position_fraction, perturbation_amplitude, perturbation_mode, &
        perturbation_phase, random_jitter_amplitude, random_seed_value, &
        transverse_periodic

    ! Physics / transport.
    character(len=32) :: mechanism_id
    logical :: molecular_diffusion_flag
    logical :: soret_diffusion_flag
    logical :: viscosity_flag
    logical :: thermal_radiation_flag

    namelist /physics_config/ mechanism_id, molecular_diffusion_flag, &
        soret_diffusion_flag, viscosity_flag, thermal_radiation_flag

    ! Low-Mach numerics.
    real(dp) :: cfl_coefficient
    real(dp) :: initial_time_step

    namelist /numerics_config/ cfl_coefficient, initial_time_step

    ! Adaptive inlet/flame anchoring.
    character(len=40) :: flame_stabilization_mode
    character(len=240) :: scientific_result_title
    character(len=80) :: flame_case_setup
    real(dp) :: flame_time_delay
    real(dp) :: flame_time_track
    real(dp) :: flame_time_control
    real(dp) :: flame_response_settle_time
    real(dp) :: flame_measurement_max_duration
    real(dp) :: flame_measurement_min_displacement_cells

    namelist /flame_stabilization_config/ flame_stabilization_mode, &
        scientific_result_title, flame_case_setup, flame_time_delay, &
        flame_time_track, flame_time_control, flame_response_settle_time, &
        flame_measurement_max_duration, flame_measurement_min_displacement_cells

    ! Chemistry backend.
    character(len=24) :: chemistry_backend
    real(dp) :: chemistry_cvode_relative_tolerance
    real(dp) :: chemistry_cvode_absolute_tolerance
    integer  :: chemistry_cvode_max_steps
    real(dp) :: chemistry_qss1_relative_change_limit
    real(dp) :: chemistry_qss1_minimum_internal_step
    real(dp) :: chemistry_qss1_active_concentration_fraction
    real(dp) :: chemistry_qss1_step_growth_factor
    integer  :: chemistry_qss1_max_steps
    real(dp) :: chemistry_qss2_error_tolerance
    real(dp) :: chemistry_qss2_minimum_internal_step
    real(dp) :: chemistry_qss2_active_concentration_fraction
    integer  :: chemistry_qss2_max_steps

    namelist /chemistry_backend_config/ chemistry_backend, &
        chemistry_cvode_relative_tolerance, &
        chemistry_cvode_absolute_tolerance, chemistry_cvode_max_steps, &
        chemistry_qss1_relative_change_limit, &
        chemistry_qss1_minimum_internal_step, &
        chemistry_qss1_active_concentration_fraction, &
        chemistry_qss1_step_growth_factor, chemistry_qss1_max_steps, &
        chemistry_qss2_error_tolerance, &
        chemistry_qss2_minimum_internal_step, &
        chemistry_qss2_active_concentration_fraction, &
        chemistry_qss2_max_steps

    ! Output.
    real(dp) :: postprocess_interval_us
    real(dp) :: field_save_interval_us
    real(dp) :: restart_check_interval_ms
    real(dp) :: wall_time_output_limit_min

    namelist /output_config/ postprocess_interval_us, field_save_interval_us, &
        restart_check_interval_ms, wall_time_output_limit_min

    ! Derived problem data.
    character(len=32) :: mech_name
    character(len=64) :: mech_file, thermo_file, transdata_file
    real(dp) :: x_h2, nu
    real(dp) :: fresh_h2, fresh_o2, fresh_n2
    real(dp) :: prod_h2, prod_o2, prod_n2, prod_h2o
    real(dp) :: front_shift, random_value
    real(dp) :: pi_local
    integer :: nx, ny, i, j
    integer :: i_h2, i_o2, i_n2, i_h2o
    integer :: front_index, base_front_index
    integer :: species_number
    integer :: log_unit, io_unit, ierr
    integer :: transducer_offset
    integer, dimension(3,2) :: utter_loop
    integer, dimension(3,2) :: observation_slice
    integer, dimension(3,2) :: summation_region
    real(dp), dimension(3) :: mesh_cell_size
    logical :: stop_flag, config_exists, case_exists

    !--------------------------------------------------------------------------
    ! Defaults: the requested 15% H2 production pilot.
    !
    ! The legacy interface initialized the bulk and inlet at 5 atm but used a
    ! 1-atm outlet.  The low-Mach formulation uses a spatially uniform
    ! thermodynamic pressure, so this replacement uses one operating pressure
    ! everywhere.  Change operating_pressure in the namelist for a 1-atm case.
    !--------------------------------------------------------------------------
    case_id = 'h2_15_qss2_prod'
    results_root = '2D_inlet_stabilized_flame_cases'
    overwrite_existing_case = .false.

    hydrogen_mole_percent = 15.0_dp
    n2_o2_molar_ratio = 3.762_dp
    operating_pressure = 5.0_dp*101325.0_dp
    fresh_temperature = 300.0_dp
    burned_temperature = 1500.0_dp
    inlet_velocity = 0.30_dp

    domain_length = 6.4e-3_dp
    domain_width = 6.4e-3_dp
    cell_size_xy = 5.0e-5_dp
    front_position_fraction = 0.65_dp
    perturbation_amplitude = 200.0e-6_dp
    perturbation_mode = 1
    perturbation_phase = 0.0_dp
    random_jitter_amplitude = 100.0e-6_dp
    random_seed_value = 20261005
    transverse_periodic = .true.

    mechanism_id = 'keromnes'
    molecular_diffusion_flag = .true.
    soret_diffusion_flag = .true.
    viscosity_flag = .true.
    thermal_radiation_flag = .false.

    cfl_coefficient = 0.25_dp
    initial_time_step = 1.0e-7_dp

    flame_stabilization_mode = 'anchor'
    scientific_result_title = &
        '2D inlet-stabilized 15 percent H2-air flame at 5 atm'
    flame_case_setup = '2D_inlet_stabilized'
    flame_time_delay = 1.0e-3_dp
    flame_time_track = 1.0e-4_dp
    flame_time_control = 5.0e-4_dp
    flame_response_settle_time = 5.0e-4_dp
    flame_measurement_max_duration = 1.0_dp
    flame_measurement_min_displacement_cells = 0.25_dp

    chemistry_backend = 'qss2'
    chemistry_cvode_relative_tolerance = 1.0e-8_dp
    chemistry_cvode_absolute_tolerance = 1.0e-12_dp
    chemistry_cvode_max_steps = 100000
    chemistry_qss1_relative_change_limit = 2.0e-3_dp
    chemistry_qss1_minimum_internal_step = 1.0e-10_dp
    chemistry_qss1_active_concentration_fraction = 1.0e-7_dp
    chemistry_qss1_step_growth_factor = 1.04_dp
    chemistry_qss1_max_steps = 100000
    chemistry_qss2_error_tolerance = 1.0e-3_dp
    chemistry_qss2_minimum_internal_step = 1.0e-10_dp
    chemistry_qss2_active_concentration_fraction = 1.0e-7_dp
    chemistry_qss2_max_steps = 100000

    postprocess_interval_us = 5.0_dp
    field_save_interval_us = 250.0_dp
    restart_check_interval_ms = 1.0_dp
    wall_time_output_limit_min = 240.0_dp

    call get_current_directory(initial_work_dir)

    config_file = '2D_inlet_stabilized_flame.nml'
    if (command_argument_count() >= 1) then
        call get_command_argument(1, config_file)
    end if

    inquire(file=trim(config_file), exist=config_exists)
    if (.not. config_exists) then
        write(*,'(A)') 'ERROR: configuration file not found: '//trim(config_file)
        error stop 1
    end if

    open(newunit=io_unit, file=trim(config_file), status='old', &
        action='read', iostat=ierr)
    if (ierr /= 0) error stop 'Unable to open configuration namelist'

    rewind(io_unit)
    read(io_unit, nml=case_config, iostat=ierr)
    if (ierr /= 0) call namelist_error('case_config', ierr)

    rewind(io_unit)
    read(io_unit, nml=mixture_config, iostat=ierr)
    if (ierr /= 0) call namelist_error('mixture_config', ierr)

    rewind(io_unit)
    read(io_unit, nml=geometry_config, iostat=ierr)
    if (ierr /= 0) call namelist_error('geometry_config', ierr)

    rewind(io_unit)
    read(io_unit, nml=physics_config, iostat=ierr)
    if (ierr /= 0) call namelist_error('physics_config', ierr)

    rewind(io_unit)
    read(io_unit, nml=numerics_config, iostat=ierr)
    if (ierr /= 0) call namelist_error('numerics_config', ierr)

    rewind(io_unit)
    read(io_unit, nml=flame_stabilization_config, iostat=ierr)
    if (ierr /= 0) call namelist_error('flame_stabilization_config', ierr)

    rewind(io_unit)
    read(io_unit, nml=chemistry_backend_config, iostat=ierr)
    if (ierr /= 0) call namelist_error('chemistry_backend_config', ierr)

    rewind(io_unit)
    read(io_unit, nml=output_config, iostat=ierr)
    if (ierr /= 0) call namelist_error('output_config', ierr)
    close(io_unit)

    call validate_configuration()

    select case (trim(to_lower(mechanism_id)))
    case ('keromnes')
        mech_name = 'KEROMNES'
        mech_file = 'KEROMNES.txt'
        thermo_file = 'KEROMNES_THERMO.txt'
        transdata_file = 'KEROMNES_TRANSDATA.txt'
    case default
        error stop 'This production interface currently supports mechanism_id=keromnes'
    end select

    nx = nint(domain_length/cell_size_xy)
    ny = nint(domain_width/cell_size_xy)

    if (abs(real(nx,dp)*cell_size_xy-domain_length) > &
        1.0e-10_dp*max(domain_length,1.0_dp)) then
        error stop 'domain_length must be an integer multiple of cell_size_xy'
    end if
    if (abs(real(ny,dp)*cell_size_xy-domain_width) > &
        1.0e-10_dp*max(domain_width,1.0_dp)) then
        error stop 'domain_width must be an integer multiple of cell_size_xy'
    end if

    work_dir = trim(results_root)//trim(fold_sep)//trim(case_id)
    inquire(file=trim(work_dir), exist=case_exists)
    if (case_exists .and. .not. overwrite_existing_case) then
        write(*,'(A)') 'ERROR: case directory already exists: '//trim(work_dir)
        write(*,'(A)') 'Set overwrite_existing_case=.true. only when intentional.'
        error stop 1
    end if
    if (case_exists .and. overwrite_existing_case) then
        call run_filesystem_command( &
            'cmake -E remove_directory "'//trim(work_dir)//'"')
    end if

    call ensure_directory(trim(work_dir))
    call copy_directory_tree(trim(task_setup_folder), &
        trim(work_dir)//trim(fold_sep)//trim(task_setup_folder))
    call run_filesystem_command('cmake -E copy_if_different "'// &
        trim(config_file)//'" "'//trim(work_dir)//trim(fold_sep)//'case_input.nml"')

    call change_directory(trim(work_dir))

    open(newunit=log_unit, file=problem_setup_log_file, status='replace', &
        form='formatted')
    call write_nrg_source_revision(log_unit, '2D inlet-stabilized flame setup')

    problem_domain = computational_domain_c( &
        dimensions=2, &
        cells_number=(/nx,ny,1/), &
        coordinate_system='cartesian', &
        lengths=reshape((/0.0_dp,0.0_dp,0.0_dp, &
            domain_length,domain_width,0.005_dp/),(/3,2/)), &
        axis_names=(/'x','y','z'/), &
        periodic=(/.false.,transverse_periodic,.false./))

    problem_chemistry = chemical_properties_c( &
        chemical_mechanism_file_name=mech_file, &
        default_enhanced_efficiencies=1.0_dp, &
        E_act_units='cal.mol')

    problem_thermophysics = thermophysical_properties_c( &
        chemistry=problem_chemistry, &
        thermo_data_file_name=thermo_file, &
        transport_data_file_name=transdata_file, &
        molar_masses_data_file_name='')

    problem_solver_options = solver_options_c( &
        solver_name='fds_low_mach', &
        hydrodynamics_flag=.true., &
        heat_transfer_flag=.true., &
        molecular_diffusion_flag=molecular_diffusion_flag, &
        soret_diffusion_flag=soret_diffusion_flag, &
        viscosity_flag=viscosity_flag, &
        thermal_radiation_flag=thermal_radiation_flag, &
        chemical_reaction_flag=.true., &
        grav_acc=(/0.0_dp,0.0_dp,0.0_dp/), &
        additional_particles_phases=0, &
        CFL_flag=.true., &
        CFL_coefficient=cfl_coefficient, &
        initial_time_step=initial_time_step, &
        chemistry_backend=trim(chemistry_backend), &
        chemistry_cvode_relative_tolerance=chemistry_cvode_relative_tolerance, &
        chemistry_cvode_absolute_tolerance=chemistry_cvode_absolute_tolerance, &
        chemistry_cvode_max_steps=chemistry_cvode_max_steps, &
        chemistry_qss1_relative_change_limit=chemistry_qss1_relative_change_limit, &
        chemistry_qss1_minimum_internal_step=chemistry_qss1_minimum_internal_step, &
        chemistry_qss1_active_concentration_fraction= &
            chemistry_qss1_active_concentration_fraction, &
        chemistry_qss1_step_growth_factor=chemistry_qss1_step_growth_factor, &
        chemistry_qss1_max_steps=chemistry_qss1_max_steps, &
        chemistry_qss2_error_tolerance=chemistry_qss2_error_tolerance, &
        chemistry_qss2_minimum_internal_step=chemistry_qss2_minimum_internal_step, &
        chemistry_qss2_active_concentration_fraction= &
            chemistry_qss2_active_concentration_fraction, &
        chemistry_qss2_max_steps=chemistry_qss2_max_steps)

    ! The configured inlet_velocity is the initial guess.  In anchor mode the
    ! shared flame_stabilization_solver subsequently updates all inlet boundary
    ! types to keep the multidimensional flame envelope inside the domain.
    problem_flame_control = flame_stabilization_control_c( &
        mode=trim(flame_stabilization_mode), &
        scientific_result_title=trim(scientific_result_title), &
        case_setup=trim(flame_case_setup), &
        time_delay=flame_time_delay, &
        time_track=flame_time_track, &
        time_control=flame_time_control, &
        response_settle_time=flame_response_settle_time, &
        measurement_max_duration=flame_measurement_max_duration, &
        measurement_min_displacement_cells= &
            flame_measurement_min_displacement_cells)
    problem_controls_setup = problem_controls_c( &
        flame_stabilization=problem_flame_control)

    problem_mpi_support = mpi_communications_c(problem_domain)
    problem_data_manager = data_manager_c( &
        problem_domain, problem_mpi_support, problem_chemistry, &
        problem_thermophysics, problem_solver_options, problem_controls_setup)

    call problem_data_manager%create_boundary_conditions( &
        problem_boundaries, number_of_boundary_types=4, default_boundary=1)
    call problem_data_manager%create_computational_mesh(problem_mesh)
    call problem_data_manager%create_scalar_field(T,'temperature','T')
    call problem_data_manager%create_scalar_field(rho,'density','rho')
    call problem_data_manager%create_scalar_field(p,'pressure','p')
    call problem_data_manager%create_vector_field(v,'velocity','v','spatial')
    call problem_data_manager%create_vector_field( &
        Y,'specie_mass_fraction','Y','chemical')

    mesh_cell_size = problem_mesh%get_cell_edges_length()
    utter_loop = problem_domain%get_global_utter_cells_bounds()
    species_number = problem_chemistry%species_number

    i_h2  = problem_chemistry%get_chemical_specie_index('H2')
    i_o2  = problem_chemistry%get_chemical_specie_index('O2')
    i_n2  = problem_chemistry%get_chemical_specie_index('N2')
    i_h2o = problem_chemistry%get_chemical_specie_index('H2O')

    if (min(i_h2,i_o2,i_n2,i_h2o) <= 0) then
        error stop 'KEROMNES mechanism is missing required H2/O2/N2/H2O species'
    end if

    ! Fresh composition on a one-mole-H2 basis.
    x_h2 = hydrogen_mole_percent/100.0_dp
    nu = (1.0_dp-x_h2)/x_h2/(1.0_dp+n2_o2_molar_ratio)
    fresh_h2 = 1.0_dp
    fresh_o2 = nu
    fresh_n2 = nu*n2_o2_molar_ratio

    ! Complete-combustion product approximation that conserves H/O/N atoms.
    ! For the requested lean 15% H2 case this retains the correct excess O2.
    prod_h2  = max(1.0_dp-2.0_dp*nu,0.0_dp)
    prod_h2o = min(1.0_dp,2.0_dp*nu)
    prod_o2  = max(nu-0.5_dp,0.0_dp)
    prod_n2  = fresh_n2

    p%cells(:,:,:) = operating_pressure
    T%cells(:,:,:) = fresh_temperature
    v%pr(1)%cells(:,:,:) = inlet_velocity
    v%pr(2)%cells(:,:,:) = 0.0_dp
    if (size(v%pr) >= 3) v%pr(3)%cells(:,:,:) = 0.0_dp

    do i = 1, species_number
        Y%pr(i)%cells(:,:,:) = 0.0_dp
    end do

    call seed_random_generator(random_seed_value)
    pi_local = acos(-1.0_dp)
    base_front_index = nint(front_position_fraction*real(nx,dp))

    do j = utter_loop(2,1), utter_loop(2,2)
        random_value = 0.5_dp
        if (random_jitter_amplitude > 0.0_dp) call random_number(random_value)

        front_shift = perturbation_amplitude*cos( &
            2.0_dp*pi_local*real(perturbation_mode,dp)* &
            (real(j,dp)-0.5_dp)/real(ny,dp) + perturbation_phase)
        front_shift = front_shift + &
            random_jitter_amplitude*(2.0_dp*random_value-1.0_dp)

        front_index = base_front_index + nint(front_shift/cell_size_xy)

        do i = utter_loop(1,1), utter_loop(1,2)
            if (i <= front_index) then
                T%cells(i,j,1) = fresh_temperature
                Y%pr(i_h2)%cells(i,j,1) = fresh_h2
                Y%pr(i_o2)%cells(i,j,1) = fresh_o2
                Y%pr(i_n2)%cells(i,j,1) = fresh_n2
                Y%pr(i_h2o)%cells(i,j,1) = 0.0_dp
            else
                T%cells(i,j,1) = burned_temperature
                Y%pr(i_h2)%cells(i,j,1) = prod_h2
                Y%pr(i_o2)%cells(i,j,1) = prod_o2
                Y%pr(i_n2)%cells(i,j,1) = prod_n2
                Y%pr(i_h2o)%cells(i,j,1) = prod_h2o
            end if
        end do
    end do

    ! Convert the mole-number ratios above to normalized mass fractions.
    call problem_thermophysics%change_field_units_mole_to_dimless(Y)

    ! Make the generated restart state itself periodic before it is written.
    ! The calls are harmless for non-periodic directions.
    if (transverse_periodic) then
        call problem_mpi_support%exchange_conservative_scalar_field(p)
        call problem_mpi_support%exchange_conservative_scalar_field(T)
        call problem_mpi_support%exchange_conservative_vector_field(v)
        call problem_mpi_support%exchange_conservative_vector_field(Y)
    end if

    ! Inlet at x-min and outlet at x-max.  The transverse y boundaries are
    ! either periodic topology or the original adiabatic slip walls.
    call problem_boundaries%create_boundary_type( &
        type_name='wall', slip=.true., conductive=.false., &
        wall_temperature=0.0_dp, wall_conductivity_ratio=0.0_dp, priority=1)

    call problem_boundaries%create_boundary_type( &
        type_name='outlet', &
        farfield_pressure=operating_pressure, &
        farfield_temperature=burned_temperature, &
        farfield_velocity=inlet_velocity, &
        farfield_species_names=(/'H2 ','O2 ','N2 ','H2O'/), &
        farfield_concentrations=(/prod_h2,prod_o2,prod_n2,prod_h2o/), &
        priority=2)

    call problem_boundaries%create_boundary_type( &
        type_name='inlet', &
        farfield_pressure=operating_pressure, &
        farfield_temperature=fresh_temperature, &
        farfield_velocity=inlet_velocity, &
        farfield_species_names=(/'H2','O2','N2'/), &
        farfield_concentrations=(/fresh_h2,fresh_o2,fresh_n2/), &
        priority=3)

    call problem_boundaries%create_boundary_type( &
        type_name='symmetry_plane', priority=4)

    problem_boundaries%bc_markers(utter_loop(1,1),:,:) = 3
    problem_boundaries%bc_markers(utter_loop(1,2),:,:) = 2
    if (.not. transverse_periodic) then
        problem_boundaries%bc_markers(:,utter_loop(2,1),:) = 1
        problem_boundaries%bc_markers(:,utter_loop(2,2),:) = 1
    end if
    call problem_mpi_support%exchange_boundary_conditions_markers(problem_boundaries)

    ! Production diagnostics.
    transducer_offset = max(1,nint(0.001_dp/mesh_cell_size(1)))
    observation_slice = utter_loop
    observation_slice(2,:) = (utter_loop(2,1)+utter_loop(2,2))/2
    summation_region = utter_loop

    problem_post_proc_manager = post_processor_manager_c( &
        problem_data_manager, number_post_processors=1)

    call problem_post_proc_manager%create_post_processor( &
        problem_data_manager, post_processor_name='proc1', &
        operations_number=8, save_time=postprocess_interval_us, &
        save_time_units='microseconds')

    call problem_post_proc_manager%create_post_processor_operation( &
        problem_data_manager,1,'energy_production_chemistry','max', &
        operation_area=observation_slice,grad_projection=1)
    call problem_post_proc_manager%create_post_processor_operation( &
        problem_data_manager,1,'specie_production_chemistry(H2)','sum', &
        operation_area=summation_region,operation_area_distance=(/0,0,0/))
    call problem_post_proc_manager%create_post_processor_operation( &
        problem_data_manager,1,'specie_mass_fraction(H2)','transducer', &
        operation_area_distance=(/transducer_offset,0,0/))
    call problem_post_proc_manager%create_post_processor_operation( &
        problem_data_manager,1,'specie_mass_fraction(H2)','transducer', &
        operation_area_distance=(/-transducer_offset,0,0/))
    call problem_post_proc_manager%create_post_processor_operation( &
        problem_data_manager,1,'density','transducer', &
        operation_area_distance=(/transducer_offset,0,0/))
    call problem_post_proc_manager%create_post_processor_operation( &
        problem_data_manager,1,'density','transducer', &
        operation_area_distance=(/-transducer_offset,0,0/))
    call problem_post_proc_manager%create_post_processor_operation( &
        problem_data_manager,1,'temperature','transducer', &
        operation_area_distance=(/transducer_offset,0,0/))
    call problem_post_proc_manager%create_post_processor_operation( &
        problem_data_manager,1,'temperature','transducer', &
        operation_area_distance=(/-transducer_offset,0,0/))

    problem_data_save = data_save_c( &
        problem_data_manager, &
        visible_fields_names=(/'pressure                       ', &
            'pressure_dynamic               ', &
            'temperature                    ', &
            'density                        ', &
            'velocity                       ', &
            'specie_mass_fraction           ', &
            'velocity_of_sound              ', &
            'velocity_production_viscosity  ', &
            'specie_production_chemistry    ', &
            'energy_production_chemistry    ', &
            'specie_production_diffusion    ', &
            'diffusivity(H2)                '/), &
        save_time=field_save_interval_us, &
        save_time_units='microseconds', &
        save_format='tecplot', &
        data_save_folder='data_save', &
        debug_flag=.false.)

    problem_data_io = data_io_c( &
        problem_data_manager, &
        check_time=restart_check_interval_ms, &
        check_time_units='milliseconds', &
        output_time=wall_time_output_limit_min, &
        data_output_folder='data_output')

    call write_case_summary(log_unit, nx, ny, nu)
    close(log_unit)

    call problem_data_io%output_all_data(0.0_dp,stop_flag,make_output=.true.)
    call problem_data_save%save_all_data(0.0_dp,stop_flag,make_save=.true.)

    call change_directory(trim(initial_work_dir))

    write(*,'(A)') 'Case generated successfully: '//trim(work_dir)
    write(*,'(A,A)') 'Chemistry backend: ',trim(chemistry_backend)
    if (trim(to_lower(chemistry_backend)) == 'qss2') then
        write(*,'(A,ES12.4)') 'QSS2 tolerance: ',chemistry_qss2_error_tolerance
    end if

contains

    subroutine validate_configuration()
        character(len=24) :: backend

        if (len_trim(case_id) == 0) error stop 'case_id cannot be empty'
        if (len_trim(results_root) == 0) error stop 'results_root cannot be empty'

        if (hydrogen_mole_percent <= 0.0_dp .or. &
            hydrogen_mole_percent >= 100.0_dp) then
            error stop 'hydrogen_mole_percent must be between 0 and 100'
        end if
        if (n2_o2_molar_ratio <= 0.0_dp) error stop 'Invalid N2/O2 ratio'
        if (operating_pressure <= 0.0_dp) error stop 'Pressure must be positive'
        if (fresh_temperature <= 0.0_dp .or. burned_temperature <= 0.0_dp) &
            error stop 'Temperatures must be positive'
        if (inlet_velocity < 0.0_dp) error stop 'inlet_velocity cannot be negative'

        if (domain_length <= 0.0_dp .or. domain_width <= 0.0_dp .or. &
            cell_size_xy <= 0.0_dp) error stop 'Invalid domain/mesh dimensions'
        if (front_position_fraction <= 0.0_dp .or. &
            front_position_fraction >= 1.0_dp) &
            error stop 'front_position_fraction must be inside (0,1)'
        if (perturbation_mode < 0) error stop 'perturbation_mode cannot be negative'
        if (perturbation_amplitude < 0.0_dp .or. &
            random_jitter_amplitude < 0.0_dp) &
            error stop 'Perturbation amplitudes cannot be negative'

        if (cfl_coefficient <= 0.0_dp .or. initial_time_step <= 0.0_dp) &
            error stop 'Invalid time-step controls'

        select case (trim(to_lower(flame_stabilization_mode)))
        case ('anchor','none')
        case ('laminar_burning_velocity')
            error stop 'Use anchor, not laminar_burning_velocity, for this 2D case'
        case default
            error stop 'flame_stabilization_mode must be anchor/none'
        end select
        if (flame_time_delay < 0.0_dp .or. flame_time_track <= 0.0_dp .or. &
            flame_time_control <= 0.0_dp .or. &
            flame_response_settle_time < 0.0_dp) then
            error stop 'Invalid flame stabilization timing controls'
        end if

        backend = trim(to_lower(chemistry_backend))
        select case (trim(backend))
        case ('slatec','cvode','qss1','qss2')
        case default
            error stop 'chemistry_backend must be slatec/cvode/qss1/qss2'
        end select

        if (chemistry_qss2_error_tolerance <= 0.0_dp) &
            error stop 'QSS2 tolerance must be positive'
        if (chemistry_qss2_minimum_internal_step <= 0.0_dp) &
            error stop 'QSS2 minimum step must be positive'
        if (chemistry_qss2_active_concentration_fraction <= 0.0_dp) &
            error stop 'QSS2 active fraction must be positive'
        if (chemistry_qss2_max_steps <= 0) &
            error stop 'QSS2 max steps must be positive'

        if (postprocess_interval_us <= 0.0_dp .or. &
            field_save_interval_us <= 0.0_dp .or. &
            restart_check_interval_ms <= 0.0_dp .or. &
            wall_time_output_limit_min <= 0.0_dp) &
            error stop 'Output intervals must be positive'
    end subroutine validate_configuration


    subroutine namelist_error(group_name, ios)
        character(len=*), intent(in) :: group_name
        integer, intent(in) :: ios
        write(*,'(A,A,A,I0)') 'ERROR reading /',trim(group_name), &
            '/ namelist, iostat=',ios
        error stop 1
    end subroutine namelist_error


    subroutine seed_random_generator(base_seed)
        integer, intent(in) :: base_seed
        integer :: seed_size, idx
        integer, allocatable :: seed_values(:)

        call random_seed(size=seed_size)
        allocate(seed_values(seed_size))
        do idx = 1, seed_size
            seed_values(idx) = abs(base_seed + 104729*(idx-1))
            if (seed_values(idx) == 0) seed_values(idx) = idx
        end do
        call random_seed(put=seed_values)
        deallocate(seed_values)
    end subroutine seed_random_generator


    character(len=32) function to_lower(text)
        character(len=*), intent(in) :: text
        integer :: idx, code
        to_lower = ''
        do idx = 1, min(len_trim(text),len(to_lower))
            code = iachar(text(idx:idx))
            if (code >= iachar('A') .and. code <= iachar('Z')) then
                to_lower(idx:idx) = achar(code+32)
            else
                to_lower(idx:idx) = text(idx:idx)
            end if
        end do
    end function to_lower


    subroutine write_case_summary(unit, nx_local, ny_local, nu_local)
        integer, intent(in) :: unit, nx_local, ny_local
        real(dp), intent(in) :: nu_local

        write(unit,'(A)') '============================================================'
        write(unit,'(A)') '2D inlet-stabilized H2-air flame production case'
        write(unit,'(A,A)') 'case_id: ',trim(case_id)
        write(unit,'(A,A)') 'mechanism: ',trim(mech_name)
        write(unit,'(A,A)') 'solver: fds_low_mach'
        write(unit,'(A,A)') 'chemistry backend: ',trim(chemistry_backend)
        write(unit,'(A,F10.4)') 'H2 mole percent: ',hydrogen_mole_percent
        write(unit,'(A,ES14.6)') 'operating pressure [Pa]: ',operating_pressure
        write(unit,'(A,F10.4)') 'fresh temperature [K]: ',fresh_temperature
        write(unit,'(A,F10.4)') 'burned seed temperature [K]: ',burned_temperature
        write(unit,'(A,F10.6)') 'fixed inlet velocity [m/s]: ',inlet_velocity
        write(unit,'(A,I0,A,I0)') 'mesh: ',nx_local,' x ',ny_local
        write(unit,'(A,ES14.6)') 'dx [m]: ',cell_size_xy
        write(unit,'(A,F10.5)') 'front position fraction: ',front_position_fraction
        write(unit,'(A,ES14.6)') 'sinusoidal amplitude [m]: ',perturbation_amplitude
        write(unit,'(A,I0)') 'sinusoidal mode: ',perturbation_mode
        write(unit,'(A,ES14.6)') 'sinusoidal phase [rad]: ',perturbation_phase
        write(unit,'(A,ES14.6)') 'random jitter amplitude [m]: ', &
            random_jitter_amplitude
        write(unit,'(A,I0)') 'random seed: ',random_seed_value
        write(unit,'(A,L1)') 'transverse y periodic: ',transverse_periodic
        write(unit,'(A,ES14.6)') 'fresh O2/H2 mole ratio: ',nu_local
        write(unit,'(A,ES14.6)') 'product residual H2 amount: ',prod_h2
        write(unit,'(A,ES14.6)') 'product residual O2 amount: ',prod_o2
        write(unit,'(A,ES14.6)') 'product H2O amount: ',prod_h2o
        write(unit,'(A,L1)') 'molecular diffusion: ',molecular_diffusion_flag
        write(unit,'(A,L1)') 'Soret diffusion: ',soret_diffusion_flag
        write(unit,'(A,A)') 'flame stabilization mode: ', &
            trim(flame_stabilization_mode)
        write(unit,'(A,ES14.6)') 'flame time delay [s]: ',flame_time_delay
        write(unit,'(A,ES14.6)') 'flame track interval [s]: ',flame_time_track
        write(unit,'(A,ES14.6)') 'flame control interval [s]: ', &
            flame_time_control
        write(unit,'(A,L1)') 'viscosity: ',viscosity_flag
        write(unit,'(A,L1)') 'thermal radiation: ',thermal_radiation_flag
        write(unit,'(A,F10.5)') 'CFL coefficient: ',cfl_coefficient
        write(unit,'(A,ES14.6)') 'initial dt [s]: ',initial_time_step
        if (trim(to_lower(chemistry_backend)) == 'qss2') then
            write(unit,'(A,ES14.6)') 'QSS2 tolerance: ', &
                chemistry_qss2_error_tolerance
            write(unit,'(A,ES14.6)') 'QSS2 minimum step [s]: ', &
                chemistry_qss2_minimum_internal_step
            write(unit,'(A,ES14.6)') 'QSS2 active fraction: ', &
                chemistry_qss2_active_concentration_fraction
            write(unit,'(A,I0)') 'QSS2 max steps: ',chemistry_qss2_max_steps
        end if
        write(unit,'(A)') 'inlet velocity is an initial guess; anchor control is adaptive'
        write(unit,'(A)') '============================================================'
    end subroutine write_case_summary


    subroutine get_current_directory(directory)
        character(len=*), intent(out) :: directory
        character(kind=c_char) :: buffer(4096)
        type(c_ptr) :: ptr
        integer :: idx, max_copy

        buffer = c_null_char
#ifdef WIN
        ptr = c_getcwd(buffer, int(size(buffer), c_int))
#else
        ptr = c_getcwd(buffer, int(size(buffer), c_size_t))
#endif
        if (.not. c_associated(ptr)) then
            error stop 'Unable to obtain current working directory.'
        end if

        directory = ''
        max_copy = min(len(directory), size(buffer))
        do idx = 1, max_copy
            if (buffer(idx) == c_null_char) exit
            directory(idx:idx) = transfer(buffer(idx), directory(idx:idx))
        end do

        if (idx > len(directory) .and. &
            buffer(min(idx,size(buffer))) /= c_null_char) then
            error stop 'Current working directory path exceeds internal buffer.'
        end if
    end subroutine get_current_directory


    subroutine change_directory(directory)
        character(len=*), intent(in) :: directory
        character(kind=c_char), allocatable :: c_path(:)
        integer(c_int) :: status
        integer :: idx, path_length

        path_length = len_trim(directory)
        allocate(c_path(path_length+1))
        do idx = 1, path_length
            c_path(idx) = transfer(directory(idx:idx), c_path(idx))
        end do
        c_path(path_length+1) = c_null_char

        status = c_chdir(c_path)
        if (status /= 0_c_int) then
            write(*,'(A)') 'Unable to change working directory to: '// &
                trim(directory)
            error stop 'Directory change failed.'
        end if
    end subroutine change_directory


    subroutine run_filesystem_command(command)
        character(len=*), intent(in) :: command
        integer :: cmdstat, exitstat

        cmdstat = 0
        exitstat = 0
        call execute_command_line(command, wait=.true., &
            exitstat=exitstat, cmdstat=cmdstat)
        if (cmdstat /= 0 .or. exitstat /= 0) then
            write(*,'(A)') 'Filesystem command failed: '//trim(command)
            write(*,'(A,I0,A,I0)') 'cmdstat=',cmdstat,', exitstat=',exitstat
            error stop 'Case-directory setup failed.'
        end if
    end subroutine run_filesystem_command


    subroutine ensure_directory(directory)
        character(len=*), intent(in) :: directory
        character(len=:), allocatable :: command

        command = 'cmake -E make_directory "'//trim(directory)//'"'
        call run_filesystem_command(command)
    end subroutine ensure_directory


    subroutine copy_directory_tree(source_directory, destination_directory)
        character(len=*), intent(in) :: source_directory
        character(len=*), intent(in) :: destination_directory
        character(len=:), allocatable :: command

        command = 'cmake -E copy_directory "'//trim(source_directory)//'" "'// &
            trim(destination_directory)//'"'
        call run_filesystem_command(command)
    end subroutine copy_directory_tree

end program package_interface
