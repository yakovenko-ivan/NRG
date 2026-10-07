!===============================================================================
! 2-D periodic-advection verification interface for the FDS low-Mach solver.
!
! Purpose:
!   Exercise y-periodic topology with a smooth composition wave that is
!   advected through the periodic seam by a uniform transverse velocity.
!
! Exact inviscid/nonreactive reference:
!   Y_w(y,t) = Yw0 + A*sin(2*pi*m*(y-v_y*t)/Ly + phase)
!   Y_c(y,t) = Yc0 - A*sin(2*pi*m*(y-v_y*t)/Ly + phase)
!
! The default NML performs exactly one traversal of the y-periodic domain.
! Physical boundaries exist only at x-min/x-max and are adiabatic slip walls.
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
    use run_control_class
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
    type(computational_domain)              :: problem_domain
    type(data_manager)                      :: problem_data_manager
    type(mpi_communications)                :: problem_mpi_support
    type(chemical_properties), target       :: problem_chemistry
    type(thermophysical_properties), target :: problem_thermophysics
    type(solver_options)                    :: problem_solver_options
    type(flame_stabilization_control)       :: problem_flame_control
    type(problem_controls)                  :: problem_controls_setup
    type(run_control)                       :: problem_run_control
    type(computational_mesh), target        :: problem_mesh
    type(boundary_conditions), target       :: problem_boundaries
    type(field_scalar_cons), target         :: p, T, rho
    type(field_vector_cons), target         :: v, Y
    type(post_processor_manager)            :: problem_post_proc_manager
    type(data_io)                           :: problem_data_io
    type(data_save)                         :: problem_data_save

    ! Case identity.
    character(len=64)   :: case_id
    character(len=512)  :: results_root
    character(len=1024) :: work_dir
    character(len=4096) :: initial_work_dir
    character(len=4096) :: config_file
    logical             :: overwrite_existing_case

    namelist /case_config/ case_id, results_root, overwrite_existing_case

    ! Reference thermodynamic/compositional state.
    character(len=32) :: mechanism_id
    character(len=16) :: wave_species
    character(len=16) :: carrier_species
    real(dp) :: operating_pressure
    real(dp) :: temperature
    real(dp) :: wave_mean_mass_fraction
    real(dp) :: carrier_mean_mass_fraction

    namelist /mixture_config/ mechanism_id, operating_pressure, temperature, &
        wave_species, carrier_species, wave_mean_mass_fraction, &
        carrier_mean_mass_fraction

    ! Domain topology.
    real(dp) :: domain_length_x
    real(dp) :: domain_length_y
    integer  :: cells_x
    integer  :: cells_y
    logical  :: periodic_x
    logical  :: periodic_y

    namelist /domain_config/ domain_length_x, domain_length_y, cells_x, cells_y, &
        periodic_x, periodic_y

    ! Uniform transport velocity.
    real(dp) :: velocity_x
    real(dp) :: velocity_y

    namelist /initial_velocity_config/ velocity_x, velocity_y

    ! Smooth periodic perturbation of the two-species mixture.
    real(dp) :: wave_amplitude
    integer  :: wave_mode
    real(dp) :: wave_phase

    namelist /species_wave_config/ wave_amplitude, wave_mode, wave_phase

    ! These switches are deliberately exposed but the verification interface
    ! requires every source/transport model to remain disabled.
    logical :: chemical_reaction_flag
    logical :: molecular_diffusion_flag
    logical :: soret_diffusion_flag
    logical :: viscosity_flag
    logical :: heat_transfer_flag
    logical :: thermal_radiation_flag

    namelist /physics_config/ chemical_reaction_flag, molecular_diffusion_flag, &
        soret_diffusion_flag, viscosity_flag, heat_transfer_flag, &
        thermal_radiation_flag

    ! Deterministic time integration.
    logical  :: cfl_flag
    real(dp) :: cfl_coefficient
    real(dp) :: initial_time_step
    real(dp) :: final_time

    namelist /numerics_config/ cfl_flag, cfl_coefficient, initial_time_step, &
        final_time

    ! Output controls.
    real(dp) :: field_save_interval_us
    real(dp) :: restart_check_interval_ms
    real(dp) :: wall_time_output_limit_min

    namelist /output_config/ field_save_interval_us, restart_check_interval_ms, &
        wall_time_output_limit_min

    ! Derived data.
    character(len=32) :: mech_name
    character(len=64) :: mech_file, thermo_file, transdata_file
    real(dp) :: dx, dy, y_center, phase_value
    real(dp) :: wave_fraction, carrier_fraction
    real(dp) :: pi_local, traversal_period
    real(dp) :: periods_requested
    integer  :: i, j, spec
    integer  :: i_wave, i_carrier
    integer  :: species_number
    integer  :: log_unit, io_unit, ierr
    integer  :: ref_unit
    integer, dimension(3,2) :: utter_loop
    logical :: stop_flag, config_exists, case_exists

    ! Defaults match the companion NML.
    case_id = 'periodic_y_species_wave'
    results_root = '2D_periodic_boundary_tests'
    overwrite_existing_case = .true.

    mechanism_id = 'keromnes'
    operating_pressure = 101325.0_dp
    temperature = 300.0_dp
    wave_species = 'O2'
    carrier_species = 'N2'
    wave_mean_mass_fraction = 0.20_dp
    carrier_mean_mass_fraction = 0.80_dp

    domain_length_x = 1.6e-3_dp
    domain_length_y = 6.4e-3_dp
    cells_x = 16
    cells_y = 64
    periodic_x = .false.
    periodic_y = .true.

    velocity_x = 0.0_dp
    velocity_y = 1.0_dp

    wave_amplitude = 0.05_dp
    wave_mode = 1
    wave_phase = 0.0_dp

    chemical_reaction_flag = .false.
    molecular_diffusion_flag = .false.
    soret_diffusion_flag = .false.
    viscosity_flag = .false.
    heat_transfer_flag = .false.
    thermal_radiation_flag = .false.

    cfl_flag = .false.
    cfl_coefficient = 0.25_dp
    initial_time_step = 1.0e-5_dp
    final_time = 6.4e-3_dp

    field_save_interval_us = 200.0_dp
    restart_check_interval_ms = 1000.0_dp
    wall_time_output_limit_min = 240.0_dp

    call get_current_directory(initial_work_dir)

    config_file = '2D_periodic_advection.nml'
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
    if (ierr /= 0) error stop 'Unable to open periodic-advection namelist'

    rewind(io_unit)
    read(io_unit, nml=case_config, iostat=ierr)
    if (ierr /= 0) call namelist_error('case_config', ierr)

    rewind(io_unit)
    read(io_unit, nml=mixture_config, iostat=ierr)
    if (ierr /= 0) call namelist_error('mixture_config', ierr)

    rewind(io_unit)
    read(io_unit, nml=domain_config, iostat=ierr)
    if (ierr /= 0) call namelist_error('domain_config', ierr)

    rewind(io_unit)
    read(io_unit, nml=initial_velocity_config, iostat=ierr)
    if (ierr /= 0) call namelist_error('initial_velocity_config', ierr)

    rewind(io_unit)
    read(io_unit, nml=species_wave_config, iostat=ierr)
    if (ierr /= 0) call namelist_error('species_wave_config', ierr)

    rewind(io_unit)
    read(io_unit, nml=physics_config, iostat=ierr)
    if (ierr /= 0) call namelist_error('physics_config', ierr)

    rewind(io_unit)
    read(io_unit, nml=numerics_config, iostat=ierr)
    if (ierr /= 0) call namelist_error('numerics_config', ierr)

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
        error stop 'Periodic-advection interface currently supports mechanism_id=keromnes'
    end select

    dx = domain_length_x/real(cells_x,dp)
    dy = domain_length_y/real(cells_y,dp)
    traversal_period = domain_length_y/abs(velocity_y)
    periods_requested = final_time/traversal_period

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
    call write_nrg_source_revision(log_unit, '2D periodic advection verification')

    problem_domain = computational_domain_c( &
        dimensions=2, &
        cells_number=(/cells_x,cells_y,1/), &
        coordinate_system='cartesian', &
        lengths=reshape((/0.0_dp,0.0_dp,0.0_dp, &
            domain_length_x,domain_length_y,1.0e-3_dp/),(/3,2/)), &
        axis_names=(/'x','y','z'/), &
        periodic=(/periodic_x,periodic_y,.false./))

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
        heat_transfer_flag=heat_transfer_flag, &
        molecular_diffusion_flag=molecular_diffusion_flag, &
        soret_diffusion_flag=soret_diffusion_flag, &
        viscosity_flag=viscosity_flag, &
        thermal_radiation_flag=thermal_radiation_flag, &
        chemical_reaction_flag=chemical_reaction_flag, &
        grav_acc=(/0.0_dp,0.0_dp,0.0_dp/), &
        additional_particles_phases=0, &
        CFL_flag=cfl_flag, &
        CFL_coefficient=cfl_coefficient, &
        initial_time_step=initial_time_step)

    problem_flame_control = flame_stabilization_control_c( &
        mode='none', &
        scientific_result_title='2D periodic species-wave advection verification', &
        case_setup='2D_periodic_advection')
    problem_controls_setup = problem_controls_c( &
        flame_stabilization=problem_flame_control)

    ! Writes run_control.dat for the computing module.
    problem_run_control = run_control_c( &
        termination_mode='simulation_time', final_time=final_time)

    problem_mpi_support = mpi_communications_c(problem_domain)
    problem_data_manager = data_manager_c( &
        problem_domain, problem_mpi_support, problem_chemistry, &
        problem_thermophysics, problem_solver_options, problem_controls_setup)

    call problem_data_manager%create_boundary_conditions( &
        problem_boundaries, number_of_boundary_types=1, default_boundary=1)
    call problem_data_manager%create_computational_mesh(problem_mesh)
    call problem_data_manager%create_scalar_field(T,'temperature','T')
    call problem_data_manager%create_scalar_field(rho,'density','rho')
    call problem_data_manager%create_scalar_field(p,'pressure','p')
    call problem_data_manager%create_vector_field(v,'velocity','v','spatial')
    call problem_data_manager%create_vector_field( &
        Y,'specie_mass_fraction','Y','chemical')

    utter_loop = problem_domain%get_global_utter_cells_bounds()
    species_number = problem_chemistry%species_number

    i_wave = problem_chemistry%get_chemical_specie_index(trim(wave_species))
    i_carrier = problem_chemistry%get_chemical_specie_index(trim(carrier_species))
    if (i_wave <= 0 .or. i_carrier <= 0) then
        error stop 'Configured wave/carrier species are missing from the mechanism'
    end if
    if (i_wave == i_carrier) then
        error stop 'wave_species and carrier_species must be different'
    end if

    p%cells(:,:,:) = operating_pressure
    T%cells(:,:,:) = temperature
    rho%cells(:,:,:) = 0.0_dp

    v%pr(1)%cells(:,:,:) = velocity_x
    v%pr(2)%cells(:,:,:) = velocity_y
    if (size(v%pr) >= 3) v%pr(3)%cells(:,:,:) = 0.0_dp

    do spec = 1, species_number
        Y%pr(spec)%cells(:,:,:) = 0.0_dp
    end do

    pi_local = acos(-1.0_dp)
    do j = utter_loop(2,1), utter_loop(2,2)
        y_center = (real(j,dp)-0.5_dp)*dy
        phase_value = 2.0_dp*pi_local*real(wave_mode,dp)* &
            y_center/domain_length_y + wave_phase
        wave_fraction = wave_mean_mass_fraction + wave_amplitude*sin(phase_value)
        carrier_fraction = carrier_mean_mass_fraction - &
            wave_amplitude*sin(phase_value)

        do i = utter_loop(1,1), utter_loop(1,2)
            Y%pr(i_wave)%cells(i,j,1) = wave_fraction
            Y%pr(i_carrier)%cells(i,j,1) = carrier_fraction
        end do
    end do

    ! Physical topology:
    !   x-min/x-max -> adiabatic slip wall;
    !   y-min/y-max -> periodic topology, no physical boundary marker.
    call problem_boundaries%create_boundary_type( &
        type_name='wall', slip=.true., conductive=.false., &
        wall_temperature=0.0_dp, wall_conductivity_ratio=0.0_dp, priority=1)

    ! Make the intended marker topology explicit even if an older boundary
    ! constructor initialized all outer ghost planes with default_boundary.
    problem_boundaries%bc_markers(:,utter_loop(2,1),:) = 0
    problem_boundaries%bc_markers(:,utter_loop(2,2),:) = 0
    problem_boundaries%bc_markers(utter_loop(1,1),:,:) = 1
    problem_boundaries%bc_markers(utter_loop(1,2),:,:) = 1

    ! No online reduction is needed: the saved field is compared with the
    ! exact translated wave after the run.
    problem_post_proc_manager = post_processor_manager_c( &
        problem_data_manager, number_post_processors=0)

    problem_data_save = data_save_c( &
        problem_data_manager, &
        visible_fields_names=(/'pressure                       ', &
            'pressure_dynamic               ', &
            'temperature                    ', &
            'density                        ', &
            'velocity                       ', &
            'specie_mass_fraction           '/), &
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

    call write_case_summary(log_unit, traversal_period, periods_requested)
    close(log_unit)

    call write_reference_profile('periodic_reference_initial.dat')

    stop_flag = .false.
    call problem_data_io%output_all_data(0.0_dp, stop_flag, make_output=.true.)
    call problem_data_save%save_all_data(0.0_dp, stop_flag, make_save=.true.)

    call change_directory(trim(initial_work_dir))

    write(*,'(A)') 'Case generated successfully: '//trim(work_dir)
    write(*,'(A,ES12.4,A)') 'Periodic traversal time: ', traversal_period, ' s'
    write(*,'(A,F10.4)') 'Requested traversals: ', periods_requested

contains

    subroutine validate_configuration()
        real(dp) :: fraction_margin

        if (len_trim(case_id) == 0) error stop 'case_id cannot be empty'
        if (len_trim(results_root) == 0) error stop 'results_root cannot be empty'

        if (operating_pressure <= 0.0_dp) error stop 'Pressure must be positive'
        if (temperature <= 0.0_dp) error stop 'Temperature must be positive'
        if (len_trim(wave_species) == 0 .or. len_trim(carrier_species) == 0) &
            error stop 'Both mixture species names are required'
        if (trim(wave_species) == trim(carrier_species)) &
            error stop 'wave_species and carrier_species must differ'

        if (wave_mean_mass_fraction < 0.0_dp .or. &
            carrier_mean_mass_fraction < 0.0_dp) &
            error stop 'Mean mass fractions cannot be negative'
        if (abs(wave_mean_mass_fraction + carrier_mean_mass_fraction - 1.0_dp) > &
            1.0e-12_dp) then
            error stop 'The two mean mass fractions must sum to one'
        end if

        fraction_margin = min(wave_mean_mass_fraction, carrier_mean_mass_fraction)
        if (wave_amplitude < 0.0_dp .or. wave_amplitude > fraction_margin) &
            error stop 'wave_amplitude would create an invalid mass fraction'
        if (wave_mode <= 0) error stop 'wave_mode must be positive'

        if (domain_length_x <= 0.0_dp .or. domain_length_y <= 0.0_dp) &
            error stop 'Domain lengths must be positive'
        if (cells_x < 2 .or. cells_y < 4) &
            error stop 'Periodic-advection mesh is too small'

        if (periodic_x) &
            error stop 'This verification case requires periodic_x=.false.'
        if (.not. periodic_y) &
            error stop 'This verification case requires periodic_y=.true.'

        if (abs(velocity_y) <= tiny(1.0_dp)) &
            error stop 'velocity_y must be non-zero to cross the periodic seam'

        if (chemical_reaction_flag .or. molecular_diffusion_flag .or. &
            soret_diffusion_flag .or. viscosity_flag .or. heat_transfer_flag .or. &
            thermal_radiation_flag) then
            error stop 'P2 periodic-advection verification requires all source/transport models OFF'
        end if

        if (initial_time_step <= 0.0_dp .or. final_time <= 0.0_dp) &
            error stop 'Time controls must be positive'
        if (cfl_coefficient <= 0.0_dp) error stop 'CFL coefficient must be positive'

        if (field_save_interval_us <= 0.0_dp .or. &
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


    subroutine write_case_summary(unit, period, requested_periods)
        integer, intent(in) :: unit
        real(dp), intent(in) :: period, requested_periods

        write(unit,'(A)') '============================================================'
        write(unit,'(A)') '2D periodic species-wave advection verification'
        write(unit,'(A,A)') 'case_id: ', trim(case_id)
        write(unit,'(A,A)') 'mechanism: ', trim(mech_name)
        write(unit,'(A,A)') 'solver: fds_low_mach'
        write(unit,'(A,ES14.6)') 'operating pressure [Pa]: ', operating_pressure
        write(unit,'(A,ES14.6)') 'temperature [K]: ', temperature
        write(unit,'(A,A)') 'wave species: ', trim(wave_species)
        write(unit,'(A,A)') 'carrier species: ', trim(carrier_species)
        write(unit,'(A,ES14.6)') 'wave mean mass fraction: ', &
            wave_mean_mass_fraction
        write(unit,'(A,ES14.6)') 'carrier mean mass fraction: ', &
            carrier_mean_mass_fraction
        write(unit,'(A,ES14.6)') 'wave amplitude: ', wave_amplitude
        write(unit,'(A,I0)') 'wave mode: ', wave_mode
        write(unit,'(A,ES14.6)') 'wave phase [rad]: ', wave_phase
        write(unit,'(A,I0,A,I0)') 'mesh: ', cells_x, ' x ', cells_y
        write(unit,'(A,ES14.6)') 'dx [m]: ', dx
        write(unit,'(A,ES14.6)') 'dy [m]: ', dy
        write(unit,'(A,L1)') 'periodic x: ', periodic_x
        write(unit,'(A,L1)') 'periodic y: ', periodic_y
        write(unit,'(A,ES14.6)') 'velocity x [m/s]: ', velocity_x
        write(unit,'(A,ES14.6)') 'velocity y [m/s]: ', velocity_y
        write(unit,'(A,ES14.6)') 'traversal period [s]: ', period
        write(unit,'(A,ES14.6)') 'final time [s]: ', final_time
        write(unit,'(A,F12.6)') 'requested traversals: ', requested_periods
        write(unit,'(A,L1)') 'CFL adaptation: ', cfl_flag
        write(unit,'(A,ES14.6)') 'initial/fixed dt [s]: ', initial_time_step
        write(unit,'(A)') 'x boundaries: adiabatic slip walls'
        write(unit,'(A)') 'y boundaries: periodic topology (bc_marker=0)'
        write(unit,'(A)') 'chemistry/diffusion/Soret/viscosity/heat/radiation: OFF'
        write(unit,'(A)') '============================================================'
    end subroutine write_case_summary


    subroutine write_reference_profile(filename)
        character(len=*), intent(in) :: filename
        real(dp) :: y_ref, arg_ref, yw_ref, yc_ref
        integer :: jj

        open(newunit=ref_unit, file=trim(filename), status='replace', &
            form='formatted')
        write(ref_unit,'(A)') '# j  y[m]  Y_wave(t=0)  Y_carrier(t=0)'
        do jj = 1, cells_y
            y_ref = (real(jj,dp)-0.5_dp)*dy
            arg_ref = 2.0_dp*pi_local*real(wave_mode,dp)* &
                y_ref/domain_length_y + wave_phase
            yw_ref = wave_mean_mass_fraction + wave_amplitude*sin(arg_ref)
            yc_ref = carrier_mean_mass_fraction - wave_amplitude*sin(arg_ref)
            write(ref_unit,'(I6,1X,ES20.12,1X,ES20.12,1X,ES20.12)') &
                jj, y_ref, yw_ref, yc_ref
        end do
        close(ref_unit)
    end subroutine write_reference_profile


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
