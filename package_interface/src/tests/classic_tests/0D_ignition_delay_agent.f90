!===============================================================================
! 0D CONSTANT-VOLUME IGNITION-DELAY INTERFACE FOR AGENTIC WORKFLOWS
!===============================================================================
!
! A spatially uniform, closed, adiabatic reactor is represented by a very small
! one-dimensional Cartesian domain. Hydrodynamics and chemistry are retained,
! while viscosity, molecular diffusion, heat conduction and radiation are
! disabled. Uniform initial conditions and two adiabatic slip walls make the
! solution equivalent to a homogeneous constant-volume reactor, subject to the
! usual numerical round-off and spatial-uniformity checks.
!
! The executable reads setup_input.nml, validates it, creates one case directory,
! copies task_setup, writes the NRG input files and configures transducer-based
! time histories of temperature, pressure and all species mass fractions.
!
! Platform notes:
!   * No USE IFPORT dependency.
!   * Standard EXECUTE_COMMAND_LINE is used for mkdir/copy operations.
!   * WIN is defined by CMake on Windows; POSIX commands are used otherwise.
!   * A small ISO_C_BINDING wrapper changes the current directory portably.
!
!===============================================================================

program package_interface

    use, intrinsic :: iso_c_binding, only : c_char, c_int, c_null_char
    use kind_parameters
    use global_data
    use computational_domain_class
    use chemical_properties_class
    use thermophysical_properties_class
    use solver_options_class
    use computational_mesh_class
    use mpi_communications_class
    use data_manager_class
    use boundary_conditions_class
    use field_scalar_class
    use field_vector_class
    use data_save_class
    use data_io_class
    use post_processor_manager_class
    use supplementary_routines

    implicit none

#ifdef WIN
    interface
        function c_chdir(path) bind(C, name='_chdir') result(status)
            import :: c_char, c_int
            character(kind=c_char), intent(in) :: path(*)
            integer(c_int) :: status
        end function c_chdir
    end interface
#else
    interface
        function c_chdir(path) bind(C, name='chdir') result(status)
            import :: c_char, c_int
            character(kind=c_char), intent(in) :: path(*)
            integer(c_int) :: status
        end function c_chdir
    end interface
#endif

    !--------------------------------------------------------------------------
    ! NRG objects
    !--------------------------------------------------------------------------
    type(computational_domain)                :: problem_domain
    type(data_manager)                        :: problem_data_manager
    type(mpi_communications)                  :: problem_mpi_support
    type(chemical_properties),       target   :: problem_chemistry
    type(thermophysical_properties), target   :: problem_thermophysics
    type(solver_options)                      :: problem_solver_options
    type(computational_mesh),        target   :: problem_mesh
    type(boundary_conditions),       target   :: problem_boundaries
    type(field_scalar_cons),         target   :: p, T, rho
    type(field_vector_cons),         target   :: v, Y
    type(post_processor_manager)              :: problem_post_proc_manager
    type(data_io)                              :: problem_data_io
    type(data_save)                            :: problem_data_save

    !--------------------------------------------------------------------------
    ! Agent-configurable namelist parameters
    !--------------------------------------------------------------------------
    character(len=80)  :: experiment_id
    character(len=32)  :: mechanism_id
    character(len=32)  :: solver_id
    real(dp)           :: initial_temperature
    real(dp)           :: initial_pressure
    real(dp)           :: hydrogen_mole_percent
    real(dp)           :: initial_time_step
    real(dp)           :: final_time
    real(dp)           :: postprocess_interval
    real(dp)           :: field_save_interval
    real(dp)           :: checkpoint_interval
    integer            :: cells_number_x
    real(dp)           :: cell_length_x
    logical            :: save_spatial_fields

    namelist /reactor_config/ experiment_id, mechanism_id, solver_id, &
        initial_temperature, initial_pressure, hydrogen_mole_percent, &
        initial_time_step, final_time, postprocess_interval,          &
        field_save_interval, checkpoint_interval, cells_number_x,    &
        cell_length_x, save_spatial_fields

    !--------------------------------------------------------------------------
    ! Derived configuration and local variables
    !--------------------------------------------------------------------------
    character(len=64)  :: mechanism_name
    character(len=64)  :: mech_file
    character(len=64)  :: thermo_file
    character(len=64)  :: transdata_file
    character(len=32)  :: solver_name
    character(len=512) :: work_dir
    character(len=1024):: command
    character(len=512) :: iomsg

    integer :: io_unit, log_unit, ierr, command_status, exit_status
    integer :: spec
    integer :: utter_loop(3,2)
    integer :: transducer_point(3)
    real(dp):: domain_length
    real(dp):: nu
    real(dp):: phi
    logical :: stop_flag

    !--------------------------------------------------------------------------
    ! Defaults: a single inexpensive, constant-volume H2-air ignition case
    !--------------------------------------------------------------------------
    experiment_id          = 'h2_air_cv_T1500_P1atm_XH2_20'
    mechanism_id           = 'keromnes'
    solver_id              = 'cpm'
    initial_temperature    = 1500.0_dp
    initial_pressure       = 101325.0_dp
    hydrogen_mole_percent  = 20.0_dp
    initial_time_step      = 1.0e-8_dp
    final_time             = 30.0e-6_dp
    postprocess_interval   = 0.1_dp       ! microseconds
    field_save_interval    = 5.0_dp       ! microseconds
    checkpoint_interval    = 25.0_dp      ! microseconds
    cells_number_x         = 2
    cell_length_x          = 1.0e-4_dp
    save_spatial_fields    = .true.

    !--------------------------------------------------------------------------
    ! Read and validate the agent-authored namelist
    !--------------------------------------------------------------------------
    open(newunit=io_unit, file='setup_input.nml', status='old', &
         action='read', iostat=ierr, iomsg=iomsg)
    if (ierr /= 0) then
        write(*,'(A)') 'ERROR: cannot open setup_input.nml: '//trim(iomsg)
        error stop 1
    end if

    read(io_unit, nml=reactor_config, iostat=ierr, iomsg=iomsg)
    close(io_unit)
    if (ierr /= 0) then
        write(*,'(A)') 'ERROR: cannot read /reactor_config/: '//trim(iomsg)
        error stop 2
    end if

    call normalize_and_validate_configuration()
    call resolve_mechanism()
    call resolve_solver()

    ! H2 + nu O2 + 3.762 nu N2, with the requested H2 mole percentage.
    nu  = (100.0_dp - hydrogen_mole_percent) / &
          hydrogen_mole_percent / 4.762_dp
    phi = 4.762_dp * 0.5_dp * hydrogen_mole_percent / 100.0_dp / &
          (1.0_dp - hydrogen_mole_percent / 100.0_dp)

    domain_length = real(cells_number_x, dp) * cell_length_x

    !--------------------------------------------------------------------------
    ! Create one isolated case directory and copy task_setup portably
    !--------------------------------------------------------------------------
    work_dir = '0D_ignition_delay_agent'//trim(fold_sep)// &
               trim(experiment_id)

#ifdef WIN
    command = 'if not exist "'//trim(work_dir)//'" mkdir "'// &
              trim(work_dir)//'"'
#else
    command = 'mkdir -p "'//trim(work_dir)//'"'
#endif
    call run_command_checked(command, 'create case directory')

#ifdef WIN
    command = 'xcopy "'//trim(task_setup_folder)//'" "'// &
              trim(work_dir)//trim(fold_sep)//trim(task_setup_folder)// &
              '" /E /I /K /Y >NUL'
#else
    command = 'cp -a "'//trim(task_setup_folder)//'" "'// &
              trim(work_dir)//trim(fold_sep)//'"'
#endif
    call run_command_checked(command, 'copy task_setup directory')

    call change_directory(trim(work_dir))

    ! Preserve the exact experiment input beside generated NRG files.
#ifdef WIN
    command = 'copy /Y "..'//trim(fold_sep)//'..'//trim(fold_sep)// &
              'setup_input.nml" "." >NUL'
#else
    command = 'cp "../../setup_input.nml" "."'
#endif
    call run_command_checked(command, 'copy setup_input.nml into case')

    !--------------------------------------------------------------------------
    ! Computational domain: minimal uniform spatial surrogate for a 0D reactor
    !--------------------------------------------------------------------------
    problem_domain = computational_domain_c(                         &
        dimensions        = 1,                                       &
        cells_number      = (/cells_number_x, 1, 1/),                &
        coordinate_system = 'cartesian',                             &
        lengths           = reshape((/0.0_dp, 0.0_dp, 0.0_dp,       &
                                      domain_length,                  &
                                      cell_length_x, cell_length_x/), &
                                    (/3,2/)),                         &
        axis_names        = (/'x','y','z'/))

    !--------------------------------------------------------------------------
    ! Chemistry, thermodynamics and solver
    !--------------------------------------------------------------------------
    problem_chemistry = chemical_properties_c(                       &
        chemical_mechanism_file_name  = trim(mech_file),             &
        default_enhanced_efficiencies = 1.0_dp,                      &
        E_act_units                   = 'cal.mol')

    problem_thermophysics = thermophysical_properties_c(             &
        chemistry                   = problem_chemistry,              &
        thermo_data_file_name       = trim(thermo_file),              &
        transport_data_file_name    = trim(transdata_file),           &
        molar_masses_data_file_name = 'molar_masses.dat')

    problem_solver_options = solver_options_c(                       &
        solver_name                  = trim(solver_name),             &
        hydrodynamics_flag           = .true.,                        &
        heat_transfer_flag           = .true.,                        &
        molecular_diffusion_flag     = .false.,                       &
        viscosity_flag               = .false.,                       &
        thermal_radiation_flag       = .false.,                       &
        chemical_reaction_flag       = .true.,                        &
        grav_acc                     = (/0.0_dp,0.0_dp,0.0_dp/),      &
        additional_particles_phases = 0,                             &
        CFL_flag                     = .false.,                       &
        CFL_coefficient              = 0.25_dp,                       &
        initial_time_step            = initial_time_step)

    problem_mpi_support  = mpi_communications_c(problem_domain)
    problem_data_manager = data_manager_c(problem_domain,            &
                                           problem_mpi_support,       &
                                           problem_chemistry,         &
                                           problem_thermophysics,     &
                                           problem_solver_options)

    call problem_data_manager%create_boundary_conditions(            &
        problem_boundaries, number_of_boundary_types=1,              &
        default_boundary=1)
    call problem_data_manager%create_computational_mesh(problem_mesh)

    call problem_data_manager%create_scalar_field(p,   'pressure',    'p')
    call problem_data_manager%create_scalar_field(T,   'temperature', 'T')
    call problem_data_manager%create_scalar_field(rho, 'density',     'rho')
    call problem_data_manager%create_vector_field(v, 'velocity',      &
                                                   'v', 'spatial')
    call problem_data_manager%create_vector_field(Y,                  &
        'specie_mass_fraction', 'Y', 'chemical')

    utter_loop = problem_domain%get_global_utter_cells_bounds()

    !--------------------------------------------------------------------------
    ! Uniform initial conditions
    !--------------------------------------------------------------------------
    T%cells(:,:,:) = initial_temperature
    p%cells(:,:,:) = initial_pressure
    v%pr(1)%cells(:,:,:) = 0.0_dp

    do spec = 1, problem_chemistry%species_number
        Y%pr(spec)%cells(:,:,:) = 0.0_dp
    end do
    Y%pr(1)%cells(:,:,:) = 1.0_dp
    Y%pr(2)%cells(:,:,:) = nu
    Y%pr(3)%cells(:,:,:) = nu * 3.762_dp

    ! Input composition above is molar. Convert it once to normalized mass
    ! fractions; the field name now correctly describes the stored quantity.
    call problem_thermophysics%change_field_units_mole_to_dimless(Y)

    ! Closed, adiabatic, slip walls on both ends => fixed volume and no fluxes.
    call problem_boundaries%create_boundary_type(                     &
        type_name               = 'wall',                             &
        slip                    = .true.,                             &
        conductive              = .false.,                            &
        wall_temperature        = 0.0_dp,                             &
        wall_conductivity_ratio = 0.0_dp,                             &
        priority                = 1)

    problem_boundaries%bc_markers(utter_loop(1,1),:,:) = 1
    problem_boundaries%bc_markers(utter_loop(1,2),:,:) = 1

    !--------------------------------------------------------------------------
    ! Time histories: one central transducer for T, p and all species Y_k.
    ! Vector-field transducer output is expected to contain all chemical
    ! components in mechanism order.
    !--------------------------------------------------------------------------
    transducer_point = (/0, 0, 0/)

    problem_post_proc_manager = post_processor_manager_c(             &
        problem_data_manager, number_post_processors=1)

    call problem_post_proc_manager%create_post_processor(              &
        problem_data_manager,                                          &
        post_processor_name = 'reactor_history',                        &
        operations_number   = 3,                                       &
        save_time           = postprocess_interval,                    &
        save_time_units     = 'microseconds')

    call problem_post_proc_manager%create_post_processor_operation(    &
        problem_data_manager, 1, 'temperature', 'transducer',          &
        operation_area_distance=transducer_point)
    call problem_post_proc_manager%create_post_processor_operation(    &
        problem_data_manager, 1, 'pressure', 'transducer',             &
        operation_area_distance=transducer_point)
    call problem_post_proc_manager%create_post_processor_operation(    &
        problem_data_manager, 1, 'specie_mass_fraction', 'transducer', &
        operation_area_distance=transducer_point)

    ! Spatial snapshots remain useful for verifying that the nominally 0D
    ! calculation stays spatially uniform. They may be disabled in the namelist.
    if (save_spatial_fields) then
        problem_data_save = data_save_c(                               &
            problem_data_manager,                                     &
            visible_fields_names=[character(len=40) ::                 &
                'pressure', 'temperature', 'density', 'velocity',      &
                'specie_mass_fraction'],                               &
            save_time        = field_save_interval,                    &
            save_time_units  = 'microseconds',                         &
            save_format      = 'tecplot',                              &
            data_save_folder = 'data_save',                            &
            debug_flag       = .false.)
    else
        problem_data_save = data_save_c(                               &
            problem_data_manager,                                     &
            visible_fields_names=[character(len=40) ::                 &
                'pressure', 'temperature', 'specie_mass_fraction'],    &
            save_time        = final_time * 1.0e6_dp,                  &
            save_time_units  = 'microseconds',                         &
            save_format      = 'tecplot',                              &
            data_save_folder = 'data_save',                            &
            debug_flag       = .false.)
    end if

    problem_data_io = data_io_c(                                      &
        problem_data_manager,                                         &
        check_time         = checkpoint_interval,                     &
        check_time_units   = 'microseconds',                          &
        data_output_folder = 'data_output')

    !--------------------------------------------------------------------------
    ! Correct problem log
    !--------------------------------------------------------------------------
    open(newunit=log_unit, file=problem_setup_log_file,                &
         status='replace', form='formatted', action='write')
    write(log_unit,'(A)') 'General description:'
    write(log_unit,'(A)') '  Homogeneous hydrogen-air ignition in a closed vessel.'
    write(log_unit,'(A)') 'Main aim:'
    write(log_unit,'(A)') '  Verify NRG constant-volume chemistry/thermodynamics against a reference solution.'
    write(log_unit,'(A)') 'Physical formulation:'
    write(log_unit,'(A)') '  Uniform Cartesian domain with closed adiabatic slip walls.'
    write(log_unit,'(A)') '  Transport, viscosity, radiation and gravity are disabled.'
    write(log_unit,'(A)') '  Fixed geometry gives an adiabatic constant-volume reactor.'
    write(log_unit,'(A,I0)') 'Computational cells: ', cells_number_x
    write(log_unit,'(A,ES16.8)') 'Cell length [m]: ', cell_length_x
    write(log_unit,'(A,ES16.8)') 'Initial temperature [K]: ', initial_temperature
    write(log_unit,'(A,ES16.8)') 'Initial pressure [Pa]: ', initial_pressure
    write(log_unit,'(A,F12.6)')  'Initial H2 mole percent: ', hydrogen_mole_percent
    write(log_unit,'(A,F12.6)')  'Equivalence ratio: ', phi
    write(log_unit,'(A,A)')      'Chemical mechanism: ', trim(mechanism_name)
    write(log_unit,'(A,A)')      'Solver: ', trim(solver_name)
    write(log_unit,'(A,ES16.8)') 'Initial time step [s]: ', initial_time_step
    write(log_unit,'(A,ES16.8)') 'Requested run time [s] (used by the laboratory runner): ', final_time
    write(log_unit,'(A)') 'Post-processing:'
    write(log_unit,'(A)') '  Central transducers: temperature, pressure and all species mass fractions.'
    write(log_unit,'(A)') 'Validation requirement:'
    write(log_unit,'(A)') '  Confirm spatial uniformity from saved fields before treating the case as 0D.'
    write(log_unit,'(A)') '----------------------------------------------------------'
    close(log_unit)

    ! Initial files used by computing_module.
    call problem_data_io%output_all_data(0.0_dp, stop_flag, make_output=.true.)
    call problem_data_save%save_all_data(0.0_dp, stop_flag, make_save=.true.)

    write(*,'(A)') 'Case generated successfully: '//trim(work_dir)

contains

    subroutine normalize_and_validate_configuration()
        character(len=:), allocatable :: id_trim

        experiment_id = adjustl(experiment_id)
        mechanism_id  = to_lower_ascii(adjustl(mechanism_id))
        solver_id     = to_lower_ascii(adjustl(solver_id))
        id_trim       = trim(experiment_id)

        if (len_trim(id_trim) == 0) then
            error stop 'ERROR: experiment_id must not be empty.'
        end if
        if (index(id_trim, '..') /= 0 .or. index(id_trim, '/') /= 0 .or. &
            index(id_trim, achar(92)) /= 0) then
            error stop 'ERROR: experiment_id must not contain path components.'
        end if
        if (initial_temperature <= 0.0_dp) then
            error stop 'ERROR: initial_temperature must be positive.'
        end if
        if (initial_pressure <= 0.0_dp) then
            error stop 'ERROR: initial_pressure must be positive.'
        end if
        if (hydrogen_mole_percent <= 0.0_dp .or. &
            hydrogen_mole_percent >= 100.0_dp) then
            error stop 'ERROR: hydrogen_mole_percent must lie strictly between 0 and 100.'
        end if
        if (initial_time_step <= 0.0_dp .or. final_time <= 0.0_dp) then
            error stop 'ERROR: time-step and final-time values must be positive.'
        end if
        if (postprocess_interval <= 0.0_dp .or. &
            field_save_interval <= 0.0_dp .or. &
            checkpoint_interval <= 0.0_dp) then
            error stop 'ERROR: all output intervals must be positive.'
        end if
        if (cells_number_x /= 2 .and. cells_number_x /= 4) then
            error stop 'ERROR: cells_number_x must be 2 or 4 for the MVP reactor.'
        end if
        if (cell_length_x <= 0.0_dp) then
            error stop 'ERROR: cell_length_x must be positive.'
        end if
    end subroutine normalize_and_validate_configuration


    subroutine resolve_mechanism()
        select case(trim(mechanism_id))
        case('agafonov')
            mechanism_name = 'TEREZA'
        case('keromnes')
            mechanism_name = 'KEROMNES'
        case('konnov')
            mechanism_name = 'KONNOV'
        case default
            write(*,'(A,A)') 'ERROR: unsupported mechanism_id: ', trim(mechanism_id)
            error stop 3
        end select

        mech_file      = trim(mechanism_name)//'.txt'
        thermo_file    = trim(mechanism_name)//'_THERMO.txt'
        transdata_file = trim(mechanism_name)//'_TRANSDATA.txt'
    end subroutine resolve_mechanism


    subroutine resolve_solver()
        select case(trim(solver_id))
        case('cpm')
            solver_name = 'cpm'
        case('cabaret')
            solver_name = 'CABARET'
        case('fds', 'fds_low_mach')
            solver_name = 'fds_low_mach'
        case default
            write(*,'(A,A)') 'ERROR: unsupported solver_id: ', trim(solver_id)
            error stop 4
        end select
    end subroutine resolve_solver


    subroutine run_command_checked(cmd, description)
        character(len=*), intent(in) :: cmd, description

        command_status = 0
        exit_status    = 0
        call execute_command_line(trim(cmd), wait=.true.,             &
                                  cmdstat=command_status,              &
                                  exitstat=exit_status, cmdmsg=iomsg)
        if (command_status /= 0 .or. exit_status /= 0) then
            write(*,'(A)') 'ERROR: failed to '//trim(description)//'.'
            write(*,'(A)') 'Command: '//trim(cmd)
            if (len_trim(iomsg) > 0) write(*,'(A)') trim(iomsg)
            write(*,'(A,I0,A,I0)') 'cmdstat=', command_status, &
                                    ', exitstat=', exit_status
            error stop 5
        end if
    end subroutine run_command_checked


    subroutine change_directory(path)
        character(len=*), intent(in) :: path
        character(kind=c_char), allocatable :: c_path(:)
        integer(c_int) :: status
        integer :: i, n

        n = len_trim(path)
        allocate(c_path(0:n))
        do i = 1, n
            c_path(i-1) = path(i:i)
        end do
        c_path(n) = c_null_char

        status = c_chdir(c_path)
        deallocate(c_path)

        if (status /= 0_c_int) then
            write(*,'(A)') 'ERROR: cannot change directory to '//trim(path)
            error stop 6
        end if
    end subroutine change_directory


    pure function to_lower_ascii(value) result(lower_value)
        character(len=*), intent(in) :: value
        character(len=len(value)) :: lower_value
        integer :: i, code

        lower_value = value
        do i = 1, len(value)
            code = iachar(lower_value(i:i))
            if (code >= iachar('A') .and. code <= iachar('Z')) then
                lower_value(i:i) = achar(code + iachar('a') - iachar('A'))
            end if
        end do
    end function to_lower_ascii

end program package_interface
