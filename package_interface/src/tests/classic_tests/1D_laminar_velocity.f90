!================================================================================
! 1D LAMINAR VELOCITY SIMULATION INTERFACE
!================================================================================
!
! PROGRAM: package_interface
!
! DESCRIPTION:
!   This program serves as the main interface for 1D laminar burning velocity
!   simulations. It sets up and executes parametric studies for different:
!   - Physical setups (counter-flow, near-wall)
!   - Coordinate systems (Cartesian, cylindrical, spherical)
!   - Solver types (FDS, CPM, CABARET)
!   - Chemical mechanisms (KEROMNES)
!   - Hydrogen concentrations (0-100%)
!   - Spatial resolutions (dx from 6.25e-06 to 4.0e-04)
!
! KEY FEATURES:
!   - Automated parameter sweep through nested loops
!   - Directory structure creation for organized output
!   - Multiple initial condition configurations
!   - Comprehensive post-processing setup
!   - Support for flamelet initialization from precomputed tables
!
! PHYSICAL SYSTEM:
!   Simulates 1D reactive flow with hydrogen-air mixtures
!   Primary species: H2, O2, N2, H2O
!   Temperature range: 300K - 2500K
!   Pressure: Atmospheric (101325 Pa)
!
!================================================================================

program package_interface

    !==========================================
    ! MODULE IMPORTS
    !==========================================
    use ifport                     ! Intel Fortran portability routines
    use kind_parameters            ! Defines precision kinds (dp, sp, etc.)
    use global_data                ! Global constants and parameters
    use computational_domain_class ! Domain definition and management
    use chemical_properties_class  ! Chemical kinetics and species data
    use thermophysical_properties_class ! Thermodynamic and transport properties
    use solver_options_class       ! Numerical solver configuration
    use computational_mesh_class   ! Mesh generation and management
    use mpi_communications_class   ! Parallel communication routines
    use data_manager_class         ! Central data management
    use boundary_conditions_class  ! Boundary condition specification
    use field_scalar_class         ! Scalar field operations (p, T, rho, etc.)
    use field_vector_class         ! Vector field operations (v, Y, etc.)
    use data_save_class            ! Data output and saving
    use data_io_class              ! Input/output operations
    use post_processor_manager_class ! Post-processing and monitoring
    use supplementary_routines     ! Utility functions and helper routines

    implicit none
    
    !==========================================
    ! PRIMARY SIMULATION OBJECTS
    !==========================================
    type(computational_domain)              :: problem_domain        ! Computational domain
    type(data_manager)                      :: problem_data_manager  ! Central data coordinator
    type(mpi_communications)                :: problem_mpi_support   ! MPI communication handler
    
    ! Physical properties objects (with TARGET attribute for pointer associations)
    type(chemical_properties)        ,target :: problem_chemistry      ! Chemical reaction data
    type(thermophysical_properties)  ,target :: problem_thermophysics  ! Thermodynamic properties
    
    ! Solver configuration
    type(solver_options)                     :: problem_solver_options ! Numerical method settings
    
    ! Solution fields (with TARGET attribute)
    type(computational_mesh)         ,target :: problem_mesh          ! Computational grid
    type(boundary_conditions)        ,target :: problem_boundaries    ! Boundary conditions
    type(field_scalar_cons)          ,target :: p, T, rho, E_f, v_s   ! Scalar fields:
                                                                      ! p - pressure
                                                                      ! T - temperature
                                                                      ! rho - density
                                                                      ! E_f - flame energy?
                                                                      ! v_s - speed of sound
    type(field_vector_cons)          ,target :: v, Y                  ! Vector fields:
                                                                      ! v - velocity (3 components)
                                                                      ! Y - species mass fractions
    
    ! Post-processing and I/O
    type(post_processor_manager)             :: problem_post_proc_manager  ! Monitoring and analysis
    type(data_io)                            :: problem_data_io        ! Data input/output
    type(data_save)                          :: problem_data_save      ! Solution file writing
    
    !==========================================
    ! GEOMETRIC AND DOMAIN PARAMETERS
    !==========================================
    real(dp)    ,dimension(3)   :: cell_size           ! Cell dimensions (dx, dy, dz)
    integer     ,dimension(3,2) :: utter_loop          ! Global computational bounds
    integer     ,dimension(3,2) :: observation_slice   ! Region for monitoring
    integer     ,dimension(3,2) :: summation_region    ! Region for integral calculations
    integer                     :: transducer_offset   ! Offset for transducer placement
    
    !==========================================
    ! CHEMICAL AND PHYSICAL PARAMETERS
    !==========================================
    integer                     :: species_number      ! Number of chemical species
    real(dp)                    :: domain_length       ! Domain length in x-direction [m]
    real(dp)                    :: domain_width        ! Domain width [m]
    real(dp)                    :: ignition_region     ! Ignition zone location/size
    real(dp)                    :: CFL_coeff           ! CFL stability coefficient
    real(dp)                    :: delta_x             ! Spatial resolution [m]
    real(dp)                    :: offset              ! General offset parameter
    real(dp)                    :: nu                  ! Stoichiometric O2/H2 ratio
    real(dp)                    :: X_H2                ! Hydrogen mole fraction [%]
    real(dp)                    :: spec_summ           ! Temporary sum for normalization
    
    !==========================================
    ! FLAMELET TABLE DATA (for counter_flow_precInc setup)
    !==========================================
    integer                     :: table_size          ! Size of flamelet lookup table
    integer                     :: scaling_factor      ! Table scaling factor
    integer                     :: n                   ! Generic counter
    real(dp), allocatable, dimension(:)   :: p_t, T_t, rho_t, vx_t, E_f_t, gamma_t, vs_t
    real(dp), allocatable, dimension(:,:) :: Y_t       ! Flamelet species data
    real(dp), allocatable, dimension(:)   :: temp      ! Temporary storage for file reading
    
    !==========================================
    ! FILE AND DIRECTORY PATHS
    !==========================================
    character(len=10)           :: string              ! Temporary string buffer
    character(len=500)          :: initial_work_dir    ! Initial working directory
    character(len=200)          :: work_dir            ! Current working directory
    character(len=200)          :: initials_file       ! Flamelet initialization file
    character(len=20)           :: solver_name         ! Solver type identifier
    character(len=20)           :: coordinate_system   ! Coordinate system type
    character(len=30)           :: setup               ! Physical setup type
    character(len=30)           :: mech_file           ! Chemical mechanism file
    character(len=30)           :: thermo_file         ! Thermodynamic data file
    character(len=30)           :: transdata_file      ! Transport data file
    
    !==========================================
    ! CONTROL AND STATUS VARIABLES
    !==========================================
    logical :: stop_flag          ! Simulation termination flag
    integer :: i, j, spec         ! Loop indices
    integer :: task1, task2, task3, task4, task5, task6  ! Parametric study indices
    integer :: flamelet_pos       ! Position of flamelet in domain
    integer :: ierr               ! Error status indicator
    integer :: io_unit            ! File I/O unit
    integer :: log_unit           ! Log file unit
    
    !==========================================
    ! PROGRAM EXECUTION BEGINS
    !==========================================
    
    ! Get initial working directory for later return
    ierr = getcwd(initial_work_dir)
    
    !================================================================
    ! PARAMETRIC STUDY LOOPS
    !================================================================
    ! Six nested loops for comprehensive parameter space exploration:
    ! task1: Physical setup (1=counter_flow, 2=counter_flow_precInc, 3=near_wall)
    ! task2: Coordinate system (1=Cartesian, 2=cylindrical, 3=spherical)
    ! task3: Solver type (1=FDS, 2=CPM, 3=CABARET)
    ! task4: Chemical mechanism (currently only KEROMNES)
    ! task5: Hydrogen concentration (9=9% H2)
    ! task6: Spatial resolution (2=dx=1.0e-04 m)
    !================================================================
    
    do task1 = 3, 3          ! Currently fixed at near_wall (3)
    do task2 = 1, 1          ! Currently fixed at Cartesian (1)
    do task3 = 1, 1          ! Currently fixed at FDS solver (1)
    do task4 = 1, 1          ! Currently fixed at KEROMNES mechanism (1)
    do task5 = 9, 9          ! Currently fixed at 9% H2
    do task6 = 2, 2          ! Currently fixed at dx=1.0e-04 (2)
        
        ! Initialize working directory structure
        work_dir = '1D_LBV_test'  ! Main results directory
        
        ! Create directory tree based on parameter choices
        ierr = system('mkdir '// work_dir)
        
        !------------------------------------------------
        ! TASK1: PHYSICAL SETUP SELECTION
        !------------------------------------------------
        select case(task1)
            case(1)
                work_dir = trim(work_dir) // trim(fold_sep) //'cf'
                setup = 'counter_flow'
            case(2)
                work_dir = trim(work_dir) // trim(fold_sep) //'cf_prcInc'
                setup = 'counter_flow_precInc'
            case(3)
                work_dir = trim(work_dir) // trim(fold_sep) //'nw'
                setup = 'near_wall'
        end select
        
        ierr = system('mkdir '// work_dir)
        
        !------------------------------------------------
        ! TASK2: COORDINATE SYSTEM SELECTION
        !------------------------------------------------
        select case(task2)
            case(1)
                work_dir = trim(work_dir) // trim(fold_sep) //'cartesian'
                coordinate_system = 'cartesian'
            case(2)
                work_dir = trim(work_dir) // trim(fold_sep) //'cylindrical'
                coordinate_system = 'cylindrical'
            case(3)
                work_dir = trim(work_dir) // trim(fold_sep) //'spherical'
                coordinate_system = 'spherical'
        end select
        
        ierr = system('mkdir '// work_dir)
        
        !------------------------------------------------
        ! TASK3: SOLVER TYPE SELECTION
        !------------------------------------------------
        select case(task3)
            case(1)
                work_dir = trim(work_dir) // trim(fold_sep) // 'FDS'
                solver_name = 'fds_low_mach'
            case(2)
                work_dir = trim(work_dir) // trim(fold_sep) // 'CPM'
                solver_name = 'cpm'
            case(3)
                work_dir = trim(work_dir) // trim(fold_sep) // 'CABARET'
                solver_name = 'CABARET'
        end select
        
        ierr = system('mkdir '// work_dir)
        
        !------------------------------------------------
        ! TASK4: CHEMICAL MECHANISM SELECTION
        !------------------------------------------------
        select case(task4)
            case(1)
                work_dir = trim(work_dir) // trim(fold_sep) // 'KEROMNES'
                mech_file      = 'KEROMNES.txt'
                thermo_file    = 'KEROMNES_THERMO.txt'
                transdata_file = 'KEROMNES_TRANSDATA.txt'
        end select
        
        ierr = system('mkdir '// work_dir)
        
        !------------------------------------------------
        ! TASK5: HYDROGEN CONCENTRATION
        !------------------------------------------------
        X_H2     = task5 * 1.0_dp  ! Convert to percentage
        work_dir = trim(work_dir) // trim(fold_sep) // trim(str_r(X_H2)) //'_pcnt'
        nu       = (100.0 - X_H2) / X_H2 / 4.762_dp  ! Calculate stoichiometric ratio
        
        ierr = system('mkdir '// work_dir)
        
        !------------------------------------------------
        ! TASK6: SPATIAL RESOLUTION
        !------------------------------------------------
        select case(task6)
            case(0)
                work_dir = trim(work_dir) // trim(fold_sep) // 'dx_4.0e-04'
                delta_x  = 4.0e-04_dp
            case(1)
                work_dir = trim(work_dir) // trim(fold_sep) // 'dx_2.0e-04'
                delta_x  = 2.0e-04_dp
            case(2)
                work_dir = trim(work_dir) // trim(fold_sep) // 'dx_1.0e-04'
                delta_x  = 1.0e-04_dp
            case(3)
                work_dir = trim(work_dir) // trim(fold_sep) // 'dx_5.0e-05'
                delta_x  = 5.0e-05_dp
            case(4)
                work_dir = trim(work_dir) // trim(fold_sep) // 'dx_2.5e-05'
                delta_x  = 2.5e-05_dp
            case(5)
                work_dir = trim(work_dir) // trim(fold_sep) // 'dx_1.25e-05'
                delta_x  = 1.25e-05_dp
            case(6)
                work_dir = trim(work_dir) // trim(fold_sep) // 'dx_6.25e-06'
                delta_x  = 6.25e-06_dp
        end select
        
        ierr = system('mkdir '// work_dir)
        
        ! Copy setup files to working directory and change to it
        ierr = system('echo d | xcopy .\' // trim(task_setup_folder) // ' .\' // trim(work_dir) // trim(fold_sep) // trim(task_setup_folder) // ' /E/K' )
        ierr = chdir(work_dir)
        
        ! Open log file for this configuration
        open(newunit = log_unit, file = problem_setup_log_file, status = 'replace', form = 'formatted')
        
        !================================================================
        ! DOMAIN DEFINITION
        !================================================================
        domain_length = 0.0256_dp  ! Fixed domain length [m]
        
        problem_domain = computational_domain_c(  &
            dimensions         = 1,                                     &
            cells_number       = (/int(domain_length/delta_x),1,1/),    &
            coordinate_system  = coordinate_system,                     &
            lengths            = reshape((/0.0_dp,0.000_dp,0.0_dp,      &
                                          domain_length,0.0025_dp,0.0025_dp/),(/3,2/)), &
            axis_names         = (/'x','y','z'/) )
        
        !================================================================
        ! CHEMICAL AND THERMOPHYSICAL PROPERTIES INITIALIZATION
        !================================================================
        problem_chemistry = chemical_properties_c( &
            chemical_mechanism_file_name     = mech_file,        &
            default_enhanced_efficiencies    = 1.0_dp,           &
            E_act_units                      = 'cal.mol')
        
        problem_thermophysics = thermophysical_properties_c( &
            chemistry                    = problem_chemistry,   &
            thermo_data_file_name        = thermo_file,        &
            transport_data_file_name     = transdata_file,     &
            molar_masses_data_file_name  = 'molar_masses.dat')
        
        !================================================================
        ! SOLVER OPTIONS CONFIGURATION
        !================================================================
        problem_solver_options = solver_options_c( &
            solver_name                 = solver_name,                              &
            hydrodynamics_flag          = .true.,      ! Solve momentum equations
            heat_transfer_flag          = .true.,      ! Solve energy equation
            molecular_diffusion_flag    = .true.,      ! Include species diffusion
            viscosity_flag              = .true.,      ! Include viscous effects
            chemical_reaction_flag      = .true.,      ! Include chemical reactions
            grav_acc                    = (/0.0_dp, 0.0_dp, 0.0_dp/),  ! No gravity
            additional_particles_phases = 0,           ! No particle phases
            CFL_flag                    = .true.,      ! Use CFL condition
            CFL_coefficient             = 0.25_dp,     ! CFL safety factor
            initial_time_step           = 1e-06_dp)    ! Initial Δt [s]
        
        !================================================================
        ! MPI AND DATA MANAGEMENT SETUP
        !================================================================
        problem_mpi_support   = mpi_communications_c(problem_domain)
        problem_data_manager  = data_manager_c(problem_domain, problem_mpi_support, &
                                               problem_chemistry, problem_thermophysics, &
                                               problem_solver_options)
        
        ! Create boundary conditions (2 types for inlet/outlet or wall/outlet)
        call problem_data_manager%create_boundary_conditions( &
            problem_boundaries, number_of_boundary_types = 2, default_boundary = 1)
        
        ! Create computational mesh
        call problem_data_manager%create_computational_mesh(problem_mesh)
        
        ! Create scalar solution fields
        call problem_data_manager%create_scalar_field(p,   'pressure',    'p')
        call problem_data_manager%create_scalar_field(T,   'temperature', 'T')
        call problem_data_manager%create_scalar_field(rho, 'density',     'rho')
        
        ! Create vector solution fields
        call problem_data_manager%create_vector_field(v, 'velocity', 'v', 'spatial')
        call problem_data_manager%create_vector_field(Y, 'specie_molar_concentration', 'Y', 'chemical')
        
        ! Get geometric information
        cell_size = problem_mesh%get_cell_edges_length()
        utter_loop = problem_domain%get_global_utter_cells_bounds()
        
        ! Define monitoring regions
        transducer_offset = 0.005 / cell_size(1)  ! 5mm offset in cell units
        
        observation_slice(:,1) = (/1, 1, 1/)
        observation_slice(:,2) = (/int(domain_length/delta_x), 1, 1/)
        
        summation_region(:,1) = (/-transducer_offset, 1, 1/)
        summation_region(:,2) = (/transducer_offset, 1, 1/)
        
        !================================================================
        ! POST-PROCESSING SETUP
        !================================================================
        problem_post_proc_manager = post_processor_manager_c(problem_data_manager, number_post_processors = 1)
        
        call problem_post_proc_manager%create_post_processor( &
            problem_data_manager,                            &
            post_processor_name = "proc1",                   &
            operations_number   = 7,                         &
            save_time           = 1.0_dp,                    &
            save_time_units     = 'microseconds')
        
        ! Define post-processing operations:
        ! 1. Minimum temperature gradient in observation slice
        ! 2-3. Pressure at ±5mm transducers
        ! 4-5. Density at ±5mm transducers
        ! 6-7. Temperature at ±5mm transducers
        call problem_post_proc_manager%create_post_processor_operation( &
            problem_data_manager, 1, 'temperature', 'min_grad', &
            operation_area = observation_slice, grad_projection = 1)
        call problem_post_proc_manager%create_post_processor_operation( &
            problem_data_manager, 1, 'pressure', 'transducer', &
            operation_area_distance = (/transducer_offset, 0, 0/))
        call problem_post_proc_manager%create_post_processor_operation( &
            problem_data_manager, 1, 'pressure', 'transducer', &
            operation_area_distance = (/-transducer_offset, 0, 0/))
        call problem_post_proc_manager%create_post_processor_operation( &
            problem_data_manager, 1, 'density', 'transducer', &
            operation_area_distance = (/transducer_offset, 0, 0/))
        call problem_post_proc_manager%create_post_processor_operation( &
            problem_data_manager, 1, 'density', 'transducer', &
            operation_area_distance = (/-transducer_offset, 0, 0/))
        call problem_post_proc_manager%create_post_processor_operation( &
            problem_data_manager, 1, 'temperature', 'transducer', &
            operation_area_distance = (/transducer_offset, 0, 0/))
        call problem_post_proc_manager%create_post_processor_operation( &
            problem_data_manager, 1, 'temperature', 'transducer', &
            operation_area_distance = (/-transducer_offset, 0, 0/))
        
        !================================================================
        ! DATA SAVING CONFIGURATION
        !================================================================
        problem_data_save = data_save_c( &
            problem_data_manager, &
            visible_fields_names = [ character(len=40) :: &
                'pressure',                       &
                'pressure_dynamic',               &
                'temperature',                    &
                'density',                        &
                'velocity',                       &
                'specie_molar_concentration',     &
                'velocity_of_sound',              &
                'velocity_production_viscosity',  &
                'mixture_cp',                     &
                'thermal_conductivity',           &
                'viscosity',                      &
                'energy_production_chemistry',    &
                'energy_production_diffusion'     &
            ], &
            save_time         = 100.0_dp,        ! Save interval
            save_time_units   = 'microseconds',  ! Time units for saving
            save_format       = 'tecplot',       ! Output format
            data_save_folder  = 'data_save',     ! Output directory
            debug_flag        = .false.)         ! Debug mode off
        
        !================================================================
        ! DATA I/O CONFIGURATION
        !================================================================
        problem_data_io = data_io_c( &
            problem_data_manager,    &
            check_time         = 5.0_dp,         ! Checkpoint interval
            check_time_units   = 'milliseconds', ! Time units for checkpoints
            data_output_folder = 'data_output')  ! Checkpoint directory
        
        !================================================================
        ! INITIAL CONDITIONS SETUP
        !================================================================
        
        ! Set uniform ambient conditions
        T%cells(:,:,:)   = 300.0_dp          ! Ambient temperature [K]
        p%cells(:,:,:)   = 1.0_dp * 101325.0_dp  ! Atmospheric pressure [Pa]
        
        ! Set initial species mole fractions (unburned mixture)
        Y%pr(1)%cells(:,:,:) = 1.0_dp        ! H2
        Y%pr(2)%cells(:,:,:) = nu            ! O2 (stoichiometric)
        Y%pr(3)%cells(:,:,:) = nu * 3.762_dp ! N2 (air composition)
        
        ! Apply setup-specific initialization
        select case(setup)
            
            !------------------------------------------------
            ! CASE: COUNTER-FLOW FLAME (Manual ignition)
            !------------------------------------------------
            case('counter_flow')
                ! Create hot ignition zone at 3/4 domain position
                do i = utter_loop(1,1), utter_loop(1,2)
                    if (i > int(3.0_dp * domain_length / delta_x / 4.0_dp) - (0.005_dp / delta_x)) then
                        T%cells(i,:,:) = 2500.0_dp  ! Ignition temperature
                        
                        ! Set products beyond ignition zone
                        if (i > int(3.0_dp * domain_length / delta_x / 4.0_dp)) then
                            T%cells(i,:,:)        = 300.0_dp          ! Cool products
                            Y%pr(1)%cells(i,:,:)  = 0.0_dp            ! No H2
                            Y%pr(2)%cells(i,:,:)  = 0.0_dp            ! No O2
                            Y%pr(3)%cells(i,:,:)  = nu * 3.762_dp     ! N2 only
                            Y%pr(7)%cells(i,:,:)  = 1.0_dp            ! H2O (assuming species 7)
                        end if
                    end if
                end do
                
                ! Convert mole fractions to dimensionless form
                call problem_thermophysics%change_field_units_mole_to_dimless(Y)
            
            !------------------------------------------------
            ! CASE: COUNTER-FLOW FLAME (Precomputed flamelet)
            !------------------------------------------------
            case('counter_flow_precInc')
                ! Load precomputed flamelet solution from file
                initials_file = trim(task_setup_folder) // trim(fold_sep) // &
                                'table_initials' // trim(fold_sep) // &
                                'Tereza_FULL' // trim(fold_sep) // &
                                'H2-Air_flamelet_' // trim(str_r(X_H2)) // &
                                '_pcnt_dx_' // trim(str_e(2*delta_x)) // '.dat'
                
                open(newunit = io_unit, file = initials_file, status = 'old', form = 'formatted')
                read(io_unit,*) string  ! Skip header
                
                ! Determine table size
                table_size = 0
                do 
                    read(io_unit,*,iostat = ierr) string
                    if(ierr /= 0) exit
                    table_size = table_size + 1
                end do
                table_size = table_size - 1  ! Adjust for header/format
                
                ! Allocate storage for flamelet data
                species_number = problem_chemistry%species_number
                allocate(Y_t(species_number, 0:table_size-1), &
                         T_t(0:table_size-1), &
                         vx_t(0:table_size-1))
                allocate(temp(species_number+4))
                
                ! Read flamelet data
                rewind(io_unit)
                read(io_unit,*) string  ! Skip header again
                
                do i = 0, table_size-1
                    read(io_unit,*) temp
                    T_t(i)   = temp(2)       ! Temperature
                    vx_t(i)  = temp(4)       ! Velocity
                    do spec = 1, species_number
                        Y_t(spec,i) = temp(spec+4)  ! Species mass fractions
                    end do
                end do
                close(io_unit)
                
                ! Interpolate flamelet onto computational grid
                flamelet_pos = int(2.0 * domain_length / delta_x / 4.0)  ! Center position
                
                do i = utter_loop(1,1), utter_loop(1,2)
                    ! Within flamelet region: interpolate
                    if ((i >= flamelet_pos - (table_size - 1)) .and. &
                        (i <= flamelet_pos + (table_size - 1))) then
                        
                        if (mod(i - flamelet_pos, 2) == 0) then 
                            ! Even indices: direct assignment
                            T%cells(i,1,1) = T_t((i - flamelet_pos + (table_size-1))/2)
                            v%pr(1)%cells(i,1,1) = vx_t((i - flamelet_pos + (table_size-1))/2)
                            do spec = 1, species_number
                                Y%pr(spec)%cells(i,1,1) = Y_t(spec,(i - flamelet_pos + (table_size-1))/2)
                            end do
                        else
                            ! Odd indices: linear interpolation
                            T%cells(i,1,1) = 0.5 * ( &
                                T_t((i - flamelet_pos + (table_size-1) + 1)/2) + &
                                T_t((i - flamelet_pos + (table_size-1) - 1)/2))
                            v%pr(1)%cells(i,1,1) = 0.5 * ( &
                                vx_t((i - flamelet_pos + (table_size-1) + 1)/2) + &
                                vx_t((i - flamelet_pos + (table_size-1) - 1)/2))
                            do spec = 1, species_number
                                Y%pr(spec)%cells(i,1,1) = 0.5 * ( &
                                    Y_t(spec,(i - flamelet_pos + (table_size-1) + 1)/2) + &
                                    Y_t(spec,(i - flamelet_pos + (table_size-1) - 1)/2))
                            end do
                        end if
                    end if
                    
                    ! Upstream of flamelet: use inlet conditions
                    if(i <= flamelet_pos - (table_size - 1)) then
                        T%cells(i,1,1) = T_t(0)
                        v%pr(1)%cells(i,1,1) = vx_t(0)
                        
                        ! Normalize species fractions
                        spec_summ = 0.0
                        do spec = 1, species_number
                            Y%pr(spec)%cells(i,1,1) = Y_t(spec,0)
                            spec_summ = spec_summ + Y%pr(spec)%cells(i,1,1)
                        end do
                        
                        do spec = 1, species_number
                            Y%pr(spec)%cells(i,1,1) = Y%pr(spec)%cells(i,1,1) / spec_summ
                        end do
                    end if
                    
                    ! Downstream of flamelet: use outlet conditions
                    if(i > flamelet_pos + (table_size - 1)) then
                        T%cells(i,1,1) = T_t(table_size-1)
                        v%pr(1)%cells(i,1,1) = vx_t(table_size-1)
                        do spec = 1, species_number
                            Y%pr(spec)%cells(i,1,1) = Y_t(spec, table_size-1)
                        end do
                    end if
                end do
            
            !------------------------------------------------
            ! CASE: NEAR-WALL FLAME (Wall at left boundary)
            !------------------------------------------------
            case('near_wall')
                ! Create hot region near wall (first quarter of domain)
                do i = 1, floor(domain_length / delta_x)
                    if (i < int(domain_length / delta_x / 4)) then
                        T%cells(i,:,:)        = 2000.0_dp      ! Hot near-wall region
                        Y%pr(1)%cells(i,:,:)  = 0.0_dp         ! No H2 (burned)
                        Y%pr(2)%cells(i,:,:)  = 0.0_dp         ! No O2 (burned)
                        Y%pr(3)%cells(i,:,:)  = nu * 3.762_dp  ! N2 only
                        Y%pr(7)%cells(i,:,:)  = 1.0_dp         ! H2O (products)
                    end if
                end do
                
                ! Convert mole fractions to dimensionless form
                call problem_thermophysics%change_field_units_mole_to_dimless(Y)
        
        end select
        
        !================================================================
        ! BOUNDARY CONDITIONS SETUP
        !================================================================
        
        select case(setup)
            
            !------------------------------------------------
            ! NEAR-WALL: Wall boundary + Outlet
            !------------------------------------------------
            case('near_wall')
                ! Left boundary: adiabatic slip wall
                call problem_boundaries%create_boundary_type( &
                    type_name               = 'wall',                &
                    slip                    = .true.,                &  ! Slip condition
                    conductive              = .false.,               &  ! Adiabatic
                    wall_temperature        = 0.0_dp,                &
                    wall_conductivity_ratio = 0.0_dp,                &
                    priority                = 1)
                
                ! Right boundary: outlet to ambient
                call problem_boundaries%create_boundary_type( &
                    type_name               = 'outlet',              &
                    farfield_pressure       = 1.0_dp * 101325.0_dp,  &  ! Ambient pressure
                    farfield_temperature    = 300.0_dp,              &  ! Ambient temperature
                    farfield_velocity       = 0.0_dp,                &  ! No mean flow
                    farfield_species_names  = [character(len=5) :: 'H2','O2','N2'], &
                    farfield_concentrations = (/1.0_dp, nu, nu * 3.762_dp/), &
                    priority                = 2)
            
            !------------------------------------------------
            ! COUNTER-FLOW: Inlet + Outlet (burned products)
            !------------------------------------------------
            case('counter_flow')
                ! Left boundary: fresh reactants inlet
                call problem_boundaries%create_boundary_type( &
                    type_name               = 'inlet',               &
                    farfield_pressure       = 1.0_dp * 101325.0_dp,  &
                    farfield_temperature    = 300.0_dp,              &  ! Cold reactants
                    farfield_velocity       = 0.0_dp,                &
                    farfield_species_names  = [character(len=5) :: 'H2','O2','N2'], &
                    farfield_concentrations = (/1.0_dp, nu, nu * 3.762_dp/), &
                    priority                = 1)
                
                ! Right boundary: hot products outlet
                call problem_boundaries%create_boundary_type( &
                    type_name               = 'outlet',              &
                    farfield_pressure       = 1.0_dp * 101325.0_dp,  &
                    farfield_temperature    = 1400.0_dp,             &  ! Hot products
                    farfield_velocity       = 0.0_dp,                &
                    farfield_species_names  = [character(len=5) :: 'H2O','N2'], &
                    farfield_concentrations = (/1.0_dp, nu * 3.762_dp/), &
                    priority                = 2)
            
            !------------------------------------------------
            ! COUNTER-FLOW PRECOMPUTED: Inlet/Outlet from flamelet
            !------------------------------------------------
            case('counter_flow_precInc')
                ! Left boundary: conditions from flamelet inlet
                call problem_boundaries%create_boundary_type( &
                    type_name               = 'inlet',               &
                    farfield_pressure       = 1.0_dp * 101325.0_dp,  &
                    farfield_temperature    = T_t(0),                &  ! From flamelet
                    farfield_velocity       = vx_t(0),               &  ! From flamelet
                    farfield_species_names  = [character(len=5) :: 'H2','O2','N2'], &
                    farfield_concentrations = (/1.0_dp, nu, nu * 3.762_dp/), &
                    priority                = 1)
                
                ! Right boundary: conditions from flamelet outlet
                call problem_boundaries%create_boundary_type( &
                    type_name               = 'outlet',              &
                    farfield_pressure       = 1.0_dp * 101325.0_dp,  &
                    farfield_temperature    = T_t(table_size-1),     &  ! From flamelet
                    farfield_velocity       = vx_t(table_size-1),    &  ! From flamelet
                    farfield_species_names  = [character(len=5) :: 'H2O','N2'], &
                    farfield_concentrations = (/1.0_dp, nu * 3.762_dp/), &
                    priority                = 2)
        
        end select
        
        ! Apply boundary markers to domain boundaries
        problem_boundaries%bc_markers(utter_loop(1,1),:,:) = 1  ! Left boundary
        problem_boundaries%bc_markers(utter_loop(1,2),:,:) = 2  ! Right boundary
        
        !================================================================
        ! LOG FILE ENTRY
        !================================================================
        write(log_unit,'(A)') 'General description: '
        write(log_unit,'(A)') 'Main aim: '
        write(log_unit,'(A)') '             '
        write(log_unit,'(A)') 'Problem setup:'
        write(log_unit,'(A)') '                 '
        write(log_unit,'(A)') '                 '
        write(log_unit,'(A)') 'Solver setup: '
        write(log_unit,'(A)') 'Validation and comparison:'
        write(log_unit,'(A)') '----------------------------------------------------------'
        
        close(log_unit)
        
        !================================================================
        ! INITIAL OUTPUT AND CLEANUP
        !================================================================
        ! Output initial solution and save data
        call problem_data_io%output_all_data(0.0_dp, stop_flag, make_output = .true.)
        call problem_data_save%save_all_data(0.0_dp, stop_flag, make_save   = .true.)
        
        ! Return to initial directory for next configuration
        ierr = chdir(initial_work_dir)
        
        ! Temporary pause/continue for debugging
        continue
        
    end do  ! task6 loop
    end do  ! task5 loop
    end do  ! task4 loop
    end do  ! task3 loop
    end do  ! task2 loop
    end do  ! task1 loop

end program package_interface
!================================================================================
! END OF PROGRAM
!================================================================================