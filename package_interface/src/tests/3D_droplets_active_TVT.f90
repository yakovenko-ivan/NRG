!===============================================================================
! NRG / research/active-matter
! TVT active-Janus-particle hydrodynamics: 3D numerical-reference configuration
!
! SOURCE: Senoshenko, Kiverin, Yakovenko et al., TVT draft supplied by user.
! Manuscript numerical case: 600 particles, dp=2 um, dD=50 dp=100 um.
!
! CAUTION: This is a NUMERICAL REFERENCE / LIQUID-LIKE SURROGATE, not yet an
! experiment-faithful oil/Janus model. FDS currently uses gas-mixture EOS;
! pressure is chosen to yield approximately 1000 kg/m^3 at 300 K.
! Mineral oil properties, Janus effective density, interfacial dynamics, laser
! heating/thermophoresis, and interparticle collisions are not resolved here.
! Do NOT confuse the 100-um numerical domain with the approximately 300-um
! particle cloud observed experimentally.
!
! Source architecture and namelist options audited against remote branch
! research/active-matter (2026-11-10). Windows/Intel ifx setup convention.
!===============================================================================
program package_interface

    use ifport
    use kind_parameters
    use global_data
    use nrg_build_info, only: write_nrg_source_revision
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

    implicit none

    type(computational_domain)     :: problem_domain
    type(data_manager)             :: problem_data_manager
    type(mpi_communications)       :: problem_mpi_support
    type(chemical_properties), target :: problem_chemistry
    type(thermophysical_properties), target :: problem_thermophysics
    type(solver_options)           :: problem_solver_options
    type(computational_mesh), target :: problem_mesh
    type(boundary_conditions), target :: problem_boundaries
    type(field_scalar_cons), target :: p, T, rho
    type(field_vector_cons), target :: v, Y
    type(particles_phase)          :: active_particles
    type(post_processor_manager)   :: problem_post_proc_manager
    type(data_io)                  :: problem_data_io
    type(data_save)                :: problem_data_save

    integer, parameter :: cells_per_side = 32
    integer, parameter :: particle_count = 600
    real(dp), parameter :: particle_diameter = 2.0e-6_dp       ! m: numerical case
    real(dp), parameter :: droplet_diameter = 50.0_dp * particle_diameter
    real(dp), parameter :: box_length = 1.2_dp * droplet_diameter
    real(dp), parameter :: reference_temperature = 300.0_dp  ! K: assumed
    real(dp), parameter :: target_fluid_density = 1000.0_dp  ! kg/m^3, surrogate
    real(dp), parameter :: reference_viscosity = 1.0e-3_dp   ! Pa s, nominal
    real(dp), parameter :: particle_density = 1500.0_dp     ! kg/m^3, PLACEHOLDER!
    real(dp), parameter :: target_particle_speed = 1.5e-6_dp ! m/s, between reported 1 and 2 um/s
    real(dp), parameter :: experimental_laser_wavelength = 532.0e-9_dp ! m, metadata
    real(dp), parameter :: experimental_laser_irradiance = 140.0e4_dp ! W/m^2, metadata
    real(dp), parameter :: boltzmann = 1.380649e-23_dp
    real(dp), parameter :: universal_R = 8.31446261815324_dp
    ! Mean molecular mass of oxygen/nitrogen mole ratio 1:3.762; GAS SURROGATE.
    real(dp), parameter :: apparent_molar_mass = &
        (0.031998_dp + 3.762_dp*0.0280134_dp)/(1.0_dp+3.762_dp)
    real(dp), parameter :: reference_pressure = &
        target_fluid_density * universal_R * reference_temperature / apparent_molar_mass
    real(dp), parameter :: estimated_active_force = &
        3.0_dp*pi*reference_viscosity*particle_diameter*target_particle_speed
    real(dp), parameter :: reference_rotational_diffusivity = &
        boltzmann*reference_temperature / &
        (8.0_dp*pi*reference_viscosity*(particle_diameter/2.0_dp)**3)

    ! This small timestep is a cautious STARTUP estimate for explicit viscous
    ! diffusion on a 32^3 mesh. Confirm using the actual FDS solver and mesh.
    real(dp), parameter :: initial_dt = 1.0e-6_dp   ! s; ~1/2 viscous stability estimate
    real(dp), parameter :: diagnostics_dt = 0.001_dp ! s; particle snapshots

    real(dp), dimension(3) :: cell_size, center
    integer, dimension(3,2) :: utter_loop
    real(dp) :: xx, yy, zz, fluid_radius, starting_radius, r2
    integer :: i, j, k, ierr, log_unit, fluid_cell_count
    logical :: stop_flag
    character(len=500) :: initial_work_dir
    character(len=100) :: work_dir
    character(len=20) :: solver_name

    ierr = getcwd(initial_work_dir)
    if (ierr /= 0) error stop 'TVT setup: cannot read starting directory'

    work_dir = 'tvt_active_3d'
    solver_name = 'fds_low_mach'

    ! Retain the existing NRG Windows interface convention. The leading
    ! 'echo d' disambiguates xcopy's directory prompt on a fresh checkout.
    ierr = system('if not exist ' // trim(work_dir) // ' mkdir ' // trim(work_dir))
    if (ierr /= 0) error stop 'TVT setup: cannot create case directory'
    ierr = system('echo d | xcopy .\' // trim(task_setup_folder) // ' .\' // &
        trim(work_dir) // trim(fold_sep) // trim(task_setup_folder) // ' /E/K/Y/I')
    if (ierr /= 0) error stop 'TVT setup: cannot copy task setup files'
    ierr = chdir(work_dir)
    if (ierr /= 0) error stop 'TVT setup: cannot enter case directory'

    open(newunit=log_unit, file=problem_setup_log_file, status='replace', form='formatted')
    call write_nrg_source_revision(log_unit, 'TVT active-matter problem setup')

    problem_domain = computational_domain_c( &
        dimensions=3, &
        cells_number=(/cells_per_side, cells_per_side, cells_per_side/), &
        coordinate_system='cartesian', &
        lengths=reshape((/0.0_dp, 0.0_dp, 0.0_dp, &
                          box_length, box_length, box_length/), (/3,2/)), &
        axis_names=(/'x','y','z'/))

    ! The present FDS carrier requires a gas thermodynamic backend. The
    ! reference oxygen/nitrogen mixture is used solely to construct the
    ! liquid-like density surrogate; it is NOT mineral oil.
    problem_chemistry = chemical_properties_c( &
        chemical_mechanism_file_name='hydrogen/KEROMNES.txt', &
        default_enhanced_efficiencies=1.0_dp, E_act_units='cal.mol')
    problem_thermophysics = thermophysical_properties_c( &
        chemistry=problem_chemistry, &
        thermo_data_file_name='KEROMNES_THERMO.txt', &
        transport_data_file_name='KEROMNES_TRANSDATA.txt', &
        molar_masses_data_file_name='molar_masses.dat')

    problem_solver_options = solver_options_c( &
        solver_name=solver_name, hydrodynamics_flag=.true., &
        heat_transfer_flag=.false., molecular_diffusion_flag=.false., &
        soret_diffusion_flag=.false., viscosity_flag=.true., &
        chemical_reaction_flag=.false., thermal_radiation_flag=.false., &
        grav_acc=(/0.0_dp, 0.0_dp, 0.0_dp/), &
        additional_particles_phases=1, &
        CFL_flag=.false., CFL_coefficient=0.25_dp, &
        initial_time_step=initial_dt)

    ! Unlike the legacy interface, each particle is an INERT SOLID, not a
    ! droplet undergoing thermal/phase conversion.
    active_particles%model = 'lagrangian'
    active_particles%initial_condition = 'random_sphere_nonoverlap'
    active_particles%boundary_interaction = 'reflect'
    active_particles%diameter = particle_diameter
    active_particles%material = 'Janus_MF_Cu'   ! label only; no chemistry in active mode
    active_particles%material_density = particle_density
    active_particles%heating = .false.
    active_particles%evaporating = .false.
    active_particles%inertial = .true.      ! enable two-way DRAG momentum coupling
    active_particles%initial_temperature = reference_temperature
    active_particles%initial_velocity = 0.0_dp
    active_particles%coordinates_are_relative = .true.
    active_particles%initial_position = (/0.5_dp,0.5_dp,0.5_dp/)
    active_particles%initial_particle_count = particle_count

    ! The initial cloud is inside the spherical region with clearance from
    ! voxel-mask boundaries. This radius is a numerical choice, NOT a measured
    ! experimental or manuscript value.
    fluid_radius = 0.5_dp*droplet_diameter
    active_particles%initial_region_radius = 0.80_dp*fluid_radius

    ! Photophoretic/thermophoretic propulsion is represented phenomenologically
    ! by constant thrust and 3D thermal rotational diffusion. No laser field
    ! is solved; the force is calibrated only to the target Stokes speed.
    active_particles%active_force_model = 'rotational_diffusion'
    active_particles%active_force_magnitude = estimated_active_force
    active_particles%active_rotational_diffusivity = reference_rotational_diffusivity
    active_particles%active_max_substep = initial_dt
    active_particles%active_random_seed = 20261110
    active_particles%output_interval = diagnostics_dt

    call problem_solver_options%create_additional_phase( &
        particles_parameters=active_particles)

    problem_mpi_support = mpi_communications_c(problem_domain)
    problem_data_manager = data_manager_c( &
        problem_domain, problem_mpi_support, problem_chemistry, &
        problem_thermophysics, problem_solver_options)

    call problem_data_manager%create_boundary_conditions( &
        problem_boundaries, number_of_boundary_types=1, default_boundary=1)
    call problem_data_manager%create_computational_mesh(problem_mesh)
    call problem_data_manager%create_scalar_field(p, 'pressure', 'p')
    call problem_data_manager%create_scalar_field(T, 'temperature', 'T')
    call problem_data_manager%create_scalar_field(rho, 'density', 'rho')
    call problem_data_manager%create_vector_field(v, 'velocity', 'v', 'spatial')
    call problem_data_manager%create_vector_field( &
        Y, 'specie_mass_fraction', 'Y', 'chemical')

    ! Active force does not create internal heat or mass source by assumption;
    ! viscous stresses and the particle -> liquid momentum acceleration remain.
    problem_post_proc_manager = post_processor_manager_c( &
        problem_data_manager, number_post_processors=0)

    problem_data_save = data_save_c( &
        problem_data_manager, &
        visible_fields_names=[character(len=40) :: &
            'pressure', 'pressure_dynamic', 'temperature', 'density', &
            'velocity', 'viscosity', 'velocity_production_particles01'], &
        save_time=10.0_dp, save_time_units='milliseconds', &
        save_format='tecplot', data_save_folder='data_save', &
        debug_flag=.false.)

    problem_data_io = data_io_c( &
        problem_data_manager, check_time=1.0_dp, check_time_units='seconds', &
        data_output_folder='data_output')

    ! Gas EOS surrogate: approximately target_fluid_density at T=300 K;
    ! DO NOT interpret reference_pressure as physical mineral-oil pressure.
    p%cells(:,:,:) = reference_pressure
    T%cells(:,:,:) = reference_temperature
    rho%cells(:,:,:) = target_fluid_density
    do i = 1, problem_chemistry%species_number
        Y%pr(i)%cells(:,:,:) = 0.0_dp
    end do
    ! Source Keromnes species numbering: 2=O2, 3=N2.
    Y%pr(2)%cells(:,:,:) = 1.0_dp
    Y%pr(3)%cells(:,:,:) = 3.762_dp
    call problem_thermophysics%change_field_units_mole_to_dimless(Y)

    call problem_boundaries%create_boundary_type( &
        type_name='wall', slip=.false., conductive=.false., &
        wall_temperature=0.0_dp, wall_conductivity_ratio=0.0_dp, &
        farfield_pressure=0.0_dp, farfield_temperature=0.0_dp, priority=1)

    cell_size = problem_mesh%get_cell_edges_length()
    utter_loop = problem_domain%get_global_utter_cells_bounds()
    center = 0.5_dp*box_length

    ! Immersed sphere represented by stair-step solid markers. Particle
    ! reflection against masked cells is approximate (not radius-aware).
    problem_boundaries%bc_markers(utter_loop(1,1),:,:) = 1
    problem_boundaries%bc_markers(utter_loop(1,2),:,:) = 1
    problem_boundaries%bc_markers(:,utter_loop(2,1),:) = 1
    problem_boundaries%bc_markers(:,utter_loop(2,2),:) = 1
    problem_boundaries%bc_markers(:,:,utter_loop(3,1)) = 1
    problem_boundaries%bc_markers(:,:,utter_loop(3,2)) = 1

    fluid_cell_count = 0
    do k = utter_loop(3,1), utter_loop(3,2)
        do j = utter_loop(2,1), utter_loop(2,2)
            do i = utter_loop(1,1), utter_loop(1,2)
                xx = problem_mesh%mesh(1,i,j,k) - center(1)
                yy = problem_mesh%mesh(2,i,j,k) - center(2)
                zz = problem_mesh%mesh(3,i,j,k) - center(3)
                r2 = xx**2 + yy**2 + zz**2
                if (r2 > fluid_radius**2) then
                    problem_boundaries%bc_markers(i,j,k) = 1
                else
                    if (problem_boundaries%bc_markers(i,j,k) == 0) &
                        fluid_cell_count = fluid_cell_count + 1
                end if
            end do
        end do
    end do
    if (fluid_cell_count == 0) error stop 'TVT setup: sphere contains no fluid cells'

    starting_radius = active_particles%initial_region_radius
    write(log_unit,'(A)') 'Case: TVT 3D active Janus particle numerical reference'
    write(log_unit,'(A,I0)') 'Number of active particles: ', particle_count
    write(log_unit,'(A,ES14.6)') 'Particle diameter [m]: ', particle_diameter
    write(log_unit,'(A,ES14.6)') 'Spherical region diameter [m]: ', droplet_diameter
    write(log_unit,'(A,ES14.6)') 'Initial cloud radius [m]: ', starting_radius
    write(log_unit,'(A,ES14.6)') 'Box side length [m]: ', box_length
    write(log_unit,'(A,I0)') 'Fluid cells (including any supported ghost layout): ', fluid_cell_count
    write(log_unit,'(A,3ES14.6)') 'Mesh spacing [m]: ', cell_size
    write(log_unit,'(A,ES14.6)') 'Assumed density [kg/m^3]: ', target_fluid_density
    write(log_unit,'(A,ES14.6)') 'Nominal viscosity [Pa s] (verify solver output): ', reference_viscosity
    write(log_unit,'(A,ES14.6)') 'Pseudo-gas EOS pressure [Pa]: ', reference_pressure
    write(log_unit,'(A,ES14.6)') 'Nominal Stokes speed [m/s]: ', target_particle_speed
    write(log_unit,'(A,ES14.6)') 'Active force amplitude [N]: ', estimated_active_force
    write(log_unit,'(A,ES14.6)') 'Rotational diffusivity [1/s]: ', reference_rotational_diffusivity
    write(log_unit,'(A,ES14.6)') 'Experimental laser wavelength [m]: ', experimental_laser_wavelength
    write(log_unit,'(A,ES14.6)') 'Experimental laser irradiance [W/m^2]: ', experimental_laser_irradiance
    write(log_unit,'(A,ES14.6)') 'Initial fluid timestep [s]: ', initial_dt
    write(log_unit,'(A)') 'LIMITATIONS: gas EOS surrogate; no true oil transport, laser field,'
    write(log_unit,'(A)') 'no Janus-specific active-slip flow, no particle collision model,'
    write(log_unit,'(A)') 'no dynamically deformable drop or surrounding liquid.'
    close(log_unit)

    call problem_data_io%output_all_data(0.0_dp,stop_flag,make_output=.true.)
    call problem_data_save%save_all_data(0.0_dp,stop_flag,make_save=.true.)

    ierr = chdir(initial_work_dir)
    if (ierr /= 0) error stop 'TVT setup: cannot return to starting directory'
end program package_interface
