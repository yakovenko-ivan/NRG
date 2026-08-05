module lagrangian_particles_solver_class

    use kind_parameters, only: dp
    use global_data, only: I_m, pi, T_ref
    use field_pointers
    use boundary_conditions_class
    use data_manager_class
    use computational_mesh_class
    use computational_domain_class
    use thermophysical_properties_class
    use chemical_properties_class
    use solver_options_class, only: particles_phase
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite

    implicit none

    private
    public :: lagrangian_particles_solver, lagrangian_particles_solver_c

    real(dp), parameter :: parcel_cfl_target = 0.5_dp
    real(dp), parameter :: particle_mass_floor = 1.0e-30_dp

    !> One computational parcel.  A parcel represents one physical particle in
    !> the present implementation; statistical parcel weights can be introduced
    !> later without changing the public facade.
    type :: lagrangian_particle
        real(dp), dimension(3) :: coordinates = 0.0_dp
        real(dp), dimension(3) :: previous_coordinates = 0.0_dp
        real(dp), dimension(3) :: initial_coordinates = 0.0_dp
        real(dp), dimension(3) :: velocity = 0.0_dp
        real(dp) :: temperature = 300.0_dp
        real(dp) :: mass = 0.0_dp
        logical :: active = .false.
    end type lagrangian_particle

    !> Lagrangian dispersed-phase backend.
    !>
    !> Coupling-field contract shared with the Eulerian backend:
    !>   gas_mass_source      [kg m^-3 s^-1]
    !>   gas_species_source   [kg m^-3 s^-1]
    !>   gas_momentum_source  [m s^-2] (gas acceleration)
    !>   gas_energy_source    [J m^-3 s^-1]
    type :: lagrangian_particles_solver
        private

        type(field_scalar_cons_pointer) :: gas_temperature
        type(field_scalar_cons_pointer) :: gas_density
        type(field_scalar_cons_pointer) :: gas_viscosity
        type(field_scalar_cons_pointer) :: gas_conductivity
        type(field_vector_flow_pointer) :: gas_face_velocity

        type(field_scalar_cons_pointer) :: gas_mass_source
        type(field_scalar_cons_pointer) :: gas_energy_source
        type(field_vector_cons_pointer) :: gas_momentum_source
        type(field_vector_cons_pointer) :: gas_species_source

        type(field_scalar_cons), pointer :: gas_mass_source_store => null()
        type(field_scalar_cons), pointer :: gas_energy_source_store => null()
        type(field_vector_cons), pointer :: gas_momentum_source_store => null()
        type(field_vector_cons), pointer :: gas_species_source_store => null()

        type(computational_domain) :: domain
        type(boundary_conditions_pointer) :: boundary
        type(computational_mesh_pointer) :: mesh
        type(chemical_properties_pointer) :: chemistry
        type(thermophysical_properties_pointer) :: thermophysics

        type(particles_phase) :: parameters
        type(lagrangian_particle), dimension(:), allocatable :: particles
        real(dp), dimension(3) :: gravity = 0.0_dp
        integer :: phase_number = 0
        integer :: material_index = 0
        integer :: output_unit = -1
        real(dp) :: next_output_time = 0.0_dp
        real(dp) :: legacy_time = 0.0_dp
    contains
        procedure :: set_initial_distributions
        procedure :: advance

        ! Legacy alias retained while gas solvers migrate to the facade API.
        procedure :: particles_solve => advance_legacy
        procedure :: apply_boundary_conditions_main

        procedure, private :: advance_legacy
        procedure, private :: reset_coupling_sources
        procedure, private :: locate_particle_cell
        procedure, private :: interpolate_gas_velocity
        procedure, private :: apply_particle_boundary
        procedure, private :: particle_mass_from_diameter
        procedure, private :: particle_diameter_from_mass
        procedure, private :: drag_relaxation_rate
        procedure, private :: cell_volume
        procedure, private :: write_snapshot
    end type lagrangian_particles_solver

    interface lagrangian_particles_solver_c
        module procedure constructor
    end interface lagrangian_particles_solver_c

contains

    type(lagrangian_particles_solver) function constructor( &
        manager, particles, phase_number)
        type(data_manager), intent(inout) :: manager
        type(particles_phase), intent(in) :: particles
        integer, intent(in) :: phase_number

        type(field_scalar_cons_pointer) :: scalar_pointer
        type(field_vector_cons_pointer) :: vector_pointer
        type(field_tensor_cons_pointer) :: tensor_pointer
        type(field_scalar_flow_pointer) :: scalar_flow_pointer
        type(field_vector_flow_pointer) :: vector_flow_pointer
        type(field_tensor_flow_pointer) :: tensor_flow_pointer
        character(len=48) :: long_name, short_name, output_name

        constructor%parameters = particles
        constructor%phase_number = phase_number
        constructor%domain = manager%domain
        constructor%mesh%mesh_ptr => manager%computational_mesh_pointer%mesh_ptr
        constructor%boundary%bc_ptr => manager%boundary_conditions_pointer%bc_ptr
        constructor%chemistry%chem_ptr => manager%chemistry%chem_ptr
        constructor%thermophysics%thermo_ptr => manager%thermophysics%thermo_ptr
        constructor%gravity = manager%solver_options%get_grav_acc()

        call manager%get_cons_field_pointer_by_name( &
            scalar_pointer, vector_pointer, tensor_pointer, 'temperature')
        constructor%gas_temperature%s_ptr => scalar_pointer%s_ptr
        call manager%get_cons_field_pointer_by_name( &
            scalar_pointer, vector_pointer, tensor_pointer, 'density')
        constructor%gas_density%s_ptr => scalar_pointer%s_ptr
        call manager%get_cons_field_pointer_by_name( &
            scalar_pointer, vector_pointer, tensor_pointer, 'viscosity')
        constructor%gas_viscosity%s_ptr => scalar_pointer%s_ptr
        call manager%get_cons_field_pointer_by_name( &
            scalar_pointer, vector_pointer, tensor_pointer, 'thermal_conductivity')
        constructor%gas_conductivity%s_ptr => scalar_pointer%s_ptr
        call manager%get_flow_field_pointer_by_name( &
            scalar_flow_pointer, vector_flow_pointer, tensor_flow_pointer, &
            'velocity_flow')
        constructor%gas_face_velocity%v_ptr => vector_flow_pointer%v_ptr

        allocate(constructor%gas_mass_source_store)
        allocate(constructor%gas_energy_source_store)
        allocate(constructor%gas_momentum_source_store)
        allocate(constructor%gas_species_source_store)

        call phase_field_names('density_production_particles', 'rho_prod_p', &
            phase_number, long_name, short_name)
        call manager%create_scalar_field( &
            constructor%gas_mass_source_store, long_name, short_name)
        constructor%gas_mass_source%s_ptr => constructor%gas_mass_source_store

        call phase_field_names('energy_production_particles', 'E_prod_p', &
            phase_number, long_name, short_name)
        call manager%create_scalar_field( &
            constructor%gas_energy_source_store, long_name, short_name)
        constructor%gas_energy_source%s_ptr => constructor%gas_energy_source_store

        call phase_field_names('velocity_production_particles', 'v_prod_p', &
            phase_number, long_name, short_name)
        call manager%create_vector_field( &
            constructor%gas_momentum_source_store, long_name, short_name, &
            'spatial')
        constructor%gas_momentum_source%v_ptr => &
            constructor%gas_momentum_source_store

        call phase_field_names('concentration_production_particles', 'Y_prod_p', &
            phase_number, long_name, short_name)
        call manager%create_vector_field( &
            constructor%gas_species_source_store, long_name, short_name, &
            'chemical')
        constructor%gas_species_source%v_ptr => &
            constructor%gas_species_source_store

        constructor%material_index = &
            constructor%chemistry%chem_ptr%get_chemical_specie_index( &
                trim(particles%material))
        if (constructor%material_index <= 0) then
            error stop 'lagrangian particles: material species was not found'
        end if

        if (particles%output_interval > 0.0_dp) then
            write(output_name, '(A,I2.2,A)') 'lagrangian_particles_phase_', &
                phase_number, '.dat'
            open(newunit=constructor%output_unit, file=trim(output_name), &
                status='replace', form='formatted')
            write(constructor%output_unit, '(A)') &
                'VARIABLES="time" "particle" "x" "y" "z" ' // &
                '"u" "v" "w" "temperature" "mass"'
        end if
    end function constructor


    !> Create parcels from interface data.  No source-code edits are required:
    !> single-particle and Cartesian-lattice initializations are selected by
    !> particles_phase%initial_condition.
    subroutine set_initial_distributions(this)
        class(lagrangian_particles_solver), intent(inout) :: this

        real(dp), dimension(:,:), allocatable :: lengths
        real(dp), dimension(3) :: lower, upper, position
        integer, dimension(3) :: lattice, index
        integer :: ix, iy, iz
        integer :: dimensions, number_of_particles, particle, dim
        character(len=24) :: initialization
        logical :: outside
        integer, dimension(3) :: cell

        dimensions = this%domain%get_domain_dimensions()
        lengths = this%domain%get_domain_lengths()
        initialization = lowercase(trim(this%parameters%initial_condition))

        select case (initialization)
        case ('none')
            number_of_particles = 0
        case ('single')
            number_of_particles = 1
        case ('lattice')
            lattice = this%parameters%initial_lattice
            if (dimensions < 3) lattice(3) = 1
            if (dimensions < 2) lattice(2) = 1
            number_of_particles = product(lattice)
        case default
            error stop 'lagrangian particles: initialization must be none, single, or lattice'
        end select

        if (allocated(this%particles)) deallocate(this%particles)
        allocate(this%particles(number_of_particles))
        if (number_of_particles == 0) then
            call this%reset_coupling_sources()
            return
        end if

        lower = 0.0_dp
        upper = 0.0_dp
        do dim = 1, dimensions
            if (this%parameters%coordinates_are_relative) then
                lower(dim) = lengths(dim,1) + &
                    this%parameters%initial_region_min(dim) * &
                    (lengths(dim,2) - lengths(dim,1))
                upper(dim) = lengths(dim,1) + &
                    this%parameters%initial_region_max(dim) * &
                    (lengths(dim,2) - lengths(dim,1))
            else
                lower(dim) = this%parameters%initial_region_min(dim)
                upper(dim) = this%parameters%initial_region_max(dim)
            end if
        end do

        lattice = this%parameters%initial_lattice
        if (dimensions < 3) lattice(3) = 1
        if (dimensions < 2) lattice(2) = 1
        particle = 0
        do iz = 1, lattice(3)
            index(3) = iz
            do iy = 1, lattice(2)
                index(2) = iy
                do ix = 1, lattice(1)
                    index(1) = ix
                    particle = particle + 1
                    if (initialization == 'single') then
                        do dim = 1, dimensions
                            if (this%parameters%coordinates_are_relative) then
                                position(dim) = lengths(dim,1) + &
                                    this%parameters%initial_position(dim) * &
                                    (lengths(dim,2) - lengths(dim,1))
                            else
                                position(dim) = &
                                    this%parameters%initial_position(dim)
                            end if
                        end do
                    else
                        do dim = 1, dimensions
                            position(dim) = lower(dim) + &
                                (real(index(dim),dp) - 0.5_dp) * &
                                (upper(dim) - lower(dim)) / real(lattice(dim),dp)
                        end do
                    end if
                    if (dimensions < 3) position(3) = 0.5_dp * &
                        (lengths(3,1) + lengths(3,2))
                    if (dimensions < 2) position(2) = 0.5_dp * &
                        (lengths(2,1) + lengths(2,2))

                    this%particles(particle)%coordinates = position
                    this%particles(particle)%previous_coordinates = position
                    this%particles(particle)%initial_coordinates = position
                    this%particles(particle)%velocity = &
                        this%parameters%initial_velocity
                    this%particles(particle)%temperature = &
                        this%parameters%initial_temperature
                    this%particles(particle)%mass = &
                        this%particle_mass_from_diameter(this%parameters%diameter)
                    cell = this%locate_particle_cell(position, outside)
                    this%particles(particle)%active = .not.outside
                    if (.not.this%particles(particle)%active) then
                        error stop 'lagrangian particles: initial parcel is outside the fluid domain'
                    end if
                end do
            end do
        end do

        this%next_output_time = 0.0_dp
        call this%reset_coupling_sources()
        if (this%output_unit > 0) call this%write_snapshot(0.0_dp)
    end subroutine set_initial_distributions


    !> Advance parcels and accumulate conservative two-way coupling sources.
    !> Drag is integrated analytically over each parcel substep.  The reaction
    !> impulse is deposited into the gas; heat and evaporated mass are deposited
    !> as volumetric rates over the full gas time step.
    subroutine advance(this, time_step, time)
        class(lagrangian_particles_solver), intent(inout) :: this
        real(dp), intent(in) :: time_step
        real(dp), intent(in) :: time

        real(dp), dimension(3) :: dx, gas_velocity, velocity_old
        real(dp), dimension(3) :: velocity_new, drag_impulse, position_trial
        real(dp) :: substep, speed_ratio, relative_speed, beta, exponential
        real(dp) :: gas_density, gas_viscosity, gas_conductivity
        real(dp) :: particle_diameter, reynolds, nusselt, heat_transfer_rate
        real(dp) :: heat, heat_to_boiling, evaporated_mass, vapor_enthalpy
        real(dp) :: volume, kinetic_energy, velocity_mid, mass_old
        integer, dimension(3) :: cell
        integer :: dimensions, particle, dim, subcycle, subcycles
        logical :: outside

        if (.not.ieee_is_finite(time_step) .or. time_step <= 0.0_dp) then
            error stop 'lagrangian particles: non-positive time step'
        end if

        dimensions = this%domain%get_domain_dimensions()
        dx = this%mesh%mesh_ptr%get_cell_edges_length()
        call this%reset_coupling_sources()

        associate(Tg => this%gas_temperature%s_ptr, &
                  rho => this%gas_density%s_ptr, &
                  mu => this%gas_viscosity%s_ptr, &
                  conductivity => this%gas_conductivity%s_ptr, &
                  mass_src => this%gas_mass_source%s_ptr, &
                  energy_src => this%gas_energy_source%s_ptr, &
                  momentum_src => this%gas_momentum_source%v_ptr, &
                  species_src => this%gas_species_source%v_ptr)
        do particle = 1, size(this%particles)
            if (.not.this%particles(particle)%active) cycle

            speed_ratio = 0.0_dp
            do dim = 1, dimensions
                speed_ratio = max(speed_ratio, &
                    abs(this%particles(particle)%velocity(dim)) / dx(dim))
            end do
            subcycles = max(1, ceiling(time_step * speed_ratio / parcel_cfl_target))
            substep = time_step / real(subcycles, dp)

            do subcycle = 1, subcycles
                cell = this%locate_particle_cell( &
                    this%particles(particle)%coordinates, outside)
                if (outside) then
                    this%particles(particle)%active = .false.
                    exit
                end if

                gas_density = rho%cells(cell(1),cell(2),cell(3))
                gas_viscosity = mu%cells(cell(1),cell(2),cell(3))
                gas_conductivity = conductivity%cells(cell(1),cell(2),cell(3))
                if (gas_density <= 0.0_dp .or. gas_viscosity <= 0.0_dp) then
                    error stop 'lagrangian particles: invalid gas transport state'
                end if

                call this%interpolate_gas_velocity( &
                    this%particles(particle)%coordinates, cell, gas_velocity)
                velocity_old = this%particles(particle)%velocity
                relative_speed = sqrt(sum( &
                    (velocity_old(1:dimensions) - gas_velocity(1:dimensions))**2))
                particle_diameter = this%particle_diameter_from_mass( &
                    this%particles(particle)%mass)
                beta = this%drag_relaxation_rate( &
                    gas_density, gas_viscosity, particle_diameter, &
                    this%particles(particle)%mass, relative_speed)

                ! Analytical solution of dv_p/dt = -beta (v_p-v_g) + g.
                ! The inertial option controls two-way drag feedback, not the
                ! parcel equation of motion; this preserves the legacy meaning.
                if (beta > sqrt(tiny(1.0_dp))) then
                    exponential = exp(-beta * substep)
                    do dim = 1, dimensions
                        velocity_new(dim) = gas_velocity(dim) + &
                            (velocity_old(dim) - gas_velocity(dim) - &
                            this%gravity(dim) / beta) * exponential + &
                            this%gravity(dim) / beta
                    end do
                else
                    velocity_new = velocity_old + this%gravity * substep
                end if

                ! Separate the gravity impulse from the particle momentum change;
                ! only the drag impulse is an interphase exchange.
                do dim = 1, dimensions
                    drag_impulse(dim) = this%particles(particle)%mass * ( &
                        velocity_new(dim) - velocity_old(dim) - &
                        this%gravity(dim) * substep)
                    position_trial(dim) = &
                        this%particles(particle)%coordinates(dim) + &
                        0.5_dp * substep * &
                        (velocity_old(dim) + velocity_new(dim))
                end do

                this%particles(particle)%previous_coordinates = &
                    this%particles(particle)%coordinates
                this%particles(particle)%coordinates = position_trial
                this%particles(particle)%velocity = velocity_new
                call this%apply_particle_boundary(particle)
                if (.not.this%particles(particle)%active) exit

                volume = this%cell_volume(cell)
                if (this%parameters%inertial) then
                    do dim = 1, dimensions
                        momentum_src%pr(dim)%cells(cell(1),cell(2),cell(3)) = &
                            momentum_src%pr(dim)%cells(cell(1),cell(2),cell(3)) - &
                            drag_impulse(dim) / (volume * time_step * gas_density)
                        ! Use the pre-boundary particle velocity for drag work.
                        ! A wall collision transfers momentum/energy to the wall,
                        ! not to the carrier gas.
                        velocity_mid = 0.5_dp * &
                            (velocity_old(dim) + velocity_new(dim))
                        energy_src%cells(cell(1),cell(2),cell(3)) = &
                            energy_src%cells(cell(1),cell(2),cell(3)) - &
                            drag_impulse(dim) * velocity_mid / (volume * time_step)
                    end do
                end if

                ! Lumped-capacitance heating with a boiling-energy limiter.
                ! Evaporated material enters the gas species source at parcel
                ! velocity and with its EOS-consistent sensible enthalpy.
                if (this%parameters%heating) then
                    reynolds = gas_density * relative_speed * &
                        particle_diameter / gas_viscosity
                    nusselt = 2.0_dp + 0.6_dp * sqrt(max(reynolds, 0.0_dp)) * &
                        0.7_dp**(1.0_dp / 3.0_dp)
                    heat_transfer_rate = nusselt * gas_conductivity * &
                        particle_surface_area(particle_diameter, dimensions) * &
                        (Tg%cells(cell(1),cell(2),cell(3)) - &
                        this%particles(particle)%temperature) / particle_diameter
                    heat = heat_transfer_rate * substep
                    mass_old = this%particles(particle)%mass
                    evaporated_mass = 0.0_dp

                    if (this%parameters%evaporating .and. heat > 0.0_dp) then
                        heat_to_boiling = mass_old * &
                            this%parameters%material_heat_capacity * max( &
                            this%parameters%material_boiling_temperature - &
                            this%particles(particle)%temperature, 0.0_dp)
                        if (heat <= heat_to_boiling) then
                            this%particles(particle)%temperature = &
                                this%particles(particle)%temperature + heat / &
                                (mass_old * this%parameters%material_heat_capacity)
                        else
                            this%particles(particle)%temperature = &
                                this%parameters%material_boiling_temperature
                            evaporated_mass = min(mass_old, &
                                (heat - heat_to_boiling) / &
                                this%parameters%material_latent_heat)
                        end if
                    else
                        this%particles(particle)%temperature = max(1.0_dp, &
                            this%particles(particle)%temperature + heat / &
                            (mass_old * this%parameters%material_heat_capacity))
                    end if

                    energy_src%cells(cell(1),cell(2),cell(3)) = &
                        energy_src%cells(cell(1),cell(2),cell(3)) - &
                        heat / (volume * time_step)

                    if (evaporated_mass > 0.0_dp) then
                        this%particles(particle)%mass = mass_old - evaporated_mass
                        mass_src%cells(cell(1),cell(2),cell(3)) = &
                            mass_src%cells(cell(1),cell(2),cell(3)) + &
                            evaporated_mass / (volume * time_step)
                        species_src%pr(this%material_index)%cells( &
                            cell(1),cell(2),cell(3)) = &
                            species_src%pr(this%material_index)%cells( &
                            cell(1),cell(2),cell(3)) + &
                            evaporated_mass / (volume * time_step)

                        vapor_enthalpy = species_sensible_enthalpy( &
                            this%thermophysics%thermo_ptr, this%material_index, &
                            this%particles(particle)%temperature)
                        kinetic_energy = 0.5_dp * sum( &
                            this%particles(particle)%velocity(1:dimensions)**2)
                        energy_src%cells(cell(1),cell(2),cell(3)) = &
                            energy_src%cells(cell(1),cell(2),cell(3)) + &
                            evaporated_mass * (vapor_enthalpy + kinetic_energy) / &
                            (volume * time_step)
                        do dim = 1, dimensions
                            momentum_src%pr(dim)%cells(cell(1),cell(2),cell(3)) = &
                                momentum_src%pr(dim)%cells(cell(1),cell(2),cell(3)) + &
                                evaporated_mass * &
                                this%particles(particle)%velocity(dim) / &
                                (volume * time_step * gas_density)
                        end do

                        if (this%particles(particle)%mass <= particle_mass_floor) then
                            this%particles(particle)%active = .false.
                            exit
                        end if
                    end if
                end if
            end do
        end do
        end associate

        if (this%output_unit > 0 .and. &
            time + time_step >= this%next_output_time) then
            call this%write_snapshot(time + time_step)
            this%next_output_time = this%next_output_time + &
                this%parameters%output_interval
        end if
    end subroutine advance


    subroutine advance_legacy(this, time_step)
        class(lagrangian_particles_solver), intent(inout) :: this
        real(dp), intent(in) :: time_step
        call this%advance(time_step, this%legacy_time)
        this%legacy_time = this%legacy_time + time_step
    end subroutine advance_legacy


    subroutine apply_boundary_conditions_main(this, time)
        class(lagrangian_particles_solver), intent(inout) :: this
        real(dp), intent(in) :: time
        ! Lagrangian boundary conditions are applied geometrically after each
        ! parcel substep.  Retaining this no-op binding avoids breaking old gas
        ! solvers during migration.
        if (.not.ieee_is_finite(time) .or. this%phase_number <= 0) then
            error stop 'lagrangian particles: non-finite boundary time'
        end if
    end subroutine apply_boundary_conditions_main


    subroutine reset_coupling_sources(this)
        class(lagrangian_particles_solver), intent(inout) :: this
        integer :: dimensions, species_number, dim, species

        dimensions = this%domain%get_domain_dimensions()
        species_number = this%chemistry%chem_ptr%species_number
        this%gas_mass_source%s_ptr%cells = 0.0_dp
        this%gas_energy_source%s_ptr%cells = 0.0_dp
        do dim = 1, dimensions
            this%gas_momentum_source%v_ptr%pr(dim)%cells = 0.0_dp
        end do
        do species = 1, species_number
            this%gas_species_source%v_ptr%pr(species)%cells = 0.0_dp
        end do
    end subroutine reset_coupling_sources


    !> Locate a parcel by physical coordinates.  The search uses actual mesh
    !> centers and therefore does not hard-code a global cell index or particle
    !> release location.
    function locate_particle_cell(this, coordinates, outside) result(cell)
        class(lagrangian_particles_solver), intent(in) :: this
        real(dp), dimension(3), intent(in) :: coordinates
        logical, intent(out) :: outside
        integer, dimension(3) :: cell

        integer, dimension(3,2) :: loop
        real(dp), dimension(3) :: dx
        integer :: dimensions, dim, index
        integer :: ii, jj, kk
        real(dp) :: center, lower_face, upper_face

        dimensions = this%domain%get_domain_dimensions()
        loop = this%domain%get_local_inner_cells_bounds()
        dx = this%mesh%mesh_ptr%get_cell_edges_length()
        cell = (/loop(1,1), loop(2,1), loop(3,1)/)
        outside = .false.

        do dim = 1, dimensions
            cell(dim) = 0
            do index = loop(dim,1), loop(dim,2)
                ii = loop(1,1)
                jj = loop(2,1)
                kk = loop(3,1)
                if (dim == 1) ii = index
                if (dim == 2) jj = index
                if (dim == 3) kk = index
                center = this%mesh%mesh_ptr%mesh(dim,ii,jj,kk)
                lower_face = center - 0.5_dp * dx(dim)
                upper_face = center + 0.5_dp * dx(dim)
                if (coordinates(dim) >= lower_face .and. &
                    coordinates(dim) < upper_face) then
                    cell(dim) = index
                    exit
                end if
                if (index == loop(dim,2) .and. &
                    abs(coordinates(dim) - upper_face) <= &
                    10.0_dp * spacing(max(abs(upper_face), 1.0_dp))) then
                    cell(dim) = index
                    exit
                end if
            end do
            if (cell(dim) == 0) then
                outside = .true.
                cell(dim) = loop(dim,1)
            end if
        end do

        if (.not.outside) then
            if (this%boundary%bc_ptr%bc_markers( &
                cell(1),cell(2),cell(3)) /= 0) outside = .true.
        end if
    end function locate_particle_cell


    subroutine interpolate_gas_velocity(this, coordinates, cell, velocity)
        class(lagrangian_particles_solver), intent(in) :: this
        real(dp), dimension(3), intent(in) :: coordinates
        integer, dimension(3), intent(in) :: cell
        real(dp), dimension(3), intent(out) :: velocity

        real(dp), dimension(3) :: dx
        integer :: dimensions, dim, ir, jr, kr
        real(dp) :: center, lower_face, fraction, left_value, right_value

        dimensions = this%domain%get_domain_dimensions()
        dx = this%mesh%mesh_ptr%get_cell_edges_length()
        velocity = 0.0_dp

        do dim = 1, dimensions
            center = this%mesh%mesh_ptr%mesh( &
                dim, cell(1), cell(2), cell(3))
            lower_face = center - 0.5_dp * dx(dim)
            fraction = min(1.0_dp, max(0.0_dp, &
                (coordinates(dim) - lower_face) / dx(dim)))
            ir = cell(1) + I_m(dim,1)
            jr = cell(2) + I_m(dim,2)
            kr = cell(3) + I_m(dim,3)
            left_value = this%gas_face_velocity%v_ptr%pr(dim)%cells( &
                dim, cell(1), cell(2), cell(3))
            right_value = this%gas_face_velocity%v_ptr%pr(dim)%cells( &
                dim, ir, jr, kr)
            velocity(dim) = (1.0_dp - fraction) * left_value + &
                fraction * right_value
        end do
    end subroutine interpolate_gas_velocity


    !> Apply the runtime-selected wall interaction: reflect, stick, or escape.
    !> Domain boundaries and cells marked as non-fluid are handled without
    !> recompiling the particle module.
    subroutine apply_particle_boundary(this, particle_number)
        class(lagrangian_particles_solver), intent(inout) :: this
        integer, intent(in) :: particle_number

        real(dp), dimension(:,:), allocatable :: lengths
        integer :: dimensions, dim
        character(len=24) :: interaction
        logical :: outside
        integer, dimension(3) :: cell

        dimensions = this%domain%get_domain_dimensions()
        lengths = this%domain%get_domain_lengths()
        interaction = lowercase(trim(this%parameters%boundary_interaction))

        do dim = 1, dimensions
            if (this%particles(particle_number)%coordinates(dim) < lengths(dim,1)) then
                select case (interaction)
                case ('reflect')
                    this%particles(particle_number)%coordinates(dim) = &
                        2.0_dp * lengths(dim,1) - &
                        this%particles(particle_number)%coordinates(dim)
                    this%particles(particle_number)%velocity(dim) = &
                        -this%particles(particle_number)%velocity(dim)
                case ('stick')
                    this%particles(particle_number)%coordinates(dim) = &
                        lengths(dim,1) + spacing(lengths(dim,1))
                    this%particles(particle_number)%velocity = 0.0_dp
                case ('escape')
                    this%particles(particle_number)%active = .false.
                    return
                end select
            else if (this%particles(particle_number)%coordinates(dim) > &
                lengths(dim,2)) then
                select case (interaction)
                case ('reflect')
                    this%particles(particle_number)%coordinates(dim) = &
                        2.0_dp * lengths(dim,2) - &
                        this%particles(particle_number)%coordinates(dim)
                    this%particles(particle_number)%velocity(dim) = &
                        -this%particles(particle_number)%velocity(dim)
                case ('stick')
                    this%particles(particle_number)%coordinates(dim) = &
                        lengths(dim,2) - spacing(lengths(dim,2))
                    this%particles(particle_number)%velocity = 0.0_dp
                case ('escape')
                    this%particles(particle_number)%active = .false.
                    return
                end select
            end if
        end do

        cell = this%locate_particle_cell( &
            this%particles(particle_number)%coordinates, outside)
        if (outside) then
            select case (interaction)
            case ('reflect')
                this%particles(particle_number)%coordinates = &
                    this%particles(particle_number)%previous_coordinates
                this%particles(particle_number)%velocity = &
                    -this%particles(particle_number)%velocity
            case ('stick')
                this%particles(particle_number)%coordinates = &
                    this%particles(particle_number)%previous_coordinates
                this%particles(particle_number)%velocity = 0.0_dp
            case ('escape')
                this%particles(particle_number)%active = .false.
            end select
        end if
    end subroutine apply_particle_boundary


    function particle_mass_from_diameter(this, diameter) result(mass)
        class(lagrangian_particles_solver), intent(in) :: this
        real(dp), intent(in) :: diameter
        real(dp) :: mass
        integer :: dimensions

        dimensions = this%domain%get_domain_dimensions()
        if (dimensions == 2) then
            mass = 0.25_dp * pi * diameter**2 * &
                this%parameters%material_density
        else
            mass = pi * diameter**3 * this%parameters%material_density / 6.0_dp
        end if
    end function particle_mass_from_diameter


    function particle_diameter_from_mass(this, mass) result(diameter)
        class(lagrangian_particles_solver), intent(in) :: this
        real(dp), intent(in) :: mass
        real(dp) :: diameter
        integer :: dimensions

        dimensions = this%domain%get_domain_dimensions()
        if (dimensions == 2) then
            diameter = sqrt(4.0_dp * max(mass, 0.0_dp) / &
                (pi * this%parameters%material_density))
        else
            diameter = (6.0_dp * max(mass, 0.0_dp) / &
                (pi * this%parameters%material_density))**(1.0_dp / 3.0_dp)
        end if
    end function particle_diameter_from_mass


    function drag_relaxation_rate( &
        this, gas_density, viscosity, diameter, mass, relative_speed) result(beta)
        class(lagrangian_particles_solver), intent(in) :: this
        real(dp), intent(in) :: gas_density, viscosity, diameter, mass
        real(dp), intent(in) :: relative_speed
        real(dp) :: beta, reynolds, drag_coefficient, projected_area
        integer :: dimensions

        dimensions = this%domain%get_domain_dimensions()
        if (relative_speed <= sqrt(tiny(1.0_dp))) then
            if (dimensions == 3) then
                beta = 3.0_dp * pi * viscosity * diameter / mass
            else
                beta = 0.0_dp
            end if
            return
        end if

        reynolds = gas_density * relative_speed * diameter / viscosity
        if (dimensions == 2) then
            if (reynolds < 1.0_dp) then
                drag_coefficient = 10.0_dp / max(reynolds, 1.0e-12_dp)**0.8_dp
            else if (reynolds < 1000.0_dp) then
                drag_coefficient = 10.0_dp * &
                    (0.6_dp + 0.4_dp * reynolds**0.8_dp) / reynolds
            else
                drag_coefficient = 1.0_dp
            end if
            projected_area = diameter
        else
            if (reynolds < 1000.0_dp) then
                drag_coefficient = 24.0_dp * &
                    (1.0_dp + 0.15_dp * reynolds**0.687_dp) / &
                    max(reynolds, 1.0e-12_dp)
            else
                drag_coefficient = 0.44_dp
            end if
            projected_area = 0.25_dp * pi * diameter**2
        end if
        beta = 0.5_dp * gas_density * drag_coefficient * projected_area * &
            relative_speed / mass
    end function drag_relaxation_rate


    function cell_volume(this, cell) result(volume)
        class(lagrangian_particles_solver), intent(in) :: this
        integer, dimension(3), intent(in) :: cell
        real(dp) :: volume, radius
        character(len=20) :: coordinate_system

        volume = this%mesh%mesh_ptr%get_cell_volume()
        coordinate_system = trim(this%domain%get_coordinate_system_name())
        select case (coordinate_system)
        case ('cartesian')
        case ('cylindrical')
            radius = this%mesh%mesh_ptr%mesh(1,cell(1),cell(2),cell(3))
            volume = volume * max(radius, 0.0_dp)
        case ('spherical')
            radius = this%mesh%mesh_ptr%mesh(1,cell(1),cell(2),cell(3))
            volume = volume * max(radius, 0.0_dp)**2
        case default
            error stop 'lagrangian particles: unknown coordinate system'
        end select
    end function cell_volume


    subroutine write_snapshot(this, time)
        class(lagrangian_particles_solver), intent(inout) :: this
        real(dp), intent(in) :: time
        integer :: particle

        if (this%output_unit <= 0) return
        do particle = 1, size(this%particles)
            if (.not.this%particles(particle)%active) cycle
            write(this%output_unit, '(ES14.6,1X,I0,1X,8ES14.6)') &
                time, particle, this%particles(particle)%coordinates, &
                this%particles(particle)%velocity, &
                this%particles(particle)%temperature, &
                this%particles(particle)%mass
        end do
        flush(this%output_unit)
    end subroutine write_snapshot


    subroutine phase_field_names(base_name, short_base, phase, long_name, short_name)
        character(len=*), intent(in) :: base_name, short_base
        integer, intent(in) :: phase
        character(len=*), intent(out) :: long_name, short_name

        write(long_name, '(A,I2.2)') trim(base_name), phase
        write(short_name, '(A,I2.2)') trim(short_base), phase
    end subroutine phase_field_names


    pure function particle_surface_area(diameter, dimensions) result(area)
        real(dp), intent(in) :: diameter
        integer, intent(in) :: dimensions
        real(dp) :: area

        if (dimensions == 2) then
            area = pi * diameter
        else
            area = pi * diameter**2
        end if
    end function particle_surface_area


    function species_sensible_enthalpy(thermo, species, temperature) result(value)
        class(thermophysical_properties), intent(in) :: thermo
        integer, intent(in) :: species
        real(dp), intent(in) :: temperature
        real(dp) :: value, reference_enthalpy

        reference_enthalpy = thermo%specie_enthalpy_molar(T_ref, species)
        value = (thermo%specie_enthalpy_molar(temperature, species) - &
            reference_enthalpy) / thermo%molar_masses(species)
    end function species_sensible_enthalpy


    pure function lowercase(text) result(lower)
        character(len=*), intent(in) :: text
        character(len=len(text)) :: lower
        integer :: i, code

        lower = text
        do i = 1, len(text)
            code = iachar(lower(i:i))
            if (code >= iachar('A') .and. code <= iachar('Z')) then
                lower(i:i) = achar(code + iachar('a') - iachar('A'))
            end if
        end do
    end function lowercase

end module lagrangian_particles_solver_class
