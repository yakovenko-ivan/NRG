module continuous_particles_solver_class

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
    use mpi_communications_class
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite

    implicit none

    private
    public :: continuous_particles_solver, continuous_particles_solver_c

    real(dp), parameter :: phase_density_floor = 1.0e-30_dp
    real(dp), parameter :: deformation_floor = 5.0e-2_dp

    type :: continuous_particles_solver
        private

        type(field_scalar_cons_pointer) :: gas_temperature
        type(field_scalar_cons_pointer) :: gas_density
        type(field_scalar_cons_pointer) :: gas_viscosity
        type(field_scalar_cons_pointer) :: gas_conductivity
        type(field_vector_cons_pointer) :: gas_velocity

        type(field_scalar_cons_pointer) :: particle_temperature
        type(field_scalar_cons_pointer) :: particle_temperature_stage
        type(field_scalar_cons_pointer) :: particle_density
        type(field_scalar_cons_pointer) :: particle_mass
        type(field_scalar_cons_pointer) :: particle_number_density
        type(field_scalar_cons_pointer) :: surface_area_density
        type(field_vector_cons_pointer) :: particle_velocity
        type(field_vector_cons_pointer) :: particle_velocity_stage
        type(field_scalar_flow_pointer) :: mass_flux
        type(field_scalar_flow_pointer) :: number_flux

        type(field_scalar_cons_pointer) :: gas_mass_source
        type(field_scalar_cons_pointer) :: gas_energy_source
        type(field_vector_cons_pointer) :: gas_momentum_source
        type(field_vector_cons_pointer) :: gas_species_source

        ! Each solver instance owns its registered fields.  This replaces the
        ! old module-global arrays and removes pre_constructor ordering rules.
        type(field_scalar_cons), pointer :: particle_temperature_store => null()
        type(field_scalar_cons), pointer :: particle_temperature_stage_store => null()
        type(field_scalar_cons), pointer :: particle_density_store => null()
        type(field_scalar_cons), pointer :: particle_mass_store => null()
        type(field_scalar_cons), pointer :: particle_number_density_store => null()
        type(field_scalar_cons), pointer :: surface_area_density_store => null()
        type(field_vector_cons), pointer :: particle_velocity_store => null()
        type(field_vector_cons), pointer :: particle_velocity_stage_store => null()
        type(field_scalar_flow), pointer :: mass_flux_store => null()
        type(field_scalar_flow), pointer :: number_flux_store => null()
        type(field_scalar_cons), pointer :: gas_mass_source_store => null()
        type(field_scalar_cons), pointer :: gas_energy_source_store => null()
        type(field_vector_cons), pointer :: gas_momentum_source_store => null()
        type(field_vector_cons), pointer :: gas_species_source_store => null()

        type(computational_domain) :: domain
        type(mpi_communications) :: mpi_support
        type(boundary_conditions_pointer) :: boundary
        type(computational_mesh_pointer) :: mesh
        type(chemical_properties_pointer) :: chemistry
        type(thermophysical_properties_pointer) :: thermophysics

        type(particles_phase) :: parameters
        real(dp), dimension(3) :: gravity = 0.0_dp
        integer :: phase_number = 0
        integer :: material_index = 0
    contains
        procedure :: set_initial_distributions
        procedure :: advance

        ! Compatibility bindings for legacy callers.  New gas solvers should
        ! call advance(), which owns the complete phase update sequence.
        procedure :: particles_euler_step_v_E => calculate_exchange_sources
        procedure :: particles_lagrange_step => calculate_transport_fluxes
        procedure :: particles_final_step => conservative_transport_update
        procedure :: apply_boundary_conditions_main
        procedure :: apply_boundary_conditions_interm_v_p

        procedure, private :: calculate_exchange_sources
        procedure, private :: calculate_transport_fluxes
        procedure, private :: conservative_transport_update
        procedure, private :: reset_coupling_sources
        procedure, private :: particle_mass_from_diameter
        procedure, private :: particle_diameter_from_mass
        procedure, private :: drag_relaxation_rate
        procedure, private :: cell_metric
        procedure, private :: face_metric
        procedure, private :: geometry_power
        procedure, private :: injection_factor
    end type continuous_particles_solver

    interface continuous_particles_solver_c
        module procedure constructor
    end interface continuous_particles_solver_c

contains

    type(continuous_particles_solver) function constructor(manager, particles, phase_number)
        type(data_manager), intent(inout) :: manager
        type(particles_phase), intent(in) :: particles
        integer, intent(in) :: phase_number

        type(field_scalar_cons_pointer) :: scalar_pointer
        type(field_vector_cons_pointer) :: vector_pointer
        type(field_tensor_cons_pointer) :: tensor_pointer
        character(len=48) :: long_name, short_name

        constructor%parameters = particles
        constructor%phase_number = phase_number
        constructor%domain = manager%domain
        constructor%mpi_support = manager%mpi_communications
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
        call manager%get_cons_field_pointer_by_name( &
            scalar_pointer, vector_pointer, tensor_pointer, 'velocity')
        constructor%gas_velocity%v_ptr => vector_pointer%v_ptr

        allocate(constructor%particle_temperature_store)
        allocate(constructor%particle_temperature_stage_store)
        allocate(constructor%particle_density_store)
        allocate(constructor%particle_mass_store)
        allocate(constructor%particle_number_density_store)
        allocate(constructor%surface_area_density_store)
        allocate(constructor%particle_velocity_store)
        allocate(constructor%particle_velocity_stage_store)
        allocate(constructor%mass_flux_store)
        allocate(constructor%number_flux_store)
        allocate(constructor%gas_mass_source_store)
        allocate(constructor%gas_energy_source_store)
        allocate(constructor%gas_momentum_source_store)
        allocate(constructor%gas_species_source_store)

        call phase_field_names('temperature_particles', 'T_d', phase_number, &
            long_name, short_name)
        call manager%create_scalar_field( &
            constructor%particle_temperature_store, long_name, short_name)
        constructor%particle_temperature%s_ptr => &
            constructor%particle_temperature_store

        call phase_field_names('temperature_particles_interm', 'T_d_int', &
            phase_number, long_name, short_name)
        call manager%create_scalar_field( &
            constructor%particle_temperature_stage_store, long_name, short_name)
        constructor%particle_temperature_stage%s_ptr => &
            constructor%particle_temperature_stage_store

        call phase_field_names('density_particles', 'rho_d', phase_number, &
            long_name, short_name)
        call manager%create_scalar_field( &
            constructor%particle_density_store, long_name, short_name)
        constructor%particle_density%s_ptr => constructor%particle_density_store

        call phase_field_names('mass_particles', 'mass_d', phase_number, &
            long_name, short_name)
        call manager%create_scalar_field( &
            constructor%particle_mass_store, long_name, short_name)
        constructor%particle_mass%s_ptr => constructor%particle_mass_store

        call phase_field_names('number_density_particles', 'numdens_d', &
            phase_number, long_name, short_name)
        call manager%create_scalar_field( &
            constructor%particle_number_density_store, long_name, short_name)
        constructor%particle_number_density%s_ptr => &
            constructor%particle_number_density_store

        call phase_field_names('Esurf', 'Esurf', &
            phase_number, long_name, short_name)
        call manager%create_scalar_field( &
            constructor%surface_area_density_store, long_name, short_name)
        constructor%surface_area_density%s_ptr => &
            constructor%surface_area_density_store

        call phase_field_names('velocity_particles', 'v_p', phase_number, &
            long_name, short_name)
        call manager%create_vector_field( &
            constructor%particle_velocity_store, long_name, short_name, 'spatial')
        constructor%particle_velocity%v_ptr => constructor%particle_velocity_store

        call phase_field_names('velocity_particles_interm', 'v_p_int', &
            phase_number, long_name, short_name)
        call manager%create_vector_field( &
            constructor%particle_velocity_stage_store, long_name, short_name, &
            'spatial')
        constructor%particle_velocity_stage%v_ptr => &
            constructor%particle_velocity_stage_store

        call phase_field_names('mass_flux_particles', 'm_flux_d', phase_number, &
            long_name, short_name)
        call manager%create_scalar_field( &
            constructor%mass_flux_store, long_name, short_name)
        constructor%mass_flux%s_ptr => constructor%mass_flux_store

        call phase_field_names('number_density_flux_particles', 'n_flux_d', &
            phase_number, long_name, short_name)
        call manager%create_scalar_field( &
            constructor%number_flux_store, long_name, short_name)
        constructor%number_flux%s_ptr => constructor%number_flux_store

        call phase_field_names('density_production_particles', 'rho_prod_d', &
            phase_number, long_name, short_name)
        call manager%create_scalar_field( &
            constructor%gas_mass_source_store, long_name, short_name)
        constructor%gas_mass_source%s_ptr => constructor%gas_mass_source_store

        call phase_field_names('energy_production_particles', 'E_prod_d', &
            phase_number, long_name, short_name)
        call manager%create_scalar_field( &
            constructor%gas_energy_source_store, long_name, short_name)
        constructor%gas_energy_source%s_ptr => constructor%gas_energy_source_store

        call phase_field_names('velocity_production_particles', 'v_prod_d', &
            phase_number, long_name, short_name)
        call manager%create_vector_field( &
            constructor%gas_momentum_source_store, long_name, short_name, &
            'spatial')
        constructor%gas_momentum_source%v_ptr => &
            constructor%gas_momentum_source_store

        call phase_field_names('concentration_production_particles', 'Y_prod_d', &
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
            error stop 'continuous particles: material species was not found'
        end if
    end function constructor


    subroutine set_initial_distributions(this)
        class(continuous_particles_solver), intent(inout) :: this

        integer, dimension(3,2) :: loop
        integer :: dimensions, i, j, k, dim
        real(dp) :: nominal_mass, density_value
        character(len=24) :: initialization

        dimensions = this%domain%get_domain_dimensions()
        loop = this%domain%get_local_inner_cells_bounds()
        nominal_mass = this%particle_mass_from_diameter(this%parameters%diameter)
        initialization = lowercase(trim(this%parameters%initial_condition))

        associate(rho_p => this%particle_density%s_ptr, &
                  mass_p => this%particle_mass%s_ptr, &
                  number_p => this%particle_number_density%s_ptr, &
                  temperature_p => this%particle_temperature%s_ptr, &
                  temperature_stage => this%particle_temperature_stage%s_ptr, &
                  velocity_p => this%particle_velocity%v_ptr, &
                  velocity_stage => this%particle_velocity_stage%v_ptr, &
                  area_density => this%surface_area_density%s_ptr, &
                  bc => this%boundary%bc_ptr)
        do k = loop(3,1), loop(3,2)
            do j = loop(2,1), loop(2,2)
                do i = loop(1,1), loop(1,2)
                    if (bc%bc_markers(i,j,k) /= 0) cycle

                    select case (initialization)
                    case ('none')
                        density_value = 0.0_dp
                        mass_p%cells(i,j,k) = nominal_mass
                        number_p%cells(i,j,k) = 0.0_dp
                        temperature_p%cells(i,j,k) = &
                            this%parameters%initial_temperature
                        do dim = 1, dimensions
                            velocity_p%pr(dim)%cells(i,j,k) = 0.0_dp
                        end do
                    case ('uniform')
                        density_value = this%parameters%initial_mass_loading
                        mass_p%cells(i,j,k) = nominal_mass
                        number_p%cells(i,j,k) = density_value / nominal_mass
                        temperature_p%cells(i,j,k) = &
                            this%parameters%initial_temperature
                        do dim = 1, dimensions
                            velocity_p%pr(dim)%cells(i,j,k) = &
                                this%parameters%initial_velocity(dim)
                        end do
                    case ('field')
                        density_value = max(rho_p%cells(i,j,k), 0.0_dp)
                        if (mass_p%cells(i,j,k) <= 0.0_dp) then
                            mass_p%cells(i,j,k) = nominal_mass
                        end if
                        if (number_p%cells(i,j,k) <= 0.0_dp) then
                            number_p%cells(i,j,k) = density_value / &
                                mass_p%cells(i,j,k)
                        end if
                        if (temperature_p%cells(i,j,k) <= 0.0_dp) then
                            temperature_p%cells(i,j,k) = &
                                this%parameters%initial_temperature
                        end if
                    case default
                        error stop 'continuous particles: invalid initialization mode'
                    end select

                    rho_p%cells(i,j,k) = density_value
                    temperature_stage%cells(i,j,k) = temperature_p%cells(i,j,k)
                    do dim = 1, dimensions
                        velocity_stage%pr(dim)%cells(i,j,k) = &
                            velocity_p%pr(dim)%cells(i,j,k)
                    end do
                    area_density%cells(i,j,k) = particle_surface_area( &
                        this%parameters%diameter, number_p%cells(i,j,k), dimensions)
                end do
            end do
        end do
        end associate

        call this%apply_boundary_conditions_main(0.0_dp)
        call this%reset_coupling_sources()
    end subroutine set_initial_distributions


    subroutine advance(this, time_step, time)
        class(continuous_particles_solver), intent(inout) :: this
        real(dp), intent(in) :: time_step
        real(dp), intent(in) :: time

        if (.not.ieee_is_finite(time_step) .or. time_step <= 0.0_dp) then
            error stop 'continuous particles: non-positive time step'
        end if

        call this%apply_boundary_conditions_main(time)
        call this%calculate_exchange_sources(time_step)
        call this%apply_boundary_conditions_interm_v_p()
        call this%calculate_transport_fluxes(time_step)
        call this%conservative_transport_update(time_step)
        call this%apply_boundary_conditions_main(time + time_step)
    end subroutine advance


    subroutine calculate_exchange_sources(this, time_step)
        class(continuous_particles_solver), intent(inout) :: this
        real(dp), intent(in) :: time_step

        integer, dimension(3,2) :: loop
        integer :: dimensions, i, j, k, dim
        real(dp) :: rho_g, rho_p_old, mass_old, mass_new, diameter
        real(dp) :: relative_speed, beta, reynolds, nusselt
        real(dp) :: heat_rate_per_mass, heat_source_rate
        real(dp) :: temperature_trial, boiling_temperature
        real(dp) :: evaporation_fraction, mass_source, vapor_enthalpy
        real(dp) :: kinetic_energy
        real(dp) :: mechanical_source, particle_velocity_mid
        real(dp), dimension(3) :: relative_velocity, drag_acceleration

        dimensions = this%domain%get_domain_dimensions()
        loop = this%domain%get_local_inner_cells_bounds()
        boiling_temperature = this%parameters%material_boiling_temperature

        call this%reset_coupling_sources()

        associate(Tg => this%gas_temperature%s_ptr, &
                  rho_gas => this%gas_density%s_ptr, &
                  mu => this%gas_viscosity%s_ptr, &
                  conductivity => this%gas_conductivity%s_ptr, &
                  vg => this%gas_velocity%v_ptr, &
                  Tp => this%particle_temperature%s_ptr, &
                  Tp_stage => this%particle_temperature_stage%s_ptr, &
                  rho_p => this%particle_density%s_ptr, &
                  mass_p => this%particle_mass%s_ptr, &
                  number_p => this%particle_number_density%s_ptr, &
                  vp => this%particle_velocity%v_ptr, &
                  vp_stage => this%particle_velocity_stage%v_ptr, &
                  area_density => this%surface_area_density%s_ptr, &
                  mass_src => this%gas_mass_source%s_ptr, &
                  energy_src => this%gas_energy_source%s_ptr, &
                  momentum_src => this%gas_momentum_source%v_ptr, &
                  species_src => this%gas_species_source%v_ptr, &
                  bc => this%boundary%bc_ptr)
        !$omp parallel do collapse(3) schedule(static) &
        !$omp& private(i,j,k,dim,rho_g,rho_p_old,mass_old,mass_new,diameter, &
        !$omp& relative_speed,beta,reynolds,nusselt,heat_rate_per_mass, &
        !$omp& heat_source_rate, &
        !$omp& temperature_trial,evaporation_fraction,mass_source, &
        !$omp& vapor_enthalpy,kinetic_energy,mechanical_source, &
        !$omp& particle_velocity_mid, &
        !$omp& relative_velocity,drag_acceleration)
        do k = loop(3,1), loop(3,2)
            do j = loop(2,1), loop(2,2)
                do i = loop(1,1), loop(1,2)
                    if (bc%bc_markers(i,j,k) /= 0) cycle

                    rho_g = rho_gas%cells(i,j,k)
                    rho_p_old = max(rho_p%cells(i,j,k), 0.0_dp)
                    mass_old = max(mass_p%cells(i,j,k), &
                        this%particle_mass_from_diameter(this%parameters%diameter) * &
                        1.0e-12_dp)

                    if (rho_p_old <= phase_density_floor) then
                        Tp_stage%cells(i,j,k) = Tp%cells(i,j,k)
                        do dim = 1, dimensions
                            vp_stage%pr(dim)%cells(i,j,k) = &
                                vp%pr(dim)%cells(i,j,k)
                        end do
                        cycle
                    end if
                    if (rho_g <= 0.0_dp .or. mu%cells(i,j,k) <= 0.0_dp) then
                        error stop 'continuous particles: invalid gas transport state'
                    end if

                    diameter = this%particle_diameter_from_mass(mass_old)
                    relative_speed = 0.0_dp
                    do dim = 1, dimensions
                        relative_velocity(dim) = &
                            vg%pr(dim)%cells(i,j,k) - vp%pr(dim)%cells(i,j,k)
                        relative_speed = relative_speed + relative_velocity(dim)**2
                    end do
                    relative_speed = sqrt(relative_speed)
                    beta = this%drag_relaxation_rate( &
                        rho_g, mu%cells(i,j,k), diameter, mass_old, relative_speed)

                    mechanical_source = 0.0_dp
                    do dim = 1, dimensions
                        drag_acceleration(dim) = beta * relative_velocity(dim)
                        vp_stage%pr(dim)%cells(i,j,k) = &
                            vp%pr(dim)%cells(i,j,k) + time_step * &
                            (drag_acceleration(dim) + this%gravity(dim))

                        if (this%parameters%inertial) then
                            momentum_src%pr(dim)%cells(i,j,k) = &
                                -rho_p_old * drag_acceleration(dim) / rho_g
                            particle_velocity_mid = 0.5_dp * ( &
                                vp%pr(dim)%cells(i,j,k) + &
                                vp_stage%pr(dim)%cells(i,j,k))
                            mechanical_source = mechanical_source - rho_p_old * &
                                drag_acceleration(dim) * particle_velocity_mid
                        else
                            momentum_src%pr(dim)%cells(i,j,k) = 0.0_dp
                        end if
                    end do
                    energy_src%cells(i,j,k) = mechanical_source

                    Tp_stage%cells(i,j,k) = Tp%cells(i,j,k)
                    evaporation_fraction = 0.0_dp
                    if (this%parameters%heating) then
                        reynolds = rho_g * relative_speed * diameter / &
                            mu%cells(i,j,k)
                        nusselt = 2.0_dp + 0.6_dp * sqrt(max(reynolds, 0.0_dp)) * &
                            0.7_dp**(1.0_dp / 3.0_dp)
                        heat_rate_per_mass = nusselt * conductivity%cells(i,j,k) * &
                            particle_surface_area(diameter, 1.0_dp, dimensions) * &
                            (Tg%cells(i,j,k) - Tp%cells(i,j,k)) / &
                            (diameter * mass_old)
                        heat_rate_per_mass = heat_rate_per_mass / &
                            this%parameters%material_heat_capacity
                        heat_source_rate = rho_p_old * &
                            this%parameters%material_heat_capacity * &
                            heat_rate_per_mass
                        energy_src%cells(i,j,k) = energy_src%cells(i,j,k) - &
                            heat_source_rate
                        temperature_trial = Tp%cells(i,j,k) + &
                            time_step * heat_rate_per_mass

                        if (this%parameters%evaporating .and. &
                            temperature_trial > boiling_temperature) then
                            evaporation_fraction = min(1.0_dp, &
                                (temperature_trial - boiling_temperature) * &
                                this%parameters%material_heat_capacity / &
                                this%parameters%material_latent_heat)
                            Tp_stage%cells(i,j,k) = boiling_temperature
                        else
                            Tp_stage%cells(i,j,k) = max(1.0_dp, temperature_trial)
                        end if
                    end if

                    evaporation_fraction = min(max(evaporation_fraction, 0.0_dp), &
                        1.0_dp - 100.0_dp * epsilon(1.0_dp))
                    mass_new = mass_old * (1.0_dp - evaporation_fraction)
                    mass_source = rho_p_old * evaporation_fraction / time_step

                    if (mass_source > 0.0_dp) then
                        mass_p%cells(i,j,k) = mass_new
                        rho_p%cells(i,j,k) = rho_p_old * &
                            (1.0_dp - evaporation_fraction)
                        mass_src%cells(i,j,k) = mass_source
                        species_src%pr(this%material_index)%cells(i,j,k) = &
                            mass_source

                        vapor_enthalpy = species_sensible_enthalpy( &
                            this%thermophysics%thermo_ptr, this%material_index, &
                            Tp_stage%cells(i,j,k))
                        kinetic_energy = 0.0_dp
                        do dim = 1, dimensions
                            kinetic_energy = kinetic_energy + 0.5_dp * &
                                vp_stage%pr(dim)%cells(i,j,k)**2
                        end do
                        energy_src%cells(i,j,k) = energy_src%cells(i,j,k) + &
                            mass_source * (vapor_enthalpy + kinetic_energy)
                        do dim = 1, dimensions
                            momentum_src%pr(dim)%cells(i,j,k) = &
                                momentum_src%pr(dim)%cells(i,j,k) + &
                                mass_source * vp_stage%pr(dim)%cells(i,j,k) / rho_g
                        end do
                    end if

                    number_p%cells(i,j,k) = rho_p%cells(i,j,k) / &
                        max(mass_p%cells(i,j,k), tiny(1.0_dp))
                    area_density%cells(i,j,k) = particle_surface_area( &
                        this%particle_diameter_from_mass(mass_p%cells(i,j,k)), &
                        number_p%cells(i,j,k), dimensions)
                end do
            end do
        end do
        !$omp end parallel do
        end associate
    end subroutine calculate_exchange_sources


    subroutine calculate_transport_fluxes(this, time_step)
        class(continuous_particles_solver), intent(inout) :: this
        real(dp), intent(in) :: time_step

        integer, dimension(3,2) :: face_loop
        real(dp), dimension(3) :: dx, base_area
        integer :: dimensions, power
        integer :: i, j, k, dim, il, jl, kl
        real(dp) :: velocity_face, deformation, area
        real(dp) :: density_upwind, number_upwind

        dimensions = this%domain%get_domain_dimensions()
        power = this%geometry_power()
        face_loop = this%domain%get_local_inner_faces_bounds()
        dx = this%mesh%mesh_ptr%get_cell_edges_length()
        base_area = this%mesh%mesh_ptr%get_cell_surface_area()

        call this%mpi_support%exchange_conservative_scalar_field( &
            this%particle_density%s_ptr)
        call this%mpi_support%exchange_conservative_scalar_field( &
            this%particle_number_density%s_ptr)
        call this%mpi_support%exchange_conservative_vector_field( &
            this%particle_velocity_stage%v_ptr)

        associate(rho_p => this%particle_density%s_ptr, &
                  number_p => this%particle_number_density%s_ptr, &
                  vp => this%particle_velocity_stage%v_ptr, &
                  mass_flow => this%mass_flux%s_ptr, &
                  number_flow => this%number_flux%s_ptr)
        !$omp parallel do collapse(3) schedule(static) &
        !$omp& private(i,j,k,dim,il,jl,kl,velocity_face,deformation,area, &
        !$omp& density_upwind,number_upwind)
        do k = face_loop(3,1), face_loop(3,2)
            do j = face_loop(2,1), face_loop(2,2)
                do i = face_loop(1,1), face_loop(1,2)
                    do dim = 1, dimensions
                        il = i - I_m(dim,1)
                        jl = j - I_m(dim,2)
                        kl = k - I_m(dim,3)

                        velocity_face = 0.5_dp * ( &
                            vp%pr(dim)%cells(il,jl,kl) + &
                            vp%pr(dim)%cells(i,j,k))
                        deformation = 1.0_dp + time_step * ( &
                            vp%pr(dim)%cells(i,j,k) - &
                            vp%pr(dim)%cells(il,jl,kl)) / dx(dim)
                        if (deformation <= deformation_floor) then
                            error stop 'continuous particles: non-positive CPM deformation'
                        end if

                        if (velocity_face >= 0.0_dp) then
                            density_upwind = max(rho_p%cells(il,jl,kl), 0.0_dp)
                            number_upwind = max(number_p%cells(il,jl,kl), 0.0_dp)
                        else
                            density_upwind = max(rho_p%cells(i,j,k), 0.0_dp)
                            number_upwind = max(number_p%cells(i,j,k), 0.0_dp)
                        end if

                        area = base_area(dim) * this%face_metric( &
                            dim, i, j, k, power)
                        mass_flow%cells(dim,i,j,k) = density_upwind * &
                            velocity_face * area * time_step / deformation
                        number_flow%cells(dim,i,j,k) = number_upwind * &
                            velocity_face * area * time_step / deformation
                    end do
                end do
            end do
        end do
        !$omp end parallel do
        end associate
    end subroutine calculate_transport_fluxes


    subroutine conservative_transport_update(this, time_step)
        class(continuous_particles_solver), intent(inout) :: this
        real(dp), intent(in) :: time_step

        integer, dimension(3,2) :: loop
        real(dp), dimension(3) :: base_area
        integer :: dimensions, power
        integer :: i, j, k, dim, il, jl, kl, ir, jr, kr
        real(dp) :: base_volume, volume, rho_old, rho_new, number_new
        real(dp) :: mass_left, mass_right, number_left, number_right
        real(dp) :: temperature_extensive, momentum_extensive
        real(dp) :: left_value, right_value, nominal_mass

        ! time_step is part of the stored face flux; retain the argument for
        ! interface symmetry and explicit validity checking.
        if (time_step <= 0.0_dp) then
            error stop 'continuous particles: invalid transport time step'
        end if

        dimensions = this%domain%get_domain_dimensions()
        power = this%geometry_power()
        loop = this%domain%get_local_inner_cells_bounds()
        base_volume = this%mesh%mesh_ptr%get_cell_volume()
        base_area = this%mesh%mesh_ptr%get_cell_surface_area()
        nominal_mass = this%particle_mass_from_diameter(this%parameters%diameter)

        call this%mpi_support%exchange_conservative_scalar_field( &
            this%particle_temperature_stage%s_ptr)
        call this%mpi_support%exchange_conservative_vector_field( &
            this%particle_velocity_stage%v_ptr)
        call this%mpi_support%exchange_flow_scalar_field(this%mass_flux%s_ptr)
        call this%mpi_support%exchange_flow_scalar_field(this%number_flux%s_ptr)

        associate(rho_p => this%particle_density%s_ptr, &
                  mass_p => this%particle_mass%s_ptr, &
                  number_p => this%particle_number_density%s_ptr, &
                  Tp => this%particle_temperature%s_ptr, &
                  Tp_stage => this%particle_temperature_stage%s_ptr, &
                  vp => this%particle_velocity%v_ptr, &
                  vp_stage => this%particle_velocity_stage%v_ptr, &
                  area_density => this%surface_area_density%s_ptr, &
                  mass_flow => this%mass_flux%s_ptr, &
                  number_flow => this%number_flux%s_ptr, &
                  bc => this%boundary%bc_ptr)
        !$omp parallel do collapse(3) schedule(static) &
        !$omp& private(i,j,k,dim,il,jl,kl,ir,jr,kr,volume,rho_old,rho_new, &
        !$omp& number_new,mass_left,mass_right,number_left,number_right, &
        !$omp& temperature_extensive,momentum_extensive,left_value,right_value)
        do k = loop(3,1), loop(3,2)
            do j = loop(2,1), loop(2,2)
                do i = loop(1,1), loop(1,2)
                    if (bc%bc_markers(i,j,k) /= 0) cycle

                    volume = base_volume * this%cell_metric(i, j, k, power)
                    rho_old = max(rho_p%cells(i,j,k), 0.0_dp)
                    rho_new = rho_old
                    number_new = max(number_p%cells(i,j,k), 0.0_dp)
                    temperature_extensive = rho_old * Tp_stage%cells(i,j,k)

                    do dim = 1, dimensions
                        il = i - I_m(dim,1)
                        jl = j - I_m(dim,2)
                        kl = k - I_m(dim,3)
                        ir = i + I_m(dim,1)
                        jr = j + I_m(dim,2)
                        kr = k + I_m(dim,3)

                        mass_left = mass_flow%cells(dim,i,j,k)
                        mass_right = mass_flow%cells(dim,ir,jr,kr)
                        number_left = number_flow%cells(dim,i,j,k)
                        number_right = number_flow%cells(dim,ir,jr,kr)

                        rho_new = rho_new + (mass_left - mass_right) / volume
                        number_new = number_new + &
                            (number_left - number_right) / volume

                        if (mass_left >= 0.0_dp) then
                            left_value = Tp_stage%cells(il,jl,kl)
                        else
                            left_value = Tp_stage%cells(i,j,k)
                        end if
                        if (mass_right >= 0.0_dp) then
                            right_value = Tp_stage%cells(i,j,k)
                        else
                            right_value = Tp_stage%cells(ir,jr,kr)
                        end if
                        temperature_extensive = temperature_extensive + &
                            (mass_left * left_value - mass_right * right_value) / volume
                    end do

                    if (rho_new < -1000.0_dp * epsilon(1.0_dp) * &
                        max(rho_old, 1.0_dp)) then
                        error stop 'continuous particles: negative phase density'
                    end if
                    rho_new = max(rho_new, 0.0_dp)
                    number_new = max(number_new, 0.0_dp)

                    do dim = 1, dimensions
                        momentum_extensive = rho_old * &
                            vp_stage%pr(dim)%cells(i,j,k)
                        do il = 1, dimensions
                            jl = i - I_m(il,1)
                            kl = j - I_m(il,2)
                            ir = k - I_m(il,3)
                            ! Reuse scalar temporaries as left indices below.
                            mass_left = mass_flow%cells(il,i,j,k)
                            mass_right = mass_flow%cells(il, &
                                i+I_m(il,1), j+I_m(il,2), k+I_m(il,3))

                            if (mass_left >= 0.0_dp) then
                                left_value = vp_stage%pr(dim)%cells(jl,kl,ir)
                            else
                                left_value = vp_stage%pr(dim)%cells(i,j,k)
                            end if
                            if (mass_right >= 0.0_dp) then
                                right_value = vp_stage%pr(dim)%cells(i,j,k)
                            else
                                right_value = vp_stage%pr(dim)%cells( &
                                    i+I_m(il,1), j+I_m(il,2), k+I_m(il,3))
                            end if
                            momentum_extensive = momentum_extensive + &
                                (mass_left * left_value - mass_right * right_value) / &
                                volume
                        end do

                        if (rho_new > phase_density_floor) then
                            vp%pr(dim)%cells(i,j,k) = momentum_extensive / rho_new
                        else
                            vp%pr(dim)%cells(i,j,k) = 0.0_dp
                        end if
                    end do

                    rho_p%cells(i,j,k) = rho_new
                    number_p%cells(i,j,k) = number_new
                    if (rho_new > phase_density_floor .and. number_new > &
                        phase_density_floor) then
                        Tp%cells(i,j,k) = temperature_extensive / rho_new
                        mass_p%cells(i,j,k) = rho_new / number_new
                    else
                        Tp%cells(i,j,k) = this%parameters%initial_temperature
                        mass_p%cells(i,j,k) = nominal_mass
                    end if
                    area_density%cells(i,j,k) = particle_surface_area( &
                        this%particle_diameter_from_mass(mass_p%cells(i,j,k)), &
                        number_new, dimensions)
                end do
            end do
        end do
        !$omp end parallel do
        end associate
    end subroutine conservative_transport_update


    subroutine apply_boundary_conditions_main(this, time)
        class(continuous_particles_solver), intent(inout) :: this
        real(dp), intent(in) :: time

        integer, dimension(3,2) :: loop
        integer :: dimensions, species_number
        integer :: i, j, k, dim, side, sign, boundary_number, component
        integer :: ig, jg, kg
        character(len=24) :: boundary_name, interaction
        real(dp) :: loading, ramp, nominal_mass

        dimensions = this%domain%get_domain_dimensions()
        species_number = this%chemistry%chem_ptr%species_number
        loop = this%domain%get_local_inner_cells_bounds()
        interaction = lowercase(trim(this%parameters%boundary_interaction))
        nominal_mass = this%particle_mass_from_diameter(this%parameters%diameter)
        ramp = this%injection_factor(time)

        associate(rho_p => this%particle_density%s_ptr, &
                  mass_p => this%particle_mass%s_ptr, &
                  number_p => this%particle_number_density%s_ptr, &
                  Tp => this%particle_temperature%s_ptr, &
                  vp => this%particle_velocity%v_ptr, &
                  area_density => this%surface_area_density%s_ptr, &
                  bc => this%boundary%bc_ptr)
        !$omp parallel do collapse(3) schedule(static) &
        !$omp& private(i,j,k,dim,side,sign,boundary_number,component,ig,jg,kg, &
        !$omp& boundary_name,loading)
        do k = loop(3,1), loop(3,2)
            do j = loop(2,1), loop(2,2)
                do i = loop(1,1), loop(1,2)
                    if (bc%bc_markers(i,j,k) /= 0) cycle
                    do dim = 1, dimensions
                        do side = 1, 2
                            sign = (-1)**side
                            ig = i + sign * I_m(dim,1)
                            jg = j + sign * I_m(dim,2)
                            kg = k + sign * I_m(dim,3)
                            boundary_number = bc%bc_markers(ig,jg,kg)
                            if (boundary_number == 0) cycle

                            boundary_name = lowercase(trim( &
                                bc%boundary_types(boundary_number)%get_type_name()))
                            rho_p%cells(ig,jg,kg) = rho_p%cells(i,j,k)
                            mass_p%cells(ig,jg,kg) = mass_p%cells(i,j,k)
                            number_p%cells(ig,jg,kg) = number_p%cells(i,j,k)
                            Tp%cells(ig,jg,kg) = Tp%cells(i,j,k)
                            area_density%cells(ig,jg,kg) = &
                                area_density%cells(i,j,k)
                            do component = 1, dimensions
                                vp%pr(component)%cells(ig,jg,kg) = &
                                    vp%pr(component)%cells(i,j,k)
                            end do

                            select case (boundary_name)
                            case ('wall', 'symmetry_plane')
                                select case (interaction)
                                case ('reflect')
                                    vp%pr(dim)%cells(ig,jg,kg) = &
                                        -vp%pr(dim)%cells(i,j,k)
                                case ('stick')
                                    do component = 1, dimensions
                                        vp%pr(component)%cells(ig,jg,kg) = &
                                            -vp%pr(component)%cells(i,j,k)
                                    end do
                                case ('escape')
                                    rho_p%cells(ig,jg,kg) = 0.0_dp
                                    number_p%cells(ig,jg,kg) = 0.0_dp
                                    area_density%cells(ig,jg,kg) = 0.0_dp
                                end select
                            case ('inlet')
                                if (this%parameters%inlet_mass_loading >= 0.0_dp) then
                                    loading = ramp * &
                                        this%parameters%inlet_mass_loading
                                    rho_p%cells(ig,jg,kg) = loading
                                    mass_p%cells(ig,jg,kg) = nominal_mass
                                    number_p%cells(ig,jg,kg) = loading / nominal_mass
                                    Tp%cells(ig,jg,kg) = &
                                        this%parameters%inlet_temperature
                                    do component = 1, dimensions
                                        vp%pr(component)%cells(ig,jg,kg) = &
                                            this%parameters%inlet_velocity(component)
                                    end do
                                    area_density%cells(ig,jg,kg) = &
                                        particle_surface_area( &
                                            this%parameters%diameter, &
                                            number_p%cells(ig,jg,kg), dimensions)
                                end if
                            case ('outlet')
                                ! Zero-gradient extrapolation already applied.
                            case default
                                error stop 'continuous particles: unsupported boundary type'
                            end select
                        end do
                    end do
                end do
            end do
        end do
        !$omp end parallel do
        end associate
    end subroutine apply_boundary_conditions_main


    subroutine apply_boundary_conditions_interm_v_p(this)
        class(continuous_particles_solver), intent(inout) :: this

        integer, dimension(3,2) :: loop
        integer :: dimensions, i, j, k, dim, side, sign, boundary_number
        integer :: component, ig, jg, kg
        character(len=24) :: boundary_name, interaction

        dimensions = this%domain%get_domain_dimensions()
        loop = this%domain%get_local_inner_cells_bounds()
        interaction = lowercase(trim(this%parameters%boundary_interaction))

        associate(vp => this%particle_velocity_stage%v_ptr, &
                  bc => this%boundary%bc_ptr)
        !$omp parallel do collapse(3) schedule(static) &
        !$omp& private(i,j,k,dim,side,sign,boundary_number,component,ig,jg,kg, &
        !$omp& boundary_name)
        do k = loop(3,1), loop(3,2)
            do j = loop(2,1), loop(2,2)
                do i = loop(1,1), loop(1,2)
                    if (bc%bc_markers(i,j,k) /= 0) cycle
                    do dim = 1, dimensions
                        do side = 1, 2
                            sign = (-1)**side
                            ig = i + sign * I_m(dim,1)
                            jg = j + sign * I_m(dim,2)
                            kg = k + sign * I_m(dim,3)
                            boundary_number = bc%bc_markers(ig,jg,kg)
                            if (boundary_number == 0) cycle

                            boundary_name = lowercase(trim( &
                                bc%boundary_types(boundary_number)%get_type_name()))
                            do component = 1, dimensions
                                vp%pr(component)%cells(ig,jg,kg) = &
                                    vp%pr(component)%cells(i,j,k)
                            end do
                            if (boundary_name == 'wall' .or. &
                                boundary_name == 'symmetry_plane') then
                                if (interaction == 'reflect') then
                                    vp%pr(dim)%cells(ig,jg,kg) = &
                                        -vp%pr(dim)%cells(i,j,k)
                                else if (interaction == 'stick') then
                                    do component = 1, dimensions
                                        vp%pr(component)%cells(ig,jg,kg) = &
                                            -vp%pr(component)%cells(i,j,k)
                                    end do
                                end if
                            end if
                        end do
                    end do
                end do
            end do
        end do
        !$omp end parallel do
        end associate
    end subroutine apply_boundary_conditions_interm_v_p


    subroutine reset_coupling_sources(this)
        class(continuous_particles_solver), intent(inout) :: this
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


    pure function particle_mass_from_diameter(this, diameter) result(mass)
        class(continuous_particles_solver), intent(in) :: this
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


    pure function particle_diameter_from_mass(this, mass) result(diameter)
        class(continuous_particles_solver), intent(in) :: this
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


    pure function drag_relaxation_rate( &
        this, gas_density, viscosity, diameter, mass, relative_speed) result(beta)
        class(continuous_particles_solver), intent(in) :: this
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


    function geometry_power(this) result(power)
        class(continuous_particles_solver), intent(in) :: this
        integer :: power
        character(len=20) :: coordinate_system

        coordinate_system = trim(this%domain%get_coordinate_system_name())
        select case (coordinate_system)
        case ('cartesian')
            power = 0
        case ('cylindrical')
            power = 1
        case ('spherical')
            power = 2
        case default
            error stop 'continuous particles: unknown coordinate system'
        end select
    end function geometry_power


    pure function cell_metric(this, i, j, k, power) result(metric)
        class(continuous_particles_solver), intent(in) :: this
        integer, intent(in) :: i, j, k, power
        real(dp) :: metric, radius

        if (power == 0) then
            metric = 1.0_dp
        else
            radius = this%mesh%mesh_ptr%mesh(1,i,j,k)
            metric = max(radius, 0.0_dp)**power
        end if
    end function cell_metric


    pure function face_metric(this, dim, i, j, k, power) result(metric)
        class(continuous_particles_solver), intent(in) :: this
        integer, intent(in) :: dim, i, j, k, power
        real(dp) :: metric, radius
        real(dp), dimension(3) :: dx

        if (power == 0 .or. dim /= 1) then
            metric = this%cell_metric(i, j, k, power)
            if (dim /= 1) metric = this%cell_metric(i, j, k, power)
        else
            dx = this%mesh%mesh_ptr%get_cell_edges_length()
            radius = this%mesh%mesh_ptr%mesh(1,i,j,k) - 0.5_dp * dx(1)
            metric = max(radius, 0.0_dp)**power
        end if
    end function face_metric


    pure function injection_factor(this, time) result(factor)
        class(continuous_particles_solver), intent(in) :: this
        real(dp), intent(in) :: time
        real(dp) :: factor

        if (time <= this%parameters%injection_start_time) then
            factor = 0.0_dp
        else if (this%parameters%injection_ramp_time <= 0.0_dp) then
            factor = 1.0_dp
        else
            factor = min(1.0_dp, (time - this%parameters%injection_start_time) / &
                this%parameters%injection_ramp_time)
        end if
    end function injection_factor


    subroutine phase_field_names(base_name, short_base, phase, long_name, short_name)
        character(len=*), intent(in) :: base_name, short_base
        integer, intent(in) :: phase
        character(len=*), intent(out) :: long_name, short_name

        write(long_name, '(A,I2.2)') trim(base_name), phase
        write(short_name, '(A,I2.2)') trim(short_base), phase
    end subroutine phase_field_names


    pure function particle_surface_area(diameter, number_density, dimensions) &
        result(area_density)
        real(dp), intent(in) :: diameter, number_density
        integer, intent(in) :: dimensions
        real(dp) :: area_density

        if (dimensions == 2) then
            area_density = pi * diameter * max(number_density, 0.0_dp)
        else
            area_density = pi * diameter**2 * max(number_density, 0.0_dp)
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

end module continuous_particles_solver_class
