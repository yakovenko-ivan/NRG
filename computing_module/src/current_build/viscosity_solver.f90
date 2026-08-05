module viscosity_solver_class

    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite

    use kind_parameters, only: dp
    use global_data, only: I_m
    use field_pointers
    use boundary_conditions_class
    use data_manager_class
    use computational_mesh_class
    use computational_domain_class
    use thermophysical_properties_class
    use chemical_properties_class
    use mpi_communications_class

    implicit none

    private
    public :: viscosity_solver, viscosity_solver_c

    type :: viscosity_solver
        private

        type(field_scalar_cons_pointer) :: energy_source
        type(field_scalar_cons_pointer) :: temperature
        type(field_scalar_cons_pointer) :: dynamic_viscosity
        type(field_scalar_cons_pointer) :: mixture_molar_mass
        type(field_vector_cons_pointer) :: velocity
        type(field_vector_cons_pointer) :: momentum_source
        type(field_vector_cons_pointer) :: mass_fraction
        type(field_vector_flow_pointer) :: stress_power_flux
        type(field_tensor_flow_pointer) :: stress

        type(field_scalar_cons), pointer :: energy_source_store => null()
        type(field_scalar_cons), pointer :: dynamic_viscosity_store => null()
        type(field_vector_cons), pointer :: momentum_source_store => null()
        type(field_vector_flow), pointer :: stress_power_flux_store => null()
        type(field_tensor_flow), pointer :: stress_store => null()

        type(computational_domain) :: domain
        type(mpi_communications) :: mpi_support
        type(boundary_conditions_pointer) :: boundary
        type(computational_mesh_pointer) :: mesh
        type(thermophysical_properties_pointer) :: thermophysics
        type(chemical_properties_pointer) :: chemistry

        real(dp), allocatable :: hoop_stress(:,:,:)
    contains
        procedure :: solve_viscosity
        procedure :: update_dynamic_viscosity
        procedure, private :: fill_viscosity_boundary_ghosts
        procedure, private :: calculate_stress_fluxes
        procedure, private :: velocity_divergence_at_cell
        procedure, private :: velocity_divergence_at_face
        procedure, private :: face_component_derivative
        procedure, private :: geometry_power
        procedure, private, nopass :: radial_face_ratio
        procedure, private, nopass :: positive_harmonic_mean
    end type viscosity_solver

    interface viscosity_solver_c
        module procedure constructor
    end interface viscosity_solver_c

contains

    type(viscosity_solver) function constructor(manager)
        type(data_manager), intent(inout) :: manager

        type(field_scalar_cons_pointer) :: scalar_pointer
        type(field_vector_cons_pointer) :: vector_pointer
        type(field_tensor_cons_pointer) :: tensor_pointer
        integer, dimension(3,2) :: cell_bounds

        constructor%domain = manager%domain
        constructor%mpi_support = manager%mpi_communications
        constructor%mesh%mesh_ptr => manager%computational_mesh_pointer%mesh_ptr
        constructor%boundary%bc_ptr => manager%boundary_conditions_pointer%bc_ptr
        constructor%thermophysics%thermo_ptr => manager%thermophysics%thermo_ptr
        constructor%chemistry%chem_ptr => manager%chemistry%chem_ptr

        call manager%get_cons_field_pointer_by_name( &
            scalar_pointer, vector_pointer, tensor_pointer, 'temperature')
        constructor%temperature%s_ptr => scalar_pointer%s_ptr
        call manager%get_cons_field_pointer_by_name( &
            scalar_pointer, vector_pointer, tensor_pointer, 'mixture_molar_mass')
        constructor%mixture_molar_mass%s_ptr => scalar_pointer%s_ptr
        call manager%get_cons_field_pointer_by_name( &
            scalar_pointer, vector_pointer, tensor_pointer, 'velocity')
        constructor%velocity%v_ptr => vector_pointer%v_ptr
        call manager%get_cons_field_pointer_by_name( &
            scalar_pointer, vector_pointer, tensor_pointer, 'specie_mass_fraction')
        constructor%mass_fraction%v_ptr => vector_pointer%v_ptr

        allocate(constructor%energy_source_store)
        allocate(constructor%dynamic_viscosity_store)
        allocate(constructor%momentum_source_store)
        allocate(constructor%stress_power_flux_store)
        allocate(constructor%stress_store)

        call manager%create_scalar_field( &
            constructor%energy_source_store, &
            'energy_production_viscosity', 'E_f_prod_visc')
        constructor%energy_source%s_ptr => constructor%energy_source_store

        call manager%create_scalar_field( &
            constructor%dynamic_viscosity_store, 'viscosity', 'nu')
        constructor%dynamic_viscosity%s_ptr => &
            constructor%dynamic_viscosity_store

        call manager%create_vector_field( &
            constructor%momentum_source_store, &
            'velocity_production_viscosity', 'v_prod_visc', 'spatial')
        constructor%momentum_source%v_ptr => constructor%momentum_source_store

        call manager%create_vector_field( &
            constructor%stress_power_flux_store, &
            'sigma_dot_v', 'sigma_dot_v', 'spatial')
        constructor%stress_power_flux%v_ptr => constructor%stress_power_flux_store

        call manager%create_tensor_field( &
            constructor%stress_store, 'stress', 'sigma')
        constructor%stress%t_ptr => constructor%stress_store

        constructor%energy_source%s_ptr%cells = 0.0_dp
        constructor%dynamic_viscosity%s_ptr%cells = 0.0_dp
        call zero_vector_cons_field(constructor%momentum_source%v_ptr)
        call zero_vector_flow_field(constructor%stress_power_flux%v_ptr)
        call zero_tensor_flow_field(constructor%stress%t_ptr)

        if (constructor%geometry_power() > 0) then
            cell_bounds = manager%domain%get_local_utter_cells_bounds()
            allocate(constructor%hoop_stress( &
                cell_bounds(1,1):cell_bounds(1,2), &
                cell_bounds(2,1):cell_bounds(2,2), &
                cell_bounds(3,1):cell_bounds(3,2)))
            constructor%hoop_stress = 0.0_dp
        end if
    end function constructor


    subroutine solve_viscosity(this, time_step)
        class(viscosity_solver), intent(inout) :: this
        real(dp), intent(in) :: time_step

        integer :: dimensions, power
        integer :: i, j, k, component, normal
        integer, dimension(3,2) :: cell_loop
        real(dp), dimension(3) :: cell_size
        real(dp) :: source_value, energy_value
        real(dp) :: left_ratio, right_ratio, radius, divergence_cell
        real(dp) :: radial_velocity_ratio

        if (.not. ieee_is_finite(time_step) .or. time_step <= 0.0_dp) then
            error stop 'Viscosity solver: time step must be finite and positive'
        end if

        dimensions = this%domain%get_domain_dimensions()
        power = this%geometry_power()
        cell_loop = this%domain%get_local_inner_cells_bounds()
        cell_size = this%mesh%mesh_ptr%get_cell_edges_length()

        this%energy_source%s_ptr%cells = 0.0_dp
        call zero_vector_cons_field(this%momentum_source%v_ptr)

        call this%update_dynamic_viscosity()
        call this%mpi_support%exchange_conservative_vector_field( &
            this%velocity%v_ptr)
        call this%calculate_stress_fluxes()

        associate( &
            energy_source => this%energy_source%s_ptr, &
            momentum_source => this%momentum_source%v_ptr, &
            viscosity => this%dynamic_viscosity%s_ptr, &
            velocity => this%velocity%v_ptr, &
            stress => this%stress%t_ptr, &
            stress_power_flux => this%stress_power_flux%v_ptr, &
            markers => this%boundary%bc_ptr%bc_markers, &
            mesh => this%mesh%mesh_ptr)

            !$omp parallel do collapse(3) schedule(static) default(shared) &
            !$omp& private(i,j,k,component,normal,source_value,energy_value, &
            !$omp& left_ratio,right_ratio,radius,divergence_cell, &
            !$omp& radial_velocity_ratio)
            do k = cell_loop(3,1), cell_loop(3,2)
                do j = cell_loop(2,1), cell_loop(2,2)
                    do i = cell_loop(1,1), cell_loop(1,2)
                        if (markers(i,j,k) /= 0) cycle

                        radius = mesh%mesh(1,i,j,k)
                        if (power > 0) then
                            if (radius <= 0.0_dp) then
                                error stop 'Viscosity solver: non-positive radial cell metric'
                            end if
                            divergence_cell = &
                                this%velocity_divergence_at_cell(i,j,k)
                            radial_velocity_ratio = velocity%pr(1)%cells(i,j,k)/radius
                            this%hoop_stress(i,j,k) = viscosity%cells(i,j,k) * &
                                (2.0_dp*radial_velocity_ratio - &
                                 (2.0_dp/3.0_dp)*divergence_cell)
                        end if

                        do component = 1, dimensions
                            source_value = 0.0_dp
                            do normal = 1, dimensions
                                left_ratio = this%radial_face_ratio( &
                                    power, normal, radius, cell_size(1), -1)
                                right_ratio = this%radial_face_ratio( &
                                    power, normal, radius, cell_size(1), 1)
                                source_value = source_value + &
                                    (right_ratio*stress%pr(component,normal)% &
                                        cells(normal, &
                                            i+I_m(normal,1), &
                                            j+I_m(normal,2), &
                                            k+I_m(normal,3)) - &
                                     left_ratio*stress%pr(component,normal)% &
                                        cells(normal,i,j,k)) / cell_size(normal)
                            end do
                            if (power > 0 .and. component == 1) then
                                source_value = source_value - &
                                    real(power,dp)*this%hoop_stress(i,j,k)/radius
                            end if
                            momentum_source%pr(component)%cells(i,j,k) = source_value
                        end do

                        energy_value = 0.0_dp
                        do normal = 1, dimensions
                            left_ratio = this%radial_face_ratio( &
                                power, normal, radius, cell_size(1), -1)
                            right_ratio = this%radial_face_ratio( &
                                power, normal, radius, cell_size(1), 1)
                            energy_value = energy_value + &
                                (right_ratio*stress_power_flux%pr(normal)% &
                                    cells(normal, &
                                        i+I_m(normal,1), &
                                        j+I_m(normal,2), &
                                        k+I_m(normal,3)) - &
                                 left_ratio*stress_power_flux%pr(normal)% &
                                    cells(normal,i,j,k)) / cell_size(normal)
                        end do
                        energy_source%cells(i,j,k) = energy_value
                    end do
                end do
            end do
            !$omp end parallel do
        end associate

        call this%mpi_support%exchange_conservative_scalar_field( &
            this%energy_source%s_ptr)
        call this%mpi_support%exchange_conservative_vector_field( &
            this%momentum_source%v_ptr)
    end subroutine solve_viscosity


    subroutine update_dynamic_viscosity(this)
        class(viscosity_solver), intent(inout) :: this

        integer :: species_number, i, j, k, specie
        integer, dimension(3,2) :: cell_loop
        real(dp), allocatable :: mass_fractions(:)
        real(dp) :: viscosity_value

        species_number = this%chemistry%chem_ptr%species_number
        cell_loop = this%domain%get_local_inner_cells_bounds()

        call this%mpi_support%exchange_conservative_scalar_field( &
            this%temperature%s_ptr)
        call this%mpi_support%exchange_conservative_scalar_field( &
            this%mixture_molar_mass%s_ptr)
        call this%mpi_support%exchange_conservative_vector_field( &
            this%mass_fraction%v_ptr)

        associate( &
            temperature => this%temperature%s_ptr, &
            viscosity => this%dynamic_viscosity%s_ptr, &
            mixture_molar_mass => this%mixture_molar_mass%s_ptr, &
            mass_fraction => this%mass_fraction%v_ptr, &
            markers => this%boundary%bc_ptr%bc_markers)

            !$omp parallel default(shared) private(i,j,k,specie,mass_fractions, &
            !$omp& viscosity_value)
            allocate(mass_fractions(species_number))
            !$omp do collapse(3) schedule(static)
            do k = cell_loop(3,1), cell_loop(3,2)
                do j = cell_loop(2,1), cell_loop(2,2)
                    do i = cell_loop(1,1), cell_loop(1,2)
                        if (markers(i,j,k) /= 0) cycle

                        do specie = 1, species_number
                            mass_fractions(specie) = &
                                mass_fraction%pr(specie)%cells(i,j,k)
                        end do
                        viscosity_value = this%thermophysics%thermo_ptr% &
                            mixture_dynamic_viscosity( &
                                temperature%cells(i,j,k), mass_fractions, &
                                mixture_molar_mass%cells(i,j,k))
                        if (.not. ieee_is_finite(viscosity_value) .or. &
                            viscosity_value <= 0.0_dp) then
                            error stop 'Viscosity solver: invalid dynamic viscosity'
                        end if
                        viscosity%cells(i,j,k) = viscosity_value
                    end do
                end do
            end do
            !$omp end do
            deallocate(mass_fractions)
            !$omp end parallel
        end associate

        call this%mpi_support%exchange_conservative_scalar_field( &
            this%dynamic_viscosity%s_ptr)
        call this%fill_viscosity_boundary_ghosts()
    end subroutine update_dynamic_viscosity


    subroutine calculate_stress_fluxes(this)
        class(viscosity_solver), intent(inout) :: this

        integer :: dimensions
        integer :: normal, component, i, j, k, dim
        integer, dimension(3,2) :: face_loop, loop
        real(dp), dimension(3) :: cell_size
        real(dp) :: viscosity_face, divergence_face
        real(dp) :: normal_derivative, transverse_derivative
        real(dp) :: stress_value, energy_flux_value

        dimensions = this%domain%get_domain_dimensions()
        face_loop = this%domain%get_local_inner_faces_bounds()
        cell_size = this%mesh%mesh_ptr%get_cell_edges_length()

        call zero_tensor_flow_field(this%stress%t_ptr)
        call zero_vector_flow_field(this%stress_power_flux%v_ptr)

        associate( &
            viscosity => this%dynamic_viscosity%s_ptr, &
            velocity => this%velocity%v_ptr, &
            stress => this%stress%t_ptr, &
            stress_power_flux => this%stress_power_flux%v_ptr)

            !$omp parallel default(shared) private(normal,component,i,j,k,dim, &
            !$omp& loop,viscosity_face,divergence_face,normal_derivative, &
            !$omp& transverse_derivative,stress_value,energy_flux_value)
            do normal = 1, dimensions
                loop = face_loop
                do dim = 1, dimensions
                    loop(dim,2) = face_loop(dim,2) - (1 - I_m(dim,normal))
                end do

                !$omp do collapse(3) schedule(static)
                do k = loop(3,1), loop(3,2)
                    do j = loop(2,1), loop(2,2)
                        do i = loop(1,1), loop(1,2)
                            viscosity_face = this%positive_harmonic_mean( &
                                viscosity%cells(i,j,k), &
                                viscosity%cells( &
                                    i-I_m(normal,1), &
                                    j-I_m(normal,2), &
                                    k-I_m(normal,3)))
                            divergence_face = &
                                this%velocity_divergence_at_face(normal,i,j,k)
                            energy_flux_value = 0.0_dp

                            do component = 1, dimensions
                                normal_derivative = &
                                    (velocity%pr(component)%cells(i,j,k) - &
                                     velocity%pr(component)%cells( &
                                        i-I_m(normal,1), &
                                        j-I_m(normal,2), &
                                        k-I_m(normal,3))) / cell_size(normal)

                                if (component == normal) then
                                    stress_value = viscosity_face * &
                                        (2.0_dp*normal_derivative - &
                                         (2.0_dp/3.0_dp)*divergence_face)
                                else
                                    transverse_derivative = &
                                        this%face_component_derivative( &
                                            normal, component, normal, i, j, k)
                                    stress_value = viscosity_face * &
                                        (normal_derivative + transverse_derivative)
                                end if

                                stress%pr(component,normal)% &
                                    cells(normal,i,j,k) = stress_value
                                energy_flux_value = energy_flux_value + &
                                    stress_value*0.5_dp*( &
                                        velocity%pr(component)%cells(i,j,k) + &
                                        velocity%pr(component)%cells( &
                                            i-I_m(normal,1), &
                                            j-I_m(normal,2), &
                                            k-I_m(normal,3)))
                            end do
                            stress_power_flux%pr(normal)% &
                                cells(normal,i,j,k) = energy_flux_value
                        end do
                    end do
                end do
                !$omp end do
            end do
            !$omp end parallel
        end associate

        do normal = 1, dimensions
            do component = 1, dimensions
                call this%mpi_support%exchange_flow_scalar_field( &
                    this%stress%t_ptr%pr(component,normal))
            end do
        end do
        call this%mpi_support%exchange_flow_vector_field( &
            this%stress_power_flux%v_ptr)
    end subroutine calculate_stress_fluxes


    real(dp) function velocity_divergence_at_cell(this, i, j, k) result(value)
        class(viscosity_solver), intent(in) :: this
        integer, intent(in) :: i, j, k

        integer :: dimensions, power, dim
        real(dp), dimension(3) :: cell_size
        real(dp) :: radius, radial_velocity_high, radial_velocity_low
        real(dp) :: radial_derivative

        dimensions = this%domain%get_domain_dimensions()
        power = this%geometry_power()
        cell_size = this%mesh%mesh_ptr%get_cell_edges_length()
        value = 0.0_dp

        if (power == 0) then
            do dim = 1, dimensions
                value = value + &
                    (this%velocity%v_ptr%pr(dim)%cells( &
                        i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3)) - &
                     this%velocity%v_ptr%pr(dim)%cells( &
                        i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3))) / &
                    (2.0_dp*cell_size(dim))
            end do
            return
        end if

        radius = this%mesh%mesh_ptr%mesh(1,i,j,k)
        radial_velocity_high = 0.5_dp*( &
            this%velocity%v_ptr%pr(1)%cells(i,j,k) + &
            this%velocity%v_ptr%pr(1)%cells(i+1,j,k))
        radial_velocity_low = 0.5_dp*( &
            this%velocity%v_ptr%pr(1)%cells(i-1,j,k) + &
            this%velocity%v_ptr%pr(1)%cells(i,j,k))
        radial_derivative = (radial_velocity_high - radial_velocity_low) / &
            cell_size(1)

        if (abs(radius) <= tiny(1.0_dp)) then
            value = real(power + 1,dp)*radial_derivative
        else
            value = ( &
                (radius + 0.5_dp*cell_size(1))**power*radial_velocity_high - &
                (radius - 0.5_dp*cell_size(1))**power*radial_velocity_low) / &
                (cell_size(1)*radius**power)
        end if

        do dim = 2, dimensions
            value = value + &
                (this%velocity%v_ptr%pr(dim)%cells( &
                    i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3)) - &
                 this%velocity%v_ptr%pr(dim)%cells( &
                    i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3))) / &
                (2.0_dp*cell_size(dim))
        end do
    end function velocity_divergence_at_cell


    real(dp) function velocity_divergence_at_face(this, normal, i, j, k) &
        result(value)
        class(viscosity_solver), intent(in) :: this
        integer, intent(in) :: normal, i, j, k

        integer :: dimensions, power, dim
        real(dp), dimension(3) :: cell_size
        real(dp) :: radius, face_radius, left_radius, right_radius
        real(dp) :: radial_velocity_left, radial_velocity_right
        real(dp) :: radial_velocity_high, radial_velocity_low

        dimensions = this%domain%get_domain_dimensions()
        power = this%geometry_power()
        cell_size = this%mesh%mesh_ptr%get_cell_edges_length()
        value = 0.0_dp

        if (power == 0) then
            do dim = 1, dimensions
                value = value + this%face_component_derivative( &
                    dim, dim, normal, i, j, k)
            end do
            return
        end if

        if (normal == 1) then
            face_radius = this%mesh%mesh_ptr%mesh(1,i,j,k) - &
                0.5_dp*cell_size(1)
            left_radius = face_radius - 0.5_dp*cell_size(1)
            right_radius = face_radius + 0.5_dp*cell_size(1)
            radial_velocity_left = this%velocity%v_ptr%pr(1)%cells(i-1,j,k)
            radial_velocity_right = this%velocity%v_ptr%pr(1)%cells(i,j,k)
            if (abs(face_radius) <= tiny(1.0_dp)) then
                value = real(power + 1,dp) * &
                    (radial_velocity_right - radial_velocity_left) / cell_size(1)
            else
                value = ( &
                    right_radius**power*radial_velocity_right - &
                    left_radius**power*radial_velocity_left) / &
                    (cell_size(1)*face_radius**power)
            end if
        else
            radius = this%mesh%mesh_ptr%mesh(1,i,j,k)
            radial_velocity_high = 0.25_dp*( &
                this%velocity%v_ptr%pr(1)%cells(i,j,k) + &
                this%velocity%v_ptr%pr(1)%cells( &
                    i-I_m(normal,1),j-I_m(normal,2),k-I_m(normal,3)) + &
                this%velocity%v_ptr%pr(1)%cells(i+1,j,k) + &
                this%velocity%v_ptr%pr(1)%cells( &
                    i+1-I_m(normal,1), &
                    j-I_m(normal,2),k-I_m(normal,3)))
            radial_velocity_low = 0.25_dp*( &
                this%velocity%v_ptr%pr(1)%cells(i,j,k) + &
                this%velocity%v_ptr%pr(1)%cells( &
                    i-I_m(normal,1),j-I_m(normal,2),k-I_m(normal,3)) + &
                this%velocity%v_ptr%pr(1)%cells(i-1,j,k) + &
                this%velocity%v_ptr%pr(1)%cells( &
                    i-1-I_m(normal,1), &
                    j-I_m(normal,2),k-I_m(normal,3)))
            if (abs(radius) <= tiny(1.0_dp)) then
                value = real(power + 1,dp) * &
                    (radial_velocity_high-radial_velocity_low)/cell_size(1)
            else
                value = ( &
                    (radius+0.5_dp*cell_size(1))**power*radial_velocity_high - &
                    (radius-0.5_dp*cell_size(1))**power*radial_velocity_low) / &
                    (cell_size(1)*radius**power)
            end if
        end if

        do dim = 2, dimensions
            value = value + this%face_component_derivative( &
                dim, dim, normal, i, j, k)
        end do
    end function velocity_divergence_at_face


    real(dp) function face_component_derivative(this, component, &
        derivative_dimension, normal, i, j, k) result(value)
        class(viscosity_solver), intent(in) :: this
        integer, intent(in) :: component, derivative_dimension, normal
        integer, intent(in) :: i, j, k

        integer :: left_i, left_j, left_k
        integer :: high_i, high_j, high_k, low_i, low_j, low_k
        real(dp), dimension(3) :: cell_size
        real(dp) :: high_value, low_value

        cell_size = this%mesh%mesh_ptr%get_cell_edges_length()
        left_i = i - I_m(normal,1)
        left_j = j - I_m(normal,2)
        left_k = k - I_m(normal,3)

        if (derivative_dimension == normal) then
            value = ( &
                this%velocity%v_ptr%pr(component)%cells(i,j,k) - &
                this%velocity%v_ptr%pr(component)%cells( &
                    left_i,left_j,left_k)) / cell_size(normal)
            return
        end if

        high_i = i + I_m(derivative_dimension,1)
        high_j = j + I_m(derivative_dimension,2)
        high_k = k + I_m(derivative_dimension,3)
        low_i = i - I_m(derivative_dimension,1)
        low_j = j - I_m(derivative_dimension,2)
        low_k = k - I_m(derivative_dimension,3)

        high_value = 0.25_dp*( &
            this%velocity%v_ptr%pr(component)%cells(i,j,k) + &
            this%velocity%v_ptr%pr(component)%cells(left_i,left_j,left_k) + &
            this%velocity%v_ptr%pr(component)%cells(high_i,high_j,high_k) + &
            this%velocity%v_ptr%pr(component)%cells( &
                high_i-I_m(normal,1), &
                high_j-I_m(normal,2),high_k-I_m(normal,3)))
        low_value = 0.25_dp*( &
            this%velocity%v_ptr%pr(component)%cells(i,j,k) + &
            this%velocity%v_ptr%pr(component)%cells(left_i,left_j,left_k) + &
            this%velocity%v_ptr%pr(component)%cells(low_i,low_j,low_k) + &
            this%velocity%v_ptr%pr(component)%cells( &
                low_i-I_m(normal,1), &
                low_j-I_m(normal,2),low_k-I_m(normal,3)))
        value = (high_value-low_value)/cell_size(derivative_dimension)
    end function face_component_derivative


    subroutine fill_viscosity_boundary_ghosts(this)
        class(viscosity_solver), intent(inout) :: this

        integer :: dimensions, i, j, k, dim, side, sign
        integer :: ghost_i, ghost_j, ghost_k
        integer, dimension(3,2) :: cell_loop

        dimensions = this%domain%get_domain_dimensions()
        cell_loop = this%domain%get_local_inner_cells_bounds()

        associate( &
            viscosity => this%dynamic_viscosity%s_ptr, &
            markers => this%boundary%bc_ptr%bc_markers)
            do k = cell_loop(3,1), cell_loop(3,2)
                do j = cell_loop(2,1), cell_loop(2,2)
                    do i = cell_loop(1,1), cell_loop(1,2)
                        if (markers(i,j,k) /= 0) cycle
                        do dim = 1, dimensions
                            do side = 1, 2
                                sign = (-1)**side
                                ghost_i = i + sign*I_m(dim,1)
                                ghost_j = j + sign*I_m(dim,2)
                                ghost_k = k + sign*I_m(dim,3)
                                if (markers(ghost_i,ghost_j,ghost_k) /= 0) then
                                    viscosity%cells(ghost_i,ghost_j,ghost_k) = &
                                        viscosity%cells(i,j,k)
                                end if
                            end do
                        end do
                    end do
                end do
            end do
        end associate
    end subroutine fill_viscosity_boundary_ghosts


    integer function geometry_power(this) result(power)
        class(viscosity_solver), intent(in) :: this

        select case(trim(this%domain%get_coordinate_system_name()))
        case('cylindrical')
            power = 1
        case('spherical')
            power = 2
        case default
            power = 0
        end select
    end function geometry_power


    real(dp) function radial_face_ratio(power, dimension, radius, &
        radial_spacing, side) result(ratio)
        integer, intent(in) :: power, dimension, side
        real(dp), intent(in) :: radius, radial_spacing
        real(dp) :: face_radius

        ratio = 1.0_dp
        if (power == 0 .or. dimension /= 1) return
        if (radius <= 0.0_dp) then
            error stop 'Viscosity solver: non-positive radial cell metric'
        end if
        face_radius = radius + 0.5_dp*real(side,dp)*radial_spacing
        if (face_radius < 0.0_dp) then
            error stop 'Viscosity solver: negative radial face coordinate'
        end if
        ratio = face_radius**power/radius**power
    end function radial_face_ratio


    real(dp) function positive_harmonic_mean(left_value, &
        right_value) result(mean_value)
        real(dp), intent(in) :: left_value, right_value

        if (left_value <= 0.0_dp .or. right_value <= 0.0_dp) then
            error stop 'Viscosity solver: non-positive dynamic viscosity'
        end if
        mean_value = 2.0_dp*left_value*right_value/(left_value+right_value)
    end function positive_harmonic_mean


    subroutine zero_vector_cons_field(field)
        type(field_vector_cons), pointer, intent(inout) :: field
        integer :: component

        do component = 1, size(field%pr)
            field%pr(component)%cells = 0.0_dp
        end do
    end subroutine zero_vector_cons_field


    subroutine zero_vector_flow_field(field)
        type(field_vector_flow), pointer, intent(inout) :: field
        integer :: component

        do component = 1, size(field%pr)
            field%pr(component)%cells = 0.0_dp
        end do
    end subroutine zero_vector_flow_field


    subroutine zero_tensor_flow_field(field)
        type(field_tensor_flow), pointer, intent(inout) :: field
        integer :: component, normal

        do normal = 1, size(field%pr,2)
            do component = 1, size(field%pr,1)
                field%pr(component,normal)%cells = 0.0_dp
            end do
        end do
    end subroutine zero_tensor_flow_field

end module viscosity_solver_class
