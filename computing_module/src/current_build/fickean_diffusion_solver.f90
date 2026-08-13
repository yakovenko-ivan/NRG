module fickean_diffusion_solver_class

    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite

    use kind_parameters, only: dp
    use global_data, only: I_m, T_ref
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
    public :: diffusion_solver, diffusion_solver_c

    real(dp), parameter :: default_soret_alpha_h = 0.895_dp
    real(dp), parameter :: default_soret_alpha_h2 = 0.910_dp

    type :: diffusion_solver
        private

        type(field_scalar_cons_pointer) :: temperature
        type(field_scalar_cons_pointer) :: pressure
        type(field_scalar_cons_pointer) :: mixture_molar_mass
        type(field_scalar_cons_pointer) :: density
        type(field_scalar_cons_pointer) :: energy_source
        type(field_vector_cons_pointer) :: mass_fraction
        type(field_vector_cons_pointer) :: species_source
        type(field_vector_cons_pointer) :: diffusivity
        type(field_vector_cons_pointer) :: thermal_diffusion_coefficient

        type(field_scalar_cons), pointer :: energy_source_store => null()
        type(field_vector_cons), pointer :: species_source_store => null()
        type(field_vector_cons), pointer :: diffusivity_store => null()
        type(field_vector_cons), pointer :: thermal_diffusion_store => null()

        type(computational_domain) :: domain
        type(mpi_communications) :: mpi_support
        type(boundary_conditions_pointer) :: boundary
        type(computational_mesh_pointer) :: mesh
        type(thermophysical_properties_pointer) :: thermophysics
        type(chemical_properties_pointer) :: chemistry

        logical :: soret_enabled = .true.
        integer, allocatable :: soret_species_indices(:)
        real(dp), allocatable :: soret_alpha(:)
        real(dp), allocatable :: reference_enthalpy_molar(:)
    contains
        procedure :: solve_diffusion
        procedure :: update_diffusion_coefficients
        procedure :: set_soret_enabled
        procedure :: configure_reduced_soret
        procedure, private :: initialize_default_soret_model
        procedure, private :: fill_coefficient_boundary_ghosts
        procedure, private :: calculate_face_flux
        procedure, private :: geometry_power
        procedure, private, nopass :: radial_face_ratio
        procedure, private, nopass :: positive_harmonic_mean
    end type diffusion_solver

    interface diffusion_solver_c
        module procedure constructor
    end interface diffusion_solver_c

contains

    type(diffusion_solver) function constructor(manager, soret_enabled)
        type(data_manager), intent(inout) :: manager
        logical, intent(in), optional :: soret_enabled

        type(field_scalar_cons_pointer) :: scalar_pointer
        type(field_vector_cons_pointer) :: vector_pointer
        type(field_tensor_cons_pointer) :: tensor_pointer
        integer :: specie, species_number

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
            scalar_pointer, vector_pointer, tensor_pointer, 'pressure')
        constructor%pressure%s_ptr => scalar_pointer%s_ptr
        call manager%get_cons_field_pointer_by_name( &
            scalar_pointer, vector_pointer, tensor_pointer, 'density')
        constructor%density%s_ptr => scalar_pointer%s_ptr
        call manager%get_cons_field_pointer_by_name( &
            scalar_pointer, vector_pointer, tensor_pointer, 'mixture_molar_mass')
        constructor%mixture_molar_mass%s_ptr => scalar_pointer%s_ptr
        call manager%get_cons_field_pointer_by_name( &
            scalar_pointer, vector_pointer, tensor_pointer, 'specie_mass_fraction')
        constructor%mass_fraction%v_ptr => vector_pointer%v_ptr

        allocate(constructor%energy_source_store)
        allocate(constructor%species_source_store)
        allocate(constructor%diffusivity_store)
        allocate(constructor%thermal_diffusion_store)

        call manager%create_scalar_field( &
            constructor%energy_source_store, &
            'energy_production_diffusion', 'E_f_prod_diff')
        constructor%energy_source%s_ptr => constructor%energy_source_store

        call manager%create_vector_field( &
            constructor%species_source_store, &
            'specie_production_diffusion', 'Y_prod_diff', 'chemical')
        constructor%species_source%v_ptr => constructor%species_source_store

        call manager%create_vector_field( &
            constructor%diffusivity_store, 'diffusivity', 'D', 'chemical')
        constructor%diffusivity%v_ptr => constructor%diffusivity_store

        call manager%create_vector_field( &
            constructor%thermal_diffusion_store, &
            'thermal_diffusion_coefficient', 'D_T', 'chemical')
        constructor%thermal_diffusion_coefficient%v_ptr => &
            constructor%thermal_diffusion_store

        constructor%energy_source%s_ptr%cells = 0.0_dp
        call zero_vector_field(constructor%species_source%v_ptr)
        call zero_vector_field(constructor%diffusivity%v_ptr)
        call zero_vector_field(constructor%thermal_diffusion_coefficient%v_ptr)

        species_number = constructor%chemistry%chem_ptr%species_number
        allocate(constructor%reference_enthalpy_molar(species_number))
        do specie = 1, species_number
            constructor%reference_enthalpy_molar(specie) = &
                constructor%thermophysics%thermo_ptr%specie_enthalpy_molar( &
                    T_ref, specie)
        end do

        call constructor%initialize_default_soret_model()
        if (present(soret_enabled)) constructor%soret_enabled = soret_enabled
    end function constructor


    subroutine initialize_default_soret_model(this)
        class(diffusion_solver), intent(inout) :: this

        allocate(this%soret_species_indices(2), this%soret_alpha(2))
        this%soret_species_indices(1) = &
            this%chemistry%chem_ptr%get_chemical_specie_index('H')
        this%soret_species_indices(2) = &
            this%chemistry%chem_ptr%get_chemical_specie_index('H2')
        this%soret_alpha = (/default_soret_alpha_h, default_soret_alpha_h2/)
    end subroutine initialize_default_soret_model


    subroutine set_soret_enabled(this, enabled)
        class(diffusion_solver), intent(inout) :: this
        logical, intent(in) :: enabled

        this%soret_enabled = enabled
    end subroutine set_soret_enabled


    subroutine configure_reduced_soret(this, species_indices, alpha)
        class(diffusion_solver), intent(inout) :: this
        integer, dimension(:), intent(in) :: species_indices
        real(dp), dimension(:), intent(in) :: alpha

        integer :: target

        if (size(species_indices) /= size(alpha)) then
            error stop 'Diffusion solver: Soret configuration size mismatch'
        end if
        if (size(species_indices) == 0) then
            error stop 'Diffusion solver: empty Soret target list'
        end if
        do target = 1, size(alpha)
            if (.not. ieee_is_finite(alpha(target))) then
                error stop 'Diffusion solver: non-finite Soret coefficient'
            end if
            if (species_indices(target) < 0 .or. &
                species_indices(target) > this%chemistry%chem_ptr%species_number) then
                error stop 'Diffusion solver: invalid Soret species index'
            end if
        end do

        if (allocated(this%soret_species_indices)) &
            deallocate(this%soret_species_indices)
        if (allocated(this%soret_alpha)) deallocate(this%soret_alpha)
        allocate(this%soret_species_indices(size(species_indices)))
        allocate(this%soret_alpha(size(alpha)))
        this%soret_species_indices = species_indices
        this%soret_alpha = alpha
    end subroutine configure_reduced_soret


    subroutine solve_diffusion(this, time_step)
        class(diffusion_solver), intent(inout) :: this
        real(dp), intent(in) :: time_step

        integer :: dimensions, species_number, power
        integer :: i, j, k, dim, specie
        integer, dimension(3,2) :: cell_loop
        real(dp), dimension(3) :: cell_size
        real(dp), allocatable :: flux_left(:), flux_right(:)
        real(dp), allocatable :: face_mass_fraction(:), mass_fraction_jump(:)
        real(dp), allocatable :: raw_mass_flux(:)
        real(dp) :: enthalpy_flux_left, enthalpy_flux_right
        real(dp) :: left_ratio, right_ratio, radius

        if (.not. ieee_is_finite(time_step) .or. time_step <= 0.0_dp) then
            error stop 'Diffusion solver: time step must be finite and positive'
        end if

        call this%update_diffusion_coefficients()

        dimensions = this%domain%get_domain_dimensions()
        species_number = this%chemistry%chem_ptr%species_number
        power = this%geometry_power()
        cell_loop = this%domain%get_local_inner_cells_bounds()
        cell_size = this%mesh%mesh_ptr%get_cell_edges_length()

        this%energy_source%s_ptr%cells = 0.0_dp
        call zero_vector_field(this%species_source%v_ptr)

        associate( &
            energy_source => this%energy_source%s_ptr, &
            species_source => this%species_source%v_ptr, &
            markers => this%boundary%bc_ptr%bc_markers, &
            mesh => this%mesh%mesh_ptr)

            !$omp parallel default(shared) private(i,j,k,dim,specie,flux_left, &
            !$omp& flux_right,face_mass_fraction,mass_fraction_jump, &
            !$omp& raw_mass_flux, &
            !$omp& enthalpy_flux_left,enthalpy_flux_right, &
            !$omp& left_ratio,right_ratio,radius)
            allocate(flux_left(species_number), flux_right(species_number), &
                face_mass_fraction(species_number), &
                mass_fraction_jump(species_number), raw_mass_flux(species_number))
            !$omp do collapse(3) schedule(static)
            do k = cell_loop(3,1), cell_loop(3,2)
                do j = cell_loop(2,1), cell_loop(2,2)
                    do i = cell_loop(1,1), cell_loop(1,2)
                        if (markers(i,j,k) /= 0) cycle

                        radius = mesh%mesh(1,i,j,k)
                        do dim = 1, dimensions
                            call this%calculate_face_flux( &
                                i, j, k, dim, -1, cell_size(dim), &
                                this%reference_enthalpy_molar, face_mass_fraction, &
                                mass_fraction_jump, raw_mass_flux, flux_left, &
                                enthalpy_flux_left)
                            call this%calculate_face_flux( &
                                i, j, k, dim, 1, cell_size(dim), &
                                this%reference_enthalpy_molar, face_mass_fraction, &
                                mass_fraction_jump, raw_mass_flux, flux_right, &
                                enthalpy_flux_right)

                            left_ratio = this%radial_face_ratio( &
                                power, dim, radius, cell_size(1), -1)
                            right_ratio = this%radial_face_ratio( &
                                power, dim, radius, cell_size(1), 1)

                            do specie = 1, species_number
                                species_source%pr(specie)%cells(i,j,k) = &
                                    species_source%pr(specie)%cells(i,j,k) + &
                                    (right_ratio*flux_right(specie) - &
                                     left_ratio*flux_left(specie)) / cell_size(dim)
                            end do
                            energy_source%cells(i,j,k) = &
                                energy_source%cells(i,j,k) + &
                                (right_ratio*enthalpy_flux_right - &
                                 left_ratio*enthalpy_flux_left) / cell_size(dim)
                        end do
                    end do
                end do
            end do
            !$omp end do
            deallocate(flux_left, flux_right, face_mass_fraction, &
                mass_fraction_jump, raw_mass_flux)
            !$omp end parallel
        end associate

        call this%mpi_support%exchange_conservative_scalar_field( &
            this%energy_source%s_ptr)
        call this%mpi_support%exchange_conservative_vector_field( &
            this%species_source%v_ptr)
    end subroutine solve_diffusion


    subroutine update_diffusion_coefficients(this)
        class(diffusion_solver), intent(inout) :: this

        integer :: species_number, i, j, k, specie
        integer, dimension(3,2) :: cell_loop
        real(dp), allocatable :: mass_fractions(:), diffusion_values(:)
        real(dp), allocatable :: thermal_diffusion_values(:)

        species_number = this%chemistry%chem_ptr%species_number
        cell_loop = this%domain%get_local_inner_cells_bounds()

        call this%mpi_support%exchange_conservative_scalar_field( &
            this%temperature%s_ptr)
        call this%mpi_support%exchange_conservative_scalar_field( &
            this%pressure%s_ptr)
        call this%mpi_support%exchange_conservative_scalar_field( &
            this%mixture_molar_mass%s_ptr)
        call this%mpi_support%exchange_conservative_vector_field( &
            this%mass_fraction%v_ptr)

        associate( &
            temperature => this%temperature%s_ptr, &
            pressure => this%pressure%s_ptr, &
            mixture_molar_mass => this%mixture_molar_mass%s_ptr, &
            mass_fraction => this%mass_fraction%v_ptr, &
            diffusivity => this%diffusivity%v_ptr, &
            thermal_diffusion => this%thermal_diffusion_coefficient%v_ptr, &
            markers => this%boundary%bc_ptr%bc_markers)

            !$omp parallel default(shared) private(i,j,k,specie,mass_fractions, &
            !$omp& diffusion_values,thermal_diffusion_values)
            allocate(mass_fractions(species_number))
            allocate(diffusion_values(species_number))
            allocate(thermal_diffusion_values(species_number))
            !$omp do collapse(3) schedule(static)
            do k = cell_loop(3,1), cell_loop(3,2)
                do j = cell_loop(2,1), cell_loop(2,2)
                    do i = cell_loop(1,1), cell_loop(1,2)
                        if (markers(i,j,k) /= 0) cycle

                        do specie = 1, species_number
                            mass_fractions(specie) = &
                                mass_fraction%pr(specie)%cells(i,j,k)
                        end do

                        call this%thermophysics%thermo_ptr% &
                            mixture_averaged_diffusion_coefficients( &
                                temperature%cells(i,j,k), pressure%cells(i,j,k), &
                                mass_fractions, &
                                mixture_molar_mass%cells(i,j,k), diffusion_values)

                        thermal_diffusion_values = 0.0_dp
                        if (this%soret_enabled) then
                            call this%thermophysics%thermo_ptr% &
                                reduced_thermal_diffusion_coefficients( &
                                    temperature%cells(i,j,k), mass_fractions, &
                                    mixture_molar_mass%cells(i,j,k), &
                                    this%soret_species_indices, this%soret_alpha, &
                                    thermal_diffusion_values)
                        end if

                        do specie = 1, species_number
                            if (.not. ieee_is_finite(diffusion_values(specie)) .or. &
                                diffusion_values(specie) < 0.0_dp) then
                                error stop 'Diffusion solver: invalid diffusivity'
                            end if
                            if (.not. ieee_is_finite( &
                                thermal_diffusion_values(specie))) then
                                error stop 'Diffusion solver: invalid Soret coefficient'
                            end if
                            diffusivity%pr(specie)%cells(i,j,k) = &
                                diffusion_values(specie)
                            thermal_diffusion%pr(specie)%cells(i,j,k) = &
                                thermal_diffusion_values(specie)
                        end do
                    end do
                end do
            end do
            !$omp end do
            deallocate(mass_fractions, diffusion_values, &
                thermal_diffusion_values)
            !$omp end parallel
        end associate

        call this%mpi_support%exchange_conservative_vector_field( &
            this%diffusivity%v_ptr)
        call this%mpi_support%exchange_conservative_vector_field( &
            this%thermal_diffusion_coefficient%v_ptr)
        call this%fill_coefficient_boundary_ghosts()
    end subroutine update_diffusion_coefficients


    subroutine calculate_face_flux(this, i, j, k, dimension, side, spacing, &
        reference_enthalpy_molar, face_mass_fraction, mass_fraction_jump, &
        raw_mass_flux, mass_flux, enthalpy_flux)
        class(diffusion_solver), intent(in) :: this
        integer, intent(in) :: i, j, k, dimension, side
        real(dp), intent(in) :: spacing
        real(dp), dimension(:), intent(in) :: reference_enthalpy_molar
        real(dp), dimension(:), intent(inout) :: face_mass_fraction
        real(dp), dimension(:), intent(inout) :: mass_fraction_jump, raw_mass_flux
        real(dp), dimension(:), intent(out) :: mass_flux
        real(dp), intent(out) :: enthalpy_flux

        integer :: neighbour_i, neighbour_j, neighbour_k
        integer :: species_number, specie
        real(dp) :: density_face, temperature_gradient
        real(dp) :: mass_fraction_sum, raw_flux_sum
        real(dp) :: diffusivity_face, thermal_diffusion_face
        real(dp) :: gradient_mass_fraction, raw_flux
        real(dp) :: sensible_enthalpy_cell, sensible_enthalpy_neighbour
        real(dp) :: temperature_cell, temperature_neighbour, side_real
        real(dp) :: mass_fraction_cell, mass_fraction_neighbour
        real(dp) :: species_molar_mass

        species_number = size(mass_flux)
        mass_flux = 0.0_dp
        enthalpy_flux = 0.0_dp

        neighbour_i = i + side*I_m(dimension,1)
        neighbour_j = j + side*I_m(dimension,2)
        neighbour_k = k + side*I_m(dimension,3)
        if (this%boundary%bc_ptr%bc_markers( &
            neighbour_i,neighbour_j,neighbour_k) /= 0) return

        density_face = 0.5_dp*( &
            this%density%s_ptr%cells(i,j,k) + &
            this%density%s_ptr%cells(neighbour_i,neighbour_j,neighbour_k))
        if (.not. ieee_is_finite(density_face) .or. density_face <= 0.0_dp) then
            error stop 'Diffusion solver: invalid face density'
        end if

        side_real = real(side,dp)
        temperature_cell = this%temperature%s_ptr%cells(i,j,k)
        temperature_neighbour = &
            this%temperature%s_ptr%cells(neighbour_i,neighbour_j,neighbour_k)
        temperature_gradient = side_real * &
            (log(max(temperature_neighbour, tiny(1.0_dp))) - &
             log(max(temperature_cell, tiny(1.0_dp)))) / spacing

        mass_fraction_sum = 0.0_dp
        do specie = 1, species_number
            mass_fraction_cell = &
                this%mass_fraction%v_ptr%pr(specie)%cells(i,j,k)
            mass_fraction_neighbour = &
                this%mass_fraction%v_ptr%pr(specie)%cells( &
                    neighbour_i,neighbour_j,neighbour_k)
            face_mass_fraction(specie) = &
                0.5_dp*(mass_fraction_cell + mass_fraction_neighbour)
            mass_fraction_jump(specie) = &
                mass_fraction_neighbour - mass_fraction_cell
            if (face_mass_fraction(specie) < -1.0e-12_dp) then
                error stop 'Diffusion solver: negative face mass fraction'
            end if
            face_mass_fraction(specie) = max(face_mass_fraction(specie), 0.0_dp)
            mass_fraction_sum = mass_fraction_sum + face_mass_fraction(specie)
        end do
        if (mass_fraction_sum <= tiny(1.0_dp)) then
            error stop 'Diffusion solver: zero face composition sum'
        end if
        face_mass_fraction = face_mass_fraction/mass_fraction_sum

        raw_flux_sum = 0.0_dp
        do specie = 1, species_number
            diffusivity_face = this%positive_harmonic_mean( &
                this%diffusivity%v_ptr%pr(specie)%cells(i,j,k), &
                this%diffusivity%v_ptr%pr(specie)%cells( &
                    neighbour_i,neighbour_j,neighbour_k))
            thermal_diffusion_face = 0.5_dp*( &
                this%thermal_diffusion_coefficient%v_ptr%pr(specie)% &
                    cells(i,j,k) + &
                this%thermal_diffusion_coefficient%v_ptr%pr(specie)% &
                    cells(neighbour_i,neighbour_j,neighbour_k))
            gradient_mass_fraction = &
                side_real*mass_fraction_jump(specie) / spacing

            raw_flux = density_face*diffusivity_face*gradient_mass_fraction + &
                thermal_diffusion_face*temperature_gradient
            raw_mass_flux(specie) = raw_flux
            raw_flux_sum = raw_flux_sum + raw_flux
        end do

        mass_flux = raw_mass_flux - face_mass_fraction*raw_flux_sum
        mass_flux(species_number) = mass_flux(species_number) - sum(mass_flux)

        do specie = 1, species_number
            species_molar_mass = &
                this%thermophysics%thermo_ptr%molar_masses(specie)
            sensible_enthalpy_cell = &
                (this%thermophysics%thermo_ptr%specie_enthalpy_molar( &
                    temperature_cell, specie) - &
                 reference_enthalpy_molar(specie)) / species_molar_mass
            sensible_enthalpy_neighbour = &
                (this%thermophysics%thermo_ptr%specie_enthalpy_molar( &
                    temperature_neighbour, specie) - &
                 reference_enthalpy_molar(specie)) / species_molar_mass
            enthalpy_flux = enthalpy_flux + mass_flux(specie) * &
                0.5_dp*(sensible_enthalpy_cell + sensible_enthalpy_neighbour)
        end do

    end subroutine calculate_face_flux


    subroutine fill_coefficient_boundary_ghosts(this)
        class(diffusion_solver), intent(inout) :: this

        integer :: dimensions, species_number
        integer :: i, j, k, dim, side, sign, specie
        integer :: ghost_i, ghost_j, ghost_k
        integer, dimension(3,2) :: cell_loop

        dimensions = this%domain%get_domain_dimensions()
        species_number = this%chemistry%chem_ptr%species_number
        cell_loop = this%domain%get_local_inner_cells_bounds()

        associate( &
            diffusivity => this%diffusivity%v_ptr, &
            thermal_diffusion => this%thermal_diffusion_coefficient%v_ptr, &
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
                                if (markers(ghost_i,ghost_j,ghost_k) == 0) cycle
                                do specie = 1, species_number
                                    diffusivity%pr(specie)%cells( &
                                        ghost_i,ghost_j,ghost_k) = &
                                        diffusivity%pr(specie)%cells(i,j,k)
                                    thermal_diffusion%pr(specie)%cells( &
                                        ghost_i,ghost_j,ghost_k) = &
                                        thermal_diffusion%pr(specie)%cells(i,j,k)
                                end do
                            end do
                        end do
                    end do
                end do
            end do
        end associate
    end subroutine fill_coefficient_boundary_ghosts


    integer function geometry_power(this) result(power)
        class(diffusion_solver), intent(in) :: this

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
            error stop 'Diffusion solver: non-positive radial cell metric'
        end if
        face_radius = radius + 0.5_dp*real(side,dp)*radial_spacing
        if (face_radius < 0.0_dp) then
            error stop 'Diffusion solver: negative radial face coordinate'
        end if
        ratio = face_radius**power/radius**power
    end function radial_face_ratio


    real(dp) function positive_harmonic_mean(left_value, &
        right_value) result(mean_value)
        real(dp), intent(in) :: left_value, right_value

        if (left_value < 0.0_dp .or. right_value < 0.0_dp) then
            error stop 'Diffusion solver: negative diffusion coefficient'
        end if
        if (left_value <= tiny(1.0_dp) .or. right_value <= tiny(1.0_dp)) then
            mean_value = 0.0_dp
        else
            mean_value = 2.0_dp*left_value*right_value / &
                (left_value + right_value)
        end if
    end function positive_harmonic_mean


    subroutine zero_vector_field(field)
        type(field_vector_cons), pointer, intent(inout) :: field

        integer :: component

        do component = 1, size(field%pr)
            field%pr(component)%cells = 0.0_dp
        end do
    end subroutine zero_vector_field

end module fickean_diffusion_solver_class
