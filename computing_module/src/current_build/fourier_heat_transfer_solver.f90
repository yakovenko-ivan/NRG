module fourier_heat_transfer_solver_class

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
    public :: heat_transfer_solver, heat_transfer_solver_c

    type :: heat_transfer_solver
        private

        type(field_scalar_cons_pointer) :: energy_source
        type(field_scalar_cons_pointer) :: temperature
        type(field_scalar_cons_pointer) :: thermal_conductivity
        type(field_scalar_cons_pointer) :: mixture_molar_mass
        type(field_vector_cons_pointer) :: mass_fraction

        type(field_scalar_cons), pointer :: energy_source_store => null()
        type(field_scalar_cons), pointer :: thermal_conductivity_store => null()

        type(computational_domain) :: domain
        type(mpi_communications) :: mpi_support
        type(boundary_conditions_pointer) :: boundary
        type(computational_mesh_pointer) :: mesh
        type(thermophysical_properties_pointer) :: thermophysics
        type(chemical_properties_pointer) :: chemistry
    contains
        procedure :: solve_heat_transfer
        procedure :: update_thermal_conductivity
        procedure, private :: fill_conductivity_boundary_ghosts
        procedure, private :: face_conductivity
        procedure, private :: geometry_power
        procedure, private, nopass :: radial_face_ratio
    end type heat_transfer_solver

    interface heat_transfer_solver_c
        module procedure constructor
    end interface heat_transfer_solver_c

contains

    type(heat_transfer_solver) function constructor(manager)
        type(data_manager), intent(inout) :: manager

        type(field_scalar_cons_pointer) :: scalar_pointer
        type(field_vector_cons_pointer) :: vector_pointer
        type(field_tensor_cons_pointer) :: tensor_pointer

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
            scalar_pointer, vector_pointer, tensor_pointer, 'specie_mass_fraction')
        constructor%mass_fraction%v_ptr => vector_pointer%v_ptr

        allocate(constructor%energy_source_store)
        allocate(constructor%thermal_conductivity_store)

        call manager%create_scalar_field( &
            constructor%energy_source_store, &
            'energy_production_heat_transfer', 'E_f_prod_heat')
        constructor%energy_source%s_ptr => constructor%energy_source_store

        call manager%create_scalar_field( &
            constructor%thermal_conductivity_store, &
            'thermal_conductivity', 'kappa')
        constructor%thermal_conductivity%s_ptr => &
            constructor%thermal_conductivity_store

        constructor%energy_source%s_ptr%cells = 0.0_dp
        constructor%thermal_conductivity%s_ptr%cells = 0.0_dp
    end function constructor


    subroutine solve_heat_transfer(this, time_step)
        class(heat_transfer_solver), intent(inout) :: this
        real(dp), intent(in) :: time_step

        integer :: dimensions, power
        integer :: i, j, k, dim
        integer :: left_i, left_j, left_k, right_i, right_j, right_k
        integer :: left_boundary, right_boundary
        integer, dimension(3,2) :: cell_loop
        real(dp), dimension(3) :: cell_size
        real(dp) :: conductivity_left, conductivity_right
        real(dp) :: temperature_gradient_left, temperature_gradient_right
        real(dp) :: left_ratio, right_ratio, source_value, radius

        if (.not. ieee_is_finite(time_step) .or. time_step <= 0.0_dp) then
            error stop 'Heat-transfer solver: time step must be finite and positive'
        end if

        call this%update_thermal_conductivity()

        dimensions = this%domain%get_domain_dimensions()
        power = this%geometry_power()
        cell_loop = this%domain%get_local_inner_cells_bounds()
        cell_size = this%mesh%mesh_ptr%get_cell_edges_length()
        this%energy_source%s_ptr%cells = 0.0_dp

        associate( &
            source => this%energy_source%s_ptr, &
            temperature => this%temperature%s_ptr, &
            conductivity => this%thermal_conductivity%s_ptr, &
            markers => this%boundary%bc_ptr%bc_markers, &
            mesh => this%mesh%mesh_ptr)

            !$omp parallel do collapse(3) schedule(static) default(shared) &
            !$omp& private(i,j,k,dim,left_i,left_j,left_k,right_i,right_j, &
            !$omp& right_k,left_boundary,right_boundary,conductivity_left, &
            !$omp& conductivity_right,temperature_gradient_left, &
            !$omp& temperature_gradient_right,left_ratio,right_ratio, &
            !$omp& source_value,radius)
            do k = cell_loop(3,1), cell_loop(3,2)
                do j = cell_loop(2,1), cell_loop(2,2)
                    do i = cell_loop(1,1), cell_loop(1,2)
                        if (markers(i,j,k) /= 0) cycle

                        source_value = 0.0_dp
                        radius = mesh%mesh(1,i,j,k)

                        do dim = 1, dimensions
                            left_i = i - I_m(dim,1)
                            left_j = j - I_m(dim,2)
                            left_k = k - I_m(dim,3)
                            right_i = i + I_m(dim,1)
                            right_j = j + I_m(dim,2)
                            right_k = k + I_m(dim,3)

                            left_boundary = markers(left_i,left_j,left_k)
                            right_boundary = markers(right_i,right_j,right_k)

                            conductivity_left = this%face_conductivity( &
                                conductivity%cells(i,j,k), &
                                conductivity%cells(left_i,left_j,left_k), &
                                left_boundary)
                            conductivity_right = this%face_conductivity( &
                                conductivity%cells(i,j,k), &
                                conductivity%cells(right_i,right_j,right_k), &
                                right_boundary)

                            temperature_gradient_left = &
                                (temperature%cells(i,j,k) - &
                                 temperature%cells(left_i,left_j,left_k)) / &
                                cell_size(dim)
                            temperature_gradient_right = &
                                (temperature%cells(right_i,right_j,right_k) - &
                                 temperature%cells(i,j,k)) / cell_size(dim)

                            left_ratio = this%radial_face_ratio( &
                                power, dim, radius, cell_size(1), -1)
                            right_ratio = this%radial_face_ratio( &
                                power, dim, radius, cell_size(1), 1)

                            source_value = source_value + &
                                (right_ratio*conductivity_right* &
                                 temperature_gradient_right - &
                                 left_ratio*conductivity_left* &
                                 temperature_gradient_left) / cell_size(dim)
                        end do

                        source%cells(i,j,k) = source_value
                    end do
                end do
            end do
            !$omp end parallel do
        end associate

        call this%mpi_support%exchange_conservative_scalar_field( &
            this%energy_source%s_ptr)
    end subroutine solve_heat_transfer


    subroutine update_thermal_conductivity(this)
        class(heat_transfer_solver), intent(inout) :: this

        integer :: species_number, i, j, k, specie
        integer, dimension(3,2) :: cell_loop
        real(dp), allocatable :: mass_fractions(:)
        real(dp) :: conductivity_value

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
            conductivity => this%thermal_conductivity%s_ptr, &
            mixture_molar_mass => this%mixture_molar_mass%s_ptr, &
            mass_fraction => this%mass_fraction%v_ptr, &
            markers => this%boundary%bc_ptr%bc_markers)

            !$omp parallel default(shared) private(i,j,k,specie,mass_fractions, &
            !$omp& conductivity_value)
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
                        conductivity_value = &
                            this%thermophysics%thermo_ptr% &
                            mixture_thermal_conductivity( &
                                temperature%cells(i,j,k), mass_fractions, &
                                mixture_molar_mass%cells(i,j,k))
                        if (.not. ieee_is_finite(conductivity_value) .or. &
                            conductivity_value <= 0.0_dp) then
                            error stop 'Heat-transfer solver: invalid conductivity'
                        end if
                        conductivity%cells(i,j,k) = conductivity_value
                    end do
                end do
            end do
            !$omp end do
            deallocate(mass_fractions)
            !$omp end parallel
        end associate

        call this%mpi_support%exchange_conservative_scalar_field( &
            this%thermal_conductivity%s_ptr)
        call this%fill_conductivity_boundary_ghosts()
    end subroutine update_thermal_conductivity


    subroutine fill_conductivity_boundary_ghosts(this)
        class(heat_transfer_solver), intent(inout) :: this

        integer :: dimensions, i, j, k, dim, side, sign
        integer :: ghost_i, ghost_j, ghost_k
        integer, dimension(3,2) :: cell_loop

        dimensions = this%domain%get_domain_dimensions()
        cell_loop = this%domain%get_local_inner_cells_bounds()

        associate( &
            conductivity => this%thermal_conductivity%s_ptr, &
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
                                    conductivity%cells(ghost_i,ghost_j,ghost_k) = &
                                        conductivity%cells(i,j,k)
                                end if
                            end do
                        end do
                    end do
                end do
            end do
        end associate
    end subroutine fill_conductivity_boundary_ghosts


    real(dp) function face_conductivity(this, cell_value, neighbour_value, &
        boundary_number) result(value)
        class(heat_transfer_solver), intent(in) :: this
        real(dp), intent(in) :: cell_value, neighbour_value
        integer, intent(in) :: boundary_number

        character(len=20) :: boundary_name
        real(dp) :: conductivity_ratio

        if (boundary_number == 0) then
            if (cell_value <= 0.0_dp .or. neighbour_value <= 0.0_dp) then
                error stop 'Heat-transfer solver: non-positive face conductivity'
            end if
            value = 2.0_dp*cell_value*neighbour_value / &
                (cell_value + neighbour_value)
            return
        end if

        boundary_name = this%boundary%bc_ptr% &
            boundary_types(boundary_number)%get_type_name()
        select case(trim(boundary_name))
        case('wall')
            if (.not. this%boundary%bc_ptr% &
                boundary_types(boundary_number)%is_conductive()) then
                value = 0.0_dp
                return
            end if
            conductivity_ratio = this%boundary%bc_ptr% &
                boundary_types(boundary_number)%get_wall_conductivity_ratio()
            if (.not. ieee_is_finite(conductivity_ratio) .or. &
                conductivity_ratio < 0.0_dp) then
                error stop 'Heat-transfer solver: invalid wall conductivity ratio'
            end if
            ! Preserve the submitted boundary model: the ghost-side material
            ! has conductivity ratio*kappa_gas and the face value is averaged.
            value = 0.5_dp*(1.0_dp + conductivity_ratio)*cell_value
        case default
            ! Inlet, outlet, symmetry, and unrecognized physical boundaries are
            ! treated as zero normal conductive-flux boundaries.
            value = 0.0_dp
        end select
    end function face_conductivity


    integer function geometry_power(this) result(power)
        class(heat_transfer_solver), intent(in) :: this

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
            error stop 'Heat-transfer solver: non-positive radial cell metric'
        end if

        face_radius = radius + 0.5_dp*real(side,dp)*radial_spacing
        if (face_radius < 0.0_dp) then
            error stop 'Heat-transfer solver: negative radial face coordinate'
        end if
        ratio = face_radius**power/radius**power
    end function radial_face_ratio

end module fourier_heat_transfer_solver_class
