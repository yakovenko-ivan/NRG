module thermal_radiation_solver_class

    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite

    use kind_parameters, only: dp
    use global_data, only: sigma_SB
    use field_pointers
    use boundary_conditions_class
    use data_manager_class
    use computational_domain_class
    use thermophysical_properties_class
    use chemical_properties_class
    use mpi_communications_class

    implicit none

    private
    public :: thermal_radiation_solver, thermal_radiation_solver_c

    real(dp), parameter :: reference_pressure = 101325.0_dp
    real(dp), parameter :: default_background_temperature = 300.0_dp
    integer, parameter :: absorption_polynomial_degree = 5

    type :: thermal_radiation_solver
        private

        type(field_scalar_cons_pointer) :: energy_source
        type(field_scalar_cons_pointer) :: temperature
        type(field_scalar_cons_pointer) :: pressure
        type(field_scalar_cons_pointer) :: mixture_molar_mass
        type(field_vector_cons_pointer) :: mass_fraction

        type(field_scalar_cons), pointer :: energy_source_store => null()

        type(computational_domain) :: domain
        type(mpi_communications) :: mpi_support
        type(boundary_conditions_pointer) :: boundary
        type(thermophysical_properties_pointer) :: thermophysics
        type(chemical_properties_pointer) :: chemistry

        real(dp), allocatable :: absorption_coefficients(:,:)
        real(dp) :: background_temperature = default_background_temperature
    contains
        procedure :: solve_radiation
        procedure :: set_background_temperature
        procedure :: get_background_temperature
        procedure, private :: initialize_absorption_model
        procedure, private :: absorption_coefficient
    end type thermal_radiation_solver

    interface thermal_radiation_solver_c
        module procedure constructor
    end interface thermal_radiation_solver_c

contains

    type(thermal_radiation_solver) function constructor(manager, background_temperature)
        type(data_manager), intent(inout) :: manager
        real(dp), intent(in), optional :: background_temperature

        type(field_scalar_cons_pointer) :: scalar_pointer
        type(field_vector_cons_pointer) :: vector_pointer
        type(field_tensor_cons_pointer) :: tensor_pointer

        constructor%domain = manager%domain
        constructor%mpi_support = manager%mpi_communications
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
            scalar_pointer, vector_pointer, tensor_pointer, 'mixture_molar_mass')
        constructor%mixture_molar_mass%s_ptr => scalar_pointer%s_ptr
        call manager%get_cons_field_pointer_by_name( &
            scalar_pointer, vector_pointer, tensor_pointer, 'specie_mass_fraction')
        constructor%mass_fraction%v_ptr => vector_pointer%v_ptr

        allocate(constructor%energy_source_store)
        call manager%create_scalar_field( &
            constructor%energy_source_store, &
            'energy_production_radiation', 'E_f_prod_radiation')
        constructor%energy_source%s_ptr => constructor%energy_source_store
        constructor%energy_source%s_ptr%cells = 0.0_dp

        call constructor%initialize_absorption_model()
        if (present(background_temperature)) then
            call constructor%set_background_temperature(background_temperature)
        end if
    end function constructor


    subroutine initialize_absorption_model(this)
        class(thermal_radiation_solver), intent(inout) :: this

        integer :: water_index, species_number

        species_number = this%chemistry%chem_ptr%species_number
        allocate(this%absorption_coefficients( &
            species_number, 0:absorption_polynomial_degree))
        this%absorption_coefficients = 0.0_dp

        water_index = this%chemistry%chem_ptr%get_chemical_specie_index('H2O')
        if (water_index <= 0) return

        ! Mean H2O absorption-coefficient fit retained from the submitted model.
        ! The polynomial variable is 1000/T and pressure is supplied in atm.
        this%absorption_coefficients(water_index,0) = -0.23093_dp
        this%absorption_coefficients(water_index,1) =  1.12390_dp
        this%absorption_coefficients(water_index,2) =  9.41530_dp
        this%absorption_coefficients(water_index,3) = -2.99880_dp
        this%absorption_coefficients(water_index,4) =  0.51382_dp
        this%absorption_coefficients(water_index,5) = -1.86840e-05_dp
    end subroutine initialize_absorption_model


    subroutine set_background_temperature(this, value)
        class(thermal_radiation_solver), intent(inout) :: this
        real(dp), intent(in) :: value

        if (.not. ieee_is_finite(value) .or. value <= 0.0_dp) then
            error stop 'Radiation solver: background temperature must be finite and positive'
        end if
        this%background_temperature = value
    end subroutine set_background_temperature


    pure real(dp) function get_background_temperature(this) result(value)
        class(thermal_radiation_solver), intent(in) :: this

        value = this%background_temperature
    end function get_background_temperature


    subroutine solve_radiation(this, time_step)
        class(thermal_radiation_solver), intent(inout) :: this
        real(dp), intent(in) :: time_step

        integer :: i, j, k, specie, species_number
        integer, dimension(3,2) :: cell_loop
        real(dp) :: temperature_cell, pressure_atm, mixture_molar_mass_cell
        real(dp) :: mole_fraction, planck_mean_absorption

        if (.not. ieee_is_finite(time_step) .or. time_step <= 0.0_dp) then
            error stop 'Radiation solver: time step must be finite and positive'
        end if

        species_number = this%chemistry%chem_ptr%species_number
        cell_loop = this%domain%get_local_inner_cells_bounds()
        this%energy_source%s_ptr%cells = 0.0_dp

        associate( &
            source => this%energy_source%s_ptr, &
            temperature => this%temperature%s_ptr, &
            pressure => this%pressure%s_ptr, &
            mixture_molar_mass => this%mixture_molar_mass%s_ptr, &
            mass_fraction => this%mass_fraction%v_ptr, &
            markers => this%boundary%bc_ptr%bc_markers, &
            molar_masses => this%thermophysics%thermo_ptr%molar_masses)

            !$omp parallel do collapse(3) schedule(static) default(shared) &
            !$omp& private(i,j,k,specie,temperature_cell,pressure_atm, &
            !$omp& mixture_molar_mass_cell,mole_fraction,planck_mean_absorption)
            do k = cell_loop(3,1), cell_loop(3,2)
                do j = cell_loop(2,1), cell_loop(2,2)
                    do i = cell_loop(1,1), cell_loop(1,2)
                        if (markers(i,j,k) /= 0) cycle

                        temperature_cell = temperature%cells(i,j,k)
                        mixture_molar_mass_cell = mixture_molar_mass%cells(i,j,k)
                        pressure_atm = pressure%cells(i,j,k)/reference_pressure

                        if (.not. ieee_is_finite(temperature_cell) .or. &
                            temperature_cell <= 0.0_dp) then
                            error stop 'Radiation solver: non-positive cell temperature'
                        end if
                        if (.not. ieee_is_finite(pressure_atm) .or. pressure_atm < 0.0_dp) then
                            error stop 'Radiation solver: invalid cell pressure'
                        end if
                        if (.not. ieee_is_finite(mixture_molar_mass_cell) .or. &
                            mixture_molar_mass_cell <= 0.0_dp) then
                            error stop 'Radiation solver: invalid mixture molar mass'
                        end if

                        planck_mean_absorption = 0.0_dp
                        do specie = 1, species_number
                            if (molar_masses(specie) <= 0.0_dp) cycle
                            if (mass_fraction%pr(specie)%cells(i,j,k) < -1.0e-12_dp) then
                                error stop 'Radiation solver: negative species mass fraction'
                            end if
                            mole_fraction = max( &
                                mass_fraction%pr(specie)%cells(i,j,k), 0.0_dp) * &
                                mixture_molar_mass_cell/molar_masses(specie)
                            planck_mean_absorption = planck_mean_absorption + &
                                mole_fraction*pressure_atm * &
                                this%absorption_coefficient(specie, temperature_cell)
                        end do

                        source%cells(i,j,k) = -4.0_dp*sigma_SB * &
                            (temperature_cell**4 - this%background_temperature**4) * &
                            planck_mean_absorption
                    end do
                end do
            end do
            !$omp end parallel do
        end associate

        call this%mpi_support%exchange_conservative_scalar_field( &
            this%energy_source%s_ptr)
    end subroutine solve_radiation


    pure real(dp) function absorption_coefficient(this, specie, temperature) result(value)
        class(thermal_radiation_solver), intent(in) :: this
        integer, intent(in) :: specie
        real(dp), intent(in) :: temperature

        integer :: coefficient
        real(dp) :: inverse_temperature

        inverse_temperature = 1000.0_dp/temperature
        value = this%absorption_coefficients(specie, absorption_polynomial_degree)
        do coefficient = absorption_polynomial_degree - 1, 0, -1
            value = value*inverse_temperature + &
                this%absorption_coefficients(specie, coefficient)
        end do

        ! A negative absorption coefficient is non-physical and can only arise
        ! from extrapolating the empirical polynomial beyond its useful range.
        value = max(value, 0.0_dp)
    end function absorption_coefficient

end module thermal_radiation_solver_class
