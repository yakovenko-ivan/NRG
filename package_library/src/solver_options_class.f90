module solver_options_class

    use kind_parameters, only: dp
    use global_data, only: solver_data_file_name

    implicit none

    private
    public :: particles_phase
    public :: solver_options, solver_options_c

    integer, parameter :: option_name_length = 24

    !> Configuration of one dispersed material phase.
    !>
    !> The model and initialization fields are deliberately part of the
    !> problem-interface data.  Selecting the Eulerian or Lagrangian backend
    !> therefore never requires editing or recompiling a gas-dynamic solver.
    type :: particles_phase
        character(len=option_name_length) :: model = 'continuous'
        character(len=option_name_length) :: initial_condition = 'field'
        character(len=option_name_length) :: boundary_interaction = 'reflect'

        real(dp) :: diameter = 0.0_dp
        real(dp) :: material_heat_capacity = 0.0_dp
        real(dp) :: material_density = 0.0_dp
        real(dp) :: material_latent_heat = 0.0_dp
        real(dp) :: material_boiling_temperature = 0.0_dp
        character(len=20) :: material = ''
        logical :: evaporating = .false.
        logical :: heating = .false.
        logical :: inertial = .true.

        ! Initial state shared by both backends.
        real(dp) :: initial_temperature = 300.0_dp
        real(dp), dimension(3) :: initial_velocity = 0.0_dp
        real(dp) :: initial_mass_loading = -1.0_dp

        ! Lagrangian parcel placement.  Coordinates and region bounds are
        ! fractions of the physical domain when coordinates_are_relative is true.
        logical :: coordinates_are_relative = .true.
        real(dp), dimension(3) :: initial_position = 0.5_dp
        real(dp), dimension(3) :: initial_region_min = 0.0_dp
        real(dp), dimension(3) :: initial_region_max = 1.0_dp
        integer, dimension(3) :: initial_lattice = 1

        ! Eulerian inflow state.  A negative loading requests zero-gradient
        ! extrapolation instead of prescribed injection.
        real(dp) :: inlet_mass_loading = -1.0_dp
        real(dp) :: inlet_temperature = 300.0_dp
        real(dp), dimension(3) :: inlet_velocity = 0.0_dp
        real(dp) :: injection_start_time = 0.0_dp
        real(dp) :: injection_ramp_time = 0.0_dp

        ! Lagrangian diagnostics.  A non-positive interval disables snapshots.
        real(dp) :: output_interval = 0.0_dp
    end type particles_phase

    type :: solver_options
        private
        character(len=20) :: solver_name = ''
        logical :: hydrodynamics_flag = .false.
        logical :: heat_transfer_flag = .false.
        logical :: molecular_diffusion_flag = .false.
        logical :: soret_diffusion_flag = .false.
        logical :: viscosity_flag = .false.
        logical :: thermal_radiation_flag = .false.
        logical :: chemical_reaction_flag = .false.
        logical :: CFL_flag = .false.
        real(dp) :: CFL_coefficient = 0.0_dp
        real(dp) :: initial_time_step = 0.0_dp
        real(dp), dimension(3) :: grav_acc = 0.0_dp
        integer :: additional_particles_phases = 0
        integer :: particles_phase_counter = 0
        type(particles_phase), dimension(:), allocatable :: particles
    contains
        procedure, private :: read_properties
        procedure, private :: write_properties
        procedure, private :: set_properties
        procedure, private :: validate_phase

        procedure :: create_additional_phase
        procedure :: get_solver_name
        procedure :: get_CFL_condition_coefficient
        procedure :: get_initial_time_step
        procedure :: get_CFL_condition_flag
        procedure :: get_hydrodynamics_flag
        procedure :: get_heat_transfer_flag
        procedure :: get_thermal_radiation_flag
        procedure :: get_molecular_diffusion_flag
        procedure :: get_soret_diffusion_flag
        procedure :: get_viscosity_flag
        procedure :: get_chemical_reaction_flag
        procedure :: get_grav_acc
        procedure :: get_additional_particles_phases_number
        procedure :: get_particles_params
        procedure :: write_log
    end type solver_options

    interface solver_options_c
        module procedure constructor
        module procedure constructor_file
    end interface solver_options_c

contains

    type(solver_options) function constructor( &
        solver_name, hydrodynamics_flag, heat_transfer_flag, &
        molecular_diffusion_flag, viscosity_flag, thermal_radiation_flag, &
        chemical_reaction_flag, grav_acc, additional_particles_phases, &
        CFL_flag, CFL_coefficient, initial_time_step, soret_diffusion_flag)

        character(len=*), intent(in) :: solver_name
        logical, intent(in) :: hydrodynamics_flag
        logical, intent(in) :: heat_transfer_flag
        logical, intent(in) :: molecular_diffusion_flag
        logical, intent(in) :: viscosity_flag
        logical, intent(in), optional :: thermal_radiation_flag
        logical, intent(in) :: chemical_reaction_flag
        real(dp), dimension(3), intent(in) :: grav_acc
        integer, intent(in), optional :: additional_particles_phases
        logical, intent(in) :: CFL_flag
        real(dp), intent(in) :: CFL_coefficient
        real(dp), intent(in) :: initial_time_step
        logical, intent(in), optional :: soret_diffusion_flag

        integer :: number_of_phases
        integer :: io_unit
        logical :: radiation_enabled, soret_enabled

        radiation_enabled = .false.
        if (present(thermal_radiation_flag)) radiation_enabled = thermal_radiation_flag

        soret_enabled = .false.
        if (present(soret_diffusion_flag)) soret_enabled = soret_diffusion_flag

        number_of_phases = 0
        if (present(additional_particles_phases)) then
            number_of_phases = additional_particles_phases
        end if

        call constructor%set_properties( &
            solver_name, hydrodynamics_flag, heat_transfer_flag, &
            molecular_diffusion_flag, soret_enabled, viscosity_flag, &
            radiation_enabled, chemical_reaction_flag, grav_acc, &
            number_of_phases, CFL_flag, CFL_coefficient, initial_time_step)

        open(newunit=io_unit, file=solver_data_file_name, status='replace', &
            form='formatted', delim='quote')
        call constructor%write_properties(io_unit)
        close(io_unit)
    end function constructor


    type(solver_options) function constructor_file()
        integer :: io_unit, phase

        character(len=option_name_length) :: particles_model
        character(len=option_name_length) :: particles_initial_condition
        character(len=option_name_length) :: particles_boundary_interaction
        character(len=20) :: particles_material
        real(dp) :: particles_diameter
        real(dp) :: particles_material_heat_capacity
        real(dp) :: particles_material_density
        real(dp) :: particles_material_latent_heat
        real(dp) :: particles_material_boiling_temperature
        logical :: particles_evaporating, particles_heating, particles_inertial
        real(dp) :: particles_initial_temperature
        real(dp), dimension(3) :: particles_initial_velocity
        real(dp) :: particles_initial_mass_loading
        logical :: particles_coordinates_are_relative
        real(dp), dimension(3) :: particles_initial_position
        real(dp), dimension(3) :: particles_initial_region_min
        real(dp), dimension(3) :: particles_initial_region_max
        integer, dimension(3) :: particles_initial_lattice
        real(dp) :: particles_inlet_mass_loading
        real(dp) :: particles_inlet_temperature
        real(dp), dimension(3) :: particles_inlet_velocity
        real(dp) :: particles_injection_start_time
        real(dp) :: particles_injection_ramp_time
        real(dp) :: particles_output_interval

        namelist /dispersed_phase/ particles_model, &
            particles_initial_condition, particles_boundary_interaction, &
            particles_diameter, particles_material_heat_capacity, &
            particles_material_density, particles_material_latent_heat, &
            particles_material_boiling_temperature, particles_material, &
            particles_evaporating, particles_heating, particles_inertial, &
            particles_initial_temperature, particles_initial_velocity, &
            particles_initial_mass_loading, particles_coordinates_are_relative, &
            particles_initial_position, particles_initial_region_min, &
            particles_initial_region_max, particles_initial_lattice, &
            particles_inlet_mass_loading, particles_inlet_temperature, &
            particles_inlet_velocity, particles_injection_start_time, &
            particles_injection_ramp_time, particles_output_interval

        open(newunit=io_unit, file=solver_data_file_name, status='old', &
            form='formatted')
        call constructor_file%read_properties(io_unit)

        do phase = 1, constructor_file%additional_particles_phases
            ! Defaults make old solver-data files backward compatible.
            particles_model = 'continuous'
            particles_initial_condition = 'field'
            particles_boundary_interaction = 'reflect'
            particles_diameter = 0.0_dp
            particles_material_heat_capacity = 0.0_dp
            particles_material_density = 0.0_dp
            particles_material_latent_heat = 0.0_dp
            particles_material_boiling_temperature = 0.0_dp
            particles_material = ''
            particles_evaporating = .false.
            particles_heating = .false.
            particles_inertial = .true.
            particles_initial_temperature = 300.0_dp
            particles_initial_velocity = 0.0_dp
            particles_initial_mass_loading = -1.0_dp
            particles_coordinates_are_relative = .true.
            particles_initial_position = 0.5_dp
            particles_initial_region_min = 0.0_dp
            particles_initial_region_max = 1.0_dp
            particles_initial_lattice = 1
            particles_inlet_mass_loading = -1.0_dp
            particles_inlet_temperature = 300.0_dp
            particles_inlet_velocity = 0.0_dp
            particles_injection_start_time = 0.0_dp
            particles_injection_ramp_time = 0.0_dp
            particles_output_interval = 0.0_dp

            read(unit=io_unit, nml=dispersed_phase)

            associate(p => constructor_file%particles(phase))
                p%model = particles_model
                p%initial_condition = particles_initial_condition
                p%boundary_interaction = particles_boundary_interaction
                p%diameter = particles_diameter
                p%material = particles_material
                p%material_heat_capacity = particles_material_heat_capacity
                p%material_density = particles_material_density
                p%material_latent_heat = particles_material_latent_heat
                p%material_boiling_temperature = &
                    particles_material_boiling_temperature
                p%evaporating = particles_evaporating
                p%heating = particles_heating
                p%inertial = particles_inertial
                p%initial_temperature = particles_initial_temperature
                p%initial_velocity = particles_initial_velocity
                p%initial_mass_loading = particles_initial_mass_loading
                p%coordinates_are_relative = particles_coordinates_are_relative
                p%initial_position = particles_initial_position
                p%initial_region_min = particles_initial_region_min
                p%initial_region_max = particles_initial_region_max
                p%initial_lattice = particles_initial_lattice
                p%inlet_mass_loading = particles_inlet_mass_loading
                p%inlet_temperature = particles_inlet_temperature
                p%inlet_velocity = particles_inlet_velocity
                p%injection_start_time = particles_injection_start_time
                p%injection_ramp_time = particles_injection_ramp_time
                p%output_interval = particles_output_interval
            end associate

            call constructor_file%validate_phase(phase)
            if (phase > 1) then
                if (is_lagrangian_model(constructor_file%particles(phase)%model) .neqv. &
                    is_lagrangian_model(constructor_file%particles(1)%model)) then
                    error stop 'solver_options: all dispersed phases must use one model'
                end if
            end if
        end do

        close(io_unit)
    end function constructor_file


    subroutine write_properties(this, solver_data_unit)
        class(solver_options), intent(in) :: this
        integer, intent(in) :: solver_data_unit

        character(len=20) :: solver_name
        logical :: hydrodynamics_flag, heat_transfer_flag
        logical :: molecular_diffusion_flag, soret_diffusion_flag
        logical :: viscosity_flag
        logical :: thermal_radiation_flag, chemical_reaction_flag, CFL_flag
        real(dp), dimension(3) :: grav_acc
        integer :: additional_particles_phases
        real(dp) :: CFL_coefficient, initial_time_step

        namelist /solver_properties/ solver_name, hydrodynamics_flag, &
            heat_transfer_flag, molecular_diffusion_flag, soret_diffusion_flag, &
            viscosity_flag, thermal_radiation_flag, chemical_reaction_flag, grav_acc, &
            additional_particles_phases, CFL_flag, CFL_coefficient, &
            initial_time_step

        solver_name = this%solver_name
        hydrodynamics_flag = this%hydrodynamics_flag
        heat_transfer_flag = this%heat_transfer_flag
        molecular_diffusion_flag = this%molecular_diffusion_flag
        soret_diffusion_flag = this%soret_diffusion_flag
        viscosity_flag = this%viscosity_flag
        thermal_radiation_flag = this%thermal_radiation_flag
        chemical_reaction_flag = this%chemical_reaction_flag
        grav_acc = this%grav_acc
        additional_particles_phases = this%additional_particles_phases
        CFL_flag = this%CFL_flag
        CFL_coefficient = this%CFL_coefficient
        initial_time_step = this%initial_time_step

        write(unit=solver_data_unit, nml=solver_properties)
    end subroutine write_properties


    subroutine read_properties(this, solver_data_unit)
        class(solver_options), intent(inout) :: this
        integer, intent(in) :: solver_data_unit

        character(len=20) :: solver_name
        logical :: hydrodynamics_flag, heat_transfer_flag
        logical :: molecular_diffusion_flag, soret_diffusion_flag
        logical :: viscosity_flag
        logical :: thermal_radiation_flag, chemical_reaction_flag, CFL_flag
        real(dp), dimension(3) :: grav_acc
        integer :: additional_particles_phases
        real(dp) :: CFL_coefficient, initial_time_step

        namelist /solver_properties/ solver_name, hydrodynamics_flag, &
            heat_transfer_flag, molecular_diffusion_flag, soret_diffusion_flag, &
            viscosity_flag, thermal_radiation_flag, chemical_reaction_flag, grav_acc, &
            additional_particles_phases, CFL_flag, CFL_coefficient, &
            initial_time_step

        ! Preserve compatibility with solver-data files written before the
        ! Soret option was introduced.
        soret_diffusion_flag = .false.

        read(unit=solver_data_unit, nml=solver_properties)
        call this%set_properties( &
            solver_name, hydrodynamics_flag, heat_transfer_flag, &
            molecular_diffusion_flag, soret_diffusion_flag, viscosity_flag, &
            thermal_radiation_flag, &
            chemical_reaction_flag, grav_acc, additional_particles_phases, &
            CFL_flag, CFL_coefficient, initial_time_step)
    end subroutine read_properties


    subroutine set_properties(this, solver_name, hydrodynamics_flag, &
        heat_transfer_flag, molecular_diffusion_flag, soret_diffusion_flag, &
        viscosity_flag, &
        thermal_radiation_flag, chemical_reaction_flag, grav_acc, &
        additional_particles_phases, CFL_flag, CFL_coefficient, &
        initial_time_step)

        class(solver_options), intent(inout) :: this
        character(len=*), intent(in) :: solver_name
        logical, intent(in) :: hydrodynamics_flag, heat_transfer_flag
        logical, intent(in) :: molecular_diffusion_flag, soret_diffusion_flag
        logical, intent(in) :: viscosity_flag
        logical, intent(in) :: thermal_radiation_flag, chemical_reaction_flag
        real(dp), dimension(3), intent(in) :: grav_acc
        integer, intent(in) :: additional_particles_phases
        logical, intent(in) :: CFL_flag
        real(dp), intent(in) :: CFL_coefficient, initial_time_step

        if (additional_particles_phases < 0) then
            error stop 'solver_options: negative number of dispersed phases'
        end if
        if (CFL_flag .and. CFL_coefficient <= 0.0_dp) then
            error stop 'solver_options: CFL coefficient must be positive'
        end if
        if (initial_time_step <= 0.0_dp) then
            error stop 'solver_options: initial time step must be positive'
        end if
        if (soret_diffusion_flag .and. .not. molecular_diffusion_flag) then
            error stop 'solver_options: Soret diffusion requires molecular diffusion'
        end if

        this%solver_name = solver_name
        this%hydrodynamics_flag = hydrodynamics_flag
        this%heat_transfer_flag = heat_transfer_flag
        this%molecular_diffusion_flag = molecular_diffusion_flag
        this%soret_diffusion_flag = soret_diffusion_flag
        this%viscosity_flag = viscosity_flag
        this%thermal_radiation_flag = thermal_radiation_flag
        this%chemical_reaction_flag = chemical_reaction_flag
        this%grav_acc = grav_acc
        this%additional_particles_phases = additional_particles_phases
        this%CFL_flag = CFL_flag
        this%CFL_coefficient = CFL_coefficient
        this%initial_time_step = initial_time_step
        this%particles_phase_counter = 0

        if (allocated(this%particles)) deallocate(this%particles)
        allocate(this%particles(additional_particles_phases))
    end subroutine set_properties


    subroutine create_additional_phase(this, particles_parameters)
        class(solver_options), intent(inout) :: this
        type(particles_phase), intent(in) :: particles_parameters

        integer :: io_unit, phase
        character(len=option_name_length) :: particles_model
        character(len=option_name_length) :: particles_initial_condition
        character(len=option_name_length) :: particles_boundary_interaction
        character(len=20) :: particles_material
        real(dp) :: particles_diameter
        real(dp) :: particles_material_heat_capacity
        real(dp) :: particles_material_density
        real(dp) :: particles_material_latent_heat
        real(dp) :: particles_material_boiling_temperature
        logical :: particles_evaporating, particles_heating, particles_inertial
        real(dp) :: particles_initial_temperature
        real(dp), dimension(3) :: particles_initial_velocity
        real(dp) :: particles_initial_mass_loading
        logical :: particles_coordinates_are_relative
        real(dp), dimension(3) :: particles_initial_position
        real(dp), dimension(3) :: particles_initial_region_min
        real(dp), dimension(3) :: particles_initial_region_max
        integer, dimension(3) :: particles_initial_lattice
        real(dp) :: particles_inlet_mass_loading
        real(dp) :: particles_inlet_temperature
        real(dp), dimension(3) :: particles_inlet_velocity
        real(dp) :: particles_injection_start_time
        real(dp) :: particles_injection_ramp_time
        real(dp) :: particles_output_interval

        namelist /dispersed_phase/ particles_model, &
            particles_initial_condition, particles_boundary_interaction, &
            particles_diameter, particles_material_heat_capacity, &
            particles_material_density, particles_material_latent_heat, &
            particles_material_boiling_temperature, particles_material, &
            particles_evaporating, particles_heating, particles_inertial, &
            particles_initial_temperature, particles_initial_velocity, &
            particles_initial_mass_loading, particles_coordinates_are_relative, &
            particles_initial_position, particles_initial_region_min, &
            particles_initial_region_max, particles_initial_lattice, &
            particles_inlet_mass_loading, particles_inlet_temperature, &
            particles_inlet_velocity, particles_injection_start_time, &
            particles_injection_ramp_time, particles_output_interval

        phase = this%particles_phase_counter + 1
        if (phase > this%additional_particles_phases) then
            error stop 'solver_options: too many dispersed phases were supplied'
        end if

        this%particles(phase) = particles_parameters
        this%particles_phase_counter = phase
        call this%validate_phase(phase)
        if (phase > 1) then
            if (is_lagrangian_model(this%particles(phase)%model) .neqv. &
                is_lagrangian_model(this%particles(1)%model)) then
                error stop 'solver_options: all dispersed phases must use one model'
            end if
        end if

        associate(p => this%particles(phase))
            particles_model = p%model
            particles_initial_condition = p%initial_condition
            particles_boundary_interaction = p%boundary_interaction
            particles_diameter = p%diameter
            particles_material = p%material
            particles_material_heat_capacity = p%material_heat_capacity
            particles_material_density = p%material_density
            particles_material_latent_heat = p%material_latent_heat
            particles_material_boiling_temperature = &
                p%material_boiling_temperature
            particles_evaporating = p%evaporating
            particles_heating = p%heating
            particles_inertial = p%inertial
            particles_initial_temperature = p%initial_temperature
            particles_initial_velocity = p%initial_velocity
            particles_initial_mass_loading = p%initial_mass_loading
            particles_coordinates_are_relative = p%coordinates_are_relative
            particles_initial_position = p%initial_position
            particles_initial_region_min = p%initial_region_min
            particles_initial_region_max = p%initial_region_max
            particles_initial_lattice = p%initial_lattice
            particles_inlet_mass_loading = p%inlet_mass_loading
            particles_inlet_temperature = p%inlet_temperature
            particles_inlet_velocity = p%inlet_velocity
            particles_injection_start_time = p%injection_start_time
            particles_injection_ramp_time = p%injection_ramp_time
            particles_output_interval = p%output_interval
        end associate

        open(newunit=io_unit, file=solver_data_file_name, status='old', &
            form='formatted', position='append')
        write(unit=io_unit, nml=dispersed_phase)
        close(io_unit)
    end subroutine create_additional_phase


    subroutine validate_phase(this, phase)
        class(solver_options), intent(in) :: this
        integer, intent(in) :: phase
        character(len=option_name_length) :: model, initialization, interaction

        if (phase < 1 .or. phase > size(this%particles)) then
            error stop 'solver_options: dispersed phase index is out of range'
        end if

        associate(p => this%particles(phase))
            model = lowercase(trim(p%model))
            initialization = lowercase(trim(p%initial_condition))
            interaction = lowercase(trim(p%boundary_interaction))

            select case (model)
            case ('continuous', 'eulerian', 'lagrangian', 'parcel')
            case default
                error stop 'solver_options: unknown dispersed-phase model'
            end select

            select case (initialization)
            case ('field', 'uniform', 'none', 'single', 'lattice')
            case default
                error stop 'solver_options: unknown dispersed-phase initialization'
            end select

            if ((model == 'continuous' .or. model == 'eulerian') .and. &
                initialization /= 'field' .and. initialization /= 'uniform' .and. &
                initialization /= 'none') then
                error stop 'solver_options: continuous phase requires field, uniform, or none'
            end if
            if ((model == 'lagrangian' .or. model == 'parcel') .and. &
                initialization /= 'single' .and. initialization /= 'lattice' .and. &
                initialization /= 'none') then
                error stop 'solver_options: Lagrangian phase requires single, lattice, or none'
            end if

            select case (interaction)
            case ('reflect', 'stick', 'escape')
            case default
                error stop 'solver_options: unknown particle boundary interaction'
            end select

            if (len_trim(p%material) == 0) then
                error stop 'solver_options: particle material must be specified'
            end if
            if (p%diameter <= 0.0_dp) then
                error stop 'solver_options: particle diameter must be positive'
            end if
            if (p%material_density <= 0.0_dp) then
                error stop 'solver_options: particle material density must be positive'
            end if
            if (p%heating .and. p%material_heat_capacity <= 0.0_dp) then
                error stop 'solver_options: particle heat capacity must be positive'
            end if
            if (p%evaporating .and. p%material_latent_heat <= 0.0_dp) then
                error stop 'solver_options: particle latent heat must be positive'
            end if
            if (p%evaporating .and. p%material_boiling_temperature <= 0.0_dp) then
                error stop 'solver_options: particle boiling temperature must be positive'
            end if
            if (any(p%initial_lattice < 1)) then
                error stop 'solver_options: initial particle lattice must be positive'
            end if
            if (any(p%initial_region_max <= p%initial_region_min)) then
                error stop 'solver_options: invalid initial particle region'
            end if
            if (initialization == 'uniform' .and. &
                p%initial_mass_loading < 0.0_dp) then
                error stop 'solver_options: uniform phase loading cannot be negative'
            end if
            if (p%coordinates_are_relative .and. &
                (any(p%initial_position < 0.0_dp) .or. &
                any(p%initial_position > 1.0_dp))) then
                error stop 'solver_options: relative particle position must lie in [0,1]'
            end if
            if (p%initial_temperature <= 0.0_dp .or. &
                p%inlet_temperature <= 0.0_dp) then
                error stop 'solver_options: particle temperature must be positive'
            end if
            if (p%injection_ramp_time < 0.0_dp) then
                error stop 'solver_options: injection ramp time cannot be negative'
            end if
        end associate
    end subroutine validate_phase


    subroutine write_log(this, log_unit)
        class(solver_options), intent(in) :: this
        integer, intent(in) :: log_unit
        integer :: phase

        write(log_unit, '(A)') repeat('*', 84)
        write(log_unit, '(A,A)') ' Solver name                 : ', trim(this%solver_name)
        write(log_unit, '(A,L1)') ' Hydrodynamics               : ', this%hydrodynamics_flag
        write(log_unit, '(A,L1)') ' Heat transfer               : ', this%heat_transfer_flag
        write(log_unit, '(A,L1)') ' Molecular diffusion         : ', &
            this%molecular_diffusion_flag
        write(log_unit, '(A,L1)') ' Soret diffusion             : ', &
            this%soret_diffusion_flag
        write(log_unit, '(A,L1)') ' Viscosity                   : ', this%viscosity_flag
        write(log_unit, '(A,L1)') ' Thermal radiation           : ', &
            this%thermal_radiation_flag
        write(log_unit, '(A,L1)') ' Chemical reaction           : ', &
            this%chemical_reaction_flag
        write(log_unit, '(A,3ES14.6)') ' Gravitational acceleration  : ', this%grav_acc
        write(log_unit, '(A,I0)') ' Dispersed phases            : ', &
            this%additional_particles_phases
        write(log_unit, '(A,L1)') ' CFL condition               : ', this%CFL_flag
        write(log_unit, '(A,ES14.6)') ' CFL coefficient             : ', &
            this%CFL_coefficient
        write(log_unit, '(A,ES14.6)') ' Initial time step           : ', &
            this%initial_time_step

        do phase = 1, this%additional_particles_phases
            write(log_unit, '(A,I0,A,A,A,A)') ' Phase ', phase, ': model=', &
                trim(this%particles(phase)%model), ', initialization=', &
                trim(this%particles(phase)%initial_condition)
        end do
        write(log_unit, '(A)') repeat('*', 84)
    end subroutine write_log


    pure function get_solver_name(this) result(value)
        class(solver_options), intent(in) :: this
        character(len=20) :: value
        value = this%solver_name
    end function get_solver_name

    pure function get_hydrodynamics_flag(this) result(value)
        class(solver_options), intent(in) :: this
        logical :: value
        value = this%hydrodynamics_flag
    end function get_hydrodynamics_flag

    pure function get_heat_transfer_flag(this) result(value)
        class(solver_options), intent(in) :: this
        logical :: value
        value = this%heat_transfer_flag
    end function get_heat_transfer_flag

    pure function get_molecular_diffusion_flag(this) result(value)
        class(solver_options), intent(in) :: this
        logical :: value
        value = this%molecular_diffusion_flag
    end function get_molecular_diffusion_flag

    pure function get_soret_diffusion_flag(this) result(value)
        class(solver_options), intent(in) :: this
        logical :: value
        value = this%soret_diffusion_flag
    end function get_soret_diffusion_flag

    pure function get_viscosity_flag(this) result(value)
        class(solver_options), intent(in) :: this
        logical :: value
        value = this%viscosity_flag
    end function get_viscosity_flag

    pure function get_thermal_radiation_flag(this) result(value)
        class(solver_options), intent(in) :: this
        logical :: value
        value = this%thermal_radiation_flag
    end function get_thermal_radiation_flag

    pure function get_chemical_reaction_flag(this) result(value)
        class(solver_options), intent(in) :: this
        logical :: value
        value = this%chemical_reaction_flag
    end function get_chemical_reaction_flag

    pure function get_additional_particles_phases_number(this) result(value)
        class(solver_options), intent(in) :: this
        integer :: value
        value = this%additional_particles_phases
    end function get_additional_particles_phases_number

    pure function get_CFL_condition_flag(this) result(value)
        class(solver_options), intent(in) :: this
        logical :: value
        value = this%CFL_flag
    end function get_CFL_condition_flag

    pure function get_CFL_condition_coefficient(this) result(value)
        class(solver_options), intent(in) :: this
        real(dp) :: value
        value = this%CFL_coefficient
    end function get_CFL_condition_coefficient

    pure function get_initial_time_step(this) result(value)
        class(solver_options), intent(in) :: this
        real(dp) :: value
        value = this%initial_time_step
    end function get_initial_time_step

    function get_particles_params(this, phase_number) result(value)
        class(solver_options), intent(in) :: this
        integer, intent(in) :: phase_number
        type(particles_phase) :: value

        if (phase_number < 1 .or. phase_number > size(this%particles)) then
            error stop 'solver_options: dispersed phase index is out of range'
        end if
        value = this%particles(phase_number)
    end function get_particles_params

    pure function get_grav_acc(this) result(value)
        class(solver_options), intent(in) :: this
        real(dp), dimension(3) :: value
        value = this%grav_acc
    end function get_grav_acc


    pure function is_lagrangian_model(model) result(value)
        character(len=*), intent(in) :: model
        logical :: value
        character(len=len(model)) :: normalized

        normalized = lowercase(trim(model))
        value = normalized == 'lagrangian' .or. normalized == 'parcel'
    end function is_lagrangian_model


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

end module solver_options_class
