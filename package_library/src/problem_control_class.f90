module problem_control_class

    use kind_parameters, only: dp
    use global_data, only: problem_controls_data_file_name

    implicit none

    private
    public :: flame_stabilization_control
    public :: flame_stabilization_control_c
    public :: energy_ignition_control
    public :: energy_ignition_control_c
    public :: turbulence_forcing_control
    public :: turbulence_forcing_control_c
    public :: problem_controls
    public :: problem_controls_c

    integer, parameter :: mode_length = 40
    integer, parameter :: scientific_title_length = 240
    integer, parameter :: case_setup_length = 80

    real(dp), parameter :: default_time_delay = 1.0e-03_dp
    real(dp), parameter :: default_time_track = 1.0e-04_dp
    real(dp), parameter :: default_time_control = 5.0e-04_dp
    real(dp), parameter :: default_response_settle_time = 5.0e-04_dp
    real(dp), parameter :: default_measurement_max_duration = 1.0_dp
    real(dp), parameter :: default_measurement_min_displacement_cells = 2.5e-01_dp

    type :: flame_stabilization_control
        private
        character(len=mode_length) :: mode = 'none'
        character(len=scientific_title_length) :: scientific_result_title = &
            'NRG laminar flame result'
        character(len=case_setup_length) :: case_setup = 'unspecified'
        real(dp) :: time_delay = default_time_delay
        real(dp) :: time_track = default_time_track
        real(dp) :: time_control = default_time_control
        real(dp) :: response_settle_time = default_response_settle_time
        real(dp) :: measurement_max_duration = default_measurement_max_duration
        real(dp) :: measurement_min_displacement_cells = &
            default_measurement_min_displacement_cells
    contains
        procedure, private :: set_properties => set_flame_stabilization_properties
        procedure :: is_enabled => flame_stabilization_is_enabled
        procedure :: is_anchor => flame_stabilization_is_anchor
        procedure :: is_laminar_burning_velocity => &
            flame_stabilization_is_laminar_burning_velocity
        procedure :: get_mode => get_flame_stabilization_mode
        procedure :: get_scientific_result_title => &
            get_flame_stabilization_scientific_result_title
        procedure :: get_case_setup => get_flame_stabilization_case_setup
        procedure :: get_time_delay => get_flame_stabilization_time_delay
        procedure :: get_time_track => get_flame_stabilization_time_track
        procedure :: get_time_control => get_flame_stabilization_time_control
        procedure :: get_response_settle_time => &
            get_flame_stabilization_response_settle_time
        procedure :: get_measurement_max_duration => &
            get_flame_stabilization_measurement_max_duration
        procedure :: get_measurement_min_displacement_cells => &
            get_flame_stabilization_measurement_min_displacement_cells
        procedure :: validate_solver_compatibility => &
            validate_flame_stabilization_solver_compatibility
        procedure :: validate_problem_compatibility => &
            validate_flame_stabilization_problem_compatibility
    end type flame_stabilization_control

    type :: energy_ignition_control
        private
        character(len=mode_length) :: mode = 'none'
        real(dp) :: start_time = 0.0_dp
        real(dp) :: duration = 0.0_dp
        real(dp) :: total_energy = 0.0_dp
        real(dp), dimension(3) :: center = 0.0_dp
        real(dp) :: radius = 0.0_dp
    contains
        procedure, private :: set_properties => set_energy_ignition_properties
        procedure :: is_enabled => energy_ignition_is_enabled
        procedure :: get_mode => get_energy_ignition_mode
        procedure :: get_start_time => get_energy_ignition_start_time
        procedure :: get_duration => get_energy_ignition_duration
        procedure :: get_total_energy => get_energy_ignition_total_energy
        procedure :: get_center => get_energy_ignition_center
        procedure :: get_radius => get_energy_ignition_radius
        procedure :: validate_solver_compatibility => &
            validate_energy_ignition_solver_compatibility
    end type energy_ignition_control

    type :: turbulence_forcing_control
        private
        character(len=mode_length) :: mode = 'none'
        real(dp) :: start_time = 0.0_dp
        real(dp) :: duration = 0.0_dp
        real(dp) :: ramp_time = 0.0_dp
        real(dp) :: acceleration_amplitude = 0.0_dp
        real(dp) :: integral_scale = 0.0_dp
        real(dp) :: update_interval = 0.0_dp
        real(dp), dimension(3) :: center = 0.0_dp
        real(dp) :: radius = 0.0_dp
        integer :: random_seed = 1
    contains
        procedure, private :: set_properties => set_turbulence_forcing_properties
        procedure :: is_enabled => turbulence_forcing_is_enabled
        procedure :: get_mode => get_turbulence_forcing_mode
        procedure :: get_start_time => get_turbulence_forcing_start_time
        procedure :: get_duration => get_turbulence_forcing_duration
        procedure :: get_ramp_time => get_turbulence_forcing_ramp_time
        procedure :: get_acceleration_amplitude => &
            get_turbulence_forcing_acceleration_amplitude
        procedure :: get_integral_scale => get_turbulence_forcing_integral_scale
        procedure :: get_update_interval => get_turbulence_forcing_update_interval
        procedure :: get_center => get_turbulence_forcing_center
        procedure :: get_radius => get_turbulence_forcing_radius
        procedure :: get_random_seed => get_turbulence_forcing_random_seed
        procedure :: validate_solver_compatibility => &
            validate_turbulence_forcing_solver_compatibility
        procedure :: validate_problem_compatibility => &
            validate_turbulence_forcing_problem_compatibility
    end type turbulence_forcing_control

    type :: problem_controls
        private
        type(flame_stabilization_control) :: flame_stabilization
        type(energy_ignition_control) :: energy_ignition
        type(turbulence_forcing_control) :: turbulence_forcing
    contains
        procedure, private :: read_properties => read_problem_controls
        procedure, private :: write_properties => write_problem_controls
        procedure :: get_flame_stabilization
        procedure :: get_energy_ignition
        procedure :: get_turbulence_forcing
        procedure :: validate_solver_compatibility => &
            validate_problem_controls_solver_compatibility
        procedure :: write_log => write_problem_controls_log
    end type problem_controls

    interface flame_stabilization_control_c
        module procedure flame_stabilization_constructor
    end interface flame_stabilization_control_c

    interface energy_ignition_control_c
        module procedure energy_ignition_constructor
    end interface energy_ignition_control_c

    interface turbulence_forcing_control_c
        module procedure turbulence_forcing_constructor
    end interface turbulence_forcing_control_c

    interface problem_controls_c
        module procedure problem_controls_constructor
        module procedure problem_controls_constructor_file
    end interface problem_controls_c

contains

    type(flame_stabilization_control) function flame_stabilization_constructor( &
            mode, time_delay, time_track, time_control, response_settle_time, &
            measurement_max_duration, measurement_min_displacement_cells, &
            scientific_result_title, case_setup)

        character(len=*), intent(in) :: mode
        character(len=*), intent(in), optional :: scientific_result_title
        character(len=*), intent(in), optional :: case_setup
        real(dp), intent(in), optional :: time_delay
        real(dp), intent(in), optional :: time_track
        real(dp), intent(in), optional :: time_control
        real(dp), intent(in), optional :: response_settle_time
        real(dp), intent(in), optional :: measurement_max_duration
        real(dp), intent(in), optional :: measurement_min_displacement_cells

        real(dp) :: time_delay_set, time_track_set, time_control_set
        real(dp) :: response_settle_time_set, measurement_max_duration_set
        real(dp) :: measurement_min_displacement_cells_set
        character(len=scientific_title_length) :: scientific_result_title_set
        character(len=case_setup_length) :: case_setup_set

        time_delay_set = default_time_delay
        time_track_set = default_time_track
        time_control_set = default_time_control
        response_settle_time_set = default_response_settle_time
        measurement_max_duration_set = default_measurement_max_duration
        measurement_min_displacement_cells_set = &
            default_measurement_min_displacement_cells
        scientific_result_title_set = 'NRG laminar flame result'
        case_setup_set = 'unspecified'

        if (present(scientific_result_title)) &
            scientific_result_title_set = trim(scientific_result_title)
        if (present(case_setup)) case_setup_set = trim(case_setup)
        if (present(time_delay)) time_delay_set = time_delay
        if (present(time_track)) time_track_set = time_track
        if (present(time_control)) time_control_set = time_control
        if (present(response_settle_time)) &
            response_settle_time_set = response_settle_time
        if (present(measurement_max_duration)) &
            measurement_max_duration_set = measurement_max_duration
        if (present(measurement_min_displacement_cells)) &
            measurement_min_displacement_cells_set = &
                measurement_min_displacement_cells

        call flame_stabilization_constructor%set_properties( &
            mode, time_delay_set, time_track_set, time_control_set, &
            response_settle_time_set, measurement_max_duration_set, &
            measurement_min_displacement_cells_set, scientific_result_title_set, &
            case_setup_set)
    end function flame_stabilization_constructor


    type(energy_ignition_control) function energy_ignition_constructor( &
            mode, start_time, duration, total_energy, center, radius)
        character(len=*), intent(in) :: mode
        real(dp), intent(in), optional :: start_time, duration, total_energy, radius
        real(dp), dimension(3), intent(in), optional :: center

        real(dp) :: start_time_set, duration_set, total_energy_set, radius_set
        real(dp), dimension(3) :: center_set

        start_time_set = 0.0_dp
        duration_set = 0.0_dp
        total_energy_set = 0.0_dp
        center_set = 0.0_dp
        radius_set = 0.0_dp

        if (present(start_time)) start_time_set = start_time
        if (present(duration)) duration_set = duration
        if (present(total_energy)) total_energy_set = total_energy
        if (present(center)) center_set = center
        if (present(radius)) radius_set = radius

        call energy_ignition_constructor%set_properties( &
            mode, start_time_set, duration_set, total_energy_set, center_set, radius_set)
    end function energy_ignition_constructor


    type(turbulence_forcing_control) function turbulence_forcing_constructor( &
            mode, start_time, duration, ramp_time, acceleration_amplitude, &
            integral_scale, update_interval, center, radius, random_seed)
        character(len=*), intent(in) :: mode
        real(dp), intent(in), optional :: start_time, duration, ramp_time
        real(dp), intent(in), optional :: acceleration_amplitude
        real(dp), intent(in), optional :: integral_scale, update_interval, radius
        real(dp), dimension(3), intent(in), optional :: center
        integer, intent(in), optional :: random_seed

        real(dp) :: start_time_set, duration_set, ramp_time_set
        real(dp) :: acceleration_amplitude_set, integral_scale_set
        real(dp) :: update_interval_set, radius_set
        real(dp), dimension(3) :: center_set
        integer :: random_seed_set

        start_time_set = 0.0_dp
        duration_set = 0.0_dp
        ramp_time_set = 0.0_dp
        acceleration_amplitude_set = 0.0_dp
        integral_scale_set = 0.0_dp
        update_interval_set = 0.0_dp
        center_set = 0.0_dp
        radius_set = 0.0_dp
        random_seed_set = 1

        if (present(start_time)) start_time_set = start_time
        if (present(duration)) duration_set = duration
        if (present(ramp_time)) ramp_time_set = ramp_time
        if (present(acceleration_amplitude)) &
            acceleration_amplitude_set = acceleration_amplitude
        if (present(integral_scale)) integral_scale_set = integral_scale
        if (present(update_interval)) update_interval_set = update_interval
        if (present(center)) center_set = center
        if (present(radius)) radius_set = radius
        if (present(random_seed)) random_seed_set = random_seed

        call turbulence_forcing_constructor%set_properties( &
            mode, start_time_set, duration_set, ramp_time_set, &
            acceleration_amplitude_set, integral_scale_set, update_interval_set, &
            center_set, radius_set, random_seed_set)
    end function turbulence_forcing_constructor


    type(problem_controls) function problem_controls_constructor( &
            flame_stabilization, energy_ignition, turbulence_forcing)
        type(flame_stabilization_control), intent(in) :: flame_stabilization
        type(energy_ignition_control), intent(in), optional :: energy_ignition
        type(turbulence_forcing_control), intent(in), optional :: turbulence_forcing

        integer :: io_unit

        problem_controls_constructor%flame_stabilization = flame_stabilization
        problem_controls_constructor%energy_ignition = &
            energy_ignition_control_c('none')
        problem_controls_constructor%turbulence_forcing = &
            turbulence_forcing_control_c('none')

        if (present(energy_ignition)) &
            problem_controls_constructor%energy_ignition = energy_ignition
        if (present(turbulence_forcing)) &
            problem_controls_constructor%turbulence_forcing = turbulence_forcing

        open(newunit=io_unit, file=problem_controls_data_file_name, &
            status='replace', form='formatted', delim='quote')
        call problem_controls_constructor%write_properties(io_unit)
        close(io_unit)
    end function problem_controls_constructor


    type(problem_controls) function problem_controls_constructor_file()
        integer :: io_unit, io_status
        logical :: file_exists

        problem_controls_constructor_file%flame_stabilization = &
            flame_stabilization_control_c('none')
        problem_controls_constructor_file%energy_ignition = &
            energy_ignition_control_c('none')
        problem_controls_constructor_file%turbulence_forcing = &
            turbulence_forcing_control_c('none')

        inquire(file=problem_controls_data_file_name, exist=file_exists)
        if (.not. file_exists) return

        open(newunit=io_unit, file=problem_controls_data_file_name, &
            status='old', form='formatted', action='read', iostat=io_status)
        if (io_status /= 0) then
            error stop 'problem_controls: unable to open problem_controls.inf'
        end if

        call problem_controls_constructor_file%read_properties(io_unit)
        close(io_unit)
    end function problem_controls_constructor_file


    subroutine set_flame_stabilization_properties(this, mode, time_delay, &
            time_track, time_control, response_settle_time, &
            measurement_max_duration, measurement_min_displacement_cells, &
            scientific_result_title, case_setup)
        class(flame_stabilization_control), intent(inout) :: this
        character(len=*), intent(in) :: mode
        character(len=*), intent(in) :: scientific_result_title
        character(len=*), intent(in) :: case_setup
        real(dp), intent(in) :: time_delay, time_track, time_control
        real(dp), intent(in) :: response_settle_time, measurement_max_duration
        real(dp), intent(in) :: measurement_min_displacement_cells
        character(len=mode_length) :: normalized_mode

        normalized_mode = normalize_mode(mode)
        select case (trim(normalized_mode))
        case ('none', 'anchor', 'laminar_burning_velocity')
        case default
            error stop 'problem_controls: unknown flame stabilization mode'
        end select

        if (time_delay < 0.0_dp) &
            error stop 'problem_controls: flame time_delay cannot be negative'
        if (time_track <= 0.0_dp) &
            error stop 'problem_controls: flame time_track must be positive'
        if (time_control <= 0.0_dp) &
            error stop 'problem_controls: flame time_control must be positive'
        if (response_settle_time < 0.0_dp) &
            error stop 'problem_controls: response_settle_time cannot be negative'
        if (measurement_max_duration <= 0.0_dp) &
            error stop 'problem_controls: measurement_max_duration must be positive'
        if (measurement_min_displacement_cells <= 0.0_dp) &
            error stop 'problem_controls: measurement displacement must be positive'

        if (len_trim(scientific_result_title) == 0) &
            error stop 'problem_controls: scientific result title cannot be empty'
        if (len_trim(case_setup) == 0) &
            error stop 'problem_controls: case setup cannot be empty'

        this%mode = normalized_mode
        this%scientific_result_title = trim(scientific_result_title)
        this%case_setup = trim(case_setup)
        this%time_delay = time_delay
        this%time_track = time_track
        this%time_control = time_control
        this%response_settle_time = response_settle_time
        this%measurement_max_duration = measurement_max_duration
        this%measurement_min_displacement_cells = &
            measurement_min_displacement_cells
    end subroutine set_flame_stabilization_properties


    subroutine set_energy_ignition_properties( &
            this, mode, start_time, duration, total_energy, center, radius)
        class(energy_ignition_control), intent(inout) :: this
        character(len=*), intent(in) :: mode
        real(dp), intent(in) :: start_time, duration, total_energy, radius
        real(dp), dimension(3), intent(in) :: center
        character(len=mode_length) :: normalized_mode

        normalized_mode = normalize_mode(mode)
        select case (trim(normalized_mode))
        case ('none')
        case ('uniform_sphere')
            if (start_time < 0.0_dp) &
                error stop 'problem_controls: ignition start_time cannot be negative'
            if (duration <= 0.0_dp) &
                error stop 'problem_controls: ignition duration must be positive'
            if (total_energy <= 0.0_dp) &
                error stop 'problem_controls: ignition total_energy must be positive'
            if (radius <= 0.0_dp) &
                error stop 'problem_controls: ignition radius must be positive'
        case default
            error stop 'problem_controls: unknown energy ignition mode'
        end select

        this%mode = normalized_mode
        this%start_time = start_time
        this%duration = duration
        this%total_energy = total_energy
        this%center = center
        this%radius = radius
    end subroutine set_energy_ignition_properties


    subroutine set_turbulence_forcing_properties( &
            this, mode, start_time, duration, ramp_time, acceleration_amplitude, &
            integral_scale, update_interval, center, radius, random_seed)
        class(turbulence_forcing_control), intent(inout) :: this
        character(len=*), intent(in) :: mode
        real(dp), intent(in) :: start_time, duration, ramp_time
        real(dp), intent(in) :: acceleration_amplitude
        real(dp), intent(in) :: integral_scale, update_interval, radius
        real(dp), dimension(3), intent(in) :: center
        integer, intent(in) :: random_seed
        character(len=mode_length) :: normalized_mode

        normalized_mode = normalize_mode(mode)
        select case (trim(normalized_mode))
        case ('none')
        case ('planar_fourier')
            if (start_time < 0.0_dp) &
                error stop 'problem_controls: turbulence start_time cannot be negative'
            if (duration <= 0.0_dp) &
                error stop 'problem_controls: turbulence duration must be positive'
            if (ramp_time < 0.0_dp) &
                error stop 'problem_controls: turbulence ramp_time cannot be negative'
            if (acceleration_amplitude < 0.0_dp) &
                error stop 'problem_controls: turbulence acceleration must be non-negative'
            if (integral_scale <= 0.0_dp) &
                error stop 'problem_controls: turbulence integral_scale must be positive'
            if (update_interval <= 0.0_dp) &
                error stop 'problem_controls: turbulence update_interval must be positive'
            if (radius < 0.0_dp) &
                error stop 'problem_controls: turbulence radius cannot be negative'
            if (random_seed <= 0) &
                error stop 'problem_controls: turbulence random_seed must be positive'
        case default
            error stop 'problem_controls: unknown turbulence forcing mode'
        end select

        this%mode = normalized_mode
        this%start_time = start_time
        this%duration = duration
        this%ramp_time = ramp_time
        this%acceleration_amplitude = acceleration_amplitude
        this%integral_scale = integral_scale
        this%update_interval = update_interval
        this%center = center
        this%radius = radius
        this%random_seed = random_seed
    end subroutine set_turbulence_forcing_properties


    subroutine write_problem_controls(this, io_unit)
        class(problem_controls), intent(in) :: this
        integer, intent(in) :: io_unit

        character(len=mode_length) :: flame_stabilization_mode
        character(len=scientific_title_length) :: scientific_result_title
        character(len=case_setup_length) :: case_setup
        real(dp) :: time_delay, time_track, time_control, response_settle_time
        real(dp) :: measurement_max_duration, measurement_min_displacement_cells
        character(len=mode_length) :: energy_ignition_mode
        real(dp) :: ignition_start_time, ignition_duration, ignition_total_energy
        real(dp), dimension(3) :: ignition_center
        real(dp) :: ignition_radius
        character(len=mode_length) :: turbulence_forcing_mode
        real(dp) :: turbulence_start_time, turbulence_duration, turbulence_ramp_time
        real(dp) :: turbulence_acceleration_amplitude, turbulence_integral_scale
        real(dp) :: turbulence_update_interval, turbulence_radius
        real(dp), dimension(3) :: turbulence_center
        integer :: turbulence_random_seed

        namelist /problem_controls_parameters/ flame_stabilization_mode, &
            scientific_result_title, case_setup, &
            time_delay, time_track, time_control, response_settle_time, &
            measurement_max_duration, measurement_min_displacement_cells
        namelist /energy_ignition_parameters/ energy_ignition_mode, &
            ignition_start_time, ignition_duration, ignition_total_energy, &
            ignition_center, ignition_radius
        namelist /turbulence_forcing_parameters/ turbulence_forcing_mode, &
            turbulence_start_time, turbulence_duration, turbulence_ramp_time, &
            turbulence_acceleration_amplitude, turbulence_integral_scale, &
            turbulence_update_interval, turbulence_center, turbulence_radius, &
            turbulence_random_seed

        flame_stabilization_mode = this%flame_stabilization%get_mode()
        scientific_result_title = &
            this%flame_stabilization%get_scientific_result_title()
        case_setup = this%flame_stabilization%get_case_setup()
        time_delay = this%flame_stabilization%get_time_delay()
        time_track = this%flame_stabilization%get_time_track()
        time_control = this%flame_stabilization%get_time_control()
        response_settle_time = this%flame_stabilization%get_response_settle_time()
        measurement_max_duration = &
            this%flame_stabilization%get_measurement_max_duration()
        measurement_min_displacement_cells = &
            this%flame_stabilization%get_measurement_min_displacement_cells()

        energy_ignition_mode = this%energy_ignition%get_mode()
        ignition_start_time = this%energy_ignition%get_start_time()
        ignition_duration = this%energy_ignition%get_duration()
        ignition_total_energy = this%energy_ignition%get_total_energy()
        ignition_center = this%energy_ignition%get_center()
        ignition_radius = this%energy_ignition%get_radius()

        turbulence_forcing_mode = this%turbulence_forcing%get_mode()
        turbulence_start_time = this%turbulence_forcing%get_start_time()
        turbulence_duration = this%turbulence_forcing%get_duration()
        turbulence_ramp_time = this%turbulence_forcing%get_ramp_time()
        turbulence_acceleration_amplitude = &
            this%turbulence_forcing%get_acceleration_amplitude()
        turbulence_integral_scale = this%turbulence_forcing%get_integral_scale()
        turbulence_update_interval = this%turbulence_forcing%get_update_interval()
        turbulence_center = this%turbulence_forcing%get_center()
        turbulence_radius = this%turbulence_forcing%get_radius()
        turbulence_random_seed = this%turbulence_forcing%get_random_seed()

        write(unit=io_unit, nml=problem_controls_parameters)
        write(unit=io_unit, nml=energy_ignition_parameters)
        write(unit=io_unit, nml=turbulence_forcing_parameters)
    end subroutine write_problem_controls


    subroutine read_problem_controls(this, io_unit)
        class(problem_controls), intent(inout) :: this
        integer, intent(in) :: io_unit

        character(len=mode_length) :: flame_stabilization_mode
        character(len=scientific_title_length) :: scientific_result_title
        character(len=case_setup_length) :: case_setup
        real(dp) :: time_delay, time_track, time_control, response_settle_time
        real(dp) :: measurement_max_duration, measurement_min_displacement_cells
        character(len=mode_length) :: energy_ignition_mode
        real(dp) :: ignition_start_time, ignition_duration, ignition_total_energy
        real(dp), dimension(3) :: ignition_center
        real(dp) :: ignition_radius
        character(len=mode_length) :: turbulence_forcing_mode
        real(dp) :: turbulence_start_time, turbulence_duration, turbulence_ramp_time
        real(dp) :: turbulence_acceleration_amplitude, turbulence_integral_scale
        real(dp) :: turbulence_update_interval, turbulence_radius
        real(dp), dimension(3) :: turbulence_center
        integer :: turbulence_random_seed
        integer :: io_status

        namelist /problem_controls_parameters/ flame_stabilization_mode, &
            scientific_result_title, case_setup, &
            time_delay, time_track, time_control, response_settle_time, &
            measurement_max_duration, measurement_min_displacement_cells
        namelist /energy_ignition_parameters/ energy_ignition_mode, &
            ignition_start_time, ignition_duration, ignition_total_energy, &
            ignition_center, ignition_radius
        namelist /turbulence_forcing_parameters/ turbulence_forcing_mode, &
            turbulence_start_time, turbulence_duration, turbulence_ramp_time, &
            turbulence_acceleration_amplitude, turbulence_integral_scale, &
            turbulence_update_interval, turbulence_center, turbulence_radius, &
            turbulence_random_seed

        flame_stabilization_mode = 'none'
        scientific_result_title = 'NRG laminar flame result'
        case_setup = 'unspecified'
        time_delay = default_time_delay
        time_track = default_time_track
        time_control = default_time_control
        response_settle_time = default_response_settle_time
        measurement_max_duration = default_measurement_max_duration
        measurement_min_displacement_cells = &
            default_measurement_min_displacement_cells

        read(unit=io_unit, nml=problem_controls_parameters, iostat=io_status)
        if (io_status /= 0) then
            error stop 'problem_controls: invalid problem_controls_parameters namelist'
        end if
        this%flame_stabilization = flame_stabilization_control_c( &
            flame_stabilization_mode, time_delay=time_delay, &
            time_track=time_track, time_control=time_control, &
            response_settle_time=response_settle_time, &
            measurement_max_duration=measurement_max_duration, &
            measurement_min_displacement_cells= &
                measurement_min_displacement_cells, &
            scientific_result_title=trim(scientific_result_title), &
            case_setup=trim(case_setup))

        energy_ignition_mode = 'none'
        ignition_start_time = 0.0_dp
        ignition_duration = 0.0_dp
        ignition_total_energy = 0.0_dp
        ignition_center = 0.0_dp
        ignition_radius = 0.0_dp
        read(unit=io_unit, nml=energy_ignition_parameters, iostat=io_status)
        if (io_status == 0) then
            this%energy_ignition = energy_ignition_control_c( &
                energy_ignition_mode, start_time=ignition_start_time, &
                duration=ignition_duration, total_energy=ignition_total_energy, &
                center=ignition_center, radius=ignition_radius)
        else if (io_status > 0) then
            error stop 'problem_controls: invalid energy_ignition_parameters namelist'
        end if

        turbulence_forcing_mode = 'none'
        turbulence_start_time = 0.0_dp
        turbulence_duration = 0.0_dp
        turbulence_ramp_time = 0.0_dp
        turbulence_acceleration_amplitude = 0.0_dp
        turbulence_integral_scale = 0.0_dp
        turbulence_update_interval = 0.0_dp
        turbulence_center = 0.0_dp
        turbulence_radius = 0.0_dp
        turbulence_random_seed = 1
        read(unit=io_unit, nml=turbulence_forcing_parameters, iostat=io_status)
        if (io_status == 0) then
            this%turbulence_forcing = turbulence_forcing_control_c( &
                turbulence_forcing_mode, start_time=turbulence_start_time, &
                duration=turbulence_duration, ramp_time=turbulence_ramp_time, &
                acceleration_amplitude=turbulence_acceleration_amplitude, &
                integral_scale=turbulence_integral_scale, &
                update_interval=turbulence_update_interval, &
                center=turbulence_center, radius=turbulence_radius, &
                random_seed=turbulence_random_seed)
        else if (io_status > 0) then
            error stop 'problem_controls: invalid turbulence_forcing_parameters namelist'
        end if
    end subroutine read_problem_controls


    function get_flame_stabilization(this) result(control)
        class(problem_controls), intent(in) :: this
        type(flame_stabilization_control) :: control
        control = this%flame_stabilization
    end function get_flame_stabilization

    function get_energy_ignition(this) result(control)
        class(problem_controls), intent(in) :: this
        type(energy_ignition_control) :: control
        control = this%energy_ignition
    end function get_energy_ignition

    function get_turbulence_forcing(this) result(control)
        class(problem_controls), intent(in) :: this
        type(turbulence_forcing_control) :: control
        control = this%turbulence_forcing
    end function get_turbulence_forcing


    subroutine validate_problem_controls_solver_compatibility(this, solver_name)
        class(problem_controls), intent(in) :: this
        character(len=*), intent(in) :: solver_name

        call this%flame_stabilization%validate_solver_compatibility(solver_name)
        call this%energy_ignition%validate_solver_compatibility(solver_name)
        call this%turbulence_forcing%validate_solver_compatibility(solver_name)
    end subroutine validate_problem_controls_solver_compatibility


    subroutine validate_flame_stabilization_solver_compatibility(this, solver_name)
        class(flame_stabilization_control), intent(in) :: this
        character(len=*), intent(in) :: solver_name

        if (.not. this%is_enabled()) return
        if (trim(solver_name) /= 'fds_low_mach') then
            error stop 'problem_controls: flame stabilization is currently supported only by fds_low_mach'
        end if
    end subroutine validate_flame_stabilization_solver_compatibility


    subroutine validate_energy_ignition_solver_compatibility(this, solver_name)
        class(energy_ignition_control), intent(in) :: this
        character(len=*), intent(in) :: solver_name

        if (.not. this%is_enabled()) return
        select case (trim(solver_name))
        case ('fds_low_mach', 'CABARET', 'CABARET_low_mach')
        case default
            error stop 'problem_controls: energy ignition is unsupported by the selected solver'
        end select
    end subroutine validate_energy_ignition_solver_compatibility


    subroutine validate_turbulence_forcing_solver_compatibility(this, solver_name)
        class(turbulence_forcing_control), intent(in) :: this
        character(len=*), intent(in) :: solver_name

        if (.not. this%is_enabled()) return
        select case (trim(solver_name))
        case ('fds_low_mach', 'CABARET', 'CABARET_low_mach')
        case default
            error stop 'problem_controls: turbulence forcing is unsupported by the selected solver'
        end select
    end subroutine validate_turbulence_forcing_solver_compatibility


    subroutine validate_flame_stabilization_problem_compatibility(this, &
            reactive, dimensions, inlet_count, outlet_count)
        class(flame_stabilization_control), intent(in) :: this
        logical, intent(in) :: reactive
        integer, intent(in) :: dimensions, inlet_count, outlet_count

        if (.not. this%is_enabled()) return
        if (.not. reactive) &
            error stop 'problem_controls: flame stabilization requires reactive flow'
        if (dimensions /= 1) &
            error stop 'problem_controls: flame stabilization currently requires a 1D problem'
        if (inlet_count /= 1) &
            error stop 'problem_controls: flame stabilization requires exactly one inlet'
        if (outlet_count < 1) &
            error stop 'problem_controls: flame stabilization requires an outlet'
    end subroutine validate_flame_stabilization_problem_compatibility


    subroutine validate_turbulence_forcing_problem_compatibility(this, dimensions)
        class(turbulence_forcing_control), intent(in) :: this
        integer, intent(in) :: dimensions

        if (.not. this%is_enabled()) return
        if (trim(this%mode) == 'planar_fourier' .and. dimensions < 2) then
            error stop 'problem_controls: planar_fourier turbulence requires at least 2D'
        end if
    end subroutine validate_turbulence_forcing_problem_compatibility


    logical function flame_stabilization_is_enabled(this)
        class(flame_stabilization_control), intent(in) :: this
        flame_stabilization_is_enabled = trim(this%mode) /= 'none'
    end function flame_stabilization_is_enabled

    logical function flame_stabilization_is_anchor(this)
        class(flame_stabilization_control), intent(in) :: this
        flame_stabilization_is_anchor = trim(this%mode) == 'anchor'
    end function flame_stabilization_is_anchor

    logical function flame_stabilization_is_laminar_burning_velocity(this)
        class(flame_stabilization_control), intent(in) :: this
        flame_stabilization_is_laminar_burning_velocity = &
            trim(this%mode) == 'laminar_burning_velocity'
    end function flame_stabilization_is_laminar_burning_velocity

    logical function energy_ignition_is_enabled(this)
        class(energy_ignition_control), intent(in) :: this
        energy_ignition_is_enabled = trim(this%mode) /= 'none'
    end function energy_ignition_is_enabled

    logical function turbulence_forcing_is_enabled(this)
        class(turbulence_forcing_control), intent(in) :: this
        turbulence_forcing_is_enabled = trim(this%mode) /= 'none'
    end function turbulence_forcing_is_enabled


    function get_flame_stabilization_mode(this) result(value)
        class(flame_stabilization_control), intent(in) :: this
        character(len=mode_length) :: value
        value = this%mode
    end function get_flame_stabilization_mode
    function get_flame_stabilization_scientific_result_title(this) result(value)
        class(flame_stabilization_control), intent(in) :: this
        character(len=scientific_title_length) :: value
        value = this%scientific_result_title
    end function get_flame_stabilization_scientific_result_title
    function get_flame_stabilization_case_setup(this) result(value)
        class(flame_stabilization_control), intent(in) :: this
        character(len=case_setup_length) :: value
        value = this%case_setup
    end function get_flame_stabilization_case_setup
    real(dp) function get_flame_stabilization_time_delay(this)
        class(flame_stabilization_control), intent(in) :: this
        get_flame_stabilization_time_delay = this%time_delay
    end function get_flame_stabilization_time_delay
    real(dp) function get_flame_stabilization_time_track(this)
        class(flame_stabilization_control), intent(in) :: this
        get_flame_stabilization_time_track = this%time_track
    end function get_flame_stabilization_time_track
    real(dp) function get_flame_stabilization_time_control(this)
        class(flame_stabilization_control), intent(in) :: this
        get_flame_stabilization_time_control = this%time_control
    end function get_flame_stabilization_time_control
    real(dp) function get_flame_stabilization_response_settle_time(this)
        class(flame_stabilization_control), intent(in) :: this
        get_flame_stabilization_response_settle_time = this%response_settle_time
    end function get_flame_stabilization_response_settle_time
    real(dp) function get_flame_stabilization_measurement_max_duration(this)
        class(flame_stabilization_control), intent(in) :: this
        get_flame_stabilization_measurement_max_duration = this%measurement_max_duration
    end function get_flame_stabilization_measurement_max_duration
    real(dp) function get_flame_stabilization_measurement_min_displacement_cells(this)
        class(flame_stabilization_control), intent(in) :: this
        get_flame_stabilization_measurement_min_displacement_cells = &
            this%measurement_min_displacement_cells
    end function get_flame_stabilization_measurement_min_displacement_cells

    function get_energy_ignition_mode(this) result(value)
        class(energy_ignition_control), intent(in) :: this
        character(len=mode_length) :: value
        value = this%mode
    end function get_energy_ignition_mode
    real(dp) function get_energy_ignition_start_time(this)
        class(energy_ignition_control), intent(in) :: this
        get_energy_ignition_start_time = this%start_time
    end function get_energy_ignition_start_time
    real(dp) function get_energy_ignition_duration(this)
        class(energy_ignition_control), intent(in) :: this
        get_energy_ignition_duration = this%duration
    end function get_energy_ignition_duration
    real(dp) function get_energy_ignition_total_energy(this)
        class(energy_ignition_control), intent(in) :: this
        get_energy_ignition_total_energy = this%total_energy
    end function get_energy_ignition_total_energy
    function get_energy_ignition_center(this) result(value)
        class(energy_ignition_control), intent(in) :: this
        real(dp), dimension(3) :: value
        value = this%center
    end function get_energy_ignition_center
    real(dp) function get_energy_ignition_radius(this)
        class(energy_ignition_control), intent(in) :: this
        get_energy_ignition_radius = this%radius
    end function get_energy_ignition_radius

    function get_turbulence_forcing_mode(this) result(value)
        class(turbulence_forcing_control), intent(in) :: this
        character(len=mode_length) :: value
        value = this%mode
    end function get_turbulence_forcing_mode
    real(dp) function get_turbulence_forcing_start_time(this)
        class(turbulence_forcing_control), intent(in) :: this
        get_turbulence_forcing_start_time = this%start_time
    end function get_turbulence_forcing_start_time
    real(dp) function get_turbulence_forcing_duration(this)
        class(turbulence_forcing_control), intent(in) :: this
        get_turbulence_forcing_duration = this%duration
    end function get_turbulence_forcing_duration
    real(dp) function get_turbulence_forcing_ramp_time(this)
        class(turbulence_forcing_control), intent(in) :: this
        get_turbulence_forcing_ramp_time = this%ramp_time
    end function get_turbulence_forcing_ramp_time
    real(dp) function get_turbulence_forcing_acceleration_amplitude(this)
        class(turbulence_forcing_control), intent(in) :: this
        get_turbulence_forcing_acceleration_amplitude = this%acceleration_amplitude
    end function get_turbulence_forcing_acceleration_amplitude
    real(dp) function get_turbulence_forcing_integral_scale(this)
        class(turbulence_forcing_control), intent(in) :: this
        get_turbulence_forcing_integral_scale = this%integral_scale
    end function get_turbulence_forcing_integral_scale
    real(dp) function get_turbulence_forcing_update_interval(this)
        class(turbulence_forcing_control), intent(in) :: this
        get_turbulence_forcing_update_interval = this%update_interval
    end function get_turbulence_forcing_update_interval
    function get_turbulence_forcing_center(this) result(value)
        class(turbulence_forcing_control), intent(in) :: this
        real(dp), dimension(3) :: value
        value = this%center
    end function get_turbulence_forcing_center
    real(dp) function get_turbulence_forcing_radius(this)
        class(turbulence_forcing_control), intent(in) :: this
        get_turbulence_forcing_radius = this%radius
    end function get_turbulence_forcing_radius
    integer function get_turbulence_forcing_random_seed(this)
        class(turbulence_forcing_control), intent(in) :: this
        get_turbulence_forcing_random_seed = this%random_seed
    end function get_turbulence_forcing_random_seed


    subroutine write_problem_controls_log(this, log_unit)
        class(problem_controls), intent(in) :: this
        integer, intent(in) :: log_unit
        real(dp), dimension(3) :: center

        write(log_unit,'(A)') &
            '*************************************************************************************'
        write(log_unit,'(A)') 'Problem controls setup:'
        write(log_unit,'(2A)') ' Flame stabilization mode: ', &
            trim(this%flame_stabilization%get_mode())
        write(log_unit,'(A,ES12.4)') '  time_delay [s]: ', &
            this%flame_stabilization%get_time_delay()
        write(log_unit,'(A,ES12.4)') '  time_track [s]: ', &
            this%flame_stabilization%get_time_track()
        write(log_unit,'(A,ES12.4)') '  time_control [s]: ', &
            this%flame_stabilization%get_time_control()
        write(log_unit,'(A,ES12.4)') '  response_settle_time [s]: ', &
            this%flame_stabilization%get_response_settle_time()
        write(log_unit,'(A,ES12.4)') '  measurement_max_duration [s]: ', &
            this%flame_stabilization%get_measurement_max_duration()
        write(log_unit,'(A,ES12.4)') '  measurement_min_displacement_cells: ', &
            this%flame_stabilization%get_measurement_min_displacement_cells()

        center = this%energy_ignition%get_center()
        write(log_unit,'(2A)') ' Energy ignition mode: ', &
            trim(this%energy_ignition%get_mode())
        if (this%energy_ignition%is_enabled()) then
            write(log_unit,'(A,ES12.4)') '  start_time [s]: ', &
                this%energy_ignition%get_start_time()
            write(log_unit,'(A,ES12.4)') '  duration [s]: ', &
                this%energy_ignition%get_duration()
            write(log_unit,'(A,ES12.4)') '  total_energy [J]: ', &
                this%energy_ignition%get_total_energy()
            write(log_unit,'(A,3ES12.4)') '  center: ', center
            write(log_unit,'(A,ES12.4)') '  radius [m]: ', &
                this%energy_ignition%get_radius()
        end if

        center = this%turbulence_forcing%get_center()
        write(log_unit,'(2A)') ' Turbulence forcing mode: ', &
            trim(this%turbulence_forcing%get_mode())
        if (this%turbulence_forcing%is_enabled()) then
            write(log_unit,'(A,ES12.4)') '  start_time [s]: ', &
                this%turbulence_forcing%get_start_time()
            write(log_unit,'(A,ES12.4)') '  duration [s]: ', &
                this%turbulence_forcing%get_duration()
            write(log_unit,'(A,ES12.4)') '  ramp_time [s]: ', &
                this%turbulence_forcing%get_ramp_time()
            write(log_unit,'(A,ES12.4)') '  acceleration_amplitude [m/s2]: ', &
                this%turbulence_forcing%get_acceleration_amplitude()
            write(log_unit,'(A,ES12.4)') '  integral_scale [m]: ', &
                this%turbulence_forcing%get_integral_scale()
            write(log_unit,'(A,ES12.4)') '  update_interval [s]: ', &
                this%turbulence_forcing%get_update_interval()
            write(log_unit,'(A,3ES12.4)') '  center: ', center
            write(log_unit,'(A,ES12.4)') '  radius [m; 0=whole domain]: ', &
                this%turbulence_forcing%get_radius()
            write(log_unit,'(A,I0)') '  random_seed: ', &
                this%turbulence_forcing%get_random_seed()
        end if
        write(log_unit,'(A)') &
            '*************************************************************************************'
    end subroutine write_problem_controls_log


    pure function normalize_mode(mode) result(normalized)
        character(len=*), intent(in) :: mode
        character(len=mode_length) :: normalized
        integer :: i, code

        normalized = adjustl(mode)
        do i = 1, len_trim(normalized)
            code = iachar(normalized(i:i))
            if (code >= iachar('A') .and. code <= iachar('Z')) then
                normalized(i:i) = achar(code + iachar('a') - iachar('A'))
            end if
        end do
    end function normalize_mode

end module problem_control_class
