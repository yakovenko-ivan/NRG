module run_control_class

    use, intrinsic :: iso_fortran_env, only: int64
    use kind_parameters, only: dp
    use global_data

    implicit none

    private
    public :: run_control, run_control_c

    integer, parameter :: mode_length = 24
    integer, parameter :: reason_length = 32

    type :: run_control
        private

        character(len=mode_length) :: termination_mode = 'none'
        real(dp) :: final_time = 1.0e300_dp
        real(dp) :: wall_time_limit = 1.0e300_dp
        real(dp) :: wall_time_reserve = 0.0_dp

        integer(int64) :: start_count = 0_int64
        integer(int64) :: clock_rate = 0_int64
        integer(int64) :: clock_max = 0_int64
        logical :: clock_started = .false.

    contains
        procedure, private :: read_properties
        procedure, private :: write_properties
        procedure, private :: set_properties

        procedure :: start
        procedure :: check_termination
        procedure :: get_elapsed_wall_time
        procedure :: get_final_time
        procedure :: get_wall_time_limit
        procedure :: get_wall_time_reserve
        procedure :: get_termination_mode
        procedure :: limits_simulation_time
        procedure :: limits_wall_time
        procedure :: cap_time_step
        procedure :: write_log
        procedure :: write_termination_log
    end type run_control

    interface run_control_c
        module procedure constructor
        module procedure constructor_file
    end interface run_control_c

contains

    type(run_control) function constructor(termination_mode, final_time, &
            wall_time_limit, wall_time_reserve)

        character(len=*), intent(in) :: termination_mode
        real(dp), intent(in), optional :: final_time
        real(dp), intent(in), optional :: wall_time_limit
        real(dp), intent(in), optional :: wall_time_reserve

        real(dp) :: final_time_set
        real(dp) :: wall_time_limit_set
        real(dp) :: wall_time_reserve_set
        integer :: io_unit

        final_time_set = 1.0e300_dp
        wall_time_limit_set = 1.0e300_dp
        wall_time_reserve_set = 0.0_dp

        if (present(final_time)) final_time_set = final_time
        if (present(wall_time_limit)) wall_time_limit_set = wall_time_limit
        if (present(wall_time_reserve)) wall_time_reserve_set = wall_time_reserve

        call constructor%set_properties(termination_mode, final_time_set, &
            wall_time_limit_set, wall_time_reserve_set)

        open(newunit=io_unit, file=run_control_data_file_name, status='replace', &
            form='formatted', delim='quote')
        call constructor%write_properties(io_unit)
        close(io_unit)
    end function constructor


    type(run_control) function constructor_file()

        integer :: io_unit, io_status
        logical :: file_exists

        inquire(file=run_control_data_file_name, exist=file_exists)

        if (.not. file_exists) then
            ! Backward compatibility: legacy problem setups do not contain a
            ! run-control file.  In that case no new termination criterion is
            ! imposed and the old data-IO/post-processing mechanisms can run.
            call constructor_file%set_properties( &
                'none', 1.0e300_dp, 1.0e300_dp, 0.0_dp)
            return
        end if

        open(newunit=io_unit, file=run_control_data_file_name, status='old', &
            form='formatted', action='read', iostat=io_status)
        if (io_status /= 0) then
            error stop 'run_control: unable to open run_control.dat'
        end if

        call constructor_file%read_properties(io_unit)
        close(io_unit)
    end function constructor_file


    subroutine set_properties(this, termination_mode, final_time, &
            wall_time_limit, wall_time_reserve)

        class(run_control), intent(inout) :: this
        character(len=*), intent(in) :: termination_mode
        real(dp), intent(in) :: final_time
        real(dp), intent(in) :: wall_time_limit
        real(dp), intent(in) :: wall_time_reserve

        character(len=mode_length) :: normalized_mode

        normalized_mode = normalize_mode(termination_mode)

        select case (trim(normalized_mode))
        case ('none')
        case ('simulation_time')
            if (final_time <= 0.0_dp) then
                error stop 'run_control: final simulation time must be positive'
            end if
        case ('wall_time')
            call validate_wall_limit(wall_time_limit, wall_time_reserve)
        case ('either')
            if (final_time <= 0.0_dp) then
                error stop 'run_control: final simulation time must be positive'
            end if
            call validate_wall_limit(wall_time_limit, wall_time_reserve)
        case default
            error stop 'run_control: unknown termination mode'
        end select

        this%termination_mode = normalized_mode
        this%final_time = final_time
        this%wall_time_limit = wall_time_limit
        this%wall_time_reserve = wall_time_reserve

        this%start_count = 0_int64
        this%clock_rate = 0_int64
        this%clock_max = 0_int64
        this%clock_started = .false.
    end subroutine set_properties


    subroutine validate_wall_limit(wall_time_limit, wall_time_reserve)
        real(dp), intent(in) :: wall_time_limit, wall_time_reserve

        if (wall_time_limit <= 0.0_dp) then
            error stop 'run_control: wall-time limit must be positive'
        end if
        if (wall_time_reserve < 0.0_dp) then
            error stop 'run_control: wall-time reserve cannot be negative'
        end if
        if (wall_time_reserve >= wall_time_limit) then
            error stop 'run_control: wall-time reserve must be smaller than limit'
        end if
    end subroutine validate_wall_limit


    subroutine write_properties(this, io_unit)
        class(run_control), intent(in) :: this
        integer, intent(in) :: io_unit

        character(len=mode_length) :: termination_mode
        real(dp) :: final_time, wall_time_limit, wall_time_reserve

        namelist /run_control_parameters/ termination_mode, final_time, &
            wall_time_limit, wall_time_reserve

        termination_mode = this%termination_mode
        final_time = this%final_time
        wall_time_limit = this%wall_time_limit
        wall_time_reserve = this%wall_time_reserve

        write(unit=io_unit, nml=run_control_parameters)
    end subroutine write_properties


    subroutine read_properties(this, io_unit)
        class(run_control), intent(inout) :: this
        integer, intent(in) :: io_unit

        character(len=mode_length) :: termination_mode
        real(dp) :: final_time, wall_time_limit, wall_time_reserve
        integer :: io_status

        namelist /run_control_parameters/ termination_mode, final_time, &
            wall_time_limit, wall_time_reserve

        ! Defaults also make hand-written files tolerant of omitted parameters
        ! that are irrelevant for the selected mode.
        termination_mode = 'none'
        final_time = 1.0e300_dp
        wall_time_limit = 1.0e300_dp
        wall_time_reserve = 0.0_dp

        read(unit=io_unit, nml=run_control_parameters, iostat=io_status)
        if (io_status /= 0) then
            error stop 'run_control: invalid run_control_parameters namelist'
        end if

        call this%set_properties(termination_mode, final_time, &
            wall_time_limit, wall_time_reserve)
    end subroutine read_properties


    subroutine start(this)
        class(run_control), intent(inout) :: this

        call system_clock(this%start_count, this%clock_rate, this%clock_max)

        if (this%clock_rate <= 0_int64) then
            error stop 'run_control: SYSTEM_CLOCK has invalid count rate'
        end if

        this%clock_started = .true.
    end subroutine start


    subroutine check_termination(this, simulation_time, terminate, reason, &
            restart_required)
        class(run_control), intent(inout) :: this
        real(dp), intent(in) :: simulation_time
        logical, intent(out) :: terminate
        character(len=*), intent(out), optional :: reason
        logical, intent(out), optional :: restart_required

        real(dp) :: elapsed_wall_time, wall_trigger_time
        real(dp) :: time_tolerance
        logical :: simulation_limit_reached, wall_limit_reached
        character(len=reason_length) :: local_reason

        if (.not. this%clock_started) call this%start()

        simulation_limit_reached = .false.
        wall_limit_reached = .false.

        if (this%limits_simulation_time()) then
            time_tolerance = 100.0_dp*epsilon(1.0_dp)* &
                max(1.0_dp, abs(this%final_time))
            simulation_limit_reached = &
                simulation_time >= this%final_time - time_tolerance
        end if

        if (this%limits_wall_time()) then
            elapsed_wall_time = this%get_elapsed_wall_time()
            wall_trigger_time = this%wall_time_limit - this%wall_time_reserve
            wall_limit_reached = elapsed_wall_time >= wall_trigger_time
        end if

        terminate = simulation_limit_reached .or. wall_limit_reached

        if (simulation_limit_reached) then
            local_reason = 'final_simulation_time'
        else if (wall_limit_reached) then
            local_reason = 'wall_time_limit'
        else
            local_reason = 'none'
        end if

        if (present(reason)) reason = trim(local_reason)

        ! A restart checkpoint is required only when a wall-clock limit stops
        ! an otherwise unfinished calculation.  If the physical final time is
        ! reached in the same step, the calculation is complete and no restart
        ! checkpoint is needed.
        if (present(restart_required)) then
            restart_required = wall_limit_reached .and. &
                .not. simulation_limit_reached
        end if
    end subroutine check_termination


    function get_elapsed_wall_time(this) result(value)
        class(run_control), intent(in) :: this
        real(dp) :: value

        integer(int64) :: current_count, elapsed_count

        if (.not. this%clock_started) then
            value = 0.0_dp
            return
        end if

        call system_clock(current_count)

        if (current_count >= this%start_count) then
            elapsed_count = current_count - this%start_count
        else if (this%clock_max > 0_int64) then
            elapsed_count = (this%clock_max - this%start_count) + &
                current_count + 1_int64
        else
            error stop 'run_control: SYSTEM_CLOCK counter moved backwards'
        end if

        value = real(elapsed_count, dp)/real(this%clock_rate, dp)
    end function get_elapsed_wall_time


    pure function limits_simulation_time(this) result(value)
        class(run_control), intent(in) :: this
        logical :: value

        value = trim(this%termination_mode) == 'simulation_time' .or. &
            trim(this%termination_mode) == 'either'
    end function limits_simulation_time


    pure function limits_wall_time(this) result(value)
        class(run_control), intent(in) :: this
        logical :: value

        value = trim(this%termination_mode) == 'wall_time' .or. &
            trim(this%termination_mode) == 'either'
    end function limits_wall_time


    pure function cap_time_step(this, current_time, proposed_time_step) result(value)
        class(run_control), intent(in) :: this
        real(dp), intent(in) :: current_time, proposed_time_step
        real(dp) :: value, remaining_time

        value = proposed_time_step
        if (.not. this%limits_simulation_time()) return

        remaining_time = max(this%final_time - current_time, 0.0_dp)
        value = min(value, remaining_time)
    end function cap_time_step


    pure function get_final_time(this) result(value)
        class(run_control), intent(in) :: this
        real(dp) :: value
        value = this%final_time
    end function get_final_time


    pure function get_wall_time_limit(this) result(value)
        class(run_control), intent(in) :: this
        real(dp) :: value
        value = this%wall_time_limit
    end function get_wall_time_limit


    pure function get_wall_time_reserve(this) result(value)
        class(run_control), intent(in) :: this
        real(dp) :: value
        value = this%wall_time_reserve
    end function get_wall_time_reserve


    pure function get_termination_mode(this) result(value)
        class(run_control), intent(in) :: this
        character(len=mode_length) :: value
        value = this%termination_mode
    end function get_termination_mode


    subroutine write_log(this, log_unit)
        class(run_control), intent(in) :: this
        integer, intent(in) :: log_unit

        write(log_unit, '(A)') repeat('*', 84)
        write(log_unit, '(A,A)') ' Termination mode            : ', &
            trim(this%termination_mode)

        if (this%limits_simulation_time()) then
            write(log_unit, '(A,ES14.6,A)') ' Final simulation time       : ', &
                this%final_time, ' s'
        else
            write(log_unit, '(A)') ' Final simulation time       : disabled'
        end if

        if (this%limits_wall_time()) then
            write(log_unit, '(A,ES14.6,A)') ' Wall-time limit             : ', &
                this%wall_time_limit, ' s'
            write(log_unit, '(A,ES14.6,A)') ' Wall-time reserve           : ', &
                this%wall_time_reserve, ' s'
        else
            write(log_unit, '(A)') ' Wall-time limit             : disabled'
        end if

        write(log_unit, '(A)') repeat('*', 84)
    end subroutine write_log


    subroutine write_termination_log(this, log_unit, simulation_time, reason)
        class(run_control), intent(in) :: this
        integer, intent(in) :: log_unit
        real(dp), intent(in) :: simulation_time
        character(len=*), intent(in) :: reason

        write(log_unit, '(A)') repeat('*', 84)
        write(log_unit, '(A)') ' Calculation terminated normally.'
        write(log_unit, '(A,A)') ' Reason                      : ', trim(reason)
        write(log_unit, '(A,ES14.6,A)') ' Final simulation time       : ', &
            simulation_time, ' s'
        write(log_unit, '(A,ES14.6,A)') ' Elapsed wall time           : ', &
            this%get_elapsed_wall_time(), ' s'
        write(log_unit, '(A)') repeat('*', 84)
    end subroutine write_termination_log


    pure function normalize_mode(mode) result(value)
        character(len=*), intent(in) :: mode
        character(len=mode_length) :: value
        character(len=len(mode)) :: lower

        lower = lowercase(adjustl(mode))

        select case (trim(lower))
        case ('none', 'disabled')
            value = 'none'
        case ('simulation_time', 'physical_time', 'final_time')
            value = 'simulation_time'
        case ('wall_time', 'wall_clock', 'wall_clock_time')
            value = 'wall_time'
        case ('either', 'both')
            value = 'either'
        case default
            value = trim(lower)
        end select
    end function normalize_mode


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

end module run_control_class
