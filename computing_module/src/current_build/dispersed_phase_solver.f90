module dispersed_phase_solver_class

    use kind_parameters, only: dp
    use data_manager_class, only: data_manager
    use solver_options_class, only: particles_phase
    use lagrangian_particles_solver_class, only: lagrangian_particles_solver, &
        lagrangian_particles_solver_c
    use continuous_particles_solver_class, only: continuous_particles_solver, &
        continuous_particles_solver_c

    implicit none

    private
    public :: dispersed_phase_solver, dispersed_phase_solver_c

    !> Runtime facade for one dispersed material phase.
    !>
    !> Both backends are compiled, but construction allocates exactly one.
    !> Parent gas solvers therefore use one stable API and contain no
    !> comment/uncomment model selection.
    type :: dispersed_phase_solver
        private
        character(len=24) :: model = ''
        type(lagrangian_particles_solver), allocatable :: lagrangian
        type(continuous_particles_solver), allocatable :: continuous
    contains
        procedure :: set_initial_distributions
        procedure :: advance
        procedure :: get_model
        procedure :: is_lagrangian
        procedure :: is_continuous
    end type dispersed_phase_solver

    interface dispersed_phase_solver_c
        module procedure constructor
    end interface dispersed_phase_solver_c

contains

    type(dispersed_phase_solver) function constructor( &
        manager, parameters, phase_number)
        type(data_manager), intent(inout) :: manager
        type(particles_phase), intent(in) :: parameters
        integer, intent(in) :: phase_number
        character(len=24) :: model

        ! Model aliases are accepted at the interface but stored canonically.
        model = lowercase(trim(parameters%model))
        select case (model)
        case ('lagrangian', 'parcel')
            constructor%model = 'lagrangian'
            allocate(constructor%lagrangian)
            constructor%lagrangian = lagrangian_particles_solver_c( &
                manager, parameters, phase_number)
        case ('continuous', 'eulerian')
            constructor%model = 'continuous'
            allocate(constructor%continuous)
            constructor%continuous = continuous_particles_solver_c( &
                manager, parameters, phase_number)
        case default
            error stop 'dispersed phase: unknown model; use lagrangian or continuous'
        end select
    end function constructor


    subroutine set_initial_distributions(this)
        class(dispersed_phase_solver), intent(inout) :: this

        select case (trim(this%model))
        case ('lagrangian')
            call this%lagrangian%set_initial_distributions()
        case ('continuous')
            call this%continuous%set_initial_distributions()
        case default
            error stop 'dispersed phase: uninitialized facade'
        end select
    end subroutine set_initial_distributions


    subroutine advance(this, time_step, time)
        class(dispersed_phase_solver), intent(inout) :: this
        real(dp), intent(in) :: time_step, time

        select case (trim(this%model))
        case ('lagrangian')
            call this%lagrangian%advance(time_step, time)
        case ('continuous')
            call this%continuous%advance(time_step, time)
        case default
            error stop 'dispersed phase: uninitialized facade'
        end select
    end subroutine advance


    pure function get_model(this) result(model)
        class(dispersed_phase_solver), intent(in) :: this
        character(len=24) :: model
        model = this%model
    end function get_model


    pure function is_lagrangian(this) result(value)
        class(dispersed_phase_solver), intent(in) :: this
        logical :: value
        value = trim(this%model) == 'lagrangian'
    end function is_lagrangian


    pure function is_continuous(this) result(value)
        class(dispersed_phase_solver), intent(in) :: this
        logical :: value
        value = trim(this%model) == 'continuous'
    end function is_continuous


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

end module dispersed_phase_solver_class
