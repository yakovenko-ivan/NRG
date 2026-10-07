program transport_mixture_diffusion_benchmark
    use, intrinsic :: iso_fortran_env, only: int64
    use kind_parameters, only: dp
    use global_data, only: P_atm
    use chemical_properties_class, only: chemical_properties, chemical_properties_c
    use thermophysical_properties_class, only: thermophysical_properties, &
        thermophysical_properties_c

    implicit none

    integer, parameter :: number_of_states = 8
    integer, parameter :: default_calls_per_repeat = 1000000
    integer, parameter :: default_repeats = 7
    integer, parameter :: default_warmup_calls = 20000

    type(chemical_properties) :: chemistry
    type(thermophysical_properties) :: thermophysics

    real(dp), allocatable :: mass_fractions(:,:)
    real(dp), allocatable :: diffusion(:)
    real(dp), dimension(number_of_states) :: temperature, pressure
    real(dp), dimension(number_of_states) :: mixture_molar_mass_state
    real(dp), allocatable :: wall_times(:), cpu_times(:)

    integer :: species_number
    integer :: calls_per_repeat, repeats, warmup_calls
    integer :: repeat, call_index, state_index
    integer :: output_unit
    integer(int64) :: clock_start, clock_finish, clock_rate
    real(dp) :: cpu_start, cpu_finish, elapsed_wall, elapsed_cpu
    real(dp) :: checksum, reference_checksum, diffusivity_weight

    calls_per_repeat = default_calls_per_repeat
    repeats = default_repeats
    warmup_calls = default_warmup_calls
    call parse_command_line(calls_per_repeat, repeats, warmup_calls)

    ! Initialize directly from the committed raw KEROMNES input files.
    ! Run this executable with package_interface/ as the working directory so
    ! that ./task_setup/chemical_mechanisms and ./task_setup/thermophysical_data
    ! resolve exactly as they do in a normal package-interface setup generator.
    chemistry = chemical_properties_c( &
        chemical_mechanism_file_name='KEROMNES.txt', &
        default_enhanced_efficiencies=1.0_dp, &
        E_act_units='cal.mol')

    thermophysics = thermophysical_properties_c( &
        chemistry=chemistry, &
        thermo_data_file_name='KEROMNES_THERMO.txt', &
        transport_data_file_name='KEROMNES_TRANSDATA.txt', &
        molar_masses_data_file_name='molar_masses.dat')

    species_number = chemistry%species_number
    if (species_number <= 1) then
        error stop 'Transport microbenchmark requires at least two species.'
    end if

    allocate(mass_fractions(species_number, number_of_states))
    allocate(diffusion(species_number))
    allocate(wall_times(repeats), cpu_times(repeats))

    call build_h2_air_states(chemistry, thermophysics, mass_fractions, &
        temperature, pressure)

    do state_index = 1, number_of_states
        mixture_molar_mass_state(state_index) = &
            thermophysics%mixture_molar_mass_from_mass_fractions( &
                mass_fractions(:,state_index))
    end do

    ! A compact numerical signature independent of the benchmark repetition
    ! count. This makes baseline/optimized output comparison straightforward.
    reference_checksum = 0.0_dp
    do state_index = 1, number_of_states
        call thermophysics%mixture_averaged_diffusion_coefficients( &
            temperature(state_index), pressure(state_index), &
            mass_fractions(:,state_index), &
            mixture_molar_mass_state(state_index), diffusion)
        do call_index = 1, species_number
            diffusivity_weight = real(call_index,dp)
            reference_checksum = reference_checksum + &
                diffusivity_weight*diffusion(call_index)
        end do
    end do

    ! Warm-up removes one-time instruction/data-cache effects and lets the CPU
    ! reach the same execution regime before timed samples begin.
    checksum = 0.0_dp
    do call_index = 1, warmup_calls
        state_index = 1 + mod(call_index-1, number_of_states)
        call thermophysics%mixture_averaged_diffusion_coefficients( &
            temperature(state_index), pressure(state_index), &
            mass_fractions(:,state_index), &
            mixture_molar_mass_state(state_index), diffusion)
        checksum = checksum + diffusion(1) + diffusion(species_number)
    end do

    call system_clock(count_rate=clock_rate)
    if (clock_rate <= 0_int64) then
        error stop 'system_clock returned a non-positive count rate.'
    end if

    open(newunit=output_unit, file='transport_microbenchmark.csv', &
        status='replace', action='write', form='formatted')
    write(output_unit,'(A)') &
        'repeat,wall_seconds,cpu_seconds,ns_per_call,reference_checksum'

    checksum = 0.0_dp
    do repeat = 1, repeats
        call cpu_time(cpu_start)
        call system_clock(clock_start)

        do call_index = 1, calls_per_repeat
            state_index = 1 + mod(call_index-1, number_of_states)
            call thermophysics%mixture_averaged_diffusion_coefficients( &
                temperature(state_index), pressure(state_index), &
                mass_fractions(:,state_index), &
                mixture_molar_mass_state(state_index), diffusion)

            ! Observable dependence on each call prevents dead-code elimination.
            checksum = checksum + diffusion(1) + diffusion(species_number)
        end do

        call system_clock(clock_finish)
        call cpu_time(cpu_finish)

        elapsed_wall = real(clock_finish-clock_start,dp)/real(clock_rate,dp)
        elapsed_cpu = cpu_finish-cpu_start

        wall_times(repeat) = elapsed_wall
        cpu_times(repeat) = elapsed_cpu

        write(output_unit,'(I0,",",ES24.16,",",ES24.16,",",ES24.16,",",ES24.16)') &
            repeat, elapsed_wall, elapsed_cpu, &
            1.0e9_dp*elapsed_wall/real(calls_per_repeat,dp), reference_checksum

        write(*,'(A,I0,A,F10.6,A,F10.6,A,F12.3)') &
            'repeat ', repeat, ': wall=', elapsed_wall, &
            ' s, cpu=', elapsed_cpu, &
            ' s, ns/call=', 1.0e9_dp*elapsed_wall/real(calls_per_repeat,dp)
    end do
    close(output_unit)

    write(*,'(/,A)') 'NRG mixture-diffusion transport microbenchmark'
    write(*,'(A,I0)') 'species:             ', species_number
    write(*,'(A,I0)') 'representative states:', number_of_states
    write(*,'(A,I0)') 'calls per repeat:    ', calls_per_repeat
    write(*,'(A,I0)') 'warm-up calls:       ', warmup_calls
    write(*,'(A,I0)') 'repeats:             ', repeats
    write(*,'(A,ES24.16)') 'reference checksum:  ', reference_checksum
    write(*,'(A,ES24.16)') 'runtime checksum:    ', checksum
    write(*,'(A,F12.3)') 'best wall ns/call:   ', &
        1.0e9_dp*minval(wall_times)/real(calls_per_repeat,dp)
    write(*,'(A,F12.3)') 'median wall ns/call: ', &
        1.0e9_dp*median(wall_times)/real(calls_per_repeat,dp)
    write(*,'(A,F12.3)') 'median CPU ns/call:  ', &
        1.0e9_dp*median(cpu_times)/real(calls_per_repeat,dp)
    write(*,'(A)') 'CSV: transport_microbenchmark.csv'

contains

    subroutine parse_command_line(calls, repeat_count, warmup)
        integer, intent(inout) :: calls, repeat_count, warmup

        integer :: arg_index, ios
        character(len=256) :: arg

        do arg_index = 1, command_argument_count()
            call get_command_argument(arg_index, arg)

            if (index(trim(arg),'--calls=') == 1) then
                read(arg(len('--calls=')+1:),*,iostat=ios) calls
                if (ios /= 0 .or. calls <= 0) &
                    error stop 'Invalid --calls=N argument.'
            else if (index(trim(arg),'--repeats=') == 1) then
                read(arg(len('--repeats=')+1:),*,iostat=ios) repeat_count
                if (ios /= 0 .or. repeat_count <= 0) &
                    error stop 'Invalid --repeats=N argument.'
            else if (index(trim(arg),'--warmup=') == 1) then
                read(arg(len('--warmup=')+1:),*,iostat=ios) warmup
                if (ios /= 0 .or. warmup < 0) &
                    error stop 'Invalid --warmup=N argument.'
            else if (trim(arg) == '--help' .or. trim(arg) == '-h') then
                write(*,'(A)') 'usage: package_interface_transport_mixture_diffusion_benchmark ' // &
                    '[--calls=N] [--repeats=N] [--warmup=N]'
                stop
            else
                write(*,'(A,A)') 'Unknown argument: ', trim(arg)
                error stop 'Unsupported transport benchmark argument.'
            end if
        end do
    end subroutine parse_command_line


    subroutine build_h2_air_states(chem, thermo, Y_states, T_states, p_states)
        type(chemical_properties), intent(inout) :: chem
        type(thermophysical_properties), intent(in) :: thermo
        real(dp), dimension(:,:), intent(out) :: Y_states
        real(dp), dimension(:), intent(out) :: T_states, p_states

        real(dp), parameter :: initial_h2 = 0.08_dp
        real(dp), dimension(number_of_states), parameter :: progress = &
            (/0.0_dp, 0.10_dp, 0.25_dp, 0.45_dp, &
              0.65_dp, 0.80_dp, 0.95_dp, 1.0_dp/)
        real(dp), dimension(number_of_states), parameter :: flame_temperature = &
            (/300.0_dp, 500.0_dp, 800.0_dp, 1100.0_dp, &
              1400.0_dp, 1700.0_dp, 2000.0_dp, 2300.0_dp/)

        real(dp), allocatable :: X(:)
        real(dp) :: x_o2_initial, x_n2_initial, normalization
        real(dp) :: mass_normalization
        integer :: state, specie
        integer :: h2_index, o2_index, n2_index, h2o_index

        if (size(Y_states,1) /= chem%species_number .or. &
            size(Y_states,2) /= number_of_states) then
            error stop 'Unexpected state-array dimensions.'
        end if

        h2_index = chem%get_chemical_specie_index('H2')
        o2_index = chem%get_chemical_specie_index('O2')
        n2_index = chem%get_chemical_specie_index('N2')
        h2o_index = chem%get_chemical_specie_index('H2O')

        if (min(h2_index,o2_index,n2_index,h2o_index) <= 0) then
            error stop 'Benchmark requires H2, O2, N2, and H2O in the mechanism.'
        end if

        allocate(X(chem%species_number))
        x_o2_initial = (1.0_dp-initial_h2)/4.762_dp
        x_n2_initial = 3.762_dp*x_o2_initial

        do state = 1, number_of_states
            X = 0.0_dp

            ! Major-species path from fresh 8% H2/air toward fully consumed H2:
            !        2 H2 + O2 -> 2 H2O
            ! It is not used as a flame solution; it simply supplies realistic
            ! composition/temperature variation to the production transport kernel.
            X(h2_index) = initial_h2*(1.0_dp-progress(state))
            X(o2_index) = x_o2_initial - 0.5_dp*initial_h2*progress(state)
            X(n2_index) = x_n2_initial
            X(h2o_index) = initial_h2*progress(state)

            normalization = sum(X)
            if (normalization <= 0.0_dp) &
                error stop 'Invalid benchmark mole-fraction state.'
            X = X/normalization

            mass_normalization = 0.0_dp
            do specie = 1, chem%species_number
                mass_normalization = mass_normalization + &
                    X(specie)*thermo%molar_masses(specie)
            end do
            if (mass_normalization <= 0.0_dp) &
                error stop 'Invalid benchmark mass normalization.'

            do specie = 1, chem%species_number
                Y_states(specie,state) = &
                    X(specie)*thermo%molar_masses(specie)/mass_normalization
            end do

            T_states(state) = flame_temperature(state)
            p_states(state) = P_atm
        end do

        deallocate(X)
    end subroutine build_h2_air_states


    real(dp) function median(values) result(value)
        real(dp), dimension(:), intent(in) :: values

        real(dp), allocatable :: sorted(:)
        real(dp) :: temporary
        integer :: i, j, count

        count = size(values)
        if (count <= 0) error stop 'median() received an empty array.'

        allocate(sorted(count))
        sorted = values

        do i = 2, count
            temporary = sorted(i)
            j = i-1
            do while (j >= 1)
                if (sorted(j) <= temporary) exit
                sorted(j+1) = sorted(j)
                j = j-1
            end do
            sorted(j+1) = temporary
        end do

        if (mod(count,2) == 1) then
            value = sorted((count+1)/2)
        else
            value = 0.5_dp*(sorted(count/2)+sorted(count/2+1))
        end if

        deallocate(sorted)
    end function median

end program transport_mixture_diffusion_benchmark
