!===============================================================================
! MESH-FREE CHEMICAL-EQUILIBRIUM / KINETICS VALIDATION INTERFACE
!===============================================================================
!
! This package interface deliberately avoids computational_domain, fields, MPI,
! and data_manager.  It instantiates the same NRG chemical/thermophysical data
! classes used by production solvers, the same mesh-free chemical-kinetics core
! used by chemical_kinetics_solver, and the Gibbs-equilibrium calculator.
!
! Modes:
!   tp        fixed-temperature / fixed-pressure ideal-gas Gibbs equilibrium
!   hp        adiabatic constant-pressure equilibrium
!   crossover H2/O2 chain-branching / HO2-stabilization crossover temperature
!   all       execute all three calculations
!
! An optional first command-line argument names a namelist file containing
! /equilibrium_validation/.  Without it, the defaults below are used.
!===============================================================================

program package_interface

    use, intrinsic :: iso_fortran_env, only: error_unit, output_unit

    use kind_parameters, only: dp
    use global_data, only: P_atm
    use chemical_properties_class, only: chemical_properties, &
        chemical_properties_c
    use thermophysical_properties_class, only: thermophysical_properties, &
        thermophysical_properties_c
    use chemical_kinetics_core_class, only: chemical_kinetics_core, &
        chemical_kinetics_core_c, reaction_rate_diagnostics
    use chemical_equilibrium_class, only: chemical_equilibrium, &
        chemical_equilibrium_c, equilibrium_state, equilibrium_diagnostics
    use chemical_kinetics_diagnostics, only: hydrogen_crossover_diagnostics, &
        hydrogen_crossover_temperature

    implicit none

    integer, parameter :: maximum_components = 64

    type(chemical_properties), target :: chemistry
    type(thermophysical_properties), target :: thermophysics
    type(chemical_kinetics_core) :: kinetics_core
    type(chemical_equilibrium) :: equilibrium
    type(equilibrium_state) :: tp_state,hp_state
    type(equilibrium_diagnostics) :: tp_diagnostics,hp_diagnostics
    type(hydrogen_crossover_diagnostics) :: crossover_diagnostics

    character(len=16) :: mode
    character(len=64) :: mechanism_file,thermo_file,transport_file
    character(len=64) :: component_species(maximum_components)
    real(dp) :: component_moles(maximum_components)
    integer :: components_number
    real(dp) :: tp_temperature,initial_temperature,pressure
    real(dp) :: crossover_temperature_min,crossover_temperature_max
    character(len=128) :: output_prefix

    namelist /equilibrium_validation/ mode,mechanism_file,thermo_file, &
        transport_file,components_number,component_species,component_moles, &
        tp_temperature,initial_temperature,pressure, &
        crossover_temperature_min,crossover_temperature_max,output_prefix

    real(dp), allocatable :: initial_moles(:),initial_x(:)
    real(dp) :: crossover_temperature
    character(len=1024) :: config_file
    character(len=512) :: io_message
    integer :: io_unit,io_status,component,specie_index
    logical :: do_tp,do_hp,do_crossover

    ! Defaults form a stoichiometric H2-air validation case for KEROMNES.
    mode = 'all'
    mechanism_file = 'KEROMNES.txt'
    thermo_file = 'KEROMNES_THERMO.txt'
    transport_file = 'KEROMNES_TRANSDATA.txt'
    components_number = 3
    component_species = ''
    component_moles = 0.0_dp
    component_species(1:3) = (/'H2                  ', &
                               'O2                  ', &
                               'N2                  '/)
    component_moles(1:3) = (/2.0_dp,1.0_dp,3.762_dp/)
    tp_temperature = 1500.0_dp
    initial_temperature = 300.0_dp
    pressure = P_atm
    crossover_temperature_min = 400.0_dp
    crossover_temperature_max = 3000.0_dp
    output_prefix = 'equilibrium_validation'

    if (command_argument_count() >= 1) then
        config_file = ''
        call get_command_argument(1,config_file)
        if (len_trim(config_file) > 0) then
            open(newunit=io_unit,file=trim(config_file),status='old', &
                action='read',iostat=io_status,iomsg=io_message)
            if (io_status /= 0) then
                write(error_unit,'(A,A)') &
                    'Unable to open validation configuration: ',trim(config_file)
                write(error_unit,'(A)') trim(io_message)
                error stop 'Chemical-equilibrium validation setup failure'
            end if
            read(io_unit,nml=equilibrium_validation,iostat=io_status, &
                iomsg=io_message)
            close(io_unit)
            if (io_status /= 0) then
                write(error_unit,'(A)') &
                    'Unable to read /equilibrium_validation/ namelist.'
                write(error_unit,'(A)') trim(io_message)
                error stop 'Chemical-equilibrium validation namelist failure'
            end if
        end if
    end if

    call validate_configuration()
    call decode_mode(do_tp,do_hp,do_crossover)

    chemistry = chemical_properties_c( &
        chemical_mechanism_file_name=trim(mechanism_file), &
        default_enhanced_efficiencies=1.0_dp,E_act_units='cal.mol')
    thermophysics = thermophysical_properties_c( &
        chemistry=chemistry,thermo_data_file_name=trim(thermo_file), &
        transport_data_file_name=trim(transport_file), &
        molar_masses_data_file_name='molar_masses.dat')

    kinetics_core = chemical_kinetics_core_c(chemistry,thermophysics)
    equilibrium = chemical_equilibrium_c(chemistry,thermophysics)

    allocate(initial_moles(chemistry%species_number))
    allocate(initial_x(chemistry%species_number))
    initial_moles = 0.0_dp
    do component = 1, components_number
        specie_index = chemistry%get_chemical_specie_index( &
            trim(component_species(component)))
        if (specie_index <= 0 .or. specie_index > chemistry%species_number) then
            write(error_unit,'(A,A)') &
                'Validation composition contains unknown species: ', &
                trim(component_species(component))
            error stop 'Chemical-equilibrium validation species lookup failure'
        end if
        initial_moles(specie_index) = initial_moles(specie_index) + &
            component_moles(component)
    end do
    if (sum(initial_moles) <= tiny(1.0_dp)) then
        error stop 'Chemical-equilibrium validation composition is empty'
    end if
    initial_x = initial_moles/sum(initial_moles)

    call write_input_summary()

    if (do_tp) then
        call equilibrium%equilibrate_tp(tp_temperature,pressure,initial_moles, &
            tp_state,tp_diagnostics)
        call write_equilibrium_state('TP',tp_state,tp_diagnostics)
        call write_equilibrium_stdout('TP',tp_state,tp_diagnostics)
    end if

    if (do_hp) then
        call equilibrium%equilibrate_hp(initial_temperature,pressure, &
            initial_moles,hp_state,hp_diagnostics)
        call write_equilibrium_state('HP',hp_state,hp_diagnostics)
        call write_equilibrium_stdout('HP',hp_state,hp_diagnostics)
    end if

    if (do_crossover) then
        call hydrogen_crossover_temperature(kinetics_core,pressure,initial_x, &
            crossover_temperature,crossover_diagnostics, &
            crossover_temperature_min,crossover_temperature_max)
        call write_crossover_output(crossover_temperature,crossover_diagnostics)
        write(output_unit,'(A,ES16.8)') 'H2/O2 crossover temperature [K]: ', &
            crossover_temperature
        write(output_unit,'(A,L1)') 'Crossover converged: ', &
            crossover_diagnostics%converged
        write(output_unit,'(A,ES16.8)') 'Crossover relative balance residual: ', &
            crossover_diagnostics%relative_balance_residual
    end if

contains

    subroutine validate_configuration()
        integer :: i

        mode = lower_ascii(adjustl(mode))
        if (components_number < 1 .or. components_number > maximum_components) then
            error stop 'components_number is outside the supported range'
        end if
        do i = 1, components_number
            if (len_trim(component_species(i)) == 0) then
                error stop 'A validation component has an empty species name'
            end if
            if (component_moles(i) < 0.0_dp) then
                error stop 'Validation component mole amounts must be nonnegative'
            end if
        end do
        if (pressure <= 0.0_dp) error stop 'Validation pressure must be positive'
        if (tp_temperature <= 0.0_dp) &
            error stop 'TP validation temperature must be positive'
        if (initial_temperature <= 0.0_dp) &
            error stop 'HP initial temperature must be positive'
        if (crossover_temperature_min <= 0.0_dp .or. &
            crossover_temperature_max <= crossover_temperature_min) then
            error stop 'Invalid crossover temperature interval'
        end if
        if (len_trim(output_prefix) == 0) then
            error stop 'Validation output_prefix may not be empty'
        end if
    end subroutine validate_configuration


    subroutine decode_mode(tp,hp,crossover)
        logical, intent(out) :: tp,hp,crossover

        tp = .false.
        hp = .false.
        crossover = .false.
        select case (trim(mode))
        case ('tp')
            tp = .true.
        case ('hp')
            hp = .true.
        case ('crossover')
            crossover = .true.
        case ('all')
            tp = .true.
            hp = .true.
            crossover = .true.
        case default
            error stop "mode must be 'tp', 'hp', 'crossover', or 'all'"
        end select
    end subroutine decode_mode


    subroutine write_input_summary()
        character(len=256) :: file_name
        integer :: unit,i

        file_name = trim(output_prefix)//'_input.csv'
        open(newunit=unit,file=trim(file_name),status='replace',action='write')
        write(unit,'(A)') 'key,value'
        write(unit,'(A,A)') 'mode,',trim(mode)
        write(unit,'(A,A)') 'mechanism_file,',trim(mechanism_file)
        write(unit,'(A,A)') 'thermo_file,',trim(thermo_file)
        write(unit,'(A,A)') 'transport_file,',trim(transport_file)
        write(unit,'(A,ES24.16E3)') 'pressure_Pa,',pressure
        write(unit,'(A,ES24.16E3)') 'tp_temperature_K,',tp_temperature
        write(unit,'(A,ES24.16E3)') 'hp_initial_temperature_K,', &
            initial_temperature
        close(unit)

        file_name = trim(output_prefix)//'_initial_composition.csv'
        open(newunit=unit,file=trim(file_name),status='replace',action='write')
        write(unit,'(A)') 'species,mole_number,mole_fraction'
        do i = 1, chemistry%species_number
            write(unit,'(A,",",ES24.16E3,",",ES24.16E3)') &
                trim(chemistry%species_names(i)),initial_moles(i),initial_x(i)
        end do
        close(unit)
    end subroutine write_input_summary


    subroutine write_equilibrium_state(label,state,diagnostics)
        character(len=*), intent(in) :: label
        type(equilibrium_state), intent(in) :: state
        type(equilibrium_diagnostics), intent(in) :: diagnostics

        character(len=256) :: state_file,diagnostics_file
        character(len=16) :: lower_label
        integer :: unit,specie

        lower_label = lower_ascii(label)
        state_file = trim(output_prefix)//'_'//trim(lower_label)//'_state.csv'
        open(newunit=unit,file=trim(state_file),status='replace',action='write')
        write(unit,'(A)') &
            'species,mole_number,mole_fraction,mass_fraction,molar_mass_kg_per_mol'
        do specie = 1, chemistry%species_number
            write(unit,'(A,4(",",ES24.16E3))') &
                trim(chemistry%species_names(specie)),state%mole_numbers(specie), &
                state%mole_fractions(specie),state%mass_fractions(specie), &
                thermophysics%molar_masses(specie)
        end do
        close(unit)

        diagnostics_file = trim(output_prefix)//'_'//trim(lower_label)// &
            '_diagnostics.csv'
        open(newunit=unit,file=trim(diagnostics_file),status='replace', &
            action='write')
        write(unit,'(A)') 'key,value'
        write(unit,'(A,L1)') 'converged,',diagnostics%converged
        write(unit,'(A,I0)') 'iterations,',diagnostics%iterations
        write(unit,'(A,ES24.16E3)') 'temperature_K,',state%temperature
        write(unit,'(A,ES24.16E3)') 'pressure_Pa,',state%pressure
        write(unit,'(A,ES24.16E3)') 'total_moles,',state%total_moles
        write(unit,'(A,ES24.16E3)') 'max_element_relative_residual,', &
            diagnostics%max_element_relative_residual
        write(unit,'(A,ES24.16E3)') 'mass_relative_residual,', &
            diagnostics%mass_relative_residual
        write(unit,'(A,ES24.16E3)') 'mole_fraction_sum_residual,', &
            diagnostics%mole_fraction_sum_residual
        write(unit,'(A,ES24.16E3)') 'stationarity_residual,', &
            diagnostics%stationarity_residual
        write(unit,'(A,ES24.16E3)') 'enthalpy_relative_residual,', &
            diagnostics%enthalpy_relative_residual
        close(unit)
    end subroutine write_equilibrium_state


    subroutine write_equilibrium_stdout(label,state,diagnostics)
        character(len=*), intent(in) :: label
        type(equilibrium_state), intent(in) :: state
        type(equilibrium_diagnostics), intent(in) :: diagnostics

        write(output_unit,'(A)') ''
        write(output_unit,'(A,A,A)') 'Equilibrium ',trim(label),' result'
        write(output_unit,'(A,L1)') '  converged: ',diagnostics%converged
        write(output_unit,'(A,ES16.8)') '  temperature [K]: ',state%temperature
        write(output_unit,'(A,ES16.8)') '  pressure [Pa]: ',state%pressure
        write(output_unit,'(A,ES16.8)') '  max element residual: ', &
            diagnostics%max_element_relative_residual
        if (trim(lower_ascii(label)) == 'hp') then
            write(output_unit,'(A,ES16.8)') '  enthalpy residual: ', &
                diagnostics%enthalpy_relative_residual
        end if
    end subroutine write_equilibrium_stdout


    subroutine write_crossover_output(temperature,diagnostics)
        real(dp), intent(in) :: temperature
        type(hydrogen_crossover_diagnostics), intent(in) :: diagnostics

        character(len=256) :: summary_file,reactions_file
        integer :: unit,i,reaction

        summary_file = trim(output_prefix)//'_crossover_diagnostics.csv'
        open(newunit=unit,file=trim(summary_file),status='replace',action='write')
        write(unit,'(A)') 'key,value'
        write(unit,'(A,L1)') 'converged,',diagnostics%converged
        write(unit,'(A,I0)') 'iterations,',diagnostics%iterations
        write(unit,'(A,ES24.16E3)') 'temperature_K,',temperature
        write(unit,'(A,ES24.16E3)') 'pressure_Pa,',pressure
        write(unit,'(A,ES24.16E3)') 'branching_rate_coefficient,', &
            diagnostics%branching_rate_coefficient
        write(unit,'(A,ES24.16E3)') 'stabilization_rate_coefficient,', &
            diagnostics%stabilization_rate_coefficient
        write(unit,'(A,ES24.16E3)') 'relative_balance_residual,', &
            diagnostics%relative_balance_residual
        close(unit)

        reactions_file = trim(output_prefix)//'_crossover_reactions.csv'
        open(newunit=unit,file=trim(reactions_file),status='replace',action='write')
        write(unit,'(A)') &
            'role,reaction_index,equation,k_inf,k_0,M_eff,Pr,F_troe,k_forward,k_reverse,reverse_factor'
        do i = 1, size(diagnostics%branching_reactions)
            reaction = diagnostics%branching_reactions(i)
            call write_reaction_diagnostic_row(unit,'branching',reaction, &
                diagnostics%branching_details(i))
        end do
        do i = 1, size(diagnostics%stabilization_reactions)
            reaction = diagnostics%stabilization_reactions(i)
            call write_reaction_diagnostic_row(unit,'stabilization',reaction, &
                diagnostics%stabilization_details(i))
        end do
        close(unit)
    end subroutine write_crossover_output


    subroutine write_reaction_diagnostic_row(unit,role,reaction,details)
        integer, intent(in) :: unit,reaction
        character(len=*), intent(in) :: role
        type(reaction_rate_diagnostics), intent(in) :: details
        character(len=256) :: equation

        equation = kinetics_core%reaction_equation_string(reaction)
        write(unit,'(A,",",I0,",",A,8(",",ES24.16E3))') &
            trim(role),reaction,trim(equation),details%high_pressure_rate, &
            details%low_pressure_rate,details%third_body_concentration, &
            details%reduced_pressure,details%troe_factor, &
            details%forward_rate_coefficient,details%reverse_rate_coefficient, &
            details%reverse_factor
    end subroutine write_reaction_diagnostic_row


    pure function lower_ascii(value) result(lower_value)
        character(len=*), intent(in) :: value
        character(len=len(value)) :: lower_value
        integer :: i,code

        lower_value = value
        do i = 1, len(value)
            code = iachar(lower_value(i:i))
            if (code >= iachar('A') .and. code <= iachar('Z')) then
                lower_value(i:i) = achar(code-iachar('A')+iachar('a'))
            end if
        end do
    end function lower_ascii

end program package_interface
