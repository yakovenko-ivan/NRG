!> Thermally perfect ideal-gas-mixture equation-of-state closures.
!!
!! The module translates between the conservative/material variables advanced by
!! the compressible, FDS low-Mach, and CABARET low-Mach solvers and the derived
!! thermodynamic fields used elsewhere in NRG.  Species thermodynamics are
!! supplied by thermophysical_properties_class (JANAF/NASA tables).
!!
!! Energy convention
!! -----------------
!! The solver stores a shifted specific internal energy
!!
!!   e_s(T,Y) = [e_m(T,X) - h_m(T_ref,X)] / M(Y),
!!
!! and the corresponding sensible enthalpy
!!
!!   h_s(T,Y) = [h_m(T,X) - h_m(T_ref,X)] / M(Y).
!!
!! The same reference enthalpy is used in the inverse-temperature routines of
!! thermophysical_properties_class.  Keeping that convention identical in cell,
!! face, initial-condition, and boundary closures is essential for avoiding
!! artificial temperature jumps at material interfaces.
module table_approximated_real_gas_class

    use kind_parameters, only: dp
    use global_data, only: r_gase_J, T_ref, tables_temperature_ceiling, I_m
    use field_scalar_class, only: field_scalar_cons, field_scalar_flow
    use field_pointers
    use data_manager_class, only: data_manager
    use computational_domain_class, only: computational_domain
    use chemical_properties_class, only: chemical_properties_pointer
    use boundary_conditions_class, only: boundary_conditions_pointer
    use thermophysical_properties_class, only: thermophysical_properties_pointer

    implicit none

    private
    public :: table_approximated_real_gas, table_approximated_real_gas_c

    real(dp), parameter :: eos_temperature_floor = 1.0_dp
    real(dp), parameter :: eos_density_floor = 100.0_dp*tiny(1.0_dp)
    real(dp), parameter :: eos_pressure_floor = 100.0_dp*tiny(1.0_dp)
    real(dp), parameter :: eos_composition_tolerance = 100.0_dp*epsilon(1.0_dp)

    ! Derived cell-centred fields owned by this module and registered with the
    ! data manager.  The objects are module targets because NRG field pointers
    ! refer to them for the lifetime of the simulation.
    type(field_scalar_cons), target :: sound_speed_field
    type(field_scalar_cons), target :: gamma_field
    type(field_scalar_cons), target :: sensible_enthalpy_field
    type(field_scalar_cons), target :: mixture_cp_field
    type(field_scalar_cons), target :: full_enthalpy_field

    ! Derived face-centred field.
    type(field_scalar_flow), target :: gamma_face_field

    type :: table_approximated_real_gas
        ! Cell fields supplied by the data manager.
        type(field_scalar_cons_pointer) :: p, p_stat, p_stat_old
        type(field_scalar_cons_pointer) :: T, e_i, E_f, rho
        type(field_scalar_cons_pointer) :: mix_mol_mass
        type(field_scalar_cons_pointer) :: dp_stat_dt
        type(field_scalar_cons_pointer) :: v_s, gamma, h_s
        type(field_scalar_cons_pointer) :: mixture_cp, h_full
        type(field_vector_cons_pointer) :: v, Y

        ! CABARET face/flux fields.
        type(field_scalar_flow_pointer) :: rho_f, p_f, e_i_f
        type(field_scalar_flow_pointer) :: v_s_f, E_f_f, T_f, gamma_f
        type(field_vector_flow_pointer) :: v_f, Y_f

        type(computational_domain) :: domain
        type(thermophysical_properties_pointer) :: thermo
        type(chemical_properties_pointer) :: chem
        type(boundary_conditions_pointer) :: boundary

    contains
        ! Compressible cell and face closures.
        procedure :: apply_state_equation
        procedure :: apply_state_equation_flow_variables

        ! Initial-condition closures.  The face routine intentionally delegates
        ! to the production face closure so both paths use identical physics.
        procedure :: apply_state_equation_for_initial_conditions
        procedure :: apply_state_equation_flow_variables_for_IC
        procedure :: apply_boundary_conditions_for_initial_conditions

        ! Legacy FDS low-Mach closure with cellwise thermodynamic pressure.
        procedure :: apply_state_equation_low_mach_fds

        ! Constant-background-pressure CABARET/RSST closures.  The first treats
        ! density as primary; the second treats sensible enthalpy as primary.
        procedure :: apply_state_equation_low_mach_rsst
        procedure :: apply_state_equation_low_mach_rsst_enthalpy

        ! Shared thermodynamic helpers.
        procedure, private :: prepare_composition
        procedure, private :: evaluate_caloric_state
        procedure, private :: assert_temperature_in_table
    end type table_approximated_real_gas

    interface table_approximated_real_gas_c
        module procedure constructor
    end interface table_approximated_real_gas_c

contains

    !> Bind EOS fields to the data manager and register derived quantities.
    type(table_approximated_real_gas) function constructor(manager) result(eos)
        type(data_manager), intent(inout) :: manager

        type(field_scalar_cons_pointer) :: scalar_cell
        type(field_vector_cons_pointer) :: vector_cell
        type(field_tensor_cons_pointer) :: tensor_cell
        type(field_scalar_flow_pointer) :: scalar_face
        type(field_vector_flow_pointer) :: vector_face
        type(field_tensor_flow_pointer) :: tensor_face

        call manager%get_cons_field_pointer_by_name( &
            scalar_cell, vector_cell, tensor_cell, 'temperature')
        eos%T%s_ptr => scalar_cell%s_ptr
        call manager%get_cons_field_pointer_by_name( &
            scalar_cell, vector_cell, tensor_cell, 'pressure')
        eos%p%s_ptr => scalar_cell%s_ptr
        call manager%get_cons_field_pointer_by_name( &
            scalar_cell, vector_cell, tensor_cell, 'pressure_static')
        eos%p_stat%s_ptr => scalar_cell%s_ptr
        call manager%get_cons_field_pointer_by_name( &
            scalar_cell, vector_cell, tensor_cell, 'pressure_static_old')
        eos%p_stat_old%s_ptr => scalar_cell%s_ptr
        call manager%get_cons_field_pointer_by_name( &
            scalar_cell, vector_cell, tensor_cell, 'pressure_static_change')
        eos%dp_stat_dt%s_ptr => scalar_cell%s_ptr
        call manager%get_cons_field_pointer_by_name( &
            scalar_cell, vector_cell, tensor_cell, 'density')
        eos%rho%s_ptr => scalar_cell%s_ptr
        call manager%get_cons_field_pointer_by_name( &
            scalar_cell, vector_cell, tensor_cell, 'internal_energy')
        eos%e_i%s_ptr => scalar_cell%s_ptr
        call manager%get_cons_field_pointer_by_name( &
            scalar_cell, vector_cell, tensor_cell, 'full_energy')
        eos%E_f%s_ptr => scalar_cell%s_ptr
        call manager%get_cons_field_pointer_by_name( &
            scalar_cell, vector_cell, tensor_cell, 'mixture_molar_mass')
        eos%mix_mol_mass%s_ptr => scalar_cell%s_ptr
        call manager%get_cons_field_pointer_by_name( &
            scalar_cell, vector_cell, tensor_cell, 'velocity')
        eos%v%v_ptr => vector_cell%v_ptr
        call manager%get_cons_field_pointer_by_name( &
            scalar_cell, vector_cell, tensor_cell, 'specie_mass_fraction')
        eos%Y%v_ptr => vector_cell%v_ptr

        call manager%get_flow_field_pointer_by_name( &
            scalar_face, vector_face, tensor_face, 'density_flow')
        eos%rho_f%s_ptr => scalar_face%s_ptr
        call manager%get_flow_field_pointer_by_name( &
            scalar_face, vector_face, tensor_face, 'pressure_flow')
        eos%p_f%s_ptr => scalar_face%s_ptr
        call manager%get_flow_field_pointer_by_name( &
            scalar_face, vector_face, tensor_face, 'internal_energy_flow')
        eos%e_i_f%s_ptr => scalar_face%s_ptr
        call manager%get_flow_field_pointer_by_name( &
            scalar_face, vector_face, tensor_face, 'full_energy_flow')
        eos%E_f_f%s_ptr => scalar_face%s_ptr
        call manager%get_flow_field_pointer_by_name( &
            scalar_face, vector_face, tensor_face, 'velocity_of_sound_flow')
        eos%v_s_f%s_ptr => scalar_face%s_ptr
        call manager%get_flow_field_pointer_by_name( &
            scalar_face, vector_face, tensor_face, 'temperature_flow')
        eos%T_f%s_ptr => scalar_face%s_ptr
        call manager%get_flow_field_pointer_by_name( &
            scalar_face, vector_face, tensor_face, 'velocity_flow')
        eos%v_f%v_ptr => vector_face%v_ptr
        call manager%get_flow_field_pointer_by_name( &
            scalar_face, vector_face, tensor_face, 'specie_mass_fraction_flow')
        eos%Y_f%v_ptr => vector_face%v_ptr

        call manager%create_scalar_field( &
            sound_speed_field, 'velocity_of_sound', 'v_s')
        eos%v_s%s_ptr => sound_speed_field
        call manager%create_scalar_field( &
            gamma_field, 'adiabatic_index', 'gamma')
        eos%gamma%s_ptr => gamma_field
        call manager%create_scalar_field( &
            sensible_enthalpy_field, 'sensible_enthalpy', 'h_s')
        eos%h_s%s_ptr => sensible_enthalpy_field
        call manager%create_scalar_field( &
            mixture_cp_field, 'mixture_cp', 'cp')
        eos%mixture_cp%s_ptr => mixture_cp_field
        call manager%create_scalar_field( &
            full_enthalpy_field, 'mixture_enthalpy', 'h_full')
        eos%h_full%s_ptr => full_enthalpy_field
        call manager%create_scalar_field( &
            gamma_face_field, 'adiabatic_index_flow', 'gamma_flow')
        eos%gamma_f%s_ptr => gamma_face_field

        eos%domain = manager%domain
        eos%thermo%thermo_ptr => manager%thermophysics%thermo_ptr
        eos%chem%chem_ptr => manager%chemistry%chem_ptr
        eos%boundary%bc_ptr => manager%boundary_conditions_pointer%bc_ptr
    end function constructor


    !> Normalize mass fractions and construct mole fractions and mixture molar mass.
    !!
    !! Small negative values generated by roundoff are clipped.  A materially
    !! negative component or a zero total composition is rejected because any
    !! silent fallback would change the gas identity and hide a transport error.
    subroutine prepare_composition(this, Y_mass, X_mole, molar_mass, context)
        class(table_approximated_real_gas), intent(in) :: this
        real(dp), intent(inout) :: Y_mass(:)
        real(dp), intent(out) :: X_mole(:)
        real(dp), intent(out) :: molar_mass
        character(len=*), intent(in) :: context

        integer :: spec, species_number
        real(dp) :: mass_sum, molar_denominator, species_molar_mass

        species_number = this%chem%chem_ptr%species_number
        if (size(Y_mass) /= species_number .or. size(X_mole) /= species_number) then
            error stop 'EOS composition helper: inconsistent species-vector size'
        end if

        mass_sum = 0.0_dp
        do spec = 1, species_number
            if (Y_mass(spec) < -eos_composition_tolerance) then
                write(*,'(A,1X,A,1X,I0,1X,ES24.16)') &
                    'EOS negative mass fraction:', trim(context), spec, Y_mass(spec)
                error stop 'EOS: materially negative mass fraction'
            end if
            Y_mass(spec) = max(Y_mass(spec), 0.0_dp)
            mass_sum = mass_sum + Y_mass(spec)
        end do
        if (mass_sum <= tiny(1.0_dp)) then
            write(*,'(A,1X,A)') 'EOS zero species sum:', trim(context)
            error stop 'EOS: zero species sum'
        end if
        Y_mass = Y_mass/mass_sum

        molar_denominator = 0.0_dp
        do spec = 1, species_number
            species_molar_mass = this%thermo%thermo_ptr%molar_masses(spec)
            if (species_molar_mass <= tiny(1.0_dp)) then
                write(*,'(A,1X,A,1X,I0)') &
                    'EOS invalid species molar mass:', trim(context), spec
                error stop 'EOS: non-positive species molar mass'
            end if
            molar_denominator = molar_denominator + &
                Y_mass(spec)/species_molar_mass
        end do
        if (molar_denominator <= tiny(1.0_dp)) then
            error stop 'EOS: invalid mixture molar denominator'
        end if

        molar_mass = 1.0_dp/molar_denominator
        do spec = 1, species_number
            X_mole(spec) = Y_mass(spec)*molar_mass/ &
                this%thermo%thermo_ptr%molar_masses(spec)
        end do
    end subroutine prepare_composition


    !> Evaluate caloric properties at a known temperature and composition.
    !!
    !! mixture_cp is returned in molar units [J/(mol K)] consistently in every
    !! closure.  Consumers that require mass-specific cp must divide by M.
    subroutine evaluate_caloric_state(this, temperature, X_mole, molar_mass, &
        cp_molar, gamma_value, sensible_enthalpy, full_enthalpy, shifted_internal_energy)
        class(table_approximated_real_gas), intent(in) :: this
        real(dp), intent(in) :: temperature, X_mole(:), molar_mass
        real(dp), intent(out) :: cp_molar, gamma_value
        real(dp), intent(out) :: sensible_enthalpy, full_enthalpy
        real(dp), intent(out) :: shifted_internal_energy

        real(dp) :: cv_molar, reference_enthalpy, enthalpy_molar
        real(dp) :: internal_energy_molar

        call this%assert_temperature_in_table(temperature, 'caloric-state evaluation')
        if (molar_mass <= tiny(1.0_dp)) then
            error stop 'EOS caloric state: non-positive mixture molar mass'
        end if

        cp_molar = this%thermo%thermo_ptr%mixture_cp_molar(temperature, X_mole)
        cv_molar = cp_molar - r_gase_J
        if (cp_molar <= 0.0_dp .or. cv_molar <= 0.0_dp) then
            error stop 'EOS caloric state: non-positive heat capacity'
        end if
        gamma_value = cp_molar/cv_molar

        reference_enthalpy = &
            this%thermo%thermo_ptr%mixture_enthalpy_molar(T_ref, X_mole)
        enthalpy_molar = &
            this%thermo%thermo_ptr%mixture_enthalpy_molar(temperature, X_mole)
        internal_energy_molar = &
            this%thermo%thermo_ptr%mixture_internal_energy_molar(temperature, X_mole)

        sensible_enthalpy = (enthalpy_molar-reference_enthalpy)/molar_mass
        full_enthalpy = enthalpy_molar/molar_mass
        shifted_internal_energy = &
            (internal_energy_molar-reference_enthalpy)/molar_mass
    end subroutine evaluate_caloric_state


    !> Reject temperatures outside the thermodynamic table interval.
    subroutine assert_temperature_in_table(this, temperature, context)
        class(table_approximated_real_gas), intent(in) :: this
        real(dp), intent(in) :: temperature
        character(len=*), intent(in) :: context

        if (.not. associated(this%thermo%thermo_ptr)) then
            error stop 'EOS: thermophysical-properties pointer is not associated'
        end if
        if (temperature <= eos_temperature_floor .or. &
            temperature >= tables_temperature_ceiling) then
            write(*,'(A,1X,A,1X,ES24.16)') &
                'EOS temperature outside table:', trim(context), temperature
            error stop 'EOS: temperature outside thermodynamic table'
        end if
    end subroutine assert_temperature_in_table


    !> Construct a complete cell state from prescribed p, T, Y, and velocity.
    !!
    !! This is the startup closure used before any conservative CABARET update.
    !! Density follows the ideal-mixture law and total energy is assembled from
    !! the shifted internal energy plus all resolved kinetic-energy components.
    subroutine apply_state_equation_for_initial_conditions(this)
        class(table_approximated_real_gas), intent(inout) :: this

        integer :: dimensions, species_number
        integer :: i, j, k, dim, spec
        integer :: cell_loop(3,2)
        real(dp) :: molar_mass, cp_molar, gamma_value
        real(dp) :: sensible_enthalpy, full_enthalpy, shifted_internal_energy
        real(dp) :: kinetic_energy
        real(dp), allocatable :: Y_mass(:), X_mole(:)

        dimensions = this%domain%get_domain_dimensions()
        species_number = this%chem%chem_ptr%species_number
        cell_loop = this%domain%get_local_inner_cells_bounds()

        associate(T => this%T%s_ptr, p => this%p%s_ptr, rho => this%rho%s_ptr, &
            e_i => this%e_i%s_ptr, E_f => this%E_f%s_ptr, h_s => this%h_s%s_ptr, &
            h_full => this%h_full%s_ptr, M => this%mix_mol_mass%s_ptr, &
            sound_speed => this%v_s%s_ptr, gamma => this%gamma%s_ptr, &
            cp => this%mixture_cp%s_ptr, Y => this%Y%v_ptr, v => this%v%v_ptr, &
            bc => this%boundary%bc_ptr)

        !$omp parallel default(shared) private(i,j,k,dim,spec,molar_mass,cp_molar, &
        !$omp& gamma_value,sensible_enthalpy,full_enthalpy,shifted_internal_energy, &
        !$omp& kinetic_energy,Y_mass,X_mole)
        allocate(Y_mass(species_number), X_mole(species_number))
        !$omp do collapse(3) schedule(static)
        do k = cell_loop(3,1), cell_loop(3,2)
            do j = cell_loop(2,1), cell_loop(2,2)
                do i = cell_loop(1,1), cell_loop(1,2)
                    if (bc%bc_markers(i,j,k) /= 0) cycle
                    if (p%cells(i,j,k) <= eos_pressure_floor) then
                        error stop 'EOS initial state: non-positive pressure'
                    end if
                    call this%assert_temperature_in_table( &
                        T%cells(i,j,k), 'initial cell state')

                    do spec = 1, species_number
                        Y_mass(spec) = Y%pr(spec)%cells(i,j,k)
                    end do
                    call this%prepare_composition( &
                        Y_mass, X_mole, molar_mass, 'initial cell state')
                    do spec = 1, species_number
                        Y%pr(spec)%cells(i,j,k) = Y_mass(spec)
                    end do

                    rho%cells(i,j,k) = p%cells(i,j,k)*molar_mass/ &
                        (r_gase_J*T%cells(i,j,k))
                    if (rho%cells(i,j,k) <= eos_density_floor) then
                        error stop 'EOS initial state: non-positive density'
                    end if

                    call this%evaluate_caloric_state( &
                        T%cells(i,j,k), X_mole, molar_mass, cp_molar, &
                        gamma_value, sensible_enthalpy, full_enthalpy, &
                        shifted_internal_energy)
                    M%cells(i,j,k) = molar_mass
                    cp%cells(i,j,k) = cp_molar
                    gamma%cells(i,j,k) = gamma_value
                    h_s%cells(i,j,k) = sensible_enthalpy
                    h_full%cells(i,j,k) = full_enthalpy
                    e_i%cells(i,j,k) = shifted_internal_energy

                    kinetic_energy = 0.0_dp
                    do dim = 1, dimensions
                        kinetic_energy = kinetic_energy + &
                            0.5_dp*v%pr(dim)%cells(i,j,k)**2
                    end do
                    E_f%cells(i,j,k) = shifted_internal_energy + kinetic_energy
                    sound_speed%cells(i,j,k) = sqrt( &
                        gamma_value*p%cells(i,j,k)/rho%cells(i,j,k))
                end do
            end do
        end do
        !$omp end do
        deallocate(Y_mass, X_mole)
        !$omp end parallel
        end associate

        call this%apply_boundary_conditions_for_initial_conditions()
    end subroutine apply_state_equation_for_initial_conditions


    !> Recover p, T, and caloric fields from compressible conservative variables.
    subroutine apply_state_equation(this)
        class(table_approximated_real_gas), intent(inout) :: this

        integer :: dimensions, species_number
        integer :: i, j, k, dim, spec
        integer :: cell_loop(3,2)
        real(dp) :: molar_mass, cp_molar, gamma_value
        real(dp) :: sensible_enthalpy, full_enthalpy, shifted_internal_energy
        real(dp) :: kinetic_energy, target_internal_energy_molar, temperature_guess
        real(dp), allocatable :: Y_mass(:), X_mole(:)

        dimensions = this%domain%get_domain_dimensions()
        species_number = this%chem%chem_ptr%species_number
        cell_loop = this%domain%get_local_inner_cells_bounds()

        associate(T => this%T%s_ptr, p => this%p%s_ptr, rho => this%rho%s_ptr, &
            e_i => this%e_i%s_ptr, E_f => this%E_f%s_ptr, h_s => this%h_s%s_ptr, &
            h_full => this%h_full%s_ptr, M => this%mix_mol_mass%s_ptr, &
            sound_speed => this%v_s%s_ptr, gamma => this%gamma%s_ptr, &
            cp => this%mixture_cp%s_ptr, Y => this%Y%v_ptr, v => this%v%v_ptr, &
            bc => this%boundary%bc_ptr)

        !$omp parallel default(shared) private(i,j,k,dim,spec,molar_mass,cp_molar, &
        !$omp& gamma_value,sensible_enthalpy,full_enthalpy,shifted_internal_energy, &
        !$omp& kinetic_energy,target_internal_energy_molar,temperature_guess, &
        !$omp& Y_mass,X_mole)
        allocate(Y_mass(species_number), X_mole(species_number))
        !$omp do collapse(3) schedule(static)
        do k = cell_loop(3,1), cell_loop(3,2)
            do j = cell_loop(2,1), cell_loop(2,2)
                do i = cell_loop(1,1), cell_loop(1,2)
                    if (bc%bc_markers(i,j,k) /= 0) cycle
                    if (rho%cells(i,j,k) <= eos_density_floor) then
                        error stop 'Compressible EOS: non-positive density'
                    end if

                    do spec = 1, species_number
                        Y_mass(spec) = Y%pr(spec)%cells(i,j,k)
                    end do
                    call this%prepare_composition( &
                        Y_mass, X_mole, molar_mass, 'compressible cell state')
                    do spec = 1, species_number
                        Y%pr(spec)%cells(i,j,k) = Y_mass(spec)
                    end do

                    kinetic_energy = 0.0_dp
                    do dim = 1, dimensions
                        kinetic_energy = kinetic_energy + &
                            0.5_dp*v%pr(dim)%cells(i,j,k)**2
                    end do
                    e_i%cells(i,j,k) = E_f%cells(i,j,k)-kinetic_energy
                    target_internal_energy_molar = &
                        e_i%cells(i,j,k)*molar_mass
                    temperature_guess = max( &
                        eos_temperature_floor, min(T%cells(i,j,k), &
                        nearest(tables_temperature_ceiling, -1.0_dp)))
                    T%cells(i,j,k) = this%thermo%thermo_ptr%calculate_temperature( &
                        temperature_guess, target_internal_energy_molar, X_mole)
                    call this%assert_temperature_in_table( &
                        T%cells(i,j,k), 'compressible cell state')

                    p%cells(i,j,k) = rho%cells(i,j,k)*r_gase_J* &
                        T%cells(i,j,k)/molar_mass
                    if (p%cells(i,j,k) <= eos_pressure_floor) then
                        error stop 'Compressible EOS: non-positive pressure'
                    end if

                    call this%evaluate_caloric_state( &
                        T%cells(i,j,k), X_mole, molar_mass, cp_molar, &
                        gamma_value, sensible_enthalpy, full_enthalpy, &
                        shifted_internal_energy)
                    M%cells(i,j,k) = molar_mass
                    cp%cells(i,j,k) = cp_molar
                    gamma%cells(i,j,k) = gamma_value
                    h_s%cells(i,j,k) = sensible_enthalpy
                    h_full%cells(i,j,k) = full_enthalpy
                    e_i%cells(i,j,k) = shifted_internal_energy
                    sound_speed%cells(i,j,k) = sqrt( &
                        gamma_value*p%cells(i,j,k)/rho%cells(i,j,k))
                end do
            end do
        end do
        !$omp end do
        deallocate(Y_mass, X_mole)
        !$omp end parallel
        end associate
    end subroutine apply_state_equation


    !> Update the legacy FDS low-Mach thermodynamic state.
    !!
    !! predictor=.true. advances p0 explicitly and stores the old pressure;
    !! predictor=.false. applies the trapezoidal correction used by the FDS
    !! pressure evolution.  Density is primary, so temperature follows directly
    !! from p0 = rho R T / M.
    subroutine apply_state_equation_low_mach_fds(this, time_step, predictor)
        class(table_approximated_real_gas), intent(inout) :: this
        real(dp), intent(in) :: time_step
        logical, intent(in) :: predictor

        integer :: dimensions, species_number
        integer :: i, j, k, dim, spec
        integer :: cell_loop(3,2)
        real(dp) :: molar_mass, cp_molar, gamma_value
        real(dp) :: sensible_enthalpy, full_enthalpy, shifted_internal_energy
        real(dp) :: kinetic_energy
        real(dp), allocatable :: Y_mass(:), X_mole(:)

        dimensions = this%domain%get_domain_dimensions()
        species_number = this%chem%chem_ptr%species_number
        cell_loop = this%domain%get_local_inner_cells_bounds()

        associate(T => this%T%s_ptr, p => this%p%s_ptr, p0 => this%p_stat%s_ptr, &
            p0_old => this%p_stat_old%s_ptr, dp0dt => this%dp_stat_dt%s_ptr, &
            rho => this%rho%s_ptr, e_i => this%e_i%s_ptr, E_f => this%E_f%s_ptr, &
            h_s => this%h_s%s_ptr, h_full => this%h_full%s_ptr, &
            M => this%mix_mol_mass%s_ptr, sound_speed => this%v_s%s_ptr, &
            gamma => this%gamma%s_ptr, cp => this%mixture_cp%s_ptr, &
            Y => this%Y%v_ptr, v => this%v%v_ptr, bc => this%boundary%bc_ptr)

        !$omp parallel default(shared) private(i,j,k,dim,spec,molar_mass,cp_molar, &
        !$omp& gamma_value,sensible_enthalpy,full_enthalpy,shifted_internal_energy, &
        !$omp& kinetic_energy,Y_mass,X_mole)
        allocate(Y_mass(species_number), X_mole(species_number))
        !$omp do collapse(3) schedule(static)
        do k = cell_loop(3,1), cell_loop(3,2)
            do j = cell_loop(2,1), cell_loop(2,2)
                do i = cell_loop(1,1), cell_loop(1,2)
                    if (bc%bc_markers(i,j,k) /= 0) cycle
                    if (rho%cells(i,j,k) <= eos_density_floor) then
                        error stop 'FDS low-Mach EOS: non-positive density'
                    end if

                    do spec = 1, species_number
                        Y_mass(spec) = Y%pr(spec)%cells(i,j,k)
                    end do
                    call this%prepare_composition( &
                        Y_mass, X_mole, molar_mass, 'FDS low-Mach cell state')
                    do spec = 1, species_number
                        Y%pr(spec)%cells(i,j,k) = Y_mass(spec)
                    end do

                    if (predictor) then
                        p0_old%cells(i,j,k) = p0%cells(i,j,k)
                        p0%cells(i,j,k) = p0%cells(i,j,k) + &
                            time_step*dp0dt%cells(i,j,k)
                    else
                        p0%cells(i,j,k) = 0.5_dp*(p0%cells(i,j,k) + &
                            p0_old%cells(i,j,k) + time_step*dp0dt%cells(i,j,k))
                    end if
                    if (p0%cells(i,j,k) <= eos_pressure_floor) then
                        error stop 'FDS low-Mach EOS: non-positive thermodynamic pressure'
                    end if

                    T%cells(i,j,k) = p0%cells(i,j,k)*molar_mass/ &
                        (rho%cells(i,j,k)*r_gase_J)
                    call this%assert_temperature_in_table( &
                        T%cells(i,j,k), 'FDS low-Mach cell state')
                    call this%evaluate_caloric_state( &
                        T%cells(i,j,k), X_mole, molar_mass, cp_molar, &
                        gamma_value, sensible_enthalpy, full_enthalpy, &
                        shifted_internal_energy)

                    M%cells(i,j,k) = molar_mass
                    cp%cells(i,j,k) = cp_molar
                    gamma%cells(i,j,k) = gamma_value
                    h_s%cells(i,j,k) = sensible_enthalpy
                    h_full%cells(i,j,k) = full_enthalpy
                    e_i%cells(i,j,k) = shifted_internal_energy
                    p%cells(i,j,k) = p0%cells(i,j,k)
                    sound_speed%cells(i,j,k) = sqrt( &
                        gamma_value*p0%cells(i,j,k)/rho%cells(i,j,k))

                    ! The submitted code reset E_f inside this loop and therefore
                    ! retained only the last velocity component in 2-D/3-D.
                    kinetic_energy = 0.0_dp
                    do dim = 1, dimensions
                        kinetic_energy = kinetic_energy + &
                            0.5_dp*v%pr(dim)%cells(i,j,k)**2
                    end do
                    E_f%cells(i,j,k) = shifted_internal_energy + kinetic_energy
                end do
            end do
        end do
        !$omp end do
        deallocate(Y_mass, X_mole)
        !$omp end parallel
        end associate
    end subroutine apply_state_equation_low_mach_fds


    !> Constant-pressure closure with density as the primary material variable.
    !!
    !! This form is used by relaxed/RSST CABARET branches after density/species
    !! transport.  Dynamic pressure is intentionally excluded: only the uniform
    !! thermodynamic pressure closes the ideal-gas relation.
    subroutine apply_state_equation_low_mach_rsst(this, background_pressure)
        class(table_approximated_real_gas), intent(inout) :: this
        real(dp), intent(in) :: background_pressure

        integer :: dimensions, species_number
        integer :: i, j, k, dim, spec
        integer :: cell_loop(3,2)
        real(dp) :: molar_mass, cp_molar, gamma_value
        real(dp) :: sensible_enthalpy, full_enthalpy, shifted_internal_energy
        real(dp) :: kinetic_energy
        real(dp), allocatable :: Y_mass(:), X_mole(:)

        if (background_pressure <= eos_pressure_floor) then
            error stop 'RSST density EOS: non-positive background pressure'
        end if
        dimensions = this%domain%get_domain_dimensions()
        species_number = this%chem%chem_ptr%species_number
        cell_loop = this%domain%get_local_inner_cells_bounds()

        associate(T => this%T%s_ptr, p => this%p%s_ptr, rho => this%rho%s_ptr, &
            e_i => this%e_i%s_ptr, E_f => this%E_f%s_ptr, h_s => this%h_s%s_ptr, &
            h_full => this%h_full%s_ptr, M => this%mix_mol_mass%s_ptr, &
            sound_speed => this%v_s%s_ptr, gamma => this%gamma%s_ptr, &
            cp => this%mixture_cp%s_ptr, Y => this%Y%v_ptr, v => this%v%v_ptr, &
            bc => this%boundary%bc_ptr)

        !$omp parallel default(shared) private(i,j,k,dim,spec,molar_mass,cp_molar, &
        !$omp& gamma_value,sensible_enthalpy,full_enthalpy,shifted_internal_energy, &
        !$omp& kinetic_energy,Y_mass,X_mole)
        allocate(Y_mass(species_number), X_mole(species_number))
        !$omp do collapse(3) schedule(static)
        do k = cell_loop(3,1), cell_loop(3,2)
            do j = cell_loop(2,1), cell_loop(2,2)
                do i = cell_loop(1,1), cell_loop(1,2)
                    if (bc%bc_markers(i,j,k) /= 0) cycle
                    if (rho%cells(i,j,k) <= eos_density_floor) then
                        error stop 'RSST density EOS: non-positive density'
                    end if

                    do spec = 1, species_number
                        Y_mass(spec) = Y%pr(spec)%cells(i,j,k)
                    end do
                    call this%prepare_composition( &
                        Y_mass, X_mole, molar_mass, 'RSST density cell state')
                    do spec = 1, species_number
                        Y%pr(spec)%cells(i,j,k) = Y_mass(spec)
                    end do

                    T%cells(i,j,k) = background_pressure*molar_mass/ &
                        (rho%cells(i,j,k)*r_gase_J)
                    call this%assert_temperature_in_table( &
                        T%cells(i,j,k), 'RSST density cell state')
                    call this%evaluate_caloric_state( &
                        T%cells(i,j,k), X_mole, molar_mass, cp_molar, &
                        gamma_value, sensible_enthalpy, full_enthalpy, &
                        shifted_internal_energy)

                    M%cells(i,j,k) = molar_mass
                    cp%cells(i,j,k) = cp_molar
                    gamma%cells(i,j,k) = gamma_value
                    h_s%cells(i,j,k) = sensible_enthalpy
                    h_full%cells(i,j,k) = full_enthalpy
                    e_i%cells(i,j,k) = shifted_internal_energy
                    p%cells(i,j,k) = background_pressure
                    sound_speed%cells(i,j,k) = sqrt( &
                        gamma_value*background_pressure/rho%cells(i,j,k))

                    kinetic_energy = 0.0_dp
                    do dim = 1, dimensions
                        kinetic_energy = kinetic_energy + &
                            0.5_dp*v%pr(dim)%cells(i,j,k)**2
                    end do
                    E_f%cells(i,j,k) = shifted_internal_energy + kinetic_energy
                end do
            end do
        end do
        !$omp end do
        deallocate(Y_mass, X_mole)
        !$omp end parallel
        end associate
    end subroutine apply_state_equation_low_mach_rsst


    !> Constant-pressure closure with sensible enthalpy as the primary variable.
    !!
    !! The low-Mach CABARET material equation advances h_s and Y.  This routine
    !! inverts h_s(T,Y), computes rho from p0, and refreshes all dependent fields.
    subroutine apply_state_equation_low_mach_rsst_enthalpy(this, background_pressure)
        class(table_approximated_real_gas), intent(inout) :: this
        real(dp), intent(in) :: background_pressure

        integer :: dimensions, species_number
        integer :: i, j, k, dim, spec
        integer :: cell_loop(3,2)
        real(dp) :: molar_mass, cp_molar, gamma_value
        real(dp) :: sensible_enthalpy, full_enthalpy, shifted_internal_energy
        real(dp) :: kinetic_energy, temperature_guess, target_enthalpy_molar
        real(dp), allocatable :: Y_mass(:), X_mole(:)

        if (background_pressure <= eos_pressure_floor) then
            error stop 'RSST enthalpy EOS: non-positive background pressure'
        end if
        dimensions = this%domain%get_domain_dimensions()
        species_number = this%chem%chem_ptr%species_number
        cell_loop = this%domain%get_local_inner_cells_bounds()

        associate(T => this%T%s_ptr, p => this%p%s_ptr, rho => this%rho%s_ptr, &
            e_i => this%e_i%s_ptr, E_f => this%E_f%s_ptr, h_s => this%h_s%s_ptr, &
            h_full => this%h_full%s_ptr, M => this%mix_mol_mass%s_ptr, &
            sound_speed => this%v_s%s_ptr, gamma => this%gamma%s_ptr, &
            cp => this%mixture_cp%s_ptr, Y => this%Y%v_ptr, v => this%v%v_ptr, &
            bc => this%boundary%bc_ptr)

        !$omp parallel default(shared) private(i,j,k,dim,spec,molar_mass,cp_molar, &
        !$omp& gamma_value,sensible_enthalpy,full_enthalpy,shifted_internal_energy, &
        !$omp& kinetic_energy,temperature_guess,target_enthalpy_molar,Y_mass,X_mole)
        allocate(Y_mass(species_number), X_mole(species_number))
        !$omp do collapse(3) schedule(static)
        do k = cell_loop(3,1), cell_loop(3,2)
            do j = cell_loop(2,1), cell_loop(2,2)
                do i = cell_loop(1,1), cell_loop(1,2)
                    if (bc%bc_markers(i,j,k) /= 0) cycle

                    do spec = 1, species_number
                        Y_mass(spec) = Y%pr(spec)%cells(i,j,k)
                    end do
                    call this%prepare_composition( &
                        Y_mass, X_mole, molar_mass, 'RSST enthalpy cell state')
                    do spec = 1, species_number
                        Y%pr(spec)%cells(i,j,k) = Y_mass(spec)
                    end do

                    temperature_guess = max( &
                        eos_temperature_floor, min(T%cells(i,j,k), &
                        nearest(tables_temperature_ceiling, -1.0_dp)))
                    target_enthalpy_molar = h_s%cells(i,j,k)*molar_mass
                    T%cells(i,j,k) = &
                        this%thermo%thermo_ptr%calculate_temperature_Pconst( &
                        temperature_guess, target_enthalpy_molar, X_mole)
                    call this%assert_temperature_in_table( &
                        T%cells(i,j,k), 'RSST enthalpy cell state')

                    rho%cells(i,j,k) = background_pressure*molar_mass/ &
                        (r_gase_J*T%cells(i,j,k))
                    if (rho%cells(i,j,k) <= eos_density_floor) then
                        error stop 'RSST enthalpy EOS: non-positive density'
                    end if
                    call this%evaluate_caloric_state( &
                        T%cells(i,j,k), X_mole, molar_mass, cp_molar, &
                        gamma_value, sensible_enthalpy, full_enthalpy, &
                        shifted_internal_energy)

                    M%cells(i,j,k) = molar_mass
                    cp%cells(i,j,k) = cp_molar
                    gamma%cells(i,j,k) = gamma_value
                    h_s%cells(i,j,k) = sensible_enthalpy
                    h_full%cells(i,j,k) = full_enthalpy
                    e_i%cells(i,j,k) = shifted_internal_energy
                    p%cells(i,j,k) = background_pressure
                    sound_speed%cells(i,j,k) = sqrt( &
                        gamma_value*background_pressure/rho%cells(i,j,k))

                    kinetic_energy = 0.0_dp
                    do dim = 1, dimensions
                        kinetic_energy = kinetic_energy + &
                            0.5_dp*v%pr(dim)%cells(i,j,k)**2
                    end do
                    E_f%cells(i,j,k) = shifted_internal_energy + kinetic_energy
                end do
            end do
        end do
        !$omp end do
        deallocate(Y_mass, X_mole)
        !$omp end parallel
        end associate
    end subroutine apply_state_equation_low_mach_rsst_enthalpy


    !> Close all active CABARET face states from p_f, rho_f, Y_f, and v_f.
    subroutine apply_state_equation_flow_variables(this)
        class(table_approximated_real_gas), intent(inout) :: this

        integer :: dimensions, species_number
        integer :: i, j, k, dim, component, spec
        integer :: face_bounds(3,2), loop(3,2)
        integer :: il, jl, kl
        real(dp) :: molar_mass, cp_molar, gamma_value
        real(dp) :: sensible_enthalpy, full_enthalpy, shifted_internal_energy
        real(dp) :: kinetic_energy
        real(dp), allocatable :: Y_mass(:), X_mole(:)

        dimensions = this%domain%get_domain_dimensions()
        species_number = this%chem%chem_ptr%species_number
        face_bounds = this%domain%get_local_inner_faces_bounds()

        associate(p_f => this%p_f%s_ptr, rho_f => this%rho_f%s_ptr, &
            e_i_f => this%e_i_f%s_ptr, E_f_f => this%E_f_f%s_ptr, &
            T_f => this%T_f%s_ptr, sound_speed_f => this%v_s_f%s_ptr, &
            gamma_f => this%gamma_f%s_ptr, Y_f => this%Y_f%v_ptr, &
            v_f => this%v_f%v_ptr, bc => this%boundary%bc_ptr)

        !$omp parallel default(shared) private(i,j,k,dim,component,spec,loop, &
        !$omp& il,jl,kl,molar_mass,cp_molar,gamma_value,sensible_enthalpy, &
        !$omp& full_enthalpy,shifted_internal_energy,kinetic_energy,Y_mass,X_mole)
        allocate(Y_mass(species_number), X_mole(species_number))
        do dim = 1, dimensions
            loop = face_bounds
            do component = 1, dimensions
                loop(component,2) = face_bounds(component,2) - &
                    (1-I_m(component,dim))
            end do

            !$omp do collapse(3) schedule(static)
            do k = loop(3,1), loop(3,2)
                do j = loop(2,1), loop(2,2)
                    do i = loop(1,1), loop(1,2)
                        il = i-I_m(dim,1)
                        jl = j-I_m(dim,2)
                        kl = k-I_m(dim,3)
                        if (bc%bc_markers(i,j,k) /= 0 .and. &
                            bc%bc_markers(il,jl,kl) /= 0) cycle
                        if (p_f%cells(dim,i,j,k) <= eos_pressure_floor) then
                            error stop 'Face EOS: non-positive pressure'
                        end if
                        if (rho_f%cells(dim,i,j,k) <= eos_density_floor) then
                            error stop 'Face EOS: non-positive density'
                        end if

                        do spec = 1, species_number
                            Y_mass(spec) = Y_f%pr(spec)%cells(dim,i,j,k)
                        end do
                        call this%prepare_composition( &
                            Y_mass, X_mole, molar_mass, 'CABARET face state')
                        do spec = 1, species_number
                            Y_f%pr(spec)%cells(dim,i,j,k) = Y_mass(spec)
                        end do

                        T_f%cells(dim,i,j,k) = p_f%cells(dim,i,j,k)*molar_mass/ &
                            (rho_f%cells(dim,i,j,k)*r_gase_J)
                        call this%assert_temperature_in_table( &
                            T_f%cells(dim,i,j,k), 'CABARET face state')
                        call this%evaluate_caloric_state( &
                            T_f%cells(dim,i,j,k), X_mole, molar_mass, cp_molar, &
                            gamma_value, sensible_enthalpy, full_enthalpy, &
                            shifted_internal_energy)

                        gamma_f%cells(dim,i,j,k) = gamma_value
                        e_i_f%cells(dim,i,j,k) = shifted_internal_energy
                        sound_speed_f%cells(dim,i,j,k) = sqrt( &
                            gamma_value*p_f%cells(dim,i,j,k)/ &
                            rho_f%cells(dim,i,j,k))

                        kinetic_energy = 0.0_dp
                        do component = 1, dimensions
                            kinetic_energy = kinetic_energy + &
                                0.5_dp*v_f%pr(component)%cells(dim,i,j,k)**2
                        end do
                        E_f_f%cells(dim,i,j,k) = &
                            shifted_internal_energy + kinetic_energy
                    end do
                end do
            end do
            !$omp end do nowait
        end do
        deallocate(Y_mass, X_mole)
        !$omp end parallel
        end associate
    end subroutine apply_state_equation_flow_variables


    !> Initial face closure retained for API compatibility.
    subroutine apply_state_equation_flow_variables_for_IC(this)
        class(table_approximated_real_gas), intent(inout) :: this
        call this%apply_state_equation_flow_variables()
    end subroutine apply_state_equation_flow_variables_for_IC


    !> Construct thermodynamically consistent physical ghost cells at startup.
    !!
    !! Hydrodynamic solvers reapply their own boundary conditions every step;
    !! this routine exists to ensure the first EOS/face reconstruction never sees
    !! stale density, energy, heat capacity, or sound speed in a boundary ghost.
    subroutine apply_boundary_conditions_for_initial_conditions(this)
        class(table_approximated_real_gas), intent(inout) :: this

        integer :: dimensions, species_number
        integer :: i, j, k, dim, side, sign, component, spec
        integer :: ghost_i, ghost_j, ghost_k, boundary_number, species_index
        integer :: cell_loop(3,2)
        character(len=20) :: boundary_name
        real(dp) :: pressure_ghost, temperature_ghost, molar_mass
        real(dp) :: cp_molar, gamma_value, sensible_enthalpy, full_enthalpy
        real(dp) :: shifted_internal_energy, kinetic_energy, farfield_velocity
        real(dp) :: farfield_pressure, farfield_temperature, farfield_density
        real(dp), allocatable :: Y_mass(:), X_mole(:)
        real(dp), allocatable :: farfield_composition(:)
        character(len=10), allocatable :: farfield_species_names(:)

        dimensions = this%domain%get_domain_dimensions()
        species_number = this%chem%chem_ptr%species_number
        cell_loop = this%domain%get_local_inner_cells_bounds()
        allocate(Y_mass(species_number), X_mole(species_number))

        associate(T => this%T%s_ptr, p => this%p%s_ptr, rho => this%rho%s_ptr, &
            e_i => this%e_i%s_ptr, E_f => this%E_f%s_ptr, h_s => this%h_s%s_ptr, &
            h_full => this%h_full%s_ptr, M => this%mix_mol_mass%s_ptr, &
            sound_speed => this%v_s%s_ptr, gamma => this%gamma%s_ptr, &
            cp => this%mixture_cp%s_ptr, Y => this%Y%v_ptr, v => this%v%v_ptr, &
            bc => this%boundary%bc_ptr)

        do k = cell_loop(3,1), cell_loop(3,2)
            do j = cell_loop(2,1), cell_loop(2,2)
                do i = cell_loop(1,1), cell_loop(1,2)
                    if (bc%bc_markers(i,j,k) /= 0) cycle
                    do dim = 1, dimensions
                        do side = 1, 2
                            sign = merge(-1, 1, side == 1)
                            ghost_i = i + sign*I_m(dim,1)
                            ghost_j = j + sign*I_m(dim,2)
                            ghost_k = k + sign*I_m(dim,3)
                            boundary_number = &
                                bc%bc_markers(ghost_i,ghost_j,ghost_k)
                            if (boundary_number == 0) cycle

                            boundary_name = &
                                bc%boundary_types(boundary_number)%get_type_name()
                            select case (boundary_name)
                            case ('wall')
                                pressure_ghost = p%cells(i,j,k)
                                if (bc%boundary_types(boundary_number)%is_conductive()) then
                                    temperature_ghost = bc%boundary_types( &
                                        boundary_number)%get_wall_temperature()
                                else
                                    temperature_ghost = T%cells(i,j,k)
                                end if
                                do spec = 1, species_number
                                    Y_mass(spec) = Y%pr(spec)%cells(i,j,k)
                                end do
                                do component = 1, dimensions
                                    if (component == dim) then
                                        v%pr(component)%cells(ghost_i,ghost_j,ghost_k) = &
                                            -v%pr(component)%cells(i,j,k)
                                    else if (bc%boundary_types( &
                                        boundary_number)%is_slip()) then
                                        v%pr(component)%cells(ghost_i,ghost_j,ghost_k) = &
                                            v%pr(component)%cells(i,j,k)
                                    else
                                        v%pr(component)%cells(ghost_i,ghost_j,ghost_k) = &
                                            -v%pr(component)%cells(i,j,k)
                                    end if
                                end do

                            case ('symmetry_plane')
                                pressure_ghost = p%cells(i,j,k)
                                temperature_ghost = T%cells(i,j,k)
                                do spec = 1, species_number
                                    Y_mass(spec) = Y%pr(spec)%cells(i,j,k)
                                end do
                                do component = 1, dimensions
                                    if (component == dim) then
                                        v%pr(component)%cells(ghost_i,ghost_j,ghost_k) = &
                                            -v%pr(component)%cells(i,j,k)
                                    else
                                        v%pr(component)%cells(ghost_i,ghost_j,ghost_k) = &
                                            v%pr(component)%cells(i,j,k)
                                    end if
                                end do

                            case ('inlet', 'outlet')
                                farfield_pressure = bc%boundary_types( &
                                    boundary_number)%get_farfield_pressure()
                                farfield_temperature = bc%boundary_types( &
                                    boundary_number)%get_farfield_temperature()
                                farfield_velocity = bc%boundary_types( &
                                    boundary_number)%get_farfield_velocity()
                                if (farfield_pressure <= eos_pressure_floor .or. &
                                    farfield_temperature <= eos_temperature_floor) then
                                    error stop 'EOS boundary initialization: invalid far field'
                                end if
                                pressure_ghost = farfield_pressure
                                temperature_ghost = farfield_temperature

                                call bc%boundary_types(boundary_number)% &
                                    get_farfield_species_names(farfield_species_names)
                                call bc%boundary_types(boundary_number)% &
                                    get_farfield_concentrations(farfield_composition)
                                Y_mass = 0.0_dp
                                if (size(farfield_species_names) /= &
                                    size(farfield_composition)) then
                                    error stop 'EOS boundary: inconsistent far-field composition'
                                end if
                                do spec = 1, size(farfield_species_names)
                                    species_index = this%chem%chem_ptr% &
                                        get_chemical_specie_index( &
                                        farfield_species_names(spec))
                                    Y_mass(species_index) = farfield_composition(spec)
                                end do
                                ! Boundary input is a mole-fraction vector.  Convert
                                ! explicitly to mass fractions before the common helper.
                                call this%thermo%thermo_ptr% &
                                    change_cell_units_mole_to_dimless(Y_mass)

                                v%pr(dim)%cells(ghost_i,ghost_j,ghost_k) = &
                                    farfield_velocity
                                do component = 1, dimensions
                                    if (component /= dim) then
                                        v%pr(component)%cells(ghost_i,ghost_j,ghost_k) = &
                                            v%pr(component)%cells(i,j,k)
                                    end if
                                end do

                            case default
                                write(*,'(A,1X,A)') &
                                    'EOS unsupported boundary type:', trim(boundary_name)
                                error stop 'EOS boundary initialization failed'
                            end select

                            call this%prepare_composition( &
                                Y_mass, X_mole, molar_mass, &
                                'initial boundary ghost state')
                            call this%assert_temperature_in_table( &
                                temperature_ghost, 'initial boundary ghost state')
                            if (pressure_ghost <= eos_pressure_floor) then
                                error stop 'EOS boundary: non-positive pressure'
                            end if
                            call this%evaluate_caloric_state( &
                                temperature_ghost, X_mole, molar_mass, cp_molar, &
                                gamma_value, sensible_enthalpy, full_enthalpy, &
                                shifted_internal_energy)

                            p%cells(ghost_i,ghost_j,ghost_k) = pressure_ghost
                            T%cells(ghost_i,ghost_j,ghost_k) = temperature_ghost
                            rho%cells(ghost_i,ghost_j,ghost_k) = &
                                pressure_ghost*molar_mass/ &
                                (r_gase_J*temperature_ghost)
                            M%cells(ghost_i,ghost_j,ghost_k) = molar_mass
                            cp%cells(ghost_i,ghost_j,ghost_k) = cp_molar
                            gamma%cells(ghost_i,ghost_j,ghost_k) = gamma_value
                            h_s%cells(ghost_i,ghost_j,ghost_k) = &
                                sensible_enthalpy
                            h_full%cells(ghost_i,ghost_j,ghost_k) = full_enthalpy
                            e_i%cells(ghost_i,ghost_j,ghost_k) = &
                                shifted_internal_energy
                            do spec = 1, species_number
                                Y%pr(spec)%cells(ghost_i,ghost_j,ghost_k) = &
                                    Y_mass(spec)
                            end do

                            kinetic_energy = 0.0_dp
                            do component = 1, dimensions
                                kinetic_energy = kinetic_energy + 0.5_dp* &
                                    v%pr(component)%cells(ghost_i,ghost_j,ghost_k)**2
                            end do
                            E_f%cells(ghost_i,ghost_j,ghost_k) = &
                                shifted_internal_energy + kinetic_energy
                            sound_speed%cells(ghost_i,ghost_j,ghost_k) = sqrt( &
                                gamma_value*pressure_ghost/ &
                                rho%cells(ghost_i,ghost_j,ghost_k))

                            if (boundary_name == 'inlet' .or. &
                                boundary_name == 'outlet') then
                                farfield_density = rho%cells(ghost_i,ghost_j,ghost_k)
                                call bc%boundary_types(boundary_number)% &
                                    set_farfield_density(farfield_density)
                                call bc%boundary_types(boundary_number)% &
                                    set_farfield_energy( &
                                    E_f%cells(ghost_i,ghost_j,ghost_k))
                            end if
                        end do
                    end do
                end do
            end do
        end do
        end associate

        if (allocated(farfield_species_names)) deallocate(farfield_species_names)
        if (allocated(farfield_composition)) deallocate(farfield_composition)
        deallocate(Y_mass, X_mole)
    end subroutine apply_boundary_conditions_for_initial_conditions

end module table_approximated_real_gas_class
