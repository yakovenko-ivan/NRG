module thermophysical_properties_class

    use kind_parameters
    use global_data

    use chemical_properties_class

    use field_pointers

    implicit none

    private
    public thermophysical_properties, thermophysical_properties_pointer, &
           thermophysical_properties_c

    type :: thermophysical_properties_pointer
        type(thermophysical_properties), pointer :: thermo_ptr
    end type

    type :: thermophysical_properties
        character(len=30) :: thermo_data_file_name
        character(len=30) :: transport_data_file_name

        ! Kept only for backward compatibility with existing package-interface
        ! calls and previously generated /thermo_properties/ namelists.  The
        ! file is no longer opened or used: species molar masses are derived
        ! from the elemental composition in the NASA/CHEMKIN thermo headers.
        character(len=30) :: molar_masses_data_file_name

        real(dp), dimension(:), allocatable :: potential_well_depth
        real(dp), dimension(:), allocatable :: collision_diameter
        real(dp), dimension(:), allocatable :: molar_masses
        real(dp), dimension(:,:), allocatable :: binary_diffusivity_constant
        real(dp), dimension(:,:,:), allocatable :: a_coeffs
        real(dp), dimension(2,2,8) :: omega_c
    contains
        procedure, private :: read_properties
        procedure, private :: write_properties
        procedure, private :: set_properties

        procedure :: specie_cp_molar
        procedure :: specie_entropy_molar
        procedure :: specie_enthalpy_molar
        procedure :: specie_internal_energy_molar

        procedure :: mixture_cp_molar
        procedure :: mixture_enthalpy_molar
        procedure :: mixture_internal_energy_molar

        procedure :: calculate_temperature
        procedure :: calculate_temperature_Pconst
        procedure :: calculate_omega
        procedure :: specie_dynamic_viscosity
        procedure :: mixture_dynamic_viscosity
        procedure :: specie_thermal_conductivity
        procedure :: mixture_thermal_conductivity
        procedure :: binary_diffusion_coefficient
        procedure :: mixture_averaged_diffusion_coefficients
        procedure :: reduced_thermal_diffusion_coefficients
        procedure :: clip_temperature

        ! CABARET/thermally-perfect ideal-gas mixture helpers.
        procedure :: mixture_molar_mass_from_mass_fractions
        procedure :: mole_fractions_from_mass_fractions
        procedure :: temperature_from_pressure_density_Y
        procedure :: mixture_specific_entropy
        procedure :: temperature_from_entropy_pressure_Y
        procedure :: pressure_derivatives_density_energy_Y
        procedure :: pressure_derivatives_density_energy_composition_Y
        procedure :: pressure_source_rate_from_conservative_sources
        procedure :: specific_internal_energy_from_temperature_Y
        procedure :: sensible_internal_energy_molar_from_temperature_X
        procedure :: sensible_energy_parameter_from_temperature_Y
        procedure :: temperature_from_sensible_energy_parameter_Y

        procedure :: change_field_units_mole_to_dimless
        procedure :: change_cell_units_mole_to_dimless

        procedure :: write_log
    end type

    interface thermophysical_properties_c
        module procedure constructor
        module procedure constructor_file
    end interface

contains

    type(thermophysical_properties) function constructor(chemistry, &
            thermo_data_file_name, transport_data_file_name, &
            molar_masses_data_file_name)

        type(chemical_properties), intent(in) :: chemistry
        character(len=*), intent(in) :: thermo_data_file_name
        character(len=*), intent(in) :: transport_data_file_name
        character(len=*), intent(in) :: molar_masses_data_file_name

        integer :: io_unit

        call constructor%set_properties(chemistry, thermo_data_file_name, &
            transport_data_file_name, molar_masses_data_file_name)

        open(newunit=io_unit, file=thermophysical_data_file_name, &
            status='replace', form='formatted', delim='quote')
        call constructor%write_properties(io_unit)
        close(io_unit)
    end function constructor


    type(thermophysical_properties) function constructor_file(chemistry)
        type(chemical_properties), intent(in) :: chemistry

        integer :: io_unit

        open(newunit=io_unit, file=thermophysical_data_file_name, &
            status='old', form='formatted')
        call constructor_file%read_properties(chemistry,io_unit)
        close(io_unit)
    end function constructor_file


    subroutine write_properties(this,thermo_data_unit)
        class(thermophysical_properties), intent(in) :: this
        integer, intent(in) :: thermo_data_unit

        character(len=30) :: thermo_data_file_name
        character(len=30) :: transport_data_file_name
        character(len=30) :: molar_masses_data_file_name

        ! molar_masses_data_file_name remains in the namelist for compatibility
        ! with existing task_setup files, but is no longer used by this class.
        namelist /thermo_properties/ thermo_data_file_name, &
            transport_data_file_name, molar_masses_data_file_name

        thermo_data_file_name        = this%thermo_data_file_name
        transport_data_file_name     = this%transport_data_file_name
        molar_masses_data_file_name  = this%molar_masses_data_file_name

        write(unit=thermo_data_unit,nml=thermo_properties)
    end subroutine write_properties


    subroutine read_properties(this,chemistry,thermo_data_unit)
        class(thermophysical_properties), intent(inout) :: this
        type(chemical_properties), intent(in) :: chemistry
        integer, intent(in) :: thermo_data_unit

        character(len=30) :: thermo_data_file_name
        character(len=30) :: transport_data_file_name
        character(len=30) :: molar_masses_data_file_name

        namelist /thermo_properties/ thermo_data_file_name, &
            transport_data_file_name, molar_masses_data_file_name

        ! Preserve compatibility with old setup files and also make a future
        ! setup without this obsolete item safe if namelist generation changes.
        thermo_data_file_name       = ''
        transport_data_file_name    = ''
        molar_masses_data_file_name = ''

        read(unit=thermo_data_unit,nml=thermo_properties)
        call this%set_properties(chemistry, thermo_data_file_name, &
            transport_data_file_name, molar_masses_data_file_name)
    end subroutine read_properties


    subroutine set_properties(this,chemistry,thermo_data_file_name, &
            transport_data_file_name,molar_masses_data_file_name)

        class(thermophysical_properties), intent(inout) :: this
        type(chemical_properties), intent(in) :: chemistry
        character(len=*), intent(in) :: thermo_data_file_name
        character(len=*), intent(in) :: transport_data_file_name
        character(len=*), intent(in) :: molar_masses_data_file_name

        integer :: species_number
        integer :: thermo_data_file_unit, transport_data_file_unit
        integer :: spec_num
        integer :: error

        species_number = chemistry%species_number

        this%thermo_data_file_name       = thermo_data_file_name
        this%transport_data_file_name    = transport_data_file_name
        this%molar_masses_data_file_name = molar_masses_data_file_name

        allocate(this%molar_masses(species_number))
        allocate(this%potential_well_depth(species_number), &
                 this%collision_diameter(species_number))
        allocate(this%binary_diffusivity_constant(species_number,species_number))
        allocate(this%a_coeffs(7,2,species_number))
        this%binary_diffusivity_constant = 0.0_dp

        open(newunit=thermo_data_file_unit, &
            file=trim(task_setup_folder)//trim(fold_sep)// &
                 trim(thermophysical_data_folder)//trim(fold_sep)// &
                 this%thermo_data_file_name, &
            status='old',iostat=error)
        if (error /= 0) then
            print *, 'Error: unable to open thermophysical data file: ', &
                trim(this%thermo_data_file_name)
            error stop 'thermophysical_properties: thermo file open failure'
        end if

        open(newunit=transport_data_file_unit, &
            file=trim(task_setup_folder)//trim(fold_sep)// &
                 trim(thermophysical_data_folder)//trim(fold_sep)// &
                 this%transport_data_file_name, &
            status='old',iostat=error)
        if (error /= 0) then
            print *, 'Error: unable to open transport data file: ', &
                trim(this%transport_data_file_name)
            error stop 'thermophysical_properties: transport file open failure'
        end if

        ! The former molar_masses.dat OPEN is intentionally absent here.
        ! M_k is now reconstructed from the fixed-width elemental-composition
        ! fields of the same NASA/CHEMKIN thermo record used for a_coeffs.
        do spec_num=1,species_number
            this%a_coeffs(:,:,spec_num) = get_a_coeffs( &
                chemistry%species_names(spec_num),thermo_data_file_unit)
            this%a_coeffs(:,:,spec_num) = &
                this%a_coeffs(:,:,spec_num)*r_gase_J

            this%molar_masses(spec_num) = &
                get_specie_molar_mass_from_thermo( &
                    chemistry%species_names(spec_num),thermo_data_file_unit)
        end do

        do spec_num=1,chemistry%species_number
            call get_transport_properties( &
                this%potential_well_depth(spec_num), &
                this%collision_diameter(spec_num), &
                chemistry%species_names(spec_num),transport_data_file_unit)
        end do

        call initialize_binary_diffusivity_constants(this)

        close(thermo_data_file_unit)
        close(transport_data_file_unit)

        !# P.D. Neufeld, A.R. Janzen, A.R. Aziz // J. Chem. Phys. 57, 1100 (1972)
        this%omega_c(1,1,1) = 1.06036_dp
        this%omega_c(1,1,2) = 0.15610_dp
        this%omega_c(1,1,3) = 0.19300_dp
        this%omega_c(1,1,4) = 0.47635_dp
        this%omega_c(1,1,5) = 1.03587_dp
        this%omega_c(1,1,6) = 1.52996_dp
        this%omega_c(1,1,7) = 1.76474_dp
        this%omega_c(1,1,8) = 3.89411_dp

        this%omega_c(1,2,1) = 1.00220_dp
        this%omega_c(1,2,2) = 0.15530_dp
        this%omega_c(1,2,3) = 0.16105_dp
        this%omega_c(1,2,4) = 0.72751_dp
        this%omega_c(1,2,5) = 0.86125_dp
        this%omega_c(1,2,6) = 2.06848_dp
        this%omega_c(1,2,7) = 1.95162_dp
        this%omega_c(1,2,8) = 4.84492_dp

        this%omega_c(2,2,1) = 1.16145_dp
        this%omega_c(2,2,2) = 0.14874_dp
        this%omega_c(2,2,3) = 0.52487_dp
        this%omega_c(2,2,4) = 0.77320_dp
        this%omega_c(2,2,5) = 2.16178_dp
        this%omega_c(2,2,6) = 2.43787_dp
        this%omega_c(2,2,7) = 0.0_dp
        this%omega_c(2,2,8) = 0.0_dp
    end subroutine set_properties


    subroutine write_log(this,log_unit)
        class(thermophysical_properties), intent(in) :: this
        integer, intent(in) :: log_unit

        write(log_unit,'(A)') &
            '************************************************************************************* '
        write(log_unit,'(A)') ' Thermophtsical data setup : '
        write(log_unit,'(A,A)') ' Thermo data file        : ', &
            this%thermo_data_file_name
        write(log_unit,'(A,A)') ' Transport data file     : ', &
            this%transport_data_file_name
        write(log_unit,'(A)') &
            ' Molar masses source     : NASA/CHEMKIN thermo elemental composition'
        write(log_unit,'(A)') &
            '************************************************************************************* '
    end subroutine write_log


    pure real(dp) function specie_cp_molar(this,temperature,specie_number) result(cp_molar)
        class(thermophysical_properties), intent(in) :: this
        real(dp), intent(in) :: temperature
        integer, intent(in) :: specie_number
        integer :: i, coeff_set
        real(dp) :: temp

        coeff_set = 1
        temp = this%clip_temperature(temperature)
        if (temp >= 1000.0_dp) coeff_set = 2

        cp_molar = this%a_coeffs(5,coeff_set,specie_number)
        do i=4,1,-1
            cp_molar = cp_molar*temp + this%a_coeffs(i,coeff_set,specie_number)
        end do
    end function specie_cp_molar


    pure real(dp) function specie_enthalpy_molar(this,temperature,specie_number) result(hs_molar)
        class(thermophysical_properties), intent(in) :: this
        real(dp), intent(in) :: temperature
        integer, intent(in) :: specie_number
        integer :: i, coeff_set
        real(dp) :: temp

        coeff_set = 1
        temp = this%clip_temperature(temperature)
        if (temp >= 1000.0_dp) coeff_set = 2

        hs_molar = this%a_coeffs(5,coeff_set,specie_number)/5.0_dp
        do i=4,1,-1
            hs_molar = hs_molar*temp + &
                this%a_coeffs(i,coeff_set,specie_number)/real(i,dp)
        end do
        hs_molar = hs_molar + this%a_coeffs(6,coeff_set,specie_number)/temp
        hs_molar = hs_molar*temp
    end function specie_enthalpy_molar


    pure real(dp) function specie_entropy_molar(this,temperature,specie_number) result(s_molar)
        class(thermophysical_properties), intent(in) :: this
        real(dp), intent(in) :: temperature
        integer, intent(in) :: specie_number
        integer :: i, coeff_set
        real(dp) :: temp

        coeff_set = 1
        temp = this%clip_temperature(temperature)
        if (temp >= 1000.0_dp) coeff_set = 2

        s_molar = this%a_coeffs(5,coeff_set,specie_number)/4.0_dp
        do i=4,2,-1
            s_molar = s_molar*temp + &
                this%a_coeffs(i,coeff_set,specie_number)/real(i-1,dp)
        end do
        s_molar = s_molar*temp + this%a_coeffs(7,coeff_set,specie_number) + &
            this%a_coeffs(1,coeff_set,specie_number)*log(temp)
    end function specie_entropy_molar


    pure real(dp) function specie_internal_energy_molar(this,temperature,specie_number) result(e_molar)
        class(thermophysical_properties), intent(in) :: this
        real(dp), intent(in) :: temperature
        integer, intent(in) :: specie_number
        real(dp) :: hs_molar

        hs_molar = this%specie_enthalpy_molar(temperature,specie_number)
        e_molar = hs_molar-r_gase_J*temperature
    end function specie_internal_energy_molar


    pure real(dp) function mixture_cp_molar(this,temperature,species_molar_fractions) result(cp_molar)
        class(thermophysical_properties), intent(in) :: this
        real(dp), intent(in) :: temperature
        real(dp), dimension(:), intent(in) :: species_molar_fractions
        real(dp) :: specie_cp
        integer :: specie_number

        cp_molar = 0.0_dp
        do specie_number=1,size(species_molar_fractions)
            specie_cp = this%specie_cp_molar(temperature,specie_number)
            cp_molar = cp_molar + specie_cp*species_molar_fractions(specie_number)
        end do
    end function mixture_cp_molar


    pure real(dp) function mixture_enthalpy_molar(this,temperature,species_molar_fractions) result(hs_molar)
        class(thermophysical_properties), intent(in) :: this
        real(dp), intent(in) :: temperature
        real(dp), dimension(:), intent(in) :: species_molar_fractions
        real(dp) :: specie_enthalpy
        integer :: specie_number

        hs_molar = 0.0_dp
        do specie_number=1,size(species_molar_fractions)
            specie_enthalpy = this%specie_enthalpy_molar(temperature,specie_number)
            hs_molar = hs_molar + specie_enthalpy*species_molar_fractions(specie_number)
        end do
    end function mixture_enthalpy_molar


    pure real(dp) function mixture_internal_energy_molar(this,temperature,species_molar_fractions) result(e_molar)
        class(thermophysical_properties), intent(in) :: this
        real(dp), intent(in) :: temperature
        real(dp), dimension(:), intent(in) :: species_molar_fractions
        real(dp) :: specie_internal_energy
        integer :: specie_number

        e_molar = 0.0_dp
        do specie_number=1,size(species_molar_fractions)
            specie_internal_energy = &
                this%specie_internal_energy_molar(temperature,specie_number)
            e_molar = e_molar + &
                specie_internal_energy*species_molar_fractions(specie_number)
        end do
    end function mixture_internal_energy_molar


    function calculate_temperature(this,temperature,e_i,molar_fractions)
        real(dp) :: calculate_temperature
        class(thermophysical_properties), intent(in) :: this
        real(dp), intent(in) :: temperature,e_i
        real(dp), dimension(:), intent(in) :: molar_fractions

        real(dp) :: temp
        integer :: iter,max_iter
        real(dp) :: eps,residual,cp,cv,e_int,h_const,e_const

        temp = this%clip_temperature(temperature)
        max_iter = 1000
        eps = 1e-08

        h_const = this%mixture_enthalpy_molar(T_ref,molar_fractions)
        e_const = this%mixture_internal_energy_molar(T_ref,molar_fractions)

        cp = this%mixture_cp_molar(temp,molar_fractions)
        cv = cp-r_gase_J

        do iter=1,max_iter
            e_int = this%mixture_internal_energy_molar(temp,molar_fractions)-h_const
            residual = (e_i-e_int)/cp
            temp = temp+residual
            if (abs(residual)<eps) exit
            if (e_int>0.0_dp) continue
        end do

        calculate_temperature = temp
    end function calculate_temperature


    recursive pure function calculate_temperature_Pconst(this,temperature,h_s,molar_fractions)
        real(dp) :: calculate_temperature_Pconst
        class(thermophysical_properties), intent(in) :: this
        real(dp), intent(in) :: temperature,h_s
        real(dp), dimension(:), intent(in) :: molar_fractions

        real(dp) :: temp,eps,residual
        integer :: iter,max_iter

        temp = this%clip_temperature(temperature)
        calculate_temperature_Pconst = temp
        max_iter = 100
        eps = 1e-08

        do iter=1,max_iter
            calculate_temperature_Pconst = calculate_temperature_Pconst + &
                (h_s-(this%mixture_enthalpy_molar(calculate_temperature_Pconst,molar_fractions)- &
                this%mixture_enthalpy_molar(T_ref,molar_fractions)))/ &
                this%mixture_cp_molar(calculate_temperature_Pconst,molar_fractions)
            residual = h_s-this%mixture_enthalpy_molar( &
                calculate_temperature_Pconst,molar_fractions)+ &
                this%mixture_enthalpy_molar(T_ref,molar_fractions)
            if (abs(residual)<eps) exit
        end do
    end function calculate_temperature_Pconst


    subroutine change_field_units_mole_to_dimless(this,Y)
        class(thermophysical_properties), intent(in) :: this
        type(field_vector_cons), intent(inout) :: Y

        real(dp) :: mass_summ
        integer, dimension(3,2) :: loop
        integer :: species_number
        integer :: i,j,k,dim,specie_number

        species_number = size(Y%pr)
        do dim=1,3
            loop(dim,1)=lbound(Y%pr(1)%cells,dim)
            loop(dim,2)=ubound(Y%pr(1)%cells,dim)
        end do

        do k=loop(3,1),loop(3,2)
        do j=loop(2,1),loop(2,2)
        do i=loop(1,1),loop(1,2)
            mass_summ=0.0_dp
            do specie_number=1,species_number
                mass_summ=mass_summ+Y%pr(specie_number)%cells(i,j,k)* &
                    this%molar_masses(specie_number)
            end do
            do specie_number=1,species_number
                Y%pr(specie_number)%cells(i,j,k)= &
                    Y%pr(specie_number)%cells(i,j,k)* &
                    this%molar_masses(specie_number)/mass_summ
            end do
        end do
        end do
        end do
    end subroutine change_field_units_mole_to_dimless


    subroutine change_cell_units_mole_to_dimless(this,X)
        class(thermophysical_properties), intent(in) :: this
        real(dp), dimension(:), intent(inout) :: X

        real(dp) :: mass_summ
        integer :: specie_number

        mass_summ=0.0_dp
        do specie_number=1,size(X)
            mass_summ=mass_summ+X(specie_number)*this%molar_masses(specie_number)
        end do
        do specie_number=1,size(X)
            X(specie_number)=X(specie_number)*this%molar_masses(specie_number)/mass_summ
        end do
    end subroutine change_cell_units_mole_to_dimless


    real(dp) function get_specie_molar_fraction(this,Y,specie_index)
        class(thermophysical_properties), intent(in) :: this
        real(dp), dimension(:), intent(in) :: Y
        integer, intent(in) :: specie_index

        real(dp) :: mass_summ
        integer :: specie_number

        mass_summ=0.0_dp
        do specie_number=1,size(Y)
            mass_summ=mass_summ+Y(specie_number)/this%molar_masses(specie_number)
        end do

        get_specie_molar_fraction=Y(specie_index)/ &
            this%molar_masses(specie_index)/mass_summ
    end function get_specie_molar_fraction


    pure real(dp) function calculate_omega(this,T_red,l,s) result(omega)
        class(thermophysical_properties), intent(in) :: this
        real(dp), intent(in) :: T_red
        integer, intent(in) :: l,s

        if (T_red<70.0_dp) then
            omega = this%omega_c(l,s,1)/(T_red**this%omega_c(l,s,2)) + &
                this%omega_c(l,s,3)/exp(this%omega_c(l,s,4)*T_red) + &
                this%omega_c(l,s,5)/exp(this%omega_c(l,s,6)*T_red) + &
                this%omega_c(l,s,7)/exp(this%omega_c(l,s,8)*T_red)
        else
            omega = this%omega_c(l,s,1)/(T_red**this%omega_c(l,s,2))
        end if
    end function calculate_omega


    pure real(dp) function specie_dynamic_viscosity(this,T,specie_number) result(mu_specie)
        class(thermophysical_properties), intent(in) :: this
        real(dp), intent(in) :: T
        integer, intent(in) :: specie_number
        real(dp) :: T_red,omega_2_2

        mu_specie=0.0_dp
        if (specie_number<1 .or. specie_number>size(this%molar_masses)) return
        if (T<=0.0_dp) return
        if (this%molar_masses(specie_number)<=0.0_dp) return
        if (this%potential_well_depth(specie_number)<=0.0_dp) return
        if (this%collision_diameter(specie_number)<=0.0_dp) return

        T_red=T/this%potential_well_depth(specie_number)
        omega_2_2=this%calculate_omega(T_red,2,2)
        if (omega_2_2<=0.0_dp) return

        mu_specie=8.4417e-07_dp*sqrt(this%molar_masses(specie_number))*sqrt(T)/ &
            (this%collision_diameter(specie_number)**2*omega_2_2)
    end function specie_dynamic_viscosity


    real(dp) function mixture_dynamic_viscosity(this,T,Y_mass,m_mix) result(mu_mix)
        class(thermophysical_properties), intent(in) :: this
        real(dp), intent(in) :: T
        real(dp), dimension(:), intent(in) :: Y_mass
        real(dp), intent(in), optional :: m_mix

        real(dp) :: mixture_molar_mass,mole_fraction,mu_specie
        real(dp) :: sum_arithmetic,sum_harmonic
        integer :: specie_number,species_number

        mu_mix=0.0_dp
        species_number=min(size(Y_mass),size(this%molar_masses))
        if (species_number<=0 .or. T<=0.0_dp) return

        if (present(m_mix)) then
            mixture_molar_mass=m_mix
        else
            mixture_molar_mass=this%mixture_molar_mass_from_mass_fractions(Y_mass)
        end if
        if (mixture_molar_mass<=0.0_dp) return

        sum_arithmetic=0.0_dp
        sum_harmonic=0.0_dp
        do specie_number=1,species_number
            if (this%molar_masses(specie_number)<=0.0_dp) cycle
            mole_fraction=Y_mass(specie_number)*mixture_molar_mass/ &
                this%molar_masses(specie_number)
            if (mole_fraction==0.0_dp) cycle
            mu_specie=this%specie_dynamic_viscosity(T,specie_number)
            if (mu_specie<=0.0_dp) cycle
            sum_arithmetic=sum_arithmetic+mole_fraction*mu_specie
            sum_harmonic=sum_harmonic+mole_fraction/mu_specie
        end do
        if (sum_harmonic>1.0e-10_dp) then
            mu_mix=0.5_dp*(sum_arithmetic+1.0_dp/sum_harmonic)
        end if
    end function mixture_dynamic_viscosity


    pure real(dp) function specie_thermal_conductivity(this,T,specie_number) result(lambda_specie)
        class(thermophysical_properties), intent(in) :: this
        real(dp), intent(in) :: T
        integer, intent(in) :: specie_number

        real(dp) :: T_red,omega_2_2,cp_molar,cv_molar,conductivity_constant

        lambda_specie=0.0_dp
        if (specie_number<1 .or. specie_number>size(this%molar_masses)) return
        if (T<=0.0_dp) return
        if (this%molar_masses(specie_number)<=0.0_dp) return
        if (this%potential_well_depth(specie_number)<=0.0_dp) return
        if (this%collision_diameter(specie_number)<=0.0_dp) return

        T_red=T/this%potential_well_depth(specie_number)
        omega_2_2=this%calculate_omega(T_red,2,2)
        if (omega_2_2<=0.0_dp) return

        cp_molar=this%specie_cp_molar(T,specie_number)
        cv_molar=cp_molar-r_gase_J
        conductivity_constant=2.6319e-05_dp* &
            sqrt(1.0_dp/this%molar_masses(specie_number))/ &
            this%collision_diameter(specie_number)**2
        lambda_specie=conductivity_constant*sqrt(T)/omega_2_2* &
            (4.0_dp*cv_molar/(15.0_dp*r_gase_J)+3.0_dp/5.0_dp)
    end function specie_thermal_conductivity


    real(dp) function mixture_thermal_conductivity(this,T,Y_mass,m_mix) result(lambda_mix)
        class(thermophysical_properties), intent(in) :: this
        real(dp), intent(in) :: T
        real(dp), dimension(:), intent(in) :: Y_mass
        real(dp), intent(in), optional :: m_mix

        real(dp) :: mixture_molar_mass,mole_fraction,lambda_specie
        real(dp) :: sum_arithmetic,sum_harmonic
        integer :: specie_number,species_number

        lambda_mix=0.0_dp
        species_number=min(size(Y_mass),size(this%molar_masses))
        if (species_number<=0 .or. T<=0.0_dp) return

        if (present(m_mix)) then
            mixture_molar_mass=m_mix
        else
            mixture_molar_mass=this%mixture_molar_mass_from_mass_fractions(Y_mass)
        end if
        if (mixture_molar_mass<=0.0_dp) return

        sum_arithmetic=0.0_dp
        sum_harmonic=0.0_dp
        do specie_number=1,species_number
            if (this%molar_masses(specie_number)<=0.0_dp) cycle
            mole_fraction=Y_mass(specie_number)*mixture_molar_mass/ &
                this%molar_masses(specie_number)
            if (mole_fraction==0.0_dp) cycle
            lambda_specie=this%specie_thermal_conductivity(T,specie_number)
            if (lambda_specie<=0.0_dp) cycle
            sum_arithmetic=sum_arithmetic+mole_fraction*lambda_specie
            sum_harmonic=sum_harmonic+mole_fraction/lambda_specie
        end do
        if (sum_harmonic>1.0e-10_dp) then
            lambda_mix=0.5_dp*(sum_arithmetic+1.0_dp/sum_harmonic)
        end if
    end function mixture_thermal_conductivity


    pure real(dp) function binary_diffusion_coefficient(this,T,p, &
            specie_number1,specie_number2) result(D_binary)
        class(thermophysical_properties), intent(in) :: this
        real(dp), intent(in) :: T,p
        integer, intent(in) :: specie_number1,specie_number2

        real(dp) :: reduced_temperature,omega_1_1
        integer :: i1,i2

        D_binary=0.0_dp
        if (T<=0.0_dp .or. p<=0.0_dp) return
        if (specie_number1==specie_number2) return
        if (specie_number1<1 .or. specie_number1>size(this%molar_masses)) return
        if (specie_number2<1 .or. specie_number2>size(this%molar_masses)) return

        i1=min(specie_number1,specie_number2)
        i2=max(specie_number1,specie_number2)
        if (this%binary_diffusivity_constant(i1,i2)<=0.0_dp) return
        if (this%potential_well_depth(i1)<=0.0_dp .or. &
            this%potential_well_depth(i2)<=0.0_dp) return

        reduced_temperature=T/sqrt(this%potential_well_depth(i1)* &
            this%potential_well_depth(i2))
        omega_1_1=this%calculate_omega(reduced_temperature,1,1)
        if (omega_1_1<=0.0_dp) return

        D_binary=this%binary_diffusivity_constant(i1,i2)*T**1.5_dp/ &
            (p*omega_1_1)
    end function binary_diffusion_coefficient


    subroutine mixture_averaged_diffusion_coefficients(this,T,p,Y_mass,m_mix,D_mix)
        class(thermophysical_properties), intent(in) :: this
        real(dp), intent(in) :: T,p,m_mix
        real(dp), dimension(:), intent(in) :: Y_mass
        real(dp), dimension(:), intent(out) :: D_mix

        real(dp), dimension(size(Y_mass),size(Y_mass)) :: D_binary
        real(dp) :: accumulation_value,Dij
        integer :: species_number,specie_number1,specie_number2

        D_mix=0.0_dp
        species_number=min(size(Y_mass),size(D_mix),size(this%molar_masses))
        if (species_number<=0 .or. T<=0.0_dp .or. p<=0.0_dp .or. &
            m_mix<=0.0_dp) return

        D_binary=0.0_dp
        do specie_number1=1,species_number-1
            do specie_number2=specie_number1+1,species_number
                D_binary(specie_number1,specie_number2)= &
                    this%binary_diffusion_coefficient( &
                        T,p,specie_number1,specie_number2)
            end do
        end do

        do specie_number1=1,species_number
            if (this%molar_masses(specie_number1)<=0.0_dp) cycle
            accumulation_value=0.0_dp
            do specie_number2=1,species_number
                if (specie_number2==specie_number1) cycle
                if (this%molar_masses(specie_number2)<=0.0_dp) cycle
                if (specie_number2>specie_number1) then
                    Dij=D_binary(specie_number1,specie_number2)
                else
                    Dij=D_binary(specie_number2,specie_number1)
                end if
                if (Dij<=0.0_dp) cycle
                accumulation_value=accumulation_value + &
                    (1.0_dp-Y_mass(specie_number1))* &
                    this%molar_masses(specie_number1)*Y_mass(specie_number2)/ &
                    (Dij*this%molar_masses(specie_number2))
                accumulation_value=accumulation_value + &
                    Y_mass(specie_number1)*Y_mass(specie_number2)/Dij
            end do
            if (accumulation_value>1.0e-10_dp) then
                D_mix(specie_number1)= &
                    (1.0_dp-Y_mass(specie_number1))* &
                    this%molar_masses(specie_number1)/(accumulation_value*m_mix)
            end if
        end do
    end subroutine mixture_averaged_diffusion_coefficients


    subroutine reduced_thermal_diffusion_coefficients(this,T,Y_mass,m_mix, &
            target_indices,alpha_soret,D_T)
        class(thermophysical_properties), intent(in) :: this
        real(dp), intent(in) :: T,m_mix
        real(dp), dimension(:), intent(in) :: Y_mass
        integer, dimension(:), intent(in) :: target_indices
        real(dp), dimension(:), intent(in) :: alpha_soret
        real(dp), dimension(:), intent(out) :: D_T

        real(dp), parameter :: schlup_temperature_scale=37.15_dp
        real(dp) :: mu_mixture,mu_specie,phi_i_m
        real(dp) :: specie_molar_mass,specie_mass_fraction,specie_mole_fraction
        real(dp) :: reduced_temperature,omega_1_1,omega_1_2,C_star
        real(dp) :: molar_mass_specie_g,molar_mass_mixture_g
        integer :: target,specie_number,number_of_targets

        D_T=0.0_dp
        if (T<=0.0_dp .or. m_mix<=0.0_dp) return

        mu_mixture=this%mixture_dynamic_viscosity(T,Y_mass,m_mix)
        if (mu_mixture<=0.0_dp) return

        number_of_targets=min(size(target_indices),size(alpha_soret))
        do target=1,number_of_targets
            specie_number=target_indices(target)
            if (specie_number<1 .or. specie_number> &
                min(size(Y_mass),size(D_T),size(this%molar_masses))) cycle

            specie_molar_mass=this%molar_masses(specie_number)
            specie_mass_fraction=Y_mass(specie_number)
            if (specie_molar_mass<=0.0_dp .or. &
                specie_mass_fraction<=0.0_dp) cycle

            specie_mole_fraction=specie_mass_fraction*m_mix/specie_molar_mass
            mu_specie=this%specie_dynamic_viscosity(T,specie_number)
            if (mu_specie<=0.0_dp) cycle

            phi_i_m=1.065_dp/(2.0_dp*sqrt(2.0_dp))* &
                (1.0_dp+sqrt(mu_specie/mu_mixture)* &
                (m_mix/specie_molar_mass)**0.25_dp)**2
            if (phi_i_m<=0.0_dp) cycle

            molar_mass_specie_g=1000.0_dp*specie_molar_mass
            molar_mass_mixture_g=1000.0_dp*m_mix
            reduced_temperature=T/schlup_temperature_scale* &
                (molar_mass_specie_g*molar_mass_mixture_g)**(-0.29_dp)
            if (reduced_temperature<=0.0_dp) cycle

            omega_1_1=this%calculate_omega(reduced_temperature,1,1)
            omega_1_2=this%calculate_omega(reduced_temperature,1,2)
            if (omega_1_1<=0.0_dp) cycle
            C_star=omega_1_2/omega_1_1

            D_T(specie_number)=-alpha_soret(target)*15.0_dp/4.0_dp* &
                specie_mole_fraction*mu_specie/phi_i_m* &
                (1.2_dp*C_star-1.0_dp)*(1.0_dp-specie_mass_fraction)
        end do
    end subroutine reduced_thermal_diffusion_coefficients


    pure real(dp) function clip_temperature(this,T) result(T_clipped)
        class(thermophysical_properties), intent(in) :: this
        real(dp), intent(in) :: T
        T_clipped=min(T,tables_temperature_ceiling)
    end function clip_temperature


    real(dp) function mixture_molar_mass_from_mass_fractions(this,Y_mass) result(m_mix)
        class(thermophysical_properties), intent(in) :: this
        real(dp), dimension(:), intent(in) :: Y_mass
        integer :: spec
        real(dp) :: denom

        denom=0.0_dp
        do spec=1,size(Y_mass)
            denom=denom+Y_mass(spec)/this%molar_masses(spec)
        end do
        m_mix=1.0_dp/denom
    end function mixture_molar_mass_from_mass_fractions


    subroutine mole_fractions_from_mass_fractions(this,Y_mass,X_mole)
        class(thermophysical_properties), intent(in) :: this
        real(dp), dimension(:), intent(in) :: Y_mass
        real(dp), dimension(:), intent(out) :: X_mole
        integer :: spec,species_number
        real(dp) :: m_mix

        species_number=size(Y_mass)
        if (size(X_mole)/=species_number) &
            error stop 'mole_fractions_from_mass_fractions: size mismatch'
        m_mix=this%mixture_molar_mass_from_mass_fractions(Y_mass)
        do spec=1,species_number
            X_mole(spec)=Y_mass(spec)*m_mix/this%molar_masses(spec)
        end do
    end subroutine mole_fractions_from_mass_fractions


    real(dp) function temperature_from_pressure_density_Y(this,pressure,density,Y_mass) result(temperature)
        class(thermophysical_properties), intent(in) :: this
        real(dp), intent(in) :: pressure,density
        real(dp), dimension(:), intent(in) :: Y_mass
        real(dp) :: m_mix

        m_mix=this%mixture_molar_mass_from_mass_fractions(Y_mass)
        temperature=pressure*m_mix/(density*r_gase_J)
    end function temperature_from_pressure_density_Y


    real(dp) function mixture_specific_entropy(this,temperature,pressure, &
            Y_mass,pressure_ref) result(s_mass)
        class(thermophysical_properties), intent(in) :: this
        real(dp), intent(in) :: temperature,pressure,pressure_ref
        real(dp), dimension(:), intent(in) :: Y_mass
        integer :: spec,species_number
        real(dp) :: m_mix,x_spec,s_molar

        species_number=size(Y_mass)
        m_mix=this%mixture_molar_mass_from_mass_fractions(Y_mass)
        s_molar=-r_gase_J*log(pressure/pressure_ref)
        do spec=1,species_number
            x_spec=Y_mass(spec)*m_mix/this%molar_masses(spec)
            if (x_spec<0.0_dp) &
                error stop 'mixture_specific_entropy: negative mole fraction'
            if (x_spec>0.0_dp) then
                s_molar=s_molar+x_spec* &
                    this%specie_entropy_molar(temperature,spec)- &
                    r_gase_J*x_spec*log(x_spec)
            end if
        end do
        s_mass=s_molar/m_mix
    end function mixture_specific_entropy


    real(dp) function temperature_from_entropy_pressure_Y(this,s_target, &
            pressure,Y_mass,T_guess,T_bracket_low,T_bracket_high, &
            pressure_ref) result(temperature)
        class(thermophysical_properties), intent(in) :: this
        real(dp), intent(in) :: s_target,pressure,T_guess
        real(dp), dimension(:), intent(in) :: Y_mass
        real(dp), intent(in) :: T_bracket_low,T_bracket_high,pressure_ref

        integer :: iter,spec,species_number
        real(dp) :: T_low,T_high,T_loc,s_low,s_high,s_cur,residual
        real(dp) :: m_mix,x_spec,cp_molar,ds_dT,dT

        species_number=size(Y_mass)
        m_mix=this%mixture_molar_mass_from_mass_fractions(Y_mass)
        T_low=T_bracket_low
        T_high=T_bracket_high
        s_low=this%mixture_specific_entropy(T_low,pressure,Y_mass,pressure_ref)
        s_high=this%mixture_specific_entropy(T_high,pressure,Y_mass,pressure_ref)

        if (s_target>s_high) then
            write(*,*) 'Thermophysical entropy inversion failed:'
            write(*,*) '  s_target, s_low, s_high = ',s_target,s_low,s_high
            write(*,*) '  pressure, T_guess = ',pressure,T_guess
            write(*,*) '  T_low, T_high = ',T_low,T_high
            print *, 'temperature_from_entropy_pressure_Y: entropy outside temperature bracket'
            temperature=T_high
            return
        end if
        if (s_target<s_low) then
            write(*,*) 'Thermophysical entropy inversion failed:'
            write(*,*) '  s_target, s_low, s_high = ',s_target,s_low,s_high
            write(*,*) '  pressure, T_guess = ',pressure,T_guess
            write(*,*) '  T_low, T_high = ',T_low,T_high
            print *, 'temperature_from_entropy_pressure_Y: entropy outside temperature bracket'
            temperature=T_low
            return
        end if
        if (T_guess<T_low .or. T_guess>T_high) then
            write(*,*) 'Thermophysical entropy inversion failed:'
            write(*,*) '  T_guess, T_low, T_high = ',T_guess,T_low,T_high
            error stop 'temperature_from_entropy_pressure_Y: initial temperature outside bracket'
        end if

        T_loc=T_guess
        do iter=1,60
            s_cur=this%mixture_specific_entropy(T_loc,pressure,Y_mass,pressure_ref)
            residual=s_cur-s_target
            if (abs(residual)<=1.0e-10_dp*max(1.0_dp,abs(s_target))) exit
            if (residual>0.0_dp) then
                T_high=T_loc
            else
                T_low=T_loc
            end if

            cp_molar=0.0_dp
            do spec=1,species_number
                x_spec=Y_mass(spec)*m_mix/this%molar_masses(spec)
                if (x_spec<0.0_dp) &
                    error stop 'temperature_from_entropy_pressure_Y: negative mole fraction'
                cp_molar=cp_molar+x_spec*this%specie_cp_molar(T_loc,spec)
            end do
            ds_dT=cp_molar/(m_mix*T_loc)
            dT=-residual/ds_dT
            if ((T_loc+dT<=T_low) .or. (T_loc+dT>=T_high) .or. &
                (abs(dT)>0.5_dp*T_loc)) then
                T_loc=0.5_dp*(T_low+T_high)
            else
                T_loc=T_loc+dT
            end if
        end do
        temperature=T_loc
    end function temperature_from_entropy_pressure_Y


    real(dp) function specific_internal_energy_from_temperature_Y(this, &
            temperature,Y_mass) result(e_mass)
        class(thermophysical_properties), intent(in) :: this
        real(dp), intent(in) :: temperature
        real(dp), dimension(:), intent(in) :: Y_mass
        integer :: spec

        e_mass=0.0_dp
        do spec=1,size(Y_mass)
            e_mass=e_mass+Y_mass(spec)* &
                this%specie_internal_energy_molar(temperature,spec)/ &
                this%molar_masses(spec)
        end do
    end function specific_internal_energy_from_temperature_Y


    real(dp) function sensible_internal_energy_molar_from_temperature_X( &
            this,temperature,X_mole) result(e_sens_molar)
        class(thermophysical_properties), intent(in) :: this
        real(dp), intent(in) :: temperature
        real(dp), dimension(:), intent(in) :: X_mole

        e_sens_molar=this%mixture_internal_energy_molar(temperature,X_mole)- &
            this%mixture_internal_energy_molar(T_ref,X_mole)
    end function sensible_internal_energy_molar_from_temperature_X


    real(dp) function sensible_energy_parameter_from_temperature_Y(this, &
            temperature,Y_mass) result(ksi_sens)
        class(thermophysical_properties), intent(in) :: this
        real(dp), intent(in) :: temperature
        real(dp), dimension(:), intent(in) :: Y_mass
        real(dp), dimension(size(Y_mass)) :: X_mole
        real(dp) :: e_sens_molar

        call this%mole_fractions_from_mass_fractions(Y_mass,X_mole)
        e_sens_molar=this%sensible_internal_energy_molar_from_temperature_X( &
            temperature,X_mole)
        ksi_sens=e_sens_molar/(r_gase_J*T_ref)
    end function sensible_energy_parameter_from_temperature_Y


    real(dp) function temperature_from_sensible_energy_parameter_Y(this, &
            ksi_sens,Y_mass,T_guess) result(temperature)
        class(thermophysical_properties), intent(in) :: this
        real(dp), intent(in) :: ksi_sens,T_guess
        real(dp), dimension(:), intent(in) :: Y_mass
        real(dp), dimension(size(Y_mass)) :: X_mole
        real(dp) :: e_target_molar

        call this%mole_fractions_from_mass_fractions(Y_mass,X_mole)
        e_target_molar=ksi_sens*r_gase_J*T_ref + &
            this%mixture_internal_energy_molar(T_ref,X_mole)- &
            this%mixture_enthalpy_molar(T_ref,X_mole)
        temperature=this%calculate_temperature(T_guess,e_target_molar,X_mole)
    end function temperature_from_sensible_energy_parameter_Y


    subroutine pressure_derivatives_density_energy_Y(this,density,pressure, &
            Y_mass,dp_deps,dp_drho,sound_speed)
        class(thermophysical_properties), intent(in) :: this
        real(dp), intent(in) :: density,pressure
        real(dp), dimension(:), intent(in) :: Y_mass
        real(dp), intent(out) :: dp_deps,dp_drho
        real(dp), intent(out), optional :: sound_speed

        integer :: spec,species_number
        real(dp) :: temperature,m_mix,cp_molar,cv_molar,x_spec,c2

        species_number=size(Y_mass)
        m_mix=this%mixture_molar_mass_from_mass_fractions(Y_mass)
        temperature=pressure*m_mix/(density*r_gase_J)
        cp_molar=0.0_dp
        do spec=1,species_number
            x_spec=Y_mass(spec)*m_mix/this%molar_masses(spec)
            cp_molar=cp_molar+x_spec*this%specie_cp_molar(temperature,spec)
        end do
        cv_molar=cp_molar-r_gase_J
        dp_deps=density*r_gase_J/cv_molar
        dp_drho=pressure/density
        if (present(sound_speed)) then
            c2=dp_drho+dp_deps*pressure/(density*density)
            sound_speed=sqrt(c2)
        end if
    end subroutine pressure_derivatives_density_energy_Y


    subroutine pressure_derivatives_density_energy_composition_Y(this, &
            density,pressure,Y_mass,dp_deps,dp_drho,dp_dY,sound_speed)
        class(thermophysical_properties), intent(in) :: this
        real(dp), intent(in) :: density,pressure
        real(dp), dimension(:), intent(in) :: Y_mass
        real(dp), intent(out) :: dp_deps,dp_drho
        real(dp), dimension(:), intent(out) :: dp_dY
        real(dp), intent(out), optional :: sound_speed

        integer :: spec,species_number
        real(dp) :: temperature,m_mix,cp_molar,cv_molar,x_spec,c2,e_mass_spec

        species_number=size(Y_mass)
        if (size(dp_dY)/=species_number) &
            error stop 'pressure_derivatives_density_energy_composition_Y: dp_dY size mismatch'

        m_mix=this%mixture_molar_mass_from_mass_fractions(Y_mass)
        temperature=pressure*m_mix/(density*r_gase_J)
        cp_molar=0.0_dp
        do spec=1,species_number
            x_spec=Y_mass(spec)*m_mix/this%molar_masses(spec)
            cp_molar=cp_molar+x_spec*this%specie_cp_molar(temperature,spec)
        end do
        cv_molar=cp_molar-r_gase_J
        dp_deps=density*r_gase_J/cv_molar
        dp_drho=pressure/density

        do spec=1,species_number
            e_mass_spec=this%specie_internal_energy_molar(temperature,spec)/ &
                this%molar_masses(spec)
            dp_dY(spec)=density*(r_gase_J*temperature/this%molar_masses(spec)- &
                (r_gase_J/cv_molar)*e_mass_spec)
        end do

        if (present(sound_speed)) then
            c2=dp_drho+dp_deps*pressure/(density*density)
            sound_speed=sqrt(c2)
        end if
    end subroutine pressure_derivatives_density_energy_composition_Y


    subroutine pressure_source_rate_from_conservative_sources(this,density, &
            pressure,E_total,velocity,Y_mass,S_rho,S_momentum,S_rhoE, &
            S_rhoY,dp_dt,de_dt,dY_dt)
        class(thermophysical_properties), intent(in) :: this
        real(dp), intent(in) :: density,pressure,E_total,S_rho,S_rhoE
        real(dp), dimension(:), intent(in) :: velocity,Y_mass,S_momentum,S_rhoY
        real(dp), intent(out) :: dp_dt,de_dt
        real(dp), dimension(:), intent(out) :: dY_dt

        integer :: spec,species_number
        real(dp) :: u2,dp_deps,dp_drho
        real(dp), dimension(size(Y_mass)) :: dp_dY

        species_number=size(Y_mass)
        if (size(dY_dt)/=species_number) &
            error stop 'pressure_source_rate_from_conservative_sources: dY size mismatch'

        u2=sum(velocity(:)*velocity(:))
        de_dt=(S_rhoE+(u2-E_total)*S_rho-sum(velocity(:)*S_momentum(:)))/density
        do spec=1,species_number
            dY_dt(spec)=(S_rhoY(spec)-Y_mass(spec)*S_rho)/density
        end do
        call this%pressure_derivatives_density_energy_composition_Y( &
            density,pressure,Y_mass,dp_deps,dp_drho,dp_dY)
        dp_dt=dp_drho*S_rho+dp_deps*de_dt+sum(dp_dY(:)*dY_dt(:))
    end subroutine pressure_source_rate_from_conservative_sources


!# Thermo I/O functions

    function get_a_coeffs(specie_name,thermo_data_file_unit) result(a_coeffs)
        ! Read a standard four-record NASA-7/CHEMKIN thermo entry.  Coefficients
        ! occupy fixed 15-character fields in columns 1:75; comments after the
        ! record marker therefore cannot interfere with conversion.
        character(len=*), intent(in) :: specie_name
        integer, intent(in) :: thermo_data_file_unit
        real(dp), dimension(7,2) :: a_coeffs

        character(len=256) :: header,line2,line3,line4,adjusted_header
        character(len=64) :: record_species
        integer :: io_status
        logical :: found

        a_coeffs=0.0_dp
        found=.false.
        rewind(thermo_data_file_unit)

        do
            read(thermo_data_file_unit,'(A)',iostat=io_status) header
            if (io_status<0) exit
            if (io_status>0) then
                error stop 'thermophysical_properties: NASA thermo header read failure'
            end if

            adjusted_header=adjustl(header)
            if (len_trim(adjusted_header)==0) cycle
            if (adjusted_header(1:1)=='!') cycle

            record_species=''
            read(adjusted_header,*,iostat=io_status) record_species
            if (io_status/=0) cycle
            if (trim(record_species)=='END') exit
            if (trim(record_species)/=trim(specie_name)) cycle

            read(thermo_data_file_unit,'(A)',iostat=io_status) line2
            if (io_status/=0) error stop 'thermophysical_properties: NASA thermo line 2 read failure'
            read(thermo_data_file_unit,'(A)',iostat=io_status) line3
            if (io_status/=0) error stop 'thermophysical_properties: NASA thermo line 3 read failure'
            read(thermo_data_file_unit,'(A)',iostat=io_status) line4
            if (io_status/=0) error stop 'thermophysical_properties: NASA thermo line 4 read failure'

            ! High-temperature coefficients a1..a5.
            call read_nasa_number(line2( 1:15),a_coeffs(1,2),specie_name)
            call read_nasa_number(line2(16:30),a_coeffs(2,2),specie_name)
            call read_nasa_number(line2(31:45),a_coeffs(3,2),specie_name)
            call read_nasa_number(line2(46:60),a_coeffs(4,2),specie_name)
            call read_nasa_number(line2(61:75),a_coeffs(5,2),specie_name)

            ! High-temperature a6,a7, then low-temperature a1..a3.
            call read_nasa_number(line3( 1:15),a_coeffs(6,2),specie_name)
            call read_nasa_number(line3(16:30),a_coeffs(7,2),specie_name)
            call read_nasa_number(line3(31:45),a_coeffs(1,1),specie_name)
            call read_nasa_number(line3(46:60),a_coeffs(2,1),specie_name)
            call read_nasa_number(line3(61:75),a_coeffs(3,1),specie_name)

            ! Remaining low-temperature a4..a7.
            call read_nasa_number(line4( 1:15),a_coeffs(4,1),specie_name)
            call read_nasa_number(line4(16:30),a_coeffs(5,1),specie_name)
            call read_nasa_number(line4(31:45),a_coeffs(6,1),specie_name)
            call read_nasa_number(line4(46:60),a_coeffs(7,1),specie_name)

            found=.true.
            exit
        end do

        rewind(thermo_data_file_unit)

        if (.not.found) then
            write(*,'(A,A)') 'Thermophysical properties: NASA thermo entry not found for ', &
                trim(specie_name)
            error stop 'thermophysical_properties: missing NASA thermo entry'
        end if
    end function get_a_coeffs


    subroutine read_nasa_number(field,value,specie_name)
        character(len=*), intent(in) :: field
        real(dp), intent(out) :: value
        character(len=*), intent(in) :: specie_name
        integer :: io_status

        read(field,*,iostat=io_status) value
        if (io_status/=0) then
            write(*,'(A,A,A,A)') 'Thermophysical properties: invalid NASA coefficient for ', &
                trim(specie_name), ': ', field
            error stop 'thermophysical_properties: NASA coefficient conversion failure'
        end if
    end subroutine read_nasa_number


    real(dp) function get_specie_molar_mass_from_thermo( &
            specie_name,thermo_data_file_unit) result(molar_mass)
        ! Derive M_k from the standard NASA/CHEMKIN elemental-composition
        ! fields in columns 25:44 of the species header.  The species token is
        ! used only as an identifier.  Therefore aliases/state labels such as
        ! TXCH2, O2X, OH*, CH2(S), etc. do not have to be interpreted as
        ! molecular formulae.
        character(len=*), intent(in) :: specie_name
        integer, intent(in) :: thermo_data_file_unit

        integer, parameter :: thermo_line_length=256
        integer, parameter :: composition_slots=4
        character(len=thermo_line_length) :: line,adjusted_line
        character(len=64) :: record_species
        character(len=5) :: composition_field
        character(len=2) :: element_symbol
        real(dp) :: atom_count,element_mass
        integer :: io_status,slot,first,last
        logical :: record_found,composition_found

        molar_mass=0.0_dp
        record_found=.false.
        composition_found=.false.

        rewind(thermo_data_file_unit)
        do
            read(thermo_data_file_unit,'(A)',iostat=io_status) line
            if (io_status<0) exit
            if (io_status>0) then
                error stop 'thermophysical_properties: thermo header read failure'
            end if

            adjusted_line=adjustl(line)
            if (len_trim(adjusted_line)==0) cycle
            if (adjusted_line(1:1)=='!') cycle

            record_species=''
            read(adjusted_line,*,iostat=io_status) record_species
            if (io_status/=0) cycle
            if (trim(record_species)/=trim(specie_name)) cycle

            record_found=.true.

            ! NASA-7/CHEMKIN species header:
            !   cols  1:24  species/date/descriptor
            !   cols 25:44  four fields (element symbol + atom count)
            !   col      45  phase
            do slot=1,composition_slots
                first=25+5*(slot-1)
                last=first+4
                composition_field=line(first:last)
                element_symbol=adjustl(composition_field(1:2))
                if (len_trim(element_symbol)==0) cycle

                atom_count=0.0_dp
                read(composition_field(3:5),*,iostat=io_status) atom_count
                if (io_status/=0) then
                    write(*,'(A,A,A,A)') &
                        'Thermophysical properties: malformed elemental field for ', &
                        trim(specie_name), ': ', composition_field
                    error stop 'thermophysical_properties: invalid NASA elemental composition'
                end if

                if (abs(atom_count)<=tiny(1.0_dp)) cycle
                element_mass=atomic_molar_mass(element_symbol)
                if (element_mass<=0.0_dp) then
                    write(*,'(A,A,A,A)') &
                        'Thermophysical properties: unsupported element ', &
                        trim(element_symbol), ' in species ', trim(specie_name)
                    error stop 'thermophysical_properties: unsupported chemical element'
                end if

                molar_mass=molar_mass+atom_count*element_mass
                composition_found=.true.
            end do

            exit
        end do

        rewind(thermo_data_file_unit)

        if (.not.record_found) then
            write(*,'(A,A)') &
                'Thermophysical properties: thermo entry not found for species ', &
                trim(specie_name)
            error stop 'thermophysical_properties: species absent from thermo file'
        end if

        if (.not.composition_found .or. molar_mass<=0.0_dp) then
            write(*,'(A,A)') &
                'Thermophysical properties: no valid elemental composition for species ', &
                trim(specie_name)
            error stop 'thermophysical_properties: empty NASA elemental composition'
        end if
    end function get_specie_molar_mass_from_thermo


    real(dp) function atomic_molar_mass(element_symbol) result(mass)
        ! Standard atomic weights in kg/mol.  Values are deliberately defined
        ! once, at elemental level, so every species is internally consistent.
        ! Conventional values are used for the elements most common in
        ! combustion mechanisms; E is the molar mass of an electron for
        ! compatibility with thermo databases that explicitly include charge.
        character(len=*), intent(in) :: element_symbol
        character(len=len(element_symbol)) :: element

        element=to_upper_ascii(adjustl(element_symbol))

        select case(trim(element))
        case('H')
            mass=1.00794e-3_dp
        case('D')
            mass=2.01410177812e-3_dp
        case('HE')
            mass=4.002602e-3_dp
        case('LI')
            mass=6.941e-3_dp
        case('BE')
            mass=9.012182e-3_dp
        case('B')
            mass=10.811e-3_dp
        case('C')
            mass=12.0107e-3_dp
        case('N')
            mass=14.0067e-3_dp
        case('O')
            mass=15.9994e-3_dp
        case('F')
            mass=18.998403163e-3_dp
        case('NE')
            mass=20.1797e-3_dp
        case('NA')
            mass=22.98976928e-3_dp
        case('MG')
            mass=24.3050e-3_dp
        case('AL')
            mass=26.9815385e-3_dp
        case('SI')
            mass=28.0855e-3_dp
        case('P')
            mass=30.973761998e-3_dp
        case('S')
            mass=32.065e-3_dp
        case('CL')
            mass=35.453e-3_dp
        case('AR')
            mass=39.948e-3_dp
        case('K')
            mass=39.0983e-3_dp
        case('CA')
            mass=40.078e-3_dp
        case('FE')
            mass=55.845e-3_dp
        case('BR')
            mass=79.904e-3_dp
        case('I')
            mass=126.90447e-3_dp
        case('E')
            mass=5.48579909065e-7_dp
        case default
            mass=-1.0_dp
        end select
    end function atomic_molar_mass


    pure function to_upper_ascii(value) result(upper_value)
        character(len=*), intent(in) :: value
        character(len=len(value)) :: upper_value
        integer :: i,code

        upper_value=value
        do i=1,len(value)
            code=iachar(upper_value(i:i))
            if (code>=iachar('a') .and. code<=iachar('z')) then
                upper_value(i:i)=achar(code-iachar('a')+iachar('A'))
            end if
        end do
    end function to_upper_ascii


    subroutine initialize_binary_diffusivity_constants(this)
        class(thermophysical_properties), intent(inout) :: this

        real(dp) :: reduced_collision_diameter,inv_reduced_molar_mass
        integer :: specie_number1,specie_number2,species_number

        this%binary_diffusivity_constant=0.0_dp
        species_number=size(this%molar_masses)

        do specie_number1=1,species_number-1
            if (this%molar_masses(specie_number1)<=0.0_dp) cycle
            if (this%collision_diameter(specie_number1)<=0.0_dp) cycle
            do specie_number2=specie_number1+1,species_number
                if (this%molar_masses(specie_number2)<=0.0_dp) cycle
                if (this%collision_diameter(specie_number2)<=0.0_dp) cycle

                reduced_collision_diameter=0.5_dp* &
                    (this%collision_diameter(specie_number1)+ &
                     this%collision_diameter(specie_number2))
                reduced_collision_diameter=reduced_collision_diameter**2
                inv_reduced_molar_mass= &
                    (this%molar_masses(specie_number1)+ &
                     this%molar_masses(specie_number2))/ &
                    (this%molar_masses(specie_number1)* &
                     this%molar_masses(specie_number2))

                this%binary_diffusivity_constant(specie_number1,specie_number2)= &
                    0.00001_dp*0.595_dp*sqrt(inv_reduced_molar_mass)/ &
                    reduced_collision_diameter
                this%binary_diffusivity_constant(specie_number2,specie_number1)= &
                    this%binary_diffusivity_constant(specie_number1,specie_number2)
            end do
        end do
    end subroutine initialize_binary_diffusivity_constants


    subroutine get_transport_properties(potential_well_depth,collision_diameter, &
            specie_name,transport_data_file_unit)
        real(dp), intent(inout) :: potential_well_depth,collision_diameter
        character(len=10), intent(in) :: specie_name
        integer, intent(in) :: transport_data_file_unit

        character(len=100) :: string
        real(dp) :: potential_well_depth_in_file,collision_diameter_in_file
        integer :: dummy,error

        potential_well_depth=0.0_dp
        collision_diameter=0.0_dp

        read(transport_data_file_unit,*,iostat=error) string
        do
            read(transport_data_file_unit,*,iostat=error) string,dummy, &
                potential_well_depth_in_file,collision_diameter_in_file
            if (error/=0) then
                print *, 'Error: unsuccessufull transport properties reading for specie:', &
                    specie_name
                error stop 'thermophysical_properties: transport property read failure'
            end if
            if (string==specie_name) then
                collision_diameter=collision_diameter_in_file
                potential_well_depth=potential_well_depth_in_file
                exit
            end if
        end do

        collision_diameter=collision_diameter*0.1_dp ! Angstrem -> nanometer
        rewind(transport_data_file_unit)
    end subroutine get_transport_properties

end module thermophysical_properties_class
