module thermophysical_properties_class

	use kind_parameters
	use global_data

	use chemical_properties_class

	use field_pointers

	implicit none

	private
	public  thermophysical_properties,thermophysical_properties_pointer ,thermophysical_properties_c

	type 	:: thermophysical_properties_pointer
		type(thermophysical_properties)	,pointer	:: thermo_ptr
    end type
    
	type    :: thermophysical_properties
		character(len=30)                            ::	thermo_data_file_name, transport_data_file_name, molar_masses_data_file_name
		real(dp) ,dimension(:)       ,allocatable    :: potential_well_depth, collision_diameter, molar_masses
		real(dp) ,dimension(:,:,:)   ,allocatable    :: a_coeffs
        real(dp) ,dimension(2,2,8)                   :: omega_c
	contains
		procedure	,private	:: read_properties
		procedure	,private	:: write_properties
		procedure	,private	:: set_properties


		procedure	:: specie_cp_molar
		procedure	:: specie_entropy_molar        
		procedure	:: specie_enthalpy_molar        
        procedure   :: specie_internal_energy_molar
        
		procedure	:: mixture_cp_molar
    	procedure	:: mixture_enthalpy_molar
		procedure	:: mixture_internal_energy_molar
        
		procedure	:: calculate_temperature
        procedure	:: calculate_temperature_Pconst
        procedure	:: calculate_omega
        procedure   :: clip_temperature

        ! CABARET/thermally-perfect ideal-gas mixture helpers.
        ! These routines use the JANAF polynomial data stored in this class and do
        ! not apply pressure/density/gamma floors. They are intended for strict
        ! gas-dynamic algorithms where invalid thermodynamic states must fail
        ! visibly rather than being silently repaired.
        procedure   :: mixture_molar_mass_from_mass_fractions
        procedure   :: mole_fractions_from_mass_fractions
        procedure   :: temperature_from_pressure_density_Y
        procedure   :: mixture_specific_entropy
        procedure   :: temperature_from_entropy_pressure_Y
        procedure   :: pressure_derivatives_density_energy_Y
        procedure   :: pressure_derivatives_density_energy_composition_Y
        procedure   :: pressure_source_rate_from_conservative_sources
        procedure   :: specific_internal_energy_from_temperature_Y

        
		procedure	:: change_field_units_mole_to_dimless
		procedure	:: change_cell_units_mole_to_dimless

		procedure	:: write_log		
	end type

	interface   thermophysical_properties_c
		module procedure   constructor
		module procedure   constructor_file
	end interface

contains

	type(thermophysical_properties) function constructor(chemistry, thermo_data_file_name, transport_data_file_name, molar_masses_data_file_name)
		type(chemical_properties)	,intent(in)	:: chemistry
		character(len=*)			,intent(in)	:: thermo_data_file_name, transport_data_file_name, molar_masses_data_file_name

		integer	:: io_unit
		
		call constructor%set_properties(chemistry,thermo_data_file_name,transport_data_file_name,molar_masses_data_file_name)
		
		open(newunit = io_unit, file = thermophysical_data_file_name, status = 'replace', form = 'formatted', delim = 'quote')
		call constructor%write_properties(io_unit)
		close(io_unit)	
	end function

	type(thermophysical_properties) function constructor_file(chemistry)
		type(chemical_properties)		,intent(in)	:: chemistry

		integer	:: io_unit

		open(newunit = io_unit, file = thermophysical_data_file_name, status = 'old', form = 'formatted')
		call constructor_file%read_properties(chemistry,io_unit)
		close(io_unit)
	end function	
	
	subroutine write_properties(this,thermo_data_unit)
		class(thermophysical_properties)	,intent(in)	:: this
		integer								,intent(in)	:: thermo_data_unit
		
		character(len=30)	:: thermo_data_file_name, transport_data_file_name, molar_masses_data_file_name
		
		namelist /thermo_properties/  thermo_data_file_name, transport_data_file_name, molar_masses_data_file_name
		
		thermo_data_file_name			= this%thermo_data_file_name
		transport_data_file_name		= this%transport_data_file_name
		molar_masses_data_file_name		= this%molar_masses_data_file_name
		
		write(unit = thermo_data_unit, nml = thermo_properties)

	end subroutine
	
	subroutine read_properties(this,chemistry,thermo_data_unit)
		class(thermophysical_properties)	,intent(inout)	:: this
		type(chemical_properties)			,intent(in)	:: chemistry
		integer								,intent(in)	:: thermo_data_unit
		
		character(len=30)	:: thermo_data_file_name, transport_data_file_name, molar_masses_data_file_name
		
		namelist /thermo_properties/  thermo_data_file_name, transport_data_file_name, molar_masses_data_file_name
		
		read(unit = thermo_data_unit, nml = thermo_properties)
		call this%set_properties(chemistry,thermo_data_file_name,transport_data_file_name, molar_masses_data_file_name)	
		
	end subroutine
	
	subroutine set_properties(this,chemistry, thermo_data_file_name, transport_data_file_name, molar_masses_data_file_name)
		class(thermophysical_properties)	,intent(inout)	:: this
		type(chemical_properties)			,intent(in)		:: chemistry
		character(len=*)					,intent(in)		:: thermo_data_file_name, transport_data_file_name, molar_masses_data_file_name
		
		integer :: species_number
		integer :: thermo_data_file_unit, transport_data_file_unit, molar_masses_data_file_unit

		integer :: spec_num

		integer :: error

		species_number	= chemistry%species_number

		this%thermo_data_file_name			=	thermo_data_file_name	
		this%transport_data_file_name		=	transport_data_file_name
		this%molar_masses_data_file_name	=	molar_masses_data_file_name
		
		allocate(this%molar_masses(species_number))
		allocate(this%potential_well_depth(species_number),this%collision_diameter(species_number))
		allocate(this%a_coeffs(7,2,species_number))

		open(newunit = thermo_data_file_unit,       file = trim(task_setup_folder) // trim(fold_sep) // trim(thermophysical_data_folder) // trim(fold_sep) // this%thermo_data_file_name,       status = 'old',    iostat = error)
		open(newunit = transport_data_file_unit,    file = trim(task_setup_folder) // trim(fold_sep) // trim(thermophysical_data_folder) // trim(fold_sep) // this%transport_data_file_name,    status = 'old',    iostat = error)
		open(newunit = molar_masses_data_file_unit, file = trim(task_setup_folder) // trim(fold_sep) // trim(thermophysical_data_folder) // trim(fold_sep) // this%molar_masses_data_file_name, status = 'old',    iostat = error)

		do spec_num = 1, species_number
			this%a_coeffs(:,:,spec_num)  = get_a_coeffs(chemistry%species_names(spec_num), thermo_data_file_unit)
			this%a_coeffs(:,:,spec_num)  = this%a_coeffs(:,:,spec_num) * r_gase_J
			this%molar_masses(spec_num)  = get_specie_molar_mass(chemistry%species_names(spec_num), molar_masses_data_file_unit)
		end do

		do spec_num = 1, chemistry%species_number
			call get_transport_properties(this%potential_well_depth(spec_num), this%collision_diameter(spec_num), chemistry%species_names(spec_num), transport_data_file_unit)
		end do

		close(thermo_data_file_unit)
		close(transport_data_file_unit)
		close(molar_masses_data_file_unit)
        
        !# P.D. Neufeld, A.R. Janzen, R.A. Aziz // J. Chem. Phys. 57, 1100 (1972), doi:10.1063/1.1678363
        
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
        
	end subroutine
	
	subroutine write_log(this,log_unit)
		class(thermophysical_properties)	,intent(in)	:: this	
		integer								,intent(in)	:: log_unit	
		
		write(log_unit,'(A)')		'************************************************************************************* '
		write(log_unit,'(A)')		' Thermophtsical data setup : '
		write(log_unit,'(A,A)')		' Thermo data file        : ',	this%thermo_data_file_name
		write(log_unit,'(A,A)')		' Transport data file     : ',	this%transport_data_file_name
		write(log_unit,'(A,A)')     ' Molar masses data file  : ',	this%molar_masses_data_file_name
		write(log_unit,'(A)')		'************************************************************************************* '
		
    end subroutine
	
    pure real(dp) function specie_cp_molar(this, temperature, specie_number) result(cp_molar)
        class(thermophysical_properties), intent(in) :: this
        real(dp), intent(in)    :: temperature
        integer, intent(in)     :: specie_number
        integer                 :: i, coeff_set
        real(dp)                :: temp
        
        coeff_set = 1
        
        temp = this%clip_temperature(temperature)
        
        if (temp >= 1000.0_dp) coeff_set = 2

        cp_molar = this%a_coeffs(5, coeff_set, specie_number)
        do i = 4, 1, -1
            cp_molar = cp_molar*temp + this%a_coeffs(i, coeff_set, specie_number)
        end do
    end function specie_cp_molar
    
    pure real(dp) function specie_enthalpy_molar(this,temperature, specie_number) result(hs_molar)
		class(thermophysical_properties),   intent(in)  :: this
		real(dp),       intent(in)  :: temperature
		integer,        intent(in)  :: specie_number
        integer                     :: i, coeff_set
        real(dp)                    :: temp

        coeff_set = 1  
        
		temp = this%clip_temperature(temperature)
		
        if (temp >= 1000.0_dp) coeff_set = 2
        
        hs_molar = this%a_coeffs(5,coeff_set,specie_number) / 5.0_dp
        do i = 4, 1, -1
            hs_molar = hs_molar*temp + this%a_coeffs(i,coeff_set,specie_number) / real(i,dp)
        end do
        hs_molar = hs_molar + this%a_coeffs(6,coeff_set,specie_number) / temp
        hs_molar = hs_molar * temp
        
    end function specie_enthalpy_molar

    pure real(dp) function specie_entropy_molar(this, temperature, specie_number) result(s_molar)
        class(thermophysical_properties), intent(in) :: this
        real(dp), intent(in)    :: temperature
        integer, intent(in)     :: specie_number
        integer                 :: i, coeff_set
		real(dp)                :: temp
        
        coeff_set = 1

        temp = this%clip_temperature(temperature)
        
        if (temp >= 1000.0_dp) coeff_set = 2

        s_molar = this%a_coeffs(5, coeff_set, specie_number)/4.0_dp
        do i = 4, 2, -1
            s_molar = s_molar*temp + this%a_coeffs(i, coeff_set, specie_number)/real(i - 1, dp)
        end do
        s_molar = s_molar*temp + this%a_coeffs(7, coeff_set, specie_number) + &
                  this%a_coeffs(1, coeff_set, specie_number)*log(temp)
        
    end function specie_entropy_molar    
    
	pure real(dp) function specie_internal_energy_molar(this, temperature, specie_number) result(e_molar)
        class(thermophysical_properties), intent(in) :: this
        real(dp), intent(in)    :: temperature
        integer, intent(in)     :: specie_number
        real(dp)                :: hs_molar
        integer                 :: i, coeff_set

        hs_molar = this%specie_enthalpy_molar(temperature, specie_number)
        
        ! Ideal-gas molar internal energy: e_k = h_k - R_u T.
        e_molar = hs_molar - r_gase_J*temperature
        
    end function specie_internal_energy_molar    
    
    
	pure real(dp) function mixture_cp_molar(this, temperature, species_molar_fractions) result(cp_molar)
		class(thermophysical_properties),   intent(in)  :: this
		real(dp)                 ,intent(in)    :: temperature
		real(dp) ,dimension(:)   ,intent(in)    :: species_molar_fractions
		real(dp)                                :: specie_cp

        integer :: species_number
		integer :: specie_number

		cp_molar = 0.0_dp
 		species_number = size(species_molar_fractions)

		do specie_number = 1,species_number
			specie_cp = this%specie_cp_molar(temperature, specie_number)
			cp_molar = cp_molar + specie_cp * species_molar_fractions(specie_number)
		end do
    end function mixture_cp_molar
    
    
	pure real(dp) function mixture_enthalpy_molar(this,temperature,species_molar_fractions) result(hs_molar)
		class(thermophysical_properties),   intent(in)  :: this
		real(dp)                ,intent(in) :: temperature        
		real(dp) ,dimension(:)  ,intent(in) :: species_molar_fractions
		real(dp)                            :: specie_enthalpy
    
		integer :: species_number
		real(dp) :: temp
		integer :: i, specie_number

        hs_molar = 0.0_dp
		species_number = size(species_molar_fractions)
		
		do specie_number = 1,species_number
			specie_enthalpy = this%specie_enthalpy_molar(temperature, specie_number)
			hs_molar = hs_molar + specie_enthalpy * species_molar_fractions(specie_number)
		end do
	end function mixture_enthalpy_molar		
	

	pure real(dp) function mixture_internal_energy_molar(this,temperature,species_molar_fractions) result(e_molar)
		class(thermophysical_properties),   intent(in) :: this
		real(dp),                           intent(in) :: temperature
		real(dp) ,dimension(:),             intent(in) :: species_molar_fractions
		real(dp)                            :: specie_internal_energy

        integer :: species_number
		real(dp) :: temp
		integer :: i, specie_number

        e_molar = 0.0_dp       
		species_number = size(species_molar_fractions)
		
		do specie_number = 1,species_number
			specie_internal_energy = this%specie_internal_energy_molar(temperature, specie_number)
			e_molar = e_molar + specie_internal_energy * species_molar_fractions(specie_number)
		end do 
    end function mixture_internal_energy_molar
	
    
	!recursive pure function calculate_temperature(this,temperature,E_full,mol_mix_conc,cv,concentrations)
 !
	!	real(dp)             :: calculate_temperature
 !
	!	class(thermophysical_properties)    ,intent(in) :: this
	!	real(dp)                         ,intent(in) :: temperature, E_full, mol_mix_conc, cv
	!	real(dp) ,dimension(:)           ,intent(in) :: concentrations
	!	!!$OMP threadprivate (calculate_temperature)
 !
	!	real(dp) :: temp		
	!	real(dp) :: specie_dcpdT_T, dcpdT_T
	!	integer :: species_number
	!	integer :: i, specie_number
 !
	!	calculate_temperature   = 0.0_dp
	!	dcpdT_T                 = 0.0_dp
 !
	!	species_number          = size(concentrations)
 !       temp = temperature
 !       if(temperature >= 4200.0_sp) temp = 4200.0_sp
 !
	!	do specie_number = 1,species_number
	!		if(temp < 1000.0) then
	!			specie_dcpdT_T = this%a_coeffs(5,1,specie_number) * 5.0_dp
	!			do i = 4, 1, -1
	!				specie_dcpdT_T = specie_dcpdT_T * temp
	!				specie_dcpdT_T = specie_dcpdT_T + this%a_coeffs(i,1,specie_number) * real(i,dp)
	!			end do
	!		else
	!			specie_dcpdT_T = this%a_coeffs(5,2,specie_number) * 5.0_dp
	!			do i=4,1,-1
	!				specie_dcpdT_T = specie_dcpdT_T * temp
	!				specie_dcpdT_T = specie_dcpdT_T + this%a_coeffs(i,2,specie_number) * real(i,dp)
	!			continue
	!			end do
	!		end if
	!		dcpdT_T = dcpdT_T + specie_dcpdT_T * concentrations(specie_number)
	!	end do
 !
	!	if (temperature >= 4200.0_sp) then
 !           calculate_temperature = E_full / cv
 !       else
	!		calculate_temperature = temp + (E_full - cv * temp) / (dcpdT_T - r_gase_J * mol_mix_conc)
	!	end if
	!	continue
	!end function

	function calculate_temperature(this,temperature,e_i,molar_fractions)

		real(dp)										:: calculate_temperature

		class(thermophysical_properties)    ,intent(in) :: this
		real(dp)                         ,intent(in) :: temperature, e_i
		real(dp) ,dimension(:)           ,intent(in) :: molar_fractions
		!!$OMP threadprivate (calculate_temperature)

		real(dp)		:: temp		
		integer			:: species_number
		integer			:: i, specie_number
		integer			:: iter, max_iter
		real(dp)		:: eps, residual, old_residual, cp, cv, e_int, h_int, h_const, e_const
		real(sp)		:: var
        
        
		species_number          = size(molar_fractions)
		temp = this%clip_temperature(temperature)

        max_iter	= 1000
		eps			= 1e-08

        h_const		= this%mixture_enthalpy_molar(T_ref, molar_fractions)
	    e_const		= this%mixture_internal_energy_molar(T_ref, molar_fractions)
        
        cp		= this%mixture_cp_molar(temp, molar_fractions) 
        cv		= cp - r_gase_J        
       
        !# Newton method
		do iter = 1, max_iter
            e_int	= this%mixture_internal_energy_molar(temp, molar_fractions) - h_const
			residual =  (e_i - e_int) / cp
			temp = temp + (residual) 
            if ( abs(residual) < eps) exit
            if ( e_int > 0.0) then
                continue
            end if
	    end do
        
        calculate_temperature = temp
        
        
    end function
    
	recursive pure function calculate_temperature_Pconst(this,temperature,h_s,molar_fractions)

		real(dp)             :: calculate_temperature_Pconst

		class(thermophysical_properties)    ,intent(in) :: this
		real(dp)                         ,intent(in) :: temperature, h_s
		real(dp) ,dimension(:)           ,intent(in) :: molar_fractions
		!!$OMP threadprivate (calculate_temperature)

		real(dp) :: temp		
		integer :: species_number
		integer :: i, specie_number
		integer	:: iter, max_iter
		real(dp)	:: eps, residual

		species_number          = size(molar_fractions)
		temp = this%clip_temperature(temperature)

		calculate_temperature_Pconst   = temp
		
		max_iter	= 100
		eps			= 1e-08
		do iter = 1, max_iter
			calculate_temperature_Pconst = calculate_temperature_Pconst + (h_s - (this%mixture_enthalpy_molar(calculate_temperature_Pconst, molar_fractions) - this%mixture_enthalpy_molar(T_ref, molar_fractions))) / (this%mixture_cp_molar(calculate_temperature_Pconst, molar_fractions))
			residual =  h_s - this%mixture_enthalpy_molar(calculate_temperature_Pconst, molar_fractions) + this%mixture_enthalpy_molar(T_ref, molar_fractions)
			if ( abs(residual) < eps) exit
        end do
        
	end function	

	subroutine change_field_units_mole_to_dimless(this,Y)
		class(thermophysical_properties)    ,intent(in) 	:: this
		type(field_vector_cons)				,intent(inout)	:: Y

		real(dp)	:: mass_summ

		integer	,dimension(3,2)		:: loop
		integer						:: species_number
		integer						:: i, j, k, dim, specie_number

		species_number	= size(Y%pr)

		do dim = 1,3
			loop(dim,1) = lbound(Y%pr(1)%cells,dim)
			loop(dim,2) = ubound(Y%pr(1)%cells,dim)
		end do

		do k = loop(3,1),loop(3,2)
		do j = loop(2,1),loop(2,2)
		do i = loop(1,1),loop(1,2)
			mass_summ = 0.0_dp
			do specie_number = 1,species_number
				mass_summ = mass_summ + Y%pr(specie_number)%cells(i,j,k) * this%molar_masses(specie_number)
			end do

			do specie_number = 1,species_number
				Y%pr(specie_number)%cells(i,j,k)	= Y%pr(specie_number)%cells(i,j,k) * this%molar_masses(specie_number) / mass_summ
			end do
		end do
		end do
		end do

	end subroutine	
	
	subroutine change_cell_units_mole_to_dimless(this,X)
		class(thermophysical_properties)    ,intent(in) 	:: this
		real(dp)	,dimension(:)			,intent(inout)	:: X

		real(dp)	:: mass_summ

		integer						:: species_number
		integer						:: specie_number

		species_number	= size(X)

		mass_summ = 0.0_dp
		do specie_number = 1,species_number
			mass_summ = mass_summ + X(specie_number) * this%molar_masses(specie_number)
		end do

		do specie_number = 1,species_number
			X(specie_number)	= X(specie_number) * this%molar_masses(specie_number) / mass_summ
		end do

    end subroutine
    
    real(dp) function get_specie_molar_fraction(this,Y,specie_index)
		class(thermophysical_properties)    ,intent(in) 	:: this
		real(dp)	,dimension(:)			,intent(in)		:: Y
        integer			                    ,intent(in)		:: specie_index

		real(dp)	:: mass_summ

		integer						:: species_number
		integer						:: specie_number

		species_number	= size(Y)

		mass_summ = 0.0_dp
		do specie_number = 1,species_number
			mass_summ = mass_summ + Y(specie_number) / this%molar_masses(specie_number)
		end do

		get_specie_molar_fraction =	Y(specie_index) / this%molar_masses(specie_index) / mass_summ
    end function
    
    pure real(dp) function calculate_omega(this,T_red,l,s) result(omega)
		class(thermophysical_properties)    ,intent(in) 	:: this
		real(dp)                            ,intent(in)		:: T_red
        integer			                    ,intent(in)		:: l,s

		integer						:: species_number
		integer						:: specie_number

		if (T_red < 70.0_dp) then
			omega   =	this%omega_c(l,s,1) / (T_red ** this%omega_c(l,s,2))    +   &
						this%omega_c(l,s,3) / exp(this%omega_c(l,s,4) * T_red)  +   &
						this%omega_c(l,s,5) / exp(this%omega_c(l,s,6) * T_red)  +   &
						this%omega_c(l,s,7) / exp(this%omega_c(l,s,8) * T_red)
		else
			omega   =	this%omega_c(l,s,1) / (T_red ** this%omega_c(l,s,2))
		end if
  
    end function
    
    pure real(dp) function clip_temperature(this,T) result(T_clipped)
        class(thermophysical_properties)    ,intent(in) 	:: this
		real(dp)                            ,intent(in)		:: T
        
        T_clipped = min(T,tables_temperature_ceiling)
    end function



    real(dp) function mixture_molar_mass_from_mass_fractions(this, Y_mass) result(m_mix)
        class(thermophysical_properties), intent(in) :: this
        real(dp), dimension(:), intent(in) :: Y_mass
        integer :: spec, species_number
        real(dp) :: denom

        species_number = size(Y_mass)
        denom = 0.0_dp
        do spec = 1, species_number
            denom = denom + Y_mass(spec)/this%molar_masses(spec)
        end do
        m_mix = 1.0_dp/denom
    end function mixture_molar_mass_from_mass_fractions


    subroutine mole_fractions_from_mass_fractions(this, Y_mass, X_mole)
        class(thermophysical_properties), intent(in) :: this
        real(dp), dimension(:), intent(in) :: Y_mass
        real(dp), dimension(:), intent(out) :: X_mole
        integer :: spec, species_number
        real(dp) :: m_mix

        species_number = size(Y_mass)
        if (size(X_mole) /= species_number) error stop 'mole_fractions_from_mass_fractions: size mismatch'

        m_mix = this%mixture_molar_mass_from_mass_fractions(Y_mass)
        do spec = 1, species_number
            X_mole(spec) = Y_mass(spec)*m_mix/this%molar_masses(spec)
        end do
    end subroutine mole_fractions_from_mass_fractions


    real(dp) function temperature_from_pressure_density_Y(this, pressure, density, Y_mass) result(temperature)
        class(thermophysical_properties), intent(in) :: this
        real(dp), intent(in) :: pressure, density
        real(dp), dimension(:), intent(in) :: Y_mass
        real(dp) :: m_mix

        m_mix = this%mixture_molar_mass_from_mass_fractions(Y_mass)
        temperature = pressure*m_mix/(density*r_gase_J)
    end function temperature_from_pressure_density_Y


    real(dp) function mixture_specific_entropy(this, temperature, pressure, Y_mass, pressure_ref) result(s_mass)
        class(thermophysical_properties), intent(in) :: this
        real(dp), intent(in) :: temperature, pressure, pressure_ref
        real(dp), dimension(:), intent(in) :: Y_mass
        integer :: spec, species_number
        real(dp) :: m_mix, x_spec, s_molar

        species_number = size(Y_mass)
        m_mix = this%mixture_molar_mass_from_mass_fractions(Y_mass)

        s_molar = -r_gase_J*log(pressure/pressure_ref)
        do spec = 1, species_number
            x_spec = Y_mass(spec)*m_mix/this%molar_masses(spec)
            if (x_spec < 0.0_dp) error stop 'mixture_specific_entropy: negative mole fraction'
            if (x_spec > 0.0_dp) then
                s_molar = s_molar + x_spec*this%specie_entropy_molar(temperature, spec) - &
                          r_gase_J*x_spec*log(x_spec)
            end if
        end do

        s_mass = s_molar/m_mix
    end function mixture_specific_entropy


    real(dp) function temperature_from_entropy_pressure_Y(this, s_target, pressure, Y_mass, T_guess, &
            T_bracket_low, T_bracket_high, pressure_ref) result(temperature)
        class(thermophysical_properties), intent(in) :: this
        real(dp), intent(in) :: s_target, pressure, T_guess
        real(dp), dimension(:), intent(in) :: Y_mass
        real(dp), intent(in) :: T_bracket_low, T_bracket_high, pressure_ref

        integer :: iter, spec, species_number
        real(dp) :: T_low, T_high, T_loc, s_low, s_high, s_cur, residual
        real(dp) :: m_mix, x_spec, cp_molar, ds_dT, dT

        species_number = size(Y_mass)
        m_mix = this%mixture_molar_mass_from_mass_fractions(Y_mass)

        T_low = T_bracket_low
        T_high = T_bracket_high
        s_low = this%mixture_specific_entropy(T_low, pressure, Y_mass, pressure_ref)
        s_high = this%mixture_specific_entropy(T_high, pressure, Y_mass, pressure_ref)

        if (s_target > s_high) then
            write(*,*) 'Thermophysical entropy inversion failed:'
            write(*,*) '  s_target, s_low, s_high = ', s_target, s_low, s_high
            write(*,*) '  pressure, T_guess = ', pressure, T_guess
            write(*,*) '  T_low, T_high = ', T_low, T_high
            print *, 'temperature_from_entropy_pressure_Y: entropy outside temperature bracket'
            temperature = T_high
            return
        end if
        
        if (s_target < s_low) then
            write(*,*) 'Thermophysical entropy inversion failed:'
            write(*,*) '  s_target, s_low, s_high = ', s_target, s_low, s_high
            write(*,*) '  pressure, T_guess = ', pressure, T_guess
            write(*,*) '  T_low, T_high = ', T_low, T_high
            print *, 'temperature_from_entropy_pressure_Y: entropy outside temperature bracket'
            temperature = T_low
            return
        end if
        
        if (T_guess < T_low .or. T_guess > T_high) then
            write(*,*) 'Thermophysical entropy inversion failed:'
            write(*,*) '  T_guess, T_low, T_high = ', T_guess, T_low, T_high
            error stop 'temperature_from_entropy_pressure_Y: initial temperature outside bracket'
        end if

        T_loc = T_guess
        do iter = 1, 60
            s_cur = this%mixture_specific_entropy(T_loc, pressure, Y_mass, pressure_ref)
            residual = s_cur - s_target
            if (abs(residual) <= 1.0e-10_dp*max(1.0_dp, abs(s_target))) exit

            if (residual > 0.0_dp) then
                T_high = T_loc
            else
                T_low = T_loc
            end if

            cp_molar = 0.0_dp
            do spec = 1, species_number
                x_spec = Y_mass(spec)*m_mix/this%molar_masses(spec)
                if (x_spec < 0.0_dp) error stop 'temperature_from_entropy_pressure_Y: negative mole fraction'
                cp_molar = cp_molar + x_spec*this%specie_cp_molar(T_loc, spec)
            end do
            ds_dT = cp_molar/(m_mix*T_loc)
            dT = -residual/ds_dT

            if ((T_loc + dT <= T_low) .or. (T_loc + dT >= T_high) .or. (abs(dT) > 0.5_dp*T_loc)) then
                T_loc = 0.5_dp*(T_low + T_high)
            else
                T_loc = T_loc + dT
            end if
        end do

        temperature = T_loc
    end function temperature_from_entropy_pressure_Y


    real(dp) function specific_internal_energy_from_temperature_Y(this, temperature, Y_mass) result(e_mass)
        class(thermophysical_properties), intent(in) :: this
        real(dp), intent(in) :: temperature
        real(dp), dimension(:), intent(in) :: Y_mass
        integer :: spec, species_number

        species_number = size(Y_mass)
        e_mass = 0.0_dp
        do spec = 1, species_number
            e_mass = e_mass + Y_mass(spec)*this%specie_internal_energy_molar(temperature, spec)/this%molar_masses(spec)
        end do
    end function specific_internal_energy_from_temperature_Y


    subroutine pressure_derivatives_density_energy_Y(this, density, pressure, Y_mass, dp_deps, dp_drho, sound_speed)
        class(thermophysical_properties), intent(in) :: this
        real(dp), intent(in) :: density, pressure
        real(dp), dimension(:), intent(in) :: Y_mass
        real(dp), intent(out) :: dp_deps, dp_drho
        real(dp), intent(out), optional :: sound_speed

        integer :: spec, species_number
        real(dp) :: temperature, m_mix, cp_molar, cv_molar, x_spec, c2

        species_number = size(Y_mass)
        m_mix = this%mixture_molar_mass_from_mass_fractions(Y_mass)
        temperature = pressure*m_mix/(density*r_gase_J)

        cp_molar = 0.0_dp
        do spec = 1, species_number
            x_spec = Y_mass(spec)*m_mix/this%molar_masses(spec)
            cp_molar = cp_molar + x_spec*this%specie_cp_molar(temperature, spec)
        end do
        cv_molar = cp_molar - r_gase_J

        ! For a thermally-perfect ideal-gas mixture at frozen composition:
        !   p = rho * R_u/M(Y) * T,
        !   eps = eps(T,Y), d eps/dT|Y = cv_molar/M(Y).
        ! Therefore dp/deps|rho,Y = rho*R_u/cv_molar and
        ! dp/drho|eps,Y = R_u*T/M(Y) = p/rho.
        dp_deps = density*r_gase_J/cv_molar
        dp_drho = pressure/density

        if (present(sound_speed)) then
            c2 = dp_drho + dp_deps*pressure/(density*density)
            sound_speed = sqrt(c2)
        end if
    end subroutine pressure_derivatives_density_energy_Y


    subroutine pressure_derivatives_density_energy_composition_Y(this, density, pressure, Y_mass, &
            dp_deps, dp_drho, dp_dY, sound_speed)
        class(thermophysical_properties), intent(in) :: this
        real(dp), intent(in) :: density, pressure
        real(dp), dimension(:), intent(in) :: Y_mass
        real(dp), intent(out) :: dp_deps, dp_drho
        real(dp), dimension(:), intent(out) :: dp_dY
        real(dp), intent(out), optional :: sound_speed

        integer :: spec, species_number
        real(dp) :: temperature, m_mix, cp_molar, cv_molar, x_spec, c2
        real(dp) :: e_mass_spec

        species_number = size(Y_mass)
        if (size(dp_dY) /= species_number) error stop 'pressure_derivatives_density_energy_composition_Y: dp_dY size mismatch'

        m_mix = this%mixture_molar_mass_from_mass_fractions(Y_mass)
        temperature = pressure*m_mix/(density*r_gase_J)

        cp_molar = 0.0_dp
        do spec = 1, species_number
            x_spec = Y_mass(spec)*m_mix/this%molar_masses(spec)
            cp_molar = cp_molar + x_spec*this%specie_cp_molar(temperature, spec)
        end do
        cv_molar = cp_molar - r_gase_J

        ! Thermally-perfect ideal-gas mixture:
        !   p = rho * R_u * T * sum_k(Y_k/M_k),
        !   eps = sum_k Y_k e_k(T).
        ! Frozen-composition derivatives used by the Euler acoustic Jacobian:
        !   dp/deps|rho,Y = rho*R_u/cv_molar,
        !   dp/drho|eps,Y = p/rho.
        dp_deps = density*r_gase_J/cv_molar
        dp_drho = pressure/density

        ! Composition derivatives at fixed rho and eps, treating Y_k as the
        ! independent composition coordinates:
        !   dp/dY_k|rho,eps,Y_{j/=k}
        !     = rho*[R_u*T/M_k - (R_u/cv_molar)*e_k^mass(T)].
        ! For physical source rates satisfying sum_k dY_k/dt = 0, only the
        ! contraction sum_k(dp/dY_k*dY_k/dt) is used, so the usual mass-fraction
        ! constraint does not introduce an ambiguity in dp/dt.
        do spec = 1, species_number
            e_mass_spec = this%specie_internal_energy_molar(temperature, spec)/this%molar_masses(spec)
            dp_dY(spec) = density*(r_gase_J*temperature/this%molar_masses(spec) - &
                         (r_gase_J/cv_molar)*e_mass_spec)
        end do

        if (present(sound_speed)) then
            c2 = dp_drho + dp_deps*pressure/(density*density)
            sound_speed = sqrt(c2)
        end if
    end subroutine pressure_derivatives_density_energy_composition_Y


    subroutine pressure_source_rate_from_conservative_sources(this, density, pressure, E_total, velocity, Y_mass, &
            S_rho, S_momentum, S_rhoE, S_rhoY, dp_dt, de_dt, dY_dt)
        class(thermophysical_properties), intent(in) :: this
        real(dp), intent(in) :: density, pressure, E_total, S_rho, S_rhoE
        real(dp), dimension(:), intent(in) :: velocity, Y_mass, S_momentum, S_rhoY
        real(dp), intent(out) :: dp_dt, de_dt
        real(dp), dimension(:), intent(out) :: dY_dt

        integer :: spec, species_number
        real(dp) :: u2, dp_deps, dp_drho
        real(dp), dimension(size(Y_mass)) :: dp_dY

        species_number = size(Y_mass)
        if (size(dY_dt) /= species_number) error stop 'pressure_source_rate_from_conservative_sources: dY size mismatch'

        u2 = sum(velocity(:)*velocity(:))

        ! Conservative source -> primitive internal-energy and composition rates.
        ! E_total is specific total energy E = eps + |u|^2/2.
        de_dt = (S_rhoE + (u2 - E_total)*S_rho - sum(velocity(:)*S_momentum(:))) / density

        do spec = 1, species_number
            dY_dt(spec) = (S_rhoY(spec) - Y_mass(spec)*S_rho)/density
        end do

        ! Exact JANAF/thermally-perfect ideal-gas mixture pressure source rate:
        !   dp/dt = p_rho*S_rho + p_eps*deps/dt + sum_k p_Yk*dY_k/dt.
        ! The composition term is essential for chemistry and diffusion; it
        ! vanishes only for single-component/frozen-composition source terms.
        call this%pressure_derivatives_density_energy_composition_Y( &
            density, pressure, Y_mass, dp_deps, dp_drho, dp_dY)

        dp_dt = dp_drho*S_rho + dp_deps*de_dt + sum(dp_dY(:)*dY_dt(:))
    end subroutine pressure_source_rate_from_conservative_sources

!# Thermo I/O functions
	
	function get_a_coeffs(specie_name,thermo_data_file_unit) result(a_coeffs)
		character(len=10)   ,intent(in)     :: specie_name
		integer             ,intent(in)     :: thermo_data_file_unit

		real(dp)         ,dimension(7,2) :: a_coeffs


		integer             :: i,i1,i2, indx, prev, beginning, string_len, string_end, spec_num, num
		character(len=100)  :: string

		integer             :: error

		a_coeffs = 0.0_dp

		do
			read(thermo_data_file_unit,*,iostat = error)    string
			if(error /= 0) exit
			if(string == 'END') exit

			if(string == specie_name) then
				num = 0

				i1 = 1
				i2 = 2

				do while (num < 16)
					read(thermo_data_file_unit,'(A)',iostat = error) string
					if (error /= 0) then
						print *, 'ERROR : Unsuccessful specific heat capacity reading', specie_name
                        pause
						stop
					end if

					prev        = -1
					beginning   = 1

					string_len  = len(trim(string))

					do i = 1,string_len

						indx = index('-+E0123456789.',string(i:i))

						if (((indx == 0).and.(prev > 0)).or.(indx == 1).and.((prev /= 3).and.(prev > 0)).or.(indx /= 0).and.(i == string_len)) then
						select case(num)
							case(5,11)
								num = num + 1
							case(16)
								exit
							case default
								if(i /= string_len) then
									num = num + 1
									read(string(beginning:i-1),*,iostat = error) a_coeffs(i1,i2)
									continue
								else
									num = num + 1
									read(string(beginning:i),*,iostat = error) a_coeffs(i1,i2)
									continue
								end if

								if (error /= 0) then
									print *, 'ERROR: Unsuccessful CHAR conversion while reading specific heat capacity for specie :', specie_name
                                    pause
									stop
								end if

								i1 = i1 + 1
								if(num == 8) then
									i1 = 1
									i2 = 1
								end if

						end select

						beginning = i

						end if

						prev = indx
					end do
				end do
			end if
		end do

		rewind(thermo_data_file_unit)

    end function get_a_coeffs
    
	function get_specie_molar_mass(specie_name,molar_masses_data_file_unit)   result(molar_mass)
		character(len=10)   ,intent(in)     :: specie_name
		integer             ,intent(in)     :: molar_masses_data_file_unit

		character(len=100)  :: string

		real(dp) :: molar_mass

		real(dp) :: molar_mass_in_file
		integer     :: error

		molar_mass = 0.0_dp

		do
			read(molar_masses_data_file_unit,*,iostat = error) string, molar_mass_in_file
			if(error /= 0) then 
                print *, 'molar mass is not defined for specie:', specie_name
                pause
                stop
            end if
			if(string == specie_name) then
				molar_mass = molar_mass_in_file
                exit
			end if
		end do

		rewind(molar_masses_data_file_unit)

	end function get_specie_molar_mass

	subroutine get_transport_properties(potential_well_depth, collision_diameter, specie_name, transport_data_file_unit)
		real(dp)         ,intent(inout)  :: potential_well_depth, collision_diameter
		character(len=10)   ,intent(in)     :: specie_name
		integer             ,intent(in)     :: transport_data_file_unit

		character(len=100)  :: string
		real(dp) :: potential_well_depth_in_file, collision_diameter_in_file

		integer     :: dummy
		integer     :: error

		potential_well_depth    = 0.0_dp
		collision_diameter      = 0.0_dp

		read(transport_data_file_unit,*,iostat=error) string

		do
			read(transport_data_file_unit,*,iostat=error) string, dummy, potential_well_depth_in_file, collision_diameter_in_file
			if(error /= 0) then
				print *, 'Error: unsuccessufull transport properties reading for specie:', specie_name
                pause
				stop
			end if
			
			if(string == specie_name) then
				collision_diameter = collision_diameter_in_file
				potential_well_depth = potential_well_depth_in_file
				exit
			end if
		end do

		collision_diameter = collision_diameter * 0.1_dp		! Angstrem -> nanometer

		rewind(transport_data_file_unit)

	end subroutine get_transport_properties    
    
    
end module
