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
	
		character(len=30)								::	thermo_data_file_name, transport_data_file_name, molar_masses_data_file_name
		real(dkind) ,dimension(:)       ,allocatable    ::  potential_well_depth, collision_diameter, molar_masses
		real(dkind) ,dimension(:,:,:)   ,allocatable    ::  a_coeffs
	contains
		procedure	,private	:: read_properties
		procedure	,private	:: write_properties
		procedure	,private	:: set_properties

		procedure	:: calculate_mixture_cp
		procedure	:: calculate_specie_cp
		procedure	:: calculate_mixture_cp_dT
		procedure	:: calculate_specie_entropy
		procedure	:: calculate_specie_enthalpy
		procedure	:: calculate_mixture_enthalpy
		procedure	:: calculate_mixture_energy
		procedure	:: calculate_temperature
        procedure	:: calculate_temperature_Pconst

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
		
		open(newunit = io_unit, file = thermophysical_data_file_name, status = 'replace', form = 'formatted')
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
	
	recursive pure function calculate_mixture_cp(this,temperature,species_concentrations)
		class(thermophysical_properties),   intent(in)  :: this
		real(dkind)                             :: calculate_mixture_cp
		real(dkind)                 ,intent(in) :: temperature
		real(dkind) ,dimension(:)   ,intent(in) :: species_concentrations
		real(dkind)                             :: specie_cp
		!!$OMP threadprivate (calculate_mixture_cp)

		real(dkind) :: temp
		integer :: species_number
		integer :: i, specie_number

		calculate_mixture_cp = 0.0_dkind
 		species_number = size(species_concentrations)

		temp = temperature
			if(temperature >= 4200.0_rkind) temp = 4200.0_rkind

		do specie_number = 1,species_number
			if(temp < 1000.0) then
				specie_cp = this%a_coeffs(5,1,specie_number)
				do i = 4, 1, -1
					specie_cp = specie_cp * temp
					specie_cp = specie_cp + this%a_coeffs(i,1,specie_number)
				end do
			else
				specie_cp = this%a_coeffs(5,2,specie_number)
				do i = 4, 1, -1
					specie_cp = specie_cp * temp
					specie_cp = specie_cp + this%a_coeffs(i,2,specie_number)
				end do
			end if
			calculate_mixture_cp = calculate_mixture_cp + specie_cp * species_concentrations(specie_number)
		end do
	end function

	recursive pure function calculate_specie_cp(this,temperature,specie_number)
		class(thermophysical_properties),   intent(in)  :: this
		real(dkind)             :: calculate_specie_cp
		real(dkind) ,intent(in) :: temperature
		integer     ,intent(in) :: specie_number
		!!$OMP threadprivate (calculate_specie_cp)

		real(dkind) :: temp
		real(dkind) :: specie_cp

		integer :: i

		specie_cp = 0.0_dkind
		temp = temperature
        if(temperature >= 4200.0_rkind) temp = 4200.0_rkind

		if (temp < 1000.0_dkind) then
			specie_cp = this%a_coeffs(5,1,specie_number)
			do i = 4, 1, -1
				specie_cp = specie_cp * temp
				specie_cp = specie_cp + this%a_coeffs(i,1,specie_number)
			end do
		else
			specie_cp = this%a_coeffs(5,2,specie_number)
			do i = 4, 1, -1
				specie_cp = specie_cp * temp
				specie_cp = specie_cp + this%a_coeffs(i,2,specie_number)
			end do
		end if

		calculate_specie_cp = specie_cp

	end function

	recursive pure function calculate_mixture_cp_dT(this,temperature,species_concentrations)
		class(thermophysical_properties),   intent(in)  :: this
		real(dkind)                             :: calculate_mixture_cp_dT
		real(dkind)                 ,intent(in) :: temperature
		real(dkind) ,dimension(:)   ,intent(in) :: species_concentrations
		real(dkind)                             :: specie_cp_dT
		!!$OMP threadprivate (calculate_mixture_cp)

		real(dkind) :: temp
		integer :: species_number
		integer :: i, specie_number

		calculate_mixture_cp_dT = 0.0_dkind
 		species_number = size(species_concentrations)

		temp = temperature
			if(temperature >= 4200.0_rkind) temp = 4200.0_rkind

		do specie_number = 1,species_number
			if(temp < 1000.0) then
				specie_cp_dT = 4.0_dkind * this%a_coeffs(5,1,specie_number)
				do i = 4, 2, -1
					specie_cp_dT = specie_cp_dT * temp
					specie_cp_dT = specie_cp_dT + real(i-1,dkind)*this%a_coeffs(i,1,specie_number)
				end do
			else
				specie_cp_dT = 4.0_dkind * this%a_coeffs(5,2,specie_number)
				do i = 4, 2, -1
					specie_cp_dT = specie_cp_dT * temp
					specie_cp_dT = specie_cp_dT + real(i-1,dkind)*this%a_coeffs(i,2,specie_number)
				end do
			end if
			calculate_mixture_cp_dT = calculate_mixture_cp_dT + specie_cp_dT * species_concentrations(specie_number)
		end do
	end function	
	
	
	recursive pure function calculate_specie_enthalpy(this,temperature, specie_number)
		class(thermophysical_properties),   intent(in)  :: this
		real(dkind)             :: calculate_specie_enthalpy
		real(dkind) ,intent(in) :: temperature
		integer     ,intent(in) :: specie_number
		real(dkind)             :: specie_enthalpy
		!!$OMP threadprivate (calculate_specie_enthalpy)

		real(dkind) :: temp
		integer :: i

		calculate_specie_enthalpy = 0.0_dkind
		temp = temperature
		if(temperature >= 4200.0_rkind) temp = 4200.0_dkind
		
		if(temp < 1000.0) then
			specie_enthalpy = this%a_coeffs(5,1,specie_number) / 5.0_dkind
			do i = 4, 1, -1
				specie_enthalpy = specie_enthalpy * temp
				specie_enthalpy = specie_enthalpy + this%a_coeffs(i,1,specie_number) / real(i,dkind)
			end do
			specie_enthalpy = specie_enthalpy + this%a_coeffs(6,1,specie_number) / temp
		else
			specie_enthalpy = this%a_coeffs(5,2,specie_number) / 5.0_dkind 
			do i = 4, 1, -1
				specie_enthalpy = specie_enthalpy * temp
				specie_enthalpy = specie_enthalpy + this%a_coeffs(i,2,specie_number) / real(i,dkind)
			end do
			specie_enthalpy = specie_enthalpy + this%a_coeffs(6,2,specie_number) / temp
		end if
		calculate_specie_enthalpy = specie_enthalpy * temp
	end function

	recursive pure function calculate_mixture_enthalpy(this,temperature,species_concentrations)
		class(thermophysical_properties),   intent(in)  :: this
		real(dkind) ,dimension(:)   ,intent(in) :: species_concentrations
		real(dkind)             :: calculate_mixture_enthalpy
		real(dkind) ,intent(in) :: temperature

		real(dkind)             :: specie_enthalpy
		!!$OMP threadprivate (calculate_specie_enthalpy)

		integer :: species_number
		real(dkind) :: temp
		integer :: i, specie_number

		species_number = size(species_concentrations)
		
		calculate_mixture_enthalpy = 0.0_dkind
		
		temp = temperature
		if(temperature >= 4200.0_rkind) temp = 4200.0_rkind
		
		do specie_number = 1,species_number
			specie_enthalpy = 0.0_dkind

			if(temp <= 1000.0) then
				specie_enthalpy = this%a_coeffs(5,1,specie_number) / 5.0_dkind
				do i = 4, 1, -1
					specie_enthalpy = specie_enthalpy * temp
					specie_enthalpy = specie_enthalpy + this%a_coeffs(i,1,specie_number) / real(i,dkind)
				end do
				specie_enthalpy = specie_enthalpy + this%a_coeffs(6,1,specie_number) / temp
			else
				specie_enthalpy = this%a_coeffs(5,2,specie_number) / 5.0_dkind
				do i = 4, 1, -1
					specie_enthalpy = specie_enthalpy * temp
					specie_enthalpy = specie_enthalpy + this%a_coeffs(i,2,specie_number) / real(i,dkind)
				end do
				specie_enthalpy = specie_enthalpy + this%a_coeffs(6,2,specie_number) / temp
			end if
			specie_enthalpy				= specie_enthalpy * temp
			calculate_mixture_enthalpy	= calculate_mixture_enthalpy + specie_enthalpy * species_concentrations(specie_number)
		end do

	end function	
	
	recursive pure function calculate_mixture_energy(this,temperature,species_concentrations)
		class(thermophysical_properties),   intent(in)  :: this
		real(dkind) ,dimension(:)   ,intent(in) :: species_concentrations
		real(dkind)             :: calculate_mixture_energy
		real(dkind) ,intent(in) :: temperature

		real(dkind)             :: specie_energy
		!!$OMP threadprivate (calculate_specie_enthalpy)

		integer :: species_number
		real(dkind) :: temp
		integer :: i, specie_number

		species_number = size(species_concentrations)
		
		calculate_mixture_energy = 0.0_dkind
		
		temp = temperature
		if(temperature >= 4200.0_rkind) temp = 4200.0_rkind
		
		do specie_number = 1,species_number
			specie_energy = 0.0_dkind

			if(temp <= 1000.0) then
				specie_energy = this%a_coeffs(5,1,specie_number) / 5.0_dkind
				do i = 4, 1, -1
					specie_energy = specie_energy * temp
					specie_energy = specie_energy + this%a_coeffs(i,1,specie_number) / real(i,dkind)
				end do
				specie_energy = specie_energy + this%a_coeffs(6,1,specie_number) / temp
			else
				specie_energy = this%a_coeffs(5,2,specie_number) / 5.0_dkind
				do i = 4, 1, -1
					specie_energy = specie_energy * temp
					specie_energy = specie_energy + this%a_coeffs(i,2,specie_number) / real(i,dkind)
				end do
				specie_energy = specie_energy + this%a_coeffs(6,2,specie_number) / temp
			end if
			specie_energy				= (specie_energy - r_gase_J )* temp
			calculate_mixture_energy	= calculate_mixture_energy + specie_energy * species_concentrations(specie_number)
		end do

	end function	
	
	recursive pure function calculate_specie_entropy(this,temperature, specie_number)
		real(dkind)             :: calculate_specie_entropy

		class(thermophysical_properties),   intent(in)  :: this
		real(dkind) ,intent(in) :: temperature
		integer     ,intent(in) :: specie_number
		real(dkind)             :: specie_entropy
		!!$OMP threadprivate (calculate_specie_entropy)

		real(dkind) :: temp
		integer :: i

		calculate_specie_entropy = 0.0_dkind
		temp = temperature
		if(temperature >= 4200.0_rkind) temp = 4200.0_rkind

		if(temp < 1000.0) then
			specie_entropy = this%a_coeffs(5,1,specie_number) / 4.0_dkind
			do i = 4, 2, -1
				specie_entropy = specie_entropy * temp
				specie_entropy = specie_entropy + this%a_coeffs(i,1,specie_number) / real(i-1,dkind)
			end do
			specie_entropy = specie_entropy * temp + this%a_coeffs(7,1,specie_number) + this%a_coeffs(1,1,specie_number) * log(temp)
		else
			specie_entropy = this%a_coeffs(5,2,specie_number) / 4.0_dkind
			do i = 4, 2, -1
				specie_entropy = specie_entropy * temp
				specie_entropy = specie_entropy + this%a_coeffs(i,2,specie_number) / real(i-1,dkind)
			end do
			specie_entropy = specie_entropy * temp + this%a_coeffs(7,2,specie_number) + this%a_coeffs(1,2,specie_number) * log(temp)
		end if
		calculate_specie_entropy = specie_entropy
	end function

	!recursive pure function calculate_temperature(this,temperature,E_full,mol_mix_conc,cv,concentrations)
 !
	!	real(dkind)             :: calculate_temperature
 !
	!	class(thermophysical_properties)    ,intent(in) :: this
	!	real(dkind)                         ,intent(in) :: temperature, E_full, mol_mix_conc, cv
	!	real(dkind) ,dimension(:)           ,intent(in) :: concentrations
	!	!!$OMP threadprivate (calculate_temperature)
 !
	!	real(dkind) :: temp		
	!	real(dkind) :: specie_dcpdT_T, dcpdT_T
	!	integer :: species_number
	!	integer :: i, specie_number
 !
	!	calculate_temperature   = 0.0_dkind
	!	dcpdT_T                 = 0.0_dkind
 !
	!	species_number          = size(concentrations)
 !       temp = temperature
 !       if(temperature >= 4200.0_rkind) temp = 4200.0_rkind
 !
	!	do specie_number = 1,species_number
	!		if(temp < 1000.0) then
	!			specie_dcpdT_T = this%a_coeffs(5,1,specie_number) * 5.0_dkind
	!			do i = 4, 1, -1
	!				specie_dcpdT_T = specie_dcpdT_T * temp
	!				specie_dcpdT_T = specie_dcpdT_T + this%a_coeffs(i,1,specie_number) * real(i,dkind)
	!			end do
	!		else
	!			specie_dcpdT_T = this%a_coeffs(5,2,specie_number) * 5.0_dkind
	!			do i=4,1,-1
	!				specie_dcpdT_T = specie_dcpdT_T * temp
	!				specie_dcpdT_T = specie_dcpdT_T + this%a_coeffs(i,2,specie_number) * real(i,dkind)
	!			continue
	!			end do
	!		end if
	!		dcpdT_T = dcpdT_T + specie_dcpdT_T * concentrations(specie_number)
	!	end do
 !
	!	if (temperature >= 4200.0_rkind) then
 !           calculate_temperature = E_full / cv
 !       else
	!		calculate_temperature = temp + (E_full - cv * temp) / (dcpdT_T - r_gase_J * mol_mix_conc)
	!	end if
	!	continue
	!end function

	recursive pure function calculate_temperature(this,temperature,e_i,concentrations)

		real(dkind)             :: calculate_temperature

		class(thermophysical_properties)    ,intent(in) :: this
		real(dkind)                         ,intent(in) :: temperature, e_i
		real(dkind) ,dimension(:)           ,intent(in) :: concentrations
		!!$OMP threadprivate (calculate_temperature)

		real(dkind) :: temp		
		integer :: species_number
		integer :: i, specie_number
		integer	:: iter, max_iter
		real(dkind)	:: eps, residual

		species_number          = size(concentrations)
        temp = temperature
        if(temperature >= 4200.0_rkind) temp = 4200.0_rkind

		calculate_temperature   = temp
		
		max_iter	= 100
		eps			= 1e-08
		do iter = 1, max_iter
			calculate_temperature = calculate_temperature + (e_i - (this%calculate_mixture_energy(calculate_temperature, concentrations) - this%calculate_mixture_enthalpy(298.15_dkind, concentrations))) / (this%calculate_mixture_cp(calculate_temperature, concentrations) - r_gase_J)
			residual =  e_i - this%calculate_mixture_energy(calculate_temperature, concentrations) + this%calculate_mixture_enthalpy(298.15_dkind, concentrations)
			if ( abs(residual) < eps) exit
		end do

    end function
    
	recursive pure function calculate_temperature_Pconst(this,temperature,h_s,concentrations)

		real(dkind)             :: calculate_temperature_Pconst

		class(thermophysical_properties)    ,intent(in) :: this
		real(dkind)                         ,intent(in) :: temperature, h_s
		real(dkind) ,dimension(:)           ,intent(in) :: concentrations
		!!$OMP threadprivate (calculate_temperature)

		real(dkind) :: temp		
		integer :: species_number
		integer :: i, specie_number
		integer	:: iter, max_iter
		real(dkind)	:: eps, residual

		species_number          = size(concentrations)
        temp = temperature
        if(temperature >= 4200.0_rkind) temp = 4200.0_rkind

		calculate_temperature_Pconst   = temp
		
		max_iter	= 100
		eps			= 1e-08
		do iter = 1, max_iter
			calculate_temperature_Pconst = calculate_temperature_Pconst + (h_s - (this%calculate_mixture_enthalpy(calculate_temperature_Pconst, concentrations) - this%calculate_mixture_enthalpy(298.15_dkind, concentrations))) / (this%calculate_mixture_cp(calculate_temperature_Pconst, concentrations))
			residual =  h_s - this%calculate_mixture_enthalpy(calculate_temperature_Pconst, concentrations) + this%calculate_mixture_enthalpy(298.15_dkind, concentrations)
			if ( abs(residual) < eps) exit
        end do
        
	end function	
	
	function get_a_coeffs(specie_name,thermo_data_file_unit) result(a_coeffs)
		character(len=10)   ,intent(in)     :: specie_name
		integer             ,intent(in)     :: thermo_data_file_unit

		real(dkind)         ,dimension(7,2) :: a_coeffs


		integer             :: i,i1,i2, indx, prev, beginning, string_len, string_end, spec_num, num
		character(len=100)  :: string

		integer             :: error

		a_coeffs = 0.0_dkind

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

	end function

	function get_specie_molar_mass(specie_name,molar_masses_data_file_unit)   result(molar_mass)
		character(len=10)   ,intent(in)     :: specie_name
		integer             ,intent(in)     :: molar_masses_data_file_unit

		character(len=100)  :: string

		real(dkind) :: molar_mass

		real(dkind) :: molar_mass_in_file
		integer     :: error

		molar_mass = 0.0_dkind

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

	end function

	subroutine get_transport_properties(potential_well_depth, collision_diameter, specie_name, transport_data_file_unit)
		real(dkind)         ,intent(inout)  :: potential_well_depth, collision_diameter
		character(len=10)   ,intent(in)     :: specie_name
		integer             ,intent(in)     :: transport_data_file_unit

		character(len=100)  :: string
		real(dkind) :: potential_well_depth_in_file, collision_diameter_in_file

		integer     :: dummy
		integer     :: error

		potential_well_depth    = 0.0_dkind
		collision_diameter      = 0.0_dkind

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

		collision_diameter = collision_diameter * 0.1_dkind		! Angstrem -> nanometer

		rewind(transport_data_file_unit)

	end subroutine

	subroutine change_field_units_mole_to_dimless(this,Y)
		class(thermophysical_properties)    ,intent(in) 	:: this
		type(field_vector_cons)				,intent(inout)	:: Y

		real(dkind)	:: mass_summ

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
			mass_summ = 0.0_dkind
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
	
	subroutine change_cell_units_mole_to_dimless(this,Y)
		class(thermophysical_properties)    ,intent(in) 	:: this
		real(dkind)	,dimension(:)			,intent(inout)	:: Y

		real(dkind)	:: mass_summ

		integer						:: species_number
		integer						:: specie_number

		species_number	= size(Y)

		mass_summ = 0.0_dkind
		do specie_number = 1,species_number
			mass_summ = mass_summ + Y(specie_number) * this%molar_masses(specie_number)
		end do

		do specie_number = 1,species_number
			Y(specie_number)	= Y(specie_number) * this%molar_masses(specie_number) / mass_summ
		end do

	end subroutine		
end module
