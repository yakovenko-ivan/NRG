module boundary_type_class

	use kind_parameters
	
	implicit none

	private
	public	:: boundary_type, boundary_type_c

	type boundary_type
		private
		character(len=20)						:: type_name
		logical									:: slip, conductive
		real(dkind)								:: wall_temperature ,wall_conductivity_ratio
		real(dkind)								:: farfield_pressure, farfield_temperature, farfield_density, farfield_velocity, farfield_energy
		real(dkind)			,dimension(:)	,allocatable	:: farfield_concentrations
		character(len=10)	,dimension(:)	,allocatable	:: farfield_species_names			
		integer									:: priority
	contains
		! IO Getters/Setters
		procedure	,private	:: read_properties
		procedure	,private	:: write_properties
		procedure	,private	:: set_properties	

		! Getters
		procedure	:: get_type_name
		procedure	:: get_wall_temperature
		procedure	:: get_wall_conductivity_ratio
		procedure	:: get_farfield_pressure
		procedure	:: get_farfield_temperature
		procedure	:: get_farfield_density
        procedure	:: get_farfield_energy
		procedure	:: get_farfield_velocity
        procedure	:: get_farfield_species_names
		procedure	:: get_farfield_concentrations
		
		procedure	:: get_condition_priority
		procedure	:: is_slip
		procedure	:: is_conductive

		!Setters
		procedure	:: set_farfield_pressure
		procedure	:: set_farfield_temperature
		procedure	:: set_farfield_density
        procedure	:: set_farfield_energy
		procedure	:: set_farfield_concentrations

		! Logger		
		procedure	:: write_log
	end type boundary_type

	interface boundary_type_c
		module procedure constructor
		module procedure constructor_file
	end interface

contains

	type(boundary_type)  function constructor(type_name,slip,conductive,wall_temperature,wall_conductivity_ratio,farfield_pressure,farfield_temperature,farfield_density,farfield_velocity,farfield_concentrations,farfield_species_names,priority,bc_data_file_unit)
		character(len=*)		,intent(in)				:: type_name
		logical					,intent(in)	,optional	:: slip, conductive
		real(dkind)				,intent(in)	,optional	:: wall_temperature, wall_conductivity_ratio
		real(dkind)				,intent(in)	,optional	:: farfield_pressure, farfield_temperature, farfield_density, farfield_velocity
		real(dkind)			,dimension(:)	,intent(in)	,optional	:: farfield_concentrations
		character(len=*)	,dimension(:)	,intent(in)	,optional	:: farfield_species_names		
		integer					,intent(in)				:: priority
		integer					,intent(in)				:: bc_data_file_unit

		call constructor%set_properties(type_name,slip,conductive,wall_temperature,wall_conductivity_ratio,farfield_pressure,farfield_temperature,farfield_density,farfield_velocity,farfield_concentrations,farfield_species_names,priority)
		call constructor%write_properties(bc_data_file_unit)
	end function

	type(boundary_type)  function constructor_file(bc_data_file_unit)
		integer	,intent(in)	:: bc_data_file_unit
		
		call constructor_file%read_properties(bc_data_file_unit)
	end function
	
	subroutine read_properties(this,bc_data_file_unit)
		class(boundary_type)	,intent(inout)	:: this
		integer					,intent(in)	:: bc_data_file_unit
		
		character(len=20)	:: type_name
		logical				:: slip, conductive
		real(dkind)			:: wall_temperature, wall_conductivity_ratio
		real(dkind)			:: farfield_pressure, farfield_temperature, farfield_density, farfield_velocity
		integer				:: farfield_species_number
		real(dkind)			,dimension(:)	,allocatable	:: farfield_concentrations
		character(len=10)	,dimension(:)	,allocatable	:: farfield_species_names
		
		integer				:: priority
		
		namelist /boundary_type/ type_name, slip, conductive,wall_temperature, wall_conductivity_ratio ,farfield_pressure, farfield_temperature, farfield_density,farfield_velocity, farfield_species_number, priority
		namelist /boundary_type_species/farfield_concentrations, farfield_species_names

		read(unit = bc_data_file_unit, nml = boundary_type)
		
		allocate(farfield_concentrations(farfield_species_number))
		allocate(farfield_species_names(farfield_species_number))
		
		read(unit = bc_data_file_unit, nml = boundary_type_species)
		
		call this%set_properties(type_name,slip,conductive,wall_temperature,wall_conductivity_ratio,farfield_pressure,farfield_temperature,farfield_density,farfield_velocity,farfield_concentrations,farfield_species_names,priority)
		
	end subroutine
	
	subroutine write_properties(this,bc_data_file_unit)
		class(boundary_type)	,intent(in)	:: this
		integer					,intent(in)	:: bc_data_file_unit
		
		character(len=20)	:: type_name
		logical				:: slip, conductive
		real(dkind)			:: wall_temperature, wall_conductivity_ratio
		real(dkind)			:: farfield_pressure, farfield_temperature, farfield_density, farfield_velocity
		integer				:: farfield_species_number
		real(dkind)			,dimension(:)	,allocatable	:: farfield_concentrations
		character(len=10)	,dimension(:)	,allocatable	:: farfield_species_names
		
		integer				:: priority
		
		namelist /boundary_type/ type_name, slip, conductive,wall_temperature, wall_conductivity_ratio ,farfield_pressure, farfield_temperature, farfield_density,farfield_velocity, farfield_species_number, priority
		namelist /boundary_type_species/farfield_concentrations, farfield_species_names

		type_name				= this%type_name
		slip					= this%slip
		conductive				= this%conductive
		wall_temperature		= this%wall_temperature
		wall_conductivity_ratio	= this%wall_conductivity_ratio
		farfield_pressure		= this%farfield_pressure
		farfield_temperature	= this%farfield_temperature
		farfield_density		= this%farfield_density
		farfield_velocity		= this%farfield_velocity
		farfield_species_number = size(this%farfield_concentrations)
		
		allocate(farfield_concentrations(farfield_species_number))
		allocate(farfield_species_names(farfield_species_number))
		
		farfield_concentrations	= this%farfield_concentrations
		farfield_species_names	= this%farfield_species_names

		priority				= this%priority
		
		write(unit = bc_data_file_unit, nml = boundary_type)
		write(unit = bc_data_file_unit, nml = boundary_type_species)
		
	end subroutine
	
	
	subroutine set_properties(this,type_name,slip,conductive,wall_temperature,wall_conductivity_ratio,farfield_pressure,farfield_temperature,farfield_density,farfield_velocity,farfield_concentrations,farfield_species_names,priority)
		class(boundary_type)	,intent(inout)			:: this
		character(len=*)		,intent(in)				:: type_name
		logical					,intent(in)	,optional	:: slip, conductive
		real(dkind)				,intent(in)	,optional	:: wall_temperature, wall_conductivity_ratio
		real(dkind)				,intent(in)	,optional	:: farfield_pressure, farfield_temperature, farfield_density, farfield_velocity
		real(dkind)			,dimension(:)	,intent(in)	,optional	:: farfield_concentrations
		character(len=*)	,dimension(:)	,intent(in)	,optional	:: farfield_species_names	
		integer					,intent(in)				:: priority

		integer	:: specie
		
		this%type_name			= type_name

		select case(type_name)
		case('wall')
			if(present(slip).and.present(conductive)) then
				this%slip						= slip
				this%conductive					= conductive
			else
				print *, ' Wall boundary with priority ', priority,' was not correctly specified. Please provide full set of wall parameters (slip flag, conductive flag).'
				stop
			end if
			if (conductive) then
				if(present(wall_temperature).and.present(wall_conductivity_ratio)) then
					this%wall_temperature			= wall_temperature
					this%wall_conductivity_ratio	= wall_conductivity_ratio
				else
					print *, ' Wall conductive boundary with priority ', priority,' was not correctly specified. Please provide full set of conductive wall parameters (wall_temperature, wall_conductivity_ratio).'
					stop
				end if
			else
				this%wall_temperature			= 0.0_dkind
				this%wall_conductivity_ratio	= 0.0_dkind	
			end if
			this%farfield_pressure			= 0.0_dkind
			this%farfield_temperature		= 0.0_dkind
			this%farfield_density			= 0.0_dkind
			this%farfield_velocity			= 0.0_dkind
			allocate(this%farfield_concentrations(0))
			allocate(this%farfield_species_names(0))
			this%priority					= priority
		case('outlet','outlet2','inlet')
			this%slip						= .false.
			this%conductive					= .false.
			this%wall_temperature			= 0.0_dkind
			this%wall_conductivity_ratio	= 0.0_dkind

			if(present(farfield_pressure).and.present(farfield_temperature).and.present(farfield_velocity).and.present(farfield_concentrations).and.present(farfield_species_names)) then
				this%farfield_pressure			= farfield_pressure
				this%farfield_temperature		= farfield_temperature
				this%farfield_velocity			= farfield_velocity
				allocate(this%farfield_concentrations(size(farfield_concentrations)))
				this%farfield_concentrations	= farfield_concentrations
				allocate(this%farfield_species_names(size(farfield_species_names)))
				this%farfield_species_names		= farfield_species_names
			else
				print *, ' Inlet/Outlet boundary with priority ', priority,' was not correctly specified. Please provide full set of inlet parameters (P_inf,T_inf,Rho_inf,v(x)_inf,Y_inf).'
				stop
			end if
			this%farfield_density			= 0.0_dkind
			this%priority					= priority
		case('symmetry_plane')
			this%slip						= .true.
			this%conductive					= .false.
			this%wall_temperature			= 0.0_dkind
			this%wall_conductivity_ratio	= 0.0_dkind
			this%farfield_pressure			= 0.0_dkind
			this%farfield_temperature		= 0.0_dkind
			this%farfield_density			= 0.0_dkind
			this%farfield_velocity			= 0.0_dkind
			allocate(this%farfield_concentrations(0))
			allocate(this%farfield_species_names(0))			
			this%priority					= priority
		end select
	end subroutine
	
	subroutine write_log(this,log_unit)
		class(boundary_type)	,intent(in)	:: this
		integer					,intent(in)	:: log_unit
		
		integer	:: bc_types_counter

		write(log_unit,'(A,A)') ' Boundary type name:', this%type_name
		write(log_unit,'(A,L)') ' Boundary slip condition:', this%slip
		write(log_unit,'(A,L)') ' Boundary conductive condition:', this%conductive
		write(log_unit,'(A,E14.7)') ' Boundary wall temperature:', this%wall_temperature
		write(log_unit,'(A,E14.7)') ' Boundary wall conductivity ratio:', this%wall_conductivity_ratio
		write(log_unit,'(A,E14.7)') ' Boundary farfield pressure:',this%farfield_pressure
		write(log_unit,'(A,E14.7)') ' Boundary farfield temperature:',this%farfield_temperature
		write(log_unit,'(A,I2)') ' Boundary priority:', this%priority
	end subroutine
	
! ************** Getters ***************

	pure function	get_type_name(this)
		class(boundary_type)	,intent(in)	:: this
		character(len=20)					:: get_type_name

		get_type_name = this%type_name
	end function

	pure function	get_wall_temperature(this)
		class(boundary_type)	,intent(in)		:: this
		real(dkind)								:: get_wall_temperature

		get_wall_temperature = this%wall_temperature
	end function

	pure function	get_wall_conductivity_ratio(this)
		class(boundary_type)	,intent(in)		:: this
		real(dkind)								:: get_wall_conductivity_ratio

		get_wall_conductivity_ratio	= this%wall_conductivity_ratio
	end function

	pure function	get_farfield_pressure(this)
		class(boundary_type)	,intent(in)		:: this
		real(dkind)								:: get_farfield_pressure

		get_farfield_pressure = this%farfield_pressure
	end function

	pure function	get_farfield_temperature(this)
		class(boundary_type)	,intent(in)		:: this
		real(dkind)								:: get_farfield_temperature

		get_farfield_temperature = this%farfield_temperature
	end function

	pure function	get_farfield_density(this)
		class(boundary_type)	,intent(in)		:: this
		real(dkind)								:: get_farfield_density

		get_farfield_density = this%farfield_density
    end function

	pure function	get_farfield_energy(this)
		class(boundary_type)	,intent(in)		:: this
		real(dkind)								:: get_farfield_energy

		get_farfield_energy = this%farfield_energy
	end function

	pure function	get_farfield_velocity(this)
		class(boundary_type)	,intent(in)		:: this
		real(dkind)								:: get_farfield_velocity

		get_farfield_velocity = this%farfield_velocity
	end function

	subroutine	get_farfield_concentrations(this, farfield_concentrations)
		class(boundary_type)		,intent(in)		:: this
		real(dkind)	,dimension(:)	,allocatable	:: farfield_concentrations

		if (allocated(farfield_concentrations)) deallocate(farfield_concentrations)
		allocate(farfield_concentrations(size(this%farfield_concentrations)))
		
		farfield_concentrations = this%farfield_concentrations
	end subroutine

	subroutine	get_farfield_species_names(this, farfield_species_names)
		class(boundary_type)				,intent(in)		:: this
		character(len=10)	,dimension(:)	,allocatable	:: farfield_species_names

		if (allocated(farfield_species_names)) deallocate(farfield_species_names)
		allocate(farfield_species_names(size(this%farfield_species_names)))
		
		farfield_species_names = this%farfield_species_names
	end subroutine	
	
	pure integer	function	get_condition_priority(this)
		class(boundary_type)	,intent(in)	:: this

		get_condition_priority = this%priority
	end function

	pure logical function	is_slip(this)
		class(boundary_type)	,intent(in)	:: this

		is_slip = this%slip
	end function

	pure logical function	is_conductive(this)
		class(boundary_type)	,intent(in)	:: this

		is_conductive = this%conductive
	end function


	pure subroutine	set_farfield_pressure(this,farfield_pressure)
		class(boundary_type)	,intent(inout)		:: this
		real(dkind)				,intent(in)			:: farfield_pressure

		this%farfield_pressure = farfield_pressure
	end subroutine

	pure subroutine	set_farfield_temperature(this,farfield_temperature)
		class(boundary_type)	,intent(inout)		:: this
		real(dkind)				,intent(in)			:: farfield_temperature

		this%farfield_temperature = farfield_temperature
	end subroutine

	pure subroutine	set_farfield_density(this,farfield_pressure)
		class(boundary_type)	,intent(inout)		:: this
		real(dkind)				,intent(in)			:: farfield_pressure

		this%farfield_density = farfield_pressure
    end subroutine

 	pure subroutine	set_farfield_energy(this,farfield_energy)
		class(boundary_type)	,intent(inout)		:: this
		real(dkind)				,intent(in)			:: farfield_energy

		this%farfield_energy = farfield_energy
	end subroutine   
    
	pure subroutine	set_farfield_velocity(this,farfield_velocity)
		class(boundary_type)	,intent(inout)		:: this
		real(dkind)				,intent(in)			:: farfield_velocity

		this%farfield_velocity = farfield_velocity
	end subroutine	
	
	pure subroutine	set_farfield_concentrations(this,farfield_concentrations)
		class(boundary_type)		,intent(inout)		:: this
		real(dkind), dimension(:)	,intent(in)			:: farfield_concentrations

		this%farfield_concentrations = farfield_concentrations
	end subroutine		

end module
