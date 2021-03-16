module solver_options_class

	use kind_parameters
	use global_data

	implicit none

	private
	public  solid_particles_phase
	public  liquid_droplets_phase
	public  solver_options ,solver_options_c
	
	type	:: solid_particles_phase
		real(dkind)	:: diameter, material_heat_capacity, material_density
	end type
	
	type	:: liquid_droplets_phase
		real(dkind)			:: diameter, material_heat_capacity, material_density, material_latent_heat, material_boiling_temperature
		character(len=20)	:: material
		logical				:: combustible
	end type
	
	type    :: solver_options
		private
		character(len=20)	:: solver_name	
		logical				:: hydrodynamics_flag, heat_transfer_flag, molecular_diffusion_flag, viscosity_flag, chemical_reaction_flag, CFL_flag
		real(dkind)			:: CFL_coefficient, initial_time_step
		integer				:: additional_particles_phases
		integer				:: additional_droplets_phases
		type(solid_particles_phase)	,dimension(:)	,allocatable	:: particles
		type(liquid_droplets_phase)	,dimension(:)	,allocatable	:: droplets
		integer				:: particles_phase_counter
		integer				:: droplets_phase_counter
	contains
		procedure	,private	:: read_properties
		procedure	,private	:: write_properties
		procedure	,private	:: set_properties

		procedure	:: create_additional_phase
		
		procedure	:: get_solver_name
		procedure	:: get_CFL_condition_coefficient
		procedure	:: get_initial_time_step
		procedure	:: get_CFL_condition_flag
		procedure	:: get_hydrodynamics_flag
		procedure	:: get_heat_transfer_flag
		procedure	:: get_molecular_diffusion_flag
		procedure	:: get_viscosity_flag
		procedure	:: get_chemical_reaction_flag
		procedure	:: get_additional_particles_phases_number
		procedure	:: get_additional_droplets_phases_number
		procedure	:: get_particles_params
		procedure	:: get_droplets_params
        
		procedure	:: write_log		
	end type
	
	interface   solver_options_c 
		module procedure   constructor
		module procedure   constructor_file
	end interface

contains

	type(solver_options) function constructor(solver_name,hydrodynamics_flag, heat_transfer_flag, molecular_diffusion_flag, viscosity_flag, chemical_reaction_flag, additional_particles_phases, additional_droplets_phases, CFL_flag, CFL_coefficient, initial_time_step)
		character(len=*)	,intent(in)	:: solver_name
		logical				,intent(in)	:: hydrodynamics_flag, heat_transfer_flag, molecular_diffusion_flag, viscosity_flag, chemical_reaction_flag, CFL_flag
		integer				,intent(in)	,optional	:: additional_particles_phases
		integer				,intent(in)	,optional	:: additional_droplets_phases
		real(dkind)			,intent(in)	:: CFL_coefficient, initial_time_step
		
		integer	:: additional_particles, additional_droplets
		
		integer	:: io_unit
		
		if(present(additional_particles_phases))	then
			additional_particles = additional_particles_phases
		else
			additional_particles = 0
		end if
		
		if(present(additional_droplets_phases))	then
			additional_droplets	= additional_droplets_phases
		else
			additional_droplets = 0
		end if
		
		call constructor%set_properties(solver_name,hydrodynamics_flag, heat_transfer_flag, molecular_diffusion_flag, viscosity_flag, chemical_reaction_flag, additional_particles, additional_droplets, CFL_flag, CFL_coefficient, initial_time_step)
		
		open(newunit = io_unit, file = solver_data_file_name, status = 'replace', form = 'formatted')
		call constructor%write_properties(io_unit)
		close(io_unit)	
	end function

	type(solver_options) function constructor_file()

		integer	:: io_unit
	
		real(dkind)	:: particles_diameter, particles_material_heat_capacity, particles_material_density
		real(dkind)	:: droplets_diameter, droplets_material_heat_capacity, droplets_material_density,droplets_material_latent_heat, droplets_material_boiling_temperature
		character(len=20)	:: droplets_material
		logical		:: droplets_combustible
		
		integer	:: particles_phase_counter, droplets_phase_counter

		namelist /particles_phase/	particles_diameter, particles_material_heat_capacity, particles_material_density
		namelist /droplets_phase/	droplets_diameter, droplets_material_heat_capacity, droplets_material_density,droplets_material_latent_heat, droplets_material_boiling_temperature, droplets_material, droplets_combustible
		
		open(newunit = io_unit, file = solver_data_file_name, status = 'old', form = 'formatted')
		call constructor_file%read_properties(io_unit)
		
		do particles_phase_counter = 1, size(constructor_file%particles)
			read(unit = io_unit, nml = particles_phase)
			constructor_file%particles(particles_phase_counter)%diameter					= particles_diameter
			constructor_file%particles(particles_phase_counter)%material_heat_capacity		= particles_material_heat_capacity
			constructor_file%particles(particles_phase_counter)%material_density			= particles_material_density
		end do
		
		rewind(io_unit)
		
		do droplets_phase_counter = 1, size(constructor_file%droplets)
			read(unit = io_unit, nml = droplets_phase)
            constructor_file%droplets(droplets_phase_counter)%diameter					= droplets_diameter
			constructor_file%droplets(droplets_phase_counter)%material					= droplets_material
			constructor_file%droplets(droplets_phase_counter)%material_heat_capacity	= droplets_material_heat_capacity
			constructor_file%droplets(droplets_phase_counter)%material_density			= droplets_material_density
            constructor_file%droplets(droplets_phase_counter)%material_latent_heat		= droplets_material_latent_heat
			constructor_file%droplets(droplets_phase_counter)%material_boiling_temperature		= droplets_material_boiling_temperature
			constructor_file%droplets(droplets_phase_counter)%combustible				= droplets_combustible
			
		end do		

		close(io_unit)
	end function	
	
	subroutine write_properties(this,solver_data_unit)
		class(solver_options)	,intent(in)	:: this
		integer					,intent(in)	:: solver_data_unit
		
		logical				:: hydrodynamics_flag, heat_transfer_flag, molecular_diffusion_flag, viscosity_flag, chemical_reaction_flag, CFL_flag
		integer				:: additional_particles_phases, additional_droplets_phases
		real(dkind)			:: CFL_coefficient, initial_time_step
		character(len=20)	:: solver_name
		
		namelist /solver_properties/  solver_name, hydrodynamics_flag, heat_transfer_flag, molecular_diffusion_flag, viscosity_flag, chemical_reaction_flag, additional_particles_phases, additional_droplets_phases, CFL_flag, CFL_coefficient, initial_time_step
		
		solver_name					=	this%solver_name
		hydrodynamics_flag			=	this%hydrodynamics_flag
		heat_transfer_flag			=	this%heat_transfer_flag
		molecular_diffusion_flag	=	this%molecular_diffusion_flag
		viscosity_flag				=	this%viscosity_flag
		chemical_reaction_flag		=	this%chemical_reaction_flag
		additional_particles_phases		=	this%additional_particles_phases	
		additional_droplets_phases		=	this%additional_droplets_phases
		CFL_flag					=	this%CFL_flag
		CFL_coefficient				=	this%CFL_coefficient
		initial_time_step			=	this%initial_time_step
		
		write(unit = solver_data_unit, nml = solver_properties)

	end subroutine
	
	subroutine read_properties(this,solver_data_unit)
		class(solver_options)	,intent(inout)	:: this
		integer					,intent(in)		:: solver_data_unit
		
		logical				:: hydrodynamics_flag, heat_transfer_flag, molecular_diffusion_flag, viscosity_flag, chemical_reaction_flag, CFL_flag
		integer				:: additional_particles_phases, additional_droplets_phases
		real(dkind)			:: CFL_coefficient, initial_time_step
		character(len=20)	:: solver_name
	
		namelist /solver_properties/  solver_name, hydrodynamics_flag, heat_transfer_flag, molecular_diffusion_flag, viscosity_flag, chemical_reaction_flag, additional_particles_phases, additional_droplets_phases, CFL_flag, CFL_coefficient, initial_time_step
		
		read(unit = solver_data_unit, nml = solver_properties)
		call this%set_properties(solver_name,hydrodynamics_flag, heat_transfer_flag, molecular_diffusion_flag, viscosity_flag, chemical_reaction_flag, additional_particles_phases, additional_droplets_phases, CFL_flag, CFL_coefficient, initial_time_step)
		
	end subroutine
	
	subroutine set_properties(this, solver_name,hydrodynamics_flag, heat_transfer_flag, molecular_diffusion_flag, viscosity_flag, chemical_reaction_flag, additional_particles_phases, additional_droplets_phases, CFL_flag, CFL_coefficient, initial_time_step)
		class(solver_options)	,intent(inout)	:: this
		logical								,intent(in)	:: hydrodynamics_flag, heat_transfer_flag, molecular_diffusion_flag, viscosity_flag, chemical_reaction_flag, CFL_flag
		integer								,intent(in)	:: additional_particles_phases, additional_droplets_phases
		real(dkind)							,intent(in)	:: CFL_coefficient, initial_time_step
		character(len=*)					,intent(in)	:: solver_name

		this%solver_name				=	solver_name
		this%hydrodynamics_flag			=	hydrodynamics_flag
		this%heat_transfer_flag			=	heat_transfer_flag
		this%molecular_diffusion_flag	=	molecular_diffusion_flag
		this%viscosity_flag				=	viscosity_flag
		this%chemical_reaction_flag		=	chemical_reaction_flag
		allocate(this%particles(additional_particles_phases))
		allocate(this%droplets(additional_droplets_phases))
		this%additional_particles_phases	= additional_particles_phases
		this%additional_droplets_phases		= additional_droplets_phases
		this%CFL_flag					=	CFL_flag
		this%CFL_coefficient			=	CFL_coefficient
		this%initial_time_step			=	initial_time_step

		this%particles_phase_counter		=	0
		this%droplets_phase_counter			=	0
		
	end subroutine
	
	subroutine write_log(this,log_unit)
		class(solver_options)	,intent(in)	:: this	
		integer					,intent(in)	:: log_unit	
		
		write(log_unit,'(A)')		'************************************************************************************* '
		write(log_unit,'(A)')		' Solver setup : '
		write(log_unit,'(A,A)')		' Solver name			: ',	this%solver_name
		write(log_unit,'(A,L)')		' Hydrodynamics			: ',	this%hydrodynamics_flag
		write(log_unit,'(A,L)')     ' Heat transfer			: ',	this%heat_transfer_flag
		write(log_unit,'(A,L)')     ' Molecular diffusion	: ',	this%molecular_diffusion_flag
		write(log_unit,'(A,L)')     ' Viscosity				: ',	this%viscosity_flag
		write(log_unit,'(A,L)')     ' Chemical reaction		: ',	this%chemical_reaction_flag
		write(log_unit,'(A,L)')     ' Solid particles		: ',	this%additional_particles_phases
		write(log_unit,'(A,L)')     ' Liquid droplets		: ',	this%additional_droplets_phases
		write(log_unit,'(A,L)')     ' Courant-Freidrichs-Lewy condition		: ',	this%CFL_flag
		write(log_unit,'(A,E14.7)')     ' Courant-Freidrichs-Lewy coefficient	: ',	this%CFL_coefficient
		write(log_unit,'(A,E14.7)')     ' Initial time step		: ',	this%initial_time_step
		write(log_unit,'(A)')		'************************************************************************************* '
		
	end subroutine
	
	pure function get_solver_name(this)
		class(solver_options)	,intent(in)	:: this
		character(len=20)					:: get_solver_name
		
		get_solver_name	= this%solver_name
	end function
		
	pure function get_hydrodynamics_flag(this)
		class(solver_options)	,intent(in)	:: this
		logical								:: get_hydrodynamics_flag
		
		get_hydrodynamics_flag	= this%hydrodynamics_flag
	end function	
	
	pure function get_heat_transfer_flag(this)
		class(solver_options)	,intent(in)	:: this
		logical								:: get_heat_transfer_flag
		
		get_heat_transfer_flag	= this%heat_transfer_flag
	end function	

	pure function get_molecular_diffusion_flag(this)
		class(solver_options)	,intent(in)	:: this
		logical								:: get_molecular_diffusion_flag
		
		get_molecular_diffusion_flag	= this%molecular_diffusion_flag
	end function	

	pure function get_viscosity_flag(this)
		class(solver_options)	,intent(in)	:: this
		logical								:: get_viscosity_flag
		
		get_viscosity_flag	= this%viscosity_flag
	end function	
	
	pure function get_chemical_reaction_flag(this)
		class(solver_options)	,intent(in)	:: this
		logical								:: get_chemical_reaction_flag
		
		get_chemical_reaction_flag	= this%chemical_reaction_flag
	end function
	
	pure function get_additional_particles_phases_number(this)
		class(solver_options)	,intent(in)	:: this
		integer								:: get_additional_particles_phases_number
		
		get_additional_particles_phases_number	= this%additional_particles_phases
	end function	

	pure function get_additional_droplets_phases_number(this)
		class(solver_options)	,intent(in)	:: this
		integer								:: get_additional_droplets_phases_number
		
		get_additional_droplets_phases_number	= this%additional_droplets_phases
	end function	
	
	pure function get_CFL_condition_flag(this)
		class(solver_options)	,intent(in)	:: this
		logical								:: get_CFL_condition_flag
		
		get_CFL_condition_flag	= this%CFL_flag
	end function	
	
	pure function get_CFL_condition_coefficient(this)
		class(solver_options)	,intent(in)	:: this
		real(dkind)							:: get_CFL_condition_coefficient
		
		get_CFL_condition_coefficient	= this%CFL_coefficient
	end function	
	
	pure function get_initial_time_step(this)
		class(solver_options)	,intent(in)	:: this
		real(dkind)							:: get_initial_time_step
		
		get_initial_time_step	= this%initial_time_step
	end function	
	
	pure function get_particles_params(this,phase_number)
		class(solver_options)	,intent(in)	:: this
		integer					,intent(in)	:: phase_number
		type(solid_particles_phase)			:: get_particles_params
		
		get_particles_params	= this%particles(phase_number)
    end function	

	pure function get_droplets_params(this,phase_number)
		class(solver_options)	,intent(in)	:: this
		integer					,intent(in)	:: phase_number
		type(liquid_droplets_phase)			:: get_droplets_params
		
		get_droplets_params	= this%droplets(phase_number)
	end function    
    
	subroutine create_additional_phase(this, phase_type, solid_particles_parameters, liquid_droplets_parameters)
		class(solver_options)	,intent(inout)	:: this
		character(len=15)		,intent(in)	:: phase_type
		
		type(solid_particles_phase)	,intent(in)	,optional	:: solid_particles_parameters
		type(liquid_droplets_phase)	,intent(in)	,optional	:: liquid_droplets_parameters

		real(dkind)	:: particles_diameter, particles_material_heat_capacity, particles_material_density
		real(dkind)	:: droplets_diameter, droplets_material_heat_capacity, droplets_material_density,droplets_material_latent_heat, droplets_material_boiling_temperature	
		character(len=20)	:: droplets_material
		logical		:: droplets_combustible
		
		integer	:: io_unit
		
		namelist /particles_phase/	particles_diameter, particles_material_heat_capacity, particles_material_density
		namelist /droplets_phase/	droplets_diameter, droplets_material_heat_capacity, droplets_material_density,droplets_material_latent_heat, droplets_material_boiling_temperature, droplets_material, droplets_combustible	
		
		select case(trim(phase_type))
			case("solid particles")
				this%particles_phase_counter = this%particles_phase_counter + 1
		
				if(this%particles_phase_counter > this%additional_particles_phases) then
					print *, 'ERROR: Trying to set too many additoinal particles phases, check additional_particles_phases argument.'
					stop
				end if			
				if (present(solid_particles_parameters)) then
					this%particles(this%particles_phase_counter) = solid_particles_parameters
					open(newunit = io_unit, file = solver_data_file_name, status = 'old', form = 'formatted', position = 'append')
					particles_diameter					= this%particles(this%particles_phase_counter)%diameter
					particles_material_heat_capacity	= this%particles(this%particles_phase_counter)%material_heat_capacity
					particles_material_density			= this%particles(this%particles_phase_counter)%material_density
					write(unit = io_unit, nml = particles_phase)
					close(io_unit)
				else
					print *, 'ERROR: Please provide parameters for solid particles phase'
				end if
			case("liquid droplets")
				this%droplets_phase_counter = this%droplets_phase_counter + 1
		
				if(this%droplets_phase_counter > this%additional_droplets_phases) then
					print *, 'ERROR: Trying to set too many additoinal droplets phases, check additional_droplets_phases argument.'
					stop
				end if				
				if (present(liquid_droplets_parameters)) then
					this%droplets(this%droplets_phase_counter) = liquid_droplets_parameters
					open(newunit = io_unit, file = solver_data_file_name, status = 'old', form = 'formatted', position = 'append')
                    droplets_diameter                   = this%droplets(this%droplets_phase_counter)%diameter
					droplets_material                   = this%droplets(this%droplets_phase_counter)%material
					droplets_material_heat_capacity		= this%droplets(this%droplets_phase_counter)%material_heat_capacity
					droplets_material_density			= this%droplets(this%droplets_phase_counter)%material_density
                    droplets_material_latent_heat       = this%droplets(this%droplets_phase_counter)%material_latent_heat
					droplets_material_boiling_temperature	= this%droplets(this%droplets_phase_counter)%material_boiling_temperature
					droplets_combustible				= this%droplets(this%droplets_phase_counter)%combustible
					write(unit = io_unit, nml = droplets_phase)
					close(io_unit)					
				else
					print *, 'ERROR: Please provide parameters for liquid droplets phase'
				end if		
		end select

	end subroutine
	
	
	
	
end module

