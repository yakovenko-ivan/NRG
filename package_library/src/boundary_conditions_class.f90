module boundary_conditions_class

	use kind_parameters
	use global_data
	use computational_domain_class
	use boundary_type_class

	implicit none

	private
	public	::	boundary_conditions, boundary_conditions_pointer, boundary_conditions_c

	type boundary_conditions_pointer
		type(boundary_conditions)	,pointer	:: bc_ptr
	end type boundary_conditions_pointer

	type boundary_conditions
		type(boundary_type)		,dimension(:)		,allocatable	:: boundary_types
		integer(ikind)			,dimension(:,:,:)	,allocatable	:: bc_markers
		integer														:: default_boundary_type
		integer														:: boundary_type_counter
	contains
		! IO Getters/Setters
		procedure	,private	:: read_properties
		procedure	,private	:: write_properties
		procedure	,private	:: set_properties

		! Logger	
		procedure	:: write_log
				
		procedure	:: create_boundary_type		
	end type boundary_conditions

	interface boundary_conditions_c
		module procedure	constructor
		module procedure	constructor_file
	end interface

contains

	type(boundary_conditions)	function constructor(domain,number_of_boundary_types,default_boundary)
		type(computational_domain)	,intent(in)	:: domain
		integer						,intent(in)	:: number_of_boundary_types
		integer						,intent(in)	:: default_boundary

		integer	:: io_unit
		
		call constructor%set_properties(domain, number_of_boundary_types, default_boundary)
		
		open(newunit = io_unit, file = boundary_data_file_name, status = 'replace', form = 'formatted')
		call constructor%write_properties(io_unit)
		close(io_unit)

	end function
	
	type(boundary_conditions)	function constructor_file(domain)
		type(computational_domain)	,intent(in)	:: domain	
		
		integer	:: io_unit
		integer	:: bc_type_counter
		
		open(newunit = io_unit, file = boundary_data_file_name, status = 'old', form = 'formatted')
		call constructor_file%read_properties(domain,io_unit)

		do bc_type_counter = 1, size(constructor_file%boundary_types)
			constructor_file%boundary_types(bc_type_counter) = boundary_type_c(io_unit)
		end do		
		close(io_unit)
		
	end function
	
	subroutine read_properties(this,domain,boundary_unit)
		class(boundary_conditions)	,intent(inout)	:: this
		type(computational_domain)	,intent(in)	:: domain
		integer						,intent(in)	:: boundary_unit
		
		integer	:: number_of_boundary_types, default_boundary_type
		
		namelist /boundary_parameters/ number_of_boundary_types, default_boundary_type
		
		read(unit = boundary_unit, nml = boundary_parameters)
		
		call this%set_properties(domain, number_of_boundary_types, default_boundary_type)		
	end subroutine
	
	subroutine write_properties(this,boundary_unit)
		class(boundary_conditions)	,intent(in)	:: this
		integer						,intent(in)	:: boundary_unit
		
		integer	:: number_of_boundary_types, default_boundary_type
		
		namelist /boundary_parameters/ number_of_boundary_types, default_boundary_type
		
		number_of_boundary_types	= size(this%boundary_types)
		default_boundary_type		= this%default_boundary_type
		
		write(unit = boundary_unit, nml = boundary_parameters)
	end subroutine
	
	subroutine set_properties(this,domain,number_of_boundary_types,default_boundary)
		class(boundary_conditions)	,intent(inout)	:: this
		type(computational_domain)	,intent(in)	:: domain
		integer						,intent(in)	:: number_of_boundary_types
		integer						,intent(in)	:: default_boundary

		integer					:: dimensions
		integer	,dimension(3,2)	:: allocation_bounds

		integer	,dimension(3)	:: processor_number	
		integer	,dimension(3)	:: processor_grid_coord 

		integer	:: bound_number ,bound
		integer	:: i,j,k
		
		allocate(this%boundary_types(number_of_boundary_types))

		this%default_boundary_type = default_boundary
		
		dimensions				= domain%get_domain_dimensions()
		allocation_bounds		= domain%get_local_utter_cells_bounds()

		processor_number		= domain%get_processor_number()
		processor_grid_coord	= domain%get_processor_grid_coord()

		allocate(this%bc_markers(	allocation_bounds(1,1):allocation_bounds(1,2)	, &
									allocation_bounds(2,1):allocation_bounds(2,2)	, &
									allocation_bounds(3,1):allocation_bounds(3,2)))

		this%bc_markers = 0

		if (processor_grid_coord(1) == 0) 						this%bc_markers(allocation_bounds(1,1),:,:)	= default_boundary
		if (processor_grid_coord(1) == processor_number(1)-1) 	this%bc_markers(allocation_bounds(1,2),:,:)	= default_boundary

		if(dimensions >= 2) then
			if (processor_grid_coord(2) == 0) 						this%bc_markers(:,allocation_bounds(2,1),:)	= default_boundary
			if (processor_grid_coord(2) == processor_number(2)-1) 	this%bc_markers(:,allocation_bounds(2,2),:)	= default_boundary

			if(dimensions == 3) then
				if (processor_grid_coord(3) == 0) 						this%bc_markers(:,:,allocation_bounds(3,1))	= default_boundary
				if (processor_grid_coord(3) == processor_number(3)-1) 	this%bc_markers(:,:,allocation_bounds(3,2))	= default_boundary				
			end if
		end if

		this%boundary_type_counter = 0
		
	end subroutine
		
	subroutine write_log(this,log_unit)
		class(boundary_conditions)	,intent(in)	:: this
		integer						,intent(in)	:: log_unit
		
		integer	:: bc_types_counter

		write(log_unit,'(A)')	'************************************************************************************* '			
		write(log_unit,'(A)')	'Boundary conditions setup: '
		write(log_unit,'(A,I3)')	' Boundary types: ',	size(this%boundary_types)
		do bc_types_counter = 1, size(this%boundary_types)
			write(log_unit,'(A,I2,A)') ' Boundary type #',	bc_types_counter, ':'
			call this%boundary_types(bc_types_counter)%write_log(log_unit)
		end do	
		write(log_unit,'(A)')	'************************************************************************************* '			
	end subroutine
	
	subroutine create_boundary_type(this,type_name,slip,conductive,wall_temperature,wall_conductivity_ratio,farfield_pressure,farfield_temperature,farfield_density,farfield_velocity,farfield_concentrations,farfield_species_names,priority)
		class(boundary_conditions)	,intent(inout)			:: this
		character(len=*)			,intent(in)				:: type_name
		logical						,intent(in)	,optional	:: slip, conductive
		real(dkind)					,intent(in)	,optional	:: wall_temperature, wall_conductivity_ratio
		real(dkind)					,intent(in)	,optional	:: farfield_pressure, farfield_temperature, farfield_density, farfield_velocity
		real(dkind)			,dimension(:)	,intent(in)	,optional	:: farfield_concentrations
		character(len=*)	,dimension(:)	,intent(in)	,optional	:: farfield_species_names
		integer						,intent(in)				:: priority		
	
		integer			:: io_unit
		
		this%boundary_type_counter = this%boundary_type_counter + 1
		
		open(newunit = io_unit, file = boundary_data_file_name, status = 'old', form = 'formatted', position = 'append')

		this%boundary_types(this%boundary_type_counter) = boundary_type_c(type_name,slip,conductive,wall_temperature,wall_conductivity_ratio,farfield_pressure,farfield_temperature,farfield_density,farfield_velocity,farfield_concentrations,farfield_species_names,priority,io_unit)
		close(io_unit)
		
	end subroutine

end module
