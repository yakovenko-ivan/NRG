module post_processor_class

	use kind_parameters
	use global_data
	use computational_domain_class	
	use data_manager_class
	use field_pointers
	use post_processor_operation_class
    use boundary_conditions_class

	implicit none

	private
	public	:: post_processor, post_processor_c

	type post_processor
		private
		type(post_processor_operation)			,dimension(:)	,allocatable	:: post_processor_operations
		real(dkind)								,dimension(:)	,allocatable	:: values
		type(computational_domain)			:: domain
		type(boundary_conditions_pointer)	:: boundaries
		character(len=30)		:: post_processor_output_file
		real(rkind)             :: save_time
		character(len=20)       :: save_time_units
		real(rkind)             :: save_time_coefficient
		character(len=10)       :: save_time_units_abbreviation
		integer					:: operations_number
		integer					:: operation_counter
		integer					:: output_counter
	contains
		! IO Getters/Setters
		procedure	,private	:: read_properties
		procedure	,private	:: write_properties
		procedure	,private	:: set_properties	

		! Logger
		procedure	:: write_log

		procedure	:: write_output_counter
		
		procedure	:: create_post_processor_operation
		procedure	:: process_data
	end type post_processor

	interface post_processor_c
		module procedure	constructor
		module procedure	constructor_file
	end interface

contains

	type(post_processor)	function constructor(manager,post_processor_output_file,operations_number,save_time,save_time_units,post_processor_data_file_name)
		type(data_manager)	,intent(in)	:: manager
		character(len=*)	,intent(in)	:: post_processor_output_file
		integer				,intent(in)	:: operations_number
		real(dkind)			,intent(in)	:: save_time
		character(len=*)	,intent(in)	:: save_time_units
		character(len=*)	,intent(in)	:: post_processor_data_file_name

		integer	:: output_counter = 0
		integer	:: io_unit
		
		call constructor%set_properties(manager,post_processor_output_file,operations_number,save_time,save_time_units,output_counter)

		open(newunit = io_unit, file = post_processor_data_file_name, status = 'replace', form = 'formatted')
		call constructor%write_properties(io_unit)
		close(io_unit)	
	end function

	type(post_processor)	function constructor_file(manager,post_processor_data_file_name)
		type(data_manager)	,intent(in)	:: manager
		character(len=50)	,intent(in)	:: post_processor_data_file_name
		
		integer	:: io_unit
	
		integer	:: operations_counter	
		
		open(newunit = io_unit, file = post_processor_data_file_name, status = 'old', form = 'formatted')
		call constructor_file%read_properties(manager,io_unit)

		do operations_counter = 1, size(constructor_file%post_processor_operations)
			constructor_file%post_processor_operations(operations_counter) = post_processor_operation_c(manager,io_unit)
		end do
		close(io_unit)			

	end function
	
	subroutine read_properties(this,manager,pproc_unit)
		class(post_processor)	,intent(inout)	:: this
		type(data_manager)		,intent(in)		:: manager
		integer					,intent(in)		:: pproc_unit
		
		character(len=30)	:: post_processor_output_file
		integer				:: operations_number
		real(dkind)			:: save_time
		character(len=20)	:: save_time_units

		integer				:: output_counter
		integer				:: ierr
		
		namelist /pproc_parameters/ post_processor_output_file, operations_number, save_time, save_time_units
		namelist /pproc_output_counter/ output_counter
		
		read(unit = pproc_unit, nml = pproc_parameters)
		read(unit = pproc_unit, nml = pproc_output_counter, iostat = ierr)
		
		if (ierr /= 0) output_counter = 0
		rewind(pproc_unit)
		call this%set_properties(manager,post_processor_output_file,operations_number,save_time,save_time_units,output_counter)
	end subroutine
	
	subroutine write_properties(this,pproc_unit)
		class(post_processor)	,intent(in)	:: this
		integer					,intent(in)	:: pproc_unit
		
		character(len=30)	:: post_processor_output_file
		integer				:: operations_number
		real(dkind)			:: save_time
		character(len=20)	:: save_time_units
		
		namelist /pproc_parameters/ post_processor_output_file, operations_number, save_time, save_time_units

		post_processor_output_file	= this%post_processor_output_file
		operations_number			= size(this%post_processor_operations)
		save_time					= this%save_time
		save_time_units				= this%save_time_units
		
		write(unit = pproc_unit, nml = pproc_parameters)
	end subroutine
	
	subroutine set_properties(this,manager,post_processor_output_file,operations_number,save_time,save_time_units,output_counter)
		class(post_processor)	,intent(inout)	:: this
		character(len=30)		,intent(in)	:: post_processor_output_file
		type(data_manager)		,intent(in)	:: manager
		integer					,intent(in)	:: operations_number
		real(dkind)				,intent(in)	:: save_time
		character(len=*)		,intent(in)	:: save_time_units	
		integer					,intent(in)	:: output_counter
		
		this%domain				= manager%domain	
		this%boundaries%bc_ptr	=> manager%boundary_conditions_pointer%bc_ptr
		
		allocate(this%post_processor_operations(operations_number))
		allocate(this%values(operations_number + 1 + this%domain%get_domain_dimensions()))

		this%post_processor_output_file	= post_processor_output_file
		this%output_counter				= output_counter
		this%save_time					= save_time
		this%save_time_units			= save_time_units
		this%output_counter				= output_counter
				
		this%operations_number			= operations_number
		this%operation_counter			= 0
		
		select case(this%save_time_units)
			 case('milliseconds')
				  this%save_time_units_abbreviation = 'ms'
				  this%save_time_coefficient        = 1e+03_dkind
			 case('microseconds')
				  this%save_time_units_abbreviation = 'us'
				  this%save_time_coefficient        = 1e+06_dkind
			 case('nanoseconds')
				  this%save_time_units_abbreviation = 'ns'
				  this%save_time_coefficient        = 1e+09_dkind
				  
		end select	
	end subroutine

	subroutine write_log(this,log_unit)
		class(post_processor)	,intent(in)	:: this
		integer					,intent(in)	:: log_unit
		
		integer	:: operations_counter

		write(log_unit,'(A,E14.7,A)')	' Write time: ',	this%save_time, this%save_time_units_abbreviation
		write(log_unit,'(A)')			' Operations: '

		do operations_counter = 1, size(this%post_processor_operations)
			write(log_unit,'(A,I2,A)') '	#',	operations_counter, ':'
			call this%post_processor_operations(operations_counter)%write_log(log_unit)
		end do		

	end subroutine
	
	subroutine write_output_counter(this,post_processor_output_file)
		class(post_processor)	,intent(in)	:: this
		character(len=*)		,intent(in)	:: post_processor_output_file

		integer				:: io_unit
		integer				:: ierr
		integer				:: output_counter
		
		namelist /pproc_output_counter/ output_counter

		open(newunit = io_unit, file = post_processor_output_file, status = 'old', form = 'formatted')
		read(io_unit, nml = pproc_output_counter, iostat = ierr)
		
		if (ierr == 0 ) then
			backspace(io_unit)
			backspace(io_unit)
			backspace(io_unit)
		end if
		
		output_counter = this%output_counter
		
		write(unit = io_unit, nml = pproc_output_counter)
		close(io_unit)
		
	end subroutine
	
	subroutine create_post_processor_operation(this,manager,field_name,operation_type,operation_area,operation_area_distance,grad_projection,post_processor_data_file_name)
		class(post_processor)			,intent(inout)	:: this
		type(data_manager)				,intent(in)	:: manager
		character(len=*)				,intent(in)	:: field_name
		character(len=*)				,intent(in)	:: operation_type
		integer	,dimension(3,2)			,intent(in)	:: operation_area	
		integer	,dimension(3)			,intent(in)	:: operation_area_distance
		
		integer							,intent(in)	:: grad_projection
		
		character(len=50)				,intent(in)	:: post_processor_data_file_name
		
!		integer	,save	:: operation_number = 0
		integer			:: io_unit
		
		this%operation_counter = this%operation_counter + 1
		
		if(this%operation_counter > this%operations_number) then
			print *, 'ERROR: Trying to set too many operations, check operations_number argument.'
			stop
		end if
		
		open(newunit = io_unit, file = post_processor_data_file_name, status = 'old', form = 'formatted', position = 'append')
		
		this%post_processor_operations(this%operation_counter)	= post_processor_operation_c	(manager					= manager			, &
																								field_name					= field_name		, &
																								operation_type				= operation_type	, &
																								operation_area				= operation_area	, &
																								operation_area_distance		= operation_area_distance		, &
																								grad_projection				= grad_projection	, &
																								post_processor_data_file_unit = io_unit)
		close(io_unit)	
	end subroutine
	
	subroutine process_data(this,time)
		class(post_processor)		,intent(inout)		:: this
		real(dkind)					,intent(in)			:: time	
		integer						:: output_file_unit		
		integer	,dimension(3)		:: leading_point_indexes, point_indexes
		
		logical	:: file_exists
		integer	:: i, dim, proj_counter, operations_counter

		integer					:: processor_rank

		processor_rank = this%domain%get_processor_rank()

		if((time*this%save_time_coefficient >= this%save_time*this%output_counter).or.(time == 0.0)) then
			if (processor_rank == 0) then
				inquire(file = this%post_processor_output_file, exist = file_exists)
				if (file_exists) then
					open(newunit = output_file_unit, file = this%post_processor_output_file, status = 'old', form = 'formatted', position = 'append')
				else
					open(newunit = output_file_unit, file = this%post_processor_output_file, status = 'new', form = 'formatted', position = 'append')
				end if
			end if

			this%values(1) = time 
			
			call this%post_processor_operations(1)%process_operation(this%boundaries,leading_point_indexes,this%values(2))
				
			do dim = 1,this%domain%get_domain_dimensions()
				this%values(2 + dim) = leading_point_indexes(dim)
			end do			
				
			do operations_counter = 2, size(this%post_processor_operations)
				call this%post_processor_operations(operations_counter)%set_point(leading_point_indexes)
				if(associated(this%post_processor_operations(operations_counter)%operation_field_scal%s_ptr)) then
					call this%post_processor_operations(operations_counter)%process_operation(this%boundaries,point_indexes,this%values(operations_counter+this%domain%get_domain_dimensions()+1))
				end if
				if(associated(this%post_processor_operations(operations_counter)%operation_field_vect%v_ptr)) then
					call this%post_processor_operations(operations_counter)%process_operation(this%boundaries,point_indexes,this%values(operations_counter+this%domain%get_domain_dimensions()+1))
				end if
			end do

			if (processor_rank == 0) then	
				write(output_file_unit,'(*(E14.6))')  (this%values(i), i=1,size(this%values))
				close(output_file_unit)
			end if
			
			this%output_counter = this%output_counter + 1

		end if
		
	end subroutine
end module
