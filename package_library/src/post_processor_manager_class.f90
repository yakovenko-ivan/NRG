module post_processor_manager_class

	use kind_parameters
	use global_data
	use computational_domain_class		
	use field_pointers
	use data_manager_class
	use post_processor_class
	use post_processor_operation_class

	implicit none

	private
	public	:: post_processor_manager, post_processor_manager_c
	
	type	:: post_processor_manager
		type(post_processor)	,dimension(:)	,allocatable	:: post_processors
		character(len=50)		,dimension(:)	,allocatable	:: post_processors_setup_file_names

		integer								:: post_processors_created
		
		type(computational_domain)			:: domain
	contains
		!IO Getters/Setters
		procedure	,private	:: read_properties
		procedure	,private	:: write_properties
		procedure	,private	:: set_properties	

		procedure	,private	:: write_output_counters
		
		! Logger
		procedure	:: write_log
		
		procedure	:: create_post_processor
		procedure	:: create_post_processor_operation
		procedure	:: process_data
	end type post_processor_manager

	interface post_processor_manager_c
		module procedure constructor
		module procedure constructor_file
	end interface
	
contains

	type(post_processor_manager)	function constructor(manager,number_post_processors) 
		type(data_manager)	,intent(in)	:: manager
		integer				,intent(in)	:: number_post_processors
		integer	:: io_unit
		
		call constructor%set_properties(manager,number_post_processors)
		
		open(newunit = io_unit, file = post_processor_manager_data_file_name, status = 'replace', form = 'formatted')
		call constructor%write_properties(io_unit)
		close(io_unit)	
	end function
	
	type(post_processor_manager)	function constructor_file(manager)
		type(data_manager)	,intent(in)	:: manager

		integer	:: io_unit
		integer	:: pproc_counter

		open(newunit = io_unit, file = post_processor_manager_data_file_name, status = 'old', form = 'formatted')
		call constructor_file%read_properties(manager,io_unit)
		close(io_unit)

		do pproc_counter = 1, size(constructor_file%post_processors)
			constructor_file%post_processors(pproc_counter) =  post_processor_c(manager, constructor_file%post_processors_setup_file_names(pproc_counter))
		end do
	end function
	
	subroutine read_properties(this,manager,pproc_manager_unit)
		class(post_processor_manager)	,intent(inout)	:: this
		type(data_manager)				,intent(in)		:: manager
		integer							,intent(in)		:: pproc_manager_unit
		
		integer	:: number_post_processors	
		
		namelist /pproc_manager_properties/ number_post_processors
		
		read(unit = pproc_manager_unit, nml = pproc_manager_properties)

		call this%set_properties(manager,number_post_processors)
	end subroutine
	
	subroutine write_properties(this,pproc_manager_unit)
		class(post_processor_manager)	,intent(in)	:: this
		integer							,intent(in)	:: pproc_manager_unit
		
		integer	:: number_post_processors	
		
		namelist /pproc_manager_properties/ number_post_processors
		
		number_post_processors = size(this%post_processors)
		write(unit = pproc_manager_unit, nml = pproc_manager_properties)
	end subroutine
	
	subroutine set_properties(this,manager,number_post_processors)
		class(post_processor_manager)	,intent(inout)	:: this
		type(data_manager)				,intent(in)		:: manager
		integer							,intent(in)		:: number_post_processors		
	
		integer	:: num_index
		integer	:: pproc_counter
		
		this%domain						= manager%domain
		this%post_processors_created	= 0
		
		allocate(this%post_processors(number_post_processors))
		allocate(this%post_processors_setup_file_names(number_post_processors))
		
		num_index = index(post_processor_data_file_name,'#')
		
		do pproc_counter = 1, number_post_processors
			write(this%post_processors_setup_file_names(pproc_counter),'(A,I2,A)') post_processor_data_file_name(1:num_index) , pproc_counter , trim(post_processor_data_file_name(num_index+1:))
		end do
		
	end subroutine

	subroutine write_log(this,log_unit)
		class(post_processor_manager)	,intent(in)	:: this
		integer							,intent(in)	:: log_unit

		integer	:: pproc		
		
		write(log_unit,'(A)')	'************************************************************************************* '			
		write(log_unit,'(A)')	' Post processors setup:'
		write(log_unit,'(A,I2)')' Amount of post processors: ', size(this%post_processors)

		do pproc = 1, size(this%post_processors)
			write(log_unit,'(A,I2)') 'Post processor #', pproc
			call this%post_processors(pproc)%write_log(log_unit)
		end do

		write(log_unit,'(A)')	'************************************************************************************* '		
		
	end subroutine
	
	subroutine write_output_counters(this)
		class(post_processor_manager)	,intent(in)	:: this

		integer	:: pproc_counter

		do pproc_counter = 1, size(this%post_processors)
			call this%post_processors(pproc_counter)%write_output_counter(this%post_processors_setup_file_names(pproc_counter))
		end do
	end subroutine
	
	subroutine create_post_processor(this,manager,post_processor_name,operations_number,save_time,save_time_units)
		class(post_processor_manager)	,intent(inout)	:: this
		type(data_manager)				,intent(in)	:: manager
		integer							,intent(in)	:: operations_number
		real(dkind)						,intent(in)	:: save_time
		character(len=*)				,intent(in)	:: save_time_units
		character(len=*)				,intent(in)	:: post_processor_name
		
		character(len=30)	:: post_processor_output_file
		
		this%post_processors_created = this%post_processors_created + 1	
		
		post_processor_output_file = post_processor_name // post_processor_data_format
		this%post_processors(this%post_processors_created)	= post_processor_c(manager,post_processor_output_file,operations_number,save_time,save_time_units,this%post_processors_setup_file_names(this%post_processors_created))

	end subroutine
	
	subroutine create_post_processor_operation(this,manager,post_processor_number,field_name,operation_type,operation_area,operation_area_distance,grad_projection)
		class(post_processor_manager)	,intent(inout)	:: this
		type(data_manager)				,intent(in)	:: manager
		integer							,intent(in)	:: post_processor_number
		character(len=*)				,intent(in)	:: field_name
		character(len=*)				,intent(in)	:: operation_type
		integer	,dimension(3,2)			,intent(in)	,optional	:: operation_area	
		integer	,dimension(3)			,intent(in)	,optional	:: operation_area_distance
		integer							,intent(in)	,optional	:: grad_projection
		
		integer	,dimension(3,2)	:: op_area 
		integer	,dimension(3)	:: op_distance
		integer					:: grad_p
		
		if(present(operation_area)) then
			op_area = operation_area
		else
			op_area = 0
		end if
		
		if(present(operation_area_distance)) then
			op_distance = operation_area_distance
		else
			op_distance = 0
		end if
		
		if(present(grad_projection)) then
			grad_p = grad_projection
		else
			grad_p = 0
		end if		
		
		call this%post_processors(post_processor_number)%create_post_processor_operation(manager,field_name,operation_type,op_area,op_distance,grad_p,this%post_processors_setup_file_names(post_processor_number))
		
	end subroutine
	
	subroutine process_data(this,time,stop_flag)
		class(post_processor_manager)	,intent(inout)	:: this
		real(dkind)						,intent(in)		:: time
		logical							,intent(in)		:: stop_flag
		integer	:: number_post_processors
		logical	:: file_exists
		integer	:: io_unit_pproc
		integer	:: pproc

		integer					:: processor_rank

		processor_rank = this%domain%get_processor_rank()

		number_post_processors		= size(this%post_processors)

		do pproc = 1, number_post_processors
			call this%post_processors(pproc)%process_data(time)
		end do
		
		if((stop_flag).and.(processor_rank == 0)) call this%write_output_counters()
	
	end subroutine
	
end module
