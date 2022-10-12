module data_io_class
	
	use kind_parameters
	use global_data
	use computational_domain_class
	use computational_mesh_class
	use boundary_conditions_class
	use field_pointers
	use data_manager_class

#ifdef mpi
    use MPI
#endif

	implicit none
	
	private
	public	::	data_io, data_io_c	
		
	type data_io
		private
		type(computational_domain)			:: domain
		type(computational_mesh_pointer)	:: mesh
		type(boundary_conditions_pointer)	:: boundaries
		
		type(field_scalar_cons_pointer)	,dimension(:)	,allocatable	:: scalar_io_fields
		type(field_vector_cons_pointer)	,dimension(:)	,allocatable	:: vector_io_fields
		type(field_tensor_cons_pointer)	,dimension(:)	,allocatable	:: tensor_io_fields
	
		type(field_scalar_flow_pointer)	,dimension(:)	,allocatable	:: scalar_flow_io_fields
		type(field_vector_flow_pointer)	,dimension(:)	,allocatable	:: vector_flow_io_fields

		character(len=35)			,dimension(:)		,allocatable	:: scalar_io_fields_names
		character(len=35)			,dimension(:)		,allocatable	:: vector_io_fields_names
		character(len=35)			,dimension(:)		,allocatable	:: tensor_io_fields_names
		
		character(len=35)			,dimension(:)		,allocatable	:: scalar_flow_io_fields_names
		character(len=35)			,dimension(:)		,allocatable	:: vector_flow_io_fields_names		

		real(dkind)					,dimension(:,:,:)	,allocatable	:: io_cons_buffer
		real(dkind)					,dimension(:,:,:)	,allocatable	:: io_flow_buffer

		integer	,dimension(3,2)	:: io_cons_loop, io_cons_loop_max, io_flow_loop, io_flow_loop_max
		integer	,dimension(3,2)	:: io_cons_buffer_bounds
		integer	,dimension(3)	:: io_cons_buffer_size
		
		real(dkind)			:: check_time
		character(len=20)	:: check_time_units		
		character(len=10)   :: check_time_units_abbreviation
		real(rkind)         :: check_time_coefficient		
		
		real(dkind)			:: start_time
		real(dkind)			:: output_time
		
		integer				:: load_counter
		character(len=20)	:: data_output_folder

		integer				:: mpi_array_bounds_type
		integer				:: mpi_array_nobounds_type
		integer				:: mpi_flow_array_bounds_type
		integer				:: mpi_flow_array_nobounds_type		
	contains
		! IO Getters/Setters
		procedure	,private	:: read_properties
		procedure	,private	:: write_properties
		procedure	,private	:: set_properties
		
		! Setters
		procedure	:: set_start_time
		procedure	:: add_io_scalar_cons_field
		procedure	:: add_io_scalar_flow_field
		procedure	:: add_io_vector_flow_field

		! Getters
		procedure	:: get_load_counter
		procedure	:: get_data_output_folder

		procedure				:: output_all_data
		! Sequential data output	
		procedure	,private	:: output_fields
		procedure	,private	:: output_mesh
		procedure	,private	:: output_bounds
		! Parallel data output
		procedure	,private	:: output_fields_mpi
		procedure	,private	:: output_mesh_mpi
		procedure	,private	:: output_bounds_mpi		

		procedure				:: input_all_data
		! Sequential data input		
		procedure	,private	:: input_fields
		procedure	,private	:: input_mesh
		procedure	,private	:: input_bounds	

		! Prallel data input
		procedure	,private	:: input_mesh_mpi
		procedure	,private	:: input_fields_mpi
		procedure	,private	:: input_bounds_mpi

		procedure	,private	:: commit_mpi_io_datatype	

		! Logger
		procedure	:: write_log
		procedure	:: write_load_counter	

	end type
	
	interface data_io_c
		module procedure constructor
		module procedure constructor_file
	end interface
	
contains

	type(data_io)	function constructor(manager,check_time,check_time_units,output_time,data_output_folder,io_fields_names)
		type(data_manager)					,intent(in)	:: manager
		real(dkind)							,intent(in)	:: check_time
		character(len=*)					,intent(in)	:: check_time_units
		real(dkind)							,intent(in)	:: output_time
		character(len=*)					,intent(in)	:: data_output_folder
		character(len=35)	,dimension(:)	,intent(in)	,optional	:: io_fields_names
					
		character(len=10)	:: field_type
		integer	:: load_counter, fields_counter
		integer	:: io_unit
				
		load_counter = 0
		
		if(present(io_fields_names)) then
			call constructor%set_properties(manager,check_time,check_time_units,output_time,data_output_folder,load_counter,io_fields_names)
		else
			call constructor%set_properties(manager,check_time,check_time_units,output_time,data_output_folder,load_counter)
		end if
		
		open(newunit = io_unit, file = data_io_data_file_name, status = 'replace', form = 'formatted', delim = 'quote')
		call constructor%write_properties(io_unit)
		close(io_unit)

	end function 
	
	type(data_io) function constructor_file(manager,calculation_time)
		type(data_manager)			,intent(in)		:: manager
		real(dkind)					,intent(out)	:: calculation_time
		
		integer	:: io_unit
		
		open(newunit = io_unit, file = data_io_data_file_name, status = 'old', form = 'formatted', delim = 'quote')
		call constructor_file%read_properties(manager,io_unit,calculation_time)
		close(io_unit)	

	end function
	
	subroutine read_properties(this,manager,data_io_file_unit,calculation_time)
		class(data_io)		,intent(inout)	:: this
		type(data_manager)	,intent(in)		:: manager
		integer				,intent(in)		:: data_io_file_unit
		real(dkind)			,intent(out)	:: calculation_time
		character(len=35)	,dimension(:)	,allocatable	:: io_fields_names
		
		integer				:: io_fields_number
		real(dkind)			:: check_time
		character(len=20)	:: check_time_units
		real(dkind)			:: output_time
		character(len=20)	:: data_output_folder		
		integer				:: load_counter	
		integer				:: ierr
		
		namelist /data_io_parameters/ io_fields_number, check_time,check_time_units, output_time, data_output_folder	
		namelist /io_fields/ io_fields_names
		namelist /data_io_load_counter/ load_counter, calculation_time
		
		read(unit = data_io_file_unit, nml = data_io_parameters)
		allocate(io_fields_names(io_fields_number))
		read(unit = data_io_file_unit, nml = io_fields)
		read(unit = data_io_file_unit, nml = data_io_load_counter, iostat = ierr)
		if(ierr /= 0) then
			load_counter = 0
			calculation_time = 0.0_dkind
		end if
		
		call this%set_properties(manager,check_time,check_time_units,output_time,data_output_folder,load_counter,io_fields_names)	
		call this%commit_mpi_io_datatype()	
	end subroutine
	
	subroutine write_properties(this,data_io_file_unit)
		class(data_io)		,intent(in)	:: this
		integer				,intent(in)	:: data_io_file_unit		
	
		character(len=35)	,dimension(:)	,allocatable	:: io_fields_names
		
		integer				:: io_fields_number
		real(dkind)			:: check_time
		character(len=20)	:: check_time_units		
		real(rkind)			:: output_time
		character(len=20)	:: data_output_folder		
		integer				:: load_counter	
		real(dkind)			:: calculation_time = 0.0_dkind
		
		namelist /data_io_parameters/ io_fields_number, check_time,check_time_units, output_time, data_output_folder	
		namelist /io_fields/ io_fields_names
		
		io_fields_number	= size(this%scalar_io_fields_names) + size(this%vector_io_fields_names) + size(this%tensor_io_fields_names)
		check_time			= this%check_time
		check_time_units	= this%check_time_units
		output_time			= this%output_time
		data_output_folder	= this%data_output_folder
		
		allocate(io_fields_names(io_fields_number))
		io_fields_names = [this%scalar_io_fields_names,this%vector_io_fields_names,this%tensor_io_fields_names]

		write(unit = data_io_file_unit, nml = data_io_parameters)
		write(unit = data_io_file_unit, nml = io_fields)

	end subroutine
	
	subroutine set_properties(this,manager,check_time,check_time_units,output_time,data_output_folder,load_counter,io_fields_names)
		class(data_io)		,intent(inout)	:: this
		type(data_manager)	,intent(in)		:: manager
		real(dkind)			,intent(in)		:: check_time
		character(len=*)	,intent(in)		:: check_time_units
		real(dkind)			,intent(in)		:: output_time
		character(len=*)	,intent(in)		:: data_output_folder		
		integer				,intent(in)		:: load_counter	

		character(len=35)	,dimension(:)	,intent(in), optional	:: io_fields_names		

		character(len=100)              :: system_command

		type(field_scalar_cons_pointer)	:: scal_ptr
		type(field_vector_cons_pointer)	:: vect_ptr
		type(field_tensor_cons_pointer)	:: tens_ptr		
		
		character(len=10)	:: field_type
		integer	:: scalar_fields = 0
		integer	:: vector_fields = 0
		integer	:: tensor_fields = 0

		integer	,dimension(3,2)	:: cons_allocation_bounds, flow_allocation_bounds
		integer	,dimension(3,2)	:: utter_loop
		integer	,dimension(3)	:: local_cells_number, grid_coord, grid_size

		integer	:: processor_rank

		integer	:: fields_counter
		integer	:: dimensions, dim

		this%check_time			= check_time
		this%check_time_units	= check_time_units
		this%output_time		= output_time
		this%data_output_folder = data_output_folder
		this%load_counter		= load_counter
		this%domain				= manager%domain
		
		this%mesh 				= manager%computational_mesh_pointer
		this%boundaries			= manager%boundary_conditions_pointer		
		
		cons_allocation_bounds	= this%domain%get_local_utter_cells_bounds()
		flow_allocation_bounds	= this%domain%get_local_utter_faces_bounds()

		allocate(this%io_cons_buffer(	cons_allocation_bounds(1,1):cons_allocation_bounds(1,2), &
										cons_allocation_bounds(2,1):cons_allocation_bounds(2,2), &
										cons_allocation_bounds(3,1):cons_allocation_bounds(3,2)))

		allocate(this%io_flow_buffer(	flow_allocation_bounds(1,1):flow_allocation_bounds(1,2), &
										flow_allocation_bounds(2,1):flow_allocation_bounds(2,2), &
										flow_allocation_bounds(3,1):flow_allocation_bounds(3,2)))									
									
		select case(this%check_time_units)
			 case('milliseconds')
				  this%check_time_units_abbreviation = 'ms'
				  this%check_time_coefficient        = 1e+03_dkind
			 case('microseconds')
				  this%check_time_units_abbreviation = 'us'
				  this%check_time_coefficient        = 1e+06_dkind
			 case('nanoseconds')
				  this%check_time_units_abbreviation = 'ns'
				  this%check_time_coefficient        = 1e+09_dkind
		end select									
									
		if(present(io_fields_names)) then
			do fields_counter = 1,size(io_fields_names)
				field_type = manager%get_field_type_by_name(io_fields_names(fields_counter))
				select case (field_type)
					case('scalar')
						scalar_fields = scalar_fields + 1
					case('vector')
						vector_fields = vector_fields + 1
					case('tensor')
						tensor_fields = tensor_fields + 1
				end select
			end do
			allocate(this%scalar_io_fields(scalar_fields))
			allocate(this%vector_io_fields(vector_fields))
			allocate(this%tensor_io_fields(tensor_fields))
			allocate(this%scalar_io_fields_names(scalar_fields))
			allocate(this%vector_io_fields_names(vector_fields))
			allocate(this%tensor_io_fields_names(tensor_fields))			
			scalar_fields = 0
			vector_fields = 0
			tensor_fields = 0
			do fields_counter = 1,size(io_fields_names)
				call manager%get_cons_field_pointer_by_name(scal_ptr,vect_ptr,tens_ptr,io_fields_names(fields_counter))
				if(associated(scal_ptr%s_ptr)) then
					scalar_fields = scalar_fields + 1 
					this%scalar_io_fields(scalar_fields)%s_ptr => scal_ptr%s_ptr
					this%scalar_io_fields_names(scalar_fields) = io_fields_names(fields_counter)
				end if
				if(associated(vect_ptr%v_ptr)) then
					vector_fields = vector_fields + 1
					this%vector_io_fields(vector_fields)%v_ptr => vect_ptr%v_ptr
					this%vector_io_fields_names(vector_fields) = io_fields_names(fields_counter)
				end if
				if(associated(tens_ptr%t_ptr)) then
					tensor_fields = tensor_fields + 1
					this%tensor_io_fields(tensor_fields)%t_ptr => tens_ptr%t_ptr
					this%tensor_io_fields_names(tensor_fields) = io_fields_names(fields_counter)
				end if				
			end do			
		else
			
			scalar_fields = manager%get_number_of_cons_scalar_fields()
			vector_fields = manager%get_number_of_cons_vector_fields()
			tensor_fields = manager%get_number_of_cons_tensor_fields()
			
			allocate(this%scalar_io_fields(scalar_fields))
			allocate(this%vector_io_fields(vector_fields))
			allocate(this%tensor_io_fields(tensor_fields))
			allocate(this%scalar_io_fields_names(scalar_fields))
			allocate(this%vector_io_fields_names(vector_fields))
			allocate(this%tensor_io_fields_names(tensor_fields))
			
			do fields_counter = 1, scalar_fields
				call manager%get_cons_field_pointer_by_number(scal_ptr,vect_ptr,tens_ptr,'scalar',fields_counter)
				this%scalar_io_fields(fields_counter)%s_ptr => scal_ptr%s_ptr
				this%scalar_io_fields_names(fields_counter) = scal_ptr%s_ptr%name_long
			end do
			do fields_counter = 1, vector_fields
				call manager%get_cons_field_pointer_by_number(scal_ptr,vect_ptr,tens_ptr,'vector',fields_counter)
				this%vector_io_fields(fields_counter)%v_ptr => vect_ptr%v_ptr
				this%vector_io_fields_names(fields_counter) = vect_ptr%v_ptr%name_long
			end do			
			do fields_counter = 1, tensor_fields
				call manager%get_cons_field_pointer_by_number(scal_ptr,vect_ptr,tens_ptr,'tensor',fields_counter)
				this%tensor_io_fields(fields_counter)%t_ptr => tens_ptr%t_ptr
				this%tensor_io_fields_names(fields_counter) = tens_ptr%t_ptr%name_long
			end do				
		end if

		dimensions = this%domain%get_domain_dimensions()

		grid_coord			= this%domain%get_processor_grid_coord()
		grid_size			= this%domain%get_processor_number()

		local_cells_number	= this%domain%get_local_cells_number()
		utter_loop					= this%domain%get_local_utter_cells_bounds()

		this%io_cons_loop			= utter_loop
		this%io_cons_loop_max		= this%domain%get_local_utter_cells_bounds_max()
		this%io_cons_buffer_bounds 	= 1
		this%io_cons_buffer_size	= local_cells_number

		do dim = 1,dimensions
			if((grid_coord(dim) == 0).and.(grid_coord(dim) /= grid_size(dim)-1)) then
				this%io_cons_loop(dim,1) 			= utter_loop(dim,1)
				this%io_cons_loop(dim,2)			= utter_loop(dim,2) - 1
				this%io_cons_buffer_size(dim)		= local_cells_number(dim) - 1				
				this%io_cons_buffer_bounds(dim,1)	= 0
				this%io_cons_buffer_bounds(dim,2)	= local_cells_number(dim) - 2 				
			end if
			if((grid_coord(dim) /= 0).and.(grid_coord(dim) == grid_size(dim)-1)) then
				this%io_cons_loop(dim,1) 			= utter_loop(dim,1) + 1
				this%io_cons_loop(dim,2)			= utter_loop(dim,2)
				this%io_cons_buffer_size(dim)		= local_cells_number(dim) - 1				
				this%io_cons_buffer_bounds(dim,1)	= 1
				this%io_cons_buffer_bounds(dim,2)	= local_cells_number(dim) - 1 				
			end if
			if((grid_coord(dim) /= 0).and.(grid_coord(dim) /= grid_size(dim)-1)) then	
				this%io_cons_loop(dim,1) 			= utter_loop(dim,1) + 1
				this%io_cons_loop(dim,2)			= utter_loop(dim,2) - 1
				this%io_cons_buffer_size(dim)  		= local_cells_number(dim) - 2				
				this%io_cons_buffer_bounds(dim,1)	= 1
				this%io_cons_buffer_bounds(dim,2)	= local_cells_number(dim) - 2 				
			end if	
			if((grid_coord(dim) == 0).and.(grid_coord(dim) == grid_size(dim)-1)) then			
				this%io_cons_loop(dim,1) 			= utter_loop(dim,1)
				this%io_cons_loop(dim,2)			= utter_loop(dim,2)
				this%io_cons_buffer_size(dim)  		= local_cells_number(dim)
				this%io_cons_buffer_bounds(dim,1)	= 0
				this%io_cons_buffer_bounds(dim,2)	= local_cells_number(dim) - 1 				
			end if
		end do

		processor_rank = this%domain%get_processor_rank()
		
		if (processor_rank == 0) then
			system_command = 'mkdir ' // trim(data_output_folder)
			call system(system_command)
		end if

	end subroutine
	
	subroutine set_start_time(this)
		class(data_io)	,intent(inout)	:: this
		
		integer	:: c(8)
		integer	:: day,h,m,s,t

		integer	:: processor_rank
		integer	:: error

		processor_rank = this%domain%get_processor_rank()

#ifdef mpi
		call MPI_BARRIER(MPI_COMM_WORLD,error)
#endif		

		call date_and_time(values=c)
		day=c(3)
		h=c(5)
		m=c(6)
		s=c(7)
		t=c(8)
		
		this%start_time = ((day*24+h)*60+m)

		if (processor_rank == 0) then
			print *, 'Calculations started.'
			print *, 'Current time = ', day,'  ',h,':',m
		end if

	end subroutine
	
	subroutine write_log(this,log_unit)
		class(data_io)	,intent(in)	:: this
		integer			,intent(in)	:: log_unit
		
		write(log_unit,'(A)')	'************************************************************************************* '	
		write(log_unit,'(A)')			' Data output setup:'
		write(log_unit,'(A,E14.7,A)')	' Data output time (minutes): ', this%output_time
		write(log_unit,'(A,A)')			' Data save folder: ', this%data_output_folder
		write(log_unit,'(A,E14.7,A)')	' Check time: ', this%check_time, this%check_time_units_abbreviation 
		write(log_unit,'(A)')	'************************************************************************************* '
	end subroutine
	
	subroutine write_load_counter(this,calculation_time)
		class(data_io)	,intent(in)	:: this
		real(dkind)		,intent(in)	:: calculation_time
		
		integer	:: load_counter
		
		integer	:: io_unit
		
		namelist /data_io_load_counter/ load_counter,calculation_time
		
		load_counter = this%load_counter
		
		open(newunit = io_unit, file = data_io_data_file_name, status = 'old', form = 'formatted', position = 'append')
		if (load_counter /= 1) then
			backspace(io_unit)
			backspace(io_unit)
			backspace(io_unit)
			backspace(io_unit)
		end if
		write(unit = io_unit, nml = data_io_load_counter)
		close(io_unit)
	end subroutine
	
	subroutine output_all_data(this,time,stop_flag,make_output)
		class(data_io)	,intent(inout)				:: this
		real(dkind)		,intent(in)					:: time
		logical			,intent(in)		,optional	:: make_output		
		logical			,intent(out)				:: stop_flag

		integer				:: processor_rank

		character(len=100)	:: system_command
		integer		,save	:: check_counter = 0
		logical		:: make_flag
		real(dkind)	:: current_time, output_time
		integer		:: c(8)
		integer		:: day,h,m,s,t
	
		integer		:: error

		processor_rank = this%domain%get_processor_rank()

!		stop_flag = .false.
		make_flag = .false.
		if(present(make_output)) make_flag = make_output

		if((time*this%check_time_coefficient >= this%check_time*check_counter).or.(make_flag)) then
		
			check_counter = check_counter + 1
		
#ifdef mpi
			call MPI_BARRIER(MPI_COMM_WORLD,error)
#endif	

			call date_and_time(values=c)
			day=c(3)
			h=c(5)
			m=c(6)
			s=c(7)
			t=c(8)
			current_time = ((day*24+h)*60+m)
			
			if (processor_rank == 0) then
				print *, 'Time left = ', this%start_time + this%output_time - current_time, ' min'
			end if

			if (((current_time - this%start_time) >= this%output_time).or.(make_flag)) then
		
				if (processor_rank == 0) then
					print *, 'Performing data output.'
					print *, 'Current time = ', day,'  ',h,':',m
					
					write(system_command,'(A,A,A,I3.3)') 'mkdir ', trim(this%data_output_folder), trim(fold_sep), this%load_counter
					call system(system_command)	
				end if			

#ifdef mpi
				call MPI_BARRIER(MPI_COMM_WORLD,error)
				! call this%output_fields_mpi()
				! call this%output_mesh_mpi()
				! call this%output_bounds_mpi()
				call this%output_fields()
				call this%output_mesh()
				call this%output_bounds()
				call MPI_BARRIER(MPI_COMM_WORLD,error)
#else			
				call this%output_fields()
				call this%output_mesh()
				call this%output_bounds()
#endif

				if (processor_rank == 0) then
					call date_and_time(values=c)
					day=c(3)
					h=c(5)
					m=c(6)
					s=c(7)
					t=c(8)
					output_time = ((day*24+h)*60+m) - current_time
					print *, 'Data output successfully finished.'
					print *, 'Output lasts ', output_time ,' min'
				end if	
		
				stop_flag = .true.
			
				if(present(make_output)) then
					if(make_output) stop_flag = .false.
				end if
			end if
		end if

		if(stop_flag.and.(processor_rank == 0)) call this%write_load_counter(time)
	end subroutine
	
	subroutine output_fields(this)
		class(data_io)	,intent(inout)	:: this
		
		integer				:: io_unit
		character(len=100)	:: file_name
		character(len=125)	:: system_command

		integer						:: dimensions
		integer						:: vector_projections_number
		integer		,dimension(2)	:: tensor_projections_number
		integer		,dimension(3,2)	:: utter_loop, inner_loop

		integer				:: processor_rank

		integer	:: i,j,k, fields_counter, dim, dim1, dim2

		integer	:: error

		utter_loop		= this%domain%get_local_utter_cells_bounds()

		inner_loop		= this%domain%get_local_inner_faces_bounds()
		
		dimensions		= this%domain%get_domain_dimensions()
		
		processor_rank = this%domain%get_processor_rank()

		do fields_counter = 1,size(this%scalar_io_fields)
#ifdef mpi
			if (processor_rank == 0) then
				write(system_command,'(A,A,A,I3.3,A,A)') 'mkdir ', trim(this%data_output_folder) , trim(fold_sep) , this%load_counter, trim(fold_sep) , trim(this%scalar_io_fields_names(fields_counter))
				call system(system_command)	
			end if	
			call MPI_BARRIER(MPI_COMM_WORLD,error)			
			write(file_name,'(A,A,I3.3,A,A,A,A,A,I4.4,A)') trim(this%data_output_folder) , trim(fold_sep) , this%load_counter, trim(fold_sep) , trim(this%scalar_io_fields_names(fields_counter)) , trim(fold_sep), trim(this%scalar_io_fields_names(fields_counter)), '_proc_', processor_rank, trim(data_io_data_format)		
#else
			write(file_name,'(A,A,I3.3,A,A,A)') trim(this%data_output_folder) , trim(fold_sep) , this%load_counter, trim(fold_sep) , trim(this%scalar_io_fields_names(fields_counter)) , trim(data_io_data_format)		
#endif
			open(newunit = io_unit, file = file_name, status = 'replace', form = 'binary')
			this%io_cons_buffer = this%scalar_io_fields(fields_counter)%s_ptr%cells(:,:,:)
			write(io_unit) (((this%io_cons_buffer(i,j,k),	i = utter_loop(1,1),utter_loop(1,2)), &
															j = utter_loop(2,1),utter_loop(2,2)), &
															k = utter_loop(3,1),utter_loop(3,2))
			close(io_unit)
		end do	
		
		do fields_counter = 1,size(this%vector_io_fields)
#ifdef mpi		
			if (processor_rank == 0) then
				write(system_command,'(A,A,A,I3.3,A,A)') 'mkdir ', trim(this%data_output_folder) , trim(fold_sep) , this%load_counter, trim(fold_sep) , trim(this%vector_io_fields_names(fields_counter))
				call system(system_command)	
			end if	
			call MPI_BARRIER(MPI_COMM_WORLD,error)	
			write(file_name,'(A,A,I3.3,A,A,A,A,A,I4.4,A)') trim(this%data_output_folder) , trim(fold_sep) , this%load_counter, trim(fold_sep) , trim(this%vector_io_fields_names(fields_counter)) , trim(fold_sep), trim(this%vector_io_fields_names(fields_counter)), '_proc_', processor_rank, trim(data_io_data_format)		
#else
			write(file_name,'(A,A,I3.3,A,A,A)') trim(this%data_output_folder) , trim(fold_sep) , this%load_counter, trim(fold_sep) , trim(this%vector_io_fields_names(fields_counter)) , trim(data_io_data_format)
#endif
			open(newunit = io_unit, file = file_name, status = 'replace', form = 'binary')
			vector_projections_number = this%vector_io_fields(fields_counter)%v_ptr%get_projections_number()
			do dim = 1, vector_projections_number
				this%io_cons_buffer = this%vector_io_fields(fields_counter)%v_ptr%pr(dim)%cells(:,:,:)
				write(io_unit) (((this%io_cons_buffer(i,j,k),	i = utter_loop(1,1),utter_loop(1,2)), &
																j = utter_loop(2,1),utter_loop(2,2)), &
																k = utter_loop(3,1),utter_loop(3,2))
			end do
			close(io_unit)
		end do		
	
		do fields_counter = 1,size(this%tensor_io_fields)
#ifdef mpi
			if (processor_rank == 0) then
				write(system_command,'(A,A,A,I3.3,A,A)') 'mkdir ', trim(this%data_output_folder) , trim(fold_sep) , this%load_counter, trim(fold_sep) , trim(this%tensor_io_fields_names(fields_counter))
				call system(system_command)	
			end if	
			call MPI_BARRIER(MPI_COMM_WORLD,error)	
			write(file_name,'(A,A,I3.3,A,A,A,A,A,I4.4,A)') trim(this%data_output_folder) , trim(fold_sep) , this%load_counter, trim(fold_sep) , trim(this%tensor_io_fields_names(fields_counter)) , trim(fold_sep), trim(this%tensor_io_fields_names(fields_counter)), '_proc_', processor_rank, trim(data_io_data_format)		
#else	
			write(file_name,'(A,A,I3.3,A,A,A)') trim(this%data_output_folder) , trim(fold_sep) , this%load_counter, trim(fold_sep) , trim(this%tensor_io_fields_names(fields_counter)) , trim(data_io_data_format)
#endif
			open(newunit = io_unit, file = file_name, status = 'replace', form = 'binary')
			tensor_projections_number = this%tensor_io_fields(fields_counter)%t_ptr%get_projections_number()
			do dim1 = 1, tensor_projections_number(1)
			do dim2 = 1, tensor_projections_number(2)
				this%io_cons_buffer = this%tensor_io_fields(fields_counter)%t_ptr%pr(dim1,dim2)%cells(:,:,:)
				write(io_unit) (((this%io_cons_buffer(i,j,k),	i = utter_loop(1,1),utter_loop(1,2)), &
																j = utter_loop(2,1),utter_loop(2,2)), &
																k = utter_loop(3,1),utter_loop(3,2))
			end do
			end do
			close(io_unit)
		end do		

		if (allocated(this%scalar_flow_io_fields)) then
			do fields_counter = 1,size(this%scalar_flow_io_fields)
#ifdef mpi		
				if (processor_rank == 0) then
					write(system_command,'(A,A,A,I3.3,A,A)') 'mkdir ', trim(this%data_output_folder) , trim(fold_sep) , this%load_counter, trim(fold_sep) , trim(this%scalar_flow_io_fields_names(fields_counter))
					call system(system_command)	
				end if	
				call MPI_BARRIER(MPI_COMM_WORLD,error)
				write(file_name,'(A,A,I3.3,A,A,A,A,A,I4.4,A)') trim(this%data_output_folder) , trim(fold_sep) , this%load_counter, trim(fold_sep) , trim(this%scalar_flow_io_fields_names(fields_counter)) , trim(fold_sep), trim(this%scalar_flow_io_fields_names(fields_counter)), '_proc_', processor_rank, trim(data_io_data_format)		
#else			
				write(file_name,'(A,A,I3.3,A,A,A)') trim(this%data_output_folder) , trim(fold_sep) , this%load_counter, trim(fold_sep) , trim(this%scalar_flow_io_fields_names(fields_counter)) , trim(data_io_data_format)		
#endif				
				open(newunit = io_unit, file = file_name, status = 'replace', form = 'binary')
				
				do dim = 1,dimensions
					this%io_flow_buffer = this%scalar_flow_io_fields(fields_counter)%s_ptr%cells(dim,:,:,:)

					write(io_unit) (((this%io_flow_buffer(i,j,k),	i = inner_loop(1,1),inner_loop(1,2)), &
																	j = inner_loop(2,1),inner_loop(2,2)), &
																	k = inner_loop(3,1),inner_loop(3,2))
				end do
				close(io_unit)
			end do
		end if

		if (allocated(this%vector_flow_io_fields)) then
			do fields_counter = 1,size(this%vector_flow_io_fields)
#ifdef mpi						
				if (processor_rank == 0) then
					write(system_command,'(A,A,A,I3.3,A,A)') 'mkdir ', trim(this%data_output_folder) , trim(fold_sep) , this%load_counter, trim(fold_sep) , trim(this%vector_flow_io_fields_names(fields_counter))
					call system(system_command)	
				end if	
				call MPI_BARRIER(MPI_COMM_WORLD,error)
				write(file_name,'(A,A,I3.3,A,A,A,A,A,I4.4,A)') trim(this%data_output_folder) , trim(fold_sep) , this%load_counter, trim(fold_sep) , trim(this%vector_flow_io_fields_names(fields_counter)) , trim(fold_sep), trim(this%vector_flow_io_fields_names(fields_counter)), '_proc_', processor_rank, trim(data_io_data_format)		
#else			
				write(file_name,'(A,A,I3.3,A,A,A)') trim(this%data_output_folder) , trim(fold_sep) , this%load_counter, trim(fold_sep) , trim(this%vector_flow_io_fields_names(fields_counter)) , trim(data_io_data_format)
#endif					
				open(newunit = io_unit, file = file_name, status = 'replace', form = 'binary')

				vector_projections_number = this%vector_flow_io_fields(fields_counter)%v_ptr%get_projections_number_flow()
				
				do dim1 = 1, vector_projections_number
				do dim2 = 1, dimensions
					this%io_flow_buffer = this%vector_flow_io_fields(fields_counter)%v_ptr%pr(dim1)%cells(dim2,:,:,:)

					write(io_unit) (((this%io_flow_buffer(i,j,k),	i = inner_loop(1,1),inner_loop(1,2)), &
																	j = inner_loop(2,1),inner_loop(2,2)), &
																	k = inner_loop(3,1),inner_loop(3,2))
				end do
				end do

				close(io_unit)
			end do
		end if		
		
		
	end subroutine
	
	subroutine output_mesh(this)
		class(data_io)	,intent(inout)	:: this
		
		integer				:: io_unit
		character(len=100)	:: file_name
		character(len=125)	:: system_command

		integer						:: dimensions
		integer		,dimension(3,2)	:: utter_loop

		integer	:: processor_rank

		integer	:: i,j,k, fields_counter, dim

		integer	:: error

		dimensions		= this%domain%get_domain_dimensions()
		utter_loop		= this%domain%get_local_utter_cells_bounds()

		processor_rank	= this%domain%get_processor_rank()

#ifdef mpi	
		if (processor_rank == 0) then
			write(system_command,'(A,A,A,I3.3,A,A)') 'mkdir ', trim(this%data_output_folder) , trim(fold_sep) , this%load_counter, trim(fold_sep) , 'mesh'
			call system(system_command)	
		end if	
		call MPI_BARRIER(MPI_COMM_WORLD,error)
		write(file_name,'(A,A,I3.3,A,A,A,A,A,I4.4,A)') trim(this%data_output_folder) , trim(fold_sep) , this%load_counter, trim(fold_sep) , 'mesh' , trim(fold_sep), 'mesh', '_proc_', processor_rank, trim(data_io_data_format)		
#else	
		write(file_name,'(A,A,I3.3,A,A,A)') trim(this%data_output_folder) , trim(fold_sep) , this%load_counter, trim(fold_sep) , 'mesh' , trim(data_io_data_format)
#endif		
		open(newunit = io_unit, file = file_name, status = 'replace', form = 'binary')
		do dim = 1, dimensions
			this%io_cons_buffer = this%mesh%mesh_ptr%mesh(dim,:,:,:)
			write(io_unit) (((this%io_cons_buffer(i,j,k),	i = utter_loop(1,1),utter_loop(1,2)), &
															j = utter_loop(2,1),utter_loop(2,2)), &
															k = utter_loop(3,1),utter_loop(3,2))
		end do
		close(io_unit)
	end subroutine	
	
	subroutine output_bounds(this)
		class(data_io)	,intent(inout)	:: this
		
		integer				:: io_unit
		character(len=100)	:: file_name
		character(len=125)	:: system_command

		integer		,dimension(3,2)	:: utter_loop

		integer	:: processor_rank

		integer	:: i,j,k, fields_counter, dim

		integer	:: error

		utter_loop		= this%domain%get_local_utter_cells_bounds()

		processor_rank	= this%domain%get_processor_rank()

#ifdef mpi
		if (processor_rank == 0) then
			write(system_command,'(A,A,A,I3.3,A,A)') 'mkdir ', trim(this%data_output_folder) , trim(fold_sep) , this%load_counter, trim(fold_sep) , 'bc'
			call system(system_command)	
		end if	
		call MPI_BARRIER(MPI_COMM_WORLD,error)	
		write(file_name,'(A,A,I3.3,A,A,A,A,A,I4.4,A)') trim(this%data_output_folder) , trim(fold_sep) , this%load_counter, trim(fold_sep) , 'bc' , trim(fold_sep), 'bc', '_proc_', processor_rank, trim(data_io_data_format)		
#else
		write(file_name,'(A,A,I3.3,A,A,A)') trim(this%data_output_folder) , trim(fold_sep) , this%load_counter, trim(fold_sep) , 'bc' , trim(data_io_data_format)
#endif		
		open(newunit = io_unit, file = file_name, status = 'replace', form = 'binary')
		this%io_cons_buffer = this%boundaries%bc_ptr%bc_markers(:,:,:)
		write(io_unit) (((this%io_cons_buffer(i,j,k),	i = utter_loop(1,1),utter_loop(1,2)), &
														j = utter_loop(2,1),utter_loop(2,2)), &
														k = utter_loop(3,1),utter_loop(3,2))

		close(io_unit)
	end subroutine	
	
	subroutine input_all_data(this)
		class(data_io)	,intent(inout)	:: this

		integer	:: processor_rank
		integer	:: io_unit, error

		call this%set_start_time()		
		
		processor_rank = this%domain%get_processor_rank()
		
		if(processor_rank == 0) then
			print *, 'Performing data input.'
		end if
		
#ifdef mpi
		call this%input_mesh_mpi()
		call this%input_fields_mpi()
		call this%input_bounds_mpi()
#else
		call this%input_fields()
		call this%input_mesh()
		call this%input_bounds()
#endif

		this%load_counter = this%load_counter + 1
	end subroutine
	
	subroutine input_fields(this)
		class(data_io)	,intent(inout)	:: this
		
		integer				:: io_unit
		character(len=65)	:: file_name
		
		integer						:: dimensions
		integer		,dimension(3,2)	:: utter_loop, inner_loop
		integer						:: vector_projections_number
		integer		,dimension(2)	:: tensor_projections_number

		integer	:: i,j,k, fields_counter, dim, dim1, dim2

		utter_loop	= this%domain%get_local_utter_cells_bounds()

		inner_loop	= this%domain%get_local_inner_faces_bounds()
		
		dimensions	= this%domain%get_domain_dimensions()
		
		do fields_counter = 1,size(this%scalar_io_fields)
			write(file_name,'(A,A,I3.3,A,A,A)') trim(this%data_output_folder) , trim(fold_sep) , this%load_counter, trim(fold_sep) , trim(this%scalar_io_fields_names(fields_counter)) , trim(data_io_data_format)		
			open(newunit = io_unit, file = file_name, status = 'old', form = 'binary')
			read(io_unit) (((this%io_cons_buffer(i,j,k),	i = utter_loop(1,1),utter_loop(1,2))	, &
															j = utter_loop(2,1),utter_loop(2,2)), &
															k = utter_loop(3,1),utter_loop(3,2))
			this%scalar_io_fields(fields_counter)%s_ptr%cells(:,:,:) = this%io_cons_buffer 
			close(io_unit)
		end do	
		
		do fields_counter = 1,size(this%vector_io_fields)
			write(file_name,'(A,A,I3.3,A,A,A)') trim(this%data_output_folder) , trim(fold_sep) , this%load_counter, trim(fold_sep) , trim(this%vector_io_fields_names(fields_counter)) , trim(data_io_data_format)
			open(newunit = io_unit, file = file_name, status = 'old', form = 'binary')
			vector_projections_number = this%vector_io_fields(fields_counter)%v_ptr%get_projections_number()
			do dim = 1, vector_projections_number
				read(io_unit) (((this%io_cons_buffer(i,j,k),	i = utter_loop(1,1),utter_loop(1,2)), &
																j = utter_loop(2,1),utter_loop(2,2)), &
																k = utter_loop(3,1),utter_loop(3,2))
				this%vector_io_fields(fields_counter)%v_ptr%pr(dim)%cells(:,:,:) = this%io_cons_buffer
			end do
			close(io_unit)
		end do		
	
		do fields_counter = 1,size(this%tensor_io_fields)
			write(file_name,'(A,A,I3.3,A,A,A)') trim(this%data_output_folder) , trim(fold_sep) , this%load_counter, trim(fold_sep) , trim(this%tensor_io_fields_names(fields_counter)) , trim(data_io_data_format)
			open(newunit = io_unit, file = file_name, status = 'old', form = 'binary')
			tensor_projections_number = this%tensor_io_fields(fields_counter)%t_ptr%get_projections_number()
			do dim1 = 1, tensor_projections_number(1)
			do dim2 = 1, tensor_projections_number(2)
				read(io_unit) (((this%io_cons_buffer(i,j,k),	i = utter_loop(1,1),utter_loop(1,2)), &
																j = utter_loop(2,1),utter_loop(2,2)), &
																k = utter_loop(3,1),utter_loop(3,2))
				this%tensor_io_fields(fields_counter)%t_ptr%pr(dim1,dim2)%cells(:,:,:) = this%io_cons_buffer									
			end do
			end do
			close(io_unit)
		end do		
		
		if (allocated(this%scalar_flow_io_fields)) then
			do fields_counter = 1,size(this%scalar_flow_io_fields)	
				write(file_name,'(A,A,I3.3,A,A,A)') trim(this%data_output_folder) , trim(fold_sep) , this%load_counter, trim(fold_sep) , trim(this%scalar_flow_io_fields_names(fields_counter)) , trim(data_io_data_format)		
				open(newunit = io_unit, file = file_name, status = 'old', form = 'binary')
				
				do dim = 1,dimensions
					read(io_unit) (((this%io_flow_buffer(i,j,k),	i = inner_loop(1,1),inner_loop(1,2)), &
																	j = inner_loop(2,1),inner_loop(2,2)), &
																	k = inner_loop(3,1),inner_loop(3,2))
					this%scalar_flow_io_fields(fields_counter)%s_ptr%cells(dim,:,:,:) = this%io_flow_buffer											
				end do
				close(io_unit)
			end do
		end if

		if (allocated(this%vector_flow_io_fields)) then
			do fields_counter = 1,size(this%vector_flow_io_fields)
				write(file_name,'(A,A,I3.3,A,A,A)') trim(this%data_output_folder) , trim(fold_sep) , this%load_counter, trim(fold_sep) , trim(this%vector_flow_io_fields_names(fields_counter)) , trim(data_io_data_format)
				open(newunit = io_unit, file = file_name, status = 'old', form = 'binary')

				vector_projections_number = this%vector_flow_io_fields(fields_counter)%v_ptr%get_projections_number_flow()
				
				do dim1 = 1, vector_projections_number
				do dim2 = 1, dimensions
					read(io_unit) (((this%io_flow_buffer(i,j,k),	i = inner_loop(1,1),inner_loop(1,2)), &
																	j = inner_loop(2,1),inner_loop(2,2)), &
																	k = inner_loop(3,1),inner_loop(3,2))
					this%vector_flow_io_fields(fields_counter)%v_ptr%pr(dim1)%cells(dim2,:,:,:)	= this%io_flow_buffer											
				end do
				end do

				close(io_unit)
			end do
		end if		
		
	end subroutine
	
	subroutine input_mesh(this)
		class(data_io)	,intent(inout)	:: this

		integer						:: dimensions
		integer		,dimension(3,2)	:: utter_loop

		integer				:: io_unit
		character(len=65)	:: file_name
		
		integer	:: i,j,k, fields_counter, dim

		dimensions	= this%domain%get_domain_dimensions()
		utter_loop	= this%domain%get_local_utter_cells_bounds()

		write(file_name,'(A,A,I3.3,A,A,A)') trim(this%data_output_folder) , trim(fold_sep) , this%load_counter, trim(fold_sep) , 'mesh' , trim(data_io_data_format)
		open(newunit = io_unit, file = file_name, status = 'old', form = 'binary')
		do dim = 1, dimensions
			read(io_unit) (((this%io_cons_buffer(i,j,k),	i = utter_loop(1,1),utter_loop(1,2)), &
															j = utter_loop(2,1),utter_loop(2,2)), &
															k = utter_loop(3,1),utter_loop(3,2))
			this%mesh%mesh_ptr%mesh(dim,:,:,:) = this%io_cons_buffer 									
		end do
		close(io_unit)
	end subroutine	
	
	subroutine input_bounds(this)
		class(data_io)	,intent(inout)	:: this
		
		integer		,dimension(3,2)	:: utter_loop

		integer				:: io_unit
		character(len=65)	:: file_name
		
		integer	:: i,j,k

		utter_loop	= this%domain%get_local_utter_cells_bounds()

		write(file_name,'(A,A,I3.3,A,A,A)') trim(this%data_output_folder) , trim(fold_sep) , this%load_counter, trim(fold_sep) , 'bc' , trim(data_io_data_format)
		open(newunit = io_unit, file = file_name, status = 'old', form = 'binary')
		read(io_unit) (((this%io_cons_buffer(i,j,k),	i = utter_loop(1,1),utter_loop(1,2)), &
														j = utter_loop(2,1),utter_loop(2,2)), &
														k = utter_loop(3,1),utter_loop(3,2))
		this%boundaries%bc_ptr%bc_markers(:,:,:) = this%io_cons_buffer
		close(io_unit)
	end subroutine	
	
	subroutine commit_mpi_io_datatype(this)
		class(data_io)	,intent(inout)	:: this

		integer					:: dimensions
		integer ,dimension(3)	:: sizes, subsizes, starts 
		integer	,dimension(3)	:: grid_coord, grid_size

		integer		:: dim, error

		this%mpi_array_bounds_type = 0

#ifdef mpi

		dimensions	= this%domain%get_domain_dimensions()

		sizes 					= this%domain%get_global_cells_number()
		subsizes 				= this%io_cons_buffer_size
		starts					= this%domain%get_global_offset()  
		starts(1:dimensions)	= starts(1:dimensions) + this%io_cons_buffer_bounds(1:dimensions,1)	

		call MPI_TYPE_CREATE_SUBARRAY(3,sizes,subsizes,starts, MPI_ORDER_FORTRAN, MPI_DOUBLE_PRECISION, this%mpi_array_bounds_type, error)
		call MPI_TYPE_COMMIT(this%mpi_array_bounds_type,error)

		sizes 		= this%domain%get_global_faces_number()
		subsizes 	= this%domain%get_local_faces_number()
		starts		= this%domain%get_global_offset()

		sizes(1:dimensions)		= sizes(1:dimensions)	 - 2
		subsizes(1:dimensions)	= subsizes(1:dimensions) - 2
		starts(1:dimensions)	= starts(1:dimensions)		

		call MPI_TYPE_CREATE_SUBARRAY(3,sizes,subsizes,starts, MPI_ORDER_FORTRAN, MPI_DOUBLE_PRECISION, this%mpi_flow_array_nobounds_type, error)
		call MPI_TYPE_COMMIT(this%mpi_flow_array_nobounds_type,error)		
#endif
	end subroutine


subroutine input_fields_mpi(this)
		class(data_io)	,intent(inout)	:: this
		
		integer					:: mpi_io_unit
		character(len=65)		:: file_name
		
		integer					:: dimensions
		integer					:: vector_projections_number
		integer	,dimension(2)	:: tensor_projections_number
		integer	,dimension(3)	:: local_cells_number, global_cells_number, local_faces_number, global_faces_number
		integer	,dimension(3,2)	:: utter_loop, utter_loop_max, inner_loop, inner_loop_max
		integer					:: buffer_size
		integer					:: i,j,k, fields_counter, dim, dim1, dim2, dim3
		integer					:: request, error, info

#ifdef mpi	
		integer(kind=MPI_OFFSET_KIND)	:: displacement
		integer 						:: status(MPI_STATUS_SIZE)
#endif

#ifdef mpi

		dimensions			= this%domain%get_domain_dimensions()

		global_cells_number	= this%domain%get_global_cells_number()

		do fields_counter = 1,size(this%scalar_io_fields)
			
			write(file_name,'(A,A,I3.3,A,A,A)') trim(this%data_output_folder) , trim(fold_sep) , this%load_counter, trim(fold_sep) , trim(this%scalar_io_fields_names(fields_counter)) , trim(data_io_data_format)		
			call MPI_FILE_OPEN(MPI_COMM_WORLD, file_name, MPI_MODE_RDONLY, MPI_INFO_NULL, mpi_io_unit, error) 

			this%io_cons_buffer	= 0.0_dkind
			do k = this%io_cons_loop_max(3,1),this%io_cons_loop_max(3,2) 
			do j = this%io_cons_loop_max(2,1),this%io_cons_loop_max(2,2) 

				displacement =  ((k-this%io_cons_loop(3,1)) * global_cells_number(2)  * global_cells_number(1) + &
								 (j-this%io_cons_loop(2,1)) * global_cells_number(1)) * sizeof(this%io_cons_buffer(1,1,1))

				call MPI_FILE_SET_VIEW(mpi_io_unit, max(displacement,0), MPI_DOUBLE_PRECISION, this%mpi_array_bounds_type, 'native', MPI_INFO_NULL, error) 

				if (((k <= this%io_cons_loop(3,2)).and.(j <= this%io_cons_loop(2,2))).and.((k >= this%io_cons_loop(3,1)).and.(j >= this%io_cons_loop(2,1)))) then
					call MPI_FILE_IREAD(mpi_io_unit, this%io_cons_buffer(this%io_cons_buffer_bounds(1,1):this%io_cons_buffer_bounds(1,2),j,k), this%io_cons_buffer_size(1), MPI_DOUBLE_PRECISION, request, error) 
					call MPI_WAIT(request,status,error)
				end if

			end do
			end do

			this%scalar_io_fields(fields_counter)%s_ptr%cells(:,:,:) = this%io_cons_buffer

			call MPI_FILE_CLOSE(mpi_io_unit,error)
		end do

		do fields_counter = 1,size(this%vector_io_fields)
			
			write(file_name,'(A,A,I3.3,A,A,A)') trim(this%data_output_folder) , trim(fold_sep) , this%load_counter, trim(fold_sep) , trim(this%vector_io_fields_names(fields_counter)) , trim(data_io_data_format)
			call MPI_FILE_OPEN(MPI_COMM_WORLD, file_name, MPI_MODE_RDONLY, MPI_INFO_NULL, mpi_io_unit, error) 

			vector_projections_number = this%vector_io_fields(fields_counter)%v_ptr%get_projections_number()

			do dim = 1, vector_projections_number
				this%io_cons_buffer = 0.0_dkind
				do k = this%io_cons_loop_max(3,1),this%io_cons_loop_max(3,2) 
				do j = this%io_cons_loop_max(2,1),this%io_cons_loop_max(2,2) 

					displacement = (dim-1) * sizeof(this%io_cons_buffer(1,1,1))
					do dim1 = 1, dimensions
						displacement = displacement * global_cells_number(dim1)
					end do

					displacement =  displacement +	((k-this%io_cons_loop(3,1)) * global_cells_number(2)  * global_cells_number(1) + &
									 				 (j-this%io_cons_loop(2,1)) * global_cells_number(1)) * sizeof(this%io_cons_buffer(1,1,1))

					call MPI_FILE_SET_VIEW(mpi_io_unit, max(displacement,0), MPI_DOUBLE_PRECISION, this%mpi_array_bounds_type, 'native', MPI_INFO_NULL, error) 

					if (((k <= this%io_cons_loop(3,2)).and.(j <= this%io_cons_loop(2,2))).and.((k >= this%io_cons_loop(3,1)).and.(j >= this%io_cons_loop(2,1)))) then
						call MPI_FILE_IREAD(mpi_io_unit, this%io_cons_buffer(this%io_cons_buffer_bounds(1,1):this%io_cons_buffer_bounds(1,2),j,k), this%io_cons_buffer_size(1), MPI_DOUBLE_PRECISION, request, error) 
						call MPI_WAIT(request,status,error)
					end if

				end do
				end do
				this%vector_io_fields(fields_counter)%v_ptr%pr(dim)%cells(:,:,:) = this%io_cons_buffer
			end do

			call MPI_FILE_CLOSE(mpi_io_unit,error)
		end do

		do fields_counter = 1,size(this%tensor_io_fields)
			write(file_name,'(A,A,I3.3,A,A,A)') trim(this%data_output_folder) , trim(fold_sep) , this%load_counter, trim(fold_sep) , trim(this%tensor_io_fields_names(fields_counter)) , trim(data_io_data_format)
			call MPI_FILE_OPEN(MPI_COMM_WORLD, file_name, MPI_MODE_RDONLY, MPI_INFO_NULL, mpi_io_unit, error) 

			tensor_projections_number = this%tensor_io_fields(fields_counter)%t_ptr%get_projections_number()

			do dim1 = 1, tensor_projections_number(1)
			do dim2 = 1, tensor_projections_number(2)
				this%io_cons_buffer = 0.0_dkind
				do k = this%io_cons_loop_max(3,1),this%io_cons_loop_max(3,2) 
				do j = this%io_cons_loop_max(2,1),this%io_cons_loop_max(2,2) 

					displacement = ((dim2-1) + (dim1-1)*tensor_projections_number(2)) * sizeof(this%io_cons_buffer(1,1,1))
					do dim3 = 1, dimensions
						displacement = displacement * global_cells_number(dim3) 
					end do

					displacement =  displacement +	((k-this%io_cons_loop(3,1)) * global_cells_number(2)  * global_cells_number(1) + &
									 				 (j-this%io_cons_loop(2,1)) * global_cells_number(1)) * sizeof(this%io_cons_buffer(1,1,1))

					call MPI_FILE_SET_VIEW(mpi_io_unit, max(displacement,0), MPI_DOUBLE_PRECISION, this%mpi_array_bounds_type, 'native', MPI_INFO_NULL, error) 

					if (((k <= this%io_cons_loop(3,2)).and.(j <= this%io_cons_loop(2,2))).and.((k >= this%io_cons_loop(3,1)).and.(j >= this%io_cons_loop(2,1)))) then
						call MPI_FILE_IREAD(mpi_io_unit, this%io_cons_buffer(this%io_cons_buffer_bounds(1,1):this%io_cons_buffer_bounds(1,2),j,k), this%io_cons_buffer_size(1), MPI_DOUBLE_PRECISION, request, error) 
						call MPI_WAIT(request,status,error)
					end if

				end do
				end do
				this%tensor_io_fields(fields_counter)%t_ptr%pr(dim1,dim2)%cells(:,:,:) = this%io_cons_buffer									
			end do
			end do
			call MPI_FILE_CLOSE(mpi_io_unit,error)
		end do

		dimensions		= this%domain%get_domain_dimensions()

		inner_loop		= this%domain%get_local_inner_faces_bounds()
		inner_loop_max	= this%domain%get_local_inner_faces_bounds_max()

		global_faces_number	= this%domain%get_global_faces_number()
		local_faces_number	= this%domain%get_local_faces_number()

		buffer_size = local_faces_number(1) - 2

		if (allocated(this%scalar_flow_io_fields)) then
			do fields_counter = 1,size(this%scalar_flow_io_fields)	
				write(file_name,'(A,A,I3.3,A,A,A)') trim(this%data_output_folder) , trim(fold_sep) , this%load_counter, trim(fold_sep) , trim(this%scalar_flow_io_fields_names(fields_counter)) , trim(data_io_data_format)		
				call MPI_FILE_OPEN(MPI_COMM_WORLD, file_name, MPI_MODE_RDONLY, MPI_INFO_NULL, mpi_io_unit, error) 

				do dim = 1,dimensions

					this%io_flow_buffer = 0.0_dkind 
					do k = inner_loop_max(3,1),inner_loop_max(3,2) 
					do j = inner_loop_max(2,1),inner_loop_max(2,2) 

						displacement = (dim-1) * sizeof(this%io_flow_buffer(1,1,1))
						do dim1 = 1, dimensions
							displacement = displacement * (global_faces_number(dim1) - 2)
						end do

						displacement =  displacement +	((k-inner_loop_max(3,1)) * (global_faces_number(2) - 2)  * (global_faces_number(1) - 2) + &
														(j-inner_loop_max(2,1)) * (global_faces_number(1) - 2)) * sizeof(this%io_flow_buffer(1,1,1))


						call MPI_FILE_SET_VIEW(mpi_io_unit, displacement, MPI_DOUBLE_PRECISION, this%mpi_flow_array_nobounds_type, 'native', MPI_INFO_NULL, error) 

						if ((k <=  inner_loop(3,2)).and.(j <= inner_loop(2,2))) then
							call MPI_FILE_IREAD(mpi_io_unit, this%io_flow_buffer(inner_loop(1,1):inner_loop(1,2),j,k), buffer_size, MPI_DOUBLE_PRECISION, request, error) 
							call MPI_WAIT(request,status,error)
						end if

					end do
					end do

					this%scalar_flow_io_fields(fields_counter)%s_ptr%cells(dim,:,:,:) = this%io_flow_buffer
				end do

				call MPI_FILE_CLOSE(mpi_io_unit,error)
			end do
		end if

		if (allocated(this%vector_flow_io_fields)) then
			do fields_counter = 1,size(this%vector_flow_io_fields)
				write(file_name,'(A,A,I3.3,A,A,A)') trim(this%data_output_folder) , trim(fold_sep) , this%load_counter, trim(fold_sep) , trim(this%vector_flow_io_fields_names(fields_counter)) , trim(data_io_data_format)
				call MPI_FILE_OPEN(MPI_COMM_WORLD, file_name, MPI_MODE_RDONLY, MPI_INFO_NULL, mpi_io_unit, error) 

				vector_projections_number = this%vector_flow_io_fields(fields_counter)%v_ptr%get_projections_number_flow()
				
				do dim1 = 1, vector_projections_number
				do dim2 = 1, dimensions
					this%io_flow_buffer = 0.0_dkind 
					do k = inner_loop_max(3,1),inner_loop_max(3,2) 
					do j = inner_loop_max(2,1),inner_loop_max(2,2) 

						displacement = ((dim2-1) + (dim1-1)*dimensions) * sizeof(this%io_flow_buffer(1,1,1))
						do dim3 = 1, dimensions
							displacement = displacement * (global_faces_number(dim3) - 2)
						end do

						displacement =  displacement +	((k-inner_loop_max(3,1)) * (global_faces_number(2) - 2)  * (global_faces_number(1) - 2) + &
														(j-inner_loop_max(2,1)) * (global_faces_number(1)-2)) * sizeof(this%io_flow_buffer(1,1,1))

						call MPI_FILE_SET_VIEW(mpi_io_unit, displacement, MPI_DOUBLE_PRECISION, this%mpi_flow_array_nobounds_type, 'native', MPI_INFO_NULL, error) 

						if ((k <=  inner_loop(3,2)).and.(j <= inner_loop(2,2))) then
							call MPI_FILE_IREAD(mpi_io_unit, this%io_flow_buffer(inner_loop(1,1):inner_loop(1,2),j,k), buffer_size, MPI_DOUBLE_PRECISION, request, error) 
							call MPI_WAIT(request,status,error)
						end if

					end do
					end do

					this%vector_flow_io_fields(fields_counter)%v_ptr%pr(dim1)%cells(dim2,:,:,:) = this%io_flow_buffer

				end do
				end do

				call MPI_FILE_CLOSE(mpi_io_unit,error)
			end do
		end if

#endif
	end subroutine


	subroutine input_mesh_mpi(this)
		class(data_io)	,intent(inout)	:: this
		
		integer				:: mpi_io_unit
		character(len=65)	:: file_name

		integer					:: dimensions
		integer	,dimension(3)	:: local_cells_number, global_cells_number
		integer	,dimension(3,2)	:: utter_loop, utter_loop_max
		integer	,dimension(3)	:: grid_coord, grid_size, buffer_bounds

		integer 				:: buffer_size
		integer					:: i,j,k, dim, dim1
		integer					:: request, error, info

#ifdef mpi	
		integer(kind=MPI_OFFSET_KIND)	:: displacement
		integer 						:: status(MPI_STATUS_SIZE)
#endif

#ifdef mpi

		dimensions		= this%domain%get_domain_dimensions()

		global_cells_number	= this%domain%get_global_cells_number()

		write(file_name,'(A,A,I3.3,A,A,A)') trim(this%data_output_folder) , trim(fold_sep) , this%load_counter, trim(fold_sep) , 'mesh' , trim(data_io_data_format)
		call MPI_FILE_OPEN(MPI_COMM_WORLD, file_name, MPI_MODE_RDONLY, MPI_INFO_NULL, mpi_io_unit, error) 

		do dim = 1, dimensions
			this%io_cons_buffer = 0.0_dkind
			do k = this%io_cons_loop_max(3,1),this%io_cons_loop_max(3,2) 
			do j = this%io_cons_loop_max(2,1),this%io_cons_loop_max(2,2) 

				displacement = (dim-1) * sizeof(this%io_cons_buffer(1,1,1))
				do dim1 = 1, dimensions
					displacement = displacement * global_cells_number(dim1) 
				end do

				displacement =  displacement +	((k-this%io_cons_loop(3,1)) * global_cells_number(2)  * global_cells_number(1) + &
												 (j-this%io_cons_loop(2,1)) * global_cells_number(1)) * sizeof(this%io_cons_buffer(1,1,1))

				call MPI_FILE_SET_VIEW(mpi_io_unit, max(displacement,0), MPI_DOUBLE_PRECISION, this%mpi_array_bounds_type, 'native', MPI_INFO_NULL, error) 

				if (((k <= this%io_cons_loop(3,2)).and.(j <= this%io_cons_loop(2,2))).and.((k >= this%io_cons_loop(3,1)).and.(j >= this%io_cons_loop(2,1)))) then
					call MPI_FILE_IREAD(mpi_io_unit, this%io_cons_buffer(this%io_cons_buffer_bounds(1,1):this%io_cons_buffer_bounds(1,2),j,k), this%io_cons_buffer_size(1), MPI_DOUBLE_PRECISION, request, error) 
					call MPI_WAIT(request,status,error)
				end if

			end do
			end do
			this%mesh%mesh_ptr%mesh(dim,:,:,:) = this%io_cons_buffer
		end do

		call MPI_FILE_CLOSE(mpi_io_unit,error)
#endif
	end subroutine

	subroutine input_bounds_mpi(this)
		class(data_io)	,intent(inout)	:: this
		
		integer				:: mpi_io_unit
		character(len=65)	:: file_name

		integer					:: dimensions
		integer	,dimension(3)	:: local_cells_number, global_cells_number
		integer	,dimension(3,2)	:: utter_loop, utter_loop_max
		integer	,dimension(3)	:: grid_coord, grid_size, buffer_bounds			

		integer :: buffer_size
		integer	:: i,j,k, dim, dim1

		integer	:: request, error, info

#ifdef mpi	
		integer(kind=MPI_OFFSET_KIND)	:: displacement
		integer 						:: status(MPI_STATUS_SIZE)
#endif

#ifdef mpi

		dimensions		= this%domain%get_domain_dimensions()

		global_cells_number	= this%domain%get_global_cells_number()

		write(file_name,'(A,A,I3.3,A,A,A)') trim(this%data_output_folder) , trim(fold_sep) , this%load_counter, trim(fold_sep) , 'bc' , trim(data_io_data_format)
		call MPI_FILE_OPEN(MPI_COMM_WORLD, file_name, MPI_MODE_RDONLY, MPI_INFO_NULL, mpi_io_unit, error) 

		do dim = 1, dimensions
			this%io_cons_buffer = 0.0_dkind
			do k = this%io_cons_loop_max(3,1),this%io_cons_loop_max(3,2) 
			do j = this%io_cons_loop_max(2,1),this%io_cons_loop_max(2,2) 

				displacement =  ((k-this%io_cons_loop(3,1)) * global_cells_number(2)  * global_cells_number(1) + &
								 (j-this%io_cons_loop(2,1)) * global_cells_number(1)) * sizeof(this%io_cons_buffer(1,1,1))

				call MPI_FILE_SET_VIEW(mpi_io_unit, max(displacement,0), MPI_DOUBLE_PRECISION, this%mpi_array_bounds_type, 'native', MPI_INFO_NULL, error) 

				if (((k <= this%io_cons_loop(3,2)).and.(j <= this%io_cons_loop(2,2))).and.((k >= this%io_cons_loop(3,1)).and.(j >= this%io_cons_loop(2,1)))) then
					call MPI_FILE_IREAD(mpi_io_unit, this%io_cons_buffer(this%io_cons_buffer_bounds(1,1):this%io_cons_buffer_bounds(1,2),j,k), this%io_cons_buffer_size(1), MPI_DOUBLE_PRECISION, request, error) 
					call MPI_WAIT(request,status,error)
				end if

			end do
			end do
			this%boundaries%bc_ptr%bc_markers(:,:,:) = this%io_cons_buffer
		end do

		call MPI_FILE_CLOSE(mpi_io_unit,error)
#endif
	end subroutine	


subroutine output_fields_mpi(this)
		class(data_io)	,intent(inout)	:: this
		
		integer					:: mpi_io_unit
		character(len=65)		:: file_name
		
		integer					:: dimensions
		integer					:: vector_projections_number
		integer	,dimension(2)	:: tensor_projections_number
		integer	,dimension(3)	:: local_cells_number, global_cells_number, local_faces_number, global_faces_number
		integer	,dimension(3,2)	:: utter_loop, utter_loop_max, inner_loop, inner_loop_max
		integer	,dimension(3)	:: grid_coord, grid_size, buffer_bounds

		integer 				:: buffer_size
		integer					:: i,j,k, fields_counter, dim, dim1, dim2, dim3
		integer					:: request, error, info

#ifdef mpi	
		integer(kind=MPI_OFFSET_KIND)	:: displacement
		integer 						:: status(MPI_STATUS_SIZE)
#endif

#ifdef mpi

		dimensions		= this%domain%get_domain_dimensions()

		global_cells_number	= this%domain%get_global_cells_number()

		do fields_counter = 1,size(this%scalar_io_fields)
			write(file_name,'(A,A,I3.3,A,A,A)') trim(this%data_output_folder) , trim(fold_sep) , this%load_counter, trim(fold_sep) , trim(this%scalar_io_fields_names(fields_counter)) , trim(data_io_data_format)		
			call MPI_FILE_OPEN(MPI_COMM_WORLD, file_name, MPI_MODE_WRONLY + MPI_MODE_CREATE, MPI_INFO_NULL, mpi_io_unit, error) 

			this%io_cons_buffer = this%scalar_io_fields(fields_counter)%s_ptr%cells(:,:,:)
			do k = this%io_cons_loop_max(3,1),this%io_cons_loop_max(3,2) 
			do j = this%io_cons_loop_max(2,1),this%io_cons_loop_max(2,2) 

				displacement =  ((k-this%io_cons_loop(3,1)) * global_cells_number(2)  * global_cells_number(1) + &
								 (j-this%io_cons_loop(2,1)) * global_cells_number(1)) * sizeof(this%io_cons_buffer(1,1,1))

				call MPI_FILE_SET_VIEW(mpi_io_unit, max(displacement,0), MPI_DOUBLE_PRECISION, this%mpi_array_bounds_type, 'native', MPI_INFO_NULL, error) 

				if (((k <= this%io_cons_loop(3,2)).and.(j <= this%io_cons_loop(2,2))).and.((k >= this%io_cons_loop(3,1)).and.(j >= this%io_cons_loop(2,1)))) then
					call MPI_FILE_IWRITE(mpi_io_unit, this%io_cons_buffer(this%io_cons_buffer_bounds(1,1):this%io_cons_buffer_bounds(1,2),j,k), this%io_cons_buffer_size(1), MPI_DOUBLE_PRECISION, request, error) 
					call MPI_WAIT(request,status,error)
				end if

			end do
			end do

			call MPI_FILE_CLOSE(mpi_io_unit,error)
		end do

		do fields_counter = 1,size(this%vector_io_fields)
			write(file_name,'(A,A,I3.3,A,A,A)') trim(this%data_output_folder) , trim(fold_sep) , this%load_counter, trim(fold_sep) , trim(this%vector_io_fields_names(fields_counter)) , trim(data_io_data_format)
			call MPI_FILE_OPEN(MPI_COMM_WORLD, file_name, MPI_MODE_WRONLY + MPI_MODE_CREATE, MPI_INFO_NULL, mpi_io_unit, error) 

			vector_projections_number = this%vector_io_fields(fields_counter)%v_ptr%get_projections_number()
			
			do dim = 1, vector_projections_number
				this%io_cons_buffer = this%vector_io_fields(fields_counter)%v_ptr%pr(dim)%cells(:,:,:)
				do k = this%io_cons_loop_max(3,1),this%io_cons_loop_max(3,2) 
				do j = this%io_cons_loop_max(2,1),this%io_cons_loop_max(2,2) 

					displacement = (dim-1) * sizeof(this%io_cons_buffer(1,1,1))
					do dim1 = 1, dimensions
						displacement = displacement * global_cells_number(dim1)
					end do

					displacement =  displacement +	((k-this%io_cons_loop(3,1)) * global_cells_number(2)  * global_cells_number(1) + &
									 				 (j-this%io_cons_loop(2,1)) * global_cells_number(1)) * sizeof(this%io_cons_buffer(1,1,1))

					call MPI_FILE_SET_VIEW(mpi_io_unit, max(displacement,0), MPI_DOUBLE_PRECISION, this%mpi_array_bounds_type, 'native', MPI_INFO_NULL, error) 

					if (((k <= this%io_cons_loop(3,2)).and.(j <= this%io_cons_loop(2,2))).and.((k >= this%io_cons_loop(3,1)).and.(j >= this%io_cons_loop(2,1)))) then
						call MPI_FILE_IWRITE(mpi_io_unit, this%io_cons_buffer(this%io_cons_buffer_bounds(1,1):this%io_cons_buffer_bounds(1,2),j,k), this%io_cons_buffer_size(1), MPI_DOUBLE_PRECISION, request, error) 
						call MPI_WAIT(request,status,error)
					end if

				end do
				end do
			end do

			call MPI_FILE_CLOSE(mpi_io_unit,error)
		end do

		do fields_counter = 1,size(this%tensor_io_fields)
			write(file_name,'(A,A,I3.3,A,A,A)') trim(this%data_output_folder) , trim(fold_sep) , this%load_counter, trim(fold_sep) , trim(this%tensor_io_fields_names(fields_counter)) , trim(data_io_data_format)
			call MPI_FILE_OPEN(MPI_COMM_WORLD, file_name, MPI_MODE_WRONLY + MPI_MODE_CREATE, MPI_INFO_NULL, mpi_io_unit, error) 

			tensor_projections_number = this%tensor_io_fields(fields_counter)%t_ptr%get_projections_number()

			do dim1 = 1, tensor_projections_number(1)
			do dim2 = 1, tensor_projections_number(2)
				this%io_cons_buffer = this%tensor_io_fields(fields_counter)%t_ptr%pr(dim1,dim2)%cells(:,:,:) 
				do k = this%io_cons_loop_max(3,1),this%io_cons_loop_max(3,2) 
				do j = this%io_cons_loop_max(2,1),this%io_cons_loop_max(2,2) 

					displacement = ((dim2-1) + (dim1-1)*tensor_projections_number(2)) * sizeof(this%io_cons_buffer(1,1,1))
					do dim3 = 1, dimensions
						displacement = displacement * global_cells_number(dim3) 
					end do

					displacement =  displacement +	((k-this%io_cons_loop(3,1)) * global_cells_number(2)  * global_cells_number(1) + &
									 				 (j-this%io_cons_loop(2,1)) * global_cells_number(1)) * sizeof(this%io_cons_buffer(1,1,1))

					call MPI_FILE_SET_VIEW(mpi_io_unit, max(displacement,0), MPI_DOUBLE_PRECISION, this%mpi_array_bounds_type, 'native', MPI_INFO_NULL, error) 

					if (((k <= this%io_cons_loop(3,2)).and.(j <= this%io_cons_loop(2,2))).and.((k >= this%io_cons_loop(3,1)).and.(j >= this%io_cons_loop(2,1)))) then
						call MPI_FILE_IWRITE(mpi_io_unit, this%io_cons_buffer(this%io_cons_buffer_bounds(1,1):this%io_cons_buffer_bounds(1,2),j,k), this%io_cons_buffer_size(1), MPI_DOUBLE_PRECISION, request, error) 
						call MPI_WAIT(request,status,error)
					end if

				end do
				end do
			end do
			end do
			call MPI_FILE_CLOSE(mpi_io_unit,error)
		end do

		dimensions		= this%domain%get_domain_dimensions()

		inner_loop		= this%domain%get_local_inner_faces_bounds()
		inner_loop_max	= this%domain%get_local_inner_faces_bounds_max()

		global_faces_number	= this%domain%get_global_faces_number()
		local_faces_number	= this%domain%get_local_faces_number()

		buffer_size = local_faces_number(1) - 2

		if (allocated(this%scalar_flow_io_fields)) then
			do fields_counter = 1,size(this%scalar_flow_io_fields)	
				write(file_name,'(A,A,I3.3,A,A,A)') trim(this%data_output_folder) , trim(fold_sep) , this%load_counter, trim(fold_sep) , trim(this%scalar_flow_io_fields_names(fields_counter)) , trim(data_io_data_format)		
				call MPI_FILE_OPEN(MPI_COMM_WORLD, file_name, MPI_MODE_WRONLY + MPI_MODE_CREATE, MPI_INFO_NULL, mpi_io_unit, error) 

				do dim = 1,dimensions

					this%io_flow_buffer = this%scalar_flow_io_fields(fields_counter)%s_ptr%cells(dim,:,:,:)
					do k = inner_loop_max(3,1),inner_loop_max(3,2) 
					do j = inner_loop_max(2,1),inner_loop_max(2,2) 

						displacement = (dim-1) * sizeof(this%io_flow_buffer(1,1,1))
						do dim1 = 1, dimensions
							displacement = displacement * (global_faces_number(dim1) - 2)
						end do

						displacement =  displacement +	((k-inner_loop_max(3,1)) * (global_faces_number(2) - 2)  * (global_faces_number(1) - 2) + &
														(j-inner_loop_max(2,1)) * (global_faces_number(1) - 2)) * sizeof(this%io_flow_buffer(1,1,1))


						call MPI_FILE_SET_VIEW(mpi_io_unit, displacement, MPI_DOUBLE_PRECISION, this%mpi_flow_array_nobounds_type, 'native', MPI_INFO_NULL, error) 

						if ((k <=  inner_loop(3,2)).and.(j <= inner_loop(2,2))) then
							call MPI_FILE_IWRITE(mpi_io_unit, this%io_flow_buffer(inner_loop(1,1):inner_loop(1,2),j,k), buffer_size, MPI_DOUBLE_PRECISION, request, error) 
							call MPI_WAIT(request,status,error)
						end if

					end do
					end do

				end do

				call MPI_FILE_CLOSE(mpi_io_unit,error)
			end do
		end if

		if (allocated(this%vector_flow_io_fields)) then
			do fields_counter = 1,size(this%vector_flow_io_fields)
				write(file_name,'(A,A,I3.3,A,A,A)') trim(this%data_output_folder) , trim(fold_sep) , this%load_counter, trim(fold_sep) , trim(this%vector_flow_io_fields_names(fields_counter)) , trim(data_io_data_format)
				call MPI_FILE_OPEN(MPI_COMM_WORLD, file_name, MPI_MODE_WRONLY + MPI_MODE_CREATE, MPI_INFO_NULL, mpi_io_unit, error) 

				vector_projections_number = this%vector_flow_io_fields(fields_counter)%v_ptr%get_projections_number_flow()
				
				do dim1 = 1, vector_projections_number
				do dim2 = 1, dimensions
					this%io_flow_buffer = this%vector_flow_io_fields(fields_counter)%v_ptr%pr(dim1)%cells(dim2,:,:,:)
					do k = inner_loop_max(3,1),inner_loop_max(3,2) 
					do j = inner_loop_max(2,1),inner_loop_max(2,2) 

						displacement = ((dim2-1) + (dim1-1)*dimensions) * sizeof(this%io_flow_buffer(1,1,1))
						do dim3 = 1, dimensions
							displacement = displacement * (global_faces_number(dim3) - 2)
						end do

						displacement =  displacement +	((k-inner_loop_max(3,1)) * (global_faces_number(2) - 2)  * (global_faces_number(1) - 2) + &
														(j-inner_loop_max(2,1)) * (global_faces_number(1)-2)) * sizeof(this%io_flow_buffer(1,1,1))

						call MPI_FILE_SET_VIEW(mpi_io_unit, displacement, MPI_DOUBLE_PRECISION, this%mpi_flow_array_nobounds_type, 'native', MPI_INFO_NULL, error) 

						if ((k <=  inner_loop(3,2)).and.(j <= inner_loop(2,2))) then
							call MPI_FILE_IWRITE(mpi_io_unit, this%io_flow_buffer(inner_loop(1,1):inner_loop(1,2),j,k), buffer_size, MPI_DOUBLE_PRECISION, request, error) 
							call MPI_WAIT(request,status,error)
						end if

					end do
					end do
				end do
				end do

				call MPI_FILE_CLOSE(mpi_io_unit,error)
			end do
		end if

#endif
	end subroutine


	subroutine output_mesh_mpi(this)
		class(data_io)	,intent(inout)	:: this
		
		integer				:: mpi_io_unit
		character(len=65)	:: file_name

		integer					:: dimensions
		integer	,dimension(3)	:: local_cells_number, global_cells_number
		integer	,dimension(3,2)	:: utter_loop, utter_loop_max
		integer	,dimension(3)	:: grid_coord, grid_size, buffer_bounds		

		integer 				:: buffer_size
		integer					:: i,j,k, dim, dim1
		integer					:: request, error, info

#ifdef mpi	
		integer(kind=MPI_OFFSET_KIND)	:: displacement
		integer 						:: status(MPI_STATUS_SIZE)
#endif

#ifdef mpi

		dimensions		= this%domain%get_domain_dimensions()

		global_cells_number	= this%domain%get_global_cells_number()

		write(file_name,'(A,A,I3.3,A,A,A)') trim(this%data_output_folder) , trim(fold_sep) , this%load_counter, trim(fold_sep) , 'mesh' , trim(data_io_data_format)
		call MPI_FILE_OPEN(MPI_COMM_WORLD, file_name, MPI_MODE_WRONLY + MPI_MODE_CREATE, MPI_INFO_NULL, mpi_io_unit, error) 

		do dim = 1, dimensions
			this%io_cons_buffer = this%mesh%mesh_ptr%mesh(dim,:,:,:)
			do k = this%io_cons_loop_max(3,1),this%io_cons_loop_max(3,2) 
			do j = this%io_cons_loop_max(2,1),this%io_cons_loop_max(2,2) 

				displacement = (dim-1) * sizeof(this%io_cons_buffer(1,1,1))
				do dim1 = 1, dimensions
					displacement = displacement * global_cells_number(dim1) 
				end do

				displacement =  displacement +	((k-this%io_cons_loop(3,1)) * global_cells_number(2)  * global_cells_number(1) + &
												 (j-this%io_cons_loop(2,1)) * global_cells_number(1)) * sizeof(this%io_cons_buffer(1,1,1))

				call MPI_FILE_SET_VIEW(mpi_io_unit, max(displacement,0), MPI_DOUBLE_PRECISION, this%mpi_array_bounds_type, 'native', MPI_INFO_NULL, error) 

				if (((k <= this%io_cons_loop(3,2)).and.(j <= this%io_cons_loop(2,2))).and.((k >= this%io_cons_loop(3,1)).and.(j >= this%io_cons_loop(2,1)))) then
					call MPI_FILE_IWRITE(mpi_io_unit, this%io_cons_buffer(this%io_cons_buffer_bounds(1,1):this%io_cons_buffer_bounds(1,2),j,k), this%io_cons_buffer_size(1), MPI_DOUBLE_PRECISION, request, error) 
					call MPI_WAIT(request,status,error)
				end if

			end do
			end do
		end do

		call MPI_FILE_CLOSE(mpi_io_unit,error)
#endif
	end subroutine

	subroutine output_bounds_mpi(this)
		class(data_io)	,intent(inout)	:: this
		
		integer				:: mpi_io_unit
		character(len=65)	:: file_name

		integer					:: dimensions
		integer	,dimension(3)	:: local_cells_number, global_cells_number
		integer	,dimension(3,2)	:: utter_loop, utter_loop_max	
		integer	,dimension(3)	:: grid_coord, grid_size, buffer_bounds	

		integer :: buffer_size
		integer	:: i,j,k, dim, dim1

		integer	:: request, error, info

#ifdef mpi	
		integer(kind=MPI_OFFSET_KIND)	:: displacement
		integer 						:: status(MPI_STATUS_SIZE)
#endif

#ifdef mpi

		dimensions		= this%domain%get_domain_dimensions()

		global_cells_number	= this%domain%get_global_cells_number()

		write(file_name,'(A,A,I3.3,A,A,A)') trim(this%data_output_folder) , trim(fold_sep) , this%load_counter, trim(fold_sep) , 'bc' , trim(data_io_data_format)
		call MPI_FILE_OPEN(MPI_COMM_WORLD, file_name, MPI_MODE_WRONLY + MPI_MODE_CREATE, MPI_INFO_NULL, mpi_io_unit, error) 

		do dim = 1, dimensions
			this%io_cons_buffer = this%boundaries%bc_ptr%bc_markers(:,:,:)
			do k = this%io_cons_loop_max(3,1),this%io_cons_loop_max(3,2) 
			do j = this%io_cons_loop_max(2,1),this%io_cons_loop_max(2,2) 

				displacement =  ((k-this%io_cons_loop(3,1)) * global_cells_number(2)  * global_cells_number(1) + &
								 (j-this%io_cons_loop(2,1)) * global_cells_number(1)) * sizeof(this%io_cons_buffer(1,1,1))

				call MPI_FILE_SET_VIEW(mpi_io_unit, max(displacement,0), MPI_DOUBLE_PRECISION, this%mpi_array_bounds_type, 'native', MPI_INFO_NULL, error) 

				if (((k <= this%io_cons_loop(3,2)).and.(j <= this%io_cons_loop(2,2))).and.((k >= this%io_cons_loop(3,1)).and.(j >= this%io_cons_loop(2,1)))) then
					call MPI_FILE_IWRITE(mpi_io_unit, this%io_cons_buffer(this%io_cons_buffer_bounds(1,1):this%io_cons_buffer_bounds(1,2),j,k), this%io_cons_buffer_size(1), MPI_DOUBLE_PRECISION, request, error) 
					call MPI_WAIT(request,status,error)
				end if

			end do
			end do
		end do

		call MPI_FILE_CLOSE(mpi_io_unit,error)
#endif
	end subroutine

	pure integer function get_load_counter(this)
		class(data_io)					,intent(in)	:: this

		get_load_counter = this%load_counter
	end function

	pure function get_data_output_folder(this)
		class(data_io)					,intent(in)	:: this
		character(len=100)				:: get_data_output_folder

		get_data_output_folder = this%data_output_folder
	end function

	subroutine add_io_scalar_cons_field(this,scal_ptr)
		class(data_io)					,intent(inout)	:: this
		type(field_scalar_cons_pointer)	,intent(in)		:: scal_ptr

		type(field_scalar_cons_pointer)	,dimension(:)	,allocatable	:: scalar_cons_io_fields_buffer
		character(len=35)				,dimension(:)	,allocatable	:: scalar_cons_io_fields_names_buffer

		if(allocated(this%scalar_io_fields)) then	
			allocate(scalar_cons_io_fields_buffer(size(this%scalar_io_fields)+1))
			allocate(scalar_cons_io_fields_names_buffer(size(this%scalar_io_fields)+1))
			scalar_cons_io_fields_buffer(:size(this%scalar_io_fields))		= this%scalar_io_fields
			scalar_cons_io_fields_buffer(size(this%scalar_io_fields)+1)	= scal_ptr
			
			scalar_cons_io_fields_names_buffer(:size(this%scalar_io_fields))	= this%scalar_io_fields_names
			scalar_cons_io_fields_names_buffer(size(this%scalar_io_fields)+1)	= scal_ptr%s_ptr%name_long
			deallocate(this%scalar_io_fields)
			deallocate(this%scalar_io_fields_names)
		else
			allocate(scalar_cons_io_fields_buffer(1))
			allocate(scalar_cons_io_fields_names_buffer(1))
			scalar_cons_io_fields_buffer(1)	= scal_ptr
			scalar_cons_io_fields_names_buffer(1)	= scal_ptr%s_ptr%name_long
		end if

		allocate(this%scalar_io_fields(size(scalar_cons_io_fields_buffer)))
		this%scalar_io_fields = scalar_cons_io_fields_buffer

		allocate(this%scalar_io_fields_names(size(scalar_cons_io_fields_buffer)))
		this%scalar_io_fields_names = scalar_cons_io_fields_names_buffer
	end subroutine	
	
	subroutine add_io_scalar_flow_field(this,scal_ptr)
		class(data_io)					,intent(inout)	:: this
		type(field_scalar_flow_pointer)	,intent(in)		:: scal_ptr

		type(field_scalar_flow_pointer)	,dimension(:)	,allocatable	:: scalar_flow_io_fields_buffer
		character(len=35)				,dimension(:)	,allocatable	:: scalar_flow_io_fields_names_buffer

		if(allocated(this%scalar_flow_io_fields)) then	
			allocate(scalar_flow_io_fields_buffer(size(this%scalar_flow_io_fields)+1))
			allocate(scalar_flow_io_fields_names_buffer(size(this%scalar_flow_io_fields)+1))
			scalar_flow_io_fields_buffer(:size(this%scalar_flow_io_fields))		= this%scalar_flow_io_fields
			scalar_flow_io_fields_buffer(size(this%scalar_flow_io_fields)+1)	= scal_ptr
			
			scalar_flow_io_fields_names_buffer(:size(this%scalar_flow_io_fields))	= this%scalar_flow_io_fields_names
			scalar_flow_io_fields_names_buffer(size(this%scalar_flow_io_fields)+1)	= scal_ptr%s_ptr%name_long
			deallocate(this%scalar_flow_io_fields)
			deallocate(this%scalar_flow_io_fields_names)
		else
			allocate(scalar_flow_io_fields_buffer(1))
			allocate(scalar_flow_io_fields_names_buffer(1))
			scalar_flow_io_fields_buffer(1)	= scal_ptr
			scalar_flow_io_fields_names_buffer(1)	= scal_ptr%s_ptr%name_long
		end if

		allocate(this%scalar_flow_io_fields(size(scalar_flow_io_fields_buffer)))
		this%scalar_flow_io_fields = scalar_flow_io_fields_buffer

		allocate(this%scalar_flow_io_fields_names(size(scalar_flow_io_fields_buffer)))
		this%scalar_flow_io_fields_names = scalar_flow_io_fields_names_buffer
	end subroutine

	subroutine add_io_vector_flow_field(this,vect_ptr)
		class(data_io)					,intent(inout)	:: this
		type(field_vector_flow_pointer)	,intent(in)		:: vect_ptr

		type(field_vector_flow_pointer)	,dimension(:)	,allocatable	:: vector_flow_io_fields_buffer
		character(len=35)				,dimension(:)	,allocatable	:: vector_flow_io_fields_names_buffer

		if(allocated(this%vector_flow_io_fields)) then	
			allocate(vector_flow_io_fields_buffer(size(this%vector_flow_io_fields)+1))
			allocate(vector_flow_io_fields_names_buffer(size(this%vector_flow_io_fields)+1))
			vector_flow_io_fields_buffer(:size(this%vector_flow_io_fields))		= this%vector_flow_io_fields
			vector_flow_io_fields_buffer(size(this%vector_flow_io_fields)+1)	= vect_ptr
			
			vector_flow_io_fields_names_buffer(:size(this%vector_flow_io_fields))	= this%vector_flow_io_fields_names
			vector_flow_io_fields_names_buffer(size(this%vector_flow_io_fields)+1)	= vect_ptr%v_ptr%name_long
			deallocate(this%vector_flow_io_fields)
			deallocate(this%vector_flow_io_fields_names)
		else
			allocate(vector_flow_io_fields_buffer(1))
			allocate(vector_flow_io_fields_names_buffer(1))
			vector_flow_io_fields_buffer(1)	= vect_ptr
			vector_flow_io_fields_names_buffer(1)	= vect_ptr%v_ptr%name_long
		end if

		allocate(this%vector_flow_io_fields(size(vector_flow_io_fields_buffer)))
		this%vector_flow_io_fields = vector_flow_io_fields_buffer

		allocate(this%vector_flow_io_fields_names(size(vector_flow_io_fields_buffer)))
		this%vector_flow_io_fields_names = vector_flow_io_fields_names_buffer
	end subroutine
	

end module
