module post_processor_operation_class

	use kind_parameters
	use computational_domain_class	
	use data_manager_class
	use field_pointers
	use boundary_conditions_class

	implicit none

	private
	public	:: post_processor_operation, post_processor_operation_c

	type ::post_processor_operation
		private
		type(field_scalar_cons_pointer)	,public	:: operation_field_scal
		type(field_vector_cons_pointer)	,public	:: operation_field_vect
		
		type(computational_domain)			:: domain
		
		character(len=40)			:: field_name
		character(len=20)			:: operation_type
		integer						:: grad_projection
		integer	,dimension(3,2)		:: operation_area
		integer	,dimension(3)		:: operation_area_center
		integer	,dimension(3)		:: operation_area_distance
	contains
		! IO Getters/Setters
		procedure	,private	:: read_properties
		procedure	,private	:: write_properties
		procedure	,private	:: set_properties	

		procedure	:: set_point
		procedure	:: process_operation
		procedure	:: get_values_number

		! Logger
		procedure	:: write_log
		
	end type

	interface post_processor_operation_c
		module procedure	constructor
		module procedure	constructor_file
	end interface
contains

	type(post_processor_operation)	function constructor(manager,field_name,operation_type,operation_area,operation_area_distance,grad_projection,post_processor_data_file_unit)
		type(data_manager)				,intent(in)	:: manager
		character(len=*)				,intent(in)	:: field_name		
		character(len=*)				,intent(in)	:: operation_type
		integer	,dimension(3,2)			,intent(in)	:: operation_area
		integer	,dimension(3)			,intent(in)	:: operation_area_distance
		integer							,intent(in)	:: grad_projection	
		integer							,intent(in)	:: post_processor_data_file_unit

		call constructor%set_properties(manager,field_name,operation_type,operation_area,operation_area_distance,grad_projection)
		call constructor%write_properties(post_processor_data_file_unit)
		
	end function

	type(post_processor_operation) function constructor_file(manager,post_processor_data_file_unit)
		type(data_manager)				,intent(in)	:: manager
		integer							,intent(in)	:: post_processor_data_file_unit
			
		call constructor_file%read_properties(manager,post_processor_data_file_unit)

	end function
	
	subroutine read_properties(this,manager,post_processor_data_file_unit)
		class(post_processor_operation)	,intent(inout)	:: this
		type(data_manager)				,intent(in)	:: manager
		integer							,intent(in)	:: post_processor_data_file_unit 
		character(len=40)		:: field_name		
		character(len=20)		:: operation_type
		integer	,dimension(3,2)	:: operation_area
		integer	,dimension(3)	:: operation_area_distance
		integer					:: grad_projection
		
		namelist /pproc_operation/ field_name, operation_type, operation_area, operation_area_distance, grad_projection

		read(unit = post_processor_data_file_unit, nml = pproc_operation)
		
		call this%set_properties(manager, field_name, operation_type, operation_area, operation_area_distance, grad_projection)
		
	end subroutine
	
	subroutine write_properties(this,post_processor_data_file_unit)
		class(post_processor_operation)	,intent(in)	:: this
		integer							,intent(in)	:: post_processor_data_file_unit 
		
		character(len=40)		:: field_name		
		character(len=20)		:: operation_type
		integer	,dimension(3,2)	:: operation_area
		integer	,dimension(3)	:: operation_area_distance
		integer					:: grad_projection	
		
		namelist /pproc_operation/ field_name, operation_type, operation_area, operation_area_distance, grad_projection
		
		field_name					= "'" // trim(this%field_name) // "'"
		operation_type				= this%operation_type
		operation_area				= this%operation_area
		operation_area_distance		= this%operation_area_distance
		grad_projection				= this%grad_projection
		
		write(unit = post_processor_data_file_unit, nml = pproc_operation)
	end subroutine
	
	subroutine set_properties(this,manager,field_name,operation_type,operation_area,operation_area_distance,grad_projection)
		class(post_processor_operation)	,intent(inout)	:: this
		type(data_manager)				,intent(in)	:: manager
		character(len=*)				,intent(in)	:: field_name		
		character(len=*)				,intent(in)	:: operation_type
		integer	,dimension(3,2)			,intent(in)	:: operation_area
		integer	,dimension(3)			,intent(in)	:: operation_area_distance
		integer							,intent(in)	:: grad_projection		
		
		integer	:: dimensions

		integer	,dimension(3,2)	:: area
		integer	,dimension(3)	:: center		
		integer	,dimension(3)	:: distance

		
		integer	:: dim
		
		type(field_scalar_cons_pointer)	:: scal_ptr
		type(field_vector_cons_pointer)	:: vect_ptr
		type(field_tensor_cons_pointer)	:: tens_ptr				
		
		dimensions = manager%domain%get_domain_dimensions()

		call manager%get_cons_field_pointer_by_name(scal_ptr,vect_ptr,tens_ptr,field_name)
			
		if(associated(scal_ptr%s_ptr)) then
			this%operation_field_scal%s_ptr => scal_ptr%s_ptr
			this%operation_field_vect%v_ptr => NULL()
		end if
		if(associated(vect_ptr%v_ptr)) then
			this%operation_field_scal%s_ptr => NULL()
			this%operation_field_vect%v_ptr => vect_ptr%v_ptr			
		end if

		area		= operation_area
		distance	= operation_area_distance
		center		= 0.5*(operation_area(:,2) - abs(operation_area(:,1)))

		area(dimensions+1:3,:)		= 1
		center(dimensions+1:3)		= 1
		distance(dimensions+1:3)	= 0
		
		this%field_name					= field_name
		this%operation_type				= operation_type
		this%operation_area				= area
		this%operation_area_distance	= distance
		this%operation_area_center		= center
		this%grad_projection			= grad_projection

		this%domain					= manager%domain		
	end subroutine
	
	subroutine write_log(this,log_unit)
		class(post_processor_operation)	,intent(in)	:: this
		integer							,intent(in)	:: log_unit
		
		write(log_unit,'(A,A)') 'Operation field name: ', this%field_name
		write(log_unit,'(A,A)') 'Operation type: ', this%operation_type
		write(log_unit,'(A,I2)') 'Gradient projection (if type is gradient): ', this%grad_projection
		write(log_unit,'(A,6I8)') 'Operational area: ', this%operation_area
		write(log_unit,'(A,3I8)') 'Operational area center distance from the main point: ', this%operation_area_distance
	end subroutine
	
	subroutine process_operation(this,bc_ptr,point_indexes,value)

#ifdef mpi
	use MPI
#endif	
	
		class(post_processor_operation)		,intent(in)		:: this
		type(boundary_conditions_pointer)	,intent(in)		:: bc_ptr
		real(dkind)							,intent(out)	:: value
		integer		,dimension(3)			,intent(out)	:: point_indexes
		
		
		integer		,dimension(:)	,allocatable	,save	:: point_indexes_array
		real(dkind)	,dimension(:)	,allocatable	,save	:: values_array
		
		real(dkind)					:: interm_value(1)
		integer						:: interm_point_indexes(3)
		
		real(dkind)					:: min_value, max_value
		
		integer		,dimension(3)	:: transducer_cell
		integer		,dimension(3,2)	:: global_inner_cells_bounds
		integer		,dimension(3,2)	:: local_operation_area
		
		integer						:: dimensions
		integer						:: processor_rank, processor_number, mpi_communicator
		
		integer						:: proc
		logical						:: flag
		integer						:: error
		
		processor_rank				= this%domain%get_processor_rank()
		mpi_communicator			= this%domain%get_mpi_communicator()	
		processor_number			= this%domain%get_mpi_communicator_size()
		global_inner_cells_bounds	= this%domain%get_global_inner_cells_bounds()

		if (.not.allocated(point_indexes_array)) then
			allocate(point_indexes_array(3*processor_number))
			allocate(values_array(processor_number))
		end if
		
		if (associated(this%operation_field_scal%s_ptr)) then
			select case(this%operation_type)
				case('max')
					call this%operation_field_scal%s_ptr%get_field_max(bc_ptr,this%operation_area,value,point_indexes)
				case('min')
					call this%operation_field_scal%s_ptr%get_field_min(bc_ptr,this%operation_area,value,point_indexes)
				case('transducer')
			
					transducer_cell = (/max(this%operation_area(1,1),global_inner_cells_bounds(1,1))	, &
										max(this%operation_area(2,1),global_inner_cells_bounds(2,1))	, &
										max(this%operation_area(3,1),global_inner_cells_bounds(3,1))/)
					
					transducer_cell = (/min(transducer_cell(1),global_inner_cells_bounds(1,2))	, &
										min(transducer_cell(2),global_inner_cells_bounds(2,2))	, &
										min(transducer_cell(3),global_inner_cells_bounds(3,2))/)
					
					if ( this%domain%cell_is_inside_local_domain(transducer_cell)) then
						interm_point_indexes = transducer_cell				
						call this%domain%convert_global_cell_to_local(transducer_cell)
						interm_value(1) = this%operation_field_scal%s_ptr%cells(transducer_cell(1),transducer_cell(2),transducer_cell(3))
					else
						interm_point_indexes	= (/0,0,0/)
						interm_value(1)			= 0.0_dkind
					end if
				case('max_grad')
					call this%operation_field_scal%s_ptr%get_field_max_grad(bc_ptr,this%operation_area,interm_value(1),this%grad_projection,interm_point_indexes)
					!if (this%domain%area_intersects_local_domain(this%operation_area)) then
					!	local_operation_area  = this%operation_area
					!	call this%domain%convert_global_cell_to_local(local_operation_area(:,1))
					!	call this%domain%convert_global_cell_to_local(local_operation_area(:,2))
					!	call this%operation_field_scal%s_ptr%get_field_max_grad(bc_ptr,local_operation_area,interm_value(1),this%grad_projection,interm_point_indexes)
					!else
					!	interm_point_indexes	= (/0,0,0/)
					!	interm_value(1)			= 0.0_dkind
					!end if
				case('min_grad')
					call this%operation_field_scal%s_ptr%get_field_min_grad(bc_ptr,this%operation_area,interm_value(1),this%grad_projection,interm_point_indexes)
				case('sum')
					if (this%domain%area_intersects_local_domain(this%operation_area)) then
						local_operation_area  = this%operation_area
						call this%domain%convert_global_cell_to_local(local_operation_area(:,1))
						call this%domain%convert_global_cell_to_local(local_operation_area(:,2))
						call this%operation_field_scal%s_ptr%get_field_sum(bc_ptr,local_operation_area,interm_value(1))
						interm_point_indexes	= (/0,0,0/)
					else
						interm_point_indexes	= (/0,0,0/)
						interm_value(1)			= 0.0_dkind
					end if					
				case('mean')
					call this%operation_field_scal%s_ptr%get_field_mean(bc_ptr,this%operation_area,interm_value(1))
				case('custom_operation_2')
			end select
		end if

		if (associated(this%operation_field_vect%v_ptr)) then
			select case(this%operation_type)
				case('max')
					!call this%operation_field_vect%v_ptr%pr(this%field_projection)%get_field_max(bc_ptr,this%operation_area,value,point_indexes)
				case('min')
					!call this%operation_field_vect%v_ptr%pr(this%field_projection)%get_field_min(bc_ptr,this%operation_area,value,point_indexes)
				case('transducer')
					!value = sqrt(sum(this%operation_field_vect%v_ptr%pr(:)%cells(this%anchor_point(1),this%anchor_point(2),this%anchor_point(3))))
				case('sum')
					!call this%operation_field_vect%v_ptr%pr(this%field_projection)%get_field_sum(bc_ptr,this%operation_area,value)
				case('mean_square')
					!call this%operation_field_vect%v_ptr%get_field_mean_square(bc_ptr,this%operation_area,value)
				case('custom_operation_2')
			end select
		end if
		
		values_array(processor_rank+1)								= interm_value(1)
		point_indexes_array(3*processor_rank+1:3*processor_rank+3)	= interm_point_indexes

#ifdef mpi					
		call mpi_gather(interm_value,1,MPI_DOUBLE_PRECISION,values_array,1,MPI_DOUBLE_PRECISION,0,mpi_communicator,error)
		call mpi_gather(interm_point_indexes,3,MPI_INTEGER,point_indexes_array,3,MPI_INTEGER,0,mpi_communicator,error)
#endif		

		if (processor_rank == 0) then
			select case(this%operation_type)
				case('transducer')
					do proc = 0, processor_number - 1
						if (.not.((point_indexes_array(3*proc+1) == 0).and.(point_indexes_array(3*proc+2) == 0).and.(point_indexes_array(3*proc+3) == 0))) then
							point_indexes	= point_indexes_array(3*proc+1:3*proc+3)
							value			= values_array(proc+1)
						end if
					end do
				case('min_grad')
					do proc = 0, processor_number - 1
						if (.not.((point_indexes_array(3*proc+1) == 0).and.(point_indexes_array(3*proc+2) == 0).and.(point_indexes_array(3*proc+3) == 0))) then
							min_value = values_array(proc+1)
							point_indexes	= point_indexes_array(3*proc+1:3*proc+3)
							exit
						end if
					end do
					do proc = 0, processor_number - 1
						if (.not.((point_indexes_array(3*proc+1) == 0).and.(point_indexes_array(3*proc+2) == 0).and.(point_indexes_array(3*proc+3) == 0))) then
							if (values_array(proc+1) <= min_value) then
								point_indexes	= point_indexes_array(3*proc+1:3*proc+3)
								value			= values_array(proc+1)
								min_value		= value
							end if
						end if
					end do	
				case('max_grad')
					do proc = 0, processor_number - 1
						if (.not.((point_indexes_array(3*proc+1) == 0).and.(point_indexes_array(3*proc+2) == 0).and.(point_indexes_array(3*proc+3) == 0))) then
							max_value = values_array(proc+1)
							point_indexes	= point_indexes_array(3*proc+1:3*proc+3)
							exit
						end if
					end do
					do proc = 0, processor_number - 1
						if (.not.((point_indexes_array(3*proc+1) == 0).and.(point_indexes_array(3*proc+2) == 0).and.(point_indexes_array(3*proc+3) == 0))) then
							if (values_array(proc+1) >= max_value) then
								point_indexes	= point_indexes_array(3*proc+1:3*proc+3)
								value			= values_array(proc+1)
								max_value		= value
							end if
						end if
					end do						
					
				case('sum')
					value = 0.0_dkind
					do proc = 0, processor_number - 1
						value			= value + values_array(proc+1)
						point_indexes	= (/0,0,0/)
					end do			
			end select
		else
			point_indexes	= (/0,0,0/)
			value			= 0.0_dkind	
		end if

	end subroutine
	
	integer	function get_values_number(this)
		class(post_processor_operation)		,intent(in)		:: this

		get_values_number = 1
	end function
	
	!subroutine write_info(this,unit)
	!	class(post_processor_operation)	,intent(in)	:: this
	!	integer						,intent(in)	:: unit
	!	integer		:: proj_counter
	!	
	!	if(associated(this%operation_field_scal%s_ptr)) then
	!			write(unit,*) trim(this%operation_type) // '(' // trim(this%operation_field_scal%s_ptr%name_short) // ')'
	!	end if
	!	if(associated(this%operation_field_vect%v_ptr)) then
	!		if (allocated(this%projections)) then			
	!			do proj_counter = 1, size(this%projections)
	!				write(unit,*) trim(this%operation_type) // '(' // trim(this%operation_field_vect%v_ptr%pr(this%projections(proj_counter))%name_short) // ')'
	!			end do
	!		else
	!			write(unit,*) trim(this%operation_type) // '(' // trim(this%operation_field_vect%v_ptr%name_short) // ')'
	!		end if
	!	end if
	!	write(unit,'(A,6I4)')	'operation area : ', this%operation_area
	!	write(unit,'(A,3I4)')	'operation point : ', this%point
	!	if(index(this%operation_type,'grad') /= 0) then	
	!		write(unit,'(A,I2)') 'operation axis : ', this%grad_axis
	!	end if
	!	
	!end subroutine
	
	subroutine set_point(this,point_indexes)
		class(post_processor_operation)	,intent(inout)	:: this
		integer		,dimension(3)		,intent(in)		:: point_indexes
		integer		,dimension(3)						:: new_center
		integer		,dimension(3)						:: shift
		
		new_center					= point_indexes + this%operation_area_distance
		shift						= new_center	- this%operation_area_center
		this%operation_area_center	= new_center
		this%operation_area(:,1)	= this%operation_area(:,1) + shift
		this%operation_area(:,2)	= this%operation_area(:,2) + shift
	end subroutine
	
end module
