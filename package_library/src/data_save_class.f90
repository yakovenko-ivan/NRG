module data_save_class

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
	public	::	data_save, data_save_c	
	
	type data_save
		private
		type(computational_domain)					:: domain
		type(computational_mesh_pointer)			:: mesh
		type(boundary_conditions_pointer)			:: boundaries
		
		type(field_scalar_cons_pointer)	,dimension(:)		,allocatable	:: visible_fields
		character(len=40)				,dimension(:)		,allocatable	:: visible_fields_names

		real(rkind)						,dimension(:,:,:)	,allocatable	:: io_buffer

		real(dkind)             :: save_time
		character(len=20)       :: save_time_units
		real(rkind)             :: save_time_coefficient
		character(len=25)       :: save_format
		character(len=10)       :: save_time_units_abbreviation
		
		integer					:: output_counter
		
		character(len=20)		:: data_save_folder
		
		logical					:: debug_flag

		integer					:: mpi_save_array_type		
	contains
		! IO Getters/Setters
		procedure	,private	:: read_properties
		procedure	,private	:: write_properties
		procedure	,private	:: set_properties	
		procedure	,private	:: write_output_counter		
			
		procedure				:: save_all_data
		procedure	,private	:: generate_header

		! Sequential data save
		procedure	,private	:: save_mesh
		procedure	,private	:: save_fields
		procedure	,private	:: save_bounds

		! Parallel data save
		procedure	,private	:: save_mesh_mpi
		procedure	,private	:: save_fields_mpi
		procedure	,private	:: save_bounds_mpi	

		procedure	,private	:: commit_mpi_save_datatype

		! Getters
		procedure	:: get_data_save_folder
		procedure	:: get_save_time
		procedure	:: get_save_time_units_abbreviation
		procedure	:: get_save_time_coefficient
		procedure	:: get_visible_fields_number

		!Setters 
		procedure	:: set_save_time
		
		! Logger
		procedure	:: write_log	
	end type

	interface data_save_c
		module procedure constructor
		module procedure constructor_file
	end interface

contains

	type(data_save)	function constructor(manager,visible_fields_names,save_time,save_time_units,save_format,data_save_folder,debug_flag)
		type(data_manager)					,intent(in)	:: manager
		character(len=40)	,dimension(:)	,intent(in)	:: visible_fields_names
		real(dkind)							,intent(in)	:: save_time
		character(len=*)					,intent(in)	:: save_time_units
		character(len=*)					,intent(in)	:: save_format
		character(len=*)					,intent(in)	:: data_save_folder
		logical								,intent(in)	:: debug_flag
		
		integer	:: output_counter = 0
		integer	:: io_unit

		call constructor%set_properties(manager,visible_fields_names,save_time,save_time_units,save_format,data_save_folder,debug_flag,output_counter)
		
		open(newunit = io_unit, file = data_save_data_file_name, status = 'replace', form = 'formatted', delim = 'quote')
		call constructor%write_properties(io_unit)
		close(io_unit)
				
	end function

	type(data_save)	function constructor_file(manager)

		type(data_manager)	,intent(in)	:: manager

		integer	:: io_unit
		
		open(newunit = io_unit, file = data_save_data_file_name, status = 'old', form = 'formatted', delim = 'quote')
		call constructor_file%read_properties(manager,io_unit)
		close(io_unit)
	end function

	subroutine read_properties(this,manager,data_save_file_unit)
		class(data_save)	,intent(inout)	:: this
		type(data_manager)	,intent(in)		:: manager
		integer				,intent(in)		:: data_save_file_unit

		character(len=40)	,dimension(:)	,allocatable	:: visible_fields_names
		
		integer				:: visible_fields_number
		real(dkind)			:: save_time
		character(len=20)	:: save_time_units
		character(len=25)	:: save_format
		character(len=20)	:: data_save_folder
		logical				:: debug_flag
		integer				:: output_counter
		
		integer				:: array_type

		integer				:: ierr

		namelist /data_save_parameters/ visible_fields_number, save_time, save_time_units, save_format, data_save_folder, debug_flag
		namelist /visible_fields/ visible_fields_names
		namelist /data_save_output_counter/ output_counter
		
		read(unit = data_save_file_unit, nml = data_save_parameters)
		
		allocate(visible_fields_names(visible_fields_number))
		
		read(unit = data_save_file_unit, nml = visible_fields)
		read(unit = data_save_file_unit, nml = data_save_output_counter, iostat = ierr)
		
		if (ierr /= 0 ) output_counter = 0
		call this%set_properties(manager,visible_fields_names,save_time,save_time_units,save_format,data_save_folder,debug_flag,output_counter)
	
	end subroutine
	
	subroutine write_properties(this,data_save_file_unit)
		class(data_save)	,intent(in)	:: this
		integer				,intent(in)	:: data_save_file_unit
		
		character(len=40)	,dimension(:)	,allocatable	:: visible_fields_names
		
		integer				:: visible_fields_number
		real(dkind)			:: save_time
		character(len=20)	:: save_time_units
		character(len=25)	:: save_format
		character(len=20)	:: data_save_folder
		logical				:: debug_flag
		
		namelist /data_save_parameters/ visible_fields_number, save_time, save_time_units, save_format, data_save_folder, debug_flag
		namelist /visible_fields/ visible_fields_names
	
		visible_fields_number = size(this%visible_fields_names)
		allocate(visible_fields_names,source = this%visible_fields_names)
		
		save_time					= this%save_time
		save_time_units				= this%save_time_units
		save_format					= this%save_format
		data_save_folder			= this%data_save_folder
		debug_flag					= this%debug_flag
		visible_fields_names		= this%visible_fields_names
		
		write(unit = data_save_file_unit, nml = data_save_parameters)
		write(unit = data_save_file_unit, nml = visible_fields)

	end subroutine
	
	subroutine set_properties(this,manager,visible_fields_names,save_time,save_time_units,save_format,data_save_folder,debug_flag, output_counter)
		class(data_save)					,intent(inout)	:: this
		type(data_manager)					,intent(in)	:: manager
		character(len=*)	,dimension(:)	,intent(in)	:: visible_fields_names
		real(dkind)							,intent(in)	:: save_time
		character(len=*)					,intent(in)	:: save_time_units
		character(len=*)					,intent(in)	:: save_format
		character(len=*)					,intent(in)	:: data_save_folder
		logical								,intent(in)	:: debug_flag	
		integer								,intent(in)	:: output_counter

		integer	:: vector_projections_number
		integer	:: number_of_cons_scalar_fields, number_of_cons_vector_fields
		integer	:: processor_rank

		type(field_scalar_cons_pointer)	:: scal_ptr
		type(field_vector_cons_pointer)	:: vect_ptr
		type(field_tensor_cons_pointer)	:: tens_ptr

		character(len=100)					:: system_command

		integer	,dimension(3,2)	:: allocation_bounds

		integer	:: number_of_visible_fields, list_count, scal_count, vect_count, proj_count, vis_s_count, vis_v_count
		
		this%domain				= manager%domain
		this%save_time			= save_time
		this%save_time_units	= save_time_units
		this%save_format		= save_format
		this%data_save_folder	= data_save_folder
		this%output_counter		= output_counter
		
		number_of_cons_scalar_fields = manager%get_number_of_cons_scalar_fields()
		number_of_cons_vector_fields = manager%get_number_of_cons_vector_fields()

		allocation_bounds		= this%domain%get_local_utter_cells_bounds()

		allocate(this%io_buffer(	allocation_bounds(1,1):allocation_bounds(1,2), &
									allocation_bounds(2,1):allocation_bounds(2,2), &
									allocation_bounds(3,1):allocation_bounds(3,2)))

		allocate(this%visible_fields_names(size(visible_fields_names)))
		
		do list_count = 1, size(visible_fields_names)
			this%visible_fields_names(list_count) = visible_fields_names(list_count)
		end do

		select case(this%save_time_units)
			case('minutes')
				this%save_time_units_abbreviation = 'm'
				this%save_time_coefficient        = 1.0_dkind/60.0_dkind			
			case('seconds')
				this%save_time_units_abbreviation = 's'
				this%save_time_coefficient        = 1.0_dkind		
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

		this%mesh 		= manager%computational_mesh_pointer
		this%boundaries	= manager%boundary_conditions_pointer
		
		this%debug_flag	= debug_flag

		number_of_visible_fields = 0
		
		if(this%debug_flag) then
			do scal_count = 1,number_of_cons_scalar_fields
				number_of_visible_fields = number_of_visible_fields + 1
			end do
			do vect_count = 1,number_of_cons_vector_fields
				call manager%get_cons_field_pointer_by_number(scal_ptr,vect_ptr,tens_ptr,'vector',vect_count)
				vector_projections_number = vect_ptr%v_ptr%get_projections_number()
				do proj_count = 1,vector_projections_number
					number_of_visible_fields = number_of_visible_fields + 1
				end do
			end do
		else
			do list_count = 1,size(visible_fields_names)
				call manager%get_cons_field_pointer_by_name(scal_ptr,vect_ptr,tens_ptr,this%visible_fields_names(list_count))
				if(associated(scal_ptr%s_ptr)) then
					number_of_visible_fields = number_of_visible_fields + 1
				end if
				if(associated(vect_ptr%v_ptr)) then
					vector_projections_number = vect_ptr%v_ptr%get_projections_number()
					do proj_count = 1, vector_projections_number
						number_of_visible_fields = number_of_visible_fields + 1
					end do
				end if
			end do
		end if

		allocate(this%visible_fields(number_of_visible_fields))
		number_of_visible_fields = 0
		
		if(this%debug_flag) then
			do scal_count = 1,number_of_cons_scalar_fields
				number_of_visible_fields = number_of_visible_fields + 1
				call manager%get_cons_field_pointer_by_number(scal_ptr,vect_ptr,tens_ptr,'scalar',scal_count)
				this%visible_fields(number_of_visible_fields)%s_ptr => scal_ptr%s_ptr
			end do
			do vect_count = 1,number_of_cons_vector_fields
				call manager%get_cons_field_pointer_by_number(scal_ptr,vect_ptr,tens_ptr,'vector',vect_count)
				vector_projections_number = vect_ptr%v_ptr%get_projections_number()
				do proj_count = 1,vector_projections_number
					number_of_visible_fields = number_of_visible_fields + 1
					this%visible_fields(number_of_visible_fields)%s_ptr => vect_ptr%v_ptr%pr(proj_count)
				end do
			end do
		else
			do list_count = 1,size(visible_fields_names)
				call manager%get_cons_field_pointer_by_name(scal_ptr,vect_ptr,tens_ptr,visible_fields_names(list_count))
				if(associated(scal_ptr%s_ptr)) then
					number_of_visible_fields = number_of_visible_fields + 1
					this%visible_fields(number_of_visible_fields)%s_ptr => scal_ptr%s_ptr
				end if
				if(associated(vect_ptr%v_ptr)) then
					vector_projections_number = vect_ptr%v_ptr%get_projections_number()
					do proj_count = 1, vector_projections_number
						number_of_visible_fields = number_of_visible_fields + 1
						this%visible_fields(number_of_visible_fields)%s_ptr => vect_ptr%v_ptr%pr(proj_count)
					end do
				end if
			end do
		end if

		processor_rank = this%domain%get_processor_rank()
		
		if (processor_rank == 0) then
			system_command = 'mkdir ' // trim(data_save_folder)
			call system(system_command)
		end if

	end subroutine
	
	pure function get_data_save_folder(this)
		class(data_save)	,intent(in)	:: this
		character(len=100)				:: get_data_save_folder

		get_data_save_folder = this%data_save_folder
	end function

	pure function get_save_time(this)
		class(data_save)	,intent(in)	:: this
		real(dkind)						:: get_save_time

		get_save_time = this%save_time
	end function

	pure function get_save_time_units_abbreviation(this)
		class(data_save)	,intent(in)	:: this
		character(len=10)				:: get_save_time_units_abbreviation

		get_save_time_units_abbreviation = this%save_time_units_abbreviation
	end function
	
	pure function get_save_time_coefficient(this)
		class(data_save)	,intent(in)	:: this
		real(rkind)			:: get_save_time_coefficient

		get_save_time_coefficient = this%save_time_coefficient
	end function

	pure function get_visible_fields_number(this)
		class(data_save)	,intent(in)	:: this
		integer							:: get_visible_fields_number

		get_visible_fields_number = size(this%visible_fields)
	end function

	pure subroutine set_save_time(this,new_save_time)
		class(data_save)	,intent(inout)	:: this
		real(dkind)			,intent(in)		:: new_save_time

		this%output_counter = this%output_counter * int(this%save_time/new_save_time)
		this%save_time = new_save_time
	end subroutine	
	
	subroutine write_log(this,log_unit)
		class(data_save)	,intent(in)	:: this
		integer				,intent(in)	:: log_unit
		
		write(log_unit,'(A)')	'************************************************************************************* '	
		write(log_unit,'(A)')			' Data save setup:'
		write(log_unit,'(A,E14.7,A)')	' Data save time: ', this%save_time, this%save_time_units_abbreviation 		
		write(log_unit,'(A,A)')			' Data save format: ', this%save_format
		write(log_unit,'(A,A)')			' Data save folder: ', this%data_save_folder
		write(log_unit,'(A)')	'************************************************************************************* '				
	end subroutine
	
	subroutine write_output_counter(this)
		class(data_save)	,intent(in)	:: this
		integer	:: output_counter
		
		integer	:: io_unit
		integer	:: ierr
		
		namelist /data_save_output_counter/ output_counter

		open(newunit = io_unit, file = data_save_data_file_name, status = 'old', form = 'formatted')
		read(unit = io_unit, nml = data_save_output_counter, iostat = ierr)
		
		if (ierr == 0 ) then
			backspace(io_unit)
			backspace(io_unit)
			backspace(io_unit)
		end if
		
		output_counter = this%output_counter		
		
		write(unit = io_unit, nml = data_save_output_counter)
		close(io_unit)
	end subroutine
	
	subroutine save_all_data(this,time,stop_flag,make_save)
		class(data_save)	,intent(inout)	:: this
		real(dkind)			,intent(in)		:: time
		logical				,intent(in)		,optional	:: make_save	
		logical				,intent(in)		:: stop_flag

		character(len=100)					:: system_command

#ifdef mpi
		integer(kind=MPI_OFFSET_KIND)	:: initial_displacement, final_displacement
#endif

		real(dkind)				:: written_time
		character(len=100)		:: file_path, file_name, proc_rank
		integer					:: unit_io, mpi_io_unit
		logical					:: debug
		integer					:: error

		logical					:: make_flag

		integer					:: processor_rank

		processor_rank = this%domain%get_processor_rank()

		make_flag = .false.
        if(stop_flag) make_flag = .true.
        
		if(present(make_save)) make_flag = make_save

		if((time*this%save_time_coefficient >= this%save_time*(this%output_counter+1)).or.(make_flag)) then

			select case (this%save_format)
				case('tecplot')
					written_time	= time*this%save_time_coefficient
					debug			= .false.
					if((make_flag).or.(this%debug_flag)) then
						debug = .true.
					end if

					write(file_path,'(I6.6,A)')  int(written_time),this%save_time_units_abbreviation
			
					if (this%debug_flag) then
						written_time = this%output_counter
						write(file_path,'(I6.6)')  int(written_time)
					end if

#ifdef mpi			

					if (processor_rank == 0) then

						system_command = 'mkdir ' // trim(this%data_save_folder) // trim(fold_sep) // trim(file_path)
						call system(system_command)

						file_name = trim(this%data_save_folder) // trim(fold_sep) // trim(file_path) // trim(fold_sep) // trim(file_path) //'_header.plt'
						open(newunit = unit_io, file = file_name, status = 'replace', form = 'binary')					
						call this%generate_header(written_time,unit_io,debug)
						close(unit_io)							
					end if
				
					call MPI_BARRIER(MPI_COMM_WORLD,error)

					write(proc_rank,'(A,I4.4)') '_proc_', processor_rank

					file_name = trim(this%data_save_folder) // trim(fold_sep) // trim(file_path) // trim(fold_sep) // trim(file_path) // trim(proc_rank) // '_data.plt'
					open(newunit = unit_io, file = file_name, status = 'replace', form = 'binary')

					call this%save_mesh(unit_io,debug)
					call this%save_fields(unit_io,debug)
					if(debug) then
						call this%save_bounds(unit_io,debug)
					end if
					close(unit_io)

					! call MPI_BARRIER(MPI_COMM_WORLD,error)

					! call this%commit_mpi_save_datatype()

					! file_name = trim(this%data_save_folder) // trim(fold_sep) // trim(file_path) //'_data.plt'
					! call MPI_FILE_OPEN(MPI_COMM_WORLD, file_name, MPI_MODE_WRONLY + MPI_MODE_CREATE, MPI_INFO_NULL, mpi_io_unit, error) 

					! initial_displacement = 0_MPI_OFFSET_KIND
					! call this%save_mesh_mpi(mpi_io_unit,initial_displacement,final_displacement)

					! initial_displacement = final_displacement

					! call this%save_fields_mpi(mpi_io_unit,initial_displacement,final_displacement)

					! initial_displacement = final_displacement

					! if(debug) then
					! 	call this%save_bounds_mpi(mpi_io_unit,initial_displacement,final_displacement)
					! end if

					! call MPI_BARRIER(MPI_COMM_WORLD,error)
					! call MPI_FILE_CLOSE(mpi_io_unit,error)
#else

					file_name = trim(this%data_save_folder) // trim(fold_sep) // trim(file_path) //'.plt'
					open(newunit = unit_io, file = file_name, status = 'replace', form = 'binary')

					call this%generate_header(written_time,unit_io,debug)

					call this%save_mesh(unit_io,debug)
					call this%save_fields(unit_io,debug)
					if(debug) then
						call this%save_bounds(unit_io,debug)
					end if
					close(unit_io)
#endif


			end select
			
			this%output_counter = this%output_counter + 1

		end if

		if(stop_flag.and.(processor_rank == 0)) call this%write_output_counter()

	end subroutine

	subroutine generate_header(this,time,tecplot_io_unit,debug)
		class(data_save)								,intent(in)				:: this
		real(dkind)									,intent(in)					:: time
		integer										,intent(in)					:: tecplot_io_unit
		logical										,intent(in)					:: debug

		integer						:: file_type = 1
		integer ,dimension(60)		:: converted_string
		character(len=40)			:: zone_name, zone_title

		integer						:: domain_dimensions
		integer		,dimension(3)	:: utter_cells_number, inner_cells_number

		character(len=5)	,dimension(:)	,allocatable	:: axis_names
		real(dkind)			,dimension(:,:)	,allocatable	:: domain_lengths

		real(dkind)				:: min_value, max_value
		integer					:: variables_number ,variables_counter
		integer					:: dim, counter

		domain_dimensions						= this%domain%get_domain_dimensions()
		utter_cells_number						= this%domain%get_global_cells_number()
		inner_cells_number						= 1
		inner_cells_number(:domain_dimensions)	= utter_cells_number(:domain_dimensions) - 2

		allocate(axis_names		, source	= this%domain%get_axis_names())
		allocate(domain_lengths	, source	= this%domain%get_domain_lengths())

		variables_number = 0
		variables_number = variables_number + domain_dimensions				! Mesh 
		if(debug) variables_number = variables_number + 1					! Boundary 
		variables_number = variables_number + size(this%visible_fields)

		!Header identification
		write(tecplot_io_unit) '#!TDV112'

		!Integer order
		write(tecplot_io_unit) 1

		!File type
		write(tecplot_io_unit) file_type

		!Title
		write(zone_title,'(A)') 'Flame_test'
		write(tecplot_io_unit) (converted_string(counter),counter=1,convert_string(converted_string,zone_title))

		!Number of variables
		write(tecplot_io_unit) variables_number

		!Variables names
		do dim = 1, domain_dimensions
			write(tecplot_io_unit)  (converted_string(counter),counter=1,convert_string(converted_string,axis_names(dim)))
		end do
		do variables_counter = 1,size(this%visible_fields)
			 write(tecplot_io_unit)  (converted_string(counter),counter=1,convert_string(converted_string,this%visible_fields(variables_counter)%s_ptr%name_short))
		end do
		if(debug) write(tecplot_io_unit) (converted_string(counter),counter=1,convert_string(converted_string,'bc_markers'))

		!Zone description
		!Zone marker
		write(tecplot_io_unit) 299.0
		!Zone name
		write(zone_name,'(A,E11.4,A)') 'zone', time, this%save_time_units_abbreviation
		write(tecplot_io_unit) (converted_string(counter),counter=1,convert_string(converted_string,zone_name))
		!Parent zone
		write(tecplot_io_unit)-1
		!StrandID
		write(tecplot_io_unit)-1
		!solution time
		write(tecplot_io_unit) time
		!Zone color
		write(tecplot_io_unit)-1
		!ZoneType
		write(tecplot_io_unit) 0
		!Data packing
		write(tecplot_io_unit) 1
		!Var location
		write(tecplot_io_unit) 0
		!Raw local neighbors
		write(tecplot_io_unit) 0
		!Miscellaneous user-defined face
		do variables_counter = 1,variables_number
			 write(tecplot_io_unit)  0
		end do
		!IMAX,JMAX,KMAX

		if(debug) then
			write(tecplot_io_unit) utter_cells_number
		else
			write(tecplot_io_unit) inner_cells_number
		end if

		!Auxilary data name
		write(tecplot_io_unit) 0
		!Data section marker
		write(tecplot_io_unit) 357.0
		!Zone implementation
		write(tecplot_io_unit) 299.0
		!Var types
		do variables_counter = 1,variables_number
			 write(tecplot_io_unit)  1
		end do
		!Passive vars
		write(tecplot_io_unit) 0
		!Variable sharing
		write(tecplot_io_unit) 0
		!Share connectivity list
		write(tecplot_io_unit) -1
		do dim = 1,domain_dimensions
			write(tecplot_io_unit) domain_lengths(dim,1)
			write(tecplot_io_unit) domain_lengths(dim,2)
		end do
		do variables_counter = 1,size(this%visible_fields)
			 !Min Value
			 call this%visible_fields(variables_counter)%s_ptr%get_field_min(bc_ptr=this%boundaries,min_value=min_value)
			 write(tecplot_io_unit) min_value
			 !Max Value
			 call this%visible_fields(variables_counter)%s_ptr%get_field_max(bc_ptr=this%boundaries,max_value=max_value)
			 write(tecplot_io_unit) max_value
		end do
		if(debug) then
			write(tecplot_io_unit) 0.0_dkind
			write(tecplot_io_unit) dble(size(this%boundaries%bc_ptr%boundary_types))
		end if

	end subroutine

	subroutine save_mesh(this,tecplot_io_unit,debug)
		class(data_save)	,intent(inout)	:: this
		integer				,intent(in)		:: tecplot_io_unit
		logical				,intent(in)		:: debug

		integer						:: dimensions

		integer		,dimension(3,2)	:: loop
		integer		,dimension(3,2)	:: inner_loop, utter_loop

		integer	:: i,j,k, dim

		dimensions				= this%domain%get_domain_dimensions()
		utter_loop				= this%domain%get_local_utter_cells_bounds()
		inner_loop				= this%domain%get_local_inner_cells_bounds()

		if (debug) then	
			loop = utter_loop
		else
			loop = inner_loop
		end if

		do dim = 1,dimensions
			this%io_buffer = this%mesh%mesh_ptr%mesh(dim,:,:,:)

			write(tecplot_io_unit) (((this%io_buffer(i,j,k),	i = loop(1,1),loop(1,2)), &
																j = loop(2,1),loop(2,2)), &
																k = loop(3,1),loop(3,2))
		end do
	end subroutine

	subroutine save_fields(this,tecplot_io_unit,debug)
		class(data_save)	,intent(inout)	:: this
		integer				,intent(in)		:: tecplot_io_unit
		logical				,intent(in)		:: debug

		integer		,dimension(3,2)	:: loop
		integer		,dimension(3,2)	:: inner_loop, utter_loop

		integer	:: i,j,k, fields_count

		utter_loop				= this%domain%get_local_utter_cells_bounds()
		inner_loop				= this%domain%get_local_inner_cells_bounds()

		if (debug) then	
			loop = utter_loop
		else
			loop = inner_loop
		end if

		do fields_count = 1,size(this%visible_fields)
			this%io_buffer = this%visible_fields(fields_count)%s_ptr%cells(:,:,:)
			write(tecplot_io_unit) (((this%io_buffer(i,j,k),	i = loop(1,1),loop(1,2)), &
																j = loop(2,1),loop(2,2)), &
																k = loop(3,1),loop(3,2))
		end do
	end subroutine

	subroutine save_bounds(this,tecplot_io_unit,debug)
		class(data_save)	,intent(inout)	:: this
		integer				,intent(in)		:: tecplot_io_unit
		logical				,intent(in)		:: debug

		integer		,dimension(3,2)	:: loop
		integer		,dimension(3,2)	:: inner_loop, utter_loop

		integer	:: i,j,k, fields_count

		utter_loop				= this%domain%get_local_utter_cells_bounds()
		inner_loop				= this%domain%get_local_inner_cells_bounds()

		if (debug) then	
			loop = utter_loop
		else
			loop = inner_loop
		end if

		this%io_buffer = this%boundaries%bc_ptr%bc_markers(:,:,:)
		write(tecplot_io_unit) (((this%io_buffer(i,j,k),	i = loop(1,1),loop(1,2)), &
															j = loop(2,1),loop(2,2)), &
															k = loop(3,1),loop(3,2))

	end subroutine

	subroutine commit_mpi_save_datatype(this)
		class(data_save)	,intent(inout)	:: this

		integer					:: dimensions
		integer ,dimension(3)	:: sizes, subsizes, starts

		integer		:: error

		this%mpi_save_array_type = 0

#ifdef mpi

		dimensions	= this%domain%get_domain_dimensions()

		sizes		= this%domain%get_global_cells_number() 
		subsizes	= this%domain%get_local_cells_number() 
		starts		= this%domain%get_global_offset() 

		sizes(1:dimensions)		= sizes(1:dimensions)	 - 2
		subsizes(1:dimensions)	= subsizes(1:dimensions) - 2
		starts(1:dimensions)	= starts(1:dimensions)		

		call MPI_TYPE_CREATE_SUBARRAY(3,sizes,subsizes,starts, MPI_ORDER_FORTRAN, MPI_REAL, this%mpi_save_array_type, error)
		call MPI_TYPE_COMMIT(this%mpi_save_array_type,error)
#endif
	end subroutine

	subroutine save_mesh_mpi(this,tecplot_io_unit, initial_displacement, final_displacement)
		class(data_save)	,intent(inout)	:: this
		integer				,intent(inout)	:: tecplot_io_unit

		integer					:: dimensions
		integer	,dimension(3)	:: local_cells_number, global_cells_number
		integer	,dimension(3,2)	:: inner_loop, inner_loop_max

		integer 				:: buffer_size
		integer					:: i,j,k, dim, dim1
		integer					:: request, error, info

#ifdef mpi	
		integer(kind=MPI_OFFSET_KIND)	,intent(inout)	,optional	:: initial_displacement, final_displacement
		integer(kind=MPI_OFFSET_KIND)	:: displacement
		integer 						:: status(MPI_STATUS_SIZE)
#else
		integer							,intent(inout)	,optional	:: initial_displacement, final_displacement
		integer							:: displacement
#endif

#ifdef mpi

		dimensions		= this%domain%get_domain_dimensions()
		
		inner_loop		= this%domain%get_local_inner_cells_bounds()
		inner_loop_max	= this%domain%get_local_inner_cells_bounds_max()

		global_cells_number	= this%domain%get_global_cells_number()
		local_cells_number	= this%domain%get_local_cells_number()

		buffer_size = local_cells_number(1) - 2
		do dim = 1,dimensions
			this%io_buffer = this%mesh%mesh_ptr%mesh(dim,:,:,:)
			do k = inner_loop_max(3,1),inner_loop_max(3,2) 
			do j = inner_loop_max(2,1),inner_loop_max(2,2) 

				displacement = (dim-1) * sizeof(this%io_buffer(1,1,1))

				do dim1 = 1, dimensions
					displacement = displacement * (global_cells_number(dim1) - 2) 
				end do
	
				displacement =  displacement + initial_displacement

				displacement =  displacement + ((k-inner_loop_max(3,1)) * (global_cells_number(2) - 2) * (global_cells_number(1) - 2) + &
												(j-inner_loop_max(2,1)) * (global_cells_number(1) - 2)) * sizeof(this%io_buffer(1,1,1))

				call MPI_FILE_SET_VIEW(tecplot_io_unit, displacement, MPI_REAL, this%mpi_save_array_type, 'native', MPI_INFO_NULL, error) 

				if ((k <=  inner_loop(3,2)).and.(j <= inner_loop(2,2))) then
					call MPI_FILE_IWRITE(tecplot_io_unit, this%io_buffer(inner_loop(1,1):inner_loop(1,2),j,k), buffer_size, MPI_REAL, request, error) 
				end if
			end do
			end do
		end do

		final_displacement = dimensions * sizeof(this%io_buffer(1,1,1))

		do dim1 = 1, dimensions
			final_displacement = final_displacement * (global_cells_number(dim1) - 2)
		end do

		final_displacement = final_displacement + initial_displacement

		call MPI_BARRIER(MPI_COMM_WORLD,error)

#endif
	end subroutine

subroutine save_fields_mpi(this,tecplot_io_unit,initial_displacement,final_displacement)
		class(data_save)	,intent(inout)	:: this
		integer				,intent(inout)	:: tecplot_io_unit
		
		integer					:: dimensions
		integer	,dimension(3)	:: local_cells_number, global_cells_number
		integer	,dimension(3,2)	:: inner_loop, inner_loop_max

		integer :: buffer_size
		integer	:: fields_count, i,j,k, dim, dim1
		integer	:: request, error, info

#ifdef mpi	
		integer(kind=MPI_OFFSET_KIND)	,intent(inout)	,optional	:: initial_displacement, final_displacement
		integer(kind=MPI_OFFSET_KIND)	:: displacement
		integer 						:: status(MPI_STATUS_SIZE)
#else
		integer							,intent(inout)	,optional	:: initial_displacement, final_displacement
		integer							:: displacement
#endif

#ifdef mpi

		dimensions		= this%domain%get_domain_dimensions()
		
		inner_loop		= this%domain%get_local_inner_cells_bounds()
		inner_loop_max	= this%domain%get_local_inner_cells_bounds_max()

		global_cells_number	= this%domain%get_global_cells_number()
		local_cells_number	= this%domain%get_local_cells_number()

		buffer_size = local_cells_number(1) - 2
		do fields_count = 1,size(this%visible_fields)
			this%io_buffer = this%visible_fields(fields_count)%s_ptr%cells(:,:,:)
			do k = inner_loop_max(3,1),inner_loop_max(3,2) 
			do j = inner_loop_max(2,1),inner_loop_max(2,2) 

				displacement = (fields_count-1) * sizeof(this%io_buffer(1,1,1))

				do dim1 = 1, dimensions
					displacement = displacement * (global_cells_number(dim1) - 2)
				end do
	
				displacement =  displacement + initial_displacement

				displacement =  displacement + ((k-inner_loop_max(3,1)) * (global_cells_number(2) - 2)  * (global_cells_number(1) - 2)  + &
												(j-inner_loop_max(2,1)) * (global_cells_number(1) - 2)) * sizeof(this%io_buffer(1,1,1))

				call MPI_FILE_SET_VIEW(tecplot_io_unit, displacement, MPI_REAL, this%mpi_save_array_type, 'native', MPI_INFO_NULL, error) 

				if ((k <=  inner_loop(3,2)).and.(j <= inner_loop(2,2))) then
					call MPI_FILE_IWRITE(tecplot_io_unit, this%io_buffer(inner_loop(1,1):inner_loop(1,2),j,k), buffer_size, MPI_REAL, request, error) 
				end if
			end do
			end do
		end do

		final_displacement = (size(this%visible_fields)) * sizeof(this%io_buffer(1,1,1))

		do dim1 = 1, dimensions
			final_displacement = final_displacement * (global_cells_number(dim1) - 2)
		end do

		final_displacement = final_displacement + initial_displacement

		call MPI_BARRIER(MPI_COMM_WORLD,error)

#endif

	end subroutine

subroutine save_bounds_mpi(this,tecplot_io_unit,initial_displacement,final_displacement)
		class(data_save)	,intent(inout)	:: this
		integer				,intent(inout)	:: tecplot_io_unit

		integer					:: dimensions
		integer	,dimension(3)	:: local_cells_number, global_cells_number
		integer	,dimension(3,2)	:: inner_loop, inner_loop_max

		integer :: buffer_size
		integer	:: i,j,k, dim, dim1
		integer	:: request, error, info

#ifdef mpi	
		integer(kind=MPI_OFFSET_KIND)	,intent(inout)	,optional	:: initial_displacement, final_displacement
		integer(kind=MPI_OFFSET_KIND)	:: displacement
		integer 						:: status(MPI_STATUS_SIZE)
#else
		integer							,intent(inout)	,optional	:: initial_displacement, final_displacement
		integer							:: displacement
#endif

#ifdef mpi

		dimensions		= this%domain%get_domain_dimensions()
		
		inner_loop		= this%domain%get_local_inner_cells_bounds()
		inner_loop_max	= this%domain%get_local_inner_cells_bounds_max()

		global_cells_number	= this%domain%get_global_cells_number()
		local_cells_number	= this%domain%get_local_cells_number()

		buffer_size = local_cells_number(1) - 2

		this%io_buffer = this%boundaries%bc_ptr%bc_markers(:,:,:)
		do k = inner_loop_max(3,1),inner_loop_max(3,2) 
		do j = inner_loop_max(2,1),inner_loop_max(2,2) 

			displacement = initial_displacement 

			displacement =  displacement + ((k-inner_loop_max(3,1)) * (global_cells_number(2) - 2) * (global_cells_number(1) - 2)  + &
											(j-inner_loop_max(2,1)) * (global_cells_number(1) - 2)) * sizeof(this%io_buffer(1,1,1))

			call MPI_FILE_SET_VIEW(tecplot_io_unit, displacement, MPI_REAL, this%mpi_save_array_type, 'native', MPI_INFO_NULL, error) 

			if ((k <=  inner_loop(3,2)).and.(j <= inner_loop(2,2))) then
				call MPI_FILE_IWRITE(tecplot_io_unit, this%io_buffer(inner_loop(1,1):inner_loop(1,2),j,k), buffer_size, MPI_REAL, request, error) 
			end if
		end do
		end do

		final_displacement = sizeof(this%io_buffer(1,1,1))

		do dim1 = 1, dimensions
			final_displacement = final_displacement * (global_cells_number(dim1) - 2)
		end do

		call MPI_BARRIER(MPI_COMM_WORLD,error)

#endif

	end subroutine

	integer function convert_string(int_str,char_str) result (num)
		 integer ,dimension(:)		::  int_str
		 character(len=*)				::  char_str

		 int_str = 0
		 do num = 1, len(char_str)
			  int_str(num) = ichar(char_str(num:num))
		 end do
		 int_str(num) = 0
		 continue
	end function

end module
