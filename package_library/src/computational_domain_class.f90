module computational_domain_class

	use kind_parameters
	use global_data
	
#ifdef mpi
    use MPI
#endif

	implicit none

	private
	public	::	computational_domain, computational_domain_c

	type computational_domain
		private
	! Global data
		integer												:: dimensions			! Problem dimensions
		integer				,dimension(3)					:: cells_number			! Cells number without ghost cells
		integer				,dimension(3)					:: faces_number			! Faces number without ghost cells
		real(dkind)			,dimension(:,:)	,allocatable	:: lengths				! Domain length in meters
		character(len=5)	,dimension(:)	,allocatable	:: axis_names			! Axis names
		character(len=20)									:: coordinate_system	! Coordinate system (cartesian/cylindrical/spherical)

	! MPI Global data
		integer												:: mpi_communicator			
		integer												:: mpi_communicator_size	! Number of processors in communicator
		integer				,dimension(3)					:: processor_number			! Cartesian processor grid dimensions

	! MPI Local data
		integer												:: processor_rank					! Processor rank in MPI_HELLO_WORLD communicator
		integer				,dimension(3)					:: processor_grid_coord				! Processor coordinate in cartesian processor grid

		integer				,dimension(3)					:: cells_number_decomposed			! Local cells number
		integer				,dimension(3)					:: cells_number_decomposed_max		! Maximum local cells number (for MPI_SET_VIEW)
		integer				,dimension(3)					:: faces_number_decomposed			! Local faces number
		integer				,dimension(3)					:: faces_number_decomposed_max

		integer				,dimension(3)					:: global_offset					! Offset from local (0,0,0) cell to global (0,0,0) cell
	contains
		! IO Getters/Setters
		procedure	,private	:: read_properties
		procedure	,private	:: write_properties
		procedure	,private	:: set_properties

		procedure				:: generate_processors_grid
		procedure				:: decompose_domain	

		! Getters
		procedure				:: get_domain_dimensions
		procedure				:: get_domain_lengths
		procedure				:: get_axis_number
		procedure				:: get_axis_names
		procedure				:: get_coordinate_system_name
		
		procedure				:: get_global_cells_number
		procedure				:: get_global_faces_number
		procedure				:: get_global_inner_cells_bounds
		procedure				:: get_global_utter_cells_bounds
		procedure				:: get_global_inner_faces_bounds
		procedure				:: get_global_utter_faces_bounds	

		procedure				:: get_processor_number
		procedure				:: get_processor_grid_coord
		procedure				:: get_processor_rank
		procedure				:: get_global_offset

		procedure				:: get_local_cells_number
		procedure				:: get_local_faces_number		
		procedure				:: get_local_cells_number_max
		procedure				:: get_local_faces_number_max
		procedure				:: get_local_inner_cells_bounds
		procedure				:: get_local_utter_cells_bounds	
		procedure				:: get_local_inner_cells_bounds_max
		procedure				:: get_local_utter_cells_bounds_max				
		procedure				:: get_local_inner_faces_bounds
		procedure				:: get_local_utter_faces_bounds
		procedure				:: get_local_inner_faces_bounds_max
		procedure				:: get_local_utter_faces_bounds_max		

		procedure				:: get_mpi_communicator	
		procedure				:: get_mpi_communicator_size	
		
		procedure				:: cell_is_inside_local_domain
		procedure				:: convert_global_cell_to_local
		procedure				:: area_intersects_local_domain

		! Setters
		procedure				:: set_mpi_communicator_size
		procedure				:: set_processor_grid_coord

		! Logger
		procedure				:: write_log
	end type computational_domain

	interface computational_domain_c
		module procedure constructor
		module procedure constructor_file
	end interface

contains

	type(computational_domain)	function constructor(dimensions,cells_number,coordinate_system,lengths,axis_names)
		integer								,intent(in)	:: dimensions
		integer				,dimension(3)	,intent(in)	:: cells_number
		character(len=*)					,intent(in)	:: coordinate_system		
		real(dkind)			,dimension(:,:)	,intent(in)	:: lengths
		character(len=*)	,dimension(:)	,intent(in)	:: axis_names

		integer	:: io_unit
		
		call constructor%set_properties(dimensions,cells_number,lengths,coordinate_system,axis_names)
		
		open(newunit = io_unit, file = domain_data_file_name, status = 'replace', form = 'formatted')
		call constructor%write_properties(io_unit)
		close(io_unit)				
	end function

	type(computational_domain)	function constructor_file()
		integer	:: io_unit
		
		open(newunit = io_unit, file = domain_data_file_name, status = 'old', form = 'formatted')
		call constructor_file%read_properties(io_unit)
		close(io_unit)	
	end function
	
	subroutine read_properties(this,domain_data_unit)
		class(computational_domain)	,intent(inout)	:: this
		integer						,intent(in)		:: domain_data_unit
	
		integer								:: dimensions
		integer				,dimension(3)	:: cells_number
		character(len=20)					:: coordinate_system
		real(dkind)			,dimension(:,:)	,allocatable	:: lengths
		character(len=5)	,dimension(:)	,allocatable	:: axis_names

		namelist /domain_properties_1/ dimensions, cells_number, coordinate_system
		namelist /domain_properties_2/ lengths, axis_names
		
		read(unit = domain_data_unit, nml = domain_properties_1)
		
		allocate(lengths(dimensions,2))
		allocate(axis_names(dimensions))
		
		read(unit = domain_data_unit, nml = domain_properties_2)
	
		call this%set_properties(dimensions,cells_number,lengths,coordinate_system,axis_names)
	end subroutine
	
	subroutine write_properties(this,domain_data_unit)
		class(computational_domain)	,intent(in)	:: this
		integer						,intent(in)	:: domain_data_unit
	
		integer								:: dimensions
		integer				,dimension(3)	:: cells_number
		character(len=20)					:: coordinate_system
		real(dkind)			,dimension(:,:)	,allocatable	:: lengths
		character(len=5)	,dimension(:)	,allocatable	:: axis_names

		namelist /domain_properties_1/ dimensions, cells_number, coordinate_system
		namelist /domain_properties_2/ lengths, axis_names
		
		dimensions			= this%dimensions
		cells_number		= this%cells_number
		coordinate_system	= this%coordinate_system
		
		allocate(axis_names(dimensions))
		allocate(lengths(dimensions,2))
		
		axis_names	= this%axis_names
		lengths		= this%lengths
		
		write(unit = domain_data_unit, nml = domain_properties_1)
		write(unit = domain_data_unit, nml = domain_properties_2)
	end subroutine
	
	subroutine set_properties(this,dimensions,cells_number,lengths,coordinate_system,axis_names)
		class(computational_domain)			,intent(inout)	:: this
		integer								,intent(in)		:: dimensions
		integer		,dimension(3)			,intent(in)		:: cells_number
		real(dkind)	,dimension(:,:)			,intent(in)		:: lengths
		character(len=*)					,intent(in)		:: coordinate_system
		character(len=*)	,dimension(:)	,intent(in)		:: axis_names
		
		integer	:: dim
		
		allocate(this%axis_names(dimensions))
		allocate(this%lengths(dimensions,2))

		this%dimensions				= dimensions
		this%cells_number			= cells_number		
		this%faces_number			= cells_number + 1
		this%lengths				= lengths(1:dimensions,:)
		this%coordinate_system		= coordinate_system
		this%axis_names				= axis_names(1:dimensions)

		this%cells_number(dimensions+1:3) = 1
		this%faces_number(dimensions+1:3) = 1

		this%mpi_communicator_size = 1
		this%processor_rank = 0

		call this%generate_processors_grid()
		call this%decompose_domain()
		
		if ((this%coordinate_system == 'spherical').and.(this%dimensions > 1))		print *, 'ERROR: dimension amount is not allowed to be greater than 1 for spherical symmetry. Check your problem setup.'
		if ((this%coordinate_system == 'cylindrical').and.(this%dimensions > 2))	print *, 'ERROR: dimension amount is not allowed to be greater than 2 for spherical symmetry. Check your problem setup.'
	end subroutine

	subroutine write_log(this,log_unit)
		class(computational_domain)	,intent(in)	:: this
		
		integer						,intent(in)	:: log_unit			
		integer									:: i

		write(log_unit,'(A)')		'************************************************************************************* '
		write(log_unit,'(A)')		' Computational domain setup : '
		write(log_unit,'(A,I3)')	' Domain dimensions         : ',	this%dimensions
		write(log_unit,'(A,3I5)')	' Domain cells number       : ',	this%cells_number
		write(log_unit,'(A,E14.7)')	' Domain lower bounds       : ',	this%lengths(1,1)
		write(log_unit,'(A,E14.7)')	' Domain upper bounds       : ',	this%lengths(1,2)
		write(log_unit,'(A,A)')		' Domain coordinate system  : ',	this%coordinate_system
		write(log_unit,'(A,3A)')	' Domain axis names         : ',	this%axis_names	
		write(log_unit,'(A)')		'************************************************************************************* '
	end subroutine
	
	subroutine generate_processors_grid(this)
		class(computational_domain)	,intent(inout)	:: this

		integer						:: comm_size
		integer						:: comm_rank
		
		integer		,dimension(3)	:: cells_number_max
		integer		,dimension(3)	:: cells_number_residual

		integer		:: min_data_transfer, data_transfer
		logical		,dimension(3)	:: is_periodic
		logical						:: reorder
		integer						:: icomm_cart
		
		integer		,dimension(3)	:: loop
		integer		:: i,j,k,dim			
		integer		:: error


#ifdef mpi
    	call MPI_COMM_SIZE(MPI_COMM_WORLD, this%mpi_communicator_size, error)
		call MPI_COMM_RANK(MPI_COMM_WORLD, this%processor_rank, error)
#endif

		this%processor_number(1) = this%mpi_communicator_size
		this%processor_number(2) = 1
		this%processor_number(3) = 1

		min_data_transfer = (this%processor_number(1)-1)*(this%cells_number(2)*this%cells_number(3)) + &
							(this%processor_number(2)-1)*(this%cells_number(1)*this%cells_number(3)) + &
							(this%processor_number(3)-1)*(this%cells_number(1)*this%cells_number(2))
		
		loop = 1
		loop(1:this%dimensions) = this%mpi_communicator_size
							
		do i = 1, loop(1)
		do j = 1, loop(2)
		do k = 1, loop(3)
			if ((i*j*k) == this%mpi_communicator_size) then
				data_transfer = (i-1)*(this%cells_number(2)*this%cells_number(3)) + &
								(j-1)*(this%cells_number(1)*this%cells_number(3)) + &
								(k-1)*(this%cells_number(1)*this%cells_number(2))
				if ( data_transfer < min_data_transfer) then
					min_data_transfer = data_transfer
					this%processor_number(1) = i
					this%processor_number(2) = j
					this%processor_number(3) = k
				end if
			end if
		end do
		end do
		end do

		is_periodic = .false.
		reorder  	= .false.

		this%processor_grid_coord	= 0
		this%global_offset			= 0

#ifdef mpi
		call MPI_CART_CREATE(MPI_COMM_WORLD, 3, this%processor_number, is_periodic,reorder, this%mpi_communicator, error)
		call MPI_CART_COORDS(this%mpi_communicator, this%processor_rank, 3, this%processor_grid_coord, error)
#endif		

	end subroutine

	subroutine decompose_domain(this)

		class(computational_domain)	,intent(inout)	:: this

		integer						:: comm_size
		integer						:: comm_rank
		
		integer		,dimension(3)	:: cells_number_max
		integer		,dimension(3)	:: cells_number_residual

		integer		:: min_data_transfer, data_transfer
		logical		,dimension(3)	:: is_periodic
		logical						:: reorder
		integer						:: icomm_cart
		
		integer		,dimension(3)	:: loop
		integer		:: i,j,k,dim			
		integer		:: error

		this%cells_number_decomposed	= this%cells_number / this%processor_number
		cells_number_residual 			= this%cells_number - this%processor_number * this%cells_number_decomposed

		this%cells_number_decomposed_max = this%cells_number_decomposed
		do dim = 1, this%dimensions
			this%cells_number_decomposed_max(dim) = this%cells_number_decomposed(dim) + 1
		end do

		do dim = 1, this%dimensions
			this%global_offset(dim) = 	this%cells_number_decomposed(dim) * this%processor_grid_coord(dim) + &
									min(this%processor_grid_coord(dim),cells_number_residual(dim))
		end do

		do dim = 1,3
			if ((this%processor_grid_coord(dim) < cells_number_residual(dim) ).and.(cells_number_residual(dim) /= 0)) then
				this%cells_number_decomposed(dim) = this%cells_number_decomposed(dim) + 1
			end if
		end do

		this%faces_number_decomposed 					= this%cells_number_decomposed
		this%faces_number_decomposed(1:this%dimensions) = this%cells_number_decomposed(1:this%dimensions) + 1

		this%faces_number_decomposed_max 					= this%cells_number_decomposed_max
		this%faces_number_decomposed_max(1:this%dimensions) = this%cells_number_decomposed_max(1:this%dimensions) + 1		

#ifdef mpi
		if (this%processor_rank == 0) then
			print *, ' Mpi domain decomposition finished.'
			print *, ' Processor grid: ', this%processor_number
		end if
#endif			
		
	end subroutine

! ************** Getters ***************

	pure integer function get_domain_dimensions(this)
		class(computational_domain)	,intent(in)	:: this

		get_domain_dimensions = this%dimensions
	end function

	pure function get_domain_lengths(this)
		class(computational_domain)	,intent(in)		:: this
		real(dkind)	,dimension(:,:)	,allocatable	:: get_domain_lengths

		
		allocate(get_domain_lengths(this%dimensions,2))
		get_domain_lengths = this%lengths
	end function

	pure function get_global_cells_number(this)
		class(computational_domain)	,intent(in)		:: this
		integer	,dimension(3)						:: get_global_cells_number

		get_global_cells_number = 1
		get_global_cells_number(1:this%dimensions) = this%cells_number(1:this%dimensions) + 2				! Global cells number WITH ghost cells
	end function

	pure function get_local_cells_number(this)
		class(computational_domain)	,intent(in)		:: this
		integer	,dimension(3)						:: get_local_cells_number

		get_local_cells_number = 1
		get_local_cells_number(1:this%dimensions) = this%cells_number_decomposed(1:this%dimensions) + 2	! Local cells number WITH ghost cells
	end function	

	pure function get_local_cells_number_max(this)
		class(computational_domain)	,intent(in)		:: this
		integer	,dimension(3)						:: get_local_cells_number_max

		get_local_cells_number_max = 1
		get_local_cells_number_max(1:this%dimensions) = this%cells_number_decomposed_max(1:this%dimensions) + 2
	end function

	pure function get_global_faces_number(this)
		class(computational_domain)	,intent(in)		:: this
		integer	,dimension(3)						:: get_global_faces_number

		get_global_faces_number = 1
		get_global_faces_number(1:this%dimensions) = this%faces_number(1:this%dimensions) + 2				! Global faces number WITH ghost cells
	end function

	pure function get_local_faces_number(this)
		class(computational_domain)	,intent(in)		:: this
		integer	,dimension(3)			:: get_local_faces_number

		get_local_faces_number = 1
		get_local_faces_number(1:this%dimensions) = this%faces_number_decomposed(1:this%dimensions) + 2	! Local faces number WITH ghost cells
	end function

	pure function get_local_faces_number_max(this)
		class(computational_domain)	,intent(in)		:: this
		integer	,dimension(3)						:: get_local_faces_number_max

		get_local_faces_number_max = 1
		get_local_faces_number_max(1:this%dimensions) = this%faces_number_decomposed_max(1:this%dimensions) + 2
	end function

	pure function get_global_offset(this)
		class(computational_domain)	,intent(in)		:: this
		integer	,dimension(3)						:: get_global_offset

		get_global_offset = this%global_offset
	end function

	pure integer function get_axis_number(this,axis_name)
		class(computational_domain)	,intent(in)	:: this
		character(len=*)			,intent(in)	:: axis_name
	
		integer	:: axis
		
		do axis = 1, this%dimensions
			if ( axis_name == this%axis_names(axis)) get_axis_number = axis		
		end do
	end function

	pure function get_axis_names(this)
		class(computational_domain)			,intent(in)		:: this
		character(len=5)	,dimension(:)	,allocatable	:: get_axis_names

		allocate(get_axis_names,source = this%axis_names)
	end function	
	
	pure function get_coordinate_system_name(this)
		class(computational_domain)			,intent(in)		:: this
		character(len=20)	:: get_coordinate_system_name
		
		get_coordinate_system_name = this%coordinate_system
	end function	

	pure integer function get_processor_rank(this)
		class(computational_domain)	,intent(in)		:: this

		get_processor_rank = this%processor_rank
	end function

	pure function get_processor_number(this)
		class(computational_domain)	,intent(in)		:: this
		integer	,dimension(3)						:: get_processor_number

		get_processor_number = this%processor_number
	end function

	pure function get_processor_grid_coord(this)
		class(computational_domain)	,intent(in)		:: this
		integer	,dimension(3)						:: get_processor_grid_coord

		get_processor_grid_coord = this%processor_grid_coord
	end function

	pure function get_global_inner_cells_bounds(this)
		class(computational_domain)	,intent(in)		:: this
		integer	,dimension(3,2)						:: get_global_inner_cells_bounds

		get_global_inner_cells_bounds		= 1
		get_global_inner_cells_bounds(:,2)	= this%cells_number
	end function

	pure function get_global_utter_cells_bounds(this)
		class(computational_domain)	,intent(in)		:: this
		integer	,dimension(3,2)						:: get_global_utter_cells_bounds

		get_global_utter_cells_bounds = 1
		get_global_utter_cells_bounds(:this%dimensions,1) = 0
		get_global_utter_cells_bounds(:this%dimensions,2) = this%cells_number(:this%dimensions) + 1
	end function

	pure function get_global_inner_faces_bounds(this)
		class(computational_domain)	,intent(in)		:: this
		integer	,dimension(3,2)						:: get_global_inner_faces_bounds

		get_global_inner_faces_bounds		= 1
		get_global_inner_faces_bounds(:,2)	= this%faces_number
	end function

	pure function get_global_utter_faces_bounds(this)
		class(computational_domain)	,intent(in)		:: this
		integer	,dimension(3,2)						:: get_global_utter_faces_bounds

		get_global_utter_faces_bounds = 1
		get_global_utter_faces_bounds(:this%dimensions,1) = 0
		get_global_utter_faces_bounds(:this%dimensions,2) = this%faces_number_decomposed(:this%dimensions) + 1
	end function

	pure function get_local_inner_cells_bounds(this)
		class(computational_domain)	,intent(in)		:: this
		integer	,dimension(3,2)						:: get_local_inner_cells_bounds

		get_local_inner_cells_bounds		= 1
		get_local_inner_cells_bounds(:,2)	= this%cells_number_decomposed
	end function


	pure function get_local_utter_cells_bounds(this)
		class(computational_domain)	,intent(in)		:: this
		integer	,dimension(3,2)						:: get_local_utter_cells_bounds

		get_local_utter_cells_bounds = 1
		get_local_utter_cells_bounds(:this%dimensions,1) = 0
		get_local_utter_cells_bounds(:this%dimensions,2) = this%cells_number_decomposed(:this%dimensions) + 1
	end function

	pure function get_local_inner_cells_bounds_max(this)
		class(computational_domain)	,intent(in)		:: this
		integer	,dimension(3,2)						:: get_local_inner_cells_bounds_max

		get_local_inner_cells_bounds_max		= 1
		get_local_inner_cells_bounds_max(:,2)	= this%cells_number_decomposed_max
	end function


	pure function get_local_utter_cells_bounds_max(this)
		class(computational_domain)	,intent(in)		:: this
		integer	,dimension(3,2)						:: get_local_utter_cells_bounds_max

		get_local_utter_cells_bounds_max = 1
		get_local_utter_cells_bounds_max(:this%dimensions,1) = 0
		get_local_utter_cells_bounds_max(:this%dimensions,2) = this%cells_number_decomposed_max(:this%dimensions) + 1
	end function

	pure function get_local_inner_faces_bounds(this)
		class(computational_domain)	,intent(in)		:: this
		integer	,dimension(3,2)						:: get_local_inner_faces_bounds

		get_local_inner_faces_bounds		= 1
		get_local_inner_faces_bounds(:,2)	= this%faces_number_decomposed
	end function

	pure function get_local_utter_faces_bounds(this)
		class(computational_domain)	,intent(in)		:: this
		integer	,dimension(3,2)						:: get_local_utter_faces_bounds

		get_local_utter_faces_bounds = 1
		get_local_utter_faces_bounds(:this%dimensions,1) = 0
		get_local_utter_faces_bounds(:this%dimensions,2) = this%faces_number_decomposed(:this%dimensions) + 1
	end function	

	pure function get_local_inner_faces_bounds_max(this)
		class(computational_domain)	,intent(in)		:: this
		integer	,dimension(3,2)						:: get_local_inner_faces_bounds_max

		get_local_inner_faces_bounds_max		= 1
		get_local_inner_faces_bounds_max(:,2)	= this%faces_number_decomposed_max
	end function


	pure function get_local_utter_faces_bounds_max(this)
		class(computational_domain)	,intent(in)		:: this
		integer	,dimension(3,2)						:: get_local_utter_faces_bounds_max

		get_local_utter_faces_bounds_max = 1
		get_local_utter_faces_bounds_max(:this%dimensions,1) = 0
		get_local_utter_faces_bounds_max(:this%dimensions,2) = this%faces_number_decomposed_max(:this%dimensions) + 1
	end function

	pure integer function get_mpi_communicator(this)
		class(computational_domain)	,intent(in)	:: this

		get_mpi_communicator = this%mpi_communicator
	end function

	pure integer function get_mpi_communicator_size(this)
		class(computational_domain)	,intent(in)	:: this

		get_mpi_communicator_size = this%mpi_communicator_size
	end function	


	subroutine set_processor_grid_coord(this,processor_grid_coord)
		class(computational_domain)	,intent(inout)		:: this
		integer	,dimension(3)		,intent(in)			:: processor_grid_coord

		this%processor_grid_coord = processor_grid_coord 
	end subroutine
	
	subroutine set_mpi_communicator_size(this,mpi_communicator_size)
		class(computational_domain)	,intent(inout)		:: this
		integer						,intent(in)			:: mpi_communicator_size

		this%mpi_communicator_size = mpi_communicator_size 
	end subroutine	

	pure logical function cell_is_inside_local_domain(this,cell)
		class(computational_domain)	,intent(in)	:: this
		integer	,dimension(3)		,intent(in)	:: cell

		integer	:: dim
		
		cell_is_inside_local_domain = .true.
		do dim = 1, this%dimensions
		
			if ((1 + this%global_offset(dim) > cell(dim)) .or. &
			    (this%cells_number_decomposed(dim) + this%global_offset(dim) < cell(dim))) then
				cell_is_inside_local_domain = .false.
			end if
			   
		end do
	end function
	
	subroutine convert_global_cell_to_local(this,cell)
		class(computational_domain)	,intent(in)			:: this
		integer	,dimension(3)		,intent(inout)		:: cell

		integer	:: dim
		
		do dim = 1, this%dimensions
			cell(dim) = cell(dim) - this%global_offset(dim)
		end do
		
	end subroutine
	
	pure logical function area_intersects_local_domain(this,area)
		class(computational_domain)	,intent(in)	:: this
		integer	,dimension(3,2)		,intent(in)	:: area

		logical	:: flag1, flag2
		integer	:: dim
		
		flag1 = .false.
		flag2 = .false.
		
		do dim = 1, this%dimensions
		
			if ((1 + this%global_offset(dim) < area(dim,1)) .and. &
			    (this%cells_number_decomposed(dim) + this%global_offset(dim) > area(dim,1))) then
				flag1 = .true.
			end if

			if ((1 + this%global_offset(dim) < area(dim,2)) .and. &
			    (this%cells_number_decomposed(dim) + this%global_offset(dim) > area(dim,2))) then
				flag2 = .true.
			end if			
			
		end do
		
		area_intersects_local_domain = .false.
		if(flag1.and.flag2) then
			area_intersects_local_domain = .true.
		end if
	end function
	
	
end module
