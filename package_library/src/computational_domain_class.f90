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
        character(len=20)	,dimension(:)   ,allocatable    :: mesh_types				! Mesh types (uniform/linear/exponential)
        real(dkind)         ,dimension(:,:) ,allocatable	:: reference_coordinates    ! Reference coordinates for cell size tendency change
        real(dkind)         ,dimension(:,:) ,allocatable	:: reference_cell_sizes		! Cell sizes in reference coordinates
        integer				,dimension(:)	,allocatable	:: reference_coordinates_number	! Reference coordinates number
        integer				,dimension(:,:) ,allocatable	:: local_mesh_cells_number	! Local mesh cells number between reference coordinates

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
        procedure				:: get_mesh_types
        procedure				:: get_reference_coordinates
        procedure				:: get_reference_cell_sizes
        procedure				:: get_reference_coordinates_number
        procedure				:: get_local_mesh_cells_number
		
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

	type(computational_domain)	function constructor(dimensions,cells_number,coordinate_system,lengths,axis_names,mesh_types,reference_coordinates,reference_coordinates_number,reference_cell_sizes)
		integer								,intent(in)				:: dimensions
		integer				,dimension(3)	,intent(in), optional	:: cells_number
		character(len=*)					,intent(in)             :: coordinate_system
        real(dkind)			,dimension(:,:)	,intent(in)				:: lengths
		character(len=*)	,dimension(:)	,intent(in)             :: axis_names
        character(len=20)	,dimension(:)	,intent(in), optional	:: mesh_types
        integer				,dimension(:)	,intent(in), optional	:: reference_coordinates_number
        real(dkind)			,dimension(:,:)	,intent(in), optional	:: reference_cell_sizes
        real(dkind)			,dimension(:,:)	,intent(in), optional	:: reference_coordinates

		integer	:: io_unit
        
        if ((present(reference_cell_sizes) .and. .not.present(mesh_types)) .or. (present(mesh_types) .and. .not.present(reference_cell_sizes))) then
            print *, "Error in problem domain set properies: mesh_types and reference_cell_sizes must be present or absent simultaneously"
            pause
        end if
        
        if  ((present(reference_cell_sizes) .and. .not.present(reference_coordinates)) .or. (present(reference_coordinates) .and. .not.present(reference_cell_sizes))) then
            print *, "Error in problem domain set properies: reference_coordinates and reference_cell_sizes must be present or absent simultaneously"
            pause
        end if
        
        if  ((present(reference_coordinates_number) .and. .not.present(reference_cell_sizes))) then
            print *, "Error in problem domain set properies: if reference_coordinates_number is present also must present reference_cell_sizes and reference_coordinates"
            pause
        end if
        
        if ((present(reference_cell_sizes) .and. present(cells_number)) .or. (.not.present(reference_cell_sizes) .and. .not.present(cells_number))) then
            print *, "Error in problem domain set properies: reference_cell_sizes and cells_number can't be present or absent simultaneously"
            pause
        end if
		
		call constructor%set_properties(dimensions,cells_number,lengths,coordinate_system,axis_names,mesh_types,reference_coordinates_number,reference_cell_sizes,reference_coordinates)
		
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
        character(len=20)	,dimension(:)	,allocatable	:: mesh_types
        integer				,dimension(:)	,allocatable	:: reference_coordinates_number
        real(dkind)			,dimension(:,:)	,allocatable	:: reference_cell_sizes
        real(dkind)			,dimension(:,:)	,allocatable	:: reference_coordinates

		namelist /domain_properties_1/ dimensions, cells_number, coordinate_system
		namelist /domain_properties_2/ lengths, axis_names, mesh_types, reference_coordinates_number
        namelist /domain_properties_3/ reference_cell_sizes, reference_coordinates
		
		read(unit = domain_data_unit, nml = domain_properties_1)
		
		allocate(lengths(dimensions,2))
		allocate(axis_names(dimensions))
        allocate(mesh_types(dimensions))
        allocate(reference_coordinates_number(dimensions))
		
		read(unit = domain_data_unit, nml = domain_properties_2)
        
        allocate(reference_cell_sizes(dimensions,maxval(reference_coordinates_number)))
        allocate(reference_coordinates(dimensions,maxval(reference_coordinates_number)))
        
        read(unit = domain_data_unit, nml = domain_properties_3)
	
		call this%set_properties(dimensions,cells_number,lengths,coordinate_system,axis_names,mesh_types,reference_coordinates_number,reference_cell_sizes,reference_coordinates)
	end subroutine
	
	subroutine write_properties(this,domain_data_unit)
		class(computational_domain)	,intent(in)	:: this
		integer						,intent(in)	:: domain_data_unit
	
		integer								:: dimensions
		integer				,dimension(3)	:: cells_number
		character(len=20)					:: coordinate_system
		real(dkind)			,dimension(:,:)	,allocatable	:: lengths
		character(len=5)	,dimension(:)	,allocatable	:: axis_names
        character(len=20)	,dimension(:)	,allocatable	:: mesh_types
        integer				,dimension(:)	,allocatable	:: reference_coordinates_number
        real(dkind)			,dimension(:,:)	,allocatable	:: reference_cell_sizes
        real(dkind)			,dimension(:,:)	,allocatable	:: reference_coordinates

		namelist /domain_properties_1/ dimensions, cells_number, coordinate_system
		namelist /domain_properties_2/ lengths, axis_names, mesh_types, reference_coordinates_number
        namelist /domain_properties_3/ reference_cell_sizes, reference_coordinates
		
		dimensions			= this%dimensions
		cells_number		= this%cells_number
		coordinate_system	= this%coordinate_system
		
		allocate(axis_names(dimensions))
		allocate(lengths(dimensions,2))
        allocate(mesh_types(dimensions))
        allocate(reference_coordinates_number(dimensions))
		
		axis_names          = this%axis_names
		lengths				= this%lengths
        mesh_types			= this%mesh_types
        reference_coordinates_number = this%reference_coordinates_number
        
        allocate(reference_cell_sizes(dimensions,maxval(reference_coordinates_number)))
        allocate(reference_coordinates(dimensions,maxval(reference_coordinates_number)))
        
        reference_cell_sizes	= this%reference_cell_sizes
        reference_coordinates	= this%reference_coordinates
		
		write(unit = domain_data_unit, nml = domain_properties_1)
		write(unit = domain_data_unit, nml = domain_properties_2)
        write(unit = domain_data_unit, nml = domain_properties_3)
	end subroutine
	
	subroutine set_properties(this,dimensions,cells_number,lengths,coordinate_system,axis_names,mesh_types,reference_coordinates_number,reference_cell_sizes,reference_coordinates)
		class(computational_domain)			,intent(inout)	:: this
		integer								,intent(in)				:: dimensions
		integer				,dimension(3)	,intent(in), optional	:: cells_number
		character(len=*)					,intent(in)             :: coordinate_system
        real(dkind)			,dimension(:,:)	,intent(in)				:: lengths
		character(len=*)	,dimension(:)	,intent(in)             :: axis_names
        character(len=20)	,dimension(:)	,intent(in), optional	:: mesh_types
        integer				,dimension(:)	,intent(in), optional	:: reference_coordinates_number
        real(dkind)			,dimension(:,:)	,intent(in), optional	:: reference_cell_sizes
        real(dkind)			,dimension(:,:)	,intent(in), optional	:: reference_coordinates
		
		real(dkind) 	:: domain_side_length, domain_local_side_length, cell_size_koef
		
		integer			:: dim, i, j
		
		allocate(this%axis_names(dimensions))
		allocate(this%lengths(dimensions,2))
        allocate(this%mesh_types(dimensions))
        allocate(this%reference_coordinates_number(dimensions))

		this%dimensions				= dimensions
		this%lengths				= lengths(1:dimensions,:)
		this%coordinate_system		= coordinate_system
		this%axis_names				= axis_names(1:dimensions)
        
        if (present(cells_number)) then
			this%cells_number			= cells_number
        else
            this%cells_number			= 0
        end if
        
        if (present(mesh_types)) then
            this%mesh_types = mesh_types
        else
            this%mesh_types = 'uniform'
        end if
        
        do dim = 1,dimensions
            if (.not.(this%mesh_types(dim) == 'uniform' .or. this%mesh_types(dim) == 'linear' .or. this%mesh_types(dim) == 'exponential')) then
                print *, "Error in problem domain set properies: unknown mesh type ", this%mesh_types(dim)
				pause
            end if
        end do
        
        if (present(reference_coordinates_number)) then
            this%reference_coordinates_number = reference_coordinates_number
        else
            this%reference_coordinates_number = 2
        end if
        
        do dim = 1,dimensions
            if (this%mesh_types(dim) == 'uniform') then
                this%reference_coordinates_number(dim) = 2
            else
                if (this%reference_coordinates_number(dim) < 2) then
                    print *, "Error in problem domain set properies: wrong reference coordinates number for dim = ", dim
					pause
                end if
            end if
        end do
        
        allocate(this%reference_cell_sizes(dimensions,maxval(this%reference_coordinates_number)))
        allocate(this%reference_coordinates(dimensions,maxval(this%reference_coordinates_number)))
        
        if (present(reference_cell_sizes)) then
            this%reference_cell_sizes(:,:)	= reference_cell_sizes(:,:maxval(this%reference_coordinates_number))
            this%reference_coordinates(:,:) = reference_coordinates(:,:maxval(this%reference_coordinates_number))
            do dim = 1,dimensions
				this%reference_coordinates(dim,1) = this%lengths(dim,1)
                this%reference_coordinates(dim,this%reference_coordinates_number(dim)) = this%lengths(dim,2)
            end do
        else
            this%reference_cell_sizes = 0.0_dkind
            do dim = 1,dimensions
                domain_side_length = this%lengths(dim,2) - this%lengths(dim,1)
				this%reference_cell_sizes(dim,:) = domain_side_length / this%cells_number(dim)
                this%reference_coordinates(dim,1) = this%lengths(dim,1)
                this%reference_coordinates(dim,2) = this%lengths(dim,2)
			end do
        end if
        
        allocate(this%local_mesh_cells_number(dimensions,maxval(this%reference_coordinates_number)-1))
        this%local_mesh_cells_number = 0
        
        do dim = 1,dimensions
            domain_side_length = this%lengths(dim,2) - this%lengths(dim,1)
            select case (this%mesh_types(dim))
            case('uniform')
                if (this%cells_number(dim) == 0) this%cells_number(dim) = ceiling(domain_side_length / this%reference_cell_sizes(dim,1))
                
            case('linear')
                this%cells_number(dim) = 0
                do i = 1,this%reference_coordinates_number(dim)-1
                    domain_local_side_length = this%reference_coordinates(dim,i+1) - this%reference_coordinates(dim,i)
                    if (domain_local_side_length > domain_side_length * 1.01_dkind) then
                        print *, "Error in problem domain set properies: check reference coordinates for dim = ", dim
						pause
					end if
					this%local_mesh_cells_number(dim,i) = ceiling(2.0_dkind*domain_local_side_length/(this%reference_cell_sizes(dim,i)+this%reference_cell_sizes(dim,i+1)))
                    this%cells_number(dim) = this%cells_number(dim) + this%local_mesh_cells_number(dim,i)
                end do
                
            case('exponential')
                this%cells_number(dim) = 0
                do i = 1,this%reference_coordinates_number(dim)-1
                    domain_local_side_length = this%reference_coordinates(dim,i+1) - this%reference_coordinates(dim,i)
                    if (domain_local_side_length > domain_side_length * 1.01_dkind) then
                        print *, "Error in problem domain set properies: check reference coordinates for dim = ", dim
						pause
                    end if
                    cell_size_koef = (domain_local_side_length - this%reference_cell_sizes(dim,i)) / (domain_local_side_length - this%reference_cell_sizes(dim,i+1))
                    if (abs(cell_size_koef - 1.0_dkind) < 1.0e-50_dkind) then
                        this%local_mesh_cells_number(dim,i) = ceiling(domain_local_side_length/this%reference_cell_sizes(dim,i))
                    else
                        this%local_mesh_cells_number(dim,i) = 1 + ceiling(log(this%reference_cell_sizes(dim,i+1) / this%reference_cell_sizes(dim,i)) / log(cell_size_koef))
                    end if
                    this%cells_number(dim) = this%cells_number(dim) + this%local_mesh_cells_number(dim,i)
                end do                
            end select
        end do
        
		this%faces_number			= this%cells_number + 1

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
    
    pure function get_mesh_types(this)
		class(computational_domain)			,intent(in)		:: this
		character(len=20)	,dimension(:)	,allocatable	:: get_mesh_types
        
        allocate(get_mesh_types,source = this%mesh_types)
	end function
    
    pure function get_reference_coordinates(this)
		class(computational_domain)			,intent(in)		:: this
        real(dkind)         ,dimension(:,:) ,allocatable	:: get_reference_coordinates
		
		allocate(get_reference_coordinates,source = this%reference_coordinates)
    end function
    
    pure function get_reference_cell_sizes(this)
		class(computational_domain)			,intent(in)		:: this
        real(dkind)         ,dimension(:,:) ,allocatable	:: get_reference_cell_sizes
		
        allocate(get_reference_cell_sizes,source = this%reference_cell_sizes)
    end function
    
    pure function get_reference_coordinates_number(this)
		class(computational_domain)			,intent(in)		:: this
        integer		       ,dimension(:)    ,allocatable	:: get_reference_coordinates_number
		
        allocate(get_reference_coordinates_number,source = this%reference_coordinates_number)
    end function
    
    pure function get_local_mesh_cells_number(this)
		class(computational_domain)			,intent(in)		:: this
        integer             ,dimension(:,:) ,allocatable	:: get_local_mesh_cells_number
		
        allocate(get_local_mesh_cells_number,source = this%local_mesh_cells_number)
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
