module mpi_communications_class

	use kind_parameters
	use global_data
	use computational_domain_class
	use field_pointers
    use boundary_conditions_class
	use computational_mesh_class
	
#ifdef mpi
	use MPI
#endif

    implicit none
    
    private
	public	:: mpi_communications, mpi_communications_c

	type mpi_communications
        private
		type(computational_domain)	:: domain

        integer ,dimension(-1:1,-1:1,-1:1)   :: processor_neigbours

        integer ,dimension(-1:1,-1:1,-1:1)    :: cons_send_subarray_type
        integer ,dimension(-1:1,-1:1,-1:1)    :: cons_recv_subarray_type    

        integer ,dimension(-1:1,-1:1,-1:1)    :: flow_send_subarray_type_longitudinal
        integer ,dimension(-1:1,-1:1,-1:1)    :: flow_recv_subarray_type_longitudinal

        integer ,dimension(-1:1,-1:1,-1:1)    :: flow_send_subarray_type_transverse
        integer ,dimension(-1:1,-1:1,-1:1)    :: flow_recv_subarray_type_transverse     
        
        real(dkind) ,dimension(:,:,:)   ,allocatable    :: cons_buffer
        real(dkind) ,dimension(:,:,:)   ,allocatable    :: flow_buffer

	contains
		procedure   :: exchange_conservative_scalar_field
        procedure   :: exchange_conservative_vector_field
        procedure   :: exchange_conservative_tensor_field
        procedure   :: exchange_boundary_conditions_markers 
		procedure   :: exchange_mesh
		procedure	:: exchange_flow_scalar_field
        procedure   :: exchange_flow_vector_field
	end type

	interface	mpi_communications_c
		module procedure	constructor
	end interface

contains

    type(mpi_communications)    function constructor(domain)
        type(computational_domain)     ,intent(in)	:: domain

        integer                 :: dim, dim1, dim2, dim3
        integer ,dimension(3)   :: neighbour_coords
        integer                 :: sign1, sign2


        integer ,dimension(3)   :: domain_size

        integer,dimension(-1:1,-1:1,-1:1,3) :: subarray_starts_recv, subarray_starts_send
        integer,dimension(-1:1,-1:1,-1:1,3) :: subarray_subsizes

        integer                 :: mpi_communicator
        integer ,dimension(3)   :: processor_number, local_cells_number, local_faces_number, processor_grid_coord
        integer ,dimension(3,2) :: cons_allocation_bounds, flow_allocation_bounds


        integer :: type
        integer :: error

        constructor%domain	= domain

        processor_grid_coord    = domain%get_processor_grid_coord()
        processor_number        = domain%get_processor_number()
        local_cells_number      = domain%get_local_cells_number()
        local_faces_number      = domain%get_local_faces_number()
        mpi_communicator        = domain%get_mpi_communicator()
		cons_allocation_bounds	    = domain%get_local_utter_cells_bounds()
        flow_allocation_bounds	    = domain%get_local_utter_faces_bounds()

#ifdef mpi

		allocate(constructor%cons_buffer(	cons_allocation_bounds(1,1):cons_allocation_bounds(1,2)	, &
										    cons_allocation_bounds(2,1):cons_allocation_bounds(2,2)	, &
										    cons_allocation_bounds(3,1):cons_allocation_bounds(3,2)))

        allocate(constructor%flow_buffer(	flow_allocation_bounds(1,1):flow_allocation_bounds(1,2)	, &
                                            flow_allocation_bounds(2,1):flow_allocation_bounds(2,2)	, &
                                            flow_allocation_bounds(3,1):flow_allocation_bounds(3,2)))                                    

        !*********************************** Get neighbours ******************************************

        constructor%processor_neigbours = -1

        do dim1 = -1,1
        do dim2 = -1,1
        do dim3 = -1,1
            neighbour_coords(1) = processor_grid_coord(1) + dim1
            neighbour_coords(2) = processor_grid_coord(2) + dim2
            neighbour_coords(3) = processor_grid_coord(3) + dim3
            if  (( all(neighbour_coords < processor_number)).and.( all(neighbour_coords >= 0)).and.((dim1 /= 0).or.(dim2 /= 0).or.(dim3 /= 0))) then
                call mpi_cart_rank(mpi_communicator, neighbour_coords, constructor%processor_neigbours(dim1,dim2,dim3),error)
            end if          
        end do
        end do
        end do

        !*********************************** Commit datatypes ******************************************

        domain_size = 1
        domain_size = local_cells_number

        do dim1 = -1,1
        do dim2 = -1,1
        do dim3 = -1,1
            ! 1 --> 1 ; -1 --> 1 ; 0 --> 
            subarray_subsizes(dim1,dim2,dim3,1) = (1 - abs(dim1))*abs(local_cells_number(1) - 2) + abs(dim1)
            subarray_subsizes(dim1,dim2,dim3,2) = (1 - abs(dim2))*abs(local_cells_number(2) - 2) + abs(dim2)
            subarray_subsizes(dim1,dim2,dim3,3) = (1 - abs(dim3))*abs(local_cells_number(3) - 2) + abs(dim3)

            if (local_cells_number(1) /= 1) then
                subarray_starts_recv(dim1,dim2,dim3,1) = 0.5*dim1*(dim1 + 1)*(local_cells_number(1) - 1)    + (1 - abs(dim1))
            else
                subarray_starts_recv(dim1,dim2,dim3,1) = 0
            end if

            if (local_cells_number(2) /= 1) then
                subarray_starts_recv(dim1,dim2,dim3,2) = 0.5*dim2*(dim2 + 1)*(local_cells_number(2) - 1)    + (1 - abs(dim2))
            else
                subarray_starts_recv(dim1,dim2,dim3,2) = 0
            end if

            if (local_cells_number(3) /= 1) then
                subarray_starts_recv(dim1,dim2,dim3,3) = 0.5*dim3*(dim3 + 1)*(local_cells_number(3) - 1)    + (1 - abs(dim3))
            else 
                subarray_starts_recv(dim1,dim2,dim3,3) = 0
            end if

            subarray_starts_send(dim1,dim2,dim3,:) = subarray_starts_recv(dim1,dim2,dim3,:) - [dim1,dim2,dim3]

            if  (constructor%processor_neigbours(dim1,dim2,dim3) >= 0) then
                call mpi_type_create_subarray(3,domain_size,subarray_subsizes(dim1,dim2,dim3,:),subarray_starts_send(dim1,dim2,dim3,:), MPI_ORDER_FORTRAN, MPI_DOUBLE_PRECISION,constructor%cons_send_subarray_type(dim1,dim2,dim3), error)
                call mpi_type_create_subarray(3,domain_size,subarray_subsizes(dim1,dim2,dim3,:),subarray_starts_recv(dim1,dim2,dim3,:), MPI_ORDER_FORTRAN, MPI_DOUBLE_PRECISION,constructor%cons_recv_subarray_type(dim1,dim2,dim3), error)
                call mpi_type_commit(constructor%cons_send_subarray_type(dim1,dim2,dim3),error)
                call mpi_type_commit(constructor%cons_recv_subarray_type(dim1,dim2,dim3),error)
            end if

        end do
        end do
        end do


        domain_size = 1
        domain_size = local_faces_number

        do dim1 = -1,1
        do dim2 = -1,1
        do dim3 = -1,1

        ! 1 --> 1 ; -1 --> 1 ; 0 --> 
           subarray_subsizes(dim1,dim2,dim3,1) = (1 - abs(dim1))*abs(local_cells_number(1) - 2) + abs(dim1)        ! сонаправленные грани
           subarray_subsizes(dim1,dim2,dim3,2) = (1 - abs(dim2))*abs(local_cells_number(2) - 2) + abs(dim2)
           subarray_subsizes(dim1,dim2,dim3,3) = (1 - abs(dim3))*abs(local_cells_number(3) - 2) + abs(dim3)

            if (local_faces_number(1) /= 1) then
                subarray_starts_recv(dim1,dim2,dim3,1) = 0.5*dim1*(dim1 + 1)*(local_faces_number(1) - 1)    + (1 - abs(dim1))   ! сонаправленные грани
            else
                subarray_starts_recv(dim1,dim2,dim3,1) = 0
            end if

            if (local_faces_number(2) /= 1) then
                subarray_starts_recv(dim1,dim2,dim3,2) = 0.5*dim2*(dim2 + 1)*(local_faces_number(2) - 1)    + (1 - abs(dim2))   ! сонаправленные грани
            else
                subarray_starts_recv(dim1,dim2,dim3,2) = 0
            end if

            if (local_faces_number(3) /= 1) then
                subarray_starts_recv(dim1,dim2,dim3,3) = 0.5*dim3*(dim3 + 1)*(local_faces_number(3) - 1)    + (1 - abs(dim3))   ! сонаправленные грани
            else 
                subarray_starts_recv(dim1,dim2,dim3,3) = 0
            end if

           subarray_starts_send(dim1,dim2,dim3,:) = subarray_starts_recv(dim1,dim2,dim3,:) - 2*[dim1,dim2,dim3]    ! сонаправленные грани

            if  ((constructor%processor_neigbours(dim1,dim2,dim3) >= 0).and.(abs(dim1+dim2+dim3) == 1)) then
                call mpi_type_create_subarray(3,domain_size,subarray_subsizes(dim1,dim2,dim3,:),subarray_starts_send(dim1,dim2,dim3,:), MPI_ORDER_FORTRAN, MPI_DOUBLE_PRECISION,constructor%flow_send_subarray_type_longitudinal(dim1,dim2,dim3), error)
                call mpi_type_create_subarray(3,domain_size,subarray_subsizes(dim1,dim2,dim3,:),subarray_starts_recv(dim1,dim2,dim3,:), MPI_ORDER_FORTRAN, MPI_DOUBLE_PRECISION,constructor%flow_recv_subarray_type_longitudinal(dim1,dim2,dim3), error)
                call mpi_type_commit(constructor%flow_send_subarray_type_longitudinal(dim1,dim2,dim3),error)
                call mpi_type_commit(constructor%flow_recv_subarray_type_longitudinal(dim1,dim2,dim3),error)
            end if

        end do
        end do
        end do

        do dim1 = -1,1
        do dim2 = -1,1
        do dim3 = -1,1

        ! 1 --> 1 ; -1 --> 1 ; 0 --> 

            subarray_subsizes(dim1,dim2,dim3,1) = (1 - abs(dim1))*abs(local_faces_number(1) - 2) + abs(dim1)        ! поперечные грани
            subarray_subsizes(dim1,dim2,dim3,2) = (1 - abs(dim2))*abs(local_faces_number(2) - 2) + abs(dim2)
            subarray_subsizes(dim1,dim2,dim3,3) = (1 - abs(dim3))*abs(local_faces_number(3) - 2) + abs(dim3)           

            if (local_faces_number(1) /= 1) then
                subarray_starts_recv(dim1,dim2,dim3,1) = 0.5*dim1*(dim1 + 1)*(local_faces_number(1) - 2)    + (1 - abs(dim1))   ! поперечные грани
            else
                subarray_starts_recv(dim1,dim2,dim3,1) = 0
            end if

            if (local_faces_number(2) /= 1) then
                subarray_starts_recv(dim1,dim2,dim3,2) = 0.5*dim2*(dim2 + 1)*(local_faces_number(2) - 2)    + (1 - abs(dim2))   ! поперечные грани
            else
                subarray_starts_recv(dim1,dim2,dim3,2) = 0
            end if

            if (local_faces_number(3) /= 1) then
                subarray_starts_recv(dim1,dim2,dim3,3) = 0.5*dim3*(dim3 + 1)*(local_faces_number(3) - 2)    + (1 - abs(dim3))   ! поперечные грани
            else 
                subarray_starts_recv(dim1,dim2,dim3,3) = 0
            end if

            subarray_starts_send(dim1,dim2,dim3,:) = subarray_starts_recv(dim1,dim2,dim3,:) - 1*[dim1,dim2,dim3]    ! поперечные грани

            if  ((constructor%processor_neigbours(dim1,dim2,dim3) >= 0).and.(abs(dim1+dim2+dim3) == 1)) then
                call mpi_type_create_subarray(3,domain_size,subarray_subsizes(dim1,dim2,dim3,:),subarray_starts_send(dim1,dim2,dim3,:), MPI_ORDER_FORTRAN, MPI_DOUBLE_PRECISION,constructor%flow_send_subarray_type_transverse(dim1,dim2,dim3), error)
                call mpi_type_create_subarray(3,domain_size,subarray_subsizes(dim1,dim2,dim3,:),subarray_starts_recv(dim1,dim2,dim3,:), MPI_ORDER_FORTRAN, MPI_DOUBLE_PRECISION,constructor%flow_recv_subarray_type_transverse(dim1,dim2,dim3), error)
                call mpi_type_commit(constructor%flow_send_subarray_type_transverse(dim1,dim2,dim3),error)
                call mpi_type_commit(constructor%flow_recv_subarray_type_transverse(dim1,dim2,dim3),error)
            end if

        end do
        end do
        end do
#endif

    end function

    subroutine exchange_conservative_scalar_field(this, scal_ptr)
        class(mpi_communications)   ,target     ,intent(inout)  :: this
        type(field_scalar_cons)	    ,target     ,intent(inout)  :: scal_ptr

        real(dkind) ,pointer    ,dimension(:) :: buffer_ptr

        integer ,dimension(3)   :: processor_number, processor_grid_coord
        
        integer                                                :: dim1,dim2,dim3
        integer ,dimension(3)                                  :: neighbour_coords
        integer                                                :: error
#ifdef mpi
        integer ,dimension(-1:1,-1:1,-1:1,mpi_status_size)     :: sr_status
#endif        

#ifdef mpi

        processor_grid_coord    = this%domain%get_processor_grid_coord()
        processor_number        = this%domain%get_processor_number()

        this%cons_buffer = scal_ptr%cells

        nullify(buffer_ptr)
        buffer_ptr(0:size(scal_ptr%cells,1)*size(scal_ptr%cells,2)*size(scal_ptr%cells,3)-1) => this%cons_buffer

        do dim1 = -1,1
        do dim2 = -1,1
        do dim3 = -1,1
            if  (this%processor_neigbours(dim1,dim2,dim3) >= 0) then
                call mpi_sendrecv(  buffer_ptr  ,1  ,this%cons_send_subarray_type(dim1,dim2,dim3)    ,this%processor_neigbours(dim1,dim2,dim3)  , 1  &
                                ,   buffer_ptr  ,1  ,this%cons_recv_subarray_type(dim1,dim2,dim3)    ,this%processor_neigbours(dim1,dim2,dim3)  , 1  &
                                ,   MPI_COMM_WORLD,MPI_STATUS_IGNORE,error)
            end if 
        end do
        end do
        end do

        scal_ptr%cells = this%cons_buffer

        call MPI_BARRIER(MPI_COMM_WORLD,error)

#endif
    end subroutine

    subroutine exchange_boundary_conditions_markers(this, bc_ptr)
        class(mpi_communications)               ,target     ,intent(inout)  :: this
        type(boundary_conditions)       	    ,target     ,intent(inout)  :: bc_ptr

        real(dkind) ,pointer    ,dimension(:) :: buffer_ptr

        integer ,dimension(3)   :: processor_number, processor_grid_coord
        
        integer                                                :: dim1,dim2,dim3
        integer ,dimension(3)                                  :: neighbour_coords
        integer                                                :: error
#ifdef mpi
        integer ,dimension(-1:1,-1:1,-1:1,mpi_status_size)     :: sr_status
#endif        

#ifdef mpi

        processor_grid_coord    = this%domain%get_processor_grid_coord()
        processor_number        = this%domain%get_processor_number()

        this%cons_buffer = bc_ptr%bc_markers

        nullify(buffer_ptr)
        buffer_ptr(0:size(bc_ptr%bc_markers,1)*size(bc_ptr%bc_markers,2)*size(bc_ptr%bc_markers,3)-1) => this%cons_buffer

        do dim1 = -1,1
        do dim2 = -1,1
        do dim3 = -1,1
            if  (this%processor_neigbours(dim1,dim2,dim3) >= 0) then
                call mpi_sendrecv(  buffer_ptr  ,1  ,this%cons_send_subarray_type(dim1,dim2,dim3)    ,this%processor_neigbours(dim1,dim2,dim3)  , 1  &
                                ,   buffer_ptr  ,1  ,this%cons_recv_subarray_type(dim1,dim2,dim3)    ,this%processor_neigbours(dim1,dim2,dim3)  , 1  &
                                ,   MPI_COMM_WORLD,MPI_STATUS_IGNORE,error)
            end if 
        end do
        end do
        end do

        bc_ptr%bc_markers = this%cons_buffer

        call MPI_BARRIER(MPI_COMM_WORLD,error)

#endif
    end subroutine
	
    subroutine exchange_mesh(this, mesh_ptr)
        class(mpi_communications)   ,target     ,intent(inout)  :: this
        type(computational_mesh)       			,target     ,intent(inout)  :: mesh_ptr

        real(dkind) ,pointer    ,dimension(:) :: buffer_ptr

        integer                 :: dimensions
        integer ,dimension(3)   :: processor_number, processor_grid_coord
        
        integer                                                :: dim,dim1,dim2,dim3
        integer ,dimension(3)                                  :: neighbour_coords
        integer                                                :: error
#ifdef mpi
        integer ,dimension(-1:1,-1:1,-1:1,mpi_status_size)     :: sr_status
#endif

#ifdef mpi

        dimensions              = this%domain%get_domain_dimensions()
        processor_grid_coord    = this%domain%get_processor_grid_coord()
        processor_number        = this%domain%get_processor_number()

        do dim = 1, dimensions

            this%cons_buffer = mesh_ptr%mesh(dim,:,:,:)

            nullify(buffer_ptr)
            buffer_ptr(0:size(mesh_ptr%mesh,2)*size(mesh_ptr%mesh,3)*size(mesh_ptr%mesh,4)-1) => this%cons_buffer

            do dim1 = -1,1
			do dim2 = -1,1
			do dim3 = -1,1
				if  (this%processor_neigbours(dim1,dim2,dim3) >= 0) then
					call mpi_sendrecv(  buffer_ptr  ,1  ,this%cons_send_subarray_type(dim1,dim2,dim3)    ,this%processor_neigbours(dim1,dim2,dim3)  , 1  &
									,   buffer_ptr  ,1  ,this%cons_recv_subarray_type(dim1,dim2,dim3)    ,this%processor_neigbours(dim1,dim2,dim3)  , 1  &
									,   MPI_COMM_WORLD,MPI_STATUS_IGNORE,error)
				end if 
			end do
			end do
			end do

            mesh_ptr%mesh(dim,:,:,:) = this%cons_buffer
        end do

        call MPI_BARRIER(MPI_COMM_WORLD,error)

#endif
    end subroutine
	
    subroutine exchange_conservative_vector_field(this, vect_ptr)
        class(mpi_communications)   ,target     ,intent(inout)  :: this
        type(field_vector_cons)	    ,target     ,intent(inout)  :: vect_ptr

        real(dkind) ,pointer    ,dimension(:) :: buffer_ptr

        integer                 :: vector_projections_number
        integer ,dimension(3)   :: processor_number, processor_grid_coord

        integer                                                :: dim,dim1,dim2,dim3
        integer ,dimension(3)                                  :: neighbour_coords
        integer                                                :: error
#ifdef mpi
        integer ,dimension(-1:1,-1:1,-1:1,mpi_status_size)     :: sr_statuses
#endif        

#ifdef mpi

        processor_grid_coord    = this%domain%get_processor_grid_coord()
        processor_number        = this%domain%get_processor_number()

        vector_projections_number = vect_ptr%get_projections_number()

        do dim = 1, vector_projections_number

            this%cons_buffer = vect_ptr%pr(dim)%cells

            nullify(buffer_ptr)
            buffer_ptr(0:size(vect_ptr%pr(dim)%cells,1)*size(vect_ptr%pr(dim)%cells,2)*size(vect_ptr%pr(dim)%cells,3)-1) => this%cons_buffer

            do dim1 = -1,1
            do dim2 = -1,1
            do dim3 = -1,1
                if  (this%processor_neigbours(dim1,dim2,dim3) >= 0) then
                    call mpi_sendrecv(  buffer_ptr  ,1  ,this%cons_send_subarray_type(dim1,dim2,dim3)    ,this%processor_neigbours(dim1,dim2,dim3)  , 1  &
                                    ,   buffer_ptr  ,1  ,this%cons_recv_subarray_type(dim1,dim2,dim3)    ,this%processor_neigbours(dim1,dim2,dim3)  , 1  &
                                    ,   MPI_COMM_WORLD,MPI_STATUS_IGNORE,error)
                end if 
            end do
            end do
            end do

            vect_ptr%pr(dim)%cells = this%cons_buffer
        end do

        call MPI_BARRIER(MPI_COMM_WORLD,error)

#endif
    end subroutine

    subroutine exchange_conservative_tensor_field(this, tens_ptr)
        class(mpi_communications)   ,target     ,intent(inout)  :: this
        type(field_tensor)	        ,target     ,intent(inout)  :: tens_ptr

        real(dkind) ,pointer    ,dimension(:) :: buffer_ptr

        integer ,dimension(2)   :: tensor_projections_number
        integer ,dimension(3)   :: processor_number, processor_grid_coord

        integer                                                :: dim,dim1,dim2,dim3,dim4,dim5
        integer ,dimension(3)                                  :: neighbour_coords
        integer                                                :: error
#ifdef mpi
        integer ,dimension(-1:1,-1:1,-1:1,mpi_status_size)     :: sr_status
#endif        

#ifdef mpi

        processor_grid_coord    = this%domain%get_processor_grid_coord()
        processor_number        = this%domain%get_processor_number()

        tensor_projections_number = tens_ptr%get_projections_number()

        do dim4 = 1, tensor_projections_number(1)
        do dim5 = 1, tensor_projections_number(2)

            this%cons_buffer = tens_ptr%pr(dim4,dim5)%cells

            nullify(buffer_ptr)
            buffer_ptr(0:size(tens_ptr%pr(dim4,dim5)%cells,1)*size(tens_ptr%pr(dim4,dim5)%cells,2)*size(tens_ptr%pr(dim4,dim5)%cells,3)-1) => this%cons_buffer

            do dim1 = -1,1
            do dim2 = -1,1
            do dim3 = -1,1
                if  (this%processor_neigbours(dim1,dim2,dim3) >= 0) then
                    call mpi_sendrecv(  buffer_ptr  ,1  ,this%cons_send_subarray_type(dim1,dim2,dim3)    ,this%processor_neigbours(dim1,dim2,dim3)  , 1  &
                                    ,   buffer_ptr  ,1  ,this%cons_recv_subarray_type(dim1,dim2,dim3)    ,this%processor_neigbours(dim1,dim2,dim3)  , 1  &
                                    ,   MPI_COMM_WORLD,MPI_STATUS_IGNORE,error)
                end if 
            end do
            end do
            end do

            tens_ptr%pr(dim4,dim5)%cells = this%cons_buffer
        end do
        end do

        call MPI_BARRIER(MPI_COMM_WORLD,error)

#endif
    end subroutine

    subroutine exchange_flow_scalar_field(this, scal_ptr)
        class(mpi_communications)   ,target     ,intent(inout)  :: this
        type(field_scalar_flow)	    ,target     ,intent(inout)  :: scal_ptr

        real(dkind) ,pointer    ,dimension(:) :: buffer_ptr

        integer                 :: dimensions
        integer ,dimension(3)   :: processor_number, processor_grid_coord
        
        integer                                                :: dim,dim1,dim2,dim3
        integer ,dimension(3)                                  :: neighbour_coords
        integer                                                :: error
#ifdef mpi
        integer ,dimension(-1:1,-1:1,-1:1,mpi_status_size)     :: sr_status
#endif

#ifdef mpi

        dimensions              = this%domain%get_domain_dimensions()
        processor_grid_coord    = this%domain%get_processor_grid_coord()
        processor_number        = this%domain%get_processor_number()

        do dim = 1, dimensions

            this%flow_buffer = scal_ptr%cells(dim,:,:,:)

            nullify(buffer_ptr)
            buffer_ptr(0:size(scal_ptr%cells,2)*size(scal_ptr%cells,3)*size(scal_ptr%cells,4)-1) => this%flow_buffer

            do dim1 = -1,1
            do dim2 = -1,1
            do dim3 = -1,1
                if  ((this%processor_neigbours(dim1,dim2,dim3) >= 0).and.(abs(dim1)+abs(dim2)+abs(dim3) == 1).and.(abs(I_m(dim,1)*dim1 + I_m(dim,2)*dim2 +I_m(dim,3)*dim3) == 1)) then  ! сонаправленные грани
                    call mpi_sendrecv(  buffer_ptr  ,1  ,this%flow_send_subarray_type_longitudinal(dim1,dim2,dim3)    ,this%processor_neigbours(dim1,dim2,dim3)  , 1  &
                                    ,   buffer_ptr  ,1  ,this%flow_recv_subarray_type_longitudinal(dim1,dim2,dim3)    ,this%processor_neigbours(dim1,dim2,dim3)  , 1  &
                                    ,   MPI_COMM_WORLD,MPI_STATUS_IGNORE,error)
                end if 
            end do
            end do
            end do

            do dim1 = -1,1
            do dim2 = -1,1
            do dim3 = -1,1
                if  ((this%processor_neigbours(dim1,dim2,dim3) >= 0).and.(abs(dim1)+abs(dim2)+abs(dim3) == 1).and.(abs(I_m(dim,1)*dim1 + I_m(dim,2)*dim2 +I_m(dim,3)*dim3) == 0)) then  ! поперечные грани    
                    call mpi_sendrecv(  buffer_ptr  ,1  ,this%flow_send_subarray_type_transverse(dim1,dim2,dim3)    ,this%processor_neigbours(dim1,dim2,dim3)  , 1  &
                                    ,   buffer_ptr  ,1  ,this%flow_recv_subarray_type_transverse(dim1,dim2,dim3)    ,this%processor_neigbours(dim1,dim2,dim3)  , 1  &
                                    ,   MPI_COMM_WORLD,MPI_STATUS_IGNORE,error)
                end if 
            end do
            end do
            end do

            scal_ptr%cells(dim,:,:,:) = this%flow_buffer
        end do

        call MPI_BARRIER(MPI_COMM_WORLD,error)

#endif
    end subroutine

    subroutine exchange_flow_vector_field(this, vect_ptr)
        class(mpi_communications)   ,target     ,intent(inout)  :: this
        type(field_vector_flow)	    ,target     ,intent(inout)  :: vect_ptr

        real(dkind) ,pointer    ,dimension(:) :: buffer_ptr

        integer                 :: dimensions, vector_projections_number
        integer ,dimension(3)   :: processor_number, processor_grid_coord

        integer                                                :: dim,dim1,dim2,dim3, dim4
        integer ,dimension(3)                                  :: neighbour_coords
        integer                                                :: error

#ifdef mpi
        integer ,dimension(-1:1,-1:1,-1:1,mpi_status_size)     :: sr_status
#endif

#ifdef mpi

        processor_grid_coord    = this%domain%get_processor_grid_coord()
        processor_number        = this%domain%get_processor_number()

        dimensions              = this%domain%get_domain_dimensions()

        vector_projections_number = vect_ptr%get_projections_number_flow()

        do dim = 1, vector_projections_number

            do dim4 = 1, dimensions

                this%flow_buffer = vect_ptr%pr(dim)%cells(dim4,:,:,:)

                nullify(buffer_ptr)
                buffer_ptr(0:size(vect_ptr%pr(dim)%cells,2)*size(vect_ptr%pr(dim)%cells,3)*size(vect_ptr%pr(dim)%cells,4)-1) => this%flow_buffer

                do dim1 = -1,1
                do dim2 = -1,1
                do dim3 = -1,1
                    if  ((this%processor_neigbours(dim1,dim2,dim3) >= 0).and.(abs(dim1)+abs(dim2)+abs(dim3) == 1).and.(abs(I_m(dim4,1)*dim1 + I_m(dim4,2)*dim2 +I_m(dim4,3)*dim3) == 1)) then       ! сонаправленные грани
                        call mpi_sendrecv(  buffer_ptr  ,1  ,this%flow_send_subarray_type_longitudinal(dim1,dim2,dim3)    ,this%processor_neigbours(dim1,dim2,dim3)  , 1  &
                                        ,   buffer_ptr  ,1  ,this%flow_recv_subarray_type_longitudinal(dim1,dim2,dim3)    ,this%processor_neigbours(dim1,dim2,dim3)  , 1  &
                                        ,   MPI_COMM_WORLD,MPI_STATUS_IGNORE,error)
                    end if 
                end do
                end do
                end do

                do dim1 = -1,1
                do dim2 = -1,1
                do dim3 = -1,1
                    if  ((this%processor_neigbours(dim1,dim2,dim3) >= 0).and.(abs(dim1)+abs(dim2)+abs(dim3) == 1).and.(abs(I_m(dim4,1)*dim1 + I_m(dim4,2)*dim2 +I_m(dim4,3)*dim3) == 0)) then       ! поперечные грани  
                        call mpi_sendrecv(  buffer_ptr  ,1  ,this%flow_send_subarray_type_transverse(dim1,dim2,dim3)    ,this%processor_neigbours(dim1,dim2,dim3)  , 1  &
                                        ,   buffer_ptr  ,1  ,this%flow_recv_subarray_type_transverse(dim1,dim2,dim3)    ,this%processor_neigbours(dim1,dim2,dim3)  , 1  &
                                        ,   MPI_COMM_WORLD,MPI_STATUS_IGNORE,error)
                    end if 
                end do
                end do
                end do

                vect_ptr%pr(dim)%cells(dim4,:,:,:) = this%flow_buffer

            end do
        end do

        call MPI_BARRIER(MPI_COMM_WORLD,error)
#endif
    end subroutine

end module mpi_communications_class
