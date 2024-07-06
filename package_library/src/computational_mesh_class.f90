module computational_mesh_class

	use kind_parameters
	use global_data
	use computational_domain_class

	implicit none

	private
	public :: computational_mesh, computational_mesh_pointer, computational_mesh_c

	type computational_mesh_pointer
		type(computational_mesh)	,pointer	:: mesh_ptr
	end type computational_mesh_pointer

	type computational_mesh
		private
		real(dkind)	,dimension(:,:,:,:)	,allocatable		:: cell_edges_length
		real(dkind) ,dimension(:,:,:)	,allocatable        :: cell_volume
		real(dkind)	,dimension(:,:,:,:)	,allocatable		:: cell_surface_area

		real(rkind)	,dimension(:,:,:,:)	,allocatable	,public	:: mesh
	contains

		procedure	,private	:: generate_uniform_mesh
		procedure	,private	:: generate_mesh

		! Getters
		procedure	:: get_cell_edges_length_loc
        procedure	:: get_cell_edges_length
		procedure	:: get_cell_volume_loc
        procedure	:: get_cell_volume
		procedure	:: get_cell_surface_area_loc
        procedure	:: get_cell_surface_area
	end type computational_mesh

	interface computational_mesh_c
		module procedure constructor
	end interface

contains

	type(computational_mesh)	function constructor(domain)
		type(computational_domain)	,intent(in)	:: domain

		integer											:: dimensions
		integer	,dimension(3,2)							:: allocation_bounds
		integer	,dimension(3)							:: cells_number, global_cells_number
        real(dkind) ,dimension(:,:) ,allocatable		:: reference_coordinates
        real(dkind) ,dimension(:,:) ,allocatable		:: reference_cell_sizes
        integer		,dimension(:)	,allocatable		:: reference_coordinates_number
        integer		,dimension(:,:) ,allocatable		:: local_mesh_cells_number
        real(dkind)	,dimension(:)	,allocatable		:: cell_sizes_koefs
		real(dkind)	,dimension(:,:)	,allocatable		:: domain_lengths
        real(dkind)										:: domain_side_length
        character(len=20) ,dimension(:)	,allocatable	:: mesh_types
        
        real(dkind)	:: cell_size_koef

		integer		:: i, j, k, dim

		dimensions				= domain%get_domain_dimensions()
		allocation_bounds		= domain%get_local_utter_cells_bounds()
		global_cells_number		= domain%get_global_cells_number()
		cells_number			= domain%get_local_cells_number()

		global_cells_number(1:dimensions)	= global_cells_number(1:dimensions) - 2

		allocate(domain_lengths,					source = domain%get_domain_lengths())
        allocate(mesh_types,						source = domain%get_mesh_types())
        allocate(reference_coordinates,				source = domain%get_reference_coordinates())
        allocate(reference_cell_sizes,				source = domain%get_reference_cell_sizes())
        allocate(reference_coordinates_number,		source = domain%get_reference_coordinates_number())
        allocate(local_mesh_cells_number,			source = domain%get_local_mesh_cells_number())
        allocate(cell_sizes_koefs(dimensions))
        
        allocate(constructor%mesh(	dimensions			, &
									allocation_bounds(1,1):allocation_bounds(1,2)	, &
									allocation_bounds(2,1):allocation_bounds(2,2)	, &
									allocation_bounds(3,1):allocation_bounds(3,2)))
        
        allocate(constructor%cell_edges_length(		dimensions			, &
													allocation_bounds(1,1):allocation_bounds(1,2)	, &
													allocation_bounds(2,1):allocation_bounds(2,2)	, &
													allocation_bounds(3,1):allocation_bounds(3,2)))
        
        allocate(constructor%cell_volume(			allocation_bounds(1,1):allocation_bounds(1,2)	, &
													allocation_bounds(2,1):allocation_bounds(2,2)	, &
													allocation_bounds(3,1):allocation_bounds(3,2)))
        
        allocate(constructor%cell_surface_area(		dimensions			, &
													allocation_bounds(1,1):allocation_bounds(1,2)	, &
													allocation_bounds(2,1):allocation_bounds(2,2)	, &
													allocation_bounds(3,1):allocation_bounds(3,2)))

		do dim = 1,dimensions
            select case (mesh_types(dim))
            case('uniform')
                domain_side_length = domain_lengths(dim,2) - domain_lengths(dim,1)
				constructor%cell_edges_length(dim, :, :, :)	= domain_side_length / global_cells_number(dim)
                
            case('linear')
                i = allocation_bounds(dim,1)
                do j = 1,reference_coordinates_number(dim)-1
                    cell_size_koef = (reference_cell_sizes(dim, j+1) - reference_cell_sizes(dim, j)) / (local_mesh_cells_number(dim, j) - 1)
                    if (dim == 1) constructor%cell_edges_length(dim,i,:,:)   = reference_cell_sizes(dim, j)
                    if (dim == 2) constructor%cell_edges_length(dim,:,i,:)   = reference_cell_sizes(dim, j)
                    if (dim == 3) constructor%cell_edges_length(dim,:,:,i)   = reference_cell_sizes(dim, j)
                    i = i + 1
                    do k = 1,local_mesh_cells_number(dim, j)
                        if (dim == 1) constructor%cell_edges_length(dim,i,:,:)	= reference_cell_sizes(dim, j) + cell_size_koef * (k - 1)
						if (dim == 2) constructor%cell_edges_length(dim,:,i,:)	= reference_cell_sizes(dim, j) + cell_size_koef * (k - 1)
						if (dim == 3) constructor%cell_edges_length(dim,:,:,i)	= reference_cell_sizes(dim, j) + cell_size_koef * (k - 1)
                        i = i + 1
                    end do
                end do
                if (dim == 1) constructor%cell_edges_length(dim,allocation_bounds(dim,2),:,:)   = reference_cell_sizes(dim, reference_coordinates_number(dim))
                if (dim == 2) constructor%cell_edges_length(dim,:,allocation_bounds(dim,2),:)   = reference_cell_sizes(dim, reference_coordinates_number(dim))
                if (dim == 3) constructor%cell_edges_length(dim,:,:,allocation_bounds(dim,2))   = reference_cell_sizes(dim, reference_coordinates_number(dim))
                
            case('exponential')
                i = allocation_bounds(dim,1)
                do j = 1,reference_coordinates_number(dim)-1
                    domain_side_length = reference_coordinates(dim, j+1) - reference_coordinates(dim, j)
                    cell_size_koef = (domain_side_length - reference_cell_sizes(dim, j)) / (domain_side_length - reference_cell_sizes(dim, j+1))
                    if (dim == 1) constructor%cell_edges_length(dim,i,:,:)   = reference_cell_sizes(dim, j)
                    if (dim == 2) constructor%cell_edges_length(dim,:,i,:)   = reference_cell_sizes(dim, j)
                    if (dim == 3) constructor%cell_edges_length(dim,:,:,i)   = reference_cell_sizes(dim, j)
                    i = i + 1
                    do k = 1,local_mesh_cells_number(dim, j)
                        if (dim == 1) constructor%cell_edges_length(dim,i,:,:)	= reference_cell_sizes(dim, j) * cell_size_koef ** (k - 1)
						if (dim == 2) constructor%cell_edges_length(dim,:,i,:)	= reference_cell_sizes(dim, j) * cell_size_koef ** (k - 1)
						if (dim == 3) constructor%cell_edges_length(dim,:,:,i)	= reference_cell_sizes(dim, j) * cell_size_koef ** (k - 1)
                        i = i + 1
                    end do
                end do
                if (dim == 1) constructor%cell_edges_length(dim,allocation_bounds(dim,2),:,:)   = reference_cell_sizes(dim, reference_coordinates_number(dim))
                if (dim == 2) constructor%cell_edges_length(dim,:,allocation_bounds(dim,2),:)   = reference_cell_sizes(dim, reference_coordinates_number(dim))
                if (dim == 3) constructor%cell_edges_length(dim,:,:,allocation_bounds(dim,2))   = reference_cell_sizes(dim, reference_coordinates_number(dim))
                
            end select
        end do
        
        constructor%cell_volume	= 1.0_dkind
        do k = allocation_bounds(3,1),allocation_bounds(3,2)
        do j = allocation_bounds(2,1),allocation_bounds(2,2)
        do i = allocation_bounds(1,1),allocation_bounds(1,2)
            do dim = 1,dimensions
                constructor%cell_volume(i,j,k)	= constructor%cell_volume(i,j,k)*constructor%cell_edges_length(dim,i,j,k)
            end do
        end do
        end do
        end do

        do dim = 1,dimensions
            do k = allocation_bounds(3,1),allocation_bounds(3,2)
            do j = allocation_bounds(2,1),allocation_bounds(2,2)
            do i = allocation_bounds(1,1),allocation_bounds(1,2)
                constructor%cell_surface_area(dim,i,j,k)	= constructor%cell_volume(i,j,k) / constructor%cell_edges_length(dim,i,j,k)
            end do
            end do
            end do
        end do
            
        call constructor%generate_mesh(domain)
	end function

	subroutine generate_uniform_mesh(this,domain)

		! Generation of the uniform structured computational grid in computational domain "domain"

		class(computational_mesh)	,intent(inout)	:: this
		type(computational_domain)	,intent(in)		:: domain

		real(dkind)	:: dimless_length
		integer		:: i, j, k, dim

		integer										:: dimensions
		integer	,dimension(3,2)						:: loop_bounds
		integer ,dimension(3)						:: decomposition_offset
		real(dkind)	,dimension(:,:)	,allocatable	:: domain_lengths

		dimensions				= domain%get_domain_dimensions()
		decomposition_offset	= domain%get_global_offset()

		loop_bounds				= domain%get_local_utter_cells_bounds()

		allocate(domain_lengths,source = domain%get_domain_lengths())		

		do dim = 1,dimensions
			do k = loop_bounds(3,1),loop_bounds(3,2)
			do j = loop_bounds(2,1),loop_bounds(2,2)
			do i = loop_bounds(1,1),loop_bounds(1,2)
				dimless_length = real((i-1)*I_m(dim,1)+(j-1)*I_m(dim,2)+(k-1)*I_m(dim,3),dkind)

				this%mesh(dim,i,j,k) = domain_lengths(dim,1) + (decomposition_offset(dim) + 0.5_dkind + dimless_length)*this%cell_edges_length(dim,i,j,k)
			end do
			end do
			end do
		end do

    end subroutine
    
    
    subroutine generate_mesh(this,domain)

		! Generation of the structured computational grid in computational domain "domain"

		class(computational_mesh)	,intent(inout)	:: this
		type(computational_domain)	,intent(in)		:: domain

		real(dkind)	:: dimless_length
		integer		:: i, j, k, dim

		integer										:: dimensions
		integer	,dimension(3,2)						:: loop_bounds
		integer ,dimension(3)						:: decomposition_offset
		real(dkind)	,dimension(:,:)	,allocatable	:: domain_lengths

		dimensions				= domain%get_domain_dimensions()
		decomposition_offset	= domain%get_global_offset()

		loop_bounds				= domain%get_local_utter_cells_bounds()

		allocate(domain_lengths,source = domain%get_domain_lengths())
        

		do dim = 1,dimensions
			do k = loop_bounds(3,1),loop_bounds(3,2)
			do j = loop_bounds(2,1),loop_bounds(2,2)
			do i = loop_bounds(1,1),loop_bounds(1,2)
                if (dim == 1) then
                    if (i == loop_bounds(1,1)) then
						dimless_length = real((i-1)*I_m(dim,1)+(j-1)*I_m(dim,2)+(k-1)*I_m(dim,3),dkind)
						this%mesh(dim,i,j,k) = domain_lengths(dim,1) + (decomposition_offset(dim) + 0.5_dkind + dimless_length)*this%cell_edges_length(dim,i,j,k)
					else
						this%mesh(dim,i,j,k) = this%mesh(dim,i-1,j,k) + (this%cell_edges_length(dim,i-1,j,k) + this%cell_edges_length(dim,i,j,k)) / 2.0_dkind
                    end if
                elseif (dim == 2) then
					if (j == loop_bounds(2,1)) then
						dimless_length = real((i-1)*I_m(dim,1)+(j-1)*I_m(dim,2)+(k-1)*I_m(dim,3),dkind)
						this%mesh(dim,i,j,k) = domain_lengths(dim,1) + (decomposition_offset(dim) + 0.5_dkind + dimless_length)*this%cell_edges_length(dim,i,j,k)
					else
						this%mesh(dim,i,j,k) = this%mesh(dim,i,j-1,k) + (this%cell_edges_length(dim,i,j-1,k) + this%cell_edges_length(dim,i,j,k)) / 2.0_dkind
                    end if
                elseif (dim == 3) then
					if (k == loop_bounds(3,1)) then
						dimless_length = real((i-1)*I_m(dim,1)+(j-1)*I_m(dim,2)+(k-1)*I_m(dim,3),dkind)
						this%mesh(dim,i,j,k) = domain_lengths(dim,1) + (decomposition_offset(dim) + 0.5_dkind + dimless_length)*this%cell_edges_length(dim,i,j,k)
					else
						this%mesh(dim,i,j,k) = this%mesh(dim,i,j,k-1) + (this%cell_edges_length(dim,i,j,k-1) + this%cell_edges_length(dim,i,j,k)) / 2.0_dkind
                    end if
                end if
                !print *, decomposition_offset(dim), dim, i, j, k, dimless_length, this%mesh(dim,i,j,k)
			end do
			end do
			end do
        end do

	end subroutine


! ************** Getters ***************

	pure function get_cell_volume_loc(this,i,j,k)
		class(computational_mesh)	,intent(in)	:: this
        integer						,intent(in)	:: i, j, k
		real(dkind)								:: get_cell_volume_loc

		get_cell_volume_loc = this%cell_volume(i,j,k)
	end function	

	pure function get_cell_edges_length_loc(this,i,j,k)
		class(computational_mesh)	,intent(in)		:: this
        integer						,intent(in)	:: i, j, k
		real(dkind)	,dimension(3)					:: get_cell_edges_length_loc

		get_cell_edges_length_loc = this%cell_edges_length(:,i,j,k)
	end function	

	pure function get_cell_surface_area_loc(this,i,j,k)
		class(computational_mesh)	,intent(in)		:: this
        integer						,intent(in)	:: i, j, k
		real(dkind)	,dimension(3)					:: get_cell_surface_area_loc

		get_cell_surface_area_loc = this%cell_surface_area(:,i,j,k)
    end function
    
    pure function get_cell_volume(this)
		class(computational_mesh)	,intent(in)	:: this
		real(dkind)								:: get_cell_volume

		get_cell_volume = this%cell_volume(1,1,1)
	end function	

	pure function get_cell_edges_length(this)
		class(computational_mesh)	,intent(in)		:: this
		real(dkind)	,dimension(3)					:: get_cell_edges_length

		get_cell_edges_length = this%cell_edges_length(:,1,1,1)
	end function	

	pure function get_cell_surface_area(this)
		class(computational_mesh)	,intent(in)		:: this
		real(dkind)	,dimension(3)					:: get_cell_surface_area

		get_cell_surface_area = this%cell_surface_area(:,1,1,1)
	end function		

end module
