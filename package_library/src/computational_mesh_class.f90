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
		real(dkind)	,dimension(3)		:: cell_edges_length
		real(dkind)						:: cell_volume
		real(dkind)	,dimension(3)		:: cell_surface_area

		real(rkind)	,dimension(:,:,:,:)	,allocatable	,public	:: mesh
	contains

		procedure	,private	:: generate_uniform_mesh

		! Getters
		procedure	:: get_cell_edges_length
		procedure	:: get_cell_volume
		procedure	:: get_cell_surface_area
	end type computational_mesh

	interface computational_mesh_c
		module procedure constructor
	end interface

contains

	type(computational_mesh)	function constructor(domain)
		type(computational_domain)	,intent(in)	:: domain

		integer					:: dimensions
		integer	,dimension(3,2)	:: allocation_bounds
		integer	,dimension(3)	:: cells_number, global_cells_number
		real(dkind)	,dimension(:,:)	,allocatable	:: domain_lengths

		integer	:: dim

		dimensions				= domain%get_domain_dimensions()
		allocation_bounds		= domain%get_local_utter_cells_bounds()
		global_cells_number		= domain%get_global_cells_number()
		cells_number			= domain%get_local_cells_number()

		global_cells_number(1:dimensions)	= global_cells_number(1:dimensions) - 2

		allocate(domain_lengths,source = domain%get_domain_lengths())

		allocate(constructor%mesh(	dimensions			, &
									allocation_bounds(1,1):allocation_bounds(1,2)	, &
									allocation_bounds(2,1):allocation_bounds(2,2)	, &
									allocation_bounds(3,1):allocation_bounds(3,2)))

		constructor%cell_volume	= 1.0_dkind
		do dim = 1,dimensions
			constructor%cell_edges_length(dim)	= (domain_lengths(dim,2) - domain_lengths(dim,1)) / global_cells_number(dim)
			constructor%cell_volume 			= constructor%cell_volume*constructor%cell_edges_length(dim)
		end do

		do dim = 1,dimensions
			constructor%cell_surface_area(dim)	= constructor%cell_volume / constructor%cell_edges_length(dim)
		end do

		call constructor%generate_uniform_mesh(domain)
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

				this%mesh(dim,i,j,k) = domain_lengths(dim,1) + (decomposition_offset(dim) + 0.5_dkind + dimless_length)*this%cell_edges_length(dim)
			end do
			end do
			end do
		end do

	end subroutine


! ************** Getters ***************

	pure function get_cell_volume(this)
		class(computational_mesh)	,intent(in)	:: this
		real(dkind)								:: get_cell_volume

		get_cell_volume = this%cell_volume
	end function	

	pure function get_cell_edges_length(this)
		class(computational_mesh)	,intent(in)		:: this
		real(dkind)	,dimension(3)					:: get_cell_edges_length

		get_cell_edges_length = this%cell_edges_length
	end function	

	pure function get_cell_surface_area(this)
		class(computational_mesh)	,intent(in)		:: this
		real(dkind)	,dimension(3)					:: get_cell_surface_area

		get_cell_surface_area = this%cell_surface_area
	end function		

end module
