module field_tensor_class

	use kind_parameters
	use field_scalar_class
	use chemical_properties_class
	use computational_domain_class

	implicit none

	private
	public	:: field_tensor, field_tensor_c

	type field_tensor
		character(len = 80)										:: name_long
		character(len = 40)										:: name_short
		type(field_scalar_cons)	,dimension(:,:)	,allocatable	:: pr
	contains
		! Getters
		procedure	:: get_projections_number
	end type field_tensor

	interface field_tensor_c
		module procedure constructor
	end interface

contains

	type(field_tensor)	function constructor(field_name_long,field_name_short,domain)
		character(len=*)							,intent(in)		:: field_name_long
		character(len=*)							,intent(in)		:: field_name_short
		type(computational_domain)				,intent(in)			:: domain

		integer								:: dimensions

		character(len=5)	,dimension(:)	,allocatable	:: projection_names
		integer	:: dim1, dim2

		constructor%name_short	= field_name_short
		constructor%name_long	= field_name_long

		dimensions			= domain%get_domain_dimensions()

		allocate(projection_names, source = domain%get_axis_names())
		allocate(constructor%pr(dimensions,dimensions))
		do dim1 = 1,dimensions
		do dim2 = 1,dimensions
			constructor%pr(dim1,dim2) = field_scalar_cons_c(	trim(field_name_long)	//'('//trim(projection_names(dim1))//trim(projection_names(dim2))//')',	&
																trim(field_name_short)	//'('//trim(projection_names(dim1))//trim(projection_names(dim2))//')',domain)
		end do
		end do

	end function

	pure function get_projections_number(this)
		class(field_tensor)	,intent(in)	:: this
		integer	,dimension(2)					:: get_projections_number
	
		get_projections_number(1) = size(this%pr,1)
		get_projections_number(2) = size(this%pr,2)
	end function

end module
