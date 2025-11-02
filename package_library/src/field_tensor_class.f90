module field_tensor_class

	use kind_parameters
	use field_scalar_class
	use chemical_properties_class
	use computational_domain_class

	implicit none

	private
	public	:: field_tensor, field_tensor_flow, field_tensor_cons, field_tensor_flow_c, field_tensor_cons_c 

	type field_tensor
		character(len = 80)										:: name_long
		character(len = 40)										:: name_short
    end type field_tensor
   
    type, extends(field_tensor) :: field_tensor_cons
        type(field_scalar_cons) ,dimension(:,:)   ,allocatable  :: pr
    contains
		! Getters
		procedure	:: get_projections_number
    end type
        
    type, extends(field_tensor) :: field_tensor_flow
        type(field_scalar_flow) ,dimension(:,:)   ,allocatable  :: pr
    contains
		! Getters
		procedure	:: get_projections_number_flow
    end type
    
	interface field_tensor_c
		module procedure constructor
	end interface
	
	interface field_tensor_cons_c
		module procedure constructor_cons
	end interface
	
	interface field_tensor_flow_c
		module procedure constructor_flow
	end interface
contains

	type(field_tensor)	function constructor(field_name_long,field_name_short)
		character(len=*)							,intent(in)		:: field_name_long
		character(len=*)							,intent(in)		:: field_name_short

		constructor%name_short	= field_name_short
		constructor%name_long	= field_name_long
    end function

	type(field_tensor_cons)	function constructor_cons(field_name_long,field_name_short,domain)
		character(len=*)					,intent(in)				:: field_name_long
		character(len=*)					,intent(in)				:: field_name_short
		type(computational_domain)			,intent(in)				:: domain
		
		integer								:: dimensions
		character(len=5)	,dimension(3)	:: projection_names

		integer	:: dim1, dim2
		
		constructor_cons%field_tensor = field_tensor_c(field_name_long,field_name_short)
		
		dimensions		= domain%get_domain_dimensions()

		projection_names = domain%get_axis_names()
        allocate(constructor_cons%pr(dimensions, dimensions))
		do dim1 = 1,dimensions
		do dim2 = 1,dimensions
			constructor_cons%pr(dim1,dim2) = field_scalar_cons_c(	trim(field_name_long)	//'('//trim(projection_names(dim1))//trim(projection_names(dim2))//')',	&
																    trim(field_name_short)	//'('//trim(projection_names(dim1))//trim(projection_names(dim2))//')',domain)
		end do
		end do
	end function	
	
	type(field_tensor_flow)	function constructor_flow(field_name_long,field_name_short,domain)
		character(len=*)					,intent(in)				:: field_name_long
		character(len=*)					,intent(in)				:: field_name_short
		type(computational_domain)			,intent(in)				:: domain
		
		integer								:: dimensions
		character(len=5)	,dimension(3)	:: projection_names

		integer	:: dim1, dim2
		
		constructor_flow%field_tensor = field_tensor_c(field_name_long,field_name_short)
		
		dimensions		= domain%get_domain_dimensions()

		allocate(constructor_flow%pr(dimensions, dimensions))
		projection_names = domain%get_axis_names()
		do dim1 = 1,dimensions
		do dim2 = 1,dimensions
			constructor_flow%pr(dim1,dim2) = field_scalar_flow_c(	trim(field_name_long)	//'('//trim(projection_names(dim1))//trim(projection_names(dim2))//')',	&
															        trim(field_name_short)	//'('//trim(projection_names(dim1))//trim(projection_names(dim2))//')',domain)
		end do
		end do
	end function	

	pure function get_projections_number(this)
		class(field_tensor_cons)	,intent(in)	:: this
		integer	,dimension(2)					:: get_projections_number
	
		get_projections_number(1) = size(this%pr,1)
		get_projections_number(2) = size(this%pr,2)
    end function

	pure function get_projections_number_flow(this)
		class(field_tensor_flow)	,intent(in)	:: this
		integer	,dimension(2)					:: get_projections_number_flow
        
		get_projections_number_flow(1) = size(this%pr,1)
		get_projections_number_flow(2) = size(this%pr,2)
	end function
end module
