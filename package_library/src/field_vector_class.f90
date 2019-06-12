module field_vector_class

	use kind_parameters
	use field_scalar_class
	use chemical_properties_class
	use computational_domain_class
	use boundary_conditions_class

	implicit none

	private
	public field_vector, field_vector_flow, field_vector_cons, field_vector_flow_c, field_vector_cons_c

	type field_vector
		character(len = 80)									:: name_long
		character(len = 40)									:: name_short
	end type field_vector

	type, extends(field_vector) ::	field_vector_cons
		type(field_scalar_cons)	,dimension(:)	,allocatable	:: pr
	contains
		! Getters
		procedure	:: get_projection_number_by_name
		procedure	:: get_projections_number
		procedure	:: get_field_mean_square
	end type
	
	type, extends(field_vector) ::	field_vector_flow
		type(field_scalar_flow)	,dimension(:)	,allocatable	:: pr
	contains
		! Getters	
		procedure	:: get_projections_number_flow
	end type

	interface field_vector_c
		module procedure constructor
	end interface
	
	interface field_vector_cons_c
		module procedure constructor_cons
	end interface
	
	interface field_vector_flow_c
		module procedure constructor_flow
	end interface

contains

	type(field_vector)	function constructor(field_name_long,field_name_short)
		character(len=*)					,intent(in)				:: field_name_long
		character(len=*)					,intent(in)				:: field_name_short

		integer	:: specie, dim

		constructor%name_short	= field_name_short
		constructor%name_long	= field_name_long
	end function
	
	type(field_vector_cons)	function constructor_cons(field_name_long,field_name_short,domain,chem)
		character(len=*)					,intent(in)				:: field_name_long
		character(len=*)					,intent(in)				:: field_name_short
		type(computational_domain)			,intent(in)				:: domain
		type(chemical_properties_pointer)	,intent(in)	,optional	:: chem
		
		integer								:: dimensions
		character(len=5)	,dimension(:)	,allocatable	:: projection_names

		integer	:: specie, dim
		
		constructor_cons%field_vector = field_vector_c(field_name_long,field_name_short)
		
		dimensions		= domain%get_domain_dimensions()

		if(present(chem)) then
			allocate(constructor_cons%pr(chem%chem_ptr%species_number))
			do specie = 1,chem%chem_ptr%species_number
				constructor_cons%pr(specie) = field_scalar_cons_c(	trim(field_name_long)	//'('//trim(chem%chem_ptr%species_names(specie))//')', &
																	trim(field_name_short)	//'('//trim(chem%chem_ptr%species_names(specie))//')',domain)
			end do
		else
			allocate(constructor_cons%pr(dimensions))
			allocate(projection_names, source = domain%get_axis_names())
			do dim = 1,dimensions
				constructor_cons%pr(dim) = field_scalar_cons_c(	trim(field_name_long)	//'('//trim(projection_names(dim))//')',	&
																trim(field_name_short)	//'('//trim(projection_names(dim))//')',domain)
			end do
		end if
	end function	
	
	type(field_vector_flow)	function constructor_flow(field_name_long,field_name_short,domain,chem)
		character(len=*)					,intent(in)				:: field_name_long
		character(len=*)					,intent(in)				:: field_name_short
		type(computational_domain)			,intent(in)				:: domain
		type(chemical_properties_pointer)	,intent(in)	,optional	:: chem
		
		integer								:: dimensions
		character(len=5)	,dimension(3)	:: projection_names

		integer	:: specie, dim
		
		constructor_flow%field_vector = field_vector_c(field_name_long,field_name_short)
		
		dimensions		= domain%get_domain_dimensions()

		if(present(chem)) then
			allocate(constructor_flow%pr(chem%chem_ptr%species_number))
			do specie = 1,chem%chem_ptr%species_number
				constructor_flow%pr(specie) = field_scalar_flow_c(	trim(field_name_long)	//'('//trim(chem%chem_ptr%species_names(specie))//')', &
																	trim(field_name_short)	//'('//trim(chem%chem_ptr%species_names(specie))//')',domain)
			end do
		else
			allocate(constructor_flow%pr(dimensions))
			projection_names = domain%get_axis_names()
			do dim = 1,dimensions
				constructor_flow%pr(dim) = field_scalar_flow_c(	trim(field_name_long)	//'('//trim(projection_names(dim))//')',	&
																trim(field_name_short)	//'('//trim(projection_names(dim))//')',domain)
			end do
		end if
	end function	
	
	subroutine get_field_mean_square(this,bc_ptr,bound,field_mean_square_value)
		class(field_vector_cons)			,intent(in)		:: this
		type(boundary_conditions_pointer)	,intent(in)		:: bc_ptr
		integer	,dimension(3,2)				,intent(in)		:: bound
		real(dkind)							,intent(out)	:: field_mean_square_value
		
		real(dkind)				:: mean_square_value
		integer	,dimension(3,2)	:: bound_sl

		integer		:: i,j,k, dim, cells_number

		do dim = 1, 3
			bound_sl(dim,1) = max(bound(dim,1),lbound(this%pr(1)%cells,dim))
			bound_sl(dim,2) = min(bound(dim,2),ubound(this%pr(1)%cells,dim))
		end do

		!$omp parallel default(none)  private(i,j,k) , &
		!$omp& shared(this,bc_ptr,bound_sl,mean_square_value,cells_number)
		mean_square_value = 0.0_dkind
		cells_number = 0
		!$omp do collapse(3) schedule(static) reduction(+:mean_square_value,cells_number)
		do k = bound_sl(3,1),bound_sl(3,2)
		do j = bound_sl(2,1),bound_sl(2,2)
		do i = bound_sl(1,1),bound_sl(1,2)
			if(bc_ptr%bc_ptr%bc_markers(i,j,k) == 0) then
				cells_number = cells_number + 1
				do dim = 1, size(this%pr)
					mean_square_value = mean_square_value + this%pr(dim)%cells(i,j,k) * this%pr(dim)%cells(i,j,k)
				end do
			end if
		end do
		end do
		end do
		!$omp end do nowait
		!$omp end parallel
		
		field_mean_square_value = sqrt (mean_square_value / cells_number)
	end subroutine
	
	integer function get_projection_number_by_name(this,projection_name)
		class(field_vector_cons)	,intent(in)	:: this
		character(len=*)			,intent(in)	:: projection_name
	
		integer	:: index_left_bracket, index_right_bracket, proj
		
		
		do proj = 1, size(this%pr)
			index_left_bracket = index(this%pr(proj)%name_short, '(')
			index_right_bracket = index(this%pr(proj)%name_short, ')') 
			if ( projection_name == this%pr(proj)%name_short(index_left_bracket+1:index_right_bracket-1)) then
				get_projection_number_by_name = proj		
			end if
		end do
		
	end function

	pure integer function get_projections_number(this)
		class(field_vector_cons)	,intent(in)	:: this
	
		get_projections_number = size(this%pr)
	end function

	pure integer function get_projections_number_flow(this)
		class(field_vector_flow)	,intent(in)	:: this
	
		get_projections_number_flow = size(this%pr)
	end function

end module
