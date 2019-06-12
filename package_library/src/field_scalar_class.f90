module field_scalar_class

	use global_data
	use kind_parameters
	use computational_domain_class
	use boundary_conditions_class

	implicit none

	private
	public field_scalar, field_scalar_flow, field_scalar_cons, field_scalar_flow_c, field_scalar_cons_c

	type :: field_scalar
		character(len = 80)			:: name_long
		character(len = 40)			:: name_short
	end type field_scalar

	type, extends(field_scalar)	:: field_scalar_cons
		real(dkind)		,dimension(:,:,:)	,allocatable	:: cells
		integer												:: dimensions
	contains

		! Data analysis 
		procedure	:: get_field_max
		procedure	:: get_field_min
		procedure	:: get_field_mean
		procedure	:: get_field_max_grad
		procedure	:: get_field_min_grad
		procedure	:: get_field_sum	
	end type
	
	type, extends(field_scalar)	:: field_scalar_flow
		real(dkind)		,dimension(:,:,:,:)	,allocatable	:: cells
	end type

	interface field_scalar_flow_c
		module procedure constructor_flow
	end interface	
	
	interface field_scalar_cons_c
		module procedure constructor_cons
	end interface
	
	interface field_scalar_c
		module procedure constructor
	end interface

contains

	type(field_scalar)	function constructor(field_name_long,field_name_short)
		character(len=*)			,intent(in)	:: field_name_long
		character(len=*)			,intent(in)	:: field_name_short

		constructor%name_short	= field_name_short
		constructor%name_long	= field_name_long
	end function

	type(field_scalar_cons)	function constructor_cons(field_name_long, field_name_short, domain)
		character(len=*)			,intent(in)	:: field_name_long
		character(len=*)			,intent(in)	:: field_name_short
		type(computational_domain)	,intent(in)	:: domain	
		
		integer					:: dimensions
		integer	,dimension(3,2)	:: allocation_bounds

		constructor_cons%field_scalar = field_scalar_c(field_name_long, field_name_short)

		dimensions			= domain%get_domain_dimensions()
		allocation_bounds	= domain%get_local_utter_cells_bounds()

		allocate(constructor_cons%cells(	allocation_bounds(1,1):allocation_bounds(1,2)	, &
											allocation_bounds(2,1):allocation_bounds(2,2)	, &
											allocation_bounds(3,1):allocation_bounds(3,2)))

		constructor_cons%cells		= 0.0_dkind
		constructor_cons%dimensions	= dimensions	
	end function

	type(field_scalar_flow)	function constructor_flow(field_name_long, field_name_short, domain)
		character(len=*)			,intent(in)	:: field_name_long
		character(len=*)			,intent(in)	:: field_name_short
		type(computational_domain)	,intent(in)	:: domain	
		
		integer					:: dimensions
		integer	,dimension(3,2)	:: allocation_bounds

		constructor_flow%field_scalar = field_scalar_c(field_name_long, field_name_short)

		dimensions			= domain%get_domain_dimensions()
		allocation_bounds	= domain%get_local_utter_faces_bounds()

		allocate(constructor_flow%cells(	dimensions			, &
											allocation_bounds(1,1):allocation_bounds(1,2)	, &
											allocation_bounds(2,1):allocation_bounds(2,2)	, &
											allocation_bounds(3,1):allocation_bounds(3,2)))

		constructor_flow%cells		= 0.0_dkind

	end function	
	
	subroutine get_field_max(this,bc_ptr,bound,max_value,max_indexes)
		class(field_scalar_cons)			,intent(in)					:: this
		type(boundary_conditions_pointer)	,intent(in)					:: bc_ptr
		integer	,dimension(3,2)				,intent(in)		,optional	:: bound
		real(dkind)							,intent(out)				:: max_value
		integer	,dimension(3)				,intent(out)	,optional	:: max_indexes
		
		integer	,dimension(3,2)	:: bound_sl

		integer		:: i,j,k ,dim

		if (present(bound)) then
			do dim = 1, 3
				bound_sl(dim,1) = max(bound(dim,1),lbound(this%cells,dim))
				bound_sl(dim,2) = min(bound(dim,2),ubound(this%cells,dim))
			end do
		else
			do dim = 1, 3
				bound_sl(dim,1) = lbound(this%cells,dim)
				bound_sl(dim,2) = ubound(this%cells,dim)
			end do
		end if

		max_value = this%cells(bound_sl(1,1),bound_sl(2,1),bound_sl(3,1))
		
		!$omp parallel default(none)  private(i,j,k,max_indexes) , &
		!$omp& shared(this,bc_ptr,bound_sl,max_value)
		!$omp do collapse(3) schedule(static) reduction(max:max_value)

		do k = bound_sl(3,1),bound_sl(3,2)
		do j = bound_sl(2,1),bound_sl(2,2)
		do i = bound_sl(1,1),bound_sl(1,2)
			if(bc_ptr%bc_ptr%bc_markers(i,j,k) == 0) then
				if (this%cells(i,j,k) >= max_value) then
					max_value = this%cells(i,j,k)
					if(present(max_indexes)) max_indexes = (/i,j,k/)
				end if
			end if
		end do
		end do
		end do
		!$omp end do nowait
		!$omp end parallel
		
	end subroutine

	subroutine get_field_min(this,bc_ptr,bound,min_value,min_indexes)
		class(field_scalar_cons)			,intent(in)					:: this
		type(boundary_conditions_pointer)	,intent(in)					:: bc_ptr
		integer	,dimension(3,2)				,intent(in)		,optional	:: bound
		real(dkind)							,intent(out)				:: min_value
		integer	,dimension(3)				,intent(out)	,optional	:: min_indexes
		integer	,dimension(3,2)	:: bound_sl

		integer		:: i,j,k, dim

		if (present(bound)) then
			do dim = 1, 3
				bound_sl(dim,1) = max(bound(dim,1),lbound(this%cells,dim))
				bound_sl(dim,2) = min(bound(dim,2),ubound(this%cells,dim))
			end do
		else
			do dim = 1, 3
				bound_sl(dim,1) = lbound(this%cells,dim)
				bound_sl(dim,2) = ubound(this%cells,dim)
			end do
		end if

		min_value = this%cells(bound_sl(1,1),bound_sl(2,1),bound_sl(3,1))

		do k = bound_sl(3,1),bound_sl(3,2)
		do j = bound_sl(2,1),bound_sl(2,2)
		do i = bound_sl(1,1),bound_sl(1,2)
			if(bc_ptr%bc_ptr%bc_markers(i,j,k) == 0) then
				if (this%cells(i,j,k) <= min_value) then
					min_value = this%cells(i,j,k)
					if(present(min_indexes))	min_indexes = (/i,j,k/)
				end if
			end if
		end do
		end do
		end do
	end subroutine

	subroutine get_field_mean(this,bc_ptr,bound,field_mean_value)
		class(field_scalar_cons)			,intent(in)			:: this
		type(boundary_conditions_pointer)	,intent(in)			:: bc_ptr
		integer	,dimension(3,2)				,intent(in)			:: bound
		real(dkind)							,intent(out)		:: field_mean_value
		
		real(dkind)				:: mean_value
		integer	,dimension(3,2)	:: bound_sl

		integer		:: i,j,k, dim, cells_number

		do dim = 1, 3
			bound_sl(dim,1) = max(bound(dim,1),lbound(this%cells,dim))
			bound_sl(dim,2) = min(bound(dim,2),ubound(this%cells,dim))
		end do

		!$omp parallel default(none)  private(i,j,k) , &
		!$omp& shared(this,bc_ptr,bound_sl,mean_value,cells_number)
			mean_value = 0.0_dkind
			cells_number = 0	
		!$omp do collapse(3) schedule(static) reduction(+:mean_value,cells_number)

		do k = bound_sl(3,1),bound_sl(3,2)
		do j = bound_sl(2,1),bound_sl(2,2)
		do i = bound_sl(1,1),bound_sl(1,2)
			if(bc_ptr%bc_ptr%bc_markers(i,j,k) == 0) then
				mean_value = mean_value + this%cells(i,j,k)
				cells_number = cells_number + 1
			end if
		end do
		end do
		end do
		!$omp end do nowait
		!$omp end parallel
		
		field_mean_value = mean_value / cells_number
		
	end subroutine

	subroutine get_field_max_grad(this,bc_ptr,bound,max_grad_value,grad_proj,max_grad_indexes)
		class(field_scalar_cons)			,intent(in)		:: this
		type(boundary_conditions_pointer)	,intent(in)		:: bc_ptr
		integer	,dimension(3,2)				,intent(in)		:: bound
		real(dkind)							,intent(out)	:: max_grad_value
		integer	,dimension(3)				,intent(out)	:: max_grad_indexes		
		integer								,intent(in)		,optional	:: grad_proj

		integer	,dimension(3,2)	:: bound_sl

		real(dkind)	:: grad, field_u, field_l, max_grad_interm
		integer		:: i,j,k, dim

		max_grad_value	= 0.0_dkind
		max_grad_interm = 0.0_dkind
		grad			= 0.0_dkind

		bound_sl(this%dimensions:,:) = 1
		
		do dim = 1, this%dimensions
			bound_sl(dim,1) = max(bound(dim,1),lbound(this%cells,dim)+1) 
			bound_sl(dim,2) = min(bound(dim,2),ubound(this%cells,dim)-1) 
		end do

		do dim = 1,this%dimensions
			if (bc_ptr%bc_ptr%bc_markers(bound_sl(1,1)+I_m(dim,1),bound_sl(2,1)+I_m(dim,2),bound_sl(3,1)+I_m(dim,3)) == 0) then
				field_u = this%cells(bound_sl(1,1)+I_m(dim,1),bound_sl(2,1)+I_m(dim,2),bound_sl(3,1)+I_m(dim,3))
			else
				field_u = this%cells(bound_sl(1,1),bound_sl(2,1),bound_sl(3,1))
			end if
			if (bc_ptr%bc_ptr%bc_markers(bound_sl(1,1)-I_m(dim,1),bound_sl(2,1)-I_m(dim,2),bound_sl(3,1)-I_m(dim,3)) == 0) then
				field_l = this%cells(bound_sl(1,1)-I_m(dim,1),bound_sl(2,1)-I_m(dim,2),bound_sl(3,1)-I_m(dim,3))
			else
				field_l = this%cells(bound_sl(1,1),bound_sl(2,1),bound_sl(3,1))
			end if
			grad	= field_u - field_l
			!grad	= grad * grad
   !
			!if (present(grad_proj)) then
			!	if(dim == grad_proj) max_grad_value = max_grad_value + grad
			!else
			!	max_grad_value = max_grad_value + grad
			!end if
			
		end do

		max_grad_value 		= grad !sqrt(max_grad_value)
		max_grad_indexes	= (/bound_sl(1,1),bound_sl(2,1),bound_sl(3,1)/)

		do k = bound_sl(3,1),bound_sl(3,2)
		do j = bound_sl(2,1),bound_sl(2,2)
		do i = bound_sl(1,1),bound_sl(1,2)
			do dim = 1,this%dimensions
				if (bc_ptr%bc_ptr%bc_markers(i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3)) == 0) then
					field_u = this%cells(i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3))
				else
					field_u = this%cells(i,j,k)
				end if
				if (bc_ptr%bc_ptr%bc_markers(i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3)) == 0) then
					field_l = this%cells(i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3))
				else
					field_l = this%cells(i,j,k)
				end if

				grad	= field_u - field_l
				!grad = grad * grad

				if (present(grad_proj)) then
					if(dim == grad_proj) max_grad_interm = grad
				else
					max_grad_interm = max_grad_interm + grad*grad
				end if
			end do
			
			if (max_grad_interm > max_grad_value) then
				max_grad_value		= max_grad_interm
				max_grad_indexes	= (/i,j,k/)
			end if
			
		end do
		end do
		end do
	end subroutine

	subroutine get_field_min_grad(this,bc_ptr,bound,min_grad_value,grad_proj,min_grad_indexes)
		class(field_scalar_cons)			,intent(in)					:: this
		type(boundary_conditions_pointer)	,intent(in)					:: bc_ptr
		integer	,dimension(3,2)				,intent(in)					:: bound
		real(dkind)							,intent(out)				:: min_grad_value
		integer	,dimension(3)				,intent(out)				:: min_grad_indexes		
		integer								,intent(in)		,optional	:: grad_proj		

		integer	,dimension(3,2)	:: bound_sl

		real(dkind)	:: grad, field_u, field_l, min_grad_interm
		integer		:: i,j,k, dim

		min_grad_value		= 0.0_dkind
		min_grad_interm		= 0.0_dkind
		grad				= 0.0_dkind

		bound_sl(this%dimensions:,:) = 1
		
		do dim = 1, this%dimensions
			bound_sl(dim,1) = max(bound(dim,1),lbound(this%cells,dim)+1) 
			bound_sl(dim,2) = min(bound(dim,2),ubound(this%cells,dim)-1) 
		end do

		do dim = 1,this%dimensions
			if (bc_ptr%bc_ptr%bc_markers(bound_sl(1,1)+I_m(dim,1),bound_sl(2,1)+I_m(dim,2),bound_sl(3,1)+I_m(dim,3)) == 0) then
				field_u = this%cells(bound_sl(1,1)+I_m(dim,1),bound_sl(2,1)+I_m(dim,2),bound_sl(3,1)+I_m(dim,3))
			else
				field_u = this%cells(bound_sl(1,1),bound_sl(2,1),bound_sl(3,1))
			end if
			if (bc_ptr%bc_ptr%bc_markers(bound_sl(1,1)-I_m(dim,1),bound_sl(2,1)-I_m(dim,2),bound_sl(3,1)-I_m(dim,3)) == 0) then
				field_l = this%cells(bound_sl(1,1)-I_m(dim,1),bound_sl(2,1)-I_m(dim,2),bound_sl(3,1)-I_m(dim,3))
			else
				field_l = this%cells(bound_sl(1,1),bound_sl(2,1),bound_sl(3,1))
			end if
			grad	= field_u - field_l
!			grad	= grad * grad

			if (present(grad_proj)) then
				if(dim == grad_proj) min_grad_value = grad
			else
				min_grad_value = min_grad_value + grad*grad
			end if

		end do

		min_grad_value = grad!sqrt(min_grad_value)
		min_grad_indexes = (/bound_sl(1,1),bound_sl(2,1),bound_sl(3,1)/)

		do k = bound_sl(3,1),bound_sl(3,2)
		do j = bound_sl(2,1),bound_sl(2,2)
		do i = bound_sl(1,1),bound_sl(1,2)
			do dim = 1,this%dimensions
				if (bc_ptr%bc_ptr%bc_markers(i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3)) == 0) then
					field_u = this%cells(i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3))
				else
					field_u = this%cells(i,j,k)
				end if
				if (bc_ptr%bc_ptr%bc_markers(i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3)) == 0) then
					field_l = this%cells(i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3))
				else
					field_l = this%cells(i,j,k)
				end if

				grad	= field_u - field_l
				!grad = grad * grad

				if (present(grad_proj)) then
					if(dim == grad_proj) min_grad_interm = grad
				else
					min_grad_interm = min_grad_interm + grad*grad
				end if				
			end do

			if (min_grad_interm < min_grad_value) then
				min_grad_value		= min_grad_interm
				min_grad_indexes	= (/i,j,k/)
			end if
		
		end do
		end do
		end do
	end subroutine

	subroutine get_field_sum(this,bc_ptr,bound,sum_value)
		class(field_scalar_cons)			,intent(in)		:: this
		type(boundary_conditions_pointer)	,intent(in)		:: bc_ptr
		integer	,dimension(3,2)				,intent(in)		:: bound
		real(dkind)							,intent(out)	:: sum_value
		integer	,dimension(3,2)	:: bound_sl

		integer		:: i,j,k, dim

		sum_value = 0.0_dkind

		bound_sl(this%dimensions:,:) = 1
		
		do dim = 1, this%dimensions
			bound_sl(dim,1) = max(bound(dim,1),lbound(this%cells,dim)+1) 
			bound_sl(dim,2) = min(bound(dim,2),ubound(this%cells,dim)-1) 
		end do

		do k = bound_sl(3,1),bound_sl(3,2)
		do j = bound_sl(2,1),bound_sl(2,2)
		do i = bound_sl(1,1),bound_sl(1,2)
			if(bc_ptr%bc_ptr%bc_markers(i,j,k) == 0) then
				sum_value = sum_value + this%cells(i,j,k)
			end if
		end do
		end do
		end do
	end subroutine

end module
