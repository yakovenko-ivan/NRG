module data_manager_class

	use kind_parameters
	use computational_domain_class
	use mpi_communications_class
	use chemical_properties_class
	use thermophysical_properties_class
	use computational_mesh_class
	use boundary_conditions_class
	use field_pointers
	use field_scalar_class
	use field_vector_class

	implicit none

	private
	public	:: data_manager, data_manager_c

	type data_manager
		private
		type(computational_domain)				,public	:: domain
		type(mpi_communications)				,public	:: mpi_communications
		type(computational_mesh_pointer)		,public	:: computational_mesh_pointer
		type(boundary_conditions_pointer)		,public	:: boundary_conditions_pointer
		type(chemical_properties_pointer)		,public	:: chemistry
		type(thermophysical_properties_pointer)	,public	:: thermophysics
		
		type(field_scalar_cons_pointer)	,dimension(100)	:: scalar_field_cons_pointers
		type(field_vector_cons_pointer)	,dimension(100)	:: vector_field_cons_pointers
		type(field_tensor_cons_pointer)	,dimension(100)	:: tensor_field_cons_pointers

		type(field_scalar_flow_pointer)	,dimension(100)	:: scalar_field_flow_pointers
		type(field_vector_flow_pointer)	,dimension(100)	:: vector_field_flow_pointers
		
		integer		:: number_of_cons_scalar_fields
		integer		:: number_of_cons_vector_fields
		integer		:: number_of_cons_tensor_fields
		
		integer		:: number_of_flow_scalar_fields
		integer		:: number_of_flow_vector_fields
		
	contains
	
		procedure	:: create_boundary_conditions
		procedure	:: create_computational_mesh
		procedure	:: create_scalar_field
		procedure	:: create_vector_field
		procedure	:: create_tensor_field

		! Getters
		procedure	:: get_cons_field_pointer_by_name
		procedure	:: get_cons_field_pointer_by_number	
		procedure	:: get_flow_field_pointer_by_name
		procedure	:: get_flow_field_pointer_by_number
		procedure	:: get_field_type_by_name
		procedure	:: get_number_of_cons_scalar_fields
		procedure	:: get_number_of_cons_vector_fields
		procedure	:: get_number_of_cons_tensor_fields
		procedure	:: get_number_of_flow_scalar_fields
		procedure	:: get_number_of_flow_vector_fields

	end type data_manager

	interface data_manager_c
		module procedure	constructor
	end interface

contains

	type(data_manager)	function constructor(domain,mpi_support,chemistry,thermophysics)
		type(computational_domain)		,intent(in)				:: domain
		type(mpi_communications)		,intent(in)				:: mpi_support
		type(chemical_properties)		,intent(in)	,target		:: chemistry
		type(thermophysical_properties)	,intent(in)	,target		:: thermophysics

		integer	:: field_number

		constructor%domain						= domain
		constructor%mpi_communications			= mpi_support
		constructor%chemistry%chem_ptr			=> chemistry
		constructor%thermophysics%thermo_ptr	=> thermophysics
		
		constructor%number_of_cons_scalar_fields = 0
		constructor%number_of_cons_vector_fields = 0
		constructor%number_of_cons_tensor_fields = 0
		
		constructor%number_of_flow_scalar_fields = 0
		constructor%number_of_flow_vector_fields = 0

	end function

	subroutine create_boundary_conditions(this,bc_instance,number_of_boundary_types, default_boundary)

		type(boundary_conditions)	,target	,intent(inout)	:: bc_instance
		class(data_manager)					,intent(inout)	:: this
		integer								,intent(in)	,optional	:: number_of_boundary_types, default_boundary

		if(present(number_of_boundary_types).and.present(default_boundary)) then
			bc_instance = boundary_conditions_c(this%domain,number_of_boundary_types,default_boundary)
		else
			bc_instance = boundary_conditions_c(this%domain)
		end if
		this%boundary_conditions_pointer%bc_ptr => bc_instance
	end subroutine

	subroutine create_computational_mesh(this,mesh_instance)

		type(computational_mesh)	,target	,intent(inout)	:: mesh_instance
		class(data_manager)					,intent(inout)	:: this

		mesh_instance = computational_mesh_c(this%domain)
		this%computational_mesh_pointer%mesh_ptr => mesh_instance
	end subroutine

	subroutine create_scalar_field(this,field_instance,field_name,field_short_name)
		class(field_scalar)	,target	,intent(inout)	:: field_instance
		class(data_manager)			,intent(inout)	:: this
		character(len=*)			,intent(in)		:: field_name
		character(len=*)			,intent(in)		:: field_short_name

		select type (field_instance)
			type is (field_scalar_cons)
				this%number_of_cons_scalar_fields = this%number_of_cons_scalar_fields + 1
				field_instance = field_scalar_cons_c(field_name, field_short_name, this%domain)
				this%scalar_field_cons_pointers(this%number_of_cons_scalar_fields)%s_ptr => field_instance
			type is (field_scalar_flow)
				this%number_of_flow_scalar_fields = this%number_of_flow_scalar_fields + 1
				field_instance = field_scalar_flow_c(field_name, field_short_name, this%domain)
				this%scalar_field_flow_pointers(this%number_of_flow_scalar_fields)%s_ptr => field_instance
		end select
	end subroutine

	subroutine create_vector_field(this,field_instance,field_name,field_short_name,field_type)
		class(field_vector)	,target	,intent(inout)	:: field_instance
		class(data_manager)			,intent(inout)	:: this
		character(len=*)			,intent(in)		:: field_name
		character(len=*)			,intent(in)		:: field_short_name
		character(len=*)			,intent(in)		:: field_type

		select type (field_instance)
			type is (field_vector_cons)	
				this%number_of_cons_vector_fields = this%number_of_cons_vector_fields + 1
				this%vector_field_cons_pointers(this%number_of_cons_vector_fields)%v_ptr => field_instance
				select case (field_type)
					case('spatial')
						field_instance = field_vector_cons_c(field_name, field_short_name, this%domain)
					case('chemical')
						field_instance = field_vector_cons_c(field_name, field_short_name, this%domain, this%chemistry)
					case default	
						print *, 'ERROR: Unknown field type'
				end select
			type is (field_vector_flow)
				this%number_of_flow_vector_fields = this%number_of_flow_vector_fields + 1
				this%vector_field_flow_pointers(this%number_of_flow_vector_fields)%v_ptr => field_instance
				select case (field_type)
					case('spatial')
						field_instance = field_vector_flow_c(field_name, field_short_name, this%domain)
					case('chemical')
						field_instance = field_vector_flow_c(field_name, field_short_name, this%domain, this%chemistry)
					case default	
						print *, 'ERROR: Unknown field type'
				end select
		end select
		
	end subroutine

	subroutine create_tensor_field(this,field_instance,field_name,field_short_name)
		type(field_tensor)	,target	,intent(inout)	:: field_instance
		class(data_manager)			,intent(inout)	:: this
		character(len=*)			,intent(in)		:: field_name
		character(len=*)			,intent(in)		:: field_short_name

		this%number_of_cons_tensor_fields = this%number_of_cons_tensor_fields + 1

		field_instance = field_tensor_c(field_name, field_short_name, this%domain)

		this%tensor_field_cons_pointers(this%number_of_cons_tensor_fields)%t_ptr => field_instance
	end subroutine

! ************** Getters ***************

	subroutine get_cons_field_pointer_by_name(this,scal_ptr,vect_ptr,tens_ptr,field_name)
		type(field_scalar_cons_pointer)	,intent(out)	:: scal_ptr
		type(field_vector_cons_pointer)	,intent(out)	:: vect_ptr
		type(field_tensor_cons_pointer)	,intent(out)	:: tens_ptr
		class(data_manager)			,intent(in)		:: this
		character(len=*)			,intent(in)		:: field_name

		integer	:: field_counter, dim, dim1, dim2

		scal_ptr%s_ptr => NULL()
		vect_ptr%v_ptr => NULL()
		tens_ptr%t_ptr => NULL()

		do field_counter = 1,this%number_of_cons_scalar_fields
			if(trim(this%scalar_field_cons_pointers(field_counter)%s_ptr%name_long) == trim(field_name)) then
				scal_ptr = this%scalar_field_cons_pointers(field_counter)
			end if
		end do

		do field_counter = 1,this%number_of_cons_vector_fields
			if(trim(this%vector_field_cons_pointers(field_counter)%v_ptr%name_long) == trim(field_name)) then
				vect_ptr = this%vector_field_cons_pointers(field_counter)
			else
				do dim = 1,size(this%vector_field_cons_pointers(field_counter)%v_ptr%pr)
					if(trim(this%vector_field_cons_pointers(field_counter)%v_ptr%pr(dim)%name_long) == trim(field_name)) then
						scal_ptr%s_ptr => this%vector_field_cons_pointers(field_counter)%v_ptr%pr(dim)
					end if
				end do
			end if
		end do

		do field_counter = 1,this%number_of_cons_tensor_fields
			if(trim(this%tensor_field_cons_pointers(field_counter)%t_ptr%name_long) == trim(field_name)) then
				tens_ptr = this%tensor_field_cons_pointers(field_counter)
			else
				do dim1 = 1,size(this%tensor_field_cons_pointers(field_counter)%t_ptr%pr,1)
				do dim2 = 1,size(this%tensor_field_cons_pointers(field_counter)%t_ptr%pr,2)
					if(trim(this%tensor_field_cons_pointers(field_counter)%t_ptr%pr(dim1,dim2)%name_long) == trim(field_name)) then
						scal_ptr%s_ptr => this%tensor_field_cons_pointers(field_counter)%t_ptr%pr(dim1,dim2)
					end if
				end do
				end do
			end if
		end do
	end subroutine

	subroutine get_cons_field_pointer_by_number(this,scal_ptr,vect_ptr,tens_ptr,field_type,field_number)
		type(field_scalar_cons_pointer)	,intent(out)	:: scal_ptr
		type(field_vector_cons_pointer)	,intent(out)	:: vect_ptr
		type(field_tensor_cons_pointer)	,intent(out)	:: tens_ptr
		class(data_manager)			,intent(in)		:: this
		character(len=*)			,intent(in)		:: field_type
		integer						,intent(in)		:: field_number

		integer	:: field_counter, dim, dim1, dim2

		scal_ptr%s_ptr => NULL()
		vect_ptr%v_ptr => NULL()
		tens_ptr%t_ptr => NULL()

		select case (field_type)
			case('scalar')
				scal_ptr = this%scalar_field_cons_pointers(field_number)
			case('vector')
				vect_ptr = this%vector_field_cons_pointers(field_number)
			case('tensor')
				tens_ptr = this%tensor_field_cons_pointers(field_number)
		end select
	end subroutine

	subroutine get_flow_field_pointer_by_name(this,scal_ptr,vect_ptr,field_name)
		type(field_scalar_flow_pointer)	,intent(out)	:: scal_ptr
		type(field_vector_flow_pointer)	,intent(out)	:: vect_ptr
		class(data_manager)			,intent(in)		:: this
		character(len=*)			,intent(in)		:: field_name

		integer	:: field_counter, dim, dim1, dim2

		scal_ptr%s_ptr => NULL()
		vect_ptr%v_ptr => NULL()

		do field_counter = 1,this%number_of_flow_scalar_fields
			if(trim(this%scalar_field_flow_pointers(field_counter)%s_ptr%name_long) == trim(field_name)) then
				scal_ptr = this%scalar_field_flow_pointers(field_counter)
			end if
		end do

		do field_counter = 1,this%number_of_flow_vector_fields
			if(trim(this%vector_field_flow_pointers(field_counter)%v_ptr%name_long) == trim(field_name)) then
				vect_ptr = this%vector_field_flow_pointers(field_counter)
			else
				do dim = 1,size(this%vector_field_flow_pointers(field_counter)%v_ptr%pr)
					if(trim(this%vector_field_flow_pointers(field_counter)%v_ptr%pr(dim)%name_long) == trim(field_name)) then
						scal_ptr%s_ptr => this%vector_field_flow_pointers(field_counter)%v_ptr%pr(dim)
					end if
				end do
			end if
		end do

	end subroutine

	subroutine get_flow_field_pointer_by_number(this,scal_ptr,vect_ptr,field_type,field_number)
		type(field_scalar_flow_pointer)	,intent(out)	:: scal_ptr
		type(field_vector_flow_pointer)	,intent(out)	:: vect_ptr
		class(data_manager)			,intent(in)		:: this
		character(len=*)			,intent(in)		:: field_type
		integer						,intent(in)		:: field_number

		integer	:: field_counter, dim, dim1, dim2

		scal_ptr%s_ptr => NULL()
		vect_ptr%v_ptr => NULL()

		select case (field_type)
			case('scalar')
				scal_ptr = this%scalar_field_flow_pointers(field_number)
			case('vector')
				vect_ptr = this%vector_field_flow_pointers(field_number)
		end select
	end subroutine


	function get_field_type_by_name(this,field_name)
		character(len=40)							:: get_field_type_by_name
		class(data_manager)			,intent(in)		:: this
		character(len=*)			,intent(in)		:: field_name
	
		integer	:: field_counter, dim, dim1, dim2

		do field_counter = 1,this%number_of_cons_scalar_fields
			if(trim(this%scalar_field_cons_pointers(field_counter)%s_ptr%name_long) == trim(field_name)) then
				get_field_type_by_name = 'scalar'
			end if
		end do

		do field_counter = 1,this%number_of_cons_vector_fields
			if(trim(this%vector_field_cons_pointers(field_counter)%v_ptr%name_long) == trim(field_name)) then
				get_field_type_by_name = 'vector'
			else
				do dim = 1,size(this%vector_field_cons_pointers(field_counter)%v_ptr%pr)
					if(trim(this%vector_field_cons_pointers(field_counter)%v_ptr%pr(dim)%name_long) == trim(field_name)) then
						get_field_type_by_name = 'vector'
					end if
				end do
			end if
		end do

		do field_counter = 1,this%number_of_cons_tensor_fields
			if(trim(this%tensor_field_cons_pointers(field_counter)%t_ptr%name_long) == trim(field_name)) then
				get_field_type_by_name = 'tensor'
			else
				do dim1 = 1,size(this%tensor_field_cons_pointers(field_counter)%t_ptr%pr,1)
				do dim2 = 1,size(this%tensor_field_cons_pointers(field_counter)%t_ptr%pr,2)
					if(trim(this%tensor_field_cons_pointers(field_counter)%t_ptr%pr(dim1,dim2)%name_long) == trim(field_name)) then
						get_field_type_by_name = 'scalar'
					end if
				end do
				end do
			end if
		end do
	
	end function
	
	pure integer function get_number_of_cons_scalar_fields(this)
		class(data_manager)			,intent(in)		:: this

		get_number_of_cons_scalar_fields = this%number_of_cons_scalar_fields
	end function

	pure integer function get_number_of_cons_vector_fields(this)
		class(data_manager)			,intent(in)		:: this

		get_number_of_cons_vector_fields = this%number_of_cons_vector_fields
	end function
	pure integer function get_number_of_cons_tensor_fields(this)
		class(data_manager)			,intent(in)		:: this

		get_number_of_cons_tensor_fields = this%number_of_cons_tensor_fields
	end function
	pure integer function get_number_of_flow_scalar_fields(this)
		class(data_manager)			,intent(in)		:: this

		get_number_of_flow_scalar_fields = this%number_of_flow_scalar_fields
	end function
	pure integer function get_number_of_flow_vector_fields(this)
		class(data_manager)			,intent(in)		:: this

		get_number_of_flow_vector_fields = this%number_of_flow_vector_fields
	end function


end module
