module field_pointers

	use kind_parameters
	use field_scalar_class
	use field_vector_class
	use field_tensor_class

	implicit none

	type 	:: field_scalar_cons_pointer
		type(field_scalar_cons)	,pointer	:: s_ptr
	end type field_scalar_cons_pointer
	
	type 	:: field_scalar_flow_pointer
		type(field_scalar_flow)	,pointer	:: s_ptr
	end type field_scalar_flow_pointer

	type 	:: field_vector_cons_pointer
		type(field_vector_cons)	,pointer	:: v_ptr
	end type field_vector_cons_pointer

	type 	:: field_vector_flow_pointer
		type(field_vector_flow)	,pointer	:: v_ptr
	end type field_vector_flow_pointer	
	
	type 	:: field_tensor_cons_pointer
		type(field_tensor)	,pointer	:: t_ptr
	end type field_tensor_cons_pointer	
	

end module
