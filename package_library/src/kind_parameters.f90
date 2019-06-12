module kind_parameters
    !    Some examples of kind parameters
    !    double precision : selected_real_kind(15,307)
    !    single precision : selected_real_kind(6,37)

	integer ,parameter	:: dkind = selected_real_kind(15,307)
	integer ,parameter	:: rkind = selected_real_kind(6,37)
	integer ,parameter	:: ikind = selected_int_kind(2)
end module
