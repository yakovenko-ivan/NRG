module implicit_fourier_heat_transfer_solver_class

	use kind_parameters
	use global_data
	use data_manager_class	
	use data_io_class
	use computational_mesh_class
	use computational_domain_class	
	use boundary_conditions_class	
	use field_pointers
	
	use thermophysical_properties_class
	use chemical_properties_class

	use mpi_communications_class

	use solver_options_class
	
	implicit none

#ifdef OMP
	include "omp_lib.h"
#endif

	private
	public	:: implicit_heat_transfer_solver, implicit_heat_transfer_solver_c

	type(field_scalar_cons)	,target	:: T,  rho, kappa, Q

	type	:: implicit_heat_transfer_solver
		
		real(dkind)		:: time, time_step, initial_time_step
	
		type(field_scalar_cons_pointer)				:: T, kappa, rho, Q, T_gas, kappa_gas
		type(computational_domain)					:: domain
		type(mpi_communications)					:: mpi_support
		type(boundary_conditions_pointer)			:: boundary
		type(computational_mesh_pointer)			:: mesh
		type(thermophysical_properties_pointer)		:: thermo
		type(chemical_properties_pointer)			:: chem

		real(dkind) ,dimension(:)		,allocatable    :: thermal_c_coeff_constant
		
		real(dkind) ,dimension(:,:,:)	,allocatable    :: B, K, C, E, V, F	
		real(dkind) ,dimension(:,:,:)	,allocatable    :: alpha, beta, gamma, delta
		real(dkind) ,dimension(:,:,:)	,allocatable    :: alpha_, beta_, gamma_, delta_
		! NOTE Variables names correspond to those used Volochinskaya J. Comp. Math & Math phys 17(2) 1977 where the description of \alpha-\beta iterative algorithm was originally given. 
		
		real(dkind)	:: phi_isoth	= 0.0_dkind, psi_isoth	= 300.0_dkind	! Uniform isothermal boundary conditions
		real(dkind)	:: phi_adiab	= 1.0_dkind, psi_adiab	= 0.0_dkind		! Symmetry plane, adiabatic boundary
		
		integer		:: mean_T_unit, energy_flux_unit
	contains
		procedure				:: solve_problem
		procedure				:: get_time_step
		procedure				:: get_time
		procedure				:: heat_transfer_solver_additional_constructor
		
		procedure	,private	::	solve_heat_transfer
		procedure	,private	::	calculate_difference_scheme_coefficients
		procedure	,private	::	alpha_beta_sweep
		procedure	,private	::	check_balance
		procedure	,private	::	get_mean_temperature
	end type

	interface	implicit_heat_transfer_solver_c
		module procedure	constructor
	end interface
contains

	type(implicit_heat_transfer_solver)	function constructor(manager)

		type(data_manager)		,intent(inout)	:: manager
		
		real(dkind)	:: calculation_time
		
		type(field_scalar_cons_pointer)	:: scal_ptr
		type(field_vector_cons_pointer)	:: vect_ptr
		type(field_tensor_cons_pointer)	:: tens_ptr

		integer	,dimension(3,2)	:: cons_allocation_bounds
		
		cons_allocation_bounds		= manager%domain%get_local_utter_cells_bounds()
		
		call manager%create_scalar_field(T	,'temperature_uran'	,'T_u')
		constructor%T%s_ptr			=> T		
		call manager%create_scalar_field(rho	,'denisty_uran'	,'rho_u')
		constructor%rho%s_ptr			=> rho		
		call manager%create_scalar_field(kappa	,'thermal_conductivity_uran'	,'kappa_u')
		constructor%kappa%s_ptr			=> kappa
		call manager%create_scalar_field(Q		,'energy_release'		,'Q')
		constructor%Q%s_ptr				=> Q		

		constructor%mesh%mesh_ptr	=> manager%computational_mesh_pointer%mesh_ptr
		constructor%boundary%bc_ptr => manager%boundary_conditions_pointer%bc_ptr

		constructor%domain				= manager%domain
		constructor%mpi_support			= manager%mpi_communications

		constructor%thermo%thermo_ptr	=> manager%thermophysics%thermo_ptr
		constructor%chem%chem_ptr		=> manager%chemistry%chem_ptr
		
		allocate(constructor%thermal_c_coeff_constant(manager%chemistry%chem_ptr%species_number))
	
		allocate(constructor%B(	cons_allocation_bounds(1,1):cons_allocation_bounds(1,2), &
								cons_allocation_bounds(2,1):cons_allocation_bounds(2,2), &
								cons_allocation_bounds(3,1):cons_allocation_bounds(3,2)))
		allocate(constructor%K(	cons_allocation_bounds(1,1):cons_allocation_bounds(1,2), &
								cons_allocation_bounds(2,1):cons_allocation_bounds(2,2), &
								cons_allocation_bounds(3,1):cons_allocation_bounds(3,2)))
		allocate(constructor%C(	cons_allocation_bounds(1,1):cons_allocation_bounds(1,2), &
								cons_allocation_bounds(2,1):cons_allocation_bounds(2,2), &
								cons_allocation_bounds(3,1):cons_allocation_bounds(3,2)))		
		allocate(constructor%E(	cons_allocation_bounds(1,1):cons_allocation_bounds(1,2), &
								cons_allocation_bounds(2,1):cons_allocation_bounds(2,2), &
								cons_allocation_bounds(3,1):cons_allocation_bounds(3,2)))
		allocate(constructor%V(	cons_allocation_bounds(1,1):cons_allocation_bounds(1,2), &
								cons_allocation_bounds(2,1):cons_allocation_bounds(2,2), &
								cons_allocation_bounds(3,1):cons_allocation_bounds(3,2)))									
		allocate(constructor%F(	cons_allocation_bounds(1,1):cons_allocation_bounds(1,2), &
								cons_allocation_bounds(2,1):cons_allocation_bounds(2,2), &
								cons_allocation_bounds(3,1):cons_allocation_bounds(3,2)))			
								
		allocate(constructor%alpha(	cons_allocation_bounds(1,1):cons_allocation_bounds(1,2), &
									cons_allocation_bounds(2,1):cons_allocation_bounds(2,2), &
									cons_allocation_bounds(3,1):cons_allocation_bounds(3,2)))
		allocate(constructor%beta(	cons_allocation_bounds(1,1):cons_allocation_bounds(1,2), &
									cons_allocation_bounds(2,1):cons_allocation_bounds(2,2), &
									cons_allocation_bounds(3,1):cons_allocation_bounds(3,2)))
		allocate(constructor%gamma(	cons_allocation_bounds(1,1):cons_allocation_bounds(1,2), &
									cons_allocation_bounds(2,1):cons_allocation_bounds(2,2), &
									cons_allocation_bounds(3,1):cons_allocation_bounds(3,2)))		
		allocate(constructor%delta(	cons_allocation_bounds(1,1):cons_allocation_bounds(1,2), &
									cons_allocation_bounds(2,1):cons_allocation_bounds(2,2), &
									cons_allocation_bounds(3,1):cons_allocation_bounds(3,2)))
		allocate(constructor%alpha_(	cons_allocation_bounds(1,1):cons_allocation_bounds(1,2), &
									cons_allocation_bounds(2,1):cons_allocation_bounds(2,2), &
									cons_allocation_bounds(3,1):cons_allocation_bounds(3,2)))
		allocate(constructor%beta_(	cons_allocation_bounds(1,1):cons_allocation_bounds(1,2), &
									cons_allocation_bounds(2,1):cons_allocation_bounds(2,2), &
									cons_allocation_bounds(3,1):cons_allocation_bounds(3,2)))
		allocate(constructor%gamma_(	cons_allocation_bounds(1,1):cons_allocation_bounds(1,2), &
									cons_allocation_bounds(2,1):cons_allocation_bounds(2,2), &
									cons_allocation_bounds(3,1):cons_allocation_bounds(3,2)))		
		allocate(constructor%delta_(	cons_allocation_bounds(1,1):cons_allocation_bounds(1,2), &
									cons_allocation_bounds(2,1):cons_allocation_bounds(2,2), &
									cons_allocation_bounds(3,1):cons_allocation_bounds(3,2)))
								
		constructor%B	= 0.0_dkind
		constructor%K	= 0.0_dkind
		constructor%C	= 0.0_dkind
		constructor%E	= 0.0_dkind
		constructor%V	= 0.0_dkind
		constructor%F	= 0.0_dkind		
		
		constructor%alpha	= 0.0_dkind
		constructor%beta	= 0.0_dkind
		constructor%gamma	= 0.0_dkind
		constructor%delta	= 0.0_dkind
										
		constructor%alpha_	= 0.0_dkind
		constructor%beta_	= 0.0_dkind
		constructor%gamma_	= 0.0_dkind
		constructor%delta_	= 0.0_dkind
		
		
		constructor%kappa%s_ptr%cells		= 0.0_dkind
		constructor%rho%s_ptr%cells			= 0.0_dkind
		constructor%T%s_ptr%cells			= 0.0_dkind		
		
		constructor%time		=	calculation_time

		open(newunit = constructor%mean_T_unit, file = 'Mean_T.dat', status = 'replace', form = 'formatted' )
		open(newunit = constructor%energy_flux_unit, file = 'energy_flux.dat', status = 'replace', form = 'formatted' )
								
	end function

	subroutine	heat_transfer_solver_additional_constructor(this, manager)
		class(implicit_heat_transfer_solver) ,intent(inout)	:: this
		type(data_manager)					,intent(inout)	:: manager
		
		integer		,dimension(3,2)	:: cons_inner_loop
		integer		:: i,j,k		
		
		type(field_scalar_cons_pointer)	:: scal_ptr
		type(field_vector_cons_pointer)	:: vect_ptr
		type(field_tensor_cons_pointer)	:: tens_ptr

		call manager%get_cons_field_pointer_by_name(scal_ptr,vect_ptr,tens_ptr,'temperature')
		this%T_gas%s_ptr			=> scal_ptr%s_ptr		
		call manager%get_cons_field_pointer_by_name(scal_ptr,vect_ptr,tens_ptr,'thermal_conductivity')
		this%kappa_gas%s_ptr		=> scal_ptr%s_ptr
		
		
		cons_inner_loop = this%domain%get_local_inner_cells_bounds()
		
		this%kappa%s_ptr%cells		= 28.0_dkind
		
		do k = cons_inner_loop(3,1),cons_inner_loop(3,2)
		do j = cons_inner_loop(2,1),cons_inner_loop(2,2)
		do i = cons_inner_loop(1,1),cons_inner_loop(1,2)
			if(this%boundary%bc_ptr%bc_markers(i,j,k) == 1) then
		
				this%kappa%s_ptr%cells(i,j,k)		= 28.0_dkind
				this%rho%s_ptr%cells(i,j,k)			= 19050.0_dkind
				this%T%s_ptr%cells(i,j,k)			= 300.0_dkind
				this%Q%s_ptr%cells(i,j,k)			= this%Q%s_ptr%cells(i,j,k)*1e-03_dkind
			!	if(j<8) this%T%s_ptr%cells(i,j,k)	= 5000.0_dkind 
					
			end if
		end do
		end do
		end do		
		
	end subroutine
	
	subroutine	solve_problem(this,time_step)
		class(implicit_heat_transfer_solver) ,intent(inout)	:: this	
		real(dkind)					,intent(in)		:: time_step
		integer	,save	:: iteration

		this%time_step = time_step
		
		if ((mod(iteration,10) == 0).or.(iteration==0)) then
			call this%check_balance()
			call this%get_mean_temperature()
		end if
		
		this%time = this%time + this%time_step
		
		call this%calculate_difference_scheme_coefficients(this%time_step)
		call this%alpha_beta_sweep()
		call this%solve_heat_transfer()

		iteration = iteration + 1
		
	end subroutine
	
	subroutine	solve_heat_transfer(this)

		class(implicit_heat_transfer_solver) ,intent(inout)	:: this

		real(dkind)	:: div_thermo_flux, thermo_flux1, thermo_flux2
		real(dkind)	:: energy_prod_summ
		real(dkind)	:: psi, phi, wall_temperature
		
		real(dkind), dimension (3,3)	:: lame_coeffs		
		
		integer						:: dimensions
		integer		,dimension(3,2)	:: cons_inner_loop
		real(dkind)	,dimension(3)	:: cell_size
		
		character(len=20)		:: boundary_type_name	
		
		integer	:: sign, bound_number
		integer :: i,j,k,plus,dim

		dimensions		= this%domain%get_domain_dimensions()

		cons_inner_loop = this%domain%get_local_inner_cells_bounds()

		cell_size		= this%mesh%mesh_ptr%get_cell_edges_length()

		associate(  rho			=> this%rho%s_ptr		, &
					T			=> this%T%s_ptr			, &
					T_gas		=> this%T_gas%s_ptr		, &
					kappa		=> this%kappa%s_ptr		, &
					kappa_gas	=> this%kappa_gas%s_ptr	, &
					bc			=> this%boundary%bc_ptr , &
					mesh		=> this%mesh%mesh_ptr)

	! Set boundary conditions 
					
			do k = cons_inner_loop(3,1),cons_inner_loop(3,2)
			do j = cons_inner_loop(2,1),cons_inner_loop(2,2)
			do i = cons_inner_loop(1,1),cons_inner_loop(1,2)
				if(bc%bc_markers(i,j,k) == 1) then
					do dim = 1,dimensions
						do plus = 1,2
							sign			= (-1)**plus
							bound_number	= bc%bc_markers(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3))
							if( bound_number /= 1 ) then
								if ( bound_number /= 0 ) then
									boundary_type_name = bc%boundary_types(bound_number)%get_type_name()
									select case(boundary_type_name)
										case('wall')
											if(bc%boundary_types(bound_number)%is_conductive()) then											
												wall_temperature = bc%boundary_types(bound_number)%get_wall_temperature()
												psi = wall_temperature
												phi = this%phi_isoth
											else
												psi = this%psi_adiab
												phi = this%phi_adiab
											end if
										case('symmetry_plane')
											psi = this%psi_adiab
											phi = this%phi_adiab
									end select
								end if
								
								if (bound_number == 0) then
									psi = T_gas%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3))
									phi = 0.0_dkind
								end if
											
								if ((dim == 1).and.(sign == -1)) then
									T%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3)) = (phi * this%delta(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3)) + psi)/(1.0_dkind - phi*this%gamma(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3)))
								end if
								if ((dim == 2).and.(sign == -1))	then
									T%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3)) = (phi * this%delta_(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3)) + psi)/(1.0_dkind - phi*this%gamma_(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3)))
								end if							
							
								if ((dim == 1).and.(sign == 1)) then
									T%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3)) = (phi * this%beta(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3)) + psi)/(1.0_dkind - phi*this%alpha(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3)))
								end if
								if ((dim == 2).and.(sign == 1))	then
									T%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3)) = (phi * this%beta_(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3)) + psi)/(1.0_dkind - phi*this%alpha_(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3)))
								end if								
							end if
						end do
					end do
				end if		
			end do
			end do
			end do						
	
			do k = cons_inner_loop(3,1),cons_inner_loop(3,2)
			do j = cons_inner_loop(2,1),cons_inner_loop(2,2)
			do i = cons_inner_loop(1,1),cons_inner_loop(1,2)
				if(bc%bc_markers(i,j,k) == 1) then
					T%cells(i,j,k) = this%gamma_(i,j-1,k)*T%cells(i,j-1,k) + this%delta_(i,j-1,k)
				end if			
			end do
			end do
			end do

		end associate
	end subroutine

	subroutine	calculate_difference_scheme_coefficients(this,time_step)

		class(implicit_heat_transfer_solver) ,intent(inout)	:: this
		real(dkind)					,intent(in)		:: time_step
		
		real(dkind)                 :: cp, kappa_eff
		
		integer						:: dimensions
		integer		,dimension(3,2)	:: cons_inner_loop
		real(dkind)	,dimension(3)	:: cell_size

		integer	:: sign, bound_number
		integer :: i,j,k,plus,dim

		dimensions		= this%domain%get_domain_dimensions()

		cons_inner_loop = this%domain%get_local_inner_cells_bounds()

		cell_size		= this%mesh%mesh_ptr%get_cell_edges_length()

		associate(  rho			=> this%rho%s_ptr		, &
					T			=> this%T%s_ptr			, &
					Q			=> this%Q%s_ptr			, &
					kappa		=> this%kappa%s_ptr		, &
					kappa_gas	=> this%kappa_gas%s_ptr	, &
					bc			=> this%boundary%bc_ptr , &
					mesh		=> this%mesh%mesh_ptr)
				
		do k = cons_inner_loop(3,1),cons_inner_loop(3,2)
		do j = cons_inner_loop(2,1),cons_inner_loop(2,2)
		do i = cons_inner_loop(1,1),cons_inner_loop(1,2)
			if(bc%bc_markers(i,j,k) == 1) then
				do dim = 1,dimensions
					do plus = 1,2
						sign			= (-1)**plus
						bound_number	= bc%bc_markers(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3))
						if (bound_number == 0) then
							kappa_eff	= 2.0_dkind * kappa%cells(i,j,k) * kappa_gas%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3))/(kappa%cells(i,j,k) + kappa_gas%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3)))
							kappa%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3)) = 2.0_dkind * kappa_eff - kappa%cells(i,j,k)
						end if						
					end do
				end do
			end if
		end do
		end do
		end do
					
					
		do k = cons_inner_loop(3,1),cons_inner_loop(3,2)
		do j = cons_inner_loop(2,1),cons_inner_loop(2,2)
		do i = cons_inner_loop(1,1),cons_inner_loop(1,2)
			if(bc%bc_markers(i,j,k) == 1) then

				cp				= 116.3_dkind	!	this%thermo%thermo_ptr%calculate_mixture_cp(T%cells(i,j,k))
						
				!this%B(i,j,k)	= 0.5_dkind*(kappa%cells(i,j,k) + kappa%cells(i,j-1,k))/cell_size(2)/cell_size(2)	- 0.5_dkind*kappa%cells(i,j,k)/cell_size(2)/mesh%mesh(2,i,j,k) 
				!
				!this%K(i,j,k)	= 0.5_dkind*(kappa%cells(i,j,k) + kappa%cells(i-1,j,k))/cell_size(1)/cell_size(1)
				!
				!this%C(i,j,k)	= 0.5_dkind*(kappa%cells(i+1,j,k) + 2.0_dkind*kappa%cells(i,j,k) + kappa%cells(i-1,j,k))/cell_size(1)/cell_size(1)	+ 0.5_dkind*(kappa%cells(i,j+1,k) + 2.0_dkind*kappa%cells(i,j,k) + kappa%cells(i,j-1,k))/cell_size(2)/cell_size(2) + (rho%cells(i,j,k) * cp)/time_step 
				!
				!this%E(i,j,k)	= 0.5_dkind*(kappa%cells(i,j,k) + kappa%cells(i+1,j,k))/cell_size(1)/cell_size(1)
				!
				!this%V(i,j,k)	= 0.5_dkind*(kappa%cells(i,j,k) + kappa%cells(i,j+1,k))/cell_size(2)/cell_size(2)	+ 0.5_dkind*kappa%cells(i,j,k)/cell_size(2)/mesh%mesh(2,i,j,k)   
				!
				!this%F(i,j,k)	= (rho%cells(i,j,k) * cp * T%cells(i,j,k))/time_step	+ Q%cells(i,j,k)
				
				
				
				this%B(i,j,k)	= 0.5_dkind*(kappa%cells(i,j,k) + kappa%cells(i,j-1,k)) * ( mesh%mesh(2,i,j,k) - 0.5_dkind*cell_size(2))/cell_size(2)/cell_size(2)/ mesh%mesh(2,i,j,k)	!- 0.5_dkind*kappa%cells(i,j,k)/cell_size(2)/mesh%mesh(2,i,j,k) 
				
				this%K(i,j,k)	= 0.5_dkind*(kappa%cells(i,j,k) + kappa%cells(i-1,j,k))/cell_size(1)/cell_size(1)
				
				this%C(i,j,k)	= 0.5_dkind*(kappa%cells(i+1,j,k) + 2.0_dkind*kappa%cells(i,j,k) + kappa%cells(i-1,j,k))/cell_size(1)/cell_size(1)	+ 0.5_dkind*(kappa%cells(i,j+1,k)*(mesh%mesh(2,i,j,k) + 0.5_dkind*cell_size(2)) + 2.0_dkind*kappa%cells(i,j,k)*mesh%mesh(2,i,j,k) + kappa%cells(i,j-1,k)*(mesh%mesh(2,i,j,k) - 0.5_dkind*cell_size(2)))/cell_size(2)/cell_size(2)/mesh%mesh(2,i,j,k) + (rho%cells(i,j,k) * cp)/time_step 
				
				this%E(i,j,k)	= 0.5_dkind*(kappa%cells(i,j,k) + kappa%cells(i+1,j,k))/cell_size(1)/cell_size(1)
				
				this%V(i,j,k)	= 0.5_dkind*(kappa%cells(i,j,k) + kappa%cells(i,j+1,k)) * ( mesh%mesh(2,i,j,k) + 0.5_dkind*cell_size(2))/cell_size(2)/cell_size(2)/ mesh%mesh(2,i,j,k)	!+ 0.5_dkind*kappa%cells(i,j,k)/cell_size(2)/mesh%mesh(2,i,j,k)   
				
				this%F(i,j,k)	= (rho%cells(i,j,k) * cp * T%cells(i,j,k))/time_step	+ Q%cells(i,j,k)				
				
				
			end if			
		end do
		end do
		end do

		end associate
	end subroutine	
	
	subroutine	alpha_beta_sweep(this)

		class(implicit_heat_transfer_solver) ,intent(inout)	:: this
		
		real(dkind)	:: phi, psi, wall_temperature
		
		integer						:: dimensions
		integer		,dimension(3,2)	:: cons_inner_loop
		real(dkind)	,dimension(3)	:: cell_size
		
		character(len=20)		:: boundary_type_name
		
		integer	:: iterations = 100

		integer	:: iter
		integer	:: sign, bound_number
		integer :: i,j,k,plus,dim


		dimensions		= this%domain%get_domain_dimensions()

		cons_inner_loop = this%domain%get_local_inner_cells_bounds()

		cell_size		= this%mesh%mesh_ptr%get_cell_edges_length()

		associate(  rho			=> this%rho%s_ptr		, &
					T			=> this%T%s_ptr			, &
					T_gas		=> this%T_gas%s_ptr		, &
					Q			=> this%Q%s_ptr			, &
					kappa		=> this%kappa%s_ptr		, &
					bc			=> this%boundary%bc_ptr , &
					mesh		=> this%mesh%mesh_ptr)
				
		! Set boundary conditions 
					
		do k = cons_inner_loop(3,1),cons_inner_loop(3,2)
		do j = cons_inner_loop(2,1),cons_inner_loop(2,2)
		do i = cons_inner_loop(1,1),cons_inner_loop(1,2)
			if(bc%bc_markers(i,j,k) == 1) then
				do dim = 1,dimensions
					do plus = 1,2
						sign			= (-1)**plus
						bound_number	= bc%bc_markers(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3))
						if( bound_number /= 1 ) then
							if ( bound_number /= 0 ) then
								boundary_type_name = bc%boundary_types(bound_number)%get_type_name()
								select case(boundary_type_name)
									case('wall')
										if(bc%boundary_types(bound_number)%is_conductive()) then											
											wall_temperature = bc%boundary_types(bound_number)%get_wall_temperature()
											psi = wall_temperature
											phi = this%phi_isoth
										else
											psi = this%psi_adiab
											phi = this%phi_adiab
										end if
									case('symmetry_plane')
										psi = this%psi_adiab
										phi = this%phi_adiab
								end select						
							end if

							if (bound_number == 0) then
								psi = T_gas%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3))
								phi = 0.0_dkind
							end if							
							
							if ((dim == 1).and.(sign == -1)) then
								this%alpha(i,j,k)	=	phi
								this%beta(i,j,k)	=	psi
							end if
							if ((dim == 2).and.(sign == -1))	then
								this%alpha_(i,j,k)	=	phi
								this%beta_(i,j,k)	=	psi
							end if							
							
							if ((dim == 1).and.(sign == 1)) then
								this%gamma(i,j,k)	=	phi
								this%delta(i,j,k)	=	psi
							end if
							if ((dim == 2).and.(sign == 1))	then
								this%gamma_(i,j,k)	=	phi
								this%delta_(i,j,k)	=	psi
							end if									
						end if	
					end do
				end do
			end if		
		end do
		end do
		end do			
		
		! First sweep for alpha and gamma
					
		do iter = 1, iterations
					
			do k = cons_inner_loop(3,1),cons_inner_loop(3,2)
			do j = cons_inner_loop(2,1),cons_inner_loop(2,2)
			do i = cons_inner_loop(1,1),cons_inner_loop(1,2)
				if(bc%bc_markers(i,j,k) == 1) then
						this%alpha(i+1,j,k) = this%E(i,j,k)/(this%C(i,j,k) - this%alpha(i,j,k)*this%K(i,j,k) - this%alpha_(i,j,k)*this%B(i,j,k) - this%gamma_(i,j,k)*this%V(i,j,k))
				end if			
			end do
			end do
			end do
			
			do k = cons_inner_loop(3,1),cons_inner_loop(3,2)
			do j = cons_inner_loop(2,1),cons_inner_loop(2,2)
			do i = cons_inner_loop(1,2),cons_inner_loop(1,1),-1
				if(bc%bc_markers(i,j,k) == 1) then
						this%gamma(i-1,j,k) = this%K(i,j,k)/(this%C(i,j,k) - this%gamma(i,j,k)*this%E(i,j,k) - this%alpha_(i,j,k)*this%B(i,j,k) - this%gamma_(i,j,k)*this%V(i,j,k))
				end if			
			end do
			end do
			end do			
		
			do k = cons_inner_loop(3,1),cons_inner_loop(3,2)
			do j = cons_inner_loop(2,1),cons_inner_loop(2,2)
			do i = cons_inner_loop(1,1),cons_inner_loop(1,2)
				if(bc%bc_markers(i,j,k) == 1) then
						this%alpha_(i,j+1,k) = this%V(i,j,k)/(this%C(i,j,k) - this%alpha_(i,j,k)*this%B(i,j,k) - this%alpha(i,j,k)*this%K(i,j,k) - this%gamma(i,j,k)*this%E(i,j,k))
				end if			
			end do
			end do
			end do
			
			do k = cons_inner_loop(3,1),cons_inner_loop(3,2)
			do j = cons_inner_loop(2,2),cons_inner_loop(2,1),-1
			do i = cons_inner_loop(1,1),cons_inner_loop(1,2)
				if(bc%bc_markers(i,j,k) == 1) then
						this%gamma_(i,j-1,k) = this%B(i,j,k)/(this%C(i,j,k) - this%gamma_(i,j,k)*this%V(i,j,k) - this%alpha(i,j,k)*this%K(i,j,k) - this%gamma(i,j,k)*this%E(i,j,k))	
				end if			
			end do
			end do
			end do			
		end do
		
		! Second sweep for beta and delta
					
		do iter = 1, iterations
					
			do k = cons_inner_loop(3,1),cons_inner_loop(3,2)
			do j = cons_inner_loop(2,1),cons_inner_loop(2,2)
			do i = cons_inner_loop(1,1),cons_inner_loop(1,2)
				if(bc%bc_markers(i,j,k) == 1) then
						this%beta(i+1,j,k)	= (this%beta(i,j,k)*this%K(i,j,k) + this%beta_(i,j,k)*this%B(i,j,k) + this%delta_(i,j,k)*this%V(i,j,k) + this%F(i,j,k))	&
											 /(this%C(i,j,k) - this%alpha(i,j,k)*this%K(i,j,k) - this%alpha_(i,j,k)*this%B(i,j,k) - this%gamma_(i,j,k)*this%V(i,j,k))
				end if			
			end do
			end do
			end do
			
			do k = cons_inner_loop(3,1),cons_inner_loop(3,2)
			do j = cons_inner_loop(2,1),cons_inner_loop(2,2)
			do i = cons_inner_loop(1,2),cons_inner_loop(1,1),-1
				if(bc%bc_markers(i,j,k) == 1) then
						this%delta(i-1,j,k)	= (this%delta(i,j,k)*this%E(i,j,k) + this%beta_(i,j,k)*this%B(i,j,k) + this%delta_(i,j,k)*this%V(i,j,k) + this%F(i,j,k))	&
											 /(this%C(i,j,k) - this%gamma(i,j,k)*this%E(i,j,k) - this%alpha_(i,j,k)*this%B(i,j,k) - this%gamma_(i,j,k)*this%V(i,j,k))
				end if			
			end do
			end do
			end do			
		
			do k = cons_inner_loop(3,1),cons_inner_loop(3,2)
			do j = cons_inner_loop(2,1),cons_inner_loop(2,2)
			do i = cons_inner_loop(1,1),cons_inner_loop(1,2)
				if(bc%bc_markers(i,j,k) == 1) then
						this%beta_(i,j+1,k) = (this%beta_(i,j,k)*this%B(i,j,k) + this%beta(i,j,k)*this%K(i,j,k) + this%delta(i,j,k)*this%E(i,j,k) + this%F(i,j,k))	&
											 /(this%C(i,j,k) - this%alpha_(i,j,k)*this%B(i,j,k) - this%alpha(i,j,k)*this%K(i,j,k) - this%gamma(i,j,k)*this%E(i,j,k))
				end if			
			end do
			end do
			end do
			
			do k = cons_inner_loop(3,1),cons_inner_loop(3,2)
			do j = cons_inner_loop(2,2),cons_inner_loop(2,1),-1
			do i = cons_inner_loop(1,1),cons_inner_loop(1,2)
				if(bc%bc_markers(i,j,k) == 1) then
						this%delta_(i,j-1,k) = (this%delta_(i,j,k)*this%V(i,j,k) + this%beta(i,j,k)*this%K(i,j,k) + this%delta(i,j,k)*this%E(i,j,k) + this%F(i,j,k))	&
											 /(this%C(i,j,k) - this%gamma_(i,j,k)*this%V(i,j,k) - this%alpha(i,j,k)*this%K(i,j,k) - this%gamma(i,j,k)*this%E(i,j,k))	
				end if			
			end do
			end do
			end do			
			
		end do		
		

		end associate
	end subroutine		

	subroutine check_balance(this)
	
		class(implicit_heat_transfer_solver) ,intent(inout)	:: this
		
		real(dkind)	:: phi, psi
		
		integer						:: dimensions
		integer		,dimension(3,2)	:: cons_inner_loop
		real(dkind)	,dimension(3)	:: cell_size
		
		real(dkind)	:: energy_input
		real(dkind)	:: energy_outflow
		real(dkind)	:: wall_temperature
		real(dkind)	:: kappa_eff
		
		character(len=20)		:: boundary_type_name

		integer	:: iter
		integer	:: sign, bound_number
		integer :: i,j,k,plus,dim

		dimensions		= this%domain%get_domain_dimensions()

		cons_inner_loop = this%domain%get_local_inner_cells_bounds()

		cell_size		= this%mesh%mesh_ptr%get_cell_edges_length()

		associate(  rho			=> this%rho%s_ptr		, &
					T			=> this%T%s_ptr			, &
					T_gas		=> this%T_gas%s_ptr		, &
					Q			=> this%Q%s_ptr			, &
					kappa		=> this%kappa%s_ptr		, &
					kappa_gas	=> this%kappa_gas%s_ptr	, &
					bc			=> this%boundary%bc_ptr , &
					mesh		=> this%mesh%mesh_ptr)
				
		! Set boundary conditions 
					
		energy_input	= 0.0_dkind			
		energy_outflow	= 0.0_dkind
		
		do k = cons_inner_loop(3,1),cons_inner_loop(3,2)
		do j = cons_inner_loop(2,1),cons_inner_loop(2,2)
		do i = cons_inner_loop(1,1),cons_inner_loop(1,2)
			if(bc%bc_markers(i,j,k) == 1) then
			
				energy_input = energy_input + 2.0_dkind * pi * mesh%mesh(2,i,j,k) * cell_size(1) * cell_size(1) * Q%cells(i,j,k)	! Q_ij * 2*pi*r*dr*dz
			
				do dim = 1,dimensions
					do plus = 1,2
						sign			= (-1)**plus
						bound_number	= bc%bc_markers(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3))
						if(( bound_number /= 1 ).and.( bound_number /= 0 )) then
							boundary_type_name = bc%boundary_types(bound_number)%get_type_name()
							select case(boundary_type_name)
								case('wall')
									if(bc%boundary_types(bound_number)%is_conductive()) then											
										wall_temperature = bc%boundary_types(bound_number)%get_wall_temperature()
										if (dim == 1) then
											energy_outflow = energy_outflow + 2.0_dkind * pi * mesh%mesh(2,i,j,k) * kappa%cells(i,j,k) * (T%cells(i,j,k) - wall_temperature)
										end if
										if (dim == 2) then
											energy_outflow = energy_outflow + 2.0_dkind * pi * mesh%mesh(2,i,j,k) * kappa%cells(i,j,k) * (T%cells(i,j,k) - wall_temperature)
										end if			
									end if
							end select	
						end if
						if(( bound_number == 0 )) then
							kappa_eff		= 2.0_dkind * kappa%cells(i,j,k) * kappa_gas%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3))/(kappa%cells(i,j,k) + kappa_gas%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3)))
						!	print *, kappa_eff, kappa%cells(i,j,k), kappa_gas%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3))
							if (dim == 1) then
								energy_outflow = energy_outflow + 2.0_dkind * pi * mesh%mesh(2,i,j,k) * kappa_eff * (T%cells(i,j,k) - T_gas%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3)))
								continue
							end if
							if (dim == 2) then
								energy_outflow = energy_outflow + 2.0_dkind * pi * mesh%mesh(2,i,j,k) * kappa_eff * (T%cells(i,j,k) - T_gas%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3)))
								continue
							end if							
						end if
					end do
				end do
			end if
							
		end do
		end do
		end do

		end associate
		
		write(this%energy_flux_unit,'(3(E14.6))') this%get_time(), energy_input, energy_outflow
		
		print *, 'Energy balance'
		print *, 'Energy input:' , energy_input
		print *, 'Energy outflow:', energy_outflow 
	!	pause
	
	end subroutine
	
	subroutine get_mean_temperature(this)
	
		class(implicit_heat_transfer_solver) ,intent(inout)	:: this
		
		real(dkind)	:: phi, psi
		real(dkind)	:: radius, height
		
		integer						:: dimensions
		integer		,dimension(3,2)	:: cons_inner_loop
		real(dkind)	,dimension(3)	:: cell_size
		
		real(dkind)	:: mean_T_radius
		real(dkind)	:: mean_T_volume
		
		character(len=20)		:: boundary_type_name

		integer	:: iter
		integer	:: sign, bound_number
		integer :: i,j,k,plus,dim

		dimensions		= this%domain%get_domain_dimensions()

		cons_inner_loop = this%domain%get_local_inner_cells_bounds()

		cell_size		= this%mesh%mesh_ptr%get_cell_edges_length()

		associate(  rho			=> this%rho%s_ptr		, &
					T			=> this%T%s_ptr			, &
					Q			=> this%Q%s_ptr			, &
					kappa		=> this%kappa%s_ptr		, &
					bc			=> this%boundary%bc_ptr , &
					mesh		=> this%mesh%mesh_ptr)
				
		! Set boundary conditions 
					
		mean_T_radius	= 0.0_dkind	
		mean_T_volume	= 0.0_dkind
		
		radius = 12.0_dkind * cell_size(1)
		height = 16.0_dkind * cell_size(1) 
		
	!	do k = cons_inner_loop(3,1),cons_inner_loop(3,2)
		do j = 1,12
		do i = 25,40
			if(bc%bc_markers(i,1,1) == 1) then
			
				mean_T_radius = mean_T_radius + cell_size(1) * T%cells(i,j,1) * cell_size(1) / (radius * height)
				mean_T_volume = mean_T_volume + (mesh%mesh(2,i,j,1)) * T%cells(i,j,1)* cell_size(1) * cell_size(1)  * (2.0_dkind/radius**2/height) 	! T_ij * r*dr*dz * (2/R**2/H)	!
			
			end if
							
		end do
		end do
	!	end do

		end associate
		
		write(this%mean_T_unit,'(3(E14.6))') this%get_time(), mean_T_radius, mean_T_volume

	end subroutine	
	
	pure function get_time_step(this)
		real(dkind)						:: get_time_step
		class(implicit_heat_transfer_solver)	,intent(in)		:: this

		get_time_step = this%time_step
	end function

	pure function get_time(this)
		real(dkind)						:: get_time
		class(implicit_heat_transfer_solver)	,intent(in)		:: this

		get_time = this%time
	end function

	
end module
