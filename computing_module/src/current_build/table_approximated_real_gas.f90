module table_approximated_real_gas_class

	use kind_parameters
	use global_data

	use field_scalar_class
	use field_vector_class
	use field_pointers
	use data_manager_class
	use computational_domain_class
	use chemical_properties_class
	use boundary_conditions_class
	use thermophysical_properties_class
	use computational_mesh_class
    	use riemann_solver_class

	implicit none

#ifdef OMP
	include "omp_lib.h"
#endif

	private
	public  table_approximated_real_gas, table_approximated_real_gas_c

	type(field_scalar_cons)			,target	:: v_s, gamma, h_s, c_v_old, e_i_old, dp_stat_dt, mixture_cp
	type(field_scalar_flow)			,target	:: gamma_f, c_v_f_old, e_i_f_old

	real(dkind)	,dimension(:)	,allocatable	:: concs
	!$omp threadprivate(concs)
    
    type(riemann_solver)						:: riemann
	!$omp threadprivate(riemann)
    
	type 	:: table_approximated_real_gas
		type(field_scalar_cons_pointer)				:: p, p_stat, p_stat_old, T, e_i, E_f, rho, mol_mix_conc, v_s, gamma, h_s,dp_stat_dt, mixture_cp
		type(field_scalar_flow_pointer)				:: rho_f, p_f, e_i_f, v_s_f, E_f_f, T_f, gamma_f
		type(field_vector_cons_pointer)				:: v, Y
		type(field_vector_flow_pointer)				:: v_f, Y_f
		type(computational_domain)					:: domain
		type(thermophysical_properties_pointer)		:: thermo
		type(chemical_properties_pointer)			:: chem
		type(computational_mesh_pointer)			:: mesh
		type(boundary_conditions_pointer)			:: boundary

		real(dkind)									:: total_energy, total_mass
		real(dkind)	,dimension(:)	,allocatable	:: total_species
		
	contains
		procedure				:: apply_state_equation
        procedure               :: apply_state_equation_low_mach	
		procedure               :: apply_state_equation_low_mach_fds
		procedure				:: apply_state_equation_flow_variables
		procedure				:: apply_state_equation_for_initial_conditions
		procedure				:: apply_state_equation_flow_variables_for_IC
		procedure				:: apply_boundary_conditions_for_initial_conditions
		procedure				:: check_conservation_laws
	end type

	interface table_approximated_real_gas_c
		module procedure constructor
	end interface

contains

	type(table_approximated_real_gas) function constructor (manager)

		type(data_manager)						,intent(inout)	:: manager

		type(field_scalar_cons_pointer)	:: scal_c_ptr
		type(field_vector_cons_pointer)	:: vect_c_ptr
		type(field_tensor_cons_pointer)	:: tens_c_ptr
		
		type(field_scalar_flow_pointer)	:: scal_f_ptr
		type(field_vector_flow_pointer)	:: vect_f_ptr		

		integer	:: species_number

		call manager%get_cons_field_pointer_by_name(scal_c_ptr,vect_c_ptr,tens_c_ptr,'temperature')
		constructor%T%s_ptr						=> scal_c_ptr%s_ptr
		call manager%get_cons_field_pointer_by_name(scal_c_ptr,vect_c_ptr,tens_c_ptr,'pressure')
		constructor%p%s_ptr						=> scal_c_ptr%s_ptr
		call manager%get_cons_field_pointer_by_name(scal_c_ptr,vect_c_ptr,tens_c_ptr,'pressure_static')
		constructor%p_stat%s_ptr				=> scal_c_ptr%s_ptr		
		call manager%get_cons_field_pointer_by_name(scal_c_ptr,vect_c_ptr,tens_c_ptr,'pressure_static_old')
		constructor%p_stat_old%s_ptr				=> scal_c_ptr%s_ptr				
		call manager%get_cons_field_pointer_by_name(scal_c_ptr,vect_c_ptr,tens_c_ptr,'density')
		constructor%rho%s_ptr					=> scal_c_ptr%s_ptr
		call manager%get_cons_field_pointer_by_name(scal_c_ptr,vect_c_ptr,tens_c_ptr,'internal_energy')
		constructor%e_i%s_ptr					=> scal_c_ptr%s_ptr
		call manager%get_cons_field_pointer_by_name(scal_c_ptr,vect_c_ptr,tens_c_ptr,'full_energy')
		constructor%E_f%s_ptr					=> scal_c_ptr%s_ptr
		call manager%get_cons_field_pointer_by_name(scal_c_ptr,vect_c_ptr,tens_c_ptr,'mixture_molar_concentration')
		constructor%mol_mix_conc%s_ptr			=> scal_c_ptr%s_ptr
		call manager%get_cons_field_pointer_by_name(scal_c_ptr,vect_c_ptr,tens_c_ptr,'pressure_static_change')
		constructor%dp_stat_dt%s_ptr			=> scal_c_ptr%s_ptr		
		
		call manager%get_flow_field_pointer_by_name(scal_f_ptr,vect_f_ptr,'density_flow')
		constructor%rho_f%s_ptr					=> scal_f_ptr%s_ptr
		call manager%get_flow_field_pointer_by_name(scal_f_ptr,vect_f_ptr,'pressure_flow')
		constructor%p_f%s_ptr					=> scal_f_ptr%s_ptr
		call manager%get_flow_field_pointer_by_name(scal_f_ptr,vect_f_ptr,'full_energy_flow')
		constructor%E_f_f%s_ptr					=> scal_f_ptr%s_ptr			
		call manager%get_flow_field_pointer_by_name(scal_f_ptr,vect_f_ptr,'internal_energy_flow')
		constructor%e_i_f%s_ptr					=> scal_f_ptr%s_ptr
		call manager%get_flow_field_pointer_by_name(scal_f_ptr,vect_f_ptr,'velocity_of_sound_flow')
		constructor%v_s_f%s_ptr					=> scal_f_ptr%s_ptr
		call manager%get_flow_field_pointer_by_name(scal_f_ptr,vect_f_ptr,'temperature_flow')
		constructor%T_f%s_ptr					=> scal_f_ptr%s_ptr		
		
		constructor%boundary%bc_ptr => manager%boundary_conditions_pointer%bc_ptr

		call manager%create_scalar_field(mixture_cp,'mixture_cp', 'cp')
		constructor%mixture_cp%s_ptr		=> mixture_cp			
		call manager%create_scalar_field(v_s	,'velocity_of_sound','v_s'	)
		constructor%v_s%s_ptr			=> v_s
		call manager%create_scalar_field(gamma	,'adiabatic_index'	,'gamma')
		constructor%gamma%s_ptr			=> gamma
		call manager%create_scalar_field(h_s	,'sensible_enthalpy'	,'h_s')
		constructor%h_s%s_ptr		=> h_s

		call manager%get_cons_field_pointer_by_name(scal_c_ptr,vect_c_ptr,tens_c_ptr,'velocity')
		constructor%v%v_ptr				=> vect_c_ptr%v_ptr
		call manager%get_cons_field_pointer_by_name(scal_c_ptr,vect_c_ptr,tens_c_ptr,'specie_molar_concentration')
		constructor%Y%v_ptr				=> vect_c_ptr%v_ptr

		call manager%get_flow_field_pointer_by_name(scal_f_ptr,vect_f_ptr,'velocity_flow')
		constructor%v_f%v_ptr			=> vect_f_ptr%v_ptr		
		call manager%get_flow_field_pointer_by_name(scal_f_ptr,vect_f_ptr,'specie_molar_concentration_flow')
		constructor%Y_f%v_ptr			=> vect_f_ptr%v_ptr		
		
		call manager%create_scalar_field(gamma_f	,'adiabatic_index_flow'	,'gamma_flow')
		constructor%gamma_f%s_ptr		=> gamma_f

		
		constructor%domain				= manager%domain
		constructor%thermo%thermo_ptr	=> manager%thermophysics%thermo_ptr
		constructor%chem%chem_ptr		=> manager%chemistry%chem_ptr
		constructor%mesh%mesh_ptr		=> manager%computational_mesh_pointer%mesh_ptr
        
        !$omp parallel
        riemann						=	riemann_solver_c(manager)
		!$omp end parallel
        
		species_number = manager%chemistry%chem_ptr%species_number
		allocate(constructor%total_species(species_number))
		constructor%total_energy	= 0.0_dkind
		constructor%total_mass		= 0.0_dkind
		constructor%total_species	= 0.0_dkind
		
		!$omp parallel
		allocate(concs(species_number))
		concs = 0.0_dkind
		!$omp end parallel

	end function

	subroutine apply_state_equation_for_initial_conditions(this)

		class(table_approximated_real_gas) ,intent(inout) :: this

		real(dkind)                     :: cp, cv

		real(dkind) ,dimension(this%chem%chem_ptr%species_number)    :: concs
		real(dkind)	:: average_molar_mass
		real(dkind)	:: farfield_mol_mix_conc

		integer	:: bound_number
		integer	:: species_number

		integer	:: dimensions
		integer	,dimension(3,2)	:: utter_loop, inner_loop
		character(len=20)		:: boundary_type_name
		real(dkind)				:: farfield_density, farfield_pressure

		integer	:: specie_number
		integer	:: sign
		integer	:: i,j,k,dim,plus

		species_number = this%chem%chem_ptr%species_number

		associate(		T               => this%T%s_ptr				, &
						mixture_cp		=> this%mixture_cp%s_ptr	, &
						p               => this%p%s_ptr				, &
						rho             => this%rho%s_ptr			, &
						e_i             => this%e_i%s_ptr			, &
						E_f             => this%E_f%s_ptr			, &
						h_s				=> this%h_s%s_ptr			, &
						mol_mix_conc    => this%mol_mix_conc%s_ptr	, &
						v_s             => this%v_s%s_ptr			, &
						Y               => this%Y%v_ptr				, &
						v               => this%v%v_ptr)

		dimensions		= this%domain%get_domain_dimensions()
		utter_loop		= this%domain%get_local_utter_cells_bounds()
		inner_loop		= this%domain%get_local_inner_cells_bounds()

		do k = inner_loop(3,1),inner_loop(3,2)
		do j = inner_loop(2,1),inner_loop(2,2)
		do i = inner_loop(1,1),inner_loop(1,2)

			if(this%boundary%bc_ptr%bc_markers(i,j,k) == 0) then
			
				average_molar_mass = 0.0_dkind
				do specie_number = 1,species_number
					if (this%thermo%thermo_ptr%molar_masses(specie_number) /= 0.0_dkind) then
						average_molar_mass = average_molar_mass + Y%pr(specie_number)%cells(i,j,k) / this%thermo%thermo_ptr%molar_masses(specie_number)
					end if
				end do
				
				average_molar_mass =  1.0_dkind / average_molar_mass
				mol_mix_conc%cells(i,j,k)	= average_molar_mass
				
				concs = 0.0_dkind
				do specie_number = 1,species_number
					if (this%thermo%thermo_ptr%molar_masses(specie_number) /= 0.0_dkind) then
						concs(specie_number)		=	Y%pr(specie_number)%cells(i,j,k) / this%thermo%thermo_ptr%molar_masses(specie_number) * mol_mix_conc%cells(i,j,k)
					end if
				end do
				
				rho%cells(i,j,k)		= p%cells(i,j,k) / (T%cells(i,j,k) * r_gase_J) * mol_mix_conc%cells(i,j,k)
				
				cp = this%thermo%thermo_ptr%calculate_mixture_cp(T%cells(i,j,k), concs)
				cv = cp - r_gase_J
				
				gamma%cells(i,j,k)		= cp / cv	
				
				e_i%cells(i,j,k)		= this%thermo%thermo_ptr%calculate_mixture_energy(T%cells(i,j,k), concs) - this%thermo%thermo_ptr%calculate_mixture_enthalpy(298.15_dkind, concs)
				h_s%cells(i,j,k)		= (this%thermo%thermo_ptr%calculate_mixture_enthalpy(T%cells(i,j,k), concs) - this%thermo%thermo_ptr%calculate_mixture_enthalpy(298.15_dkind, concs)) / mol_mix_conc%cells(i,j,k)
				
				E_f%cells(i,j,k)		= e_i%cells(i,j,k) / mol_mix_conc%cells(i,j,k)			

                
        			!# P = const 
               			!E_f%cells(i,j,k)		= h_s%cells(i,j,k)

				do dim = 1,dimensions
					E_f%cells(i,j,k)	= E_f%cells(i,j,k) + 0.5_dkind * ( v%pr(dim)%cells(i,j,k) * v%pr(dim)%cells(i,j,k) )
				end do

				v_s%cells(i,j,k)		= sqrt(gamma%cells(i,j,k)*p%cells(i,j,k)/rho%cells(i,j,k))				
				
				mixture_cp%cells(i,j,k)	= cp
				
				this%total_energy	= this%total_energy	+ E_f%cells(i,j,k)
				this%total_mass		= this%total_mass	+ rho%cells(i,j,k)
				do specie_number = 1,species_number
					this%total_species(specie_number)	=	this%total_species(specie_number) + Y%pr(specie_number)%cells(i,j,k) * rho%cells(i,j,k)
				end do

			end if
		end do
		end do
		end do

		call this%apply_boundary_conditions_for_initial_conditions()
		
		end associate

	end subroutine

	subroutine apply_state_equation_flow_variables_for_IC(this)

		class(table_approximated_real_gas) ,intent(inout) :: this

		real(dkind)	:: velocity, T_old, e_i_old, t_initial, t_final, h_s, e_internal, mol_mix_conc, cp, cv, gamma
		real(dkind)	:: average_molar_mass
		real(dkind)	:: rho_Y
		
		real(dkind)	:: old_Mach
		
		integer	:: species_number
		integer	:: dimensions
		integer	,dimension(3,2)	:: flow_inner_loop, loop

		integer	:: specie_number, dim, dim1
		integer	:: i,j,k

		species_number = this%chem%chem_ptr%species_number
		
		dimensions		= this%domain%get_domain_dimensions()
		flow_inner_loop	= this%domain%get_local_inner_faces_bounds()		
		
		associate( 	p_f     => this%p_f%s_ptr		, &
					rho_f	=> this%rho_f%s_ptr		, &
					v_s_f	=> this%v_s_f%s_ptr		, &
					v_f		=> this%v_f%v_ptr		, &
					Y_f		=> this%Y_f%v_ptr		, &
					E_f_f	=> this%E_f_f%s_ptr		, &
					e_i_f	=> this%e_i_f%s_ptr		, &
					T_f		=> this%T_f%s_ptr		, &
					gamma_f	=> this%gamma_f%s_ptr	, &
					bc		=> this%boundary%bc_ptr)

	!$omp parallel default(none) private(i,j,k,dim,dim1,specie_number,loop,average_molar_mass,t_initial,cp,cv,h_s) , &
	!$omp& firstprivate(this)	,&
	!$omp& shared(T_f,p_f,rho_f,e_i_f,e_i_f_old,v_s_f,E_f_f,v_f,Y_f,c_v_f_old,gamma_f,flow_inner_loop,dimensions,species_number,bc)
	
		do dim = 1, dimensions

			loop = flow_inner_loop

			do dim1 = 1, dimensions
				loop(dim1,2) = flow_inner_loop(dim1,2) - (1 - I_m(dim1,dim))	
			end do		

		!$omp do collapse(3) schedule(guided)

			do k = loop(3,1),loop(3,2)
			do j = loop(2,1),loop(2,2)
			do i = loop(1,1),loop(1,2)
				if ((bc%bc_markers(i,j,k) == 0).or.(bc%bc_markers(i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3)) == 0)) then

					average_molar_mass = 0.0_dkind
					do specie_number = 1,species_number
						average_molar_mass = average_molar_mass + Y_f%pr(specie_number)%cells(dim,i,j,k) / this%thermo%thermo_ptr%molar_masses(specie_number)
					end do
				
					average_molar_mass =  1.0_dkind / average_molar_mass		
					
					do specie_number = 1,species_number
						concs(specie_number)				= Y_f%pr(specie_number)%cells(dim,i,j,k) * average_molar_mass / this%thermo%thermo_ptr%molar_masses(specie_number)
					end do						

					T_f%cells(dim,i,j,k)	= p_f%cells(dim,i,j,k) / rho_f%cells(dim,i,j,k) / r_gase_J * average_molar_mass
					t_initial				= T_f%cells(dim,i,j,k)
				
					cp = this%thermo%thermo_ptr%calculate_mixture_cp(t_initial, concs)
					cv = cp - r_gase_J 

					gamma_f%cells(dim,i,j,k)	= cp / cv

					h_s						= (this%thermo%thermo_ptr%calculate_mixture_enthalpy(T_f%cells(dim,i,j,k), concs) - this%thermo%thermo_ptr%calculate_mixture_enthalpy(298.15_dkind, concs)) 
					e_i_f%cells(dim,i,j,k)	= this%thermo%thermo_ptr%calculate_mixture_energy(T_f%cells(dim,i,j,k), concs) - this%thermo%thermo_ptr%calculate_mixture_enthalpy(298.15_dkind, concs)
				
					v_s_f%cells(dim,i,j,k) = sqrt(gamma_f%cells(dim,i,j,k)*p_f%cells(dim,i,j,k)/rho_f%cells(dim,i,j,k))
					
					E_f_f%cells(dim,i,j,k) = e_i_f%cells(dim,i,j,k) / average_molar_mass

					do dim1 = 1, dimensions
						E_f_f%cells(dim,i,j,k) = E_f_f%cells(dim,i,j,k) + 0.5_dkind*v_f%pr(dim1)%cells(dim,i,j,k)*v_f%pr(dim1)%cells(dim,i,j,k)
					end do
				end if	
			end do
			end do
			end do

			!$omp end do nowait
		end do
			
		!$omp end parallel

		continue

		end associate

	end subroutine	

	subroutine apply_state_equation(this)
		class(table_approximated_real_gas) ,intent(inout) :: this

		real(dkind)	:: velocity, t_initial, t_final, e_internal, cp, cv
		real(dkind)	:: average_molar_mass

		integer	:: species_number
		integer	:: dimensions
		integer	,dimension(3,2)	:: cons_inner_loop
		
		integer	:: specie_number, dim, dim2, T_iter
		integer	:: i,j,k

		species_number	= this%chem%chem_ptr%species_number
		
		dimensions		= this%domain%get_domain_dimensions()
		cons_inner_loop	= this%domain%get_local_inner_cells_bounds()		
		
		associate(  T               => this%T%s_ptr				, &
					mixture_cp		=> this%mixture_cp%s_ptr	, &
					p               => this%p%s_ptr				, &
					rho             => this%rho%s_ptr          , &
					E_f				=> this%E_f%s_ptr		   , &
					e_i             => this%e_i%s_ptr          , &
					h_s             => this%h_s%s_ptr          , &
					v_s             => this%v_s%s_ptr          , &
					mol_mix_conc    => this%mol_mix_conc%s_ptr , &
					v				=> this%v%v_ptr			   , &
					Y				=> this%Y%v_ptr				, &
					mesh			=> this%mesh%mesh_ptr)

	!$omp parallel default(none) private(i,j,k,dim,specie_number,t_initial,t_final,e_internal,cp,cv,average_molar_mass,T_iter) , &
	!$omp& firstprivate(this)	,&
	!$omp& shared(T,mixture_cp,p,rho,e_i,h_s,gamma,v_s,mol_mix_conc,E_f,v,Y,dimensions,cons_inner_loop,species_number)

	!$omp do collapse(3) schedule(guided)
		do k = cons_inner_loop(3,1),cons_inner_loop(3,2)
		do j = cons_inner_loop(2,1),cons_inner_loop(2,2)
		do i = cons_inner_loop(1,1),cons_inner_loop(1,2)

			if(this%boundary%bc_ptr%bc_markers(i,j,k) == 0) then

				e_i%cells(i,j,k) = E_f%cells(i,j,k)
				
				do dim = 1,dimensions
					 e_i%cells(i,j,k) = e_i%cells(i,j,k) - 0.5_dkind * ( v%pr(dim)%cells(i,j,k) * v%pr(dim)%cells(i,j,k) )
				end do

				average_molar_mass = 0.0_dkind
				do specie_number = 1,species_number
					if (this%thermo%thermo_ptr%molar_masses(specie_number) /= 0.0_dkind) then
						average_molar_mass = average_molar_mass + Y%pr(specie_number)%cells(i,j,k) / this%thermo%thermo_ptr%molar_masses(specie_number)
					end if
				end do
				
				average_molar_mass =  1.0_dkind / average_molar_mass
				mol_mix_conc%cells(i,j,k)	= average_molar_mass

				concs = 0.0_dkind
				do specie_number = 1,species_number
					if (this%thermo%thermo_ptr%molar_masses(specie_number) /= 0.0_dkind) then
						concs(specie_number)	= Y%pr(specie_number)%cells(i,j,k) / this%thermo%thermo_ptr%molar_masses(specie_number) * mol_mix_conc%cells(i,j,k)
					end if
				end do				
				
				e_i%cells(i,j,k) = e_i%cells(i,j,k) * mol_mix_conc%cells(i,j,k)

			!# New de = cv*dT equation of state		
				T%cells(i,j,k)			= this%thermo%thermo_ptr%calculate_temperature(T%cells(i,j,k),e_i%cells(i,j,k),concs)
				p%cells(i,j,k)			= T%cells(i,j,k) * rho%cells(i,j,k) * r_gase_J / mol_mix_conc%cells(i,j,k)

    			!# P = const
            		!	T%cells(i,j,k)			= this%thermo%thermo_ptr%calculate_temperature_Pconst(T%cells(i,j,k),e_i%cells(i,j,k),concs)
			!	rho%cells(i,j,k)		= p%cells(i,j,k) * mol_mix_conc%cells(i,j,k) / T%cells(i,j,k) / r_gase_J
                
				cp						= this%thermo%thermo_ptr%calculate_mixture_cp(T%cells(i,j,k), concs)
				cv						= cp - r_gase_J 
				gamma%cells(i,j,k)		= cp / cv

				h_s%cells(i,j,k)		= (this%thermo%thermo_ptr%calculate_mixture_enthalpy(T%cells(i,j,k), concs) - this%thermo%thermo_ptr%calculate_mixture_enthalpy(298.15_dkind, concs)) / mol_mix_conc%cells(i,j,k)
				mixture_cp%cells(i,j,k) = cp / mol_mix_conc%cells(i,j,k)
	
				if((gamma%cells(i,j,k) < 0.0).or.(p%cells(i,j,k) < 0.0).or.(rho%cells(i,j,k) < 0.0)) then
					print *, 'Cons EOS exception: ', dim, i,j,k
					do dim2 = 1, dimensions
						print *, 'Coordinates', dim2, mesh%mesh(dim2,i,j,k)
					end do
					print *, gamma%cells(i,j,k)
					print *, p%cells(i,j,k)
					print *, rho%cells(i,j,k)
					stop
				end if	
					
				v_s%cells(i,j,k)	= sqrt(gamma%cells(i,j,k)*p%cells(i,j,k)/rho%cells(i,j,k))

			end if

		end do
		end do
		end do

	!$omp end do nowait
	!$omp end parallel

		continue

		end associate

		continue
	end subroutine

	subroutine apply_state_equation_low_mach(this)
		class(table_approximated_real_gas) ,intent(inout) :: this

		real(dkind)	:: velocity, t_initial, t_final, e_internal, cp, cv
		real(dkind)	:: average_molar_mass

		integer	:: species_number
		integer	:: dimensions
		integer	,dimension(3,2)	:: cons_inner_loop
		
		integer	:: specie_number, dim
		integer	:: i,j,k

		species_number = this%chem%chem_ptr%species_number
		
		dimensions		= this%domain%get_domain_dimensions()
		cons_inner_loop	= this%domain%get_local_inner_cells_bounds()		
		
		associate(  T               => this%T%s_ptr            , &
					p               => this%p%s_ptr            , &
					p_stat			=> this%p_stat%s_ptr	   , &
					rho             => this%rho%s_ptr          , &
					E_f				=> this%E_f%s_ptr		   , &
					e_i             => this%e_i%s_ptr          , &
					v_s             => this%v_s%s_ptr          , &
					mol_mix_conc    => this%mol_mix_conc%s_ptr , &
					v				=> this%v%v_ptr			   , &
					Y				=> this%Y%v_ptr			   , &
					bc				=> this%boundary%bc_ptr)

	!$omp parallel default(none) private(i,j,k,dim,specie_number,t_initial,t_final,e_internal,cp,cv,average_molar_mass) , &
	!$omp& firstprivate(this)	,&
	!$omp& shared(T,p,p_stat,rho,e_i,gamma,v_s,mol_mix_conc,E_f,v,Y,dimensions,cons_inner_loop,species_number,bc)

	!$omp do collapse(3) schedule(static)
		do k = cons_inner_loop(3,1),cons_inner_loop(3,2)
		do j = cons_inner_loop(2,1),cons_inner_loop(2,2)
		do i = cons_inner_loop(1,1),cons_inner_loop(1,2)

			if(bc%bc_markers(i,j,k) == 0) then

				e_i%cells(i,j,k) = E_f%cells(i,j,k)
			
				average_molar_mass = 0.0_dkind
				do specie_number = 1,species_number
					average_molar_mass = average_molar_mass + Y%pr(specie_number)%cells(i,j,k) / this%thermo%thermo_ptr%molar_masses(specie_number)
				end do
				
				average_molar_mass =  1.0_dkind / average_molar_mass
				mol_mix_conc%cells(i,j,k)	= average_molar_mass

				do specie_number = 1,species_number
					concs(specie_number)				= Y%pr(specie_number)%cells(i,j,k) * mol_mix_conc%cells(i,j,k) / this%thermo%thermo_ptr%molar_masses(specie_number)
				end do				
				
				t_initial           = T%cells(i,j,k)

				cp = this%thermo%thermo_ptr%calculate_mixture_cp(t_initial, concs)
				cv = cp - r_gase_J 

				gamma%cells(i,j,k)	= cp / cv				
				
				T%cells(i,j,k)      = e_i%cells(i,j,k) * (gamma%cells(i,j,k) - 1.0_dkind) / r_gase_J * mol_mix_conc%cells(i,j,k)
				v_s%cells(i,j,k)	= sqrt(gamma%cells(i,j,k)*p_stat%cells(i,j,k)/rho%cells(i,j,k))
                
			    e_i%cells(i,j,k) = p_stat%cells(i,j,k) / rho%cells(i,j,k) / (gamma%cells(i,j,k) - 1.0_dkind)
            end if

		end do
		end do
		end do

	!$omp end do nowait
	!$omp end parallel

		continue

		end associate

		continue
    end subroutine	

	subroutine apply_state_equation_low_mach_fds(this,time_step,predictor)
		class(table_approximated_real_gas) ,intent(inout) :: this
		real(dkind)			,intent(in)		:: time_step
		logical				,intent(in)		:: predictor
		
		real(dkind)	:: velocity,  e_internal, cp, cv
		real(dkind)	:: average_molar_mass

		integer	:: species_number
		integer	:: dimensions
		integer	,dimension(3,2)	:: cons_inner_loop
		
		integer	:: specie_number, dim
		integer	:: i,j,k

		species_number = this%chem%chem_ptr%species_number
		
		dimensions		= this%domain%get_domain_dimensions()
		cons_inner_loop	= this%domain%get_local_inner_cells_bounds()		
		
		associate(  T               => this%T%s_ptr            , &
					mixture_cp		=> this%mixture_cp%s_ptr	, &
					p               => this%p%s_ptr            , &
					p_stat			=> this%p_stat%s_ptr	   , &
					p_stat_old		=> this%p_stat_old%s_ptr   , &
					dp_stat_dt		=> this%dp_stat_dt%s_ptr   , &
					rho             => this%rho%s_ptr          , &
					E_f				=> this%E_f%s_ptr		   , &
					e_i             => this%e_i%s_ptr          , &
					h_s				=> this%h_s%s_ptr		   , &
					v_s             => this%v_s%s_ptr          , &
					mol_mix_conc    => this%mol_mix_conc%s_ptr , &
					v				=> this%v%v_ptr			   , &
					Y				=> this%Y%v_ptr			   , &
					bc				=> this%boundary%bc_ptr)

	!$omp parallel default(none) private(i,j,k,dim,specie_number,cp,cv,average_molar_mass) , &
	!$omp& firstprivate(this)	,&
	!$omp& shared(T,p,h_s,p_stat,p_stat_old,dp_stat_dt,rho,e_i,gamma,v_s,mixture_cp,mol_mix_conc,E_f,v,Y,dimensions,cons_inner_loop,species_number,predictor,time_step,bc)

	!$omp do collapse(3) schedule(static)					
					
		do k = cons_inner_loop(3,1),cons_inner_loop(3,2)
		do j = cons_inner_loop(2,1),cons_inner_loop(2,2)
		do i = cons_inner_loop(1,1),cons_inner_loop(1,2)

			if(bc%bc_markers(i,j,k) == 0) then
			
				e_i%cells(i,j,k) = E_f%cells(i,j,k)
			
				average_molar_mass = 0.0_dkind
				do specie_number = 1,species_number
					average_molar_mass = average_molar_mass + Y%pr(specie_number)%cells(i,j,k) / this%thermo%thermo_ptr%molar_masses(specie_number)
				end do
				
				average_molar_mass =  1.0_dkind / average_molar_mass
				mol_mix_conc%cells(i,j,k)	= average_molar_mass

				do specie_number = 1,species_number
					concs(specie_number)				= Y%pr(specie_number)%cells(i,j,k) * mol_mix_conc%cells(i,j,k) / this%thermo%thermo_ptr%molar_masses(specie_number)
				end do				

				if(predictor) then
					p_stat_old%cells(i,j,k)	= p_stat%cells(i,j,k)
					p_stat%cells(i,j,k)		= p_stat%cells(i,j,k) + dp_stat_dt%cells(i,j,k) * time_step
				else 
					p_stat%cells(i,j,k)	= 0.5_dkind * (p_stat%cells(i,j,k) + p_stat_old%cells(i,j,k) + dp_stat_dt%cells(i,j,k) * time_step)
				end if					

				T%cells(i,j,k)      = p_stat%cells(i,j,k) / rho%cells(i,j,k) / r_gase_J * mol_mix_conc%cells(i,j,k)
				
				cp = this%thermo%thermo_ptr%calculate_mixture_cp(T%cells(i,j,k), concs)
				cv = cp - r_gase_J 

				gamma%cells(i,j,k)	= cp / cv
		
				h_s%cells(i,j,k)	= (this%thermo%thermo_ptr%calculate_mixture_enthalpy(T%cells(i,j,k), concs) - this%thermo%thermo_ptr%calculate_mixture_enthalpy(298.15_dkind, concs)) / mol_mix_conc%cells(i,j,k)

				mixture_cp%cells(i,j,k)	= cp
				
				p%cells(i,j,k)      = p_stat%cells(i,j,k)
				
			!	print*, p_stat%cells(i,j,k), rho%cells(i,j,k), gamma%cells(i,j,k)
				v_s%cells(i,j,k)	= sqrt(gamma%cells(i,j,k)*p_stat%cells(i,j,k)/rho%cells(i,j,k))
			    e_i%cells(i,j,k)	= p_stat%cells(i,j,k) / rho%cells(i,j,k) / (gamma%cells(i,j,k) - 1.0_dkind)

				do dim = 1,dimensions
					 E_f%cells(i,j,k) = e_i%cells(i,j,k) + 0.5_dkind * ( v%pr(dim)%cells(i,j,k) * v%pr(dim)%cells(i,j,k) )
				end do
				
				
            end if

		end do
		end do
		end do

	!$omp end do nowait
	!$omp end parallel		
		
		continue
		
		end associate		
    end subroutine	

	subroutine apply_state_equation_flow_variables(this)
		class(table_approximated_real_gas) ,intent(inout) :: this

		real(dkind)	:: velocity, T_old, e_i_old, t_initial, t_final, h_s, e_internal, mol_mix_conc, cp, cv, gamma
        real(dkind)	:: rho_l,rho_r,p_l,p_r,gamma_l,gamma_r,v_l,v_r
		real(dkind)	:: average_molar_mass
		real(dkind)	:: rho_Y
		
		real(dkind)	:: old_Mach
		
		integer	:: species_number
		integer	:: dimensions
		integer	,dimension(3,2)	:: flow_inner_loop, loop

		integer	:: specie_number, dim, dim1, dim2
		integer	:: i,j,k

		species_number = this%chem%chem_ptr%species_number
		
		dimensions		= this%domain%get_domain_dimensions()
		flow_inner_loop	= this%domain%get_local_inner_faces_bounds()		
		
		associate( 	p_f     => this%p_f%s_ptr		, &
					rho_f	=> this%rho_f%s_ptr		, &
					v_s_f	=> this%v_s_f%s_ptr		, &
					v_f		=> this%v_f%v_ptr		, &
					Y_f		=> this%Y_f%v_ptr		, &
					E_f_f	=> this%E_f_f%s_ptr		, &
					e_i_f	=> this%e_i_f%s_ptr		, &
					T_f		=> this%T_f%s_ptr		, &
					gamma_f	=> this%gamma_f%s_ptr	, &
					p       => this%p%s_ptr         , &
					rho     => this%rho%s_ptr       , &
					gamma   => this%gamma%s_ptr     , &
					v       => this%v%v_ptr         , &
					mesh	=> this%mesh%mesh_ptr	, &
					bc		=> this%boundary%bc_ptr)

	!$omp parallel default(none) private(i,j,k,dim,dim1,specie_number,loop,average_molar_mass,t_initial,t_final,e_internal,cp,cv,h_s,T_old,rho_l,rho_r,p_l,p_r,gamma_l,gamma_r,v_l,v_r) , &
	!$omp& firstprivate(this)	,&
	!$omp& shared(T_f,p_f,rho_f,e_i_f,e_i_f_old,v_s_f,E_f_f,v_f,Y_f,c_v_f_old,gamma_f,p,rho,gamma,v,flow_inner_loop,dimensions,species_number,bc)
	
		do dim = 1, dimensions

			loop = flow_inner_loop

			do dim1 = 1, dimensions
				loop(dim1,2) = flow_inner_loop(dim1,2) - (1 - I_m(dim1,dim))	
			end do		

		!$omp do collapse(3) schedule(guided)

			do k = loop(3,1),loop(3,2)
			do j = loop(2,1),loop(2,2)
			do i = loop(1,1),loop(1,2)
				if ((bc%bc_markers(i,j,k) == 0).or.(bc%bc_markers(i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3)) == 0)) then

					average_molar_mass = 0.0_dkind
					do specie_number = 1,species_number
						average_molar_mass = average_molar_mass + Y_f%pr(specie_number)%cells(dim,i,j,k) / this%thermo%thermo_ptr%molar_masses(specie_number)
					end do
				
					average_molar_mass =  1.0_dkind / average_molar_mass		
					
					do specie_number = 1,species_number
						concs(specie_number)				= Y_f%pr(specie_number)%cells(dim,i,j,k) * average_molar_mass / this%thermo%thermo_ptr%molar_masses(specie_number)
                    end do						
					
                    
					if((p_f%cells(dim,i,j,k) < 0.0).or.(rho_f%cells(dim,i,j,k) < 0.0)) then
                    !if ((bc%bc_markers(i,j,k) == 0).and.(bc%bc_markers(i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3)) == 0)) then
						!print *, 'Flow EOS exception: ', dim, i,j,k
						!do dim2 = 1, dimensions
						!	print *, 'Coordinates', dim2, mesh%mesh(dim2,i,j,k)
						!end do						
						!print *, gamma_f%cells(dim,i,j,k)
						!print *, p_f%cells(dim,i,j,k)
						!print *, rho_f%cells(dim,i,j,k)
						!print *, v_s_f%cells(dim,i,j,k)
                        !gamma_f%cells(dim,i,j,k) = abs(gamma_f%cells(dim,i,j,k))
						!p_f%cells(dim,i,j,k) = abs(p_f%cells(dim,i,j,k))
						!rho_f%cells(dim,i,j,k) = abs(rho_f%cells(dim,i,j,k))
                        
						!stop
                        
                        rho_l	= rho%cells(i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3))
                        rho_r	= rho%cells(i,j,k)
                        p_l		= p%cells(i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3))
                        p_r		= p%cells(i,j,k)
                        gamma_l = gamma%cells(i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3))
                        gamma_r	= gamma%cells(i,j,k)
                        v_l		= v%pr(dim)%cells(i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3))
                        v_r		= v%pr(dim)%cells(i,j,k)
                        
                        print *, p_f%cells(dim,i,j,k), rho_f%cells(dim,i,j,k), v_f%pr(dim)%cells(dim,i,j,k)
                        
                        call riemann%set_parameters(	rho_l	, rho_r		,	&
														p_l		, p_r		,   &
														gamma_l	, gamma_r	,	&
														v_l		, v_r)
                       
                        call riemann%solve()
                        
                        rho_f%cells(dim,i,j,k)			= riemann%get_density()
                        p_f%cells(dim,i,j,k)			= riemann%get_pressure()
                        v_f%pr(dim)%cells(dim,i,j,k)	= riemann%get_velocity()
                        
                        call riemann%clear()
                        
                        print *, p_f%cells(dim,i,j,k), rho_f%cells(dim,i,j,k), v_f%pr(dim)%cells(dim,i,j,k)
                        
                        print *, 'Riemann activation count: ', riemann%get_activation_count()
!                        pause
					end if                   
                    

					T_old	= T_f%cells(dim,i,j,k)
					T_f%cells(dim,i,j,k)	= p_f%cells(dim,i,j,k) / rho_f%cells(dim,i,j,k) / r_gase_J * average_molar_mass
					t_initial				= T_f%cells(dim,i,j,k)
				
					cp = this%thermo%thermo_ptr%calculate_mixture_cp(t_initial, concs)
					cv = cp - r_gase_J 

					gamma_f%cells(dim,i,j,k)	= cp / cv
					
					e_i_f%cells(dim,i,j,k)	= this%thermo%thermo_ptr%calculate_mixture_energy(T_f%cells(dim,i,j,k), concs) - this%thermo%thermo_ptr%calculate_mixture_enthalpy(298.15_dkind, concs)
					
				!# Experimental iterative procedure. Due to discontnuity in T=1000K for sensible enthalpy this approach provides poor results		
					!h_s					= this%thermo%thermo_ptr%calculate_mixture_enthalpy(T_f%cells(dim,i,j,k), concs)
					!e_i_f%cells(dim,i,j,k)	= (h_s - r_gase_J*T_f%cells(dim,i,j,k)) / average_molar_mass !p%cells(i,j,k)/rho%cells(i,j,k)
					
				!#  Old e = cv*T equation of state.	
				!	e_i_f%cells(dim,i,j,k)	= p_f%cells(dim,i,j,k) / rho_f%cells(dim,i,j,k) / (gamma_f%cells(dim,i,j,k) - 1.0_dkind) 
					
				!# New de = cv*dT equation of state		
				!	e_i_f%cells(dim,i,j,k)	= (e_i_f_old%cells(dim,i,j,k) + (T_f%cells(dim,i,j,k) - T_old)*c_v_f_old%cells(dim,i,j,k)) /average_molar_mass 

					v_s_f%cells(dim,i,j,k)	= sqrt(gamma_f%cells(dim,i,j,k)*p_f%cells(dim,i,j,k)/rho_f%cells(dim,i,j,k))
                    
					E_f_f%cells(dim,i,j,k)	= e_i_f%cells(dim,i,j,k) / average_molar_mass

					do dim1 = 1, dimensions
						E_f_f%cells(dim,i,j,k) = E_f_f%cells(dim,i,j,k) + 0.5_dkind*v_f%pr(dim1)%cells(dim,i,j,k)*v_f%pr(dim1)%cells(dim,i,j,k)
					end do
				end if	
			end do
			end do
			end do

			!$omp end do nowait
		end do
			
		!$omp end parallel

		continue

		end associate

		continue
	end subroutine	

	subroutine check_conservation_laws(this)
		class(table_approximated_real_gas) ,intent(in) :: this

		real(dkind)	:: total_mass_error, total_energy_error
		real(dkind)	,dimension(:), allocatable	:: total_species_error
		
		integer		:: species_number
		integer	:: dimensions
		integer	,dimension(3,2)	:: inner_loop
		integer	:: specie_number
		integer	:: i,j,k

		species_number = this%chem%chem_ptr%species_number
		
		allocate(total_species_error(species_number))
		
		total_energy_error	=  this%total_energy
		total_mass_error	=  this%total_mass
		do specie_number = 1,species_number
			total_species_error(specie_number)	=	this%total_species(specie_number) 
		end do
				
		associate(  rho     => this%rho%s_ptr   , &
					E_f		=> this%E_f%s_ptr	, &
					Y		=> this%Y%v_ptr)	

			dimensions		= this%domain%get_domain_dimensions()
			inner_loop		= this%domain%get_local_inner_cells_bounds()				
			
			do k = inner_loop(3,1),inner_loop(3,2)
			do j = inner_loop(2,1),inner_loop(2,2)
			do i = inner_loop(1,1),inner_loop(1,2)					
					
				if(this%boundary%bc_ptr%bc_markers(i,j,k) == 0) then

					total_energy_error	= total_energy_error	- E_f%cells(i,j,k)
					total_mass_error	= total_mass_error		- rho%cells(i,j,k)
					
					do specie_number = 1,species_number
						total_species_error(specie_number)	=	total_species_error(specie_number) - Y%pr(specie_number)%cells(i,j,k) * rho%cells(i,j,k)
					end do
				
				end if

			end do
			end do
			end do	
			
			print *, total_energy_error / this%total_energy
			print *, total_mass_error /  this%total_mass
		!	print *, total_species_error
			
			continue
		end associate
	end subroutine
	
	subroutine apply_boundary_conditions_for_initial_conditions(this)
		class(table_approximated_real_gas) ,intent(inout) :: this

		real(dkind)                     :: cp, cv

		real(dkind) ,dimension(this%chem%chem_ptr%species_number)    :: concs
		real(dkind)	:: average_molar_mass, mol_mix_conc

		integer	:: bound_number
		integer	:: species_number, specie_index

		integer	:: dimensions
		integer	,dimension(3,2)	:: utter_loop, inner_loop
		character(len=20)		:: boundary_type_name
		real(dkind)				:: wall_temperature
		real(dkind)				:: farfield_density, farfield_pressure, farfield_temperature, farfield_velocity     !		farfield_gamma, farfield_v_s, farfield_mol_mix_conc, farfield_E_f, farfield_e_i, farfield_v
		real(dkind)				:: farfield_E_f, farfield_e_i, farfield_gamma, farfield_v_s
		real(dkind)			,dimension(:)	,allocatable	:: farfield_concentrations
		character(len=10)	,dimension(:)	,allocatable	:: farfield_species_names

		integer	:: specie_number
		integer	:: sign
		integer	:: i,j,k,dim,plus

		species_number = this%chem%chem_ptr%species_number

		associate(		T               => this%T%s_ptr				, &
						p               => this%p%s_ptr				, &
						rho             => this%rho%s_ptr			, &
						e_i             => this%e_i%s_ptr			, &
						E_f             => this%E_f%s_ptr			, &
						h_s				=> this%h_s%s_ptr			, &
						v_s             => this%v_s%s_ptr			, &
						Y               => this%Y%v_ptr				, &
						v               => this%v%v_ptr)

		dimensions		= this%domain%get_domain_dimensions()
		utter_loop		= this%domain%get_local_utter_cells_bounds()
		inner_loop		= this%domain%get_local_inner_cells_bounds()

		do k = inner_loop(3,1),inner_loop(3,2)
		do j = inner_loop(2,1),inner_loop(2,2)
		do i = inner_loop(1,1),inner_loop(1,2)
			if(this%boundary%bc_ptr%bc_markers(i,j,k) == 0) then

				do dim = 1,dimensions
					if(	((I_m(dim,1)*i + I_m(dim,2)*j + I_m(dim,3)*k) /= utter_loop(dim,1)).and.	&
						((I_m(dim,1)*i + I_m(dim,2)*j + I_m(dim,3)*k) /= utter_loop(dim,2)))  then
						do plus	= 1, 2
							sign		= (-1)**plus
							bound_number	= this%boundary%bc_ptr%bc_markers(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3))
							if(bound_number /= 0) then

								v_s%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3)) =  v_s%cells(i,j,k)
								rho%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3)) =  rho%cells(i,j,k)
								E_f%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3)) =  E_f%cells(i,j,k)

								boundary_type_name = this%boundary%bc_ptr%boundary_types(bound_number)%get_type_name()
								select case(boundary_type_name)
									case('wall')
										if(this%boundary%bc_ptr%boundary_types(bound_number)%is_conductive()) then
											wall_temperature = this%boundary%bc_ptr%boundary_types(bound_number)%get_wall_temperature()
											T%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3)) = wall_temperature
										else
											T%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3)) = T%cells(i,j,k)
										end if

										do specie_number = 1, species_number
											Y%pr(specie_number)%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3))	=	Y%pr(specie_number)%cells(i,j,k)
										end do
										
									case('outlet','inlet')
								
										farfield_pressure		= this%boundary%bc_ptr%boundary_types(bound_number)%get_farfield_pressure()
										farfield_temperature	= this%boundary%bc_ptr%boundary_types(bound_number)%get_farfield_temperature()
									!	farfield_velocity		= this%boundary%bc_ptr%boundary_types(bound_number)%get_farfield_velocity()

										call this%boundary%bc_ptr%boundary_types(bound_number)%get_farfield_species_names(farfield_species_names)
										call this%boundary%bc_ptr%boundary_types(bound_number)%get_farfield_concentrations(farfield_concentrations)
								
										concs = 0.0_dkind
										do specie_number = 1, size(farfield_species_names)
											specie_index		= this%chem%chem_ptr%get_chemical_specie_index(farfield_species_names(specie_number))
											concs(specie_index) = farfield_concentrations(specie_number)
										end do

										call this%thermo%thermo_ptr%change_cell_units_mole_to_dimless(concs)									

										do specie_number = 1, size(concs)
											Y%pr(specie_number)%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3))	= 0.0_dkind
										end do
										
										do specie_number = 1, size(farfield_species_names)
											specie_index		= this%chem%chem_ptr%get_chemical_specie_index(farfield_species_names(specie_number))
											Y%pr(specie_index)%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3))	=	concs(specie_index)
										end do

										average_molar_mass = 0.0_dkind
										do specie_number = 1,size(farfield_species_names)
											specie_index		= this%chem%chem_ptr%get_chemical_specie_index(farfield_species_names(specie_number))
											average_molar_mass = average_molar_mass + concs(specie_index) / this%thermo%thermo_ptr%molar_masses(specie_index)
										end do
				
										average_molar_mass	=  1.0_dkind / average_molar_mass
										mol_mix_conc		=  average_molar_mass
				
										farfield_density			= farfield_pressure / (farfield_temperature * r_gase_J) * mol_mix_conc
										call this%boundary%bc_ptr%boundary_types(bound_number)%set_farfield_density(farfield_density)

										do specie_number = 1,size(farfield_species_names)
											specie_index			= this%chem%chem_ptr%get_chemical_specie_index(farfield_species_names(specie_number))
											concs(specie_index)		= concs(specie_index) * mol_mix_conc / this%thermo%thermo_ptr%molar_masses(specie_index)
										end do
									
										cp = this%thermo%thermo_ptr%calculate_mixture_cp(farfield_temperature, concs)
										cv = cp - r_gase_J

										farfield_gamma	= cp / cv	
								
										!h_s%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3))		=	this%thermo%thermo_ptr%calculate_mixture_enthalpy(farfield_temperature, concs)/ mol_mix_conc
										farfield_e_i	= (this%thermo%thermo_ptr%calculate_mixture_energy(farfield_temperature, concs) - this%thermo%thermo_ptr%calculate_mixture_enthalpy(298.15_dkind, concs)) / mol_mix_conc
										
                                        farfield_velocity	=  sqrt(abs((p%cells(i,j,k) - farfield_pressure)*(rho%cells(i,j,k) - farfield_density)/farfield_density/rho%cells(i,j,k)))
                                        
										!farfield_e_i	= farfield_pressure / farfield_density / (farfield_gamma - 1.0_dkind) 
										farfield_E_f	= farfield_e_i + 0.5_dkind * ( farfield_velocity * farfield_velocity )
										farfield_v_s	= sqrt(farfield_gamma*farfield_pressure/farfield_density)	
				
										p%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3))			=	farfield_pressure 
										T%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3))			=	farfield_temperature
										v_s%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3))		=	farfield_v_s 
										rho%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3))		=	farfield_density
										v%pr(dim)%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3))	=	farfield_velocity

										E_f%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3))		=	farfield_E_f 
										call this%boundary%bc_ptr%boundary_types(bound_number)%set_farfield_energy(farfield_E_f)
                                        
										h_s%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3))		=	(this%thermo%thermo_ptr%calculate_mixture_enthalpy(farfield_temperature, concs) - this%thermo%thermo_ptr%calculate_mixture_enthalpy(298.15_dkind, concs)) / mol_mix_conc
										!h_s%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3))		=	cp*farfield_temperature/mol_mix_conc
										
								end select
							end if
						end do
					end if
				end do
			end if
		end do
		end do
		end do

		end associate		
	end subroutine
	
end module
