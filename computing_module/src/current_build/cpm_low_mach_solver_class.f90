module cpm_low_mach_solver_class

	use kind_parameters
	use global_data
	use data_manager_class
	use data_io_class	
	use computational_domain_class
	use computational_mesh_class
	use boundary_conditions_class
	use field_pointers
	use table_approximated_real_gas_class
	use chemical_kinetics_solver_class
	use chemical_properties_class
	
	use viscosity_solver_class
	use coarse_particles_method
	use fourier_heat_transfer_solver_class
	use fickean_diffusion_solver_class

	use mpi_communications_class

	use solver_options_class	
	
	implicit none

#ifdef OMP	
	include "omp_lib.h"
#endif
	
	private
	public	:: cpm_low_mach_solver, cpm_low_mach_solver_c

	type(field_scalar_cons)	,target	:: E_f_int
	type(field_vector_cons)	,target	:: v_int, Y_int

	type cpm_low_mach_solver
		logical			:: diffusion_flag, viscosity_flag, heat_trans_flag, reactive_flag, sources_flag, hydrodynamics_flag, CFL_condition_flag
		real(dkind)		:: courant_fraction
		real(dkind)		:: time, time_step, initial_time_step
		type(viscosity_solver)				:: visc_solver
		type(coarse_particles)				:: gas_dynamics_solver
		type(heat_transfer_solver)			:: heat_trans_solver
		type(diffusion_solver)				:: diff_solver
		type(chemical_kinetics_solver)		:: chem_kin_solver
		type(table_approximated_real_gas)	:: state_eq

		type(computational_domain)					:: domain
		type(mpi_communications)					:: mpi_support
		type(chemical_properties_pointer)			:: chem
		type(computational_mesh_pointer)			:: mesh
		type(boundary_conditions_pointer)			:: boundary

		type(field_scalar_cons_pointer)	:: rho	, T				, p				, v_s			, mol_mix_conc
		type(field_scalar_cons_pointer)	:: E_f	, E_f_prod_chem	, E_f_prod_heat	, E_f_prod_gd	, E_f_prod_visc	, E_f_prod_diff , E_f_int   ,e_i
		type(field_scalar_cons_pointer)	:: p_stat, p_dyn, div_v
		
		type(field_vector_cons_pointer)	:: v	, v_prod_gd		, v_prod_visc	, v_prod_source	, v_int
		type(field_vector_cons_pointer)	:: Y	, Y_prod_diff	, Y_prod_chem	, Y_int
      
	contains
		procedure	,private	:: apply_boundary_conditions_main
		procedure	,private	:: apply_boundary_conditions_interm_v
		procedure	,private	:: apply_boundary_conditions_interm_E_Y		
		procedure	,private	:: calculate_interm_v
		procedure	,private	:: calculate_interm_E_Y
 
		procedure				:: solve_problem
		procedure				:: calculate_time_step
		procedure				:: get_time_step
		procedure				:: get_time
	end type

	interface	cpm_low_mach_solver_c
		module procedure	constructor
	end interface

contains

	type(cpm_low_mach_solver)	function constructor(manager,problem_data_io, problem_solver_options)
		type(data_manager)						,intent(inout)	:: manager
		type(data_io)							,intent(inout)	:: problem_data_io
		type(solver_options)					,intent(in)		:: problem_solver_options

		real(dkind)	:: calculation_time	

		type(field_scalar_cons_pointer)	:: scal_ptr
		type(field_vector_cons_pointer)	:: vect_ptr
		type(field_tensor_cons_pointer)	:: tens_ptr
        
		constructor%diffusion_flag		= problem_solver_options%get_molecular_diffusion_flag()
		constructor%viscosity_flag		= problem_solver_options%get_viscosity_flag()
		constructor%heat_trans_flag		= problem_solver_options%get_heat_transfer_flag()
		constructor%reactive_flag		= problem_solver_options%get_chemical_reaction_flag()
		constructor%hydrodynamics_flag	= problem_solver_options%get_hydrodynamics_flag()
		constructor%courant_fraction	= problem_solver_options%get_CFL_condition_coefficient()
		constructor%CFL_condition_flag	= problem_solver_options%get_CFL_condition_flag()
		constructor%sources_flag		= .false.
		
		constructor%domain				= manager%domain
		constructor%mpi_support			= manager%mpi_communications
		constructor%chem%chem_ptr		=> manager%chemistry%chem_ptr
		constructor%boundary%bc_ptr	=> manager%boundary_conditions_pointer%bc_ptr
		constructor%mesh%mesh_ptr		=> manager%computational_mesh_pointer%mesh_ptr

		call manager%get_cons_field_pointer_by_name(scal_ptr,vect_ptr,tens_ptr,'density')
		constructor%rho%s_ptr					=> scal_ptr%s_ptr
		call manager%get_cons_field_pointer_by_name(scal_ptr,vect_ptr,tens_ptr,'temperature')
		constructor%T%s_ptr					=> scal_ptr%s_ptr
		call manager%get_cons_field_pointer_by_name(scal_ptr,vect_ptr,tens_ptr,'pressure')
		constructor%p%s_ptr					=> scal_ptr%s_ptr
		call manager%get_cons_field_pointer_by_name(scal_ptr,vect_ptr,tens_ptr,'full_energy')
		constructor%E_f%s_ptr					=> scal_ptr%s_ptr
		call manager%get_cons_field_pointer_by_name(scal_ptr,vect_ptr,tens_ptr,'internal_energy')
		constructor%e_i%s_ptr					=> scal_ptr%s_ptr        
		call manager%get_cons_field_pointer_by_name(scal_ptr,vect_ptr,tens_ptr,'velocity')
		constructor%v%v_ptr						=> vect_ptr%v_ptr
		call manager%get_cons_field_pointer_by_name(scal_ptr,vect_ptr,tens_ptr,'specie_molar_concentration')
		constructor%Y%v_ptr						=> vect_ptr%v_ptr
		call manager%get_cons_field_pointer_by_name(scal_ptr,vect_ptr,tens_ptr,'mixture_molar_concentration')
		constructor%mol_mix_conc%s_ptr		=> scal_ptr%s_ptr	

		call manager%create_scalar_field(E_f_int	,'full_energy_interm'					,'E_f_int')
		constructor%E_f_int%s_ptr				=> E_f_int
		
		call manager%create_vector_field(v_int		,'velocity_interm'						,'v_int'	,'spatial')
		constructor%v_int%v_ptr					=> v_int
		call manager%create_vector_field(Y_int		,'specie_molar_concentration_interm'	,'Y_int'	,'chemical')
		constructor%Y_int%v_ptr					=> Y_int
		
		constructor%gas_dynamics_solver	= coarse_particles_c(manager,low_mach_flag = .true.)
		
		
		call manager%get_cons_field_pointer_by_name(scal_ptr,vect_ptr,tens_ptr,'velocity_production_gas_dynamics')
		constructor%v_prod_gd%v_ptr				=> vect_ptr%v_ptr
		call manager%get_cons_field_pointer_by_name(scal_ptr,vect_ptr,tens_ptr,'energy_production_gas_dynamics')
		constructor%E_f_prod_gd%s_ptr			=> scal_ptr%s_ptr

		call manager%get_cons_field_pointer_by_name(scal_ptr,vect_ptr,tens_ptr,'pressure_dynamic')
		constructor%p_dyn%s_ptr				=> scal_ptr%s_ptr		
		call manager%get_cons_field_pointer_by_name(scal_ptr,vect_ptr,tens_ptr,'pressure_static')
		constructor%p_stat%s_ptr			=> scal_ptr%s_ptr
		call manager%get_cons_field_pointer_by_name(scal_ptr,vect_ptr,tens_ptr,'velocity_divergence')
		constructor%div_v%s_ptr				=> scal_ptr%s_ptr		
		
		if(constructor%viscosity_flag) then
			constructor%visc_solver			= viscosity_solver_c(manager)
			call manager%get_cons_field_pointer_by_name(scal_ptr,vect_ptr,tens_ptr,'energy_production_viscosity')
			constructor%E_f_prod_visc%s_ptr			=> scal_ptr%s_ptr
			call manager%get_cons_field_pointer_by_name(scal_ptr,vect_ptr,tens_ptr,'velocity_production_viscosity')
			constructor%v_prod_visc%v_ptr			=> vect_ptr%v_ptr
		end if

		if (constructor%heat_trans_flag) then
			constructor%heat_trans_solver	= heat_transfer_solver_c(manager)
			call manager%get_cons_field_pointer_by_name(scal_ptr,vect_ptr,tens_ptr,'energy_production_heat_transfer')
			constructor%E_f_prod_heat%s_ptr			=> scal_ptr%s_ptr
		end if

		if (constructor%diffusion_flag) then
			constructor%diff_solver			= diffusion_solver_c(manager)
			call manager%get_cons_field_pointer_by_name(scal_ptr,vect_ptr,tens_ptr,'energy_production_diffusion')
			constructor%E_f_prod_diff%s_ptr			=> scal_ptr%s_ptr			
			call manager%get_cons_field_pointer_by_name(scal_ptr,vect_ptr,tens_ptr,'specie_production_diffusion')
			constructor%Y_prod_diff%v_ptr			=> vect_ptr%v_ptr
		end if

		if (constructor%reactive_flag) then
			constructor%chem_kin_solver		= chemical_kinetics_solver_c(manager)
			call manager%get_cons_field_pointer_by_name(scal_ptr,vect_ptr,tens_ptr,'energy_production_chemistry')
			constructor%E_f_prod_chem%s_ptr			=> scal_ptr%s_ptr
			call manager%get_cons_field_pointer_by_name(scal_ptr,vect_ptr,tens_ptr,'specie_production_chemistry')
			constructor%Y_prod_chem%v_ptr			=> vect_ptr%v_ptr
        end if

		!constructor%sources	= sources_c(manager)
		!call manager%get_field_pointer(scal_ptr,vect_ptr,tens_ptr,'velocity_production_sources')
		!constructor%v_prod_source%v_ptr			=> vect_ptr%v_ptr

 		constructor%state_eq	=	table_approximated_real_gas_c(manager)
		call manager%get_cons_field_pointer_by_name(scal_ptr,vect_ptr,tens_ptr,'velocity_of_sound')
		constructor%v_s%s_ptr			=> scal_ptr%s_ptr										

		problem_data_io		= data_io_c(manager,calculation_time)
		
		if(problem_data_io%get_load_counter() /= 0) then
			call problem_data_io%add_io_scalar_cons_field(constructor%p_stat)
			call problem_data_io%add_io_scalar_cons_field(constructor%E_f)
			call problem_data_io%add_io_scalar_cons_field(constructor%e_i)
			call problem_data_io%add_io_scalar_cons_field(constructor%mol_mix_conc)
			call problem_data_io%add_io_scalar_cons_field(constructor%div_v)
		end if
		
		call problem_data_io%input_all_data()
			
		if(problem_data_io%get_load_counter() == 1) then
			call problem_data_io%add_io_scalar_cons_field(constructor%p_stat)
			call problem_data_io%add_io_scalar_cons_field(constructor%E_f)
			call problem_data_io%add_io_scalar_cons_field(constructor%e_i)
			call problem_data_io%add_io_scalar_cons_field(constructor%mol_mix_conc)
			call problem_data_io%add_io_scalar_cons_field(constructor%div_v)
			
			constructor%p_dyn%s_ptr%cells	= 0.0_dkind      
			constructor%p_stat%s_ptr%cells	= constructor%p%s_ptr%cells	
		end if				

		if(problem_data_io%get_load_counter() == 1) then
			call constructor%state_eq%apply_state_equation_for_initial_conditions()
		!else
		!	call constructor%state_eq%apply_state_equation_low_mach()
		end if

		call constructor%mpi_support%exchange_conservative_scalar_field(constructor%p%s_ptr)
		call constructor%mpi_support%exchange_conservative_scalar_field(constructor%p_stat%s_ptr)
		call constructor%mpi_support%exchange_conservative_scalar_field(constructor%rho%s_ptr)
		call constructor%mpi_support%exchange_conservative_scalar_field(constructor%T%s_ptr)

		call constructor%mpi_support%exchange_conservative_vector_field(constructor%Y%v_ptr)
		call constructor%mpi_support%exchange_conservative_vector_field(constructor%v%v_ptr)

		call constructor%mpi_support%exchange_boundary_conditions_markers(constructor%boundary%bc_ptr)
		call constructor%mpi_support%exchange_mesh(constructor%mesh%mesh_ptr)

		constructor%p_dyn%s_ptr%cells				= constructor%p%s_ptr%cells - constructor%p_stat%s_ptr%cells
	!	constructor%p_dyn%s_ptr%cells(90:111,1,1)	= 10.0_dkind  	
		
		constructor%time				=	calculation_time
		constructor%time_step			=	problem_solver_options%get_initial_time_step()
		constructor%initial_time_step	=	problem_solver_options%get_initial_time_step()
		
	end function

	subroutine solve_problem(this)
		class(cpm_low_mach_solver)	,intent(inout)	:: this

		integer	:: specie

		this%time = this%time + this%time_step		
		
		call this%apply_boundary_conditions_main()

		if (this%viscosity_flag)	call this%visc_solver%solve_viscosity(this%time_step)
		if (this%heat_trans_flag)	call this%heat_trans_solver%solve_heat_transfer(this%time_step)
		if (this%diffusion_flag)	call this%diff_solver%solve_diffusion(this%time_step)
		if (this%reactive_flag)		call this%chem_kin_solver%solve_chemical_kinetics(this%time_step)
		!call this%sources%calculate_sources(this%time,this%time_step)
        
		call this%gas_dynamics_solver%euler_step_v(this%time_step)
		call this%calculate_interm_v(this%time_step)
		call this%apply_boundary_conditions_interm_v()

		call this%calculate_interm_E_Y(this%time_step)
		call this%apply_boundary_conditions_interm_E_Y()
		
		call this%gas_dynamics_solver%lagrange_step(this%time_step)
		call this%gas_dynamics_solver%final_step(this%time_step)

		call this%state_eq%apply_state_equation_low_mach()
        
        call this%gas_dynamics_solver%calculate_p_dyn(this%time_step)
		!call this%state_eq%check_conservation_laws()

		if (this%CFL_condition_flag) then
			call this%calculate_time_step()
		end if		
		
	end subroutine

	subroutine calculate_interm_v(this,time_step)
		class(cpm_low_mach_solver)	,intent(inout)	:: this
		real(dkind)					,intent(in)		:: time_step
		
		integer	:: dimensions
		integer	,dimension(3,2)	:: cons_inner_loop

		integer	:: bound_number
		integer :: i,j,k,dim,specie_number
		
		dimensions		= this%domain%get_domain_dimensions()
		
		cons_inner_loop	= this%domain%get_local_inner_cells_bounds()
				
		associate (	rho				=> this%rho%s_ptr			, &
					v				=> this%v%v_ptr				, &
					v_int			=> this%v_int%v_ptr			, &
					v_prod_gd		=> this%v_prod_gd%v_ptr		, &
					v_prod_visc		=> this%v_prod_visc%v_ptr	, &
					v_prod_source	=> this%v_prod_source%v_ptr , &
					bc				=> this%boundary%bc_ptr)
					
		!$omp parallel default(none)  private(i,j,k,dim) , &
		!$omp& firstprivate(this)	,&
		!$omp& shared(rho,v,v_int,v_prod_gd,v_prod_visc,v_prod_source,bc,cons_inner_loop,dimensions,time_step)
		!$omp do collapse(3) schedule(guided)					
					
			do k = cons_inner_loop(3,1),cons_inner_loop(3,2)
			do j = cons_inner_loop(2,1),cons_inner_loop(2,2)
			do i = cons_inner_loop(1,1),cons_inner_loop(1,2)
				if(bc%bc_markers(i,j,k) == 0) then					
					do dim = 1, dimensions
						v_int%pr(dim)%cells(i,j,k) = v%pr(dim)%cells(i,j,k) + v_prod_gd%pr(dim)%cells(i,j,k)! + v_prod_source%pr(dim)%cells(i,j,k)
						if (this%viscosity_flag)	v_int%pr(dim)%cells(i,j,k) = v_int%pr(dim)%cells(i,j,k) + v_prod_visc%pr(dim)%cells(i,j,k) / rho%cells(i,j,k) * time_step
					end do
				end if
			end do
			end do
			end do
		!$omp end do nowait
		!$omp end parallel
	
			
		end associate

	end subroutine
	
	
	subroutine calculate_interm_E_Y(this,time_step)
		class(cpm_low_mach_solver)	,intent(inout)	:: this
		real(dkind)					,intent(in)		:: time_step
		
		integer	:: dimensions, species_number
		integer	,dimension(3,2)	:: cons_inner_loop

		integer	:: bound_number
		integer :: i,j,k,dim,spec

		real(dkind)		:: spec_summ
		real(dkind)		:: energy_source = 1.0_dkind
		real(dkind)	,save	:: energy_output = 0.0_dkind
				
		dimensions		= this%domain%get_domain_dimensions()
		species_number	= this%chem%chem_ptr%species_number
				
		cons_inner_loop	= this%domain%get_local_inner_cells_bounds()		
		
		associate (	rho				=> this%rho%s_ptr			, &
					E_f				=> this%E_f%s_ptr			, &
					E_f_int 		=> this%E_f_int%s_ptr		, &
					E_f_prod_chem 	=> this%E_f_prod_chem%s_ptr	, &
					E_f_prod_heat	=> this%E_f_prod_heat%s_ptr	, &
					E_f_prod_gd 	=> this%E_f_prod_gd%s_ptr	, &
					E_f_prod_visc	=> this%E_f_prod_visc%s_ptr	, &
					E_f_prod_diff	=> this%E_f_prod_diff%s_ptr	, &
					Y				=> this%Y%v_ptr				, &
					Y_int			=> this%Y_int%v_ptr			, &
					Y_prod_diff		=> this%Y_prod_diff%v_ptr	, &
					Y_prod_chem		=> this%Y_prod_chem%v_ptr	, &
					T				=> this%T%s_ptr				, &
					bc				=> this%boundary%bc_ptr)

		!$omp parallel default(none)  private(i,j,k,dim,spec,spec_summ,energy_output) , &
		!$omp& firstprivate(this)	,&
		!$omp& shared(rho,E_f,E_f_int,E_f_prod_chem,E_f_prod_heat,E_f_prod_gd,E_f_prod_visc,E_f_prod_diff,Y,Y_int,Y_prod_diff,Y_prod_chem,dimensions,species_number,bc,cons_inner_loop,energy_source,time_step)
		!$omp do collapse(3) schedule(guided)					
					
			do k = cons_inner_loop(3,1),cons_inner_loop(3,2)
			do j = cons_inner_loop(2,1),cons_inner_loop(2,2)
			do i = cons_inner_loop(1,1),cons_inner_loop(1,2)
			
				if(bc%bc_markers(i,j,k) == 0) then					
					E_f_int%cells(i,j,k) = E_f%cells(i,j,k) 
                    
					if (this%viscosity_flag)	E_f_int%cells(i,j,k) = E_f_int%cells(i,j,k) + E_f_prod_visc%cells(i,j,k) * time_step / rho%cells(i,j,k) 
					if (this%heat_trans_flag)	E_f_int%cells(i,j,k) = E_f_int%cells(i,j,k) + E_f_prod_heat%cells(i,j,k) * time_step / rho%cells(i,j,k)
					if (this%diffusion_flag)	E_f_int%cells(i,j,k) = E_f_int%cells(i,j,k) + E_f_prod_diff%cells(i,j,k) * time_step / rho%cells(i,j,k)
					if (this%reactive_flag)		E_f_int%cells(i,j,k) = E_f_int%cells(i,j,k) + E_f_prod_chem%cells(i,j,k) * time_step / rho%cells(i,j,k)

					! ************************* Energy release ******************
					!if(this%time <= 1.0E-07_dkind) then
						if (((i == 25).or.(i == 26)).and.((j == 1).or.(j == 1)).and.(k == 1)) then
						!	if (this%sources_flag)	then
								E_f_int%cells(i,j,k)	= E_f_int%cells(i,j,k)	+ energy_source * this%time_step * 1.0E06_dkind 
								energy_output			= energy_output			+ energy_source * this%time_step * 1.0E06_dkind *  rho%cells(i,1,1)
						!	end if
						end if
					!end if
					! ***********************************************************	

					spec_summ = 0.0_dkind
					do spec = 1, species_number
						Y_int%pr(spec)%cells(i,j,k) = Y%pr(spec)%cells(i,j,k)
					
						if (this%diffusion_flag)	Y_int%pr(spec)%cells(i,j,k) = Y_int%pr(spec)%cells(i,j,k) + Y_prod_diff%pr(spec)%cells(i,j,k)  * time_step / rho%cells(i,j,k)
						if (this%reactive_flag)		Y_int%pr(spec)%cells(i,j,k) = Y_int%pr(spec)%cells(i,j,k) + Y_prod_chem%pr(spec)%cells(i,j,k)  * time_step / rho%cells(i,j,k)
						
						spec_summ = spec_summ + max(Y_int%pr(spec)%cells(i,j,k), 0.0_dkind)
						
					end do
					
					do spec = 1,species_number
						Y_int%pr(spec)%cells(i,j,k) = max(Y_int%pr(spec)%cells(i,j,k), 0.0_dkind) / spec_summ
					end do					
				
				end if
		
			end do
			end do
			end do
		!$omp end do nowait
		!$omp end parallel	
			
			
		if((this%time <= 1.0E-06_dkind).and.(this%sources_flag)) then		
			print *, energy_output
		end if				
			
		end associate

    end subroutine

    
	subroutine apply_boundary_conditions_main(this)

		class(cpm_low_mach_solver)		,intent(inout)		:: this

		character(len=20)		:: boundary_type_name
		real(dkind)				:: farfield_density, farfield_pressure, wall_temperature
		
		integer					:: dimensions, species_number
		integer	,dimension(3,2)	:: cons_inner_loop

		integer	:: sign, bound_number
		integer :: i,j,k,plus,dim,dim1,spec
						
		dimensions		= this%domain%get_domain_dimensions()
		species_number	= this%chem%chem_ptr%species_number
				
		cons_inner_loop	= this%domain%get_local_inner_cells_bounds()			
		
		associate(  T				=> this%T%s_ptr				, &
					mol_mix_conc	=> this%mol_mix_conc%s_ptr	, &
					p				=> this%p%s_ptr				, &
					p_dyn			=> this%p_dyn%s_ptr			, &
					p_stat			=> this%p_stat%s_ptr		, &
					rho				=> this%rho%s_ptr			, &
					v				=> this%v%v_ptr				, &
					Y				=> this%Y%v_ptr				, &
					bc				=> this%boundary%bc_ptr		, &
					mesh			=> this%mesh%mesh_ptr)

		!$omp parallel default(none)  private(i,j,k,plus,dim,dim1,spec,sign,bound_number,boundary_type_name,farfield_pressure,farfield_density,wall_temperature) , &
		!$omp& firstprivate(this)	,&
		!$omp& shared(T,p,p_dyn,p_stat,rho,v,mol_mix_conc,Y,mesh,bc,cons_inner_loop,dimensions,species_number)
		!$omp do collapse(3) schedule(guided)					
					
			do k = cons_inner_loop(3,1),cons_inner_loop(3,2)
			do j = cons_inner_loop(2,1),cons_inner_loop(2,2)
			do i = cons_inner_loop(1,1),cons_inner_loop(1,2)
				if(bc%bc_markers(i,j,k) == 0) then
					do dim = 1,dimensions
						do plus = 1,2
							sign			= (-1)**plus
							bound_number	= bc%bc_markers(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3))
							if( bound_number /= 0 ) then

	
								boundary_type_name = bc%boundary_types(bound_number)%get_type_name()
								select case(boundary_type_name)
									case('wall')
									
										p%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3))				= p%cells(i,j,k)
										p_dyn%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3))			= p_dyn%cells(i,j,k)
										p_stat%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3))			= p_stat%cells(i,j,k)
										rho%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3))			= rho%cells(i,j,k)
										T%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3))				= T%cells(i,j,k)
										mol_mix_conc%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3))	= mol_mix_conc%cells(i,j,k)
								
										do spec = 1, species_number
											Y%pr(spec)%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3))	=	Y%pr(spec)%cells(i,j,k)
										end do

										do dim1 = 1, dimensions
											if(dim1 == dim) then
												v%pr(dim1)%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3)) = - v%pr(dim1)%cells(i,j,k)
											else
												v%pr(dim1)%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3)) = v%pr(dim1)%cells(i,j,k)
											end if
										end do												
									
										if(bc%boundary_types(bound_number)%is_conductive()) then 
											wall_temperature = bc%boundary_types(bound_number)%get_wall_temperature()
											T%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3)) = wall_temperature
										end if
										if(.not.bc%boundary_types(bound_number)%is_slip()) then
											do dim1 = 1, dimensions
												if (dim1 /= dim) then	
													v%pr(dim1)%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3)) = v%pr(dim1)%cells(i-sign*I_m(dim,1),j-sign*I_m(dim,2),k-sign*I_m(dim,3)) - 6.0_dkind * v%pr(dim1)%cells(i,j,k)
												end if
											end do
										end if
									case ('outlet')
									
										p%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3))				= p%cells(i,j,k)
										p_dyn%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3))			= p_dyn%cells(i,j,k)
										p_stat%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3))			= p_stat%cells(i,j,k)
										rho%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3))			= rho%cells(i,j,k)
										T%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3))				= T%cells(i,j,k)
										mol_mix_conc%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3))	= mol_mix_conc%cells(i,j,k)
								
										do spec = 1, species_number
											Y%pr(spec)%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3))	=	Y%pr(spec)%cells(i,j,k)
										end do

										do dim1 = 1, dimensions
											if(dim1 == dim) then
												v%pr(dim1)%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3)) = - v%pr(dim1)%cells(i,j,k)
											else
												v%pr(dim1)%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3)) = v%pr(dim1)%cells(i,j,k)
											end if
										end do										
									
										farfield_pressure	= bc%boundary_types(bound_number)%get_farfield_pressure()
										farfield_density	= bc%boundary_types(bound_number)%get_farfield_density()
										v%pr(dim)%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3)) = sign*sqrt(abs((p%cells(i,j,k) - farfield_pressure)*(rho%cells(i,j,k) - farfield_density)/farfield_density/rho%cells(i,j,k)))
									case ('outlet2')
									
										p%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3))				= p%cells(i,j,k)
										p_dyn%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3))			= p_dyn%cells(i,j,k)
										p_stat%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3))			= p_stat%cells(i,j,k)
										rho%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3))			= rho%cells(i,j,k)
										T%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3))				= T%cells(i,j,k)
										mol_mix_conc%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3))	= mol_mix_conc%cells(i,j,k)
								
										do spec = 1, species_number
											Y%pr(spec)%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3))	=	Y%pr(spec)%cells(i,j,k)
										end do
									
										!farfield_pressure	= bc%boundary_types(bound_number)%get_farfield_pressure()
										!farfield_density	= bc%boundary_types(bound_number)%get_farfield_density()
										v%pr(dim)%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3)) = v%pr(dim)%cells(i,j,k)
									case ('inlet')
										farfield_pressure	= bc%boundary_types(bound_number)%get_farfield_pressure()
										p_dyn%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3))			= p_dyn%cells(i,j,k)!farfield_pressure - p_stat%cells(i,j,k)!p_dyn%cells(i,j,k)
										p_stat%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3))			= p_stat%cells(i,j,k)
								end select
								
							end if
						end do
					end do
				end if
			end do
			end do
			end do

		!$omp end do nowait
		!$omp end parallel			
			
		end associate

	end subroutine


	subroutine apply_boundary_conditions_interm_E_Y(this)

		class(cpm_low_mach_solver)		,intent(inout)		:: this

		integer	:: dimensions, species_number
		integer	,dimension(3,2)	:: cons_inner_loop

		character(len=20)		:: boundary_type_name	
		
		integer	:: sign, bound_number
		integer :: i,j,k,plus,dim,dim1,spec
								
		dimensions		= this%domain%get_domain_dimensions()
		species_number	= this%chem%chem_ptr%species_number
				
		cons_inner_loop	= this%domain%get_local_inner_cells_bounds()			
		
		associate(  E_f				=> this%E_f%s_ptr		, &
					E_f_int			=> this%E_f_int%s_ptr	, &
					Y				=> this%Y%v_ptr			, &
					Y_int			=> this%Y_int%v_ptr		, &
					bc				=> this%boundary%bc_ptr	, &
					mesh			=> this%mesh%mesh_ptr)

		!$omp parallel default(none)  private(i,j,k,plus,dim,dim1,spec,sign,bound_number,boundary_type_name) , &
		!$omp& firstprivate(this)	,&
		!$omp& shared(E_f,Y,E_f_int,Y_int,bc,cons_inner_loop,species_number,dimensions)
		!$omp do collapse(3) schedule(guided)					
					
			do k = cons_inner_loop(3,1),cons_inner_loop(3,2)
			do j = cons_inner_loop(2,1),cons_inner_loop(2,2)
			do i = cons_inner_loop(1,1),cons_inner_loop(1,2)
				if(bc%bc_markers(i,j,k) == 0) then
					do dim = 1,dimensions
						do plus = 1,2
							sign			= (-1)**plus
							bound_number	= bc%bc_markers(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3))
							if( bound_number /= 0 ) then
							
								boundary_type_name = bc%boundary_types(bound_number)%get_type_name()
								select case(boundary_type_name)
									case ('wall')
										do dim1 = 1, dimensions
											E_f_int%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3))		= E_f_int%cells(i,j,k)
								
											do spec = 1, species_number
												Y_int%pr(spec)%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3))	=	Y_int%pr(spec)%cells(i,j,k)
											end do
										end do									
									case ('outlet')
										E_f_int%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3))		= E_f_int%cells(i,j,k)
								
										do spec = 1, species_number
											Y_int%pr(spec)%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3))	=	Y_int%pr(spec)%cells(i,j,k)
										end do										
									case ('outlet2')
										E_f_int%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3))		= E_f_int%cells(i,j,k)
								
										do spec = 1, species_number
											Y_int%pr(spec)%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3))	=	Y_int%pr(spec)%cells(i,j,k)
										end do										
									case ('inlet')
										E_f_int%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3))		= E_f%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3))
								
										do spec = 1, species_number
											Y_int%pr(spec)%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3))	=	Y%pr(spec)%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3))
										end do										
								end select							

							end if
						end do
					end do
				end if
			end do
			end do
			end do

		!$omp end do nowait
		!$omp end parallel
			
		end associate

	end subroutine

	subroutine apply_boundary_conditions_interm_v(this)

		class(cpm_low_mach_solver)		,intent(inout)		:: this

		integer	:: dimensions, species_number
		integer	,dimension(3,2)	:: cons_inner_loop
		
		character(len=20)		:: boundary_type_name		

		integer	:: sign, bound_number
		integer :: i,j,k,plus,dim,dim1,specie_number
								
		dimensions		= this%domain%get_domain_dimensions()
		species_number	= this%chem%chem_ptr%species_number
				
		cons_inner_loop	= this%domain%get_local_inner_cells_bounds()			
		
		associate(  v				=> this%v%v_ptr			, &
					v_int			=> this%v_int%v_ptr		, &
					bc				=> this%boundary%bc_ptr	, &
					mesh			=> this%mesh%mesh_ptr)

		!$omp parallel default(none)  private(i,j,k,plus,dim,dim1,sign,bound_number,boundary_type_name) , &
		!$omp& firstprivate(this)	,&
		!$omp& shared(v_int,bc,cons_inner_loop,dimensions)
		!$omp do collapse(3) schedule(guided)					
					
			do k = cons_inner_loop(3,1),cons_inner_loop(3,2)
			do j = cons_inner_loop(2,1),cons_inner_loop(2,2)
			do i = cons_inner_loop(1,1),cons_inner_loop(1,2)
				if(bc%bc_markers(i,j,k) == 0) then
					do dim = 1,dimensions
						do plus = 1,2
							sign			= (-1)**plus
							bound_number	= bc%bc_markers(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3))
							if( bound_number /= 0 ) then

								boundary_type_name = bc%boundary_types(bound_number)%get_type_name()
								select case(boundary_type_name)
									case ('wall')
										do dim1 = 1, dimensions
											if(dim1 == dim) then
												v_int%pr(dim1)%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3)) = - v_int%pr(dim1)%cells(i,j,k)
											else
												v_int%pr(dim1)%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3)) = v_int%pr(dim1)%cells(i,j,k)
											end if
										end do									
									case ('outlet')
										do dim1 = 1, dimensions
											v_int%pr(dim1)%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3)) = v_int%pr(dim1)%cells(i,j,k)
										end do											
									case ('outlet2')
										do dim1 = 1, dimensions
											v_int%pr(dim1)%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3)) = v_int%pr(dim1)%cells(i,j,k)
										end do											
									case ('inlet')
										do dim1 = 1, dimensions
											if(dim1 == dim)	v_int%pr(dim1)%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3)) = v%pr(dim1)%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3))
										end do												
								end select

							end if
						end do
					end do
				end if
			end do
			end do
			end do

		!$omp end do nowait
		!$omp end parallel			
			
		end associate

	end subroutine

	subroutine calculate_time_step(this)
		class(cpm_low_mach_solver)	,intent(inout)	:: this
		
		real(dkind)	:: delta_t_interm, time_step, velocity_value

		real(dkind)	,dimension(3)	:: cell_size

		integer	:: dimensions, species_number
		integer	,dimension(3,2)	:: cons_inner_loop
		
		integer	:: sign
		integer :: i,j,k,dim

		time_step		= 2.0e-06_dkind !this%initial_time_step
								
		dimensions		= this%domain%get_domain_dimensions()
		species_number	= this%chem%chem_ptr%species_number
				
		cell_size		= this%mesh%mesh_ptr%get_cell_edges_length()
		
		
		cons_inner_loop	= this%domain%get_local_inner_cells_bounds()			
		
		
		associate(  v				=> this%v_int%v_ptr		, &
					v_s				=> this%v_s%s_ptr		, &
					bc				=> this%boundary%bc_ptr	, &
					mesh			=> this%mesh%mesh_ptr)
		
		!$omp parallel default(none)  private(i,j,k,dim,delta_t_interm,velocity_value) , &
		!$omp& firstprivate(this)	,&
		!$omp& shared(v,v_s,bc,time_step,cons_inner_loop,dimensions,cell_size)
		!$omp do collapse(3) schedule(guided) reduction(min:time_step)						
						
		do k = cons_inner_loop(3,1),cons_inner_loop(3,2)
		do j = cons_inner_loop(2,1),cons_inner_loop(2,2)
		do i = cons_inner_loop(1,1),cons_inner_loop(1,2)
			if(bc%bc_markers(i,j,k) == 0) then
				velocity_value		= 0.0_dkind
				do dim = 1,dimensions
					velocity_value = velocity_value + v%pr(dim)%cells(i,j,k)*v%pr(dim)%cells(i,j,k)
				end do
!				delta_t_interm = minval(cell_size,cell_size > 0.0_dkind) / (sqrt(velocity_value) + v_s%cells(i,j,k))
				delta_t_interm = minval(cell_size,cell_size > 0.0_dkind) / (sqrt(velocity_value) + 10.0_dkind)
				if (delta_t_interm < time_step) then
					time_step = delta_t_interm
				end if
			end if
		end do
		end do
		end do

		!$omp end do nowait
		!$omp end parallel			
			
		this%time_step = 2.0e-06_dkind !this%initial_time_step			
			
		end associate			

	end subroutine
		
	
	pure function get_time_step(this)
		real(dkind)						:: get_time_step
		class(cpm_low_mach_solver)	,intent(in)		:: this

		get_time_step = this%time_step
	end function

	pure function get_time(this)
		real(dkind)						:: get_time
		class(cpm_low_mach_solver)	,intent(in)		:: this

		get_time = this%time
	end function

end module
