module cpm_solver_class

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
	use thermophysical_properties_class
	
	use viscosity_solver_class
	use coarse_particles_method
	use fourier_heat_transfer_solver_class
	use fickean_diffusion_solver_class
	
	use lagrangian_particles_solver_class
	use lagrangian_droplets_solver_class	
	
	use droplets_solver_class
	use particles_solver_class	
    
	use mpi_communications_class

	use solver_options_class
	
	implicit none

#ifdef OMP	
	include "omp_lib.h"
#endif
	
	private
	public	:: cpm_solver, cpm_solver_c

	type(field_scalar_cons)	,target	:: E_f_int
	type(field_vector_cons)	,target	:: v_int, Y_int
	
    real(dkind)	,dimension(:)	,allocatable	:: flame_front_coords
    integer	:: flame_loc_unit    

	type cpm_solver
		logical			:: diffusion_flag, viscosity_flag, heat_trans_flag, reactive_flag, sources_flag, hydrodynamics_flag, multiphase_flag, CFL_condition_flag
		real(dkind)		:: courant_fraction
		real(dkind)		:: time, time_step, initial_time_step
		integer			:: additional_particles_phases_number, additional_droplets_phases_number
		
		type(viscosity_solver)				:: visc_solver
		type(coarse_particles)				:: gas_dynamics_solver
		type(heat_transfer_solver)			:: heat_trans_solver
		type(diffusion_solver)				:: diff_solver
		type(chemical_kinetics_solver)		:: chem_kin_solver
		type(table_approximated_real_gas)	:: state_eq
		
!		type(lagrangian_droplets_solver), dimension(:)	    ,allocatable	:: droplets_solver			!# Lagrangian droplets solver
		type(droplets_solver)			, dimension(:)	    ,allocatable	:: droplets_solver			!# Continuum droplets solver

!		type(lagrangian_particles_solver), dimension(:)	    ,allocatable	:: particles_solver			!# Lagrangian particles solver
		type(particles_solver)			, dimension(:)	    ,allocatable	:: particles_solver			!# Continuum particles solver		

		type(computational_domain)					:: domain
		type(mpi_communications)					:: mpi_support		
		type(chemical_properties_pointer)			:: chem
		type(thermophysical_properties_pointer)		:: thermo
		type(computational_mesh_pointer)			:: mesh
		type(boundary_conditions_pointer)			:: boundary

		type(field_scalar_cons_pointer)	:: rho	, T				, p				, v_s			,mol_mix_conc
		type(field_scalar_cons_pointer)	:: E_f	, E_f_prod_chem	, E_f_prod_heat	, E_f_prod_gd	, E_f_prod_visc	, E_f_prod_diff, E_f_int
		type(field_vector_cons_pointer)	:: v	, v_prod_gd		, v_prod_visc	, v_prod_source	, v_int
		type(field_vector_cons_pointer)	:: Y	, Y_prod_diff	, Y_prod_chem	, Y_int
		type(field_vector_flow_pointer)	:: v_f
		
		type(field_scalar_cons_pointer)	,dimension(:)	,allocatable	:: rho_prod_droplets, E_f_prod_droplets
		type(field_scalar_cons_pointer)	,dimension(:)	,allocatable	:: E_f_prod_particles 
		type(field_vector_cons_pointer)	,dimension(:)	,allocatable	:: Y_prod_droplets
		
!		type(field_vector_flow_pointer)	,dimension(:)	,allocatable	::  v_prod_droplets				!# Lagrangian droplets solver
		type(field_vector_cons_pointer)	,dimension(:)	,allocatable	::  v_prod_droplets				!# Continuum droplets solver
		
!		type(field_vector_flow_pointer)	,dimension(:)	,allocatable	::  v_prod_particles			!# Lagrangian particles solver
		type(field_vector_cons_pointer)	,dimension(:)	,allocatable	::  v_prod_particles			!# Continuum particles solver		
	contains
		procedure	,private	:: apply_boundary_conditions_main
		procedure	,private	:: apply_boundary_conditions_interm_v
		procedure	,private	:: apply_boundary_conditions_interm_E_Y		
		procedure	,private	:: calculate_interm_v
		procedure	,private	:: calculate_interm_E_Y
		procedure				:: solve_problem
		procedure				:: calculate_time_step
		procedure				:: set_CFL_coefficient	
		procedure				:: get_CFL_coefficient		
		procedure				:: get_time_step
		procedure				:: get_time
        procedure	,private	:: if_stabilized
	end type

	interface	cpm_solver_c
		module procedure	constructor
	end interface

contains

	type(cpm_solver)	function constructor(manager,problem_data_io, problem_solver_options)
		type(data_manager)						,intent(inout)	:: manager
		type(data_io)							,intent(inout)	:: problem_data_io
		type(solver_options)					,intent(in)		:: problem_solver_options

		real(dkind)	:: calculation_time
		
		type(field_scalar_cons_pointer)	:: scal_ptr
		type(field_vector_cons_pointer)	:: vect_ptr
		type(field_tensor_cons_pointer)	:: tens_ptr
		
		type(field_scalar_flow_pointer)	:: scal_f_ptr		
		type(field_vector_flow_pointer)	:: vect_f_ptr		
		
		type(solid_particles_phase)		:: particles_params
        type(liquid_droplets_phase)     :: droplets_params
		integer				:: particles_phase_counter, droplets_phase_counter
		character(len=40)	:: var_name
		
        integer	,dimension(3,2)	:: cons_allocation_bounds

		constructor%diffusion_flag		= problem_solver_options%get_molecular_diffusion_flag()
		constructor%viscosity_flag		= problem_solver_options%get_viscosity_flag()
		constructor%heat_trans_flag		= problem_solver_options%get_heat_transfer_flag()
		constructor%reactive_flag		= problem_solver_options%get_chemical_reaction_flag()
		constructor%hydrodynamics_flag	= problem_solver_options%get_hydrodynamics_flag()
		constructor%courant_fraction	= problem_solver_options%get_CFL_condition_coefficient()
		constructor%CFL_condition_flag	= problem_solver_options%get_CFL_condition_flag()
		constructor%sources_flag		= .false.
		constructor%additional_particles_phases_number	= problem_solver_options%get_additional_particles_phases_number()
		constructor%additional_droplets_phases_number	= problem_solver_options%get_additional_droplets_phases_number()
		
		
		constructor%domain				= manager%domain
		constructor%mpi_support			= manager%mpi_communications		
		constructor%chem%chem_ptr		=> manager%chemistry%chem_ptr
		constructor%thermo%thermo_ptr	=> manager%thermophysics%thermo_ptr
		constructor%boundary%bc_ptr	=> manager%boundary_conditions_pointer%bc_ptr
		constructor%mesh%mesh_ptr		=> manager%computational_mesh_pointer%mesh_ptr
        
        cons_allocation_bounds		= manager%domain%get_local_utter_cells_bounds()

		call manager%get_cons_field_pointer_by_name(scal_ptr,vect_ptr,tens_ptr,'density')
		constructor%rho%s_ptr					=> scal_ptr%s_ptr
		call manager%get_cons_field_pointer_by_name(scal_ptr,vect_ptr,tens_ptr,'temperature')
		constructor%T%s_ptr					=> scal_ptr%s_ptr
		call manager%get_cons_field_pointer_by_name(scal_ptr,vect_ptr,tens_ptr,'pressure')
		constructor%p%s_ptr					=> scal_ptr%s_ptr
		call manager%get_cons_field_pointer_by_name(scal_ptr,vect_ptr,tens_ptr,'full_energy')
		constructor%E_f%s_ptr					=> scal_ptr%s_ptr
		call manager%get_cons_field_pointer_by_name(scal_ptr,vect_ptr,tens_ptr,'velocity')
		constructor%v%v_ptr						=> vect_ptr%v_ptr
		call manager%get_cons_field_pointer_by_name(scal_ptr,vect_ptr,tens_ptr,'specie_molar_concentration')
		constructor%Y%v_ptr						=> vect_ptr%v_ptr
		call manager%get_cons_field_pointer_by_name(scal_ptr,vect_ptr,tens_ptr,'mixture_molar_concentration')
		constructor%mol_mix_conc%s_ptr		=> scal_ptr%s_ptr	
		
		call manager%get_flow_field_pointer_by_name(scal_f_ptr,vect_f_ptr,'velocity_flow')
		constructor%v_f%v_ptr				=> vect_f_ptr%v_ptr		

		call manager%create_scalar_field(E_f_int	,'full_energy_interm'					,'E_f_int')
		constructor%E_f_int%s_ptr				=> E_f_int
		call manager%create_vector_field(v_int		,'velocity_interm'						,'v_int'	,'spatial')
		constructor%v_int%v_ptr					=> v_int
		call manager%create_vector_field(Y_int		,'specie_molar_concentration_interm'	,'Y_int'	,'chemical')
		constructor%Y_int%v_ptr					=> Y_int
		constructor%gas_dynamics_solver	= coarse_particles_c(manager)
		
		call manager%get_cons_field_pointer_by_name(scal_ptr,vect_ptr,tens_ptr,'velocity_production_gas_dynamics')
		constructor%v_prod_gd%v_ptr				=> vect_ptr%v_ptr
		call manager%get_cons_field_pointer_by_name(scal_ptr,vect_ptr,tens_ptr,'energy_production_gas_dynamics')
		constructor%E_f_prod_gd%s_ptr			=> scal_ptr%s_ptr

		!constructor%sources	= sources_c(manager)
		!call manager%get_field_pointer(scal_ptr,vect_ptr,tens_ptr,'velocity_production_sources')
		!constructor%v_prod_source%v_ptr			=> vect_ptr%v_ptr

		constructor%state_eq	=	table_approximated_real_gas_c(manager)
		call manager%get_cons_field_pointer_by_name(scal_ptr,vect_ptr,tens_ptr,'velocity_of_sound')
		constructor%v_s%s_ptr			=> scal_ptr%s_ptr		

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

		if(constructor%additional_droplets_phases_number /= 0) then
			allocate(constructor%droplets_solver(constructor%additional_droplets_phases_number))
			call constructor%droplets_solver(1)%pre_constructor(constructor%additional_droplets_phases_number)
			allocate(constructor%rho_prod_droplets(constructor%additional_droplets_phases_number))
			allocate(constructor%E_f_prod_droplets(constructor%additional_droplets_phases_number))
			allocate(constructor%v_prod_droplets(constructor%additional_droplets_phases_number))
			allocate(constructor%Y_prod_droplets(constructor%additional_droplets_phases_number))
			do droplets_phase_counter = 1, constructor%additional_droplets_phases_number
				droplets_params = problem_solver_options%get_droplets_params(droplets_phase_counter)
!				constructor%droplets_solver(droplets_phase_counter)	= lagrangian_droplets_solver_c(manager, droplets_params, droplets_phase_counter)		!# Lagrangian droplets solver
				constructor%droplets_solver(droplets_phase_counter)	= droplets_solver_c(manager, droplets_params, droplets_phase_counter)					!# Continuum droplets solver
				write(var_name,'(A,I2.2)') 'energy_production_droplets', droplets_phase_counter
				call manager%get_cons_field_pointer_by_name(scal_ptr,vect_ptr,tens_ptr,var_name)
				constructor%E_f_prod_droplets(droplets_phase_counter)%s_ptr	=> scal_ptr%s_ptr
				write(var_name,'(A,I2.2)') 'density_production_droplets', droplets_phase_counter
				call manager%get_cons_field_pointer_by_name(scal_ptr,vect_ptr,tens_ptr,var_name)
				constructor%rho_prod_droplets(droplets_phase_counter)%s_ptr	=> scal_ptr%s_ptr                
				write(var_name,'(A,I2.2)') 'velocity_production_droplets', droplets_phase_counter						
!				call manager%get_flow_field_pointer_by_name(scal_f_ptr,vect_f_ptr,var_name)								!# Lagrangian droplets solver
!				constructor%v_prod_droplets(droplets_phase_counter)%v_ptr		=> vect_f_ptr%v_ptr						!# Lagrangian droplets solver
				call manager%get_cons_field_pointer_by_name(scal_ptr,vect_ptr,tens_ptr,var_name)						!# Continuum droplets solver								
				constructor%v_prod_droplets(droplets_phase_counter)%v_ptr		=> vect_ptr%v_ptr						!# Continuum droplets solver
				write(var_name,'(A,I2.2)') 'concentration_production_droplets', droplets_phase_counter
				call manager%get_cons_field_pointer_by_name(scal_ptr,vect_ptr,tens_ptr,var_name)
				constructor%Y_prod_droplets(droplets_phase_counter)%v_ptr		=> vect_ptr%v_ptr                
			end do		
		end if 	
		
		if(constructor%additional_particles_phases_number /= 0) then
			allocate(constructor%particles_solver(constructor%additional_particles_phases_number))
			call constructor%particles_solver(1)%pre_constructor(constructor%additional_particles_phases_number)			
			allocate(constructor%E_f_prod_particles(constructor%additional_particles_phases_number))
			allocate(constructor%v_prod_particles(constructor%additional_particles_phases_number))
			do particles_phase_counter = 1, constructor%additional_particles_phases_number
				particles_params = problem_solver_options%get_particles_params(particles_phase_counter)
!				constructor%particles_solver(particles_phase_counter)	= lagrangian_particles_solver_c(manager, particles_params, particles_phase_counter)		!# Lagrangian particles solver
				constructor%particles_solver(particles_phase_counter)	= particles_solver_c(manager, particles_params, particles_phase_counter)				!# Continuum particles solver
				write(var_name,'(A,I2.2)') 'energy_production_particles', particles_phase_counter
				call manager%get_cons_field_pointer_by_name(scal_ptr,vect_ptr,tens_ptr,var_name)
				constructor%E_f_prod_particles(particles_phase_counter)%s_ptr	=> scal_ptr%s_ptr
				write(var_name,'(A,I2.2)') 'velocity_production_particles', particles_phase_counter
!				call manager%get_flow_field_pointer_by_name(scal_f_ptr,vect_f_ptr,var_name)								!# Lagrangian particles solver
!				constructor%v_prod_particles(particles_phase_counter)%v_ptr		=> vect_f_ptr%v_ptr						!# Lagrangian particles solver				
				call manager%get_cons_field_pointer_by_name(scal_ptr,vect_ptr,tens_ptr,var_name)						!# Continuum particles solver
				constructor%v_prod_particles(particles_phase_counter)%v_ptr		=> vect_ptr%v_ptr						!# Continuum particles solver
			end do		
        end if		
		
		
		if (constructor%reactive_flag) then
			constructor%chem_kin_solver		= chemical_kinetics_solver_c(manager)
			call manager%get_cons_field_pointer_by_name(scal_ptr,vect_ptr,tens_ptr,'energy_production_chemistry')
			constructor%E_f_prod_chem%s_ptr			=> scal_ptr%s_ptr
			call manager%get_cons_field_pointer_by_name(scal_ptr,vect_ptr,tens_ptr,'specie_production_chemistry')
			constructor%Y_prod_chem%v_ptr			=> vect_ptr%v_ptr
		end if		
		
		problem_data_io		= data_io_c(manager,calculation_time)
		
		if(problem_data_io%get_load_counter() /= 0) then
			call problem_data_io%add_io_scalar_cons_field(constructor%E_f)
		end if
		
		call problem_data_io%input_all_data()
			
		if(problem_data_io%get_load_counter() == 1) then
			call problem_data_io%add_io_scalar_cons_field(constructor%E_f)
		end if
		
		if(problem_data_io%get_load_counter() == 1) then
			call constructor%state_eq%apply_state_equation_for_initial_conditions()
			if(constructor%additional_particles_phases_number /= 0) then
				do particles_phase_counter = 1, constructor%additional_particles_phases_number
					call constructor%particles_solver(particles_phase_counter)%set_initial_distributions()
				end do
            end if
			if(constructor%additional_droplets_phases_number /= 0) then
				do droplets_phase_counter = 1, constructor%additional_droplets_phases_number
					call constructor%droplets_solver(droplets_phase_counter)%set_initial_distributions()
				end do
			end if            
		else
			call constructor%state_eq%apply_state_equation()
			call constructor%state_eq%apply_boundary_conditions_for_initial_conditions()
		end if		
		
		
		call constructor%mpi_support%exchange_conservative_scalar_field(constructor%p%s_ptr)
		call constructor%mpi_support%exchange_conservative_scalar_field(constructor%rho%s_ptr)
		call constructor%mpi_support%exchange_conservative_scalar_field(constructor%T%s_ptr)

		call constructor%mpi_support%exchange_conservative_vector_field(constructor%Y%v_ptr)
		call constructor%mpi_support%exchange_conservative_vector_field(constructor%v%v_ptr)

		call constructor%mpi_support%exchange_boundary_conditions_markers(constructor%boundary%bc_ptr)
		call constructor%mpi_support%exchange_mesh(constructor%mesh%mesh_ptr)

		constructor%time		=	calculation_time
		constructor%time_step	=	problem_solver_options%get_initial_time_step()
		constructor%initial_time_step = problem_solver_options%get_initial_time_step()

        allocate(flame_front_coords(cons_allocation_bounds(2,1):cons_allocation_bounds(2,2)))
        
        open(newunit = flame_loc_unit, file = 'av_flame_data.dat', status = 'replace', form = 'formatted')
        
	end function

	subroutine solve_problem(this)
		class(cpm_solver)	,intent(inout)	:: this

		integer	:: particles_phase_counter, droplets_phase_counter
		integer	:: specie
        
        logical	:: stabilized

		this%time = this%time + this%time_step		
		
		call this%apply_boundary_conditions_main()
		if(this%additional_particles_phases_number /= 0) then
			do particles_phase_counter = 1, this%additional_particles_phases_number
    			call this%particles_solver(particles_phase_counter)%apply_boundary_conditions_main()
            end do		
        end if
        if(this%additional_droplets_phases_number /= 0) then
			do droplets_phase_counter = 1, this%additional_droplets_phases_number
    			call this%droplets_solver(droplets_phase_counter)%apply_boundary_conditions_main(this%time)
			end do	            
		end if
        
		if (this%viscosity_flag)	call this%visc_solver%solve_viscosity(this%time_step)
		if (this%heat_trans_flag)	call this%heat_trans_solver%solve_heat_transfer(this%time_step)
		if (this%diffusion_flag)	call this%diff_solver%solve_diffusion(this%time_step)
		if (this%reactive_flag)		call this%chem_kin_solver%solve_chemical_kinetics(this%time_step)
		!call this%sources%calculate_sources(this%time,this%time_step)
		
		if(this%additional_droplets_phases_number /= 0) then
			do droplets_phase_counter = 1, this%additional_droplets_phases_number
!				call this%droplets_solver(droplets_phase_counter)%droplets_solve(this%time_step)				!# Lagrangian droplets solver
				call this%droplets_solver(droplets_phase_counter)%droplets_euler_step_v_E(this%time_step)		!# Continuum droplets solver
				call this%droplets_solver(droplets_phase_counter)%apply_boundary_conditions_interm_v_d()		!# Continuum droplets solver
				call this%droplets_solver(droplets_phase_counter)%droplets_lagrange_step(this%time_step)		!# Continuum droplets solver
				call this%droplets_solver(droplets_phase_counter)%droplets_final_step(this%time_step)			!# Continuum droplets solver		
			end do		
		end if 		
		
		call this%gas_dynamics_solver%euler_step_v(this%time_step)
        
		if(this%additional_particles_phases_number /= 0) then
			do particles_phase_counter = 1, this%additional_particles_phases_number
!				call this%particles_solver(particles_phase_counter)%particles_solve(this%time_step)				!# Lagrangian particles solver
				call this%particles_solver(particles_phase_counter)%particles_euler_step_v_E(this%time_step)	!# Continuum particles solver
				call this%particles_solver(particles_phase_counter)%apply_boundary_conditions_interm_v_p()		!# Continuum particles solver
				call this%particles_solver(particles_phase_counter)%particles_lagrange_step(this%time_step)		!# Continuum particles solver
				call this%particles_solver(particles_phase_counter)%particles_final_step(this%time_step)		!# Continuum particles solver	
			end do		
        end if

		call this%calculate_interm_v(this%time_step)
		call this%apply_boundary_conditions_interm_v()

		call this%gas_dynamics_solver%euler_step_E(this%time_step)
		call this%calculate_interm_E_Y(this%time_step)
		call this%apply_boundary_conditions_interm_E_Y()
		
		call this%gas_dynamics_solver%lagrange_step(this%time_step)
		call this%gas_dynamics_solver%final_step(this%time_step)

		call this%state_eq%apply_state_equation()

		if (this%CFL_condition_flag) then
			call this%calculate_time_step()
		end if		
		
        call this%if_stabilized(this%time, stabilized)
        if (stabilized) stop       
        
		!call this%state_eq%check_conservation_laws()

	end subroutine

	subroutine calculate_interm_v(this,time_step)
		class(cpm_solver)	,intent(inout)	:: this
		real(dkind)			,intent(in)		:: time_step

		integer	:: particles_phase_counter, droplets_phase_counter
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
					v_prod_particles	=> this%v_prod_particles	, &
                    v_prod_droplets	=> this%v_prod_droplets	, &
					bc				=> this%boundary%bc_ptr)

		!$omp parallel default(none)  private(i,j,k,dim,particles_phase_counter) , &
		!$omp& firstprivate(this)	,&
		!$omp& shared(rho,v,v_int,v_prod_gd,v_prod_visc,v_prod_source,v_prod_particles,bc,cons_inner_loop,dimensions,time_step)
		!$omp do collapse(3) schedule(guided)

			do k = cons_inner_loop(3,1),cons_inner_loop(3,2)
			do j = cons_inner_loop(2,1),cons_inner_loop(2,2)
			do i = cons_inner_loop(1,1),cons_inner_loop(1,2)
				if(bc%bc_markers(i,j,k) == 0) then					
					do dim = 1, dimensions
						v_int%pr(dim)%cells(i,j,k) = v%pr(dim)%cells(i,j,k)
						
						if (this%hydrodynamics_flag)	v_int%pr(dim)%cells(i,j,k) = v_int%pr(dim)%cells(i,j,k) + v_prod_gd%pr(dim)%cells(i,j,k)! + v_prod_source%pr(dim)%cells(i,j,k)
						
						if (this%viscosity_flag)		v_int%pr(dim)%cells(i,j,k) = v_int%pr(dim)%cells(i,j,k) + v_prod_visc%pr(dim)%cells(i,j,k) / rho%cells(i,j,k) * time_step
						
						if (this%additional_droplets_phases_number /= 0) then
							do droplets_phase_counter = 1, this%additional_droplets_phases_number
								v_int%pr(dim)%cells(i,j,k) = v_int%pr(dim)%cells(i,j,k) + v_prod_droplets(droplets_phase_counter)%v_ptr%pr(dim)%cells(i,j,k) * time_step																		!# Continuum droplets solver
!								v_int%pr(dim)%cells(i,j,k) =  v_int%pr(dim)%cells(i,j,k)	!+ 0.5_dkind*(v_prod_droplets(droplets_phase_counter)%v_ptr%pr(dim)%cells(dim,i,j,k) &
																							!			+ v_prod_droplets(droplets_phase_counter)%v_ptr%pr(dim)%cells(dim,i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3))) * time_step		!# Lagrangian droplets solver
							end do		
						end if						
						
						if (this%additional_particles_phases_number /= 0) then
							do particles_phase_counter = 1, this%additional_particles_phases_number
!								v_int%pr(dim)%cells(i,j,k) = v_int%pr(dim)%cells(i,j,k) + v_prod_particles(particles_phase_counter)%v_ptr%pr(dim)%cells(i,j,k)																		!# Continuum particles solver
								v_int%pr(dim)%cells(i,j,k) =  v_int%pr(dim)%cells(i,j,k)	!+ 0.5_dkind*(v_prod_particles(particles_phase_counter)%v_ptr%pr(dim)%cells(dim,i,j,k) & 
																							!			+ v_prod_particles(particles_phase_counter)%v_ptr%pr(dim)%cells(dim,i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3))) * time_step	!# Lagrangian particles solver
							end do		
                        end if
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
		class(cpm_solver)	,intent(inout)	:: this
		real(dkind)			,intent(in)		:: time_step
		
		integer	:: particles_phase_counter, droplets_phase_counter
		integer	:: dimensions, species_number
		integer	,dimension(3,2)	:: cons_inner_loop

		integer	:: bound_number
		integer :: i,j,k,dim,spec

		real(dkind)		:: spec_summ
		real(dkind)		:: energy_source = 397500000.0_dkind
		real(dkind)	,save	:: energy_output = 0.0_dkind
				
		dimensions		= this%domain%get_domain_dimensions()
		species_number	= this%chem%chem_ptr%species_number
				
		cons_inner_loop	= this%domain%get_local_inner_cells_bounds()		
		
		associate (	rho				=> this%rho%s_ptr			, &
                    rho_prod_droplets   => this%rho_prod_droplets	, &
					E_f				=> this%E_f%s_ptr			, &
					E_f_int 		=> this%E_f_int%s_ptr		, &
					E_f_prod_chem 	=> this%E_f_prod_chem%s_ptr	, &
					E_f_prod_heat	=> this%E_f_prod_heat%s_ptr	, &
					E_f_prod_gd 	=> this%E_f_prod_gd%s_ptr	, &
					E_f_prod_visc	=> this%E_f_prod_visc%s_ptr	, &
					E_f_prod_diff	=> this%E_f_prod_diff%s_ptr	, &
					E_f_prod_particles => this%E_f_prod_particles	,&
                    E_f_prod_droplets => this%E_f_prod_droplets	,&
					Y				=> this%Y%v_ptr				, &
					Y_int			=> this%Y_int%v_ptr			, &
					Y_prod_diff		=> this%Y_prod_diff%v_ptr	, &
					Y_prod_chem		=> this%Y_prod_chem%v_ptr	, &
                    Y_prod_droplets	=> this%Y_prod_droplets	, &
					T				=> this%T%s_ptr				, &
					bc				=> this%boundary%bc_ptr)
					
		!$omp parallel default(none)  private(i,j,k,dim,particles_phase_counter,spec,spec_summ,energy_output) , &
		!$omp& firstprivate(this)	,&
		!$omp& shared(rho,E_f,E_f_int,E_f_prod_chem,E_f_prod_heat,E_f_prod_gd,E_f_prod_visc,E_f_prod_diff,E_f_prod_particles,Y,Y_int,Y_prod_diff,Y_prod_chem,dimensions,species_number,bc,cons_inner_loop,energy_source,time_step)
		!$omp do collapse(3) schedule(guided)

			do k = cons_inner_loop(3,1),cons_inner_loop(3,2)
			do j = cons_inner_loop(2,1),cons_inner_loop(2,2)
			do i = cons_inner_loop(1,1),cons_inner_loop(1,2)
			
				if(bc%bc_markers(i,j,k) == 0) then					
					E_f_int%cells(i,j,k) = E_f%cells(i,j,k) 

					if (this%hydrodynamics_flag)	E_f_int%cells(i,j,k) = E_f_int%cells(i,j,k) + E_f_prod_gd%cells(i,j,k)
					
					if (this%viscosity_flag)	E_f_int%cells(i,j,k) = E_f_int%cells(i,j,k) + E_f_prod_visc%cells(i,j,k) * time_step / rho%cells(i,j,k) 
					if (this%heat_trans_flag)	E_f_int%cells(i,j,k) = E_f_int%cells(i,j,k) + E_f_prod_heat%cells(i,j,k) * time_step / rho%cells(i,j,k) 
					if (this%diffusion_flag)	E_f_int%cells(i,j,k) = E_f_int%cells(i,j,k) + E_f_prod_diff%cells(i,j,k) * time_step / rho%cells(i,j,k) 
					if (this%reactive_flag)		E_f_int%cells(i,j,k) = E_f_int%cells(i,j,k) + E_f_prod_chem%cells(i,j,k) * time_step / rho%cells(i,j,k)

					if (this%additional_particles_phases_number /= 0) then
						do particles_phase_counter = 1, this%additional_particles_phases_number
							E_f_int%cells(i,j,k) = E_f_int%cells(i,j,k) !+ E_f_prod_particles(particles_phase_counter)%s_ptr%cells(i,j,k)
						end do		
                    end if					
					
					if (this%additional_droplets_phases_number /= 0) then
						do droplets_phase_counter = 1, this%additional_droplets_phases_number
							E_f_int%cells(i,j,k) = E_f_int%cells(i,j,k) + E_f_prod_droplets(droplets_phase_counter)%s_ptr%cells(i,j,k) * time_step / rho%cells(i,j,k)
!                            rho%cells(i,j,k)     = rho%cells(i,j,k) + rho_prod_droplets(droplets_phase_counter)%s_ptr%cells(i,j,k)
						end do		
					end if	                    
                    
					! ************************* Energy release ******************
					!if(this%time <= 1.0E-07_dkind) then
					!	if(i <= 80) then
					!		if (this%sources_flag)	then
					!			E_f_int%cells(i,j,k)	= E_f_int%cells(i,j,k)	+ energy_source * this%time_step * 1.0E06_dkind 
					!			energy_output			= energy_output			+ energy_source * this%time_step * 1.0E06_dkind *  rho%cells(i,1,1)
					!		end if
					!	end if
					!end if
					! ***********************************************************	

					spec_summ = 0.0
					do spec = 1, species_number

						Y_int%pr(spec)%cells(i,j,k) = Y%pr(spec)%cells(i,j,k)
					
						if (this%diffusion_flag)	Y_int%pr(spec)%cells(i,j,k) = Y_int%pr(spec)%cells(i,j,k) + Y_prod_diff%pr(spec)%cells(i,j,k)  * time_step / rho%cells(i,j,k)
						if (this%reactive_flag)		Y_int%pr(spec)%cells(i,j,k) = Y_int%pr(spec)%cells(i,j,k) + Y_prod_chem%pr(spec)%cells(i,j,k)  * time_step / rho%cells(i,j,k)
						
						if (this%additional_droplets_phases_number /= 0) then
							do droplets_phase_counter = 1, this%additional_droplets_phases_number
								Y_int%pr(spec)%cells(i,j,k) = Y_int%pr(spec)%cells(i,j,k) + Y_prod_droplets(droplets_phase_counter)%v_ptr%pr(spec)%cells(i,j,k) * time_step / rho%cells(i,j,k)
							end do		
						end if
							
						if (this%thermo%thermo_ptr%molar_masses(spec) /= 0.0_dkind) then
							spec_summ = spec_summ + max(Y_int%pr(spec)%cells(i,j,k), 0.0_dkind)
						end if
					end do
					
					do spec = 1,species_number
						if (this%thermo%thermo_ptr%molar_masses(spec) /= 0.0_dkind) then
					!		Y_int%pr(spec)%cells(i,j,k) = max(Y_int%pr(spec)%cells(i,j,k), 0.0_dkind) / spec_summ
						end if
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

		class(cpm_solver)		,intent(inout)		:: this

		character(len=20)		:: boundary_type_name
		real(dkind)				:: farfield_density, farfield_pressure, farfield_velocity, wall_temperature
		
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
					rho				=> this%rho%s_ptr			, &
					v				=> this%v%v_ptr				, &
					Y				=> this%Y%v_ptr				, &
					bc				=> this%boundary%bc_ptr		, &
					mesh			=> this%mesh%mesh_ptr)

		!$omp parallel default(none)  private(i,j,k,plus,dim,dim1,spec,sign,bound_number,boundary_type_name,farfield_pressure,farfield_density,farfield_velocity,wall_temperature) , &
		!$omp& firstprivate(this)	,&
		!$omp& shared(T,p,rho,v,mol_mix_conc,Y,mesh,bc,cons_inner_loop,dimensions,species_number)
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
									case('wall','symmetry_plane')
									
										p%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3))				= p%cells(i,j,k)
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
													v%pr(dim1)%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3)) = 0.0_dkind !v%pr(dim1)%cells(i-sign*I_m(dim,1),j-sign*I_m(dim,2),k-sign*I_m(dim,3)) - 6.0_dkind * v%pr(dim1)%cells(i,j,k)
												!	v%pr(dim1)%cells(i,j,k) = 0.0_dkind
												!	v%pr(dim1)%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3)) = 0.0_dkind
												end if
											end do
										end if
									case ('outlet')

										farfield_pressure	= bc%boundary_types(bound_number)%get_farfield_pressure()
										farfield_density	= bc%boundary_types(bound_number)%get_farfield_density()
									!	print *, farfield_pressure
									!	print *, farfield_density
									
										p%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3))				= p%cells(i,j,k)
										rho%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3))			= rho%cells(i,j,k)
										T%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3))				= T%cells(i,j,k)
										mol_mix_conc%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3))	= mol_mix_conc%cells(i,j,k)									
									
										v%pr(dim)%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3)) = sign*sqrt(abs((p%cells(i,j,k) - farfield_pressure)*(rho%cells(i,j,k) - farfield_density)/farfield_density/rho%cells(i,j,k)))
										do spec = 1, species_number
											Y%pr(spec)%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3))	=	Y%pr(spec)%cells(i,j,k)
										end do

										do dim1 = 1, dimensions
											if(dim1 == dim) then
												v%pr(dim1)%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3)) = v%pr(dim1)%cells(i,j,k) !max(0.0,sign*v%pr(dim1)%cells(i,j,k))
											else
												v%pr(dim1)%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3)) = v%pr(dim1)%cells(i,j,k)
											end if
										end do	
										
									case ('inlet')

										!**** Relaxing inlet ****
								!		farfield_pressure	= 101325.0_dkind + (bc%boundary_types(bound_number)%get_farfield_pressure() - 101325.0_dkind) * min(this%time/10e-06_dkind,1.0_dkind)
								!		farfield_density	= 1.17195723916649_dkind + (bc%boundary_types(bound_number)%get_farfield_density() - 1.17195723916649_dkind) * min(this%time/10e-06_dkind,1.0_dkind)
                                       
								!		farfield_pressure	= bc%boundary_types(bound_number)%get_farfield_pressure()
								!		farfield_density	= bc%boundary_types(bound_number)%get_farfield_density()	
								!		farfield_velocity	= bc%boundary_types(bound_number)%get_farfield_velocity()	
								!		p%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3))			= farfield_pressure
								!		rho%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3))		= farfield_density
                                !        E_f%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3))		= farfield_energy
                                        
								!		v%pr(dim)%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3))	= -sign*sqrt(abs((p%cells(i,j,k) - farfield_pressure)*(rho%cells(i,j,k) - farfield_density)/farfield_density/rho%cells(i,j,k)))
										continue
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

		class(cpm_solver)		,intent(inout)		:: this

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
		!$omp& shared(E_f,E_f_int,Y,Y_int,bc,cons_inner_loop,species_number,dimensions)
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
									case('wall','outlet','symmetry_plane')
										E_f_int%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3))		= E_f_int%cells(i,j,k)
								
										do spec = 1, species_number
											Y_int%pr(spec)%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3))	=	Y_int%pr(spec)%cells(i,j,k)
										end do
									case('inlet')
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

		class(cpm_solver)		,intent(inout)		:: this

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
		!$omp& shared(v,v_int,bc,cons_inner_loop,dimensions)
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
									case ('wall','symmetry_plane')
										do dim1 = 1, dimensions
											if(dim1 == dim) then
												v_int%pr(dim1)%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3)) = - v_int%pr(dim1)%cells(i,j,k)
											else
												v_int%pr(dim1)%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3)) = v_int%pr(dim1)%cells(i,j,k)
											end if
										end do
									case ('outlet')
									
										do dim1 = 1, dimensions
											if(dim1 == dim) then
												v_int%pr(dim1)%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3)) = v_int%pr(dim1)%cells(i,j,k)
											else
												v_int%pr(dim1)%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3)) = v_int%pr(dim1)%cells(i,j,k)
											end if
										end do	
										
										!if (sign == 1) then
										!	v_int%pr(dim)%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3)) = max(0.0_dkind,v_int%pr(dim)%cells(i,j,k))
										!else
										!	v_int%pr(dim)%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3)) = min(0.0_dkind,v_int%pr(dim)%cells(i,j,k))
										!end if
									case ('inlet')
										v_int%pr(dim)%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3)) = v%pr(dim)%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3))
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

#ifdef mpi
	use MPI
#endif

		class(cpm_solver)	,intent(inout)	:: this
		
		real(dkind)	:: delta_t_interm, time_step(1), velocity_value
		real(dkind)	,dimension(:)	,allocatable	,save	:: time_step_array

		integer						:: dimensions
		integer						:: processor_rank, processor_number, mpi_communicator

		integer		,dimension(3,2)	:: cons_inner_loop
		real(dkind)	,dimension(3)	:: cell_size
		integer	:: sign
		integer :: i,j,k,dim,error

		processor_rank		= this%domain%get_processor_rank()
		mpi_communicator	= this%domain%get_mpi_communicator()

		if (.not.allocated(time_step_array)) then
			processor_number = this%domain%get_mpi_communicator_size()
			allocate(time_step_array(processor_number))
		end if

		time_step(1)	= 10.0_dkind !this%initial_time_step 

		associate(  v				=> this%v%v_ptr		, &
					v_s				=> this%v_s%s_ptr		, &
					bc				=> this%boundary%bc_ptr	, &
					mesh			=> this%mesh%mesh_ptr)
		
		dimensions			= this%domain%get_domain_dimensions()
		cons_inner_loop		= this%domain%get_local_inner_cells_bounds()
		cell_size			= mesh%get_cell_edges_length()					
					
		!!$omp parallel default(shared)  private(i,j,k,dim,delta_t_interm,velocity_value) , &
		!!$omp& firstprivate(this)	,&
		!!$omp& shared(v,v_s,mesh,bc,time_step)
		!!$omp do collapse(3) schedule(static) reduction(min:time_step)
					
		do k = cons_inner_loop(3,1),cons_inner_loop(3,2)
		do j = cons_inner_loop(2,1),cons_inner_loop(2,2)
		do i = cons_inner_loop(1,1),cons_inner_loop(1,2)
			if(bc%bc_markers(i,j,k) == 0) then
				velocity_value		= 0.0_dkind
				do dim = 1,dimensions
					velocity_value = velocity_value + v%pr(dim)%cells(i,j,k)*v%pr(dim)%cells(i,j,k)
				end do
				delta_t_interm = minval(cell_size,cell_size > 0.0_dkind) / (sqrt(velocity_value) + v_s%cells(i,j,k)) !
				if (delta_t_interm < time_step(1)) then
					time_step(1) = delta_t_interm
				end if
			end if
		end do
		end do
		end do
	
		!!$omp end do nowait
		!!$omp end parallel

		time_step_array(processor_rank+1) = time_step(1) 
		
#ifdef mpi
		call mpi_gather(time_step,1,MPI_DOUBLE_PRECISION,time_step_array,1,MPI_DOUBLE_PRECISION,0,mpi_communicator,error)
#endif
		
		if (processor_rank == 0) then
			do i = 0,size(time_step_array) - 1 
				if (time_step_array(i+1) < time_step(1)) time_step(1) = time_step_array(i+1)
			end do
		end if

#ifdef mpi	
		call mpi_bcast(time_step,1,MPI_DOUBLE_PRECISION,0,mpi_communicator,error)
#endif

		this%time_step = this%courant_fraction * time_step(1)
		
		!# Implantation
		!if(this%time_step < 3.0e-08) then
		!	this%time_step = 0.025_dkind*this%time_step
		!end if
  !
		!if(this%time_step < 6.5e-10) then
		!	this%time_step = 0.5_dkind*this%time_step
		!end if
		

		end associate			
	
	end subroutine
		
	subroutine set_CFL_coefficient(this,coefficient)
		class(cpm_solver)	,intent(inout)		:: this
		real(dkind)				,intent(in)		:: coefficient
	
		this%courant_fraction = coefficient
		
	end subroutine
	
	pure function get_CFL_coefficient(this)
		real(dkind)						:: get_CFL_coefficient
		class(cpm_solver)	,intent(in)		:: this

		get_CFL_coefficient = this%courant_fraction
	end function	
	
	pure function get_time_step(this)
		real(dkind)						:: get_time_step
		class(cpm_solver)	,intent(in)		:: this

		get_time_step = this%time_step
	end function

	pure function get_time(this)
		real(dkind)						:: get_time
		class(cpm_solver)	,intent(in)		:: this

		get_time = this%time
	end function

    subroutine if_stabilized(this,time,stabilized)
		class(cpm_solver)	,intent(inout)	:: this
		real(dkind)			,intent(in)     :: time    

		logical				,intent(out)	:: stabilized
		
        logical						:: boundary 
		real(dkind)	,dimension(3)	:: cell_size		
		
		real(dkind)					:: max_val, left_val, right_val, flame_velocity, flame_surface_length, surface_factor
		real(dkind)					:: a, b 
		real(dkind)					:: time_diff, time_delay, time_stabilization
		real(dkind), save			:: previous_flame_location = 0.0_dkind, current_flame_location = 0.0_dkind, farfield_velocity = 0.0_dkind
		real(dkind), save			:: previous_time = 0.0_dkind, current_time = 0.0_dkind
        real(dkind), save			:: av_flame_velocity = 0.0_dkind, previous_av_flame_velocity = 0.0_dkind
		real(dkind), dimension(20),	save	:: flame_velocity_array = 0.0_dkind
		integer		,save			:: correction = 0, counter = 0
		integer						:: flame_front_index
		character(len=200)			:: file_name
		
		
		integer	:: dimensions, species_number
		integer	,dimension(3,2)	:: cons_inner_loop

		real(dkind)				:: tip_coord, side_coord_x, side_coord_y, T_flame, lp_dist, x_f, min_dist, min_y, max_x
		integer					:: lp_number, lp_neighbour, lp_copies, lp_tip, lp_bound, lp_start, lp_number2 
		
		integer :: CO_index, H2O2_index, HO2_index
		integer	:: bound_number,sign
		integer :: i,j,k,plus,dim,dim1,spec, lp_index,lp_index2,lp_index3
		
		
		character(len=20)		:: boundary_type_name
		
		dimensions		= this%domain%get_domain_dimensions()
		species_number	= this%chem%chem_ptr%species_number
		
		cons_inner_loop	= this%domain%get_local_inner_cells_bounds()
				
		cell_size		= this%mesh%mesh_ptr%get_cell_edges_length()

!		CO_index	= this%chem%chem_ptr%get_chemical_specie_index('CO')
!		HO2_index	= this%chem%chem_ptr%get_chemical_specie_index('HO2')
		HO2_index		= this%chem%chem_ptr%get_chemical_specie_index('HO2')
		
		stabilized = .false.
		
		associate (	v				=> this%v%v_ptr				, &
					T				=> this%T%s_ptr				, &
					Y				=> this%Y%v_ptr				, &
					bc				=> this%boundary%bc_ptr)

	!	time_delay			= 1e-04_dkind			
	!	time_diff			= 5e-05_dkind
	!	time_stabilization	= 5e-04_dkind			
					
		time_delay			= 1e-05_dkind			
		time_diff			= 1e-04_dkind
		time_stabilization	= 5e-06_dkind			
		
		if ( time > (correction+1)*(time_diff) + time_delay) then			
					
			current_time = time
		
			!# 1D front tracer
			do k = cons_inner_loop(3,1),cons_inner_loop(3,2)
			do j = cons_inner_loop(2,1),cons_inner_loop(2,2)
				max_val = 0.0_dkind
				do i = cons_inner_loop(1,1),cons_inner_loop(1,2)-1
					if(bc%bc_markers(i,j,k) == 0) then	
					
						!! Grad temp
						!if (abs(T%cells(i+1,j,k)-T%cells(i-1,j,k)) > max_val) then
						!	max_val = abs(T%cells(i+1,j,k)-T%cells(i-1,j,k))
						!	flame_front_coords(j) = (i - 0.5_dkind)*cell_size(1) 
						!	flame_front_index = i
						!end if
					
						! max H
						if (abs(Y%pr(HO2_index)%cells(i,j,k)) > max_val) then
							max_val = Y%pr(HO2_index)%cells(i,j,k)
							flame_front_coords(j) = (i - 0.5_dkind)*cell_size(1) 
							flame_front_index = i
						end if					
					
						! max CO
						!if (abs(Y%pr(CO_index)%cells(i,j,k)) > max_val) then
						!	max_val = Y%pr(CO_index)%cells(i,j,k)
						!	flame_front_coords(j) = (i - 0.5_dkind)*cell_size(1) 
						!	flame_front_index = i
						!end if
				
					end if
				end do
			end do
			end do
	
			!left_val	= T%cells(flame_front_index,1,1) - T%cells(flame_front_index-2,1,1)
			!right_val	= T%cells(flame_front_index+2,1,1) - T%cells(flame_front_index,1,1)
			
			!left_val	= Y%pr(CO_index)%cells(flame_front_index-1,1,1)
			!right_val	= Y%pr(CO_index)%cells(flame_front_index+1,1,1)	
			 
			left_val	= Y%pr(HO2_index)%cells(flame_front_index-1,1,1)
			right_val	= Y%pr(HO2_index)%cells(flame_front_index+1,1,1)				 
			
			a = (right_val + left_val - 2.0_dkind * max_val)/2.0_dkind/cell_size(1)**2
			b = (max_val - left_val)/cell_size(1) - a*(2.0_dkind*flame_front_coords(1) - cell_size(1))
			
			current_flame_location = -b/2.0_dkind/a

			if(correction == 0) then
				previous_flame_location = current_flame_location
            end if
			
            boundary = .false.
            if(flame_front_index > cons_inner_loop(1,2) - 10) boundary = .true.
            
			if( (correction /= 0).and.(current_flame_location /=  previous_flame_location) )then 
                
				flame_velocity = (current_flame_location - previous_flame_location)/(current_time - previous_time)

				previous_flame_location = current_flame_location
				previous_time = current_time
				
                previous_av_flame_velocity = sum(flame_velocity_array)
                
                do i = 19, 1, -1
					flame_velocity_array(i+1) = flame_velocity_array(i)  
                end do
                flame_velocity_array(1) = flame_velocity
                
                av_flame_velocity = sum(flame_velocity_array)
                
				if( (correction /= 0).and.(abs(av_flame_velocity - previous_av_flame_velocity) < 1e-03))then 
					counter = counter + 1
                end if
                
                write (flame_loc_unit,'(5E14.6)') time, current_flame_location, flame_velocity, av_flame_velocity, abs(av_flame_velocity - previous_av_flame_velocity)

			end if

			if ((counter > 10).or.boundary) then
				stabilized = .true.
			end if
			
			correction = correction + 1
			
		end if	
			
		end associate
    end subroutine	    
    
end module
