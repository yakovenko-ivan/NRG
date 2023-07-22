module particles_solver_class

	use kind_parameters
	use global_data
	use field_pointers
	use boundary_conditions_class
	use data_manager_class
	use computational_mesh_class
	use computational_domain_class
	use thermophysical_properties_class
	use chemical_properties_class
	use solver_options_class
	
	use mpi_communications_class

	implicit none

#ifdef OMP	
	include "omp_lib.h"
#endif

	private
	public	:: particles_solver, particles_solver_c

	type(field_scalar_cons)	,dimension(:)	,allocatable ,target	:: T_p, T_p_int, E_f_prod_p, rho_p, alpha_p
	type(field_scalar_flow)	,dimension(:)	,allocatable ,target	:: m_flux_p
	type(field_vector_cons)	,dimension(:)	,allocatable ,target	:: v_prod_p, v_p, v_p_int	

	type	:: particles_solver
		type(field_scalar_cons_pointer)		:: T_p, T_p_int, E_f_prod, rho_p, alpha_p, T, rho, nu, kappa
		type(field_scalar_flow_pointer)		:: m_flux_p
		type(field_vector_cons_pointer)		:: v_prod, v_p, v_p_int, v
		type(computational_domain)			:: domain
		type(mpi_communications)			:: mpi_support
		type(boundary_conditions_pointer)	:: boundary
		type(computational_mesh_pointer)	:: mesh
		type(chemical_properties_pointer)	:: chem
		
		type(solid_particles_phase)			:: particles_params

		real(dkind)							:: particle_mass
	contains
		procedure				::  set_initial_distributions
		procedure				::	particles_euler_step_v_E
		procedure				::	particles_lagrange_step
		procedure				::	particles_final_step
		procedure				::  apply_boundary_conditions_main
		procedure				::  apply_boundary_conditions_interm_v_p
		procedure				::	pre_constructor
	end type

	interface	particles_solver_c
		module procedure	constructor
	end interface

contains

	subroutine pre_constructor(this,number_of_phases)
		class(particles_solver)	,intent(inout)	:: this	
		integer					,intent(in)		:: number_of_phases
		
		allocate(	T_p(number_of_phases), T_p_int(number_of_phases), rho_p(number_of_phases), E_f_prod_p(number_of_phases), &
					m_flux_p(number_of_phases), alpha_p(number_of_phases), v_p(number_of_phases), v_p_int(number_of_phases), v_prod_p(number_of_phases))	
	
	end subroutine

	type(particles_solver)	function constructor(manager,solid_particles,phase_number)

		type(data_manager)			, intent(inout)	:: manager
		type(solid_particles_phase)	, intent(in)	:: solid_particles
		integer						, intent(in)	:: phase_number

		type(field_scalar_cons_pointer)	:: scal_ptr
		type(field_vector_cons_pointer)	:: vect_ptr
		type(field_tensor_cons_pointer)	:: tens_ptr
		
 		integer					:: dimensions       
 		integer	,dimension(3,2)	:: cons_allocation_bounds, flow_allocation_bounds   
        
		character(len=40)		:: var_name, var_short_name
		
		cons_allocation_bounds		= manager%domain%get_local_utter_cells_bounds()
		flow_allocation_bounds		= manager%domain%get_local_utter_faces_bounds()
		dimensions					= manager%domain%get_domain_dimensions()   		
		
		constructor%particles_params	= solid_particles
		
		call manager%get_cons_field_pointer_by_name(scal_ptr,vect_ptr,tens_ptr,'density')
		constructor%rho%s_ptr			=> scal_ptr%s_ptr
		call manager%get_cons_field_pointer_by_name(scal_ptr,vect_ptr,tens_ptr,'temperature')
		constructor%T%s_ptr			=> scal_ptr%s_ptr
		call manager%get_cons_field_pointer_by_name(scal_ptr,vect_ptr,tens_ptr,'viscosity')
		constructor%nu%s_ptr			=> scal_ptr%s_ptr
		call manager%get_cons_field_pointer_by_name(scal_ptr,vect_ptr,tens_ptr,'thermal_conductivity')
		constructor%kappa%s_ptr			=> scal_ptr%s_ptr
		
		call manager%get_cons_field_pointer_by_name(scal_ptr,vect_ptr,tens_ptr,'velocity')
		constructor%v%v_ptr				=> vect_ptr%v_ptr		
		
		write(var_name,'(A,I2.2)')		'temperature_particles', phase_number
		write(var_short_name,'(A,I2.2)')	'T_p', phase_number
		call manager%create_scalar_field(T_p(phase_number),	var_name,	var_short_name)
		constructor%T_p%s_ptr		=> T_p(phase_number)
		
		write(var_name,'(A,I2.2)')		'temperature_particles_interm', phase_number
		write(var_short_name,'(A,I2.2)')	'T_p_int', phase_number
		call manager%create_scalar_field(T_p_int(phase_number),	var_name,	var_short_name)
		constructor%T_p_int%s_ptr		=> T_p_int(phase_number)
		
		write(var_name,'(A,I2.2)')		'density_particles', phase_number
		write(var_short_name,'(A,I2.2)')	'rho_p', phase_number
		call manager%create_scalar_field(rho_p(phase_number),	var_name,	var_short_name)
		constructor%rho_p%s_ptr		=> rho_p(phase_number)
		
		write(var_name,'(A,I2.2)')		'volume_fraction_particles', phase_number
		write(var_short_name,'(A,I2.2)')	'alpha_p', phase_number
		call manager%create_scalar_field(alpha_p(phase_number),	var_name,	var_short_name)
		constructor%alpha_p%s_ptr		=> alpha_p(phase_number)
				
		write(var_name,'(A,I2.2)')		'velocity_particles', phase_number
		write(var_short_name,'(A,I2.2)')	'v_p', phase_number		
		call manager%create_vector_field(v_p(phase_number),	var_name,	var_short_name,	'spatial')
		constructor%v_p%v_ptr		=> v_p(phase_number)		
		
		write(var_name,'(A,I2.2)')		'velocity_particles_interm', phase_number
		write(var_short_name,'(A,I2.2)')	'v_p_int', phase_number		
		call manager%create_vector_field(v_p_int(phase_number),	var_name,	var_short_name,	'spatial')
		constructor%v_p_int%v_ptr		=> v_p_int(phase_number)			
		
		write(var_name,'(A,I2.2)')		'mass_flux_particles', phase_number
		write(var_short_name,'(A,I2.2)')	'm_flux_p', phase_number		
		call manager%create_scalar_field(m_flux_p(phase_number),		var_name,	var_short_name)
		constructor%m_flux_p%s_ptr		=> m_flux_p(phase_number)		
		
		write(var_name,'(A,I2.2)')		'energy_production_particles', phase_number
		write(var_short_name,'(A,I2.2)')	'E_f_prod_p', phase_number
		call manager%create_scalar_field(E_f_prod_p(phase_number),	var_name,	var_short_name)
		constructor%E_f_prod%s_ptr		=> E_f_prod_p(phase_number)

		write(var_name,'(A,I2.2)')		'velocity_production_particles', phase_number
		write(var_short_name,'(A,I2.2)')	'v_prod_p', phase_number		
		call manager%create_vector_field(v_prod_p(phase_number),	var_name,	var_short_name,	'spatial')
		constructor%v_prod%v_ptr		=> v_prod_p(phase_number)
		
		constructor%rho_p%s_ptr%cells	= constructor%particles_params%material_density * constructor%alpha_p%s_ptr%cells
		
		constructor%mesh%mesh_ptr	=> manager%computational_mesh_pointer%mesh_ptr
		constructor%boundary%bc_ptr => manager%boundary_conditions_pointer%bc_ptr
		constructor%chem%chem_ptr	=> manager%chemistry%chem_ptr
		constructor%domain			= manager%domain
		constructor%mpi_support		= manager%mpi_communications
		
	end function

	subroutine set_initial_distributions(this)
		class(particles_solver)	,intent(inout)	:: this

		integer		,dimension(3,2)	:: cons_inner_loop

		integer :: i,j,k

		cons_inner_loop	= this%domain%get_local_inner_cells_bounds()		

		associate(  alpha_p		=> this%alpha_p%s_ptr		, &
					rho_p		=> this%rho_p%s_ptr			, &	
					T_p			=> this%T_p%s_ptr			, &
					T			=> this%T%s_ptr				, &
					particle	=> this%particles_params	, &
					mesh		=> this%mesh%mesh_ptr		, &
					bc			=> this%boundary%bc_ptr)	
		
					
		do k = cons_inner_loop(3,1),cons_inner_loop(3,2)
		do j = cons_inner_loop(2,1),cons_inner_loop(2,2)
		do i = cons_inner_loop(1,1),cons_inner_loop(1,2)	
			if(bc%bc_markers(i,j,k) == 0) then
				rho_p%cells(i,j,k)	= alpha_p%cells(i,j,k) * particle%material_density
				T_p%cells(i,j,k)	= 300.0_dkind !T%cells(i,j,k)
				this%particle_mass	= 1.0_dkind/6.0_dkind*particle%diameter**3 * particle%material_density
			end if
		end do
		end do
		end do

		end associate		
	end subroutine
	
	
	subroutine particles_euler_step_v_E(this,time_step)

		class(particles_solver)	,intent(inout)	:: this
		real(dkind)				,intent(in)		:: time_step
		real(dkind)								:: F_stokes, Q_stokes, Nusselt
		
		integer		,dimension(3,2)	:: cons_inner_loop
		real(dkind)	,dimension(3)	:: cell_size
		
		integer	:: dimensions
		integer	:: sign
		integer :: i,j,k,plus,dim

		dimensions		= this%domain%get_domain_dimensions()
		
		cons_inner_loop	= this%domain%get_local_inner_cells_bounds()		

		cell_size		= this%mesh%mesh_ptr%get_cell_edges_length()

		Nusselt			= 2.0_dkind
		
		associate(  T			=> this%T%s_ptr			, &
					T_p			=> this%T_p%s_ptr		, &	
					T_p_int		=> this%T_p_int%s_ptr	, &
					rho			=> this%rho%s_ptr		, &
					rho_p		=> this%rho_p%s_ptr		, &	
					v			=> this%v%v_ptr			, &
					v_p			=> this%v_p%v_ptr		, &
					v_p_int		=> this%v_p_int%v_ptr	, &
					kappa		=> this%kappa%s_ptr		, &
					nu			=> this%nu%s_ptr		, &
					E_f_prod	=> this%E_f_prod%s_ptr	, &
					v_prod		=> this%v_prod%v_ptr	, &					
					particle	=> this%particles_params	, &
					mesh		=> this%mesh%mesh_ptr	, &
					bc			=> this%boundary%bc_ptr)

		!$omp parallel default(none)  private(i,j,k,dim,F_stokes,Q_stokes) , &
		!$omp& firstprivate(this)	,&
		!$omp& shared(T,T_p,T_p_int,E_f_prod,rho,rho_p,v_p,v_p_int,v,v_prod,nu,kappa,particle,Nusselt,time_step,mesh,bc,cons_inner_loop,dimensions)
		!$omp do collapse(3) schedule(guided)					
					
		do k = cons_inner_loop(3,1),cons_inner_loop(3,2)
		do j = cons_inner_loop(2,1),cons_inner_loop(2,2)
		do i = cons_inner_loop(1,1),cons_inner_loop(1,2)
			if(bc%bc_markers(i,j,k) == 0) then
				do dim = 1,dimensions
					F_stokes						= 3.0_dkind * Pi * particle%diameter * nu%cells(i,j,k) / this%particle_mass * ( v%pr(dim)%cells(i,j,k) - v_p%pr(dim)%cells(i,j,k))
					v_prod%pr(dim)%cells(i,j,k)		= - F_stokes * rho_p%cells(i,j,k) / rho%cells(i,j,k) * time_step
					v_p_int%pr(dim)%cells(i,j,k)	= v_p%pr(dim)%cells(i,j,k) + F_stokes * time_step

					E_f_prod%cells(i,j,k)			= - F_stokes * v%pr(dim)%cells(i,j,k) * rho_p%cells(i,j,k) / rho%cells(i,j,k) * time_step
				end do		
					
				Q_stokes					= 6.0_dkind * kappa%cells(i,j,k) * Nusselt / (particle%diameter ** 2 * particle%material_heat_capacity * particle%material_density) * (T%cells(i,j,k) - T_p%cells(i,j,k))
				E_f_prod%cells(i,j,k)		= E_f_prod%cells(i,j,k) - Q_stokes * rho_p%cells(i,j,k) / rho%cells(i,j,k) * particle%material_heat_capacity * time_step
				T_p_int%cells(i,j,k)		= T_p%cells(i,j,k) + Q_stokes * time_step
			end if
		end do
		end do
		end do

		!$omp end do nowait
		!$omp end parallel		
		
		end associate
	end subroutine

	
	subroutine particles_lagrange_step(this,time_step)
 
		class(particles_solver)	,intent(inout)	:: this
		real(dkind)				,intent(in)		:: time_step
 
		real(dkind)	:: av_velocity, dif_velocity
		
		integer		,dimension(3,2)	:: flow_inner_loop, cons_inner_loop
		
		character(len=20)	:: coordinate_system
 
		character(len=20)	:: boundary_type_name		
		
		real(dkind)	,dimension(3)	:: cell_size	, cell_surface_area	
		integer	:: dimensions
		integer	:: sign, bound_number
		integer :: i,j,k,plus,dim
 
		dimensions		= this%domain%get_domain_dimensions()
		
		cell_size			= this%mesh%mesh_ptr%get_cell_edges_length()
		coordinate_system	= this%domain%get_coordinate_system_name()
		
		flow_inner_loop	= this%domain%get_local_inner_faces_bounds()
		cons_inner_loop	= this%domain%get_local_inner_cells_bounds()
		
		associate(  m_flux_p	=> this%m_flux_p%s_ptr	, &
					rho_p		=> this%rho_p%s_ptr		, &
					v_p_int		=> this%v_p_int%v_ptr	, &
					mesh		=> this%mesh%mesh_ptr	, &
					bc			=> this%boundary%bc_ptr)
					
		call this%mpi_support%exchange_conservative_vector_field(v_p_int)
		call this%mpi_support%exchange_conservative_scalar_field(rho_p)					
					
		!$omp parallel default(none)  private(i,j,k,dim,plus,sign,bound_number,boundary_type_name,av_velocity,dif_velocity,cell_surface_area) , &
		!$omp& firstprivate(this)	,&
		!$omp& shared(bc,mesh,m_flux_p,rho_p,v_p_int,dimensions,cell_size,time_step,flow_inner_loop,cons_inner_loop,coordinate_system)
		!$omp do collapse(3) schedule(guided)
					
		do k = flow_inner_loop(3,1),flow_inner_loop(3,2)
		do j = flow_inner_loop(2,1),flow_inner_loop(2,2)
		do i = flow_inner_loop(1,1),flow_inner_loop(1,2)
	!		if(bc%bc_markers(i,j,k) == 0) then
                
				cell_surface_area	= this%mesh%mesh_ptr%get_cell_surface_area()
	 
				do dim = 1,dimensions
				
					select case(coordinate_system)
						case ('cartesian')	
							cell_surface_area	= cell_surface_area
						case ('cylindrical')
							! x -> r, y -> z
							if(dim==1) cell_surface_area(dim) = cell_surface_area(dim) * (mesh%mesh(1,i,j,k) - 0.5_dkind*cell_size(1))									
							if(dim==2) cell_surface_area(dim) = cell_surface_area(dim) * (mesh%mesh(1,i,j,k))		! - 0.5_dkind*cell_size(1)
						case ('spherical')
							! x -> r
							if(dim==1) cell_surface_area(dim) = cell_surface_area(dim) * (mesh%mesh(1,i,j,k))**2	! - 0.5_dkind*cell_size(1)
					end select				
		
                    av_velocity     = 0.5_dkind *(v_p_int%pr(dim)%cells(i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3)) + v_p_int%pr(dim)%cells(i,j,k))
                    dif_velocity    = time_step *(v_p_int%pr(dim)%cells(i,j,k) - v_p_int%pr(dim)%cells(i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3)))/ cell_size(dim)
                    if( av_velocity >= 0 ) then
                        m_flux_p%cells(dim,i,j,k) = rho_p%cells(i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3))		* av_velocity * (cell_surface_area(dim) ) * time_step / (1.0_dkind + dif_velocity)
                    else
                        m_flux_p%cells(dim,i,j,k) = rho_p%cells(i,j,k)										* av_velocity * (cell_surface_area(dim) ) * time_step / (1.0_dkind + dif_velocity)
                    end if
					
					if((bc%bc_markers(i,j,k) == 0).and.(i <= cons_inner_loop(1,2)).and.(j <= cons_inner_loop(2,2)).and.(k <= cons_inner_loop(3,2))) then
						do plus = 1,2
							sign			= (-1)**plus
							bound_number	= bc%bc_markers(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3))
							if( bound_number /= 0 ) then
								boundary_type_name = bc%boundary_types(bound_number)%get_type_name()
								select case(boundary_type_name)
									case ('wall')		! Particles stay near top wall
										if ((dim == 2).and.(sign == 1)) then
											m_flux_p%cells(dim,i,j,k) = max(0.0_dkind,m_flux_p%cells(dim,i,j,k))							
										end if
								end select
							end if					
						end do
					end if
                end do
	!		end if
		end do
		end do
		end do
 
		!$omp end do nowait
		!$omp end parallel
		
		continue
		
		end associate
 
	end subroutine

	subroutine particles_final_step(this,time_step)
 
		class(particles_solver)	,intent(inout)	:: this
		real(dkind)				,intent(in)		:: time_step
 
		real(dkind)	:: D11, D12, D21, D22, rho_p_old, av_velocity1, av_velocity2
		
		integer		,dimension(3,2)	:: cons_inner_loop
		real(dkind)	,dimension(3)	:: cell_size
		real(dkind)					:: cell_volume
		integer	:: dimensions, species_number
		character(len=20)	:: coordinate_system
		integer	:: sign
		integer :: i,j,k,plus,dim,dim1,dim2,spec,i_ind1,i_ind2,j_ind1,j_ind2,k_ind1,k_ind2
 
		dimensions			= this%domain%get_domain_dimensions()
		species_number		= this%chem%chem_ptr%species_number
		coordinate_system	= this%domain%get_coordinate_system_name()
		
		cons_inner_loop	= this%domain%get_local_inner_cells_bounds()
	
		cell_size			= this%mesh%mesh_ptr%get_cell_edges_length()
		
		associate(  m_flux_p	=> this%m_flux_p%s_ptr	, &
					rho_p		=> this%rho_p%s_ptr		, &
					alpha_p		=> this%alpha_p%s_ptr	, &
					T_p			=> this%T_p%s_ptr		, &
					T_p_int		=> this%T_p_int%s_ptr	, &
					v_p			=> this%v_p%v_ptr		, &
					v_p_int		=> this%v_p_int%v_ptr	, &
					mesh		=> this%mesh%mesh_ptr	, &
					bc			=> this%boundary%bc_ptr)
 
		call this%mpi_support%exchange_conservative_scalar_field(T_p_int)

		call this%mpi_support%exchange_flow_scalar_field(m_flux_p)					
					
		!$omp parallel default(none)  private(i,j,k,dim,dim1,dim2,spec,rho_p_old,av_velocity1,av_velocity2,cell_volume,D11,D21,D12,D22,i_ind1,i_ind2,j_ind1,j_ind2,k_ind1,k_ind2) , &
		!$omp& firstprivate(this)	,&
		!$omp& shared(bc,mesh,m_flux_p,rho_p,T_p_int,T_p,v_p,v_p_int,alpha_p,dimensions,cons_inner_loop,coordinate_system)
		!$omp do collapse(3) schedule(guided)
					
		do k = cons_inner_loop(3,1),cons_inner_loop(3,2)
		do j = cons_inner_loop(2,1),cons_inner_loop(2,2)
		do i = cons_inner_loop(1,1),cons_inner_loop(1,2)
			if(bc%bc_markers(i,j,k) == 0) then
			
				cell_volume	= this%mesh%mesh_ptr%get_cell_volume()
				select case(coordinate_system)
					case ('cartesian')
						cell_volume			= cell_volume
					case ('cylindrical')
						cell_volume			= cell_volume * mesh%mesh(1,i,j,k)
					case ('spherical')
						cell_volume			= cell_volume * mesh%mesh(1,i,j,k)**2
				end select				
 
				rho_p_old			= rho_p%cells(i,j,k)
				do dim = 1,dimensions
					rho_p%cells(i,j,k) = rho_p%cells(i,j,k) + (m_flux_p%cells(dim,i,j,k)-m_flux_p%cells(dim,i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3))) / (cell_volume)
				end do
 
				T_p%cells(i,j,k)		= T_p_int%cells(i,j,k)  * rho_p_old / rho_p%cells(i,j,k)
				do dim = 1,dimensions
					v_p%pr(dim)%cells(i,j,k)	= v_p_int%pr(dim)%cells(i,j,k) * rho_p_old / rho_p%cells(i,j,k)
				end do

 
                do dim2 = 1,dimensions
                    D11 = 0.0_dkind
					D12 = 0.0_dkind
					D21 = 0.0_dkind
					D22 = 0.0_dkind
 
					i_ind1 = i - I_m(dim2,1)
					i_ind2 = i + I_m(dim2,1)
					j_ind1 = j - I_m(dim2,2)
					j_ind2 = j + I_m(dim2,2)
					k_ind1 = k - I_m(dim2,3)
					k_ind2 = k + I_m(dim2,3)
 
                    av_velocity1 = v_p_int%pr(dim2)%cells(i_ind1,j_ind1,k_ind1) + v_p_int%pr(dim2)%cells(i,j,k)
                    av_velocity2 = v_p_int%pr(dim2)%cells(i_ind2,j_ind2,k_ind2) + v_p_int%pr(dim2)%cells(i,j,k)
 
                    if(av_velocity1 > 0.0) then
						D11 = 1.0_dkind
					else
						D12 = 1.0_dkind
					end if
                    if(av_velocity2 < 0.0) then
						D21 = 1.0_dkind
					else
						D22 = 1.0_dkind
					end if
 
 
					T_p%cells(i,j,k) = T_p%cells(i,j,k)		+ (D11 * T_p_int%cells(i_ind1,j_ind1,k_ind1) * abs(m_flux_p%cells(dim2,i,j,k)) &
															+  D21 * T_p_int%cells(i_ind2,j_ind2,k_ind2) * abs(m_flux_p%cells(dim2,i_ind2,j_ind2,k_ind2)) &
															-  D12 * T_p_int%cells(i,j,k)				 * abs(m_flux_p%cells(dim2,i,j,k)) &
															-  D22 * T_p_int%cells(i,j,k)				 * abs(m_flux_p%cells(dim2,i_ind2,j_ind2,k_ind2))) /  rho_p%cells(i,j,k)  /(cell_volume)
 
					do dim1 = 1,dimensions
						v_p%pr(dim1)%cells(i,j,k) = v_p%pr(dim1)%cells(i,j,k)	+ (D11 * v_p_int%pr(dim1)%cells(i_ind1,j_ind1,k_ind1)	* abs(m_flux_p%cells(dim2,i,j,k)) &
																				+  D21 * v_p_int%pr(dim1)%cells(i_ind2,j_ind2,k_ind2)	* abs(m_flux_p%cells(dim2,i_ind2,j_ind2,k_ind2)) &
																				-  D12 * v_p_int%pr(dim1)%cells(i,j,k)				* abs(m_flux_p%cells(dim2,i,j,k)) &
																				-  D22 * v_p_int%pr(dim1)%cells(i,j,k)				* abs(m_flux_p%cells(dim2,i_ind2,j_ind2,k_ind2))) / rho_p%cells(i,j,k) / (cell_volume)
					end do
 
					alpha_p%cells(i,j,k) = alpha_p%cells(i,j,k) * rho_p%cells(i,j,k) / rho_p_old
					
                end do
			end if
		end do
		end do
		end do
 
		!$omp end do nowait
		!$omp end parallel
		
		end associate
	end subroutine

	subroutine apply_boundary_conditions_interm_v_p(this)

		class(particles_solver)		,intent(inout)		:: this

		integer	:: dimensions
		integer	,dimension(3,2)	:: cons_inner_loop
		
		character(len=20)		:: boundary_type_name		

		integer	:: sign, bound_number
		integer :: i,j,k,plus,dim,dim1,specie_number
								
		dimensions		= this%domain%get_domain_dimensions()
				
		cons_inner_loop	= this%domain%get_local_inner_cells_bounds()			
		
		associate(  v_p_int			=> this%v_p_int%v_ptr		, &
					bc				=> this%boundary%bc_ptr	, &
					mesh			=> this%mesh%mesh_ptr)

		!$omp parallel default(none)  private(i,j,k,plus,dim,dim1,sign,bound_number,boundary_type_name) , &
		!$omp& firstprivate(this)	,&
		!$omp& shared(v_p_int,bc,cons_inner_loop,dimensions)
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

								do dim1 = 1, dimensions
									if(dim1 == dim) then
										v_p_int%pr(dim1)%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3)) = - v_p_int%pr(dim1)%cells(i,j,k)
									else
										v_p_int%pr(dim1)%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3)) = v_p_int%pr(dim1)%cells(i,j,k)
									end if
								end do
								
								boundary_type_name = bc%boundary_types(bound_number)%get_type_name()
								select case(boundary_type_name)
									case ('outlet')
										do dim1 = 1, dimensions
											if(dim1 == dim)	v_p_int%pr(dim1)%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3)) = v_p_int%pr(dim1)%cells(i,j,k)
										end do
										
									case ('wall')		! Particles stay near top wall
										if ((dim == 2).and.(sign == 1)) then
											do dim1 = 1, dimensions
												v_p_int%pr(dim1)%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3)) = 0.0_dkind
												v_p_int%pr(dim1)%cells(i,j,k)	= 0.0_dkind
											end do										
										end if
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
	
	subroutine apply_boundary_conditions_main(this)

		class(particles_solver)		,intent(inout)		:: this

		character(len=20)		:: boundary_type_name
		real(dkind)				:: farfield_density, farfield_pressure, wall_temperature
		
		integer					:: dimensions
		integer	,dimension(3,2)	:: cons_inner_loop

		integer	:: sign, bound_number
		integer :: i,j,k,plus,dim,dim1
						
		dimensions		= this%domain%get_domain_dimensions()
			
		cons_inner_loop	= this%domain%get_local_inner_cells_bounds()			
		
		associate(  T_p				=> this%T_p%s_ptr			, &
					rho_p			=> this%rho_p%s_ptr			, &
					v_p				=> this%v_p%v_ptr			, &
					bc				=> this%boundary%bc_ptr		, &
					mesh			=> this%mesh%mesh_ptr)

		!$omp parallel default(none)  private(i,j,k,plus,dim,dim1,sign,bound_number,boundary_type_name,wall_temperature) , &
		!$omp& firstprivate(this)	,&
		!$omp& shared(T_p,rho_p,v_p,mesh,bc,cons_inner_loop,dimensions)
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

								rho_p%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3))	= rho_p%cells(i,j,k)
								T_p%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3))	= T_p%cells(i,j,k)

								do dim1 = 1, dimensions
									if(dim1 == dim) then
										v_p%pr(dim1)%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3)) = - v_p%pr(dim1)%cells(i,j,k)
									else
										v_p%pr(dim1)%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3)) = v_p%pr(dim1)%cells(i,j,k)
									end if
								end do								
								
								boundary_type_name = bc%boundary_types(bound_number)%get_type_name()
								select case(boundary_type_name)
									case('wall')
										if(bc%boundary_types(bound_number)%is_conductive()) then 
											wall_temperature = bc%boundary_types(bound_number)%get_wall_temperature()
											T_p%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3)) = wall_temperature
										end if
										if(.not.bc%boundary_types(bound_number)%is_slip()) then
											do dim1 = 1, dimensions
												if (dim1 /= dim) then	
													v_p%pr(dim1)%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3)) = 0.0_dkind
												end if
											end do
										end if
										if ((dim == 2).and.(sign == 1)) then
											do dim1 = 1, dimensions
												v_p%pr(dim1)%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3)) = 0.0_dkind
												v_p%pr(dim1)%cells(i,j,k) = 0.0_dkind
											end do										
										end if
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
	
	
	
end module
