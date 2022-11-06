module lagrangian_droplets_solver_class

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
	public	:: lagrangian_droplets_solver, lagrangian_droplets_solver_c

	type(field_scalar_cons)	,dimension(:)	,allocatable	,target	:: E_f_prod_d, rho_prod_d
	type(field_vector_cons)	,dimension(:)	,allocatable	,target	:: Y_prod_d
	type(field_vector_flow)	,dimension(:)	,allocatable	,target	:: v_prod_d

	integer	:: average_io_unit
	
	type	:: lagrangian_droplet
		real(dkind)	,dimension(3)	:: coords, velocity
		integer		,dimension(3)	:: cell
		real(dkind)					:: temperature
		real(dkind)					:: mass	
		real(dkind)					:: dm
		logical						:: outside_domain
	end type
	
	type	:: lagrangian_droplets_solver
		type(field_scalar_cons_pointer)		:: T, rho, nu, kappa, E_f_prod, rho_prod, p, mol_mix_conc, h_s
		type(field_vector_cons_pointer)		:: Y, Y_prod, D
		type(field_vector_flow_pointer)		:: v_f, v_prod
		type(computational_domain)			:: domain

		type(boundary_conditions_pointer)	:: boundary
		type(computational_mesh_pointer)	:: mesh
		type(chemical_properties_pointer)	:: chem
		type(thermophysical_properties_pointer)		:: thermo
		
		type(liquid_droplets_phase)			:: droplets_params

		real(dkind)							:: droplet_mass
        		
		type(lagrangian_droplet)	,dimension(:)	,allocatable	:: droplets	
		
		integer								:: droplets_number
        integer                             :: phase_number
	contains
		procedure				::  set_initial_distributions
		procedure				::	droplets_solve
		procedure				::  apply_boundary_conditions_main
		procedure				::  apply_boundary_conditions_interm_v_d
		procedure				::	pre_constructor
		
		procedure	,private	:: release_droplet
		procedure	,private	:: get_droplet_cell
		procedure	,private	:: count_droplets_in_cell
	end type

	interface	lagrangian_droplets_solver_c
		module procedure	constructor
	end interface

contains

	subroutine pre_constructor(this,number_of_phases)
		class(lagrangian_droplets_solver)	,intent(inout)	:: this	
		integer					,intent(in)		:: number_of_phases
		
		allocate(	E_f_prod_d(number_of_phases), rho_prod_d(number_of_phases), &
					v_prod_d(number_of_phases), Y_prod_d(number_of_phases))	
	
	end subroutine

	type(lagrangian_droplets_solver)	function constructor(manager,liquid_droplets,phase_number)

		type(data_manager)			, intent(inout)	:: manager
		type(liquid_droplets_phase)	, intent(in)	:: liquid_droplets
		integer						, intent(in)	:: phase_number

		type(field_scalar_cons_pointer)	:: scal_ptr
		type(field_vector_cons_pointer)	:: vect_ptr
		type(field_tensor_cons_pointer)	:: tens_ptr

		type(field_scalar_flow_pointer)	:: scal_f_ptr
		type(field_vector_flow_pointer)	:: vect_f_ptr			
		
 		integer					:: dimensions       
 		integer	,dimension(3,2)	:: cons_allocation_bounds, flow_allocation_bounds   
        
		character(len=100)		:: system_command
		character(len=40)		:: var_name, var_short_name
		character(len=100)		:: file_name
		
		integer	:: drop, dim
		
		cons_allocation_bounds		= manager%domain%get_local_utter_cells_bounds()
		flow_allocation_bounds		= manager%domain%get_local_utter_faces_bounds()
		dimensions					= manager%domain%get_domain_dimensions()   		
		
		constructor%droplets_params	= liquid_droplets
		
		call manager%get_cons_field_pointer_by_name(scal_ptr,vect_ptr,tens_ptr,'density')
		constructor%rho%s_ptr				=> scal_ptr%s_ptr
		call manager%get_cons_field_pointer_by_name(scal_ptr,vect_ptr,tens_ptr,'temperature')
		constructor%T%s_ptr					=> scal_ptr%s_ptr
		call manager%get_cons_field_pointer_by_name(scal_ptr,vect_ptr,tens_ptr,'pressure')
		constructor%p%s_ptr					=> scal_ptr%s_ptr
		call manager%get_cons_field_pointer_by_name(scal_ptr,vect_ptr,tens_ptr,'specie_molar_concentration')
		constructor%Y%v_ptr					=> vect_ptr%v_ptr		
		call manager%get_cons_field_pointer_by_name(scal_ptr,vect_ptr,tens_ptr,'mixture_molar_concentration')
		constructor%mol_mix_conc%s_ptr		=> scal_ptr%s_ptr		
		call manager%get_cons_field_pointer_by_name(scal_ptr,vect_ptr,tens_ptr,'sensible_enthalpy')
		constructor%h_s%s_ptr				=> scal_ptr%s_ptr		
		
		call manager%get_cons_field_pointer_by_name(scal_ptr,vect_ptr,tens_ptr,'viscosity')
		constructor%nu%s_ptr			=> scal_ptr%s_ptr
		call manager%get_cons_field_pointer_by_name(scal_ptr,vect_ptr,tens_ptr,'thermal_conductivity')
		constructor%kappa%s_ptr			=> scal_ptr%s_ptr
		call manager%get_cons_field_pointer_by_name(scal_ptr,vect_ptr,tens_ptr,'diffusivity')
		constructor%D%v_ptr				=> vect_ptr%v_ptr		
		
		
		call manager%get_flow_field_pointer_by_name(scal_f_ptr,vect_f_ptr,'velocity_flow')
		constructor%v_f%v_ptr			=> vect_f_ptr%v_ptr		

		
		write(var_name,'(A,I2.2)')		'density_production_droplets', phase_number
		write(var_short_name,'(A,I2.2)')	'rho_prod_d', phase_number
		call manager%create_scalar_field(rho_prod_d(phase_number),	var_name,	var_short_name)
		constructor%rho_prod%s_ptr		=> rho_prod_d(phase_number)        
        
		write(var_name,'(A,I2.2)')		'energy_production_droplets', phase_number
		write(var_short_name,'(A,I2.2)')	'E_f_prod_d', phase_number
		call manager%create_scalar_field(E_f_prod_d(phase_number),	var_name,	var_short_name)
		constructor%E_f_prod%s_ptr		=> E_f_prod_d(phase_number)

		write(var_name,'(A,I2.2)')		'velocity_production_droplets', phase_number
		write(var_short_name,'(A,I2.2)')	'v_prod_d', phase_number		
		call manager%create_vector_field(v_prod_d(phase_number),	var_name,	var_short_name,	'spatial')
		constructor%v_prod%v_ptr		=> v_prod_d(phase_number)
        
		write(var_name,'(A,I2.2)')		'concentration_production_droplets', phase_number
		write(var_short_name,'(A,I2.2)')	'Y_prod_d', phase_number		
		call manager%create_vector_field(Y_prod_d(phase_number),	var_name,	var_short_name,	'chemical')
		constructor%Y_prod%v_ptr		=> Y_prod_d(phase_number)        

		constructor%phase_number    = phase_number

		constructor%mesh%mesh_ptr	=> manager%computational_mesh_pointer%mesh_ptr
		constructor%boundary%bc_ptr => manager%boundary_conditions_pointer%bc_ptr
		constructor%thermo%thermo_ptr	=> manager%thermophysics%thermo_ptr
		constructor%chem%chem_ptr	=> manager%chemistry%chem_ptr
		constructor%domain			= manager%domain
		
		system_command = 'mkdir data_save_droplets'
		call system(system_command)		
		
		file_name = 'data_save_droplets' // trim(fold_sep) // 'average_droplets_data.dat'
		open(newunit = average_io_unit, file = file_name, status = 'replace', form = 'formatted')

	end function

	subroutine set_initial_distributions(this)
		class(lagrangian_droplets_solver)	,intent(inout)	:: this

		integer		,dimension(3,2)	:: cons_inner_loop
		real(dkind)	,dimension(3)	:: cell_size
		
		real(dkind)	:: delta, drop_number_max, a ,b, dimless_length
		integer		,dimension(3)	:: drop_number
		real(dkind)	,dimension(3)	:: coords
		real(dkind)	,dimension(:,:)	,allocatable	:: lengths
		
		integer	:: dimensions
		integer :: i,j,k, drop, drop_x, drop_y, drop_z, dim, droplet

		integer, dimension(3)	:: cell, initial_cell
		logical	:: out_flag
		
		dimensions		= this%domain%get_domain_dimensions()
		cons_inner_loop	= this%domain%get_local_inner_cells_bounds()		

		allocate(lengths(dimensions,2))
		lengths			= this%domain%get_domain_lengths()
		
		cell_size		= this%mesh%mesh_ptr%get_cell_edges_length()
		
!		delta = 0.0009_dkind

!		lengths(1,1) = 0.01_dkind
!		lengths(1,2) = 0.02_dkind
		
		drop_number = 1
		!do dim = 1,dimensions
		!	drop_number(dim) =  (lengths(dim,2)-lengths(dim,1)-cell_size(1)) /delta
		!end do
		
		drop_number(1)  = 900
		drop_number(2)  = 15
        
        
		
!		delta = (lengths(1,2)-lengths(1,1) - 3 * cell_size(1)) / (drop_number(1) - 1) 
		delta = (lengths(2,2)-lengths(2,1) - 3 * cell_size(1)) / (drop_number(2) - 1) 
        
		this%droplets_number = drop_number(1)*drop_number(2)*drop_number(3)
		
		allocate(this%droplets(this%droplets_number))	
		
		do drop = 1, this%droplets_number
			this%droplets(drop)%outside_domain	= .true.
			this%droplets(drop)%coords			= 0.0_dkind
			this%droplets(drop)%velocity		= 0.0_dkind
			this%droplets(drop)%dm				= 0.0_dkind 
		end do 

		do drop_x = 1, drop_number(1)
		do drop_y = 1, drop_number(2)
		do drop_z = 1, drop_number(3)
		
			droplet = drop_x + (drop_y-1)*drop_number(1) + (drop_z-1)*drop_number(1)*drop_number(2)
		
			this%droplets(droplet)%coords			= 0.0_dkind
			this%droplets(droplet)%velocity			= 0.0_dkind
			this%droplets(droplet)%dm				= 0.0_dkind  
			this%droplets(droplet)%outside_domain = .false.
				
			this%droplets(droplet)%temperature	= 300.0_dkind
			this%droplets(droplet)%mass			= Pi*this%droplets_params%diameter**3 / 6.0_dkind * this%droplets_params%material_density
			do dim = 1, dimensions
				dimless_length = real((drop_x-1)*I_m(dim,1)+(drop_y-1)*I_m(dim,2)+(drop_z-1)*I_m(dim,3),dkind)
			
				this%droplets(droplet)%velocity(dim)	= 0.0_dkind
				this%droplets(droplet)%coords(dim)		= lengths(dim,1) + cell_size(1)*1.5 + dimless_length*delta
			end do
			
			initial_cell = this%get_droplet_cell(this%droplets(droplet)%coords,out_flag)
			this%droplets(droplet)%cell				= initial_cell
			this%droplets(droplet)%outside_domain	= out_flag				
		end do
		end do
		end do
				
		continue
		
	end subroutine
	
	
	subroutine droplets_solve(this,time_step)

		class(lagrangian_droplets_solver)	,intent(inout)	:: this
		real(dkind)				,intent(in)		:: time_step
		real(dkind)								:: droplet_time_step
		integer									:: droplet_iterations
		
		integer		,dimension(3,2)	:: cons_inner_loop
		real(dkind)	,dimension(3)	:: cell_size
		real(dkind)					:: cell_volume
		real(dkind)	,dimension(3)	:: gas_velocity, relative_velocity, old_velocity, F_a
		real(dkind)					:: abs_relative_velocity
		real(dkind)					:: lower_face, higher_face, a, b
		real(dkind)					:: Re_p, Nu_p, Sc_p, Sh_p, Pr_p, C_drag, A_p, alpha_p, beta_p, D_air_water
		real(dkind)					:: H_m, H_v, H_l, H_h
		real(dkind)					:: X_vap, Y_vap, W_ratio, dY_dT
		real(dkind)					:: specie_enthalpy_gas, specie_enthalpy_liquid, mixture_cp
		real(dkind)					:: dmp, Q,	Tg_new, Tp_new, Tg_old, Tp_old, rhog_old, rhog_new, M_gas_old, M_gas_new, m_liq, mol_mix_w, Hg_old, Hg_new, Hp_old, Hp_new
		real(dkind)					:: delta_H
		real(dkind)	,dimension(:)	, allocatable	:: Yg_old, Yg_new
		real(dkind)					:: DTOP, DTOG, DTGOG, DTGOP, AGHRHO, DADYDTHVHL, DADYDTHV
		real(dkind)	,dimension(2)	:: A_col, B_col, D_vec
		
		integer						:: droplets_in_cell
		character(len=100)			:: file_path, file_name
		
		real(dkind)                 :: local_diameter, evaporation_rate, velocity_2, specie_enthalpy
		real(dkind)					:: average_diameter, average_temperature, droplets_inside
				
		integer	:: dimensions, species_number
		integer, dimension(3)	:: cell, initial_cell
		
        integer :: droplet_material_index
		integer :: C7H16_index
		integer :: i,j,k,dim, drop, drop_iter, spec
		logical	:: out_flag
		
		real(dkind)	,save	:: time = 0.0_dkind
		integer		,save	:: output_counter = 0
		integer		,save	:: droplet_release_counter = 0

		integer				:: lagrangian_droplets_io_unit
		
		dimensions		= this%domain%get_domain_dimensions()
		species_number	= this%chem%chem_ptr%species_number
		
		allocate(Yg_old(species_number),Yg_new(species_number))
		
		cons_inner_loop	= this%domain%get_local_inner_cells_bounds()		

		cell_size		= this%mesh%mesh_ptr%get_cell_edges_length()
		cell_volume		= this%mesh%mesh_ptr%get_cell_volume()
		
		droplet_material_index		= this%chem%chem_ptr%get_chemical_specie_index(this%droplets_params%material)
		
		associate(  T			=> this%T%s_ptr			, &
					rho			=> this%rho%s_ptr		, &
					p			=> this%p%s_ptr			, &
					h_s			=> this%h_s%s_ptr		, &
					Y			=> this%Y%v_ptr			, &
                    v_f			=> this%v_f%v_ptr		, &
					D			=> this%D%v_ptr			, &
					kappa		=> this%kappa%s_ptr		, &
					nu			=> this%nu%s_ptr		, &
                    rho_prod    => this%rho_prod%s_ptr	, &
					E_f_prod	=> this%E_f_prod%s_ptr	, &
					v_prod		=> this%v_prod%v_ptr	, &
                    Y_prod		=> this%Y_prod%v_ptr	, &		
					droplet	    => this%droplets_params	, &
                    chem        => this%chem%chem_ptr   , &
					mesh		=> this%mesh%mesh_ptr	, &
					bc			=> this%boundary%bc_ptr)

		!!$omp parallel default(none)  private(i,j,k,dim,F_stokes,Q_stokes,local_diameter,evaporation_rate,velabs,temp_cr) , &
		!!$omp& firstprivate(this)	,&
		!!$omp& shared(T,T_d,T_d_int,E_f_prod,rho,rho_d,mass_d,numdens_d,v_d,v_d_int,v,v_prod,Y_prod,nu,kappa,foam_marker,time_boil,droplet,Nusselt,droplet_material_index,time_step,mesh,bc,cons_inner_loop,dimensions,marker1,specie_enthalpy)
		!!$omp do collapse(3) schedule(guided)					
			
		E_f_prod%cells = 0.0_dkind			
		Y_prod%pr(droplet_material_index)%cells = 0.0_dkind	
		
		do dim = 1, dimensions
			v_prod%pr(dim)%cells = 0.0_dkind
		end do

		average_temperature = 0.0_dkind
		average_diameter	= 0.0_dkind
		droplets_inside		= 0
		
		do drop = 1, this%droplets_number

			if (.not.this%droplets(drop)%outside_domain) then
		
				average_diameter	= average_diameter		+ (6.0_dkind * this%droplets(drop)%mass / Pi /  droplet%material_density) ** (1.0_dkind / 3.0_dkind)
				average_temperature	= average_temperature	+ this%droplets(drop)%temperature
				droplets_inside		= droplets_inside + 1
				
				!# Move droplet, account momentum transfer

				droplet_time_step = time_step
				do dim = 1, dimensions
					if (abs(this%droplets(drop)%velocity(dim)) > 1e-10_dkind) then
						if (abs(cell_size(dim) / this%droplets(drop)%velocity(dim)) < droplet_time_step) droplet_time_step = abs(cell_size(dim) / this%droplets(drop)%velocity(dim))
					end if
				end do
				droplet_iterations	= ceiling(time_step/ droplet_time_step)
				droplet_time_step	= time_step / real(droplet_iterations,dkind)

				F_a =	0.0_dkind			
			
				initial_cell = this%get_droplet_cell(this%droplets(drop)%coords,out_flag)
				i = initial_cell(1)
				j = initial_cell(2)
				k = initial_cell(3)
				this%droplets(drop)%cell = initial_cell

				Tg_old		= T%cells(i,j,k)
				rhog_old	= rho%cells(i,j,k)
				
				do spec = 1, species_number
					Yg_old(spec)	= Y%pr(spec)%cells(i,j,k)
				end do	
				M_gas_old	= rhog_old * cell_volume 
				Hg_old		= h_s%cells(i,j,k) * M_gas_old
			
				do drop_iter = 1, droplet_iterations
			
					Tp_old	= this%droplets(drop)%temperature
					Hp_old	= this%droplets(drop)%temperature * droplet%material_heat_capacity * this%droplets(drop)%mass
			
					local_diameter = (6.0_dkind * this%droplets(drop)%mass / Pi /  droplet%material_density) ** (1.0_dkind / 3.0_dkind)
				
					cell = this%get_droplet_cell(this%droplets(drop)%coords,out_flag)
					i = cell(1)
					j = cell(2)
					k = cell(3)				
				
					if ((this%droplets(drop)%cell(1) /= cell(1)).and.(this%droplets(drop)%cell(2) /= cell(2)).and.(this%droplets(drop)%cell(3) /= cell(3))) then
						Tg_old		= T%cells(i,j,k)
						rhog_old	= rho%cells(i,j,k)
						M_gas_old	= rhog_old * cell_volume 
						Hg_old		= h_s%cells(i,j,k) * M_gas_old
						do spec = 1, species_number
							Yg_old(spec)	= Y%pr(spec)%cells(i,j,k)
						end do
					
						E_f_prod%cells(i,j,k) = 0.0_dkind
						Y_prod%pr(droplet_material_index)%cells(i,j,k) = 0.0_dkind					
					
						this%droplets(drop)%cell = cell
					end if

					old_velocity = this%droplets(drop)%velocity
			
				!# Interpolate gas velocity to the droplet location using linear interpolation
					gas_velocity = 0.0_dkind
					do dim = 1, dimensions
						lower_face	= mesh%mesh(dim,i,j,k) - 0.5_dkind*cell_size(dim)
						higher_face	= mesh%mesh(dim,i,j,k) + 0.5_dkind*cell_size(dim)
					
						b = (v_f%pr(dim)%cells(dim,i,j,k)*higher_face - v_f%pr(dim)%cells(dim,i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3))*lower_face) / (higher_face - lower_face) 
						a = (v_f%pr(dim)%cells(dim,i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3)) - v_f%pr(dim)%cells(dim,i,j,k)) / (higher_face - lower_face) 
					
						gas_velocity(dim) = a * this%droplets(drop)%coords(dim) + b
					end do
				
					abs_relative_velocity = 0.0_dkind
					do dim = 1, dimensions
						relative_velocity(dim)	= this%droplets(drop)%velocity(dim) - gas_velocity(dim)
						abs_relative_velocity	= abs_relative_velocity + relative_velocity(dim)*relative_velocity(dim)
					end do
				
					abs_relative_velocity = sqrt(abs_relative_velocity)
				
					!# Calculate droplet mass and cross section
					m_liq	= Pi*local_diameter**3 / 6.0_dkind * droplet%material_density
					A_p		= Pi*local_diameter**2
					!# Calculate droplet Reynolds number
					Re_p	= rhog_old * abs_relative_velocity * local_diameter / nu%cells(i,j,k)
				
					if ( abs_relative_velocity > 1.0e-10_dkind) then	
					!# Calculate drag coefficient for spherical droplet
						if (Re_p < 1e-010) then
							C_drag = 100
						elseif (Re_p < 1.0_dkind) then
								C_drag = 24.0_dkind / Re_p
						elseif (Re_p < 1000.0_dkind) then
								C_drag = 24.0_dkind * ( 0.85_dkind + 0.15 * Re_p**0.687) / Re_p
						elseif (Re_p >= 1000.0_dkind) then
								C_drag = 0.44_dkind
						end if
				
						droplets_in_cell = this%count_droplets_in_cell(cell)
				
					!# Calculate new droplet coordinates and velocities
						alpha_p	= M_gas_old / m_liq
						beta_p	= 0.5_dkind * rhog_old * C_drag * A_p * ( 1/m_liq + 1/M_gas_old) * abs_relative_velocity
					
						do dim = 1, dimensions
							this%droplets(drop)%coords(dim) =	this%droplets(drop)%coords(dim) + (this%droplets(drop)%velocity(dim) + alpha_p * gas_velocity(dim)) * droplet_time_step / (1.0_dkind + alpha_p) + &
																alpha_p * log(1.0_dkind + beta_p*droplet_time_step) / beta_p / (1.0_dkind + alpha_p) * relative_velocity(dim)
														
							this%droplets(drop)%velocity(dim)	= (this%droplets(drop)%velocity(dim) + (this%droplets(drop)%velocity(dim) + alpha_p*gas_velocity(dim)) * beta_p * droplet_time_step / (1.0_dkind + alpha_p)) / (1.0_dkind + beta_p*droplet_time_step)	
					 
!							droplet_acceleration(dim) = droplet_acceleration(dim) + 1.0_dkind/(cell_volume * rhog_old) * (m_liq*(old_velocity(dim) - this%droplets(drop)%velocity(dim))/time_step + this%droplets(drop)%dm * relative_velocity(dim)) 
                            F_a(dim) = F_a(dim) + 1.0_dkind/(cell_volume * rhog_old) * (m_liq*(old_velocity(dim) - this%droplets(drop)%velocity(dim))/time_step) ! 1/rho_g*f_b (4.43 FDS guide)

						end do 
					end if
					
					cell = this%get_droplet_cell(this%droplets(drop)%coords,out_flag)
					if (out_flag) then
						this%droplets(drop)%outside_domain = out_flag
						exit
					end if
			
					do dim = 1, dimensions
                        
						lower_face	= mesh%mesh(dim,i,j,k) - 0.5_dkind*cell_size(dim)
						higher_face	= mesh%mesh(dim,i,j,k) + 0.5_dkind*cell_size(dim)
						
						a = this%droplets(drop)%coords(dim) - lower_face
						
						b = higher_face - this%droplets(drop)%coords(dim) 
					
						v_prod%pr(dim)%cells(dim,i,j,k)										= v_prod%pr(dim)%cells(dim,i,j,k)									- (1.0_dkind - a/cell_size(dim)) * F_a(dim) 
						v_prod%pr(dim)%cells(dim,i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3))	= v_prod%pr(dim)%cells(dim,i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3))	- (1.0_dkind - b/cell_size(dim)) * F_a(dim) 
						
                        
					!	v_prod%pr(dim)%cells(dim,i,j,k)										= v_prod%pr(dim)%cells(dim,i,j,k)									- (1.0_dkind - this%droplets(drop)%coords(dim)) * droplet_acceleration(dim) 
					!	v_prod%pr(dim)%cells(dim,i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3))	= v_prod%pr(dim)%cells(dim,i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3))	- (this%droplets(drop)%coords(dim)) * droplet_acceleration(dim) 
					end do
			
			
					!# Heat and evaporate droplet, account energy transfer
					Pr_p = 0.7_dkind
					Sc_p = 0.6_dkind
			
					Nu_p = 2.0_dkind + 0.6_dkind * sqrt(Re_p) * Pr_p ** (1.0_dkind / 3.0_dkind)
					Sh_p = 2.0_dkind + 0.6_dkind * sqrt(Re_p) * Sc_p ** (1.0_dkind / 3.0_dkind)
			
					D_air_water = D%pr(droplet_material_index)%cells(i,j,k)

					H_h	= Nu_p * kappa%cells(i,j,k)	/ local_diameter
					H_m	= Sh_p * D_air_water / local_diameter
			
					H_v	= 2454208.88831280_dkind ! droplet%material_latent_heat
			
					mol_mix_w = 0.0_dkind
					do spec = 1, species_number
						mol_mix_w = mol_mix_w + Yg_old(spec)/this%thermo%thermo_ptr%molar_masses(spec)
					end do
					mol_mix_w = 1.0_dkind / mol_mix_w
			
					W_ratio	= mol_mix_w*(1.0_dkind - Yg_old(droplet_material_index)) / this%thermo%thermo_ptr%molar_masses(droplet_material_index)
			
					X_vap = min(1.0_dkind,  exp(H_v * this%thermo%thermo_ptr%molar_masses(droplet_material_index) / r_gase_J * (1.0_dkind / droplet%material_boiling_temperature - 1.0_dkind / this%droplets(drop)%temperature)))
					Y_vap = X_vap / (X_vap * (1.0_dkind - W_ratio) + W_ratio)
			
					dY_dT = (W_ratio/(X_vap*(1.0_dkind - W_ratio) + W_ratio)**2) * H_v * this%thermo%thermo_ptr%molar_masses(droplet_material_index) / r_gase_J * X_vap / this%droplets(drop)%temperature ** 2
				
					specie_enthalpy_gas		= (this%thermo%thermo_ptr%calculate_specie_cp(Tg_old,droplet_material_index))*Tg_old / this%thermo%thermo_ptr%molar_masses(droplet_material_index)
					specie_enthalpy_liquid	= (this%thermo%thermo_ptr%calculate_specie_cp(this%droplets(drop)%temperature,droplet_material_index))*this%droplets(drop)%temperature / this%thermo%thermo_ptr%molar_masses(droplet_material_index)
			
					mixture_cp = 0.0
					do spec = 1, species_number
						mixture_cp		= mixture_cp + this%thermo%thermo_ptr%calculate_specie_cp(Tg_old,spec)	* Yg_old(spec) / this%thermo%thermo_ptr%molar_masses(spec)
					end do			

					DTOG = droplet_time_step / (2.0_dkind * mixture_cp * M_gas_old)
					DTOP = droplet_time_step / (2.0_dkind * droplet%material_heat_capacity	* m_liq)
					DTGOG = droplet_time_step * A_p * H_h / (2.0_dkind * mixture_cp * M_gas_old)
					DTGOP = droplet_time_step * A_p * H_h / (2.0_dkind * droplet%material_heat_capacity	* m_liq)
				
					AGHRHO = A_p * H_m * rhog_old / ( 1.0_dkind + 0.5_dkind * droplet_time_step * A_p * H_m / cell_volume)
			
					DADYDTHVHL	= DTOG * AGHRHO	* (specie_enthalpy_liquid - specie_enthalpy_gas)
					DADYDTHV	= DTOP * AGHRHO * droplet%material_latent_heat
			
					A_col(1)	= 1.0_dkind + DTGOG
					B_col(1)	= - (DTGOG + DADYDTHVHL * dY_dT)
					A_col(2)	= -DTGOP
					B_col(2)	= 1.0_dkind + DTGOP + DADYDTHV * dY_dT
					D_vec(1)	= (1.0_dkind - DTGOG) * Tg_old + (DTGOG - DADYDTHVHL * dY_dT) * this%droplets(drop)%temperature + 2.0_dkind * DADYDTHVHL * (Y_vap - Yg_old(droplet_material_index))
					D_vec(2)	= DTGOP * Tg_old + (1.0_dkind - DTGOP + DADYDTHV * dY_dT) * this%droplets(drop)%temperature - 2.0_dkind * DADYDTHV * (Y_vap - Yg_old(droplet_material_index))
			
					Tp_new	= -(A_col(2) * D_vec(1) - A_col(1) * D_vec(2))/(A_col(1) * B_col(2) - A_col(2) * B_col(1))
					Tg_new	= (D_vec(1) - B_col(1) * Tp_new) / A_col(1)
			
					dmp = max(0.0_dkind, min(this%droplets(drop)%mass, droplet_time_step * AGHRHO * ( Y_vap -  Yg_old(droplet_material_index) + 0.5_dkind * dY_dT * (Tp_new - this%droplets(drop)%temperature))))
			
					Q	= droplet_time_step * A_p * H_h * (Tg_new - 0.5_dkind * (Tp_old + Tp_new))!(Tg_old + Tg_new - Tp_old - Tp_new)
					if (Q > this%droplets(drop)%mass * H_v) dmp = this%droplets(drop)%mass
			
					if ( dmp < this%droplets(drop)%mass) then
						Tp_new = Tp_old + (Q - dmp * H_v) / (droplet%material_heat_capacity * ( this%droplets(drop)%mass - dmp))
						if ( Tp_new > droplet%material_boiling_temperature) then
							dmp = min(this%droplets(drop)%mass, (Q - this%droplets(drop)%mass * droplet%material_heat_capacity * (droplet%material_boiling_temperature - Tp_old))/(H_v - droplet%material_heat_capacity * (droplet%material_boiling_temperature - Tp_old)))
							if ( dmp == this%droplets(drop)%mass) then
								Q = dmp * H_v
							end if
							Tp_new = droplet%material_boiling_temperature
						end if
					else
						Q		= dmp * H_v
						Tp_new	= droplet%material_boiling_temperature
						this%droplets(drop)%outside_domain = .true.
					end if
			
					this%droplets(drop)%mass			= this%droplets(drop)%mass - dmp
					this%droplets(drop)%temperature		= Tp_new
					this%droplets(drop)%dm				= this%droplets(drop)%dm + dmp / droplet_time_step
			
					M_gas_new	= M_gas_old + dmp
					Yg_new		= Yg_old*M_gas_old/M_gas_new
					Yg_new(droplet_material_index)		=	Yg_new(droplet_material_index) + dmp/M_gas_new
			
					Hg_new		= Hg_old	+ (Hp_old - this%droplets(drop)%mass*droplet%material_heat_capacity*Tp_new)
					Tg_new		= Tg_old	+ (Hg_new	- mixture_cp * Tg_old * M_gas_new) / M_gas_new / mixture_cp		!# Assuming cp = const
			
					rhog_new	= M_gas_new / cell_volume
			
					delta_H		= ((this%thermo%thermo_ptr%calculate_specie_cp(Tp_old,droplet_material_index))*Tp_old - (this%thermo%thermo_ptr%calculate_specie_cp(Tg_old,droplet_material_index))*Tg_old) / this%thermo%thermo_ptr%molar_masses(droplet_material_index)
			
					E_f_prod%cells(i,j,k)	=	E_f_prod%cells(i,j,k) + (W_ratio * dmp / M_gas_old + (dmp * delta_H - Q) / Hg_old) / time_step
					
					do dim = 1, dimensions
						dmp = dmp / cell_size(dim)
					end do
					
					Y_prod%pr(droplet_material_index)%cells(i,j,k)	= Y_prod%pr(droplet_material_index)%cells(i,j,k) + dmp / droplet_time_step

					Yg_old		= Yg_new
					Hg_old		= Hg_new
					rhog_old	= rhog_new
					Tg_old		= Tg_new
				end do
			else
			!	print *, 'droplet', drop, ' with coodinates', this%droplets(drop)%coords, '  is outside the domain'
			!	pause
			end if
		end do

		if ((time*1e06 >= 1.0_dkind*(output_counter+1)).or.(time==0.0_dkind)) then
			write (average_io_unit,'(3E14.6)')	time, average_temperature/droplets_inside, average_diameter/droplets_inside
		end if		
		
		if ((time*1e06 >= 1.0_dkind*(output_counter+1)).or.(time==0.0_dkind)) then
			write(file_path,'(I6.6,A)')  int(time*1e06),'us'
			file_name = 'data_save_droplets' // trim(fold_sep) // 'droplets_' // trim(file_path) // '.plt'
			open(newunit = lagrangian_droplets_io_unit, file = file_name, status = 'replace', form = 'formatted')
			do drop = 1, this%droplets_number
				if(.not.this%droplets(drop)%outside_domain) then
					write (lagrangian_droplets_io_unit,'(8E14.6)')	this%droplets(drop)%coords, this%droplets(drop)%velocity, this%droplets(drop)%temperature, this%droplets(drop)%mass			
				end if
			end do
			
			output_counter = output_counter + 1
			close(lagrangian_droplets_io_unit)
		end if
		
		if ((time*1e06 >= 50.0_dkind*(droplet_release_counter+1)).or.(time==0.0_dkind)) then
			droplet_release_counter = droplet_release_counter + 1
		!	call this%release_droplet()
		end if	


		time = time + time_step

		
		!!$omp end do nowait
		!!$omp end parallel		
		
		end associate
	end subroutine

	function get_droplet_cell(this,droplet_coords, out_flag)
		
		class(lagrangian_droplets_solver)			,intent(inout)					:: this
		real(dkind)		,dimension(3)				,intent(in)						:: droplet_coords
		logical										,intent(inout)	,optional		:: out_flag
		
		integer			,dimension(3)	:: get_droplet_cell
		integer		,dimension(3,2)	:: cons_inner_loop
		real(dkind)	,dimension(3)	:: cell_size
		
		integer	:: dimensions
		
		integer	:: sign, bound_number
		integer :: i,j,k,plus,dim,dim1,specie_number
		
		dimensions		= this%domain%get_domain_dimensions()
				
		cons_inner_loop	= this%domain%get_local_inner_cells_bounds()		
		
		cell_size		= this%mesh%mesh_ptr%get_cell_edges_length()
		
		get_droplet_cell	= 1
		if(present(out_flag))	out_flag	= .false.
		
		get_droplet_cell(:dimensions) = 0
		associate(	bc		=> this%boundary%bc_ptr	,&
					mesh	=> this%mesh%mesh_ptr)		
		
		if (dimensions	== 3) then
			do k = cons_inner_loop(3,1),cons_inner_loop(3,2)
				if (( 0.5_dkind*(mesh%mesh(3,1,1,k-1) + mesh%mesh(3,1,1,k)) <= droplet_coords(3)).and. &
					( 0.5_dkind*(mesh%mesh(3,1,1,k+1) + mesh%mesh(3,1,1,k)) > droplet_coords(3))) then
					get_droplet_cell(3) = k
					exit
				end if		
			end do
		end if
		
		
		if (dimensions	== 2) then
			do j = cons_inner_loop(2,1),cons_inner_loop(2,2)
				if (( 0.5_dkind*(mesh%mesh(2,1,j-1,get_droplet_cell(3)) + mesh%mesh(2,1,j,get_droplet_cell(3))) <= droplet_coords(2)).and. &
					( 0.5_dkind*(mesh%mesh(2,1,j+1,get_droplet_cell(3)) + mesh%mesh(2,1,j,get_droplet_cell(3))) > droplet_coords(2))) then
					get_droplet_cell(2) = j
					exit
				end if		
			end do
		end if
		
		do i = cons_inner_loop(1,1),cons_inner_loop(1,2)

			if ((  0.5_dkind*(mesh%mesh(1,i-1,get_droplet_cell(2),get_droplet_cell(3)) + mesh%mesh(1,i,get_droplet_cell(2),get_droplet_cell(3))) <= droplet_coords(1)).and. &
				(  0.5_dkind*(mesh%mesh(1,i+1,get_droplet_cell(2),get_droplet_cell(3)) + mesh%mesh(1,i,get_droplet_cell(2),get_droplet_cell(3))) > droplet_coords(1))) then				
				get_droplet_cell(1) = i
				exit
			end if		
		end do			
		
		if(present(out_flag)) then
			do dim = 1, dimensions 
				if (get_droplet_cell(dim) == 0 ) then 
					out_flag = .true.
					exit
				end if 
			end do
			if (bc%bc_markers(get_droplet_cell(1),get_droplet_cell(2),get_droplet_cell(3)) /= 0) then
				out_flag = .true.
			end if
		end if
		
		end associate
	
	end function
	
	integer	function count_droplets_in_cell(this,cell_indexes)
	
		class(lagrangian_droplets_solver)			,intent(inout)		:: this
		integer			,dimension(3)	,intent(in)			:: cell_indexes
		
		integer			,dimension(3)	:: cell
		logical							:: out_flag
		
		integer	:: dimensions
		
		integer	:: sign, bound_number
		integer :: drop
		
		count_droplets_in_cell = 0
		do drop = 1, this%droplets_number
			
			if (( this%droplets(drop)%cell(1) == cell_indexes(1)).and.(this%droplets(drop)%cell(2) == cell_indexes(2)).and.(this%droplets(drop)%cell(3) == cell_indexes(3)).and.(.not.out_flag))	count_droplets_in_cell = count_droplets_in_cell + 1
		end do
		
	end function
	
	subroutine apply_boundary_conditions_interm_v_d(this)
		class(lagrangian_droplets_solver)		,intent(inout)		:: this
		integer	:: dimensions
		integer	,dimension(3,2)	:: cons_inner_loop
		
		character(len=20)		:: boundary_type_name		
		integer	:: sign, bound_number
		integer :: i,j,k,plus,dim,dim1,specie_number
								
		dimensions		= this%domain%get_domain_dimensions()
				
		cons_inner_loop	= this%domain%get_local_inner_cells_bounds()			
		
		!associate(  v_d_int			=> this%v_d_int%v_ptr		, &
		!			bc				=> this%boundary%bc_ptr	, &
		!			mesh			=> this%mesh%mesh_ptr)
  !
		!!$omp parallel default(none)  private(i,j,k,plus,dim,dim1,sign,bound_number,boundary_type_name) , &
		!!$omp& firstprivate(this)	,&
		!!$omp& shared(v_d_int,bc,cons_inner_loop,dimensions)
		!!$omp do collapse(3) schedule(guided)
  !
		!	do k = cons_inner_loop(3,1),cons_inner_loop(3,2)
		!	do j = cons_inner_loop(2,1),cons_inner_loop(2,2)
		!	do i = cons_inner_loop(1,1),cons_inner_loop(1,2)
		!		if(bc%bc_markers(i,j,k) == 0) then
		!			do dim = 1,dimensions
		!				do plus = 1,2
		!					sign			= (-1)**plus
		!					bound_number	= bc%bc_markers(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3))
		!					if( bound_number /= 0 ) then
  !
		!						do dim1 = 1, dimensions
		!							if(dim1 == dim) then
		!								v_d_int%pr(dim1)%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3)) = - v_d_int%pr(dim1)%cells(i,j,k)
		!							else
		!								v_d_int%pr(dim1)%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3)) = v_d_int%pr(dim1)%cells(i,j,k)
		!							end if
		!						end do
		!						
		!						boundary_type_name = bc%boundary_types(bound_number)%get_type_name()
		!						select case(boundary_type_name)
		!							case ('outlet')
		!								do dim1 = 1, dimensions
		!									if(dim1 == dim)	v_d_int%pr(dim1)%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3)) = v_d_int%pr(dim1)%cells(i,j,k)
		!								end do
		!						end select
  !
		!					end if
		!				end do
		!			end do
		!		end if
		!	end do
		!	end do
		!	end do
  !
		!!$omp end do nowait
		!!$omp end parallel
  !
		!end associate
	end subroutine	
	
	subroutine apply_boundary_conditions_main(this)
		class(lagrangian_droplets_solver)		,intent(inout)		:: this
		character(len=20)		:: boundary_type_name
		real(dkind)				:: farfield_density, farfield_pressure, wall_temperature
		
		integer					:: dimensions
		integer	,dimension(3,2)	:: cons_inner_loop
		integer	:: sign, bound_number
		integer :: i,j,k,plus,dim,dim1
						
		dimensions		= this%domain%get_domain_dimensions()
			
		cons_inner_loop	= this%domain%get_local_inner_cells_bounds()			
		
		!associate(  T_d				=> this%T_d%s_ptr			, &
		!			rho_d			=> this%rho_d%s_ptr			, &
		!			mass_d          => this%mass_d%s_ptr        , &
		!			v_d				=> this%v_d%v_ptr			, &
		!			bc				=> this%boundary%bc_ptr		, &
		!			mesh			=> this%mesh%mesh_ptr)
  !
		!!$omp parallel default(none)  private(i,j,k,plus,dim,dim1,sign,bound_number,boundary_type_name,wall_temperature) , &
		!!$omp& firstprivate(this)	,&
		!!$omp& shared(T_d,rho_d,v_d,mass_d,mesh,bc,cons_inner_loop,dimensions)
		!!$omp do collapse(3) schedule(guided)
  !
		!	do k = cons_inner_loop(3,1),cons_inner_loop(3,2)
		!	do j = cons_inner_loop(2,1),cons_inner_loop(2,2)
		!	do i = cons_inner_loop(1,1),cons_inner_loop(1,2)
		!		if(bc%bc_markers(i,j,k) == 0) then
		!			do dim = 1,dimensions
		!				do plus = 1,2
		!					sign			= (-1)**plus
		!					bound_number	= bc%bc_markers(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3))
		!					if( bound_number /= 0 ) then
  !
		!						rho_d%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3))	= rho_d%cells(i,j,k)
		!						T_d%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3))	= T_d%cells(i,j,k)
  !                              mass_d%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3)) = mass_d%cells(i,j,k)
  !
		!						do dim1 = 1, dimensions
		!							if(dim1 == dim) then
		!								v_d%pr(dim1)%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3)) = - v_d%pr(dim1)%cells(i,j,k)
		!							else
		!								v_d%pr(dim1)%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3)) = v_d%pr(dim1)%cells(i,j,k)
		!							end if
		!						end do								
		!						
		!						boundary_type_name = bc%boundary_types(bound_number)%get_type_name()
		!						select case(boundary_type_name)
		!							case('wall')
		!								if(bc%boundary_types(bound_number)%is_conductive()) then 
		!									wall_temperature = bc%boundary_types(bound_number)%get_wall_temperature()
		!									T_d%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3)) = wall_temperature
		!								end if
		!								if(.not.bc%boundary_types(bound_number)%is_slip()) then
		!									do dim1 = 1, dimensions
		!										if (dim1 /= dim) then	
		!											v_d%pr(dim1)%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3)) = 0.0_dkind
		!										end if
		!									end do
		!								end if
		!						end select
		!						
		!					end if
		!				end do
		!			end do
		!		end if
		!	end do
		!	end do
		!	end do
  !
		!!$omp end do nowait
		!!$omp end parallel
  !
		!end associate
	end subroutine
	
	subroutine release_droplet(this)	
		class(lagrangian_droplets_solver)	,intent(inout)	:: this
		
		integer	:: dimensions
		integer :: i,j,k, drop_x, drop_y, drop_z, dim, drop

		integer, dimension(3)	:: cell, initial_cell
		logical	:: out_flag
		
		dimensions		= this%domain%get_domain_dimensions()
		
		associate(  v_f			=> this%v_f%v_ptr		, &
			droplet	    => this%droplets_params	, &
            chem        => this%chem%chem_ptr   , &
			mesh		=> this%mesh%mesh_ptr	, &
			bc			=> this%boundary%bc_ptr)
		
			do drop = 1, this%droplets_number
				if (this%droplets(drop)%outside_domain) then
				
					this%droplets(drop)%coords(1)		= 0.001_dkind

					this%droplets(drop)%dm				= 0.0_dkind  
					
					this%droplets(drop)%outside_domain = .false.
				
					this%droplets(drop)%temperature		= 300.0_dkind
					this%droplets(drop)%mass			= Pi*this%droplets_params%diameter**3 / 6.0_dkind * this%droplets_params%material_density
					
					initial_cell = this%get_droplet_cell(this%droplets(drop)%coords,out_flag)
					this%droplets(drop)%cell = initial_cell
					
					this%droplets(drop)%velocity(1)	= 5.0_dkind !v_f%pr(1)%cells(1,initial_cell(1),initial_cell(2),initial_cell(3))
					
					this%droplets(drop)%outside_domain = out_flag	
					
					exit
				end if
			end do	
	
		end associate		
		
	end subroutine
	
	
	
end module
