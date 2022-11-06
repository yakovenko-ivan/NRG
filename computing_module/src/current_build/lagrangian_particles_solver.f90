module lagrangian_particles_solver_class

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
	public	:: lagrangian_particles_solver, lagrangian_particles_solver_c

	type(field_scalar_cons)	,dimension(:)	,allocatable	,target	:: E_f_prod_p
	type(field_vector_flow)	,dimension(:)	,allocatable	,target	:: v_prod_p

	type	:: lagrangian_particle
		real(dkind)	,dimension(3)	:: coords, velocity
		integer		,dimension(3)	:: cell
		real(dkind)					:: temperature
		real(dkind)					:: mass
		logical						:: outside_domain
		integer						:: particle_number
	end type
	
	type	:: lagrangian_particles_solver
		type(field_scalar_cons_pointer)		:: T, rho, nu, kappa, E_f_prod, rho_prod, p, mol_mix_conc, h_s
		type(field_vector_cons_pointer)		:: Y, Y_prod, D
		type(field_vector_flow_pointer)		:: v_f, v_prod
		type(computational_domain)			:: domain

		type(boundary_conditions_pointer)	:: boundary
		type(computational_mesh_pointer)	:: mesh
		type(chemical_properties_pointer)	:: chem
		type(thermophysical_properties_pointer)		:: thermo
		
		type(solid_particles_phase)			:: particles_params

		real(dkind)							:: particle_mass
        		
		type(lagrangian_particle)	,dimension(:)	,allocatable	:: particles	
		
		integer								:: particles_number
        integer                             :: phase_number
	contains
		procedure				::  set_initial_distributions
		procedure				::	particles_solve
		procedure				::  apply_boundary_conditions_main
		procedure				::  apply_boundary_conditions_interm_v_d
		procedure				::	pre_constructor
		
		procedure	,private	:: release_particle
		procedure	,private	:: get_particle_cell
		procedure	,private	:: count_particles_in_cell
	end type

	interface	lagrangian_particles_solver_c
		module procedure	constructor
	end interface

contains

	subroutine pre_constructor(this,number_of_phases)
		class(lagrangian_particles_solver)	,intent(inout)	:: this	
		integer					,intent(in)		:: number_of_phases
		
		allocate(	E_f_prod_p(number_of_phases), v_prod_p(number_of_phases))	
	
	end subroutine

	type(lagrangian_particles_solver)	function constructor(manager,solid_particles,phase_number)

		type(data_manager)			, intent(inout)	:: manager
		type(solid_particles_phase)	, intent(in)	:: solid_particles
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
		
		integer	:: part, dim
		
		cons_allocation_bounds		= manager%domain%get_local_utter_cells_bounds()
		flow_allocation_bounds		= manager%domain%get_local_utter_faces_bounds()
		dimensions					= manager%domain%get_domain_dimensions()   		
		
		constructor%particles_params	= solid_particles
		
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

		write(var_name,'(A,I2.2)')		'energy_production_particles', phase_number
		write(var_short_name,'(A,I2.2)')	'E_f_prod_p', phase_number
		call manager%create_scalar_field(E_f_prod_p(phase_number),	var_name,	var_short_name)
		constructor%E_f_prod%s_ptr		=> E_f_prod_p(phase_number)

		write(var_name,'(A,I2.2)')		'velocity_production_particles', phase_number
		write(var_short_name,'(A,I2.2)')	'v_prod_p', phase_number		
		call manager%create_vector_field(v_prod_p(phase_number),	var_name,	var_short_name,	'spatial')
		constructor%v_prod%v_ptr		=> v_prod_p(phase_number)

		constructor%phase_number    = phase_number

		constructor%mesh%mesh_ptr	=> manager%computational_mesh_pointer%mesh_ptr
		constructor%boundary%bc_ptr => manager%boundary_conditions_pointer%bc_ptr
		constructor%thermo%thermo_ptr	=> manager%thermophysics%thermo_ptr
		constructor%chem%chem_ptr	=> manager%chemistry%chem_ptr
		constructor%domain			= manager%domain
		
		system_command = 'mkdir data_save_particles'
		call system(system_command)		
		
	end function

	subroutine set_initial_distributions(this)
		class(lagrangian_particles_solver)	,intent(inout)	:: this

		integer		,dimension(3,2)	:: cons_inner_loop
		real(dkind)	,dimension(3)	:: cell_size
		
		real(dkind)	:: delta, part_number_max, a ,b, dimless_length
		integer		,dimension(3)	:: part_number
		real(dkind)	,dimension(3)	:: coords
		real(dkind)	,dimension(:,:)	,allocatable	:: lengths
		real(dkind)	:: gamma1, gamma2, d
		real(dkind)	,dimension(2)	:: AA
		
		integer	:: dimensions
		integer :: i,j,k, part_x, part_y, part_z, dim, particle

		integer, dimension(3)	:: cell, initial_cell
		logical	:: out_flag
		
		dimensions		= this%domain%get_domain_dimensions()
		cons_inner_loop	= this%domain%get_local_inner_cells_bounds()		

		allocate(lengths(dimensions,2))
		lengths			= this%domain%get_domain_lengths()
		
		cell_size		= this%mesh%mesh_ptr%get_cell_edges_length()
		
!		delta = 0.0006_dkind
		
		
!		lengths(1,1) = 0.01_dkind
!		lengths(1,2) = 0.02_dkind
		
		part_number = 1
!		do dim = 1,dimensions
!			part_number(dim) =  (lengths(dim,2)-lengths(dim,1)) /delta
!		end do
		
		part_number(1) = 100
		part_number(2) = 100
		
		delta = (lengths(1,2)-lengths(1,1)) / part_number(1)
		
		this%particles_number = part_number(1)*part_number(2)*part_number(3)
		
		allocate(this%particles(this%particles_number))	
		
		do particle = 1, this%particles_number
			this%particles(particle)%outside_domain	= .true.
			this%particles(particle)%coords			= 0.0_dkind
			this%particles(particle)%velocity		= 0.0_dkind
		end do 		
		
		!do part_x = 1, part_number(1)
		!do part_y = 1, part_number(2)
		!do part_z = 1, part_number(3)
		!
		!	particle = part_x + (part_y-1)*part_number(1) + (part_z-1)*part_number(1)*part_number(2)
		!
		!	this%particles(particle)%coords			= 0.0_dkind
		!	this%particles(particle)%velocity		= 0.0_dkind
  !
		!	this%particles(particle)%outside_domain = .false.
		!		
		!	this%particles(particle)%temperature	= 300.0_dkind
		!	this%particles(particle)%mass			= Pi*this%particles_params%diameter**3 / 6.0_dkind * this%particles_params%material_density
		!	
		!	
		!	call sleepqq(1)
		!	
		!	d = 2
		!	do while (d > 1)
		!	
		!		call RANDOM_SEED()
		!		call RANDOM_NUMBER(gamma1)
		!		call RANDOM_NUMBER(gamma2)
		!		
		!		AA(1) = 1.0_dkind - 2.0_dkind*gamma1
		!		AA(2) = 1.0_dkind - 2.0_dkind*gamma2
		!	
		!		d = AA(1)*AA(1) + AA(2)*AA(2)
		!	end do
		!	
		!	do dim = 1, dimensions
		!		dimless_length = real((part_x-1)*I_m(dim,1)+(part_y-1)*I_m(dim,2)+(part_z-1)*I_m(dim,3),dkind)
		!	
		!		this%particles(particle)%velocity(dim)	= AA(dim)/sqrt(d)
		!		this%particles(particle)%coords(dim)	= lengths(dim,1) + (dimless_length+0.5)*delta
		!	end do
		!	
		!	initial_cell = this%get_particle_cell(this%particles(particle)%coords,out_flag)
		!	this%particles(particle)%cell = initial_cell
		!	this%particles(particle)%outside_domain = out_flag	
		!	this%particles(particle)%particle_number = particle
		!end do
		!end do
		!end do
				
		continue
		
	end subroutine
	
	
	subroutine particles_solve(this,time_step)

		class(lagrangian_particles_solver)	,intent(inout)	:: this
		real(dkind)				,intent(in)		:: time_step
		real(dkind)								:: particle_time_step
		integer									:: particle_iterations
		
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
		real(dkind)					:: dmp, Q,	Tg_new, Tp_new, Tg_old, Tp_old, rhog_old, rhog_new, M_gas_old, M_gas_new, m_p, mol_mix_w, Hg_old, Hg_new, Hp_old, Hp_new
		real(dkind)					:: delta_H
		real(dkind)	,dimension(:)	, allocatable	:: Yg_old, Yg_new
		real(dkind)					:: DTOP, DTOG, DTGOG, DTGOP, AGHRHO, DADYDTHVHL, DADYDTHV
		real(dkind)	,dimension(2)	:: A_col, B_col, D_vec
		
		integer						:: particles_in_cell
		character(len=100)			:: file_path, file_name
		
		real(dkind)                 :: local_diameter, evaporation_rate, velocity_2, specie_enthalpy
		real(dkind)					:: v_prod_p1, E_prod_p1
		
		integer	:: dimensions, species_number
		integer, dimension(3)	:: cell, initial_cell

		integer :: i,j,k,dim, part, part_iter, spec
		logical	:: out_flag
		
		real(dkind)	,save	:: time = 0.0_dkind
		integer		,save	:: output_counter = 0
		integer		,save	:: particle_release_counter = 0
		
		integer				:: lagrangian_particles_io_unit
		
		dimensions		= this%domain%get_domain_dimensions()
		species_number	= this%chem%chem_ptr%species_number
		
		allocate(Yg_old(species_number),Yg_new(species_number))
		
		cons_inner_loop	= this%domain%get_local_inner_cells_bounds()		

		cell_size		= this%mesh%mesh_ptr%get_cell_edges_length()
		cell_volume		= this%mesh%mesh_ptr%get_cell_volume()

		associate(  T			=> this%T%s_ptr				, &
					rho			=> this%rho%s_ptr			, &
					p			=> this%p%s_ptr				, &
					h_s			=> this%h_s%s_ptr			, &
					Y			=> this%Y%v_ptr				, &
                    v_f			=> this%v_f%v_ptr			, &
					D			=> this%D%v_ptr				, &
					kappa		=> this%kappa%s_ptr			, &
					nu			=> this%nu%s_ptr			, &
					E_f_prod	=> this%E_f_prod%s_ptr		, &
					v_prod		=> this%v_prod%v_ptr		, &
					particle	=> this%particles_params	, &
                    chem        => this%chem%chem_ptr		, &
					mesh		=> this%mesh%mesh_ptr		, &
					bc			=> this%boundary%bc_ptr)

		!!$omp parallel default(none)  private() , &
		!!$omp& firstprivate(this)	,&
		!!$omp& shared()
		!!$omp do collapse(3) schedule(guided)					
			
		E_f_prod%cells = 0.0_dkind		
		
		do dim = 1, dimensions
			v_prod%pr(dim)%cells = 0.0_dkind
		end do
		
		do part = 1, this%particles_number
		
			if (.not.this%particles(part)%outside_domain) then
		
				!# Move particle, account momentum transfer
		
				particle_time_step = time_step
				do dim = 1, dimensions
					if (abs(this%particles(part)%velocity(dim)) > 1e-10_dkind) then
						if (abs(cell_size(dim) / this%particles(part)%velocity(dim)) < particle_time_step) particle_time_step = abs(cell_size(dim) / this%particles(part)%velocity(dim))
					end if
				end do
				particle_iterations	= ceiling(time_step/ particle_time_step)
				particle_time_step	= time_step / real(particle_iterations,dkind)

				F_a =	0.0_dkind			
			
				initial_cell = this%get_particle_cell(this%particles(part)%coords,out_flag)
				i = initial_cell(1)
				j = initial_cell(2)
				k = initial_cell(3)
				this%particles(part)%cell = initial_cell

				Tg_old		= T%cells(i,j,k)
				rhog_old	= rho%cells(i,j,k)
				
				do spec = 1, species_number
					Yg_old(spec)	= Y%pr(spec)%cells(i,j,k)
				end do	
				M_gas_old	= rhog_old * cell_volume 
				Hg_old		= h_s%cells(i,j,k) * M_gas_old
			
				do part_iter = 1, particle_iterations
			
					Tp_old	= this%particles(part)%temperature
					Hp_old	= this%particles(part)%temperature * particle%material_heat_capacity * this%particles(part)%mass
			
					local_diameter = (6.0_dkind * this%particles(part)%mass / Pi /  particle%material_density) ** (1.0_dkind / 3.0_dkind)
				
					cell = this%get_particle_cell(this%particles(part)%coords,out_flag)
					i = cell(1)
					j = cell(2)
					k = cell(3)				
				
					if ((this%particles(part)%cell(1) /= cell(1)).and.(this%particles(part)%cell(2) /= cell(2)).and.(this%particles(part)%cell(3) /= cell(3))) then
						Tg_old		= T%cells(i,j,k)
						rhog_old	= rho%cells(i,j,k)
						M_gas_old	= rhog_old * cell_volume 
						Hg_old		= h_s%cells(i,j,k) * M_gas_old
						do spec = 1, species_number
							Yg_old(spec)	= Y%pr(spec)%cells(i,j,k)
						end do
					
						E_f_prod%cells(i,j,k) = 0.0_dkind
				
						this%particles(part)%cell = cell
					end if

					old_velocity = this%particles(part)%velocity
			
				!# Interpolate gas velocity to the particle location using linear interpolation
					gas_velocity = 0.0_dkind
					do dim = 1, dimensions
						lower_face	= mesh%mesh(dim,i,j,k) - 0.5_dkind*cell_size(dim)
						higher_face	= mesh%mesh(dim,i,j,k) + 0.5_dkind*cell_size(dim)
					
						b = (v_f%pr(dim)%cells(dim,i,j,k)*higher_face - v_f%pr(dim)%cells(dim,i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3))*lower_face) / (higher_face - lower_face) 
						a = (v_f%pr(dim)%cells(dim,i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3)) - v_f%pr(dim)%cells(dim,i,j,k)) / (higher_face - lower_face) 
					
						gas_velocity(dim) = a * this%particles(part)%coords(dim) + b
					end do
				
					abs_relative_velocity = 0.0_dkind
					do dim = 1, dimensions
						relative_velocity(dim)	= this%particles(part)%velocity(dim) - gas_velocity(dim)
						abs_relative_velocity	= abs_relative_velocity + relative_velocity(dim)*relative_velocity(dim)
					end do
				
					abs_relative_velocity = sqrt(abs_relative_velocity)
				
					!# Calculate particle mass and cross section
					m_p	= Pi*local_diameter**3 / 6.0_dkind * particle%material_density
					A_p		= Pi*local_diameter**2
					!# Calculate particle Reynolds number
					Re_p	= rhog_old * abs_relative_velocity * local_diameter / nu%cells(i,j,k)
				
					if ( abs_relative_velocity > 1.0e-10_dkind) then	
					!# Calculate drag coefficient for spherical particle
						if (Re_p < 1e-010) then
							C_drag = 100
						elseif (Re_p < 1.0_dkind) then
								C_drag = 24.0_dkind / Re_p
						elseif (Re_p < 1000.0_dkind) then
								C_drag = 24.0_dkind * ( 0.85_dkind + 0.15 * Re_p**0.687) / Re_p
						elseif (Re_p >= 1000.0_dkind) then
								C_drag = 0.44_dkind
						end if
				
						particles_in_cell = this%count_particles_in_cell(cell)
				
					!# Calculate new particle coordinates and velocities
						alpha_p	= M_gas_old / m_p
						beta_p	= 0.5_dkind * rhog_old * C_drag * A_p * ( 1/m_p + 1/M_gas_old) * abs_relative_velocity
					
						do dim = 1, dimensions
							this%particles(part)%coords(dim) =	this%particles(part)%coords(dim) + (this%particles(part)%velocity(dim) + alpha_p * gas_velocity(dim)) * particle_time_step / (1.0_dkind + alpha_p) + &
																alpha_p * log(1.0_dkind + beta_p*particle_time_step) / beta_p / (1.0_dkind + alpha_p) * relative_velocity(dim)

							this%particles(part)%velocity(dim)	= (this%particles(part)%velocity(dim) + (this%particles(part)%velocity(dim) + alpha_p*gas_velocity(dim)) * beta_p * particle_time_step / (1.0_dkind + alpha_p)) / (1.0_dkind + beta_p*particle_time_step) !+ !... 	
														
							F_a(dim) = F_a(dim) + 1.0_dkind/(cell_volume * rhog_old) * (m_p*(old_velocity(dim) - this%particles(part)%velocity(dim))/time_step) ! 1/rho_g*f_b (4.43 FDS guide)
						end do 
					end if
					
					cell = this%get_particle_cell(this%particles(part)%coords,out_flag)
					if (out_flag) then
						this%particles(part)%outside_domain = out_flag
						exit
					end if
			
					do dim = 1, dimensions
					
						lower_face	= mesh%mesh(dim,i,j,k) - 0.5_dkind*cell_size(dim)
						higher_face	= mesh%mesh(dim,i,j,k) + 0.5_dkind*cell_size(dim)
						
						a = this%particles(part)%coords(dim) - lower_face
						
						b = higher_face - this%particles(part)%coords(dim) 
					
						v_prod%pr(dim)%cells(dim,i,j,k)										= v_prod%pr(dim)%cells(dim,i,j,k)									- (1.0_dkind - a/cell_size(dim)) * F_a(dim) 
						v_prod%pr(dim)%cells(dim,i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3))	= v_prod%pr(dim)%cells(dim,i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3))	- (1.0_dkind - b/cell_size(dim)) * F_a(dim) 
						
						v_prod_p1 =  F_a(dim)
					end do
			
			
					!# Heat and evaporate particle, account energy transfer
					Pr_p = 0.7_dkind
					Sc_p = 0.6_dkind
			
					Nu_p = 2.0_dkind + 0.6_dkind * sqrt(Re_p) * Pr_p ** (1.0_dkind / 3.0_dkind)
					Sh_p = 2.0_dkind + 0.6_dkind * sqrt(Re_p) * Sc_p ** (1.0_dkind / 3.0_dkind)

					H_h	= Nu_p * kappa%cells(i,j,k)	/ local_diameter
			
					mol_mix_w = 0.0_dkind
					do spec = 1, species_number
						mol_mix_w = mol_mix_w + Yg_old(spec)/this%thermo%thermo_ptr%molar_masses(spec)
					end do
					mol_mix_w = 1.0_dkind / mol_mix_w
			
					mixture_cp = 0.0
					do spec = 1, species_number
						mixture_cp		= mixture_cp + this%thermo%thermo_ptr%calculate_specie_cp(Tg_old,spec)	* Yg_old(spec) / this%thermo%thermo_ptr%molar_masses(spec)
					end do			

					DTOG = particle_time_step / (2.0_dkind * mixture_cp * M_gas_old)
					DTOP = particle_time_step / (2.0_dkind * particle%material_heat_capacity	* m_p)
					DTGOG = particle_time_step * A_p * H_h / (2.0_dkind * mixture_cp * M_gas_old)
					DTGOP = particle_time_step * A_p * H_h / (2.0_dkind * particle%material_heat_capacity	* m_p)
				
					AGHRHO = A_p * H_m * rhog_old / ( 1.0_dkind + 0.5_dkind * particle_time_step * A_p * H_m / cell_volume)
			
					A_col(1)	= 1.0_dkind + DTGOG
					B_col(1)	= - (DTGOG)
					A_col(2)	= -DTGOP
					B_col(2)	= 1.0_dkind + DTGOP
					D_vec(1)	= (1.0_dkind - DTGOG) * Tg_old + (DTGOG ) * this%particles(part)%temperature
					D_vec(2)	= DTGOP * Tg_old + (1.0_dkind - DTGOP) * this%particles(part)%temperature
			
					Tp_new	= -(A_col(2) * D_vec(1) - A_col(1) * D_vec(2))/(A_col(1) * B_col(2) - A_col(2) * B_col(1))
					Tg_new	= (D_vec(1) - B_col(1) * Tp_new) / A_col(1)

					Q	= particle_time_step * A_p * H_h * (Tg_new - 0.5_dkind * (Tp_old + Tp_new))!(Tg_old + Tg_new - Tp_old - Tp_new)
		
					Tp_new = Tp_old + Q / (particle%material_heat_capacity * this%particles(part)%mass)
			
					this%particles(part)%temperature	= Tp_new
			
					M_gas_new	= M_gas_old
			
					Hg_new		= Hg_old	+ (Hp_old - this%particles(part)%mass*particle%material_heat_capacity*Tp_new)
					Tg_new		= Tg_old	+ (Hg_new	- mixture_cp * Tg_old * M_gas_new) / M_gas_new / mixture_cp		!# Assuming cp = const
			
					rhog_new	= M_gas_new / cell_volume

				!	E_f_prod%cells(i,j,k)	=	E_f_prod%cells(i,j,k) - (Q / Hg_old) / time_step

					E_prod_p1 = E_f_prod%cells(i,j,k)
					
					Hg_old		= Hg_new
					rhog_old	= rhog_new
					Tg_old		= Tg_new
				end do
			else
			!	print *, 'Particle', part, ' with coodinates', this%particles(part)%coords, '  is outside the domain'
			!	pause
			end if
		end do

		if ((time*1e06 >= 1.0_dkind*(output_counter+1)).or.(time==0.0_dkind)) then
 
	!		write(file_path,'(A,I6.6,A)') 'particle', int(time*1e06),'us'
	!		write(file_path,'(A,I3.3)')  'particle', this%particles(part)%particle_number
	!		file_name = 'data_save_particles' // trim(fold_sep) // 'particles_' // trim(file_path) // '.plt'
	!		open(newunit = lagrangian_particles_io_unit, file = file_name, status = 'replace', position = 'append', form = 'formatted')	
	!		
	!		do part = 1, this%particles_number
	!			if(.not.this%particles(part)%outside_domain) then
	!	
	!				write (lagrangian_particles_io_unit,'(7E14.6)')	time, this%particles(part)%coords(1:dimensions), this%particles(part)%velocity(1:dimensions), this%particles(part)%temperature, this%particles(part)%mass !, v_prod_p1, E_prod_p1	
	!				
	!			end if
	!		end do
	!		
	!		close(lagrangian_particles_io_unit)
	!		output_counter = output_counter + 1
	!	end if
		
			do part = 1, this%particles_number
				if(.not.this%particles(part)%outside_domain) then
			!		print *, 'writing data for part ', part
					write(file_path,'(A,I4.4)')  'particle', this%particles(part)%particle_number
					file_name = 'data_save_particles' // trim(fold_sep) // 'particles_' // trim(file_path) // '.plt'
					open(newunit = lagrangian_particles_io_unit, file = file_name, status = 'unknown', position = 'append', form = 'formatted')	
					write (lagrangian_particles_io_unit,'(7E14.6)')	time, this%particles(part)%coords(1:dimensions), this%particles(part)%velocity(1:dimensions), this%particles(part)%temperature, this%particles(part)%mass !, v_prod_p1, E_prod_p1	
					close(lagrangian_particles_io_unit)
				end if
			end do		
		
			
			output_counter = output_counter + 1
		end if		

		if ((time*1e06 >= 1.0_dkind*(particle_release_counter+1))) then
			particle_release_counter = particle_release_counter + 1
			call this%release_particle()
		end if			
		
		time = time + time_step

		
		!!$omp end do nowait
		!!$omp end parallel		
		
		end associate
	end subroutine

	function get_particle_cell(this,particle_coords, out_flag)
		
		class(lagrangian_particles_solver)			,intent(inout)					:: this
		real(dkind)		,dimension(3)				,intent(in)						:: particle_coords
		logical										,intent(inout)	,optional		:: out_flag
		
		integer			,dimension(3)	:: get_particle_cell
		integer		,dimension(3,2)	:: cons_inner_loop
		real(dkind)	,dimension(3)	:: cell_size
		
		integer	:: dimensions
		
		integer	:: sign, bound_number
		integer :: i,j,k,plus,dim,dim1,specie_number
		
		dimensions		= this%domain%get_domain_dimensions()
				
		cons_inner_loop	= this%domain%get_local_inner_cells_bounds()		
		
		cell_size		= this%mesh%mesh_ptr%get_cell_edges_length()
		
		get_particle_cell	= 1
		if(present(out_flag))	out_flag	= .false.
		
		get_particle_cell(:dimensions) = 0
		associate(	bc		=> this%boundary%bc_ptr	,&
					mesh	=> this%mesh%mesh_ptr)		
		
		if (dimensions	== 3) then
			do k = cons_inner_loop(3,1),cons_inner_loop(3,2)
				if (( 0.5_dkind*(mesh%mesh(3,1,1,k-1) + mesh%mesh(3,1,1,k)) <= particle_coords(3)).and. &
					( 0.5_dkind*(mesh%mesh(3,1,1,k+1) + mesh%mesh(3,1,1,k)) > particle_coords(3))) then
					get_particle_cell(3) = k
					exit
				end if		
			end do
		end if
		
		
		if (dimensions	== 2) then
			do j = cons_inner_loop(2,1),cons_inner_loop(2,2)
				if (( 0.5_dkind*(mesh%mesh(2,1,j-1,get_particle_cell(3)) + mesh%mesh(2,1,j,get_particle_cell(3))) <= particle_coords(2)).and. &
					( 0.5_dkind*(mesh%mesh(2,1,j+1,get_particle_cell(3)) + mesh%mesh(2,1,j,get_particle_cell(3))) > particle_coords(2))) then
					get_particle_cell(2) = j
					exit
				end if		
			end do
		end if
		
		do i = cons_inner_loop(1,1),cons_inner_loop(1,2)

			if ((  0.5_dkind*(mesh%mesh(1,i-1,get_particle_cell(2),get_particle_cell(3)) + mesh%mesh(1,i,get_particle_cell(2),get_particle_cell(3))) <= particle_coords(1)).and. &
				(  0.5_dkind*(mesh%mesh(1,i+1,get_particle_cell(2),get_particle_cell(3)) + mesh%mesh(1,i,get_particle_cell(2),get_particle_cell(3))) > particle_coords(1))) then				
				get_particle_cell(1) = i
				exit
			end if		
		end do			
		
		if(present(out_flag)) then
			do dim = 1, dimensions 
				if (get_particle_cell(dim) == 0 ) then 
					out_flag = .true.
					exit
				end if 
			end do
			if (bc%bc_markers(get_particle_cell(1),get_particle_cell(2),get_particle_cell(3)) /= 0) then
				out_flag = .true.
			end if
		end if
		
		end associate
	
	end function
	
	integer	function count_particles_in_cell(this,cell_indexes)
	
		class(lagrangian_particles_solver)			,intent(inout)		:: this
		integer			,dimension(3)	,intent(in)			:: cell_indexes
		
		integer			,dimension(3)	:: cell
		logical							:: out_flag
		
		integer	:: dimensions
		
		integer	:: sign, bound_number
		integer :: part
		
		count_particles_in_cell = 0
		do part = 1, this%particles_number
			
			if (( this%particles(part)%cell(1) == cell_indexes(1)).and.(this%particles(part)%cell(2) == cell_indexes(2)).and.(this%particles(part)%cell(3) == cell_indexes(3)).and.(.not.out_flag))	count_particles_in_cell = count_particles_in_cell + 1
		end do
		
	end function
	
	subroutine apply_boundary_conditions_interm_v_d(this)
		class(lagrangian_particles_solver)		,intent(inout)		:: this
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
		class(lagrangian_particles_solver)		,intent(inout)		:: this
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
	
	subroutine release_particle(this)	
		class(lagrangian_particles_solver)	,intent(inout)	:: this
		
		integer	:: dimensions
		integer :: i,j,k, part_x, part_y, part_z, dim, part

		integer, dimension(3)	:: cell, initial_cell
		logical	:: out_flag
		integer, save	:: part_number = 0
		
		dimensions		= this%domain%get_domain_dimensions()
		
		associate(  v_f			=> this%v_f%v_ptr			, &
				particle	    => this%particles_params	, &
					chem        => this%chem%chem_ptr		, &
					mesh		=> this%mesh%mesh_ptr		, &
					bc			=> this%boundary%bc_ptr)
		
		part_number = part_number + 1
		if (part_number <= 999) then
			do part = 1, this%particles_number
				if (this%particles(part)%outside_domain) then
				
			!		this%particles(part)%coords(1)		= 0.0005_dkind
			!		this%particles(part)%coords(2)		= 0.015_dkind
                    
            !# alpha = 4
            !        this%particles(part)%coords(1)		= 0.00295_dkind
            !        this%particles(part)%coords(2)		= 0.04_dkind
            !# alpha = 5
            !        this%particles(part)%coords(1)		= 0.00325_dkind
            !        this%particles(part)%coords(2)		= 0.04_dkind     
            !# alpha = 6
                    this%particles(part)%coords(1)		= 0.00345_dkind
                    this%particles(part)%coords(2)		= 0.04_dkind           
					
					this%particles(part)%outside_domain = .false.
				
					this%particles(part)%temperature	= 300.0_dkind
					this%particles(part)%mass			= Pi*this%particles_params%diameter**3 / 6.0_dkind * this%particles_params%material_density
					
					initial_cell = this%get_particle_cell(this%particles(part)%coords,out_flag)
					this%particles(part)%cell = initial_cell
					
                    this%particles(part)%velocity(1)	= 0.0_dkind
                    this%particles(part)%velocity(2)	= 50.0_dkind
                    
            !# Normal        
			!		this%particles(part)%velocity(1)	= -50.0_dkind * 10.0_dkind / 101**0.5!0.0_dkind
			!		this%particles(part)%velocity(2)	=  50.0_dkind * 1.0_dkind / 101**0.5!0.0_dkind! 0.99_dkind * v_f%pr(2)%cells(2,initial_cell(1),initial_cell(2),initial_cell(3))
				
            !# 45 deg
            !		this%particles(part)%velocity(1)	= -50.0_dkind * 19.0_dkind / (19.0**2 + 21.0**2)**0.5!0.0_dkind
			!		this%particles(part)%velocity(2)	=  50.0_dkind * 21.0_dkind / (19.0**2 + 21.0**2)**0.5!0.0_dkind!
            		
            !# 15 deg
            		this%particles(part)%velocity(1)	= -50.0_dkind * sin(pi*15/180)!0.0_dkind
					this%particles(part)%velocity(2)	=  50.0_dkind * cos(pi*15/180)!0.0_dkind!
                    
					this%particles(part)%outside_domain = out_flag	
					
					this%particles(part)%particle_number = part_number
					exit
				end if
			end do	
        end if
        
		end associate		
		
	end subroutine	
	
	
end module
