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

	type(field_scalar_cons)	,dimension(:)	,allocatable	,target	:: E_f_prod_p, rho_prod_p
	type(field_vector_cons)	,dimension(:)	,allocatable	,target	:: Y_prod_p, v_prod_p

	integer	            :: average_io_unit, trajectory_io_unit
	character(len=500)  :: av_header, grid_header, scatter_header, trajectory_header
	integer	:: tracked_particle
    
    real(dp)	:: streams_distance, release_time, particles_velocity
    	
	type	:: lagrangian_particle
		real(dp)	,dimension(3)	:: coords, velocity, coords_prev, coords0
		integer		,dimension(3)	:: cell
		real(dp)					:: temperature
		real(dp)					:: mass	
		real(dp)					:: dm
        logical                     :: evaporating, heating, inertial
        
		logical						:: outside_domain
	end type
	
	type	:: part_grid
		real(dp)	,dimension(:,:,:,:)	,	allocatable	:: coords, velocity
        real(dp)	,dimension(:,:,:)		,	allocatable	:: density, temperature
        integer		,dimension(3)						:: part_grid_size
        real(dp)                                        :: grid_cell_size
	end type    
    
	type	:: lagrangian_particles_solver
		type(field_scalar_cons_pointer)		:: T, rho, nu, kappa, E_f_prod, rho_prod, p, mol_mix_conc, h_s
		type(field_vector_cons_pointer)		:: Y, Y_prod, D, v_prod
		type(field_vector_flow_pointer)		:: v_f
		type(computational_domain)			:: domain

		type(boundary_conditions_pointer)	:: boundary
		type(computational_mesh_pointer)	:: mesh
		type(chemical_properties_pointer)	:: chem
		type(thermophysical_properties_pointer)		:: thermo
		
		type(particles_phase)               :: particles_params

		real(dp)							:: particle_mass
        real(dp)    , dimension(3)          :: g
                
		type(lagrangian_particle)	,dimension(:)	,allocatable	:: particles	
		
        type(part_grid)						:: part_grid
		
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
		integer					            ,intent(in)		:: number_of_phases
		
		allocate(	E_f_prod_p(number_of_phases), rho_prod_p(number_of_phases), &
					v_prod_p(number_of_phases), Y_prod_p(number_of_phases))	
	
	end subroutine

	type(lagrangian_particles_solver)	function constructor(manager, particles, phase_number)

		type(data_manager)		, intent(inout)	:: manager
		type(particles_phase)	, intent(in)	:: particles
		integer					, intent(in)	:: phase_number

		type(field_scalar_cons_pointer)	:: scal_c_ptr
		type(field_vector_cons_pointer)	:: vect_c_ptr
		type(field_tensor_cons_pointer)	:: tens_c_ptr		

		type(field_scalar_flow_pointer)	:: scal_f_ptr		
		type(field_vector_flow_pointer)	:: vect_f_ptr
		type(field_tensor_flow_pointer)	:: tens_f_ptr			
		
 		integer					:: dimensions       
 		integer	,dimension(3,2)	:: cons_allocation_bounds, flow_allocation_bounds   
        
		character(len=100)		:: system_command
		character(len=40)		:: var_name, var_short_name
		character(len=100)		:: file_name
        
		integer	:: part, dim
		
		cons_allocation_bounds		= manager%domain%get_local_utter_cells_bounds()
		flow_allocation_bounds		= manager%domain%get_local_utter_faces_bounds()
		dimensions					= manager%domain%get_domain_dimensions()   		
		
		constructor%particles_params	= particles
		
		call manager%get_cons_field_pointer_by_name(scal_c_ptr,vect_c_ptr,tens_c_ptr,'density')
		constructor%rho%s_ptr				=> scal_c_ptr%s_ptr
		call manager%get_cons_field_pointer_by_name(scal_c_ptr,vect_c_ptr,tens_c_ptr,'temperature')
		constructor%T%s_ptr					=> scal_c_ptr%s_ptr
		call manager%get_cons_field_pointer_by_name(scal_c_ptr,vect_c_ptr,tens_c_ptr,'pressure')
		constructor%p%s_ptr					=> scal_c_ptr%s_ptr
		call manager%get_cons_field_pointer_by_name(scal_c_ptr,vect_c_ptr,tens_c_ptr,'specie_molar_concentration')
		constructor%Y%v_ptr					=> vect_c_ptr%v_ptr		
		call manager%get_cons_field_pointer_by_name(scal_c_ptr,vect_c_ptr,tens_c_ptr,'mixture_molar_concentration')
		constructor%mol_mix_conc%s_ptr		=> scal_c_ptr%s_ptr		
		call manager%get_cons_field_pointer_by_name(scal_c_ptr,vect_c_ptr,tens_c_ptr,'sensible_enthalpy')
		constructor%h_s%s_ptr				=> scal_c_ptr%s_ptr		
		
		call manager%get_cons_field_pointer_by_name(scal_c_ptr,vect_c_ptr,tens_c_ptr,'viscosity')
		constructor%nu%s_ptr			=> scal_c_ptr%s_ptr
		call manager%get_cons_field_pointer_by_name(scal_c_ptr,vect_c_ptr,tens_c_ptr,'thermal_conductivity')
		constructor%kappa%s_ptr			=> scal_c_ptr%s_ptr
		call manager%get_cons_field_pointer_by_name(scal_c_ptr,vect_c_ptr,tens_c_ptr,'diffusivity')
		constructor%D%v_ptr				=> vect_c_ptr%v_ptr		
		
		call manager%get_flow_field_pointer_by_name(scal_f_ptr,vect_f_ptr,tens_f_ptr,'velocity_flow')
		constructor%v_f%v_ptr			=> vect_f_ptr%v_ptr		
		
		write(var_name,'(A,I2.2)')		'density_production_particles', phase_number
		write(var_short_name,'(A,I2.2)')	'rho_prod_p', phase_number
		call manager%create_scalar_field(rho_prod_p(phase_number),	var_name,	var_short_name)
		constructor%rho_prod%s_ptr		=> rho_prod_p(phase_number)        
        
		write(var_name,'(A,I2.2)')		'energy_production_particles', phase_number
		write(var_short_name,'(A,I2.2)')	'E_f_prod_p', phase_number
		call manager%create_scalar_field(E_f_prod_p(phase_number),	var_name,	var_short_name)
		constructor%E_f_prod%s_ptr		=> E_f_prod_p(phase_number)

		write(var_name,'(A,I2.2)')		'velocity_production_particles', phase_number
		write(var_short_name,'(A,I2.2)')	'v_prod_p', phase_number
		call manager%create_vector_field(v_prod_p(phase_number),	var_name,	var_short_name,	'spatial')
		constructor%v_prod%v_ptr		=> v_prod_p(phase_number)
        
		write(var_name,'(A,I2.2)')		'concentration_production_particles', phase_number
		write(var_short_name,'(A,I2.2)')	'Y_prod_p', phase_number		
		call manager%create_vector_field(Y_prod_p(phase_number),	var_name,	var_short_name,	'chemical')
		constructor%Y_prod%v_ptr		=> Y_prod_p(phase_number)        

		constructor%phase_number        = phase_number

		constructor%mesh%mesh_ptr	    => manager%computational_mesh_pointer%mesh_ptr
		constructor%boundary%bc_ptr     => manager%boundary_conditions_pointer%bc_ptr
		constructor%thermo%thermo_ptr	=> manager%thermophysics%thermo_ptr
		constructor%chem%chem_ptr	    => manager%chemistry%chem_ptr
		constructor%domain			    = manager%domain
		
        constructor%g                   = manager%solver_options%get_grav_acc()
        
		system_command = 'mkdir data_save_particles'
		call system(system_command)		
		
		system_command = 'mkdir particle_grid'
		call system(system_command)	

		file_name = 'average_particles_data.dat'
		open(newunit = average_io_unit, file = file_name, status = 'replace', form = 'formatted')
        
		file_name = 'single_particle_trajectory.dat'
		open(newunit = trajectory_io_unit, file = file_name, status = 'replace', form = 'formatted')
        
	end function

	subroutine set_initial_distributions(this)
		class(lagrangian_particles_solver)	,intent(inout)	:: this

		integer		,dimension(3,2)	:: cons_inner_loop
		real(dp)	,dimension(3)	:: cell_size
		
		real(dp)	:: delta, part_number_max, a ,b, dimless_length
		integer		,dimension(3)	:: part_number
		real(dp)	,dimension(3)	:: coords
		real(dp)	,dimension(:,:)	,allocatable	:: lengths
        character(len=5)	,dimension(3)		    :: axis_names
        
		integer	:: dimensions
		integer :: i,j,k, part, part_x, part_y, part_z, dim, particle, unit, unit_h
        integer     ,dimension(3)   :: bins
		integer		,dimension(1)   :: gamma_hist

		integer, dimension(3)	:: cell, initial_cell
		logical	:: out_flag
		
		dimensions		= this%domain%get_domain_dimensions()
		cons_inner_loop	= this%domain%get_local_inner_cells_bounds()		

		allocate(lengths(dimensions,2))
		lengths			= this%domain%get_domain_lengths()
		
		cell_size		= this%mesh%mesh_ptr%get_cell_edges_length()
        axis_names      = this%domain%get_axis_names()
        !# Single particle test
        
		part_number = 1
		
		tracked_particle = 1
        
		this%particles_number = part_number(1)*part_number(2)*part_number(3)
		
		allocate(this%particles(this%particles_number))	
		
		do part = 1, this%particles_number
			this%particles(part)%outside_domain	= .true.
			this%particles(part)%coords			= 0.0_dp
			this%particles(part)%velocity		= 0.0_dp
			this%particles(part)%dm				= 0.0_dp 
            
            this%particles(part)%evaporating    = this%particles_params%evaporating
            this%particles(part)%heating        = this%particles_params%heating
            this%particles(part)%inertial       = this%particles_params%inertial
        end do 

        this%particles(1)%coords(1)	        = 0.5_dp * (lengths(1,2) + cell_size(1))
        this%particles(1)%coords(2)	        = 0.5_dp * (lengths(2,2) + cell_size(1))
        this%particles(1)%coords(3)	        = 0.5_dp * (cell_size(1)) !0.5_dp * (lengths(3,2) + cell_size(1))
        
        this%particles(1)%coords_prev	    = this%particles(1)%coords
		this%particles(1)%coords0		    = this%particles(1)%coords       

        this%particles(1)%outside_domain	= .false. 
        this%particles(1)%temperature		= 300.0_dp
        
		this%particles(1)%mass	= Pi*this%particles_params%diameter**3 / 6.0_dp * this%particles_params%material_density
        continue
        
        !# Two particle streams setup
        
  !      this%particles_number = 1000
		!allocate(this%particles(this%particles_number))	
  !      
  !      release_time = 2.5e-05
  !      particles_velocity = 20.0_dp
  !      streams_distance = 2e-03
  !      
  !      delta = release_time * particles_velocity
  !      
  !      part_number(1) = 25
  !      
  !      do part = 1,  this%particles_number
		!	this%particles(part)%outside_domain		= .true.
  !          this%particles(part)%coords(1)			= 0.0
		!	this%particles(part)%velocity(1)			= 0.0
  !      end do
  !      
		!do part = 1, part_number(1)
		!	this%particles(part)%outside_domain		= .false.
		!	this%particles(part)%coords(1)			= delta * (part + 1)
  !          this%particles(part)%coords(2)			= lengths(2,2) / 2 + streams_distance / 2
		!	this%particles(part)%velocity(1)			= 1.0_dp
		!	this%particles(part)%dm					= 0.0_dp
  !      	this%particles(part)%temperature			= 300.0_dp
		!	this%particles(part)%mass				= Pi*this%particles_params%diameter**3 / 6.0_dp * this%particles_params%material_density
  !      end do 
  !      
  !      do part = part_number(1)+1, 2*part_number(1)
		!	this%particles(part)%outside_domain		= .false.
		!	this%particles(part)%coords(1)			= delta * (part + 1 - part_number(1))
  !          this%particles(part)%coords(2)			= lengths(2,2) / 2 - streams_distance / 2
		!	this%particles(part)%velocity(1)			= 1.0_dp
		!	this%particles(part)%dm					= 0.0_dp
  !          this%particles(part)%temperature			= 300.0_dp
		!	this%particles(part)%mass				= Pi*this%particles_params%diameter**3 / 6.0_dp * this%particles_params%material_density
  !      end do 
                
        
        lengths			                    = this%domain%get_domain_lengths()
        
        this%part_grid%grid_cell_size       = 2*cell_size(1)
        do dim = 1, dimensions
        this%part_grid%part_grid_size(1)    = nint((lengths(1,2)-lengths(1,1))/(this%part_grid%grid_cell_size))	!8
        this%part_grid%part_grid_size(2)    = nint((lengths(2,2)-lengths(2,1))/(this%part_grid%grid_cell_size))	!8
            this%part_grid%part_grid_size(3) = nint((lengths(3,2)-lengths(3,1))/(this%part_grid%grid_cell_size))
        
            bins(dim)	= this%part_grid%part_grid_size(dim) 
        end do 
        
        allocate(this%part_grid%coords(bins(1), bins(2), bins(3), 3))
        allocate(this%part_grid%velocity(bins(1), bins(2), bins(3), 3))
        allocate(this%part_grid%density(bins(1), bins(2), bins(3)))
        allocate(this%part_grid%temperature(bins(1), bins(2), bins(3)))
		
        this%part_grid%coords		= 0.0_dp
        this%part_grid%velocity		= 0.0_dp
        this%part_grid%density		= 0.0_dp
        this%part_grid%temperature	= 0.0_dp
        
        do i = 1, bins(1)
        do j = 1, bins(2)
        do k = 1, bins(3)
            
            do dim = 1, dimensions
			    this%part_grid%coords(i,j,k,dim) = lengths(dim,1) + (i * I_m(dim,1) + j * I_m(dim,2) + k * I_m(dim,3) - 0.5_dp) * this%part_grid%grid_cell_size
            end do
            
            do particle = 1, this%particles_number
                
                if ((abs(this%particles(particle)%coords(1) - this%part_grid%coords(i,j,k,1)) < 0.5_dp * this%part_grid%grid_cell_size).and. &
                    (abs(this%particles(particle)%coords(2) - this%part_grid%coords(i,j,k,2)) < 0.5_dp * this%part_grid%grid_cell_size).and. &
                    (abs(this%particles(particle)%coords(3) - this%part_grid%coords(i,j,k,3)) < 0.5_dp * this%part_grid%grid_cell_size)) then
					this%part_grid%density(i,j,k)		= this%part_grid%density(i,j,k) + 1
                    this%part_grid%temperature(i,j,k)	= this%part_grid%temperature(i,j,k) + this%particles(particle)%temperature
                    this%part_grid%velocity(i,j,k,1)	= this%part_grid%velocity(i,j,k,1) + this%particles(particle)%velocity(1)
                    this%part_grid%velocity(i,j,k,2)	= this%part_grid%velocity(i,j,k,2) + this%particles(particle)%velocity(2)
                    this%part_grid%velocity(i,j,k,3)	= this%part_grid%velocity(i,j,k,3) + this%particles(particle)%velocity(3)
                end if
            end do
        end do
        end do
        end do

        av_header = 'VARIABLES="time" "av_path" "av_velo"'
        trajectory_header = 'VARIABLES="time"'
        grid_header = 'VARIABLES='
        scatter_header = 'VARIABLES='
        do dim = 1, dimensions
            grid_header = trim(grid_header) // '"' // trim(axis_names(dim)) // '"'
            scatter_header = trim(scatter_header) // '"' // trim(axis_names(dim)) // '"'
            av_header =  trim(av_header) // '"av_' // trim(axis_names(dim)) // '"'
            trajectory_header = trim(trajectory_header) // '"' // trim(axis_names(dim)) // '"'
        end do
        do dim = 1, dimensions
            grid_header = trim(grid_header) // '"v(' // trim(axis_names(dim)) // ')"'
            scatter_header = trim(scatter_header) // '"v(' // trim(axis_names(dim)) // ')"'
            trajectory_header = trim(trajectory_header) // '"v(' // trim(axis_names(dim)) // ')"'
        end do
        
        do dim = 1, dimensions
            scatter_header = trim(scatter_header) // '"p_orient(' // trim(axis_names(dim)) // ')"'
        end do   
        do dim = 1, dimensions
            scatter_header = trim(scatter_header) // '"F_a(' // trim(axis_names(dim)) // ')"'
        end do   
        do dim = 1, dimensions
            scatter_header = trim(scatter_header) // '"F_St(' // trim(axis_names(dim)) // ')"'
        end do         
        
        grid_header = trim(grid_header) // '"rho_d" "T_d"'
        scatter_header = trim(scatter_header) // '"T_d" "m_d"'
        av_header = trim(av_header) // '"av_temp" "av_diameter"'
        
        write (average_io_unit,'(A)')	trim(av_header)

        write (trajectory_io_unit,'(A)')	trim(trajectory_header)
        
        continue 
        
	end subroutine
	
	
	subroutine particles_solve(this,time_step)

		class(lagrangian_particles_solver)	,intent(inout)	:: this
		real(dp)				,intent(in)		:: time_step
		real(dp)								:: particle_time_step
		integer									:: particle_iterations
		
		integer		,dimension(3,2)	:: cons_inner_loop
		real(dp)	,dimension(3)	:: cell_size
		real(dp)					:: cell_volume
		real(dp)	,dimension(3)	:: gas_velocity, relative_velocity, old_coords, old_velocity, F_a, particle_acceleration
		real(dp)					:: abs_relative_velocity, dot_prod1, dot_prod2, dot_prod3, dot_prod4, dot_prod5
		real(dp)					:: lower_face, higher_face, a, b
		real(dp)					:: Re_p, Nu_p, Sc_p, Sh_p, Pr_p, C_drag, A_pc, A_ps, alpha_p, beta_p, D_air_water, Bm
		real(dp)					:: H_mass, H_heat, H_v, H_l, dH_v_dT
		real(dp)					:: X_vap, Y_vap, W_ratio, dY_dT
		real(dp)					:: specie_enthalpy_gas, specie_enthalpy_liquid, mixture_cp, h_s_Tref
		real(dp)					:: dmp, Q,	Tg_new, Tp_new, Tg_old, Tp_old, rhog_old, rhog_new, M_gas_old, M_gas_new, m_liq, mol_mix_w, Hg_old, Hg_new, Hp_old, Hp_new
		real(dp)					:: delta_H
		real(dp)	,dimension(:)	, allocatable	:: Yg_old, Yg_new
        real(dp)	,dimension(:,:)	,allocatable	:: lengths
		real(dp)					:: DTOP, DTOG, DTGOG, DTGOP, AGHRHO, DADYDTHVHL, DADYDTHV, DAHVHLDY, DADYHV
		real(dp)	,dimension(2)	:: A_col, B_col, D_vec
		
		integer						:: particles_in_cell
		character(len=100)			:: file_path, file_name
        character(len=20)           :: output_fmt
		
		real(dp)                    :: particle_diameter, evaporation_rate, velocity_2, specie_enthalpy
		real(dp)					:: average_path, average_velocity, average_diameter, average_temperature, particles_inside
		real(dp)	,dimension(3)	:: average_coordinates

        real(dp)					:: d_ij
        real(dp)					:: velocity_mag, delta
		real(dp)	,dimension(2)	:: r_ij, v_ij, abs_r1, abs_r2, wall_distance

		integer	:: dimensions, species_number
		integer, dimension(3)	:: cell, initial_cell
		
        integer :: particle_material_index
		integer :: C7H16_index
		integer :: i,j,k,dim, part, part_iter, spec, part1, part2
		logical	:: out_flag
		
		real(dp)	,save	:: time = 0.0_dp, CFL_p
		integer		,save	:: output_counter = 0
        integer		,save	:: av_output_counter = 0
        integer		,save	:: traj_output_counter = 0
		integer		,save	:: particle_release_counter = 0

		integer				:: lagrangian_particles_io_unit
		
        lengths			= this%domain%get_domain_lengths()

		dimensions		= this%domain%get_domain_dimensions()
		species_number	= this%chem%chem_ptr%species_number
		
		allocate(Yg_old(species_number),Yg_new(species_number))
		
		cons_inner_loop	= this%domain%get_local_inner_cells_bounds()		

		cell_size		= this%mesh%mesh_ptr%get_cell_edges_length()
		cell_volume		= this%mesh%mesh_ptr%get_cell_volume()
		
		particle_material_index		= this%chem%chem_ptr%get_chemical_specie_index(this%particles_params%material)
		
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
					particle	=> this%particles_params	, &
                    chem        => this%chem%chem_ptr   , &
					mesh		=> this%mesh%mesh_ptr	, &
					bc			=> this%boundary%bc_ptr)

		!!$omp parallel default(none)  private(i,j,k,dim,F_stokes,Q_stokes,particle_diameter,evaporation_rate,velabs,temp_cr) , &
		!!$omp& firstprivate(this)	,&
		!!$omp& shared(T,T_d,T_d_int,E_f_prod,rho,rho_d,mass_d,numdens_d,v_d,v_d_int,v,v_prod,Y_prod,nu,kappa,foam_marker,time_boil,particle,Nusselt,particle_material_index,time_step,mesh,bc,cons_inner_loop,dimensions,marker1,specie_enthalpy)
		!!$omp do collapse(3) schedule(guided)					
			
		E_f_prod%cells = 0.0_dp			
		Y_prod%pr(particle_material_index)%cells = 0.0_dp	
		
		do dim = 1, dimensions
			v_prod%pr(dim)%cells = 0.0_dp
		end do

		average_temperature = 0.0_dp
		average_diameter	= 0.0_dp
		average_path		= 0.0_dp
		average_velocity	= 0.0_dp
        average_coordinates	= 0.0_dp
		particles_inside    = 0
		
        this%part_grid%velocity		= 0.0_dp
        this%part_grid%density		= 0.0_dp
        this%part_grid%temperature	= 0.0_dp        
                
        !# Dynamics
		do part = 1, this%particles_number

			if (.not.this%particles(part)%outside_domain) then

				velocity_mag = 0.0_dp
		
				average_diameter	= average_diameter		+ (6.0_dp * this%particles(part)%mass / Pi /  particle%material_density) ** (1.0_dp / 3.0_dp)
				average_temperature	= average_temperature	+ this%particles(part)%temperature
                do dim = 1, dimensions
                    velocity_mag	= velocity_mag + this%particles(part)%velocity(dim)**2.0_dp
                end do
                
                do dim = 1, dimensions
					average_coordinates(dim)    = average_coordinates(dim)  + this%particles(part)%coords(dim)
                    average_path		        = average_path		        + (this%particles(part)%coords(dim)-this%particles(part)%coords0(dim)) ** 2.0_dp 
                end do
                
				average_velocity	= average_velocity  + velocity_mag

				particles_inside    = particles_inside + 1
				
				!# Move particle, account momentum transfer

				!particle_time_step = time_step
                CFL_p = 0.0
                if (velocity_mag > 1e-10_dp) then
				    do dim = 1, dimensions
                        !if (abs(this%particles(part)%velocity(dim)/(0.5*this%particles_params%diameter)) > CFL_p) CFL_p = abs(this%particles(part)%velocity(dim)/(0.5*this%particles_params%diameter))
                        if (abs(this%particles(part)%velocity(dim))/cell_size(dim) > CFL_p) CFL_p = abs(this%particles(part)%velocity(dim))/cell_size(dim)
                    end do
                else
                    CFL_p = 1.0_dp/0.9_dp
                end if
                CFL_p = time_step*CFL_p

                particle_time_step = time_step / ceiling(0.9*CFL_p) / 10.0_dp			!Time step for particles
                particle_iterations	= time_step/ particle_time_step
                
				F_a =	0.0_dp			
			
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
			
                    old_velocity	= this%particles(part)%velocity
                    old_coords		= this%particles(part)%coords

                    this%particles(part)%coords_prev    = this%particles(part)%coords
                    
					particle_diameter = (6.0_dp * this%particles(part)%mass / Pi /  particle%material_density) ** (1.0_dp / 3.0_dp)
					
					cell = this%get_particle_cell(this%particles(part)%coords,out_flag)
					i = cell(1)
					j = cell(2)
					k = cell(3)				
				
                    !# particle moved to a different cell during the sub time step
					if ((this%particles(part)%cell(1) /= cell(1)).and.(this%particles(part)%cell(2) /= cell(2)).and.(this%particles(part)%cell(3) /= cell(3))) then
						Tg_old		= T%cells(i,j,k)
						rhog_old	= rho%cells(i,j,k)
						M_gas_old	= rhog_old * cell_volume 
						Hg_old		= h_s%cells(i,j,k) * M_gas_old
						do spec = 1, species_number
							Yg_old(spec)	= Y%pr(spec)%cells(i,j,k)
						end do
					
						E_f_prod%cells(i,j,k) = 0.0_dp
						Y_prod%pr(particle_material_index)%cells(i,j,k) = 0.0_dp					
					
						this%particles(part)%cell = cell
					end if

				!# Interpolate gas velocity to the particle location using linear interpolation
					gas_velocity = 0.0_dp
					do dim = 1, dimensions
 						lower_face	= mesh%mesh(dim,i,j,k) - 0.5_dp*cell_size(dim)
						higher_face	= mesh%mesh(dim,i,j,k) + 0.5_dp*cell_size(dim)
					
                        a = (v_f%pr(dim)%cells(dim,i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3)) - v_f%pr(dim)%cells(dim,i,j,k)) / cell_size(dim)
						b = v_f%pr(dim)%cells(dim,i,j,k)

						gas_velocity(dim) = a * (this%particles(part)%coords(dim) - lower_face) + b
					end do
				
                    !# Calculate particle-gas relative velocity
					abs_relative_velocity = 0.0_dp
					do dim = 1, dimensions
						relative_velocity(dim)	= this%particles(part)%velocity(dim) - gas_velocity(dim)
						abs_relative_velocity	= abs_relative_velocity + relative_velocity(dim)*relative_velocity(dim)
					end do
				
					abs_relative_velocity = sqrt(abs_relative_velocity)
				
					!# Calculate particle mass and cross section
					m_liq	= this%particles(part)%mass
					A_pc	= Pi*particle_diameter**2 / 4.0_dp
					A_ps	= Pi*particle_diameter**2

					!# Calculate particle Reynolds number
					Re_p	= rhog_old * abs_relative_velocity * particle_diameter / nu%cells(i,j,k)
				
					if ( abs_relative_velocity > 1.0e-10_dp) then	
					!# Calculate drag coefficient for spherical particle
						if (Re_p < 1e-010) then
							C_drag = 100
						elseif (Re_p < 1.0_dp) then
								C_drag = 24.0_dp / Re_p
						elseif (Re_p < 1000.0_dp) then
								C_drag = 24.0_dp * ( 0.85_dp + 0.15 * Re_p**0.687) / Re_p
						elseif (Re_p >= 1000.0_dp) then
								C_drag = 0.44_dp
						end if
				
						particles_in_cell = this%count_particles_in_cell(cell)
				
					!# Calculate new particle coordinates and velocities
						!alpha_p	= M_gas_old / m_liq
						!beta_p	= 0.5_dp * rhog_old * C_drag * A_p * ( 1/m_liq + 1/M_gas_old) * abs_relative_velocity
					
                        beta_p = 0.5_dp * rhog_old * C_drag * A_pc * abs_relative_velocity / m_liq
                        
						do dim = 1, dimensions
                            if (beta_p > 1e-10_dp) then							
                                this%particles(part)%velocity(dim)	= gas_velocity(dim) + (this%particles(part)%velocity(dim) - gas_velocity(dim))*exp(-beta_p * particle_time_step) 
                                
                                this%particles(part)%velocity(dim)	= this%particles(part)%velocity(dim)	+ this%g(dim)/beta_p * ( 1.0_dp - exp(-beta_p * particle_time_step))
                            else
                                this%particles(part)%velocity(dim)	= gas_velocity(dim) + this%g(dim) * particle_time_step
                            end if

				            if(this%particles(part)%velocity(dim) > 0.0) then
                                wall_distance(dim) = (mesh%mesh(dim,(1-I_m(dim,1))+I_m(dim,1)*cons_inner_loop(dim,2),(1-I_m(dim,2))+I_m(dim,2)*cons_inner_loop(dim,2),(1-I_m(dim,3))+I_m(dim,3)*cons_inner_loop(dim,2)) + 0.5_dp*cell_size(dim) )
								wall_distance(dim) = wall_distance(dim) - this%particles(part)%coords(dim)
                            else
								wall_distance(dim) = (mesh%mesh(dim,(1-I_m(dim,1))+I_m(dim,1)*cons_inner_loop(dim,1),(1-I_m(dim,2))+I_m(dim,2)*cons_inner_loop(dim,1),(1-I_m(dim,3))+I_m(dim,3)*cons_inner_loop(dim,1)) - 0.5_dp*cell_size(dim) ) 
                                wall_distance(dim) = wall_distance(dim) - this%particles(part)%coords(dim)
                            end if
                            
                            if ( abs(0.5_dp * particle_time_step*(this%particles(part)%velocity(dim) + old_velocity(dim))) < abs(wall_distance(dim))) then
	                            this%particles(part)%coords(dim)		= this%particles(part)%coords_prev(dim) + 0.5_dp * particle_time_step*(this%particles(part)%velocity(dim) + old_velocity(dim))
                            else
                                !# Closed boundaries    
	                            this%particles(part)%coords(dim)		= this%particles(part)%coords_prev(dim) + 2.0_dp * wall_distance(dim) - 0.5_dp * particle_time_step*(this%particles(part)%velocity(dim) + old_velocity(dim)) 
                                this%particles(part)%velocity(dim)	    = -this%particles(part)%velocity(dim)
                                !print *, "reflection!", particle, dim, this%particles(part)%coords(dim), old_coords(dim)
                                !pause
                            end if
                            
                            F_a(dim) = 1.0_dp / (cell_volume) * m_liq * ( - this%g(dim) - (this%particles(part)%velocity(dim) - old_velocity(dim)) / particle_time_step - this%particles(part)%dm * (this%particles(part)%velocity(dim) - gas_velocity(dim)))
						end do 
					end if
					
					cell = this%get_particle_cell(this%particles(part)%coords,out_flag)
					if (out_flag) then
						this%particles(part)%outside_domain = out_flag
						exit
					end if 
			
					do dim = 1, dimensions
                        
						lower_face	= mesh%mesh(dim,i,j,k) - 0.5_dp*cell_size(dim)
						higher_face	= mesh%mesh(dim,i,j,k) + 0.5_dp*cell_size(dim)
						
						a = this%particles(part)%coords(dim) - lower_face

						b = higher_face - this%particles(part)%coords(dim) 
					
                        if (this%particles(part)%inertial) then
						    v_prod%pr(dim)%cells(i,j,k)	= v_prod%pr(dim)%cells(i,j,k) - F_a(dim) 
						end if
                    end do
                   
                    if (this%particles(part)%heating) then
                        
                        !# Heat and evaporate particle, account energy transfer
                        !# Temperature dependence of latent heat of water
					        
			            dH_v_dT = -2316.0_dp
                        H_v	    = 3.133e+06_dp + dH_v_dT * Tp_old 
                        
					    mol_mix_w = 0.0_dp
					    do spec = 1, species_number
						    mol_mix_w = mol_mix_w + Yg_old(spec)/this%thermo%thermo_ptr%molar_masses(spec)
					    end do
					    mol_mix_w = 1.0_dp / mol_mix_w
			
					    W_ratio	= mol_mix_w*(1.0_dp - Yg_old(particle_material_index)) / this%thermo%thermo_ptr%molar_masses(particle_material_index)
			
					    X_vap = min(1.0_dp,  exp(H_v * this%thermo%thermo_ptr%molar_masses(particle_material_index) / r_gase_J * (1.0_dp / particle%material_boiling_temperature - 1.0_dp / this%particles(part)%temperature)))
					    Y_vap = X_vap / (X_vap * (1.0_dp - W_ratio) + W_ratio)
			
                        !# dY_dT with an account of hv(T) from FDS part.90 ln 3940. 
					    dY_dT = (W_ratio / (X_vap*(1.0_dp - W_ratio) + W_ratio)**2) * &
                                (H_v * X_vap * this%thermo%thermo_ptr%molar_masses(particle_material_index) / r_gase_J / this%particles(part)%temperature ** 2 + &
                                (1.0_dp / particle%material_boiling_temperature - 1.0_dp / this%particles(part)%temperature)*(dH_v_dT) * this%thermo%thermo_ptr%molar_masses(particle_material_index) / r_gase_J)
				
					    Pr_p = 0.7_dp
					    Sc_p = 0.6_dp
			
                        !# BM is the Spalding mass transfer number
                        Bm = (Y_vap - Yg_old(particle_material_index)) / (1.0_dp - Y_vap)
                    
                        !# FDS TRG 6ed pg. 99 (func.f90 particle_H_MASS_H_HEAT_GAS)
					    Nu_p = log(1.0_dp + Bm)/Bm * (2.0_dp + 0.6_dp * sqrt(Re_p) * Pr_p ** (1.0_dp / 3.0_dp))
					    Sh_p = log(1.0_dp + Bm)/Bm * (2.0_dp + 0.6_dp * sqrt(Re_p) * Sc_p ** (1.0_dp / 3.0_dp))
			
					    D_air_water = D%pr(particle_material_index)%cells(i,j,k)

    !					H_g	= 814.814366872002 *1e-04/ particle_diameter       ! FDS Nu_p * kappa%cells(i,j,k)	/ particle_diameter !814.814366872002 ! FDS
    !					H_m	= 1.35020265065071 *1e-04/ particle_diameter       ! FDS Sh_p * D_air_water / particle_diameter          !1.35020265065071 ! FDS
                    
                        !# FDS TRG 6ed pg. 99 (func.f90 particle_H_MASS_H_HEAT_GAS)
					    H_heat	= Nu_p * kappa%cells(i,j,k)	/ particle_diameter ! * 1.0_dp                                                              ! Nu_p * kappa%cells(i,j,k)	/ particle_diameter   !814.814366872002 ! FDS
					    H_mass	= Sh_p * D_air_water * Bm / (Y_vap - Yg_old(particle_material_index)) / particle_diameter !/ 50.0_dp                     ! Sh_p * D_air_water / particle_diameter          !1.35020265065071 ! FDS
                    
                    
					    specie_enthalpy_gas		= (this%thermo%thermo_ptr%calculate_specie_cp(Tg_old,particle_material_index))*Tg_old / this%thermo%thermo_ptr%molar_masses(particle_material_index)
					    specie_enthalpy_liquid	= (this%thermo%thermo_ptr%calculate_specie_cp(this%particles(part)%temperature,particle_material_index))*this%particles(part)%temperature / this%thermo%thermo_ptr%molar_masses(particle_material_index)

                        h_s_Tref			        = this%thermo%thermo_ptr%calculate_specie_enthalpy(T_ref, particle_material_index)
                        specie_enthalpy_gas		    = (this%thermo%thermo_ptr%calculate_specie_enthalpy(Tg_old, particle_material_index) - h_s_Tref) / this%thermo%thermo_ptr%molar_masses(particle_material_index)
                        specie_enthalpy_liquid		= (this%thermo%thermo_ptr%calculate_specie_enthalpy(this%particles(part)%temperature, particle_material_index) - h_s_Tref) / this%thermo%thermo_ptr%molar_masses(particle_material_index)
                    
					    mixture_cp = 0.0
					    do spec = 1, species_number
						    mixture_cp		= mixture_cp + this%thermo%thermo_ptr%calculate_specie_cp(Tg_old,spec)	* Yg_old(spec) / this%thermo%thermo_ptr%molar_masses(spec)
                        end do			

    				    AGHRHO      = A_ps * H_mass * (rhog_old) / ( 1.0_dp + 0.5_dp * particle_time_step * A_ps * H_mass * (1.0_dp - Yg_old(particle_material_index)) / cell_volume / rhog_old)
                                        
                        DTOG        = particle_time_step / (mixture_cp * M_gas_old)
                        DTGOG       = 0.5_dp * DTOG * A_ps * H_heat                    
                        DADYDTHVHL  = 0.5_dp * DTOG * AGHRHO	* (specie_enthalpy_liquid - specie_enthalpy_gas) * dY_dT
                        DAHVHLDY    = DTOG * AGHRHO * (specie_enthalpy_liquid - specie_enthalpy_gas) * (Y_vap - Yg_old(particle_material_index))
                    
                        DTOP        = particle_time_step / (particle%material_heat_capacity * m_liq)
                        DTGOP       = 0.5_dp * DTOP * A_ps * H_heat
                        DADYDTHV    = 0.5_dp * DTOP * AGHRHO  * dY_dT * H_v
                        DADYHV      = DTOP*AGHRHO*H_v*(Y_vap - Yg_old(particle_material_index))
               
					    A_col(1)	= 1.0_dp + DTGOG
					    B_col(1)	= - (DTGOG + DADYDTHVHL)
					    A_col(2)	= -DTGOP
					    B_col(2)	= 1.0_dp + DTGOP + DADYDTHV
					    D_vec(1)	= (1.0_dp - DTGOG) * Tg_old + (DTGOG - DADYDTHVHL) * this%particles(part)%temperature + DAHVHLDY
					    D_vec(2)	= DTGOP * Tg_old + (1.0_dp - DTGOP + DADYDTHV) * this%particles(part)%temperature - DADYHV
			
					    Tp_new	= -(A_col(2) * D_vec(1) - A_col(1) * D_vec(2))/(A_col(1) * B_col(2) - A_col(2) * B_col(1))
					    Tg_new	= (D_vec(1) - B_col(1) * Tp_new) / A_col(1)
			
					    Q	= particle_time_step * A_ps * H_heat * 0.5_dp * (Tg_old + Tg_new - Tp_old - Tp_new)
					    if (Q > this%particles(part)%mass * H_v) dmp = this%particles(part)%mass
                    
                        if (this%particles(part)%evaporating) then
					        dmp = max(0.0_dp, min(this%particles(part)%mass, particle_time_step * AGHRHO * ( Y_vap -  Yg_old(particle_material_index) + 0.5_dp * dY_dT * (Tp_new - Tp_old))))
                        else 
                            dmp = 0.0_dp
                        end if
			
					    if ( dmp < this%particles(part)%mass) then
                            !print *, 'Q', Q
                            !print *, 'Tp_new', Tp_new
                            !print *, 'Tg_new', Tg_new
                            !print *, 'dmp', dmp
                            !print *, 'H_v', H_v
                            !print *, 'Q-dmp*H_v', Q - dmp * H_v
                            !print *, 'Cpd*(m-dm)', (particle%material_heat_capacity * ( this%particles(part)%mass - dmp))
                        
						    Tp_new = Tp_old + (Q - dmp * H_v) / (particle%material_heat_capacity * ( this%particles(part)%mass - dmp))
                        
                            !print *, 'Tp_new_adj', Tp_new

						    if ( Tp_new > particle%material_boiling_temperature) then
							    dmp = min(this%particles(part)%mass, (Q - this%particles(part)%mass * particle%material_heat_capacity * (particle%material_boiling_temperature - Tp_old))/(H_v - particle%material_heat_capacity * (particle%material_boiling_temperature - Tp_old)))
							    if ( dmp == this%particles(part)%mass) then
								    Q = dmp * H_v
							    end if
							    Tp_new = particle%material_boiling_temperature
						    end if
					    else
						    Q		= dmp * H_v
						    Tp_new	= particle%material_boiling_temperature
						    this%particles(part)%outside_domain = .true.
					    end if
			
					    this%particles(part)%mass			= this%particles(part)%mass - dmp
					    this%particles(part)%temperature	= Tp_new
					    this%particles(part)%dm				= this%particles(part)%dm + dmp / particle_time_step
			
					    M_gas_new	= M_gas_old + dmp
					    Yg_new		= Yg_old*M_gas_old/M_gas_new
					    Yg_new(particle_material_index)		=	Yg_new(particle_material_index) + dmp/M_gas_new
			
					    Hg_new		= Hg_old	+ (Hp_old - this%particles(part)%mass*particle%material_heat_capacity*Tp_new)
    !					Tg_new		= Tg_old	+ (Hg_new - mixture_cp * Tg_old * M_gas_new) / M_gas_new / mixture_cp		!# Assuming cp = const
			
                        Tg_new		= Tg_old	+ (Hg_new - Hg_old) / M_gas_new / mixture_cp
                    
					    rhog_new	= M_gas_new / cell_volume
			
					    delta_H		= ((this%thermo%thermo_ptr%calculate_specie_cp(Tp_old,particle_material_index))*Tp_old - (this%thermo%thermo_ptr%calculate_specie_cp(Tg_old,particle_material_index))*Tg_old) / this%thermo%thermo_ptr%molar_masses(particle_material_index)
			
					    E_f_prod%cells(i,j,k)	=	E_f_prod%cells(i,j,k) + (W_ratio * dmp / M_gas_old + (dmp * delta_H - Q) / Hg_old) / time_step
					
					    do dim = 1, dimensions
						    dmp = dmp / cell_size(dim)
					    end do
					
					    Y_prod%pr(particle_material_index)%cells(i,j,k)	= Y_prod%pr(particle_material_index)%cells(i,j,k) + dmp / particle_time_step

					    Yg_old		= Yg_new
					    Hg_old		= Hg_new
					    rhog_old	= rhog_new
					    Tg_old		= Tg_new
                    end if
				end do
			else
				!print *, 'particle', part, ' with coodinates', this%particles(part)%coords, '  is outside the domain'
				!pause
			end if
        end do
        
        !# Collisions
        
        do part1 = 1, this%particles_number
            this%particles(part1)%coords_prev	= this%particles(part1)%coords
        end do
        
        do i = 1, this%part_grid%part_grid_size(1)
        do j = 1, this%part_grid%part_grid_size(2)
        do k = 1, this%part_grid%part_grid_size(3)    
            do part = 1, this%particles_number
				if ((abs(this%particles(part)%coords(1) - this%part_grid%coords(i,j,k,1)) < 0.5_dp * this%part_grid%grid_cell_size).and. &
					(abs(this%particles(part)%coords(2) - this%part_grid%coords(i,j,k,2)) < 0.5_dp * this%part_grid%grid_cell_size).and. &
					(abs(this%particles(part)%coords(3) - this%part_grid%coords(i,j,k,3)) < 0.5_dp * this%part_grid%grid_cell_size)) then
					this%part_grid%density(i,j,k)         = this%part_grid%density(i,j,k) + 1
					this%part_grid%temperature(i,j,k)     = this%part_grid%temperature(i,j,k) + this%particles(part)%temperature
					this%part_grid%velocity(i,j,k,1)      = this%part_grid%velocity(i,j,k,1) + this%particles(part)%velocity(1)
					this%part_grid%velocity(i,j,k,2)      = this%part_grid%velocity(i,j,k,2) + this%particles(part)%velocity(2)
                    this%part_grid%velocity(i,j,k,3)      = this%part_grid%velocity(i,j,k,3) + this%particles(part)%velocity(3)
				end if
            end do
        end do
        end do
        end do
        
		if (((time/(1e-04) >= (av_output_counter)).or.(time==0.0_dp)).and.(particles_inside /= 0)) then
            write (output_fmt,'("(",I1,"E14.6)")') dimensions + 5
		    write (average_io_unit,output_fmt)	time, average_path/particles_inside     , &
                                                average_velocity/particles_inside       , &
                                                average_coordinates/particles_inside    , &
                                                average_temperature/particles_inside    , &
                                                average_diameter/particles_inside
            
            av_output_counter = av_output_counter + 1
		end if			
			
		if (((time/(50.0e-06) >= (traj_output_counter)).or.(time==0.0_dp)).and.(.not.this%particles(tracked_particle)%outside_domain)) then
            
            write (output_fmt,'("(",I1,"E14.6)")') 2*dimensions + 1
		    write (trajectory_io_unit,output_fmt)	time                        , & 
                                                    this%particles(tracked_particle)%coords   , &
                                                    this%particles(tracked_particle)%velocity               

            traj_output_counter = traj_output_counter + 1
		end if		
		
		if ((time/(1.0e-03) >= (output_counter)).or.(time==0.0_dp)) then
			write(file_path,'(I6.6,A)')  int(time*1e03),'ms'
			file_name = 'data_save_particles' // trim(fold_sep) // 'particles_' // trim(file_path) // '.plt'
			open(newunit = lagrangian_particles_io_unit, file = file_name, status = 'replace', form = 'formatted')
            write (lagrangian_particles_io_unit,'(A)') trim(scatter_header)
			write (output_fmt,'("(",I2,"E14.6)")') 2*dimensions + 2
			do part = 1, this%particles_number
				if(.not.this%particles(part)%outside_domain) then
					write (lagrangian_particles_io_unit,output_fmt)	this%particles(part)%coords(1:dimensions), this%particles(part)%velocity(1:dimensions), this%particles(part)%temperature, this%particles(part)%mass			
				end if
			end do

			close(lagrangian_particles_io_unit)
            
            write(file_path,'(I6.6,A)')  int(time*1e03),'ms'
			file_name =  'particle_grid' // trim(fold_sep) // 'particles_' // trim(file_path) // '_grid.plt'
			open(newunit = lagrangian_particles_io_unit, file = file_name, status = 'replace', form = 'formatted')
			write (lagrangian_particles_io_unit,'(A)') trim(grid_header)
            write (output_fmt,'("(",I1,"E14.6)")') 2*dimensions + 2
            do i = 1, this%part_grid%part_grid_size(1)
            do j = 1, this%part_grid%part_grid_size(2)
            do k = 1, this%part_grid%part_grid_size(3)
				write (lagrangian_particles_io_unit,output_fmt)	this%part_grid%coords(i,j,k,:), this%part_grid%velocity(i,j,k,:), this%part_grid%density(i,j,k), this%part_grid%temperature(i,j,k)	
            end do
            end do
            end do 
            
			close(lagrangian_particles_io_unit)
            
            output_counter = output_counter + 1
		end if
		
		!if ((time/release_time >= 1.0_dp*(particle_release_counter)).or.(time==0.0_dp)) then
		!	particle_release_counter = particle_release_counter + 1
		!	call this%release_particle()
		!end if


		time = time + time_step

		
		!!$omp end do nowait
		!!$omp end parallel		
		
		end associate
	end subroutine

	function get_particle_cell(this,particle_coords, out_flag)
		
		class(lagrangian_particles_solver)			,intent(inout)					:: this
		real(dp)		,dimension(3)				,intent(in)						:: particle_coords
		logical										,intent(inout)	,optional		:: out_flag
		
		integer			,dimension(3)	:: get_particle_cell
		integer		,dimension(3,2)	:: cons_inner_loop
		real(dp)	,dimension(3)	:: cell_size
		
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
				if (( 0.5_dp*(mesh%mesh(3,1,1,k-1) + mesh%mesh(3,1,1,k)) <= particle_coords(3)).and. &
					( 0.5_dp*(mesh%mesh(3,1,1,k+1) + mesh%mesh(3,1,1,k)) > particle_coords(3))) then
					get_particle_cell(3) = k
					exit
				end if		
			end do
		end if
		
		
		if (dimensions	>= 2) then
			do j = cons_inner_loop(2,1),cons_inner_loop(2,2)
				if (( 0.5_dp*(mesh%mesh(2,1,j-1,get_particle_cell(3)) + mesh%mesh(2,1,j,get_particle_cell(3))) <= particle_coords(2)).and. &
					( 0.5_dp*(mesh%mesh(2,1,j+1,get_particle_cell(3)) + mesh%mesh(2,1,j,get_particle_cell(3))) > particle_coords(2))) then
					get_particle_cell(2) = j
					exit
				end if		
			end do
		end if
		
		do i = cons_inner_loop(1,1),cons_inner_loop(1,2)

			if ((  0.5_dp*(mesh%mesh(1,i-1,get_particle_cell(2),get_particle_cell(3)) + mesh%mesh(1,i,get_particle_cell(2),get_particle_cell(3))) <= particle_coords(1)).and. &
				(  0.5_dp*(mesh%mesh(1,i+1,get_particle_cell(2),get_particle_cell(3)) + mesh%mesh(1,i,get_particle_cell(2),get_particle_cell(3))) > particle_coords(1))) then				
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
	
	subroutine apply_boundary_conditions_main(this, time)
		class(lagrangian_particles_solver)		,intent(inout)		:: this
		real(dp)					            ,intent(in)			:: time
        
		character(len=20)		:: boundary_type_name
		real(dp)				:: farfield_density, farfield_pressure, wall_temperature
		
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
		!											v_d%pr(dim1)%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3)) = 0.0_dp
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

        real(dp)	,dimension(:,:)	,allocatable	:: lengths
		integer, dimension(3)	:: cell, initial_cell
		logical	:: out_flag
		
		dimensions		= this%domain%get_domain_dimensions()
		
        
        allocate(lengths(dimensions,2))
		lengths			= this%domain%get_domain_lengths()
 
		associate(  v_f			=> this%v_f%v_ptr		, &
					particle	    => this%particles_params	, &
					chem        => this%chem%chem_ptr   , &
					mesh		=> this%mesh%mesh_ptr	, &
					bc			=> this%boundary%bc_ptr)
		
			do part = 1, this%particles_number
				if (this%particles(part)%outside_domain) then
				
					this%particles(part)%coords(1)		= particles_velocity * release_time !0.001_dp
                    this%particles(part)%coords(2)		= lengths(2,2) / 2 + streams_distance / 2

					this%particles(part)%dm				= 0.0_dp  
					
					this%particles(part)%outside_domain = .false.
				
					this%particles(part)%temperature		= 300.0_dp
					this%particles(part)%mass			= Pi*this%particles_params%diameter**3 / 6.0_dp * this%particles_params%material_density
					
					initial_cell = this%get_particle_cell(this%particles(part)%coords,out_flag)
					this%particles(part)%cell = initial_cell
					
					this%particles(part)%velocity(1)	= particles_velocity !1.0_dp!v_f%pr(1)%cells(1,initial_cell(1),initial_cell(2),initial_cell(3))
					
					this%particles(part)%outside_domain = out_flag	
					
					exit
				end if
            end do	
            
            do part = 1, this%particles_number
				if (this%particles(part)%outside_domain) then
				
					this%particles(part)%coords(1)		= particles_velocity * release_time !0.001_dp
                    this%particles(part)%coords(2)		= lengths(2,2) / 2 - streams_distance / 2

					this%particles(part)%dm				= 0.0_dp  
					
					this%particles(part)%outside_domain = .false.
				
					this%particles(part)%temperature		= 300.0_dp
					this%particles(part)%mass			= Pi*this%particles_params%diameter**3 / 6.0_dp * this%particles_params%material_density
					
					initial_cell = this%get_particle_cell(this%particles(part)%coords,out_flag)
					this%particles(part)%cell = initial_cell
					
					this%particles(part)%velocity(1)	= particles_velocity !1.0_dp!v_f%pr(1)%cells(1,initial_cell(1),initial_cell(2),initial_cell(3))
					
					this%particles(part)%outside_domain = out_flag	
					
					exit
				end if
			end do	
	
		end associate		
		
	end subroutine
	
	
	
end module
