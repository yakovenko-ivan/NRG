module droplets_solver_class

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
	public	:: droplets_solver, droplets_solver_c

	type(field_scalar_cons)	,dimension(:)	,allocatable	,target	:: T_d, T_d_int, rho_d, mass_d, numdens_d, E_f_prod_d, rho_prod_d, Esurf
	type(field_scalar_flow)	,dimension(:)	,allocatable	,target	:: m_flux_d, numdens_flux_d
	type(field_vector_cons)	,dimension(:)	,allocatable	,target	:: v_d, v_d_int, v_prod_d, Y_prod_d	

	type	:: droplets_solver
		type(field_scalar_cons_pointer)		:: T_d, T_d_int, rho_d, mass_d, numdens_d, T, rho, nu, kappa, E_f_prod, rho_prod, Esurf, foam_marker, p, time_boil, time_decay, werhop, werhop2
		type(field_scalar_flow_pointer)		:: m_flux_d, numdens_flux_d
		type(field_vector_cons_pointer)		:: v_d, v_d_int, v, v_prod, Y_prod	
		type(computational_domain)			:: domain

		type(boundary_conditions_pointer)	:: boundary
		type(computational_mesh_pointer)	:: mesh
		type(chemical_properties_pointer)	:: chem
		type(thermophysical_properties_pointer)		:: thermo
		
		type(liquid_droplets_phase)			:: droplets_params

		real(dkind)							:: droplet_mass
        
        integer                             :: phase_number
	contains
		procedure				::  set_initial_distributions
		procedure				::	droplets_euler_step_v_E
		procedure				::	droplets_lagrange_step
		procedure				::	droplets_final_step
		procedure				::  apply_boundary_conditions_main
		procedure				::  apply_boundary_conditions_interm_v_d
		procedure				::	pre_constructor
	end type

	interface	droplets_solver_c
		module procedure	constructor
	end interface

contains

	subroutine pre_constructor(this,number_of_phases)
		class(droplets_solver)	,intent(inout)	:: this	
		integer					,intent(in)		:: number_of_phases
		
		allocate(	T_d(number_of_phases), T_d_int(number_of_phases), rho_d(number_of_phases), mass_d(number_of_phases), numdens_d(number_of_phases), E_f_prod_d(number_of_phases), rho_prod_d(number_of_phases), &
					m_flux_d(number_of_phases), numdens_flux_d(number_of_phases), v_d(number_of_phases), v_d_int(number_of_phases), v_prod_d(number_of_phases), Y_prod_d(number_of_phases), Esurf(number_of_phases))	
	
	end subroutine

	type(droplets_solver)	function constructor(manager,liquid_droplets,phase_number)

		type(data_manager)			, intent(inout)	:: manager
		type(liquid_droplets_phase)	, intent(in)	:: liquid_droplets
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
		
		constructor%droplets_params	= liquid_droplets
		
		call manager%get_cons_field_pointer_by_name(scal_ptr,vect_ptr,tens_ptr,'density')
		constructor%rho%s_ptr			=> scal_ptr%s_ptr
		call manager%get_cons_field_pointer_by_name(scal_ptr,vect_ptr,tens_ptr,'temperature')
		constructor%T%s_ptr			=> scal_ptr%s_ptr
		call manager%get_cons_field_pointer_by_name(scal_ptr,vect_ptr,tens_ptr,'viscosity')
		constructor%nu%s_ptr			=> scal_ptr%s_ptr
		call manager%get_cons_field_pointer_by_name(scal_ptr,vect_ptr,tens_ptr,'thermal_conductivity')
		constructor%kappa%s_ptr			=> scal_ptr%s_ptr
		call manager%get_cons_field_pointer_by_name(scal_ptr,vect_ptr,tens_ptr,'foam_marker')
		constructor%foam_marker%s_ptr			=> scal_ptr%s_ptr
		call manager%get_cons_field_pointer_by_name(scal_ptr,vect_ptr,tens_ptr,'pressure')
		constructor%p%s_ptr			=> scal_ptr%s_ptr
		call manager%get_cons_field_pointer_by_name(scal_ptr,vect_ptr,tens_ptr,'time_boil')
		constructor%time_boil%s_ptr			=> scal_ptr%s_ptr
		!call manager%get_cons_field_pointer_by_name(scal_ptr,vect_ptr,tens_ptr,'time_decay')
		!constructor%time_decay%s_ptr			=> scal_ptr%s_ptr
		!call manager%get_cons_field_pointer_by_name(scal_ptr,vect_ptr,tens_ptr,'werhop')
		!constructor%werhop%s_ptr			=> scal_ptr%s_ptr
		!call manager%get_cons_field_pointer_by_name(scal_ptr,vect_ptr,tens_ptr,'werhop2')
		!constructor%werhop2%s_ptr			=> scal_ptr%s_ptr
		
		call manager%get_cons_field_pointer_by_name(scal_ptr,vect_ptr,tens_ptr,'velocity')
		constructor%v%v_ptr				=> vect_ptr%v_ptr		
		
		write(var_name,'(A,I2.2)')		'temperature_droplets', phase_number
		write(var_short_name,'(A,I2.2)')	'T_d', phase_number
		call manager%create_scalar_field(T_d(phase_number),	var_name,	var_short_name)
		constructor%T_d%s_ptr		=> T_d(phase_number)
		
		write(var_name,'(A,I2.2)')		'temperature_droplets_interm', phase_number
		write(var_short_name,'(A,I2.2)')	'T_d_int', phase_number
		call manager%create_scalar_field(T_d_int(phase_number),	var_name,	var_short_name)
		constructor%T_d_int%s_ptr		=> T_d_int(phase_number)
		
		write(var_name,'(A,I2.2)')		'density_droplets', phase_number
		write(var_short_name,'(A,I2.2)')	'rho_d', phase_number
		call manager%create_scalar_field(rho_d(phase_number),	var_name,	var_short_name)
		constructor%rho_d%s_ptr		=> rho_d(phase_number)

		write(var_name,'(A,I2.2)')		'mass_droplets', phase_number
		write(var_short_name,'(A,I2.2)')	'mass_d', phase_number
		call manager%create_scalar_field(mass_d(phase_number),	var_name,	var_short_name)
		constructor%mass_d%s_ptr		=> mass_d(phase_number)        
    
		write(var_name,'(A,I2.2)')		'number_density_droplets', phase_number
		write(var_short_name,'(A,I2.2)')	'numdens_d', phase_number
		call manager%create_scalar_field(numdens_d(phase_number),	var_name,	var_short_name)
		constructor%numdens_d%s_ptr		=> numdens_d(phase_number)    
        
		write(var_name,'(A,I2.2)')		'velocity_droplets', phase_number
		write(var_short_name,'(A,I2.2)')	'v_d', phase_number		
		call manager%create_vector_field(v_d(phase_number),	var_name,	var_short_name,	'spatial')
		constructor%v_d%v_ptr		=> v_d(phase_number)		
		
		write(var_name,'(A,I2.2)')		'velocity_droplets_interm', phase_number
		write(var_short_name,'(A,I2.2)')	'v_d_int', phase_number		
		call manager%create_vector_field(v_d_int(phase_number),	var_name,	var_short_name,	'spatial')
		constructor%v_d_int%v_ptr		=> v_d_int(phase_number)			
		
		write(var_name,'(A,I2.2)')		'mass_flux_droplets', phase_number
		write(var_short_name,'(A,I2.2)')	'm_flux_d', phase_number		
		call manager%create_scalar_field(m_flux_d(phase_number),		var_name,	var_short_name)
		constructor%m_flux_d%s_ptr		=> m_flux_d(phase_number)		
	
		write(var_name,'(A,I2.2)')		'number_density_flux_droplets', phase_number
		write(var_short_name,'(A,I2.2)')	'm_flux_d', phase_number		
		call manager%create_scalar_field(numdens_flux_d(phase_number),		var_name,	var_short_name)
		constructor%numdens_flux_d%s_ptr		=> numdens_flux_d(phase_number)	        
 
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

		write(var_name,'(A,I2.2)')		'Esurf', phase_number
		write(var_short_name,'(A,I2.2)')	'Esurf', phase_number
		call manager%create_scalar_field(Esurf(phase_number),	var_name,	var_short_name)
		constructor%Esurf%s_ptr		=> Esurf(phase_number) 

		constructor%phase_number    = phase_number

		constructor%mesh%mesh_ptr	=> manager%computational_mesh_pointer%mesh_ptr
		constructor%boundary%bc_ptr => manager%boundary_conditions_pointer%bc_ptr
		constructor%thermo%thermo_ptr	=> manager%thermophysics%thermo_ptr
		constructor%chem%chem_ptr	=> manager%chemistry%chem_ptr
		constructor%domain			= manager%domain
		
	end function

	subroutine set_initial_distributions(this)
		class(droplets_solver)	,intent(inout)	:: this

		integer		,dimension(3,2)	:: cons_inner_loop

		integer :: i,j,k

		cons_inner_loop	= this%domain%get_local_inner_cells_bounds()		

		associate( 	rho_d		=> this%rho_d%s_ptr			, &
                    mass_d      => this%mass_d%s_ptr        , &
                    numdens_d   => this%numdens_d%s_ptr     , &
                    Esurf       => this%Esurf%s_ptr         , &
            		T_d			=> this%T_d%s_ptr			, &
					T			=> this%T%s_ptr				, &
					droplet 	=> this%droplets_params	    , &
					mesh		=> this%mesh%mesh_ptr		, &
					bc			=> this%boundary%bc_ptr)	
		
					
		do k = cons_inner_loop(3,1),cons_inner_loop(3,2)
		do j = cons_inner_loop(2,1),cons_inner_loop(2,2)
		do i = cons_inner_loop(1,1),cons_inner_loop(1,2)	
			if(bc%bc_markers(i,j,k) == 0) then
                mass_d%cells(i,j,k)     = Pi*droplet%diameter**3 / 6.0_dkind * droplet%material_density
                numdens_d%cells(i,j,k)  = rho_d%cells(i,j,k)/mass_d%cells(i,j,k)
				T_d%cells(i,j,k)	    = 300.0_dkind 
                Esurf%cells(i,j,k)      = 0.03_dkind * 4.0 * Pi * 0.0001_dkind * 0.0001_dkind * numdens_d%cells(i,j,k) 
			end if
		end do
		end do
		end do

		end associate		
	end subroutine
	
	
	subroutine droplets_euler_step_v_E(this,time_step)

		class(droplets_solver)	,intent(inout)	:: this
		real(dkind)				,intent(in)		:: time_step
		real(dkind)								:: F_stokes, Q_stokes, Nusselt
		
		integer		,dimension(3,2)	:: cons_inner_loop
		real(dkind)	,dimension(3)	:: cell_size
		real(dkind)                 :: local_diameter, evaporation_rate, local_Esurf, source_Esurf, p_0, T_0, velocity_2, marker1, specie_enthalpy, velabs, temp_cr
        
        
		integer	:: dimensions
        integer :: H2O_index
		integer :: C7H16_index
		integer	:: sign
		integer :: i,j,k,plus,dim

		dimensions		= this%domain%get_domain_dimensions()
		
		cons_inner_loop	= this%domain%get_local_inner_cells_bounds()		

		cell_size		= this%mesh%mesh_ptr%get_cell_edges_length()

		Nusselt			= 2.0_dkind
		
		H2O_index		= this%chem%chem_ptr%get_chemical_specie_index('H2O')
		C7H16_index		= this%chem%chem_ptr%get_chemical_specie_index('C7H16')	
		
		associate(  T			=> this%T%s_ptr			, &
					T_d			=> this%T_d%s_ptr		, &	
					T_d_int		=> this%T_d_int%s_ptr	, &
					rho			=> this%rho%s_ptr		, &
					rho_d		=> this%rho_d%s_ptr		, &	
                    mass_d      => this%mass_d%s_ptr    , &
                    Esurf       => this%Esurf%s_ptr     , &
                    foam_marker => this%foam_marker%s_ptr, &
                    time_boil   => this%time_boil%s_ptr , &
					p			=> this%p%s_ptr			, &
                    numdens_d   => this%numdens_d%s_ptr , &
                    v			=> this%v%v_ptr			, &
					v_d			=> this%v_d%v_ptr		, &
					v_d_int		=> this%v_d_int%v_ptr	, &
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

		!$omp parallel default(none)  private(i,j,k,dim,F_stokes,Q_stokes,local_diameter,evaporation_rate,velabs,temp_cr) , &
		!$omp& firstprivate(this)	,&
		!$omp& shared(T,T_d,T_d_int,E_f_prod,rho,rho_d,mass_d,numdens_d,v_d,v_d_int,v,v_prod,Y_prod,nu,kappa,foam_marker,time_boil,droplet,Nusselt,H2O_index,C7H16_index,time_step,mesh,bc,cons_inner_loop,dimensions,marker1,specie_enthalpy)
		!$omp do collapse(3) schedule(guided)					
					
		do k = cons_inner_loop(3,1),cons_inner_loop(3,2)
		do j = cons_inner_loop(2,1),cons_inner_loop(2,2)
		do i = cons_inner_loop(1,1),cons_inner_loop(1,2)
			if(bc%bc_markers(i,j,k) == 0) then
                
                local_diameter = 6.0_dkind * mass_d%cells(i,j,k) / Pi / droplet%material_density 
                local_diameter = local_diameter ** 0.3333333_dkind

				F_stokes = 0.0_dkind
				
				do dim = 1,dimensions
					F_stokes			= 3.0_dkind * Pi * local_diameter * nu%cells(i,j,k) / mass_d%cells(i,j,k) * ( v%pr(dim)%cells(i,j,k) - v_d%pr(dim)%cells(i,j,k)) ! [m/s^2]
			
					v_prod%pr(dim)%cells(i,j,k)		= - F_stokes * rho_d%cells(i,j,k) / rho%cells(i,j,k)					! [m/s^2]
					v_d_int%pr(dim)%cells(i,j,k)	= v_d%pr(dim)%cells(i,j,k) + F_stokes * time_step						

					E_f_prod%cells(i,j,k)			= - F_stokes * v%pr(dim)%cells(i,j,k) * rho_d%cells(i,j,k)				! [J/m^3/s]
				end do		
					
				Q_stokes					= 3.0_dkind * kappa%cells(i,j,k) * Nusselt / (0.5_dkind * local_diameter ** 2 * droplet%material_heat_capacity * droplet%material_density) * (T%cells(i,j,k) - T_d%cells(i,j,k))

				E_f_prod%cells(i,j,k)		= E_f_prod%cells(i,j,k) - Q_stokes * rho_d%cells(i,j,k) * droplet%material_heat_capacity ! [J/m^3/s]
				T_d_int%cells(i,j,k)		= T_d%cells(i,j,k) + Q_stokes * time_step

                if (droplet%combustible == .false.) then
					temp_cr = 373.15
                else
					temp_cr = 371.55    
                end if
				
				if (T_d_int%cells(i,j,k) > temp_cr) T_d_int%cells(i,j,k) = temp_cr
                
                if (T%cells(i,j,k) >= temp_cr) then
                    evaporation_rate = 2.0_dkind * local_diameter * Pi * kappa%cells(i,j,k)/droplet%material_heat_capacity  &
                                     * log(1.0_dkind+droplet%material_heat_capacity*(T%cells(i,j,k) - T_d%cells(i,j,k))/droplet%material_latent_heat)/mass_d%cells(i,j,k)	! [1/s]
                    
                    if (j > cons_inner_loop(2,2) - 2)  evaporation_rate = 0.0_dkind
                    
					if (mass_d%cells(i,j,k)*(1.0_dkind - evaporation_rate*time_step) >= 0.0_dkind) then
						mass_d%cells(i,j,k) =	mass_d%cells(i,j,k)*(1.0_dkind - evaporation_rate*time_step)
                        
						if(mass_d%cells(i,j,k) <= 1.0E-15_dkind) then
							mass_d%cells(i,j,k) = 1.0E-15_dkind
							evaporation_rate = 0.0_dkind
						end if
                        
						rho_d%cells(i,j,k)  =   rho_d%cells(i,j,k)*(1.0_dkind - evaporation_rate*time_step)

						if(rho_d%cells(i,j,k) <= 1.0E-16_dkind) then
							rho_d%cells(i,j,k) = 1.0E-16_dkind
							evaporation_rate = 0.0_dkind
						end if

!						rho%cells(i,j,k)    =   rho%cells(i,j,k) + evaporation_rate*rho_d%cells(i,j,k)*time_step
                        
						if (droplet%combustible == .false.) then
							Y_prod%pr(H2O_index)%cells(i,j,k)	= evaporation_rate*rho_d%cells(i,j,k)	! [kg/m^3/s]
						else
							Y_prod%pr(C7H16_index)%cells(i,j,k)	= evaporation_rate*rho_d%cells(i,j,k)
						end if

						do dim = 1,dimensions
							v_prod%pr(dim)%cells(i,j,k)     = v_prod%pr(dim)%cells(i,j,k)	+ evaporation_rate*rho_d%cells(i,j,k)/rho%cells(i,j,k)*v_d%pr(dim)%cells(i,j,k)	! [m/s^2]
						end do

						E_f_prod%cells(i,j,k)			    = E_f_prod%cells(i,j,k) - evaporation_rate*droplet%material_latent_heat*rho_d%cells(i,j,k)		! [J/m^3/s]
                        
					else
						evaporation_rate    = 0.0_dkind
					end if 
                end if
			end if
		end do
		end do
		end do

		!$omp end do nowait
		!$omp end parallel		
		
		end associate
	end subroutine

	
	subroutine droplets_lagrange_step(this,time_step)
 
		class(droplets_solver)	,intent(inout)	:: this
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
		
		associate(  m_flux_d	=> this%m_flux_d%s_ptr	, &
					numdens_flux_d	=> this%numdens_flux_d%s_ptr	, &
                    mass_d      => this%mass_d%s_ptr    , &
					numdens_d	=> this%numdens_d%s_ptr	, &		
					rho_d		=> this%rho_d%s_ptr		, &
					v_d_int		=> this%v_d_int%v_ptr	, &
					mesh		=> this%mesh%mesh_ptr	, &
					bc			=> this%boundary%bc_ptr)
					
		!$omp parallel default(none)  private(i,j,k,dim,plus,sign,bound_number,boundary_type_name,av_velocity,dif_velocity,cell_surface_area) , &
		!$omp& firstprivate(this)	,&
		!$omp& shared(bc,mesh,m_flux_d,numdens_flux_d,mass_d,rho_d,v_d_int,dimensions,cell_size,time_step,flow_inner_loop,cons_inner_loop,coordinate_system)
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
							if(dim==2) cell_surface_area(dim) = cell_surface_area(dim) * (mesh%mesh(1,i,j,k))		
						case ('spherical')
							! x -> r
							if(dim==1) cell_surface_area(dim) = cell_surface_area(dim) * (mesh%mesh(1,i,j,k) - 0.5_dkind*cell_size(1))**2	
					end select				
		
                    av_velocity     = 0.5_dkind *(v_d_int%pr(dim)%cells(i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3)) + v_d_int%pr(dim)%cells(i,j,k))
                    dif_velocity    = time_step *(v_d_int%pr(dim)%cells(i,j,k) - v_d_int%pr(dim)%cells(i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3)))/ cell_size(dim)
                    if( av_velocity >= 0 ) then
                        m_flux_d%cells(dim,i,j,k)       = rho_d%cells(i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3))	* av_velocity * (cell_surface_area(dim) ) * time_step / (1.0_dkind + dif_velocity)
                        if (mass_d%cells(i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3)) /= 0.0_dkind) then 
                            numdens_flux_d%cells(dim,i,j,k) = m_flux_d%cells(dim,i,j,k) / mass_d%cells(i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3))
                        else
                            numdens_flux_d%cells(dim,i,j,k) = 0.0_dkind
                        end if
                    else
                        m_flux_d%cells(dim,i,j,k)       = rho_d%cells(i,j,k)									* av_velocity * (cell_surface_area(dim) ) * time_step / (1.0_dkind + dif_velocity)
                        if (mass_d%cells(i,j,k) /= 0.0_dkind) then
                            numdens_flux_d%cells(dim,i,j,k) = m_flux_d%cells(dim,i,j,k) / mass_d%cells(i,j,k)
                        else
                            numdens_flux_d%cells(dim,i,j,k) = 0.0_dkind
                        end if
                    end if
                    if((bc%bc_markers(i,j,k) == 0).and.(i <= cons_inner_loop(1,2)).and.(j <= cons_inner_loop(2,2)-2).and.(k <= cons_inner_loop(3,2))) then
						do plus = 1,2
							sign			= (-1)**plus
							bound_number	= bc%bc_markers(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3))
							if( bound_number /= 0 ) then
								boundary_type_name = bc%boundary_types(bound_number)%get_type_name()
								select case(boundary_type_name)
									case ('wall')		! Particles stay near top wall
										if ((dim == 2).and.(sign == 1)) then
											m_flux_d%cells(dim,i,j,k) = max(0.0_dkind,m_flux_d%cells(dim,i,j,k))							
										end if
								end select
							end if					
						end do
					end if
                    
					continue
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

	subroutine droplets_final_step(this,time_step)
 
		class(droplets_solver)	,intent(inout)	:: this
		real(dkind)				,intent(in)		:: time_step
 
		real(dkind)	:: D11, D12, D21, D22, rho_d_old, av_velocity1, av_velocity2
		
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
		
		associate(  m_flux_d	=> this%m_flux_d%s_ptr	, &
                    numdens_flux_d	=> this%numdens_flux_d%s_ptr	, &
                    mass_d      => this%mass_d%s_ptr    , &
					numdens_d	=> this%numdens_d%s_ptr	, &
					rho_d		=> this%rho_d%s_ptr		, &
					T_d			=> this%T_d%s_ptr		, &
					T_d_int		=> this%T_d_int%s_ptr	, &
					v_d			=> this%v_d%v_ptr		, &
					v_d_int		=> this%v_d_int%v_ptr	, &
					mesh		=> this%mesh%mesh_ptr	, &
					bc			=> this%boundary%bc_ptr)
 
		!$omp parallel default(none)  private(i,j,k,dim,dim1,dim2,spec,rho_d_old,av_velocity1,av_velocity2,cell_volume,D11,D21,D12,D22,i_ind1,i_ind2,j_ind1,j_ind2,k_ind1,k_ind2) , &
		!$omp& firstprivate(this)	,&
		!$omp& shared(bc,mesh,m_flux_d,rho_d,T_d_int,T_d,v_d,v_d_int,numdens_d,numdens_flux_d,mass_d,dimensions,cons_inner_loop,coordinate_system)
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
 
				rho_d_old			    = rho_d%cells(i,j,k)
				do dim = 1,dimensions
					rho_d%cells(i,j,k)  = rho_d%cells(i,j,k) + (m_flux_d%cells(dim,i,j,k)-m_flux_d%cells(dim,i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3))) / (cell_volume)
                end do

				do dim = 1,dimensions
					numdens_d%cells(i,j,k)  = numdens_d%cells(i,j,k) + (numdens_flux_d%cells(dim,i,j,k)-numdens_flux_d%cells(dim,i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3))) / (cell_volume)
                end do
                
                mass_d%cells(i,j,k) = rho_d%cells(i,j,k) / numdens_d%cells(i,j,k)					
                
				T_d%cells(i,j,k)	= T_d_int%cells(i,j,k)  * rho_d_old / rho_d%cells(i,j,k)

				do dim = 1,dimensions
					v_d%pr(dim)%cells(i,j,k)	= v_d_int%pr(dim)%cells(i,j,k) * rho_d_old / rho_d%cells(i,j,k)
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
 
                    av_velocity1 = v_d_int%pr(dim2)%cells(i_ind1,j_ind1,k_ind1) + v_d_int%pr(dim2)%cells(i,j,k)
                    av_velocity2 = v_d_int%pr(dim2)%cells(i_ind2,j_ind2,k_ind2) + v_d_int%pr(dim2)%cells(i,j,k)
 
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
 
 
					T_d%cells(i,j,k) = T_d%cells(i,j,k)		+ (D11 * T_d_int%cells(i_ind1,j_ind1,k_ind1) * abs(m_flux_d%cells(dim2,i,j,k)) &
															+  D21 * T_d_int%cells(i_ind2,j_ind2,k_ind2) * abs(m_flux_d%cells(dim2,i_ind2,j_ind2,k_ind2)) &
															-  D12 * T_d_int%cells(i,j,k)				 * abs(m_flux_d%cells(dim2,i,j,k)) &
															-  D22 * T_d_int%cells(i,j,k)				 * abs(m_flux_d%cells(dim2,i_ind2,j_ind2,k_ind2))) /  rho_d%cells(i,j,k)  /(cell_volume)
 
					do dim1 = 1,dimensions
						v_d%pr(dim1)%cells(i,j,k) = v_d%pr(dim1)%cells(i,j,k)	+ (D11 * v_d_int%pr(dim1)%cells(i_ind1,j_ind1,k_ind1)	* abs(m_flux_d%cells(dim2,i,j,k)) &
																				+  D21 * v_d_int%pr(dim1)%cells(i_ind2,j_ind2,k_ind2)	* abs(m_flux_d%cells(dim2,i_ind2,j_ind2,k_ind2)) &
																				-  D12 * v_d_int%pr(dim1)%cells(i,j,k)				* abs(m_flux_d%cells(dim2,i,j,k)) &
																				-  D22 * v_d_int%pr(dim1)%cells(i,j,k)				* abs(m_flux_d%cells(dim2,i_ind2,j_ind2,k_ind2))) / rho_d%cells(i,j,k) / (cell_volume)
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

	subroutine apply_boundary_conditions_interm_v_d(this)

		class(droplets_solver)		,intent(inout)		:: this

		integer	:: dimensions
		integer	,dimension(3,2)	:: cons_inner_loop
		
		character(len=20)		:: boundary_type_name		

		integer	:: sign, bound_number
		integer :: i,j,k,plus,dim,dim1,specie_number
								
		dimensions		= this%domain%get_domain_dimensions()
				
		cons_inner_loop	= this%domain%get_local_inner_cells_bounds()			
		
		associate(  v_d_int			=> this%v_d_int%v_ptr		, &
					bc				=> this%boundary%bc_ptr	, &
					mesh			=> this%mesh%mesh_ptr)

		!$omp parallel default(none)  private(i,j,k,plus,dim,dim1,sign,bound_number,boundary_type_name) , &
		!$omp& firstprivate(this)	,&
		!$omp& shared(v_d_int,bc,cons_inner_loop,dimensions)
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
										v_d_int%pr(dim1)%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3)) = - v_d_int%pr(dim1)%cells(i,j,k)
									else
										v_d_int%pr(dim1)%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3)) = v_d_int%pr(dim1)%cells(i,j,k)
									end if
								end do
								
								boundary_type_name = bc%boundary_types(bound_number)%get_type_name()
								select case(boundary_type_name)
									case ('inlet')
										do dim1 = 1, dimensions
											if(dim1 == dim)	v_d_int%pr(dim1)%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3)) = v_d_int%pr(dim1)%cells(i,j,k)
										end do
									case ('outlet')
										do dim1 = 1, dimensions
											if(dim1 == dim)	v_d_int%pr(dim1)%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3)) = v_d_int%pr(dim1)%cells(i,j,k)
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
	
	subroutine apply_boundary_conditions_main(this, time)

		class(droplets_solver)		,intent(inout)		:: this
		real(dkind)					,intent(in)			:: time

		character(len=20)		:: boundary_type_name
		real(dkind)				:: farfield_density, farfield_pressure, farfield_rhod, wall_temperature
		real(dkind)				:: delay
		
		
		integer					:: dimensions
		integer	,dimension(3,2)	:: cons_inner_loop

		integer	:: sign, bound_number
		integer :: i,j,k,plus,dim,dim1
						
		dimensions		= this%domain%get_domain_dimensions()
			
		cons_inner_loop	= this%domain%get_local_inner_cells_bounds()			
		
		delay = 0.2_dkind
		farfield_rhod = 1.0e-02_dkind

		associate(  T_d				=> this%T_d%s_ptr			, &
					rho_d			=> this%rho_d%s_ptr			, &
					mass_d          => this%mass_d%s_ptr        , &
					v_d				=> this%v_d%v_ptr			, &
					droplet			=> this%droplets_params	    , &
					bc				=> this%boundary%bc_ptr		, &
					mesh			=> this%mesh%mesh_ptr)

		!$omp parallel default(none)  private(i,j,k,plus,dim,dim1,sign,bound_number,boundary_type_name,wall_temperature) , &
		!$omp& firstprivate(this)	,&
		!$omp& shared(T_d,rho_d,v_d,mass_d,droplet,mesh,bc,cons_inner_loop,dimensions,time, farfield_rhod, delay)
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

								rho_d%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3))	= rho_d%cells(i,j,k)
								T_d%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3))	= T_d%cells(i,j,k)
                                mass_d%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3)) = mass_d%cells(i,j,k)

								do dim1 = 1, dimensions
									if(dim1 == dim) then
										v_d%pr(dim1)%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3)) = - v_d%pr(dim1)%cells(i,j,k)
									else
										v_d%pr(dim1)%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3)) = v_d%pr(dim1)%cells(i,j,k)
									end if
								end do								
								
								boundary_type_name = bc%boundary_types(bound_number)%get_type_name()
								select case(boundary_type_name)
									case('wall')
										if(bc%boundary_types(bound_number)%is_conductive()) then 
											wall_temperature = bc%boundary_types(bound_number)%get_wall_temperature()
											T_d%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3)) = wall_temperature
										end if
										if(.not.bc%boundary_types(bound_number)%is_slip()) then
											do dim1 = 1, dimensions
												if (dim1 /= dim) then	
													v_d%pr(dim1)%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3)) = 0.0_dkind
												end if
											end do
										end if
									case('inlet')
									
										if (time > 0.2) then
											rho_d%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3)) = farfield_rhod * (time - 0.2) / delay + 1.0e-05_dkind
											if (time > 0.2 + delay) then
												rho_d%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3)) = farfield_rhod
											end if
										else
											rho_d%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3)) = 1.0e-05_dkind
										end if
										
										T_d%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3))		= 300.0_dkind
										mass_d%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3))     = Pi*droplet%diameter**3 / 6.0_dkind * droplet%material_density
										do dim1 = 1, dimensions
											if(dim1 == dim) then
												v_d%pr(dim1)%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3)) = v_d%pr(dim1)%cells(i,j,k)
											else
												v_d%pr(dim1)%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3)) = v_d%pr(dim1)%cells(i,j,k)
											end if
										end do									
									case ('outlet')
										do dim1 = 1, dimensions
											if(dim1 == dim)	v_d%pr(dim1)%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3)) = v_d%pr(dim1)%cells(i,j,k)
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
	
	
	
end module
