module coarse_particles_method

	use kind_parameters
	use global_data
	use field_pointers
	use boundary_conditions_class
	use data_manager_class
	use computational_mesh_class
	use computational_domain_class
	use thermophysical_properties_class
	use chemical_properties_class

	use mpi_communications_class

	implicit none

#ifdef OMP	
	include "omp_lib.h"
#endif

	private
	public	:: coarse_particles, coarse_particles_c

	type(field_scalar_cons)	,target	:: p_dyn, p_stat, p_int	
	type(field_scalar_cons)	,target	:: div_v
	type(field_scalar_cons)	,target	:: E_f_prod_gd
	type(field_scalar_flow)	,target	:: m_flux
	type(field_vector_cons)	,target	:: v_prod_gd	
    type(field_vector_cons)	,target	:: v_prod_sources
	type(field_vector_flow)	,target	:: v_f

	type	:: coarse_particles
		type(field_scalar_cons_pointer)		:: p, E_f_prod, rho, E_f, e_i, E_f_int, p_int, foam_marker
		type(field_scalar_flow_pointer)		:: m_flux
		type(field_vector_cons_pointer)		:: v_prod_gd, v_prod_sources, v, v_int, Y_int, Y
		type(field_scalar_cons_pointer)		:: p_dyn, p_stat, div_v		
		type(field_vector_flow_pointer)		:: v_f
		
		type(computational_domain)			:: domain
		type(mpi_communications)			:: mpi_support		
		type(boundary_conditions_pointer)	:: boundary
		type(computational_mesh_pointer)	:: mesh
		type(chemical_properties_pointer)	:: chem

		real(dp),   dimension(:,:,:)	,allocatable    :: div_v_old
		real(dp),   dimension(3)                        :: g
	contains
		procedure				::	euler_step_v
		procedure				::	euler_step_E
		procedure				::	lagrange_step
		procedure				::	final_step
        procedure				::	perturb_velocity_field
	end type

	interface	coarse_particles_c
		module procedure	constructor
	end interface

contains

	type(coarse_particles)	function constructor(manager,g,low_mach_flag)

		type(data_manager)	    , intent(inout)				:: manager
        real(dp), dimension(3)  , intent(in)                :: g
		logical				    , intent(in)	,optional	:: low_mach_flag

		type(field_scalar_cons_pointer)	:: scal_ptr
		type(field_vector_cons_pointer)	:: vect_ptr
		type(field_tensor_cons_pointer)	:: tens_ptr
		
 		integer					:: dimensions       
 		integer	,dimension(3,2)	:: cons_allocation_bounds, flow_allocation_bounds   
        
		cons_allocation_bounds		= manager%domain%get_local_utter_cells_bounds()
		flow_allocation_bounds		= manager%domain%get_local_utter_faces_bounds()
		dimensions					= manager%domain%get_domain_dimensions()   		
		
		call manager%get_cons_field_pointer_by_name(scal_ptr,vect_ptr,tens_ptr,'pressure')
		constructor%p%s_ptr				=> scal_ptr%s_ptr
		call manager%get_cons_field_pointer_by_name(scal_ptr,vect_ptr,tens_ptr,'density')
		constructor%rho%s_ptr			=> scal_ptr%s_ptr
		call manager%get_cons_field_pointer_by_name(scal_ptr,vect_ptr,tens_ptr,'full_energy')
		constructor%E_f%s_ptr			=> scal_ptr%s_ptr
		call manager%get_cons_field_pointer_by_name(scal_ptr,vect_ptr,tens_ptr,'internal_energy')
		constructor%e_i%s_ptr			=> scal_ptr%s_ptr 		
		call manager%get_cons_field_pointer_by_name(scal_ptr,vect_ptr,tens_ptr,'full_energy_interm')
		constructor%E_f_int%s_ptr		=> scal_ptr%s_ptr

	!	call manager%get_cons_field_pointer_by_name(scal_ptr,vect_ptr,tens_ptr,'foam_marker')
	!	constructor%foam_marker%s_ptr	=> scal_ptr%s_ptr        
        
		call manager%create_scalar_field(p_dyn		,'pressure_dynamic'						,'p_dyn')
		constructor%p_dyn%s_ptr						=> p_dyn	
		call manager%create_scalar_field(p_stat		,'pressure_static'						,'p_stat')
		constructor%p_stat%s_ptr					=> p_stat		
		call manager%create_scalar_field(div_v		,'velocity_divergence'					,'div_v')
		constructor%div_v%s_ptr						=> div_v
		
		call manager%create_scalar_field(p_int		,'pressure_interm'						,'p_int')
		constructor%p_int%s_ptr						=> p_int		
		
		call manager%create_scalar_field(E_f_prod_gd,'energy_production_gas_dynamics'	,'E_f_prod_gd')
		constructor%E_f_prod%s_ptr		=> E_f_prod_gd
		call manager%create_scalar_field(m_flux	,'mass_flux'							,'m_flux')
		constructor%m_flux%s_ptr		=> m_flux
		call manager%create_vector_field(v_f	,'velocity_flow'						,'v_f'		,'spatial')
		constructor%v_f%v_ptr			=> v_f
		
		call manager%get_cons_field_pointer_by_name(scal_ptr,vect_ptr,tens_ptr,'velocity')
		constructor%v%v_ptr				=> vect_ptr%v_ptr
		call manager%get_cons_field_pointer_by_name(scal_ptr,vect_ptr,tens_ptr,'velocity_interm')
		constructor%v_int%v_ptr			=> vect_ptr%v_ptr
		call manager%get_cons_field_pointer_by_name(scal_ptr,vect_ptr,tens_ptr,'specie_molar_concentration')
		constructor%Y%v_ptr				=> vect_ptr%v_ptr
		call manager%get_cons_field_pointer_by_name(scal_ptr,vect_ptr,tens_ptr,'specie_molar_concentration_interm')
		constructor%Y_int%v_ptr			=> vect_ptr%v_ptr		
		
		call manager%create_vector_field(v_prod_gd	,'velocity_production_gas_dynamics'		,'v_prod_gd'	,'spatial')
		constructor%v_prod_gd%v_ptr		=> v_prod_gd
		
		call manager%create_vector_field(v_prod_sources	,'velocity_production_sources'		,'v_prod_sources'	,'spatial')
		constructor%v_prod_sources%v_ptr	=> v_prod_sources
				

		constructor%mesh%mesh_ptr	=> manager%computational_mesh_pointer%mesh_ptr
		constructor%boundary%bc_ptr => manager%boundary_conditions_pointer%bc_ptr
		constructor%chem%chem_ptr	=> manager%chemistry%chem_ptr
		constructor%domain			= manager%domain
		constructor%mpi_support		= manager%mpi_communications
	
		if(present(low_mach_flag)) then
			allocate(constructor%div_v_old(	cons_allocation_bounds(1,1):cons_allocation_bounds(1,2), &
											cons_allocation_bounds(2,1):cons_allocation_bounds(2,2), &
											cons_allocation_bounds(3,1):cons_allocation_bounds(3,2)))
        
			constructor%div_v_old			= 0.0_dp  
 		end if
										
	end function

	subroutine euler_step_v(this,time_step)

		class(coarse_particles)	,intent(inout)	:: this
		real(dp)				,intent(in)		:: time_step
		real(dp)	:: av_pressure1, av_velocity1, av_pressure2, av_velocity2, div_p_v	

		integer		,dimension(3,2)	:: cons_inner_loop
		real(dp)	,dimension(3)	:: cell_size
		
		integer	:: dimensions
		integer	:: sign
		integer :: i,j,k,plus,dim

		dimensions		= this%domain%get_domain_dimensions()
		
		cons_inner_loop	= this%domain%get_local_inner_cells_bounds()		

		cell_size		= this%mesh%mesh_ptr%get_cell_edges_length()

		
        associate(  p			=> this%p%s_ptr)
		call this%mpi_support%exchange_conservative_scalar_field(p)
        end associate
 
		associate(  p			=> this%p%s_ptr			, &
					v_prod_gd	=> this%v_prod_gd%v_ptr	, &
					v			=> this%v%v_ptr			, &
					rho			=> this%rho%s_ptr		, &
					mesh		=> this%mesh%mesh_ptr	, &
					bc			=> this%boundary%bc_ptr)

	!$omp parallel default(shared)  private(i,j,k,dim) !, &
	!!$omp& shared(this,dimensions,cell_size,time_step,cons_inner_loop)
 
	!$omp do collapse(3) schedule(guided)

		do k = cons_inner_loop(3,1),cons_inner_loop(3,2)
		do j = cons_inner_loop(2,1),cons_inner_loop(2,2)
		do i = cons_inner_loop(1,1),cons_inner_loop(1,2)
			if(bc%bc_markers(i,j,k) == 0) then
				do dim = 1,dimensions
					v_prod_gd%pr(dim)%cells(i,j,k)  = - 0.5_dp*(p%cells(i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3)) - p%cells(i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3))) * time_step / cell_size(dim) / rho%cells(i,j,k) 
					v_prod_gd%pr(dim)%cells(i,j,k) = v_prod_gd%pr(dim)%cells(i,j,k)  + this%g(dim) * (rho%cells(1,1,1) - rho%cells(i,j,k)) * time_step / rho%cells(i,j,k) 
				end do
			end if
		end do
		end do
		end do

	!$omp end do nowait
	!$omp end parallel

		end associate

	end subroutine

	subroutine euler_step_E(this,time_step)
		class(coarse_particles)	,intent(inout)	:: this
		real(dp)				,intent(in)		:: time_step
		real(dp)	:: av_pressure1, av_velocity1, av_pressure2, av_velocity2, div_p_v	
		
		real(dp), dimension (3,3)			:: lame_coeffs
		
		integer		,dimension(3,2)	:: cons_inner_loop
		
		character(len=20)	:: coordinate_system
		
		real(dp)	,dimension(3)	:: cell_size
		integer	:: dimensions
		integer	:: sign
		integer :: i,j,k,plus,dim

		dimensions		= this%domain%get_domain_dimensions()
				
		cons_inner_loop	= this%domain%get_local_inner_cells_bounds()

		cell_size		= this%mesh%mesh_ptr%get_cell_edges_length()
		
		coordinate_system	= this%domain%get_coordinate_system_name()
		
        associate(		v_int		=> this%v_int%v_ptr)
		    call this%mpi_support%exchange_conservative_vector_field(v_int)
        end associate
 
		associate(		p			=> this%p%s_ptr			, &
						v_int		=> this%v_int%v_ptr		, &
						E_f_prod	=> this%E_f_prod%s_ptr	, &
						rho			=> this%rho%s_ptr		, &
						mesh		=> this%mesh%mesh_ptr	, &
						bc			=> this%boundary%bc_ptr)
						
	!$omp parallel default(shared)  private(i,j,k,dim,div_p_v,av_pressure1,av_velocity1,av_pressure2,av_velocity2,lame_coeffs) !, &
	!!$omp& shared(this,dimensions,cell_size,time_step,cons_inner_loop,coordinate_system)
	
	!$omp do collapse(3) schedule(guided)
		do k = cons_inner_loop(3,1),cons_inner_loop(3,2)
		do j = cons_inner_loop(2,1),cons_inner_loop(2,2)
		do i = cons_inner_loop(1,1),cons_inner_loop(1,2)
			if(bc%bc_markers(i,j,k) == 0) then
				div_p_v = 0.0_dp

				lame_coeffs		= 1.0_dp				
				
				select case(coordinate_system)
					case ('cartesian')	
						lame_coeffs			= 1.0_dp
					case ('cylindrical')
						! x -> z, y -> r
						lame_coeffs(1,1)	=  mesh%mesh(1,i,j,k) - 0.5_dp*cell_size(1)			
						lame_coeffs(1,2)	=  mesh%mesh(1,i,j,k)
						lame_coeffs(1,3)	=  mesh%mesh(1,i,j,k) + 0.5_dp*cell_size(1)	
					case ('spherical')
						! x -> r
						lame_coeffs(1,1)	=  (mesh%mesh(1,i,j,k) - 0.5_dp*cell_size(1))**2
						lame_coeffs(1,2)	=  (mesh%mesh(1,i,j,k))**2
						lame_coeffs(1,3)	=  (mesh%mesh(1,i,j,k) + 0.5_dp*cell_size(1))**2
				end select				
				
				do dim = 1,dimensions
		
                    av_velocity1	= 0.5_dp * (v_int%pr(dim)%cells(i,j,k)	+ v_int%pr(dim)	%cells(i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3)))
                    av_pressure1	= 0.5_dp * (p			%cells(i,j,k)	+ p				%cells(i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3)))
					
                    av_velocity2	= 0.5_dp * (v_int%pr(dim)%cells(i,j,k)	+ v_int%pr(dim)	%cells(i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3)))
                    av_pressure2	= 0.5_dp * (p			%cells(i,j,k)	+ p				%cells(i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3)))
					
					div_p_v = div_p_v + (av_velocity2*av_pressure2*lame_coeffs(dim,3) - av_velocity1*av_pressure1*lame_coeffs(dim,1)) * time_step / lame_coeffs(dim,2) / cell_size(dim) / rho%cells(i,j,k) 
				end do

				E_f_prod%cells(i,j,k) = - (div_p_v) 
			end if
		end do
		end do
		end do
	
	!$omp end do nowait
	!$omp end parallel


		end associate
	end subroutine
	
	
	subroutine lagrange_step(this,time_step)

		class(coarse_particles)	,intent(inout)	:: this
		real(dp)				,intent(in)		:: time_step

		real(dp)	:: av_velocity, dif_velocity
		
		integer		,dimension(3,2)	:: flow_inner_loop
		
		character(len=20)	:: coordinate_system

		real(dp)	,dimension(3)	:: cell_size	, cell_surface_area	
		integer	:: dimensions
		integer	:: sign
		integer :: i,j,k,plus,dim

		dimensions		= this%domain%get_domain_dimensions()
		
		cell_size			= this%mesh%mesh_ptr%get_cell_edges_length()
		coordinate_system	= this%domain%get_coordinate_system_name()
		
		flow_inner_loop	= this%domain%get_local_inner_faces_bounds()
		
        associate(	rho			=> this%rho%s_ptr		, &
					v_int		=> this%v_int%v_ptr)
            
		call this%mpi_support%exchange_conservative_vector_field(v_int)
		call this%mpi_support%exchange_conservative_scalar_field(rho)

    end associate
    
		associate(  m_flux		=> this%m_flux%s_ptr	, &
        !            foam_marker => this%foam_marker%s_ptr, &
					rho			=> this%rho%s_ptr		, &
					v_int		=> this%v_int%v_ptr		, &
					mesh		=> this%mesh%mesh_ptr	, &
					bc			=> this%boundary%bc_ptr)
	!$omp parallel default(shared)  private(i,j,k,dim,av_velocity,dif_velocity,cell_surface_area) !, &
	!!$omp& shared(this,dimensions,cell_size,time_step,flow_inner_loop,coordinate_system)

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
							if(dim==1) cell_surface_area(dim) = cell_surface_area(dim) * (mesh%mesh(1,i,j,k) - 0.5_dp*cell_size(1))									
							if(dim==2) cell_surface_area(dim) = cell_surface_area(dim) * (mesh%mesh(1,i,j,k))		
						case ('spherical')
							! x -> r
							if(dim==1) cell_surface_area(dim) = cell_surface_area(dim) * (mesh%mesh(1,i,j,k) - 0.5_dp*cell_size(1))**2	
					end select		
	
                    av_velocity     = 0.5_dp *(v_int%pr(dim)%cells(i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3)) + v_int%pr(dim)%cells(i,j,k))
                    dif_velocity    = time_step *(v_int%pr(dim)%cells(i,j,k) - v_int%pr(dim)%cells(i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3))) / cell_size(dim)
                    if( av_velocity >= 0 ) then
                        m_flux%cells(dim,i,j,k) = rho%cells(i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3))		* av_velocity * (cell_surface_area(dim) ) * time_step / (1.0_dp + dif_velocity)
                    else
                        m_flux%cells(dim,i,j,k) = rho%cells(i,j,k)										* av_velocity * (cell_surface_area(dim) ) * time_step / (1.0_dp + dif_velocity)
                    end if
                    
                  !  if(abs(foam_marker%cells(i,j,k) - foam_marker%cells(i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3))) == 1.0_dp) then
                  !      m_flux%cells(dim,i,j,k) = 0.0_dp
                  !  end if    

                end do
	!		end if
		end do
		end do
		end do

	!$omp end do nowait
	!$omp end parallel

		end associate

	end subroutine

	subroutine final_step(this,time_step)

		class(coarse_particles)	,intent(inout)	:: this
		real(dp)				,intent(in)		:: time_step

		real(dp)	:: D11, D12, D21, D22, rho_old, av_velocity1, av_velocity2
		
		integer		,dimension(3,2)	:: cons_inner_loop
		real(dp)	,dimension(3)	:: cell_size
		real(dp)					:: cell_volume
		integer	:: dimensions, species_number
		character(len=20)	:: coordinate_system
		integer	:: sign
		integer :: i,j,k,plus,dim,dim1,dim2,spec,i_ind1,i_ind2,j_ind1,j_ind2,k_ind1,k_ind2

		dimensions			= this%domain%get_domain_dimensions()
		species_number		= this%chem%chem_ptr%species_number
		coordinate_system	= this%domain%get_coordinate_system_name()
		
		cons_inner_loop	= this%domain%get_local_inner_cells_bounds()
	
		cell_size			= this%mesh%mesh_ptr%get_cell_edges_length()
		

    associate(  m_flux		=> this%m_flux%s_ptr	, &
			    E_f_int		=> this%E_f_int%s_ptr	, &
			    Y_int		=> this%Y_int%v_ptr)            

	call this%mpi_support%exchange_conservative_scalar_field(E_f_int)
	call this%mpi_support%exchange_conservative_vector_field(Y_int)

	call this%mpi_support%exchange_flow_scalar_field(m_flux)
    
    end associate

                
		associate(  m_flux		=> this%m_flux%s_ptr	, &
					rho			=> this%rho%s_ptr		, &
					E_f			=> this%E_f%s_ptr		, &
					E_f_int		=> this%E_f_int%s_ptr	, &
					v			=> this%v%v_ptr			, &
					v_int		=> this%v_int%v_ptr		, &
					v_f			=> this%v_f%v_ptr		, &
					Y			=> this%Y%v_ptr			, &
					Y_int		=> this%Y_int%v_ptr		, &
					mesh		=> this%mesh%mesh_ptr	, &
					bc			=> this%boundary%bc_ptr)

	!$omp parallel default(shared)  private(i,j,k,dim,dim1,dim2,spec,rho_old,av_velocity1,av_velocity2,cell_volume,D11,D21,D12,D22,i_ind1,i_ind2,j_ind1,j_ind2,k_ind1,k_ind2) !, &
	!!$omp& shared(this,dimensions,species_number,cons_inner_loop,coordinate_system)
    
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
				
				rho_old			= rho%cells(i,j,k)
				do dim = 1,dimensions
					rho%cells(i,j,k) = rho%cells(i,j,k) + (m_flux%cells(dim,i,j,k)-m_flux%cells(dim,i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3))) / (cell_volume)
				end do

				E_f%cells(i,j,k)		= E_f_int%cells(i,j,k)  * rho_old / rho%cells(i,j,k)
				do dim = 1,dimensions
					v%pr(dim)%cells(i,j,k)	= v_int%pr(dim)%cells(i,j,k) * rho_old / rho%cells(i,j,k)
				end do
				do spec = 1,species_number
					Y%pr(spec)%cells(i,j,k) = Y_int%pr(spec)%cells(i,j,k) * rho_old / rho%cells(i,j,k)
				end do

                do dim2 = 1,dimensions
                    D11 = 0.0_dp
					D12 = 0.0_dp
					D21 = 0.0_dp
					D22 = 0.0_dp

					i_ind1 = i - I_m(dim2,1)
					i_ind2 = i + I_m(dim2,1)
					j_ind1 = j - I_m(dim2,2)
					j_ind2 = j + I_m(dim2,2)
					k_ind1 = k - I_m(dim2,3)
					k_ind2 = k + I_m(dim2,3)

                    av_velocity1 = v_int%pr(dim2)%cells(i_ind1,j_ind1,k_ind1) + v_int%pr(dim2)%cells(i,j,k)
                    av_velocity2 = v_int%pr(dim2)%cells(i_ind2,j_ind2,k_ind2) + v_int%pr(dim2)%cells(i,j,k)

                    if(av_velocity1 > 0.0) then
						D11 = 1.0_dp
					else
						D12 = 1.0_dp
					end if
                    if(av_velocity2 < 0.0) then
						D21 = 1.0_dp
					else
						D22 = 1.0_dp
					end if


					E_f%cells(i,j,k) = E_f%cells(i,j,k)		+ (D11 * E_f_int%cells(i_ind1,j_ind1,k_ind1) * abs(m_flux%cells(dim2,i,j,k)) &
															+  D21 * E_f_int%cells(i_ind2,j_ind2,k_ind2) * abs(m_flux%cells(dim2,i_ind2,j_ind2,k_ind2)) &
															-  D12 * E_f_int%cells(i,j,k)				 * abs(m_flux%cells(dim2,i,j,k)) &
															-  D22 * E_f_int%cells(i,j,k)				 * abs(m_flux%cells(dim2,i_ind2,j_ind2,k_ind2))) /  rho%cells(i,j,k)  /(cell_volume)

					do dim1 = 1,dimensions
						v%pr(dim1)%cells(i,j,k) = v%pr(dim1)%cells(i,j,k)	+ (D11 * v_int%pr(dim1)%cells(i_ind1,j_ind1,k_ind1)	* abs(m_flux%cells(dim2,i,j,k)) &
																			+  D21 * v_int%pr(dim1)%cells(i_ind2,j_ind2,k_ind2)	* abs(m_flux%cells(dim2,i_ind2,j_ind2,k_ind2)) &
																			-  D12 * v_int%pr(dim1)%cells(i,j,k)				* abs(m_flux%cells(dim2,i,j,k)) &
																			-  D22 * v_int%pr(dim1)%cells(i,j,k)				* abs(m_flux%cells(dim2,i_ind2,j_ind2,k_ind2))) / rho%cells(i,j,k) / (cell_volume)
					end do

					do spec = 1,species_number
						Y%pr(spec)%cells(i,j,k) = Y%pr(spec)%cells(i,j,k)		+ (D11 * Y_int%pr(spec)%cells(i_ind1,j_ind1,k_ind1) * abs(m_flux%cells(dim2,i,j,k)) &
																				+  D21 * Y_int%pr(spec)%cells(i_ind2,j_ind2,k_ind2) * abs(m_flux%cells(dim2,i_ind2,j_ind2,k_ind2)) &
																				-  D12 * Y_int%pr(spec)%cells(i,j,k)				 * abs(m_flux%cells(dim2,i,j,k)) &
																				-  D22 * Y_int%pr(spec)%cells(i,j,k)				 * abs(m_flux%cells(dim2,i_ind2,j_ind2,k_ind2))) / rho%cells(i,j,k) / (cell_volume) 
					end do
	
                end do
			end if
		end do
		end do
		end do
	!$omp end do 
		
	!$omp do collapse(3) schedule(guided)	
		do k = cons_inner_loop(3,1),cons_inner_loop(3,2)
		do j = cons_inner_loop(2,1),cons_inner_loop(2,2)
		do i = cons_inner_loop(1,1),cons_inner_loop(1,2)
			if(bc%bc_markers(i,j,k) == 0) then
				do dim = 1,dimensions
					v_f%pr(dim)%cells(dim,i,j,k) = 0.5_dp * ( v%pr(dim)%cells(i,j,k) + v%pr(dim)%cells(i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3)))
				end do
			end if
		end do
		end do
		end do		
    !$omp end do 		
	!$omp end parallel

		end associate


	end subroutine


    subroutine perturb_velocity_field(this,time_step)
		class(coarse_particles)	,intent(inout)	:: this
		real(dp)				,intent(in)		:: time_step

		real(dp)	,dimension(3)	:: cell_size		
		real(dp) ,save			:: time = 0.0_dp
        
		integer	:: dimensions, iterations
		integer	,dimension(3,2)	:: cons_inner_loop, cons_utter_loop, flow_inner_loop
		integer	,dimension(3,2)	:: loop
		
		character(len=20)		:: boundary_type_name
		
		integer	:: sign, bound_number
		integer :: i,j,k,dim,dim1,dim2,spec,plus

		real(dp)	:: kappa_turb, gamma_turb, kx_turb, ky_turb, kappa_turb_max, time_scale
		
		dimensions		= this%domain%get_domain_dimensions()
		
		cons_inner_loop	= this%domain%get_local_inner_cells_bounds()
		cons_utter_loop = this%domain%get_local_utter_cells_bounds()
		
		flow_inner_loop	= this%domain%get_local_inner_faces_bounds()
		
		cell_size		= this%mesh%mesh_ptr%get_cell_edges_length()
		
		associate (	v_prod_source	=> this%v_prod_sources%v_ptr	, &
					bc				=> this%boundary%bc_ptr)

            time_scale		= 1.0e-03_dp
            kappa_turb_max	= 100000.0_dp
            
            kappa_turb	= min(kappa_turb_max, kappa_turb_max/time_scale * time)

			call RANDOM_SEED()
			call RANDOM_NUMBER(gamma_turb)
			
			gamma_turb = 2.0_dp * pi * gamma_turb   !2.0_dp*gamma_turb - 1.0_dp   !cos(alpha) sin(alpha)=sqrt(1-gamma_turb**2.0)
			
			do k = cons_inner_loop(3,1),cons_inner_loop(3,2)
			do j = cons_inner_loop(2,1),cons_inner_loop(2,2)
			do i = cons_inner_loop(1,1),cons_inner_loop(1,2)
                do dim = 1, dimensions
					if((bc%bc_markers(i,j,k) == 0).and.(bc%bc_markers(i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3)) == 0)) then	

						!call RANDOM_NUMBER(gamma_turb)
						!gamma_turb = 2.0_dp*gamma_turb - 1.0_dp
						!v_f%pr(dim)%cells(dim,i,j,k) = v_f%pr(dim)%cells(dim,i,j,k) + kappa_turb*gamma_turb*sqrt(this%time_step)
							
						kx_turb = 2.0_dp*pi*sin(gamma_turb)/0.005_dp   !0.002   !sqrt(1.0_dp-gamma_turb**2.0)/0.002  
						ky_turb = 2.0_dp*pi*cos(gamma_turb)/0.005_dp   !0.002   !gamma_turb/0.002
							
						if(dim == 1)then
							v_prod_sources%pr(dim)%cells(i,j,k) = kappa_turb*cos(gamma_turb)*sqrt(time_step)*cos(kx_turb*(i-0.5)*cell_size(1)+ky_turb*(j-0.5)*cell_size(1))
						endif
						if(dim == 2)then
							v_prod_sources%pr(dim)%cells(i,j,k) = -kappa_turb*sin(gamma_turb)*sqrt(time_step)*cos(kx_turb*(i-0.5)*cell_size(1)+ky_turb*(j-0.5)*cell_size(1))
						endif                        
					end if
				end do
			end do
			end do
            end do

            time = time + time_step
            
			continue
		end associate
	end subroutine
	
end module
