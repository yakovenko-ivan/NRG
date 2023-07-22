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
	type(field_vector_flow)	,target	:: v_f

	type	:: coarse_particles
		type(field_scalar_cons_pointer)		:: p, E_f_prod, rho, E_f, e_i, E_f_int, p_int, foam_marker
		type(field_scalar_flow_pointer)		:: m_flux
		type(field_vector_cons_pointer)		:: v_prod, v, v_int, Y_int, Y
		type(field_scalar_cons_pointer)		:: p_dyn, p_stat, div_v		
		type(field_vector_flow_pointer)		:: v_f
		
		type(computational_domain)			:: domain
		type(mpi_communications)			:: mpi_support		
		type(boundary_conditions_pointer)	:: boundary
		type(computational_mesh_pointer)	:: mesh
		type(chemical_properties_pointer)	:: chem

		real(dkind) ,dimension(:,:,:)	,allocatable    :: div_v_old
		
	contains
		procedure				::	euler_step_v
		procedure				::	euler_step_E
		procedure				::	lagrange_step
		procedure				:: calculate_p_dyn		
		procedure				::	final_step
	end type

	interface	coarse_particles_c
		module procedure	constructor
	end interface

contains

	type(coarse_particles)	function constructor(manager,low_mach_flag)

		type(data_manager)	, intent(inout)				:: manager
		logical				, intent(in)	,optional	:: low_mach_flag

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
		constructor%v_prod%v_ptr		=> v_prod_gd
		
		

		constructor%mesh%mesh_ptr	=> manager%computational_mesh_pointer%mesh_ptr
		constructor%boundary%bc_ptr => manager%boundary_conditions_pointer%bc_ptr
		constructor%chem%chem_ptr	=> manager%chemistry%chem_ptr
		constructor%domain			= manager%domain
		constructor%mpi_support		= manager%mpi_communications
	
		if(present(low_mach_flag)) then
			allocate(constructor%div_v_old(	cons_allocation_bounds(1,1):cons_allocation_bounds(1,2), &
											cons_allocation_bounds(2,1):cons_allocation_bounds(2,2), &
											cons_allocation_bounds(3,1):cons_allocation_bounds(3,2)))
        
			constructor%div_v_old			= 0.0_dkind  
 		end if
										
	end function

	subroutine euler_step_v(this,time_step)

		class(coarse_particles)	,intent(inout)	:: this
		real(dkind)				,intent(in)		:: time_step
		real(dkind)	:: av_pressure1, av_velocity1, av_pressure2, av_velocity2, div_p_v	

		integer		,dimension(3,2)	:: cons_inner_loop
		real(dkind)	,dimension(3)	:: cell_size
		
		integer	:: dimensions
		integer	:: sign
		integer :: i,j,k,plus,dim

		dimensions		= this%domain%get_domain_dimensions()
		
		cons_inner_loop	= this%domain%get_local_inner_cells_bounds()		

		cell_size		= this%mesh%mesh_ptr%get_cell_edges_length()

		associate(  p			=> this%p%s_ptr			, &
					v_prod		=> this%v_prod%v_ptr	, &
					v			=> this%v%v_ptr			, &
					rho			=> this%rho%s_ptr		, &
					mesh		=> this%mesh%mesh_ptr	, &
					bc			=> this%boundary%bc_ptr)

		call this%mpi_support%exchange_conservative_scalar_field(p)

	!$omp parallel default(none)  private(i,j,k,dim) , &
	!$omp& firstprivate(this)	,&
	!$omp& shared(bc,p,v_prod,rho,dimensions,cell_size,time_step,cons_inner_loop)
	!$omp do collapse(3) schedule(guided)

		do k = cons_inner_loop(3,1),cons_inner_loop(3,2)
		do j = cons_inner_loop(2,1),cons_inner_loop(2,2)
		do i = cons_inner_loop(1,1),cons_inner_loop(1,2)
			if(bc%bc_markers(i,j,k) == 0) then
				do dim = 1,dimensions
					v_prod%pr(dim)%cells(i,j,k)  = - 0.5_dkind*(p%cells(i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3)) - p%cells(i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3))) * time_step / cell_size(dim) / rho%cells(i,j,k) 
					v_prod%pr(dim)%cells(i,j,k) = v_prod%pr(dim)%cells(i,j,k)  + g(dim) * (rho%cells(1,1,1) - rho%cells(i,j,k)) * time_step / rho%cells(i,j,k) 
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
		real(dkind)				,intent(in)		:: time_step
		real(dkind)	:: av_pressure1, av_velocity1, av_pressure2, av_velocity2, div_p_v	
		
		real(dkind), dimension (3,3)			:: lame_coeffs
		
		integer		,dimension(3,2)	:: cons_inner_loop
		
		character(len=20)	:: coordinate_system
		
		real(dkind)	,dimension(3)	:: cell_size
		integer	:: dimensions
		integer	:: sign
		integer :: i,j,k,plus,dim

		dimensions		= this%domain%get_domain_dimensions()
				
		cons_inner_loop	= this%domain%get_local_inner_cells_bounds()

		cell_size		= this%mesh%mesh_ptr%get_cell_edges_length()
		
		coordinate_system	= this%domain%get_coordinate_system_name()
		
		associate(		p			=> this%p%s_ptr			, &
						v_int		=> this%v_int%v_ptr		, &
						E_f_prod	=> this%E_f_prod%s_ptr	, &
						rho			=> this%rho%s_ptr		, &
						mesh		=> this%mesh%mesh_ptr	, &
						bc			=> this%boundary%bc_ptr)
						
		call this%mpi_support%exchange_conservative_vector_field(v_int)
			
	!$omp parallel default(none)  private(i,j,k,dim,div_p_v,av_pressure1,av_velocity1,av_pressure2,av_velocity2,lame_coeffs) , &
	!$omp& firstprivate(this)	,&
	!$omp& shared(bc,mesh,p,v_int,rho,E_f_prod,dimensions,cell_size,time_step,cons_inner_loop,coordinate_system)
	!$omp do collapse(3) schedule(guided)

		do k = cons_inner_loop(3,1),cons_inner_loop(3,2)
		do j = cons_inner_loop(2,1),cons_inner_loop(2,2)
		do i = cons_inner_loop(1,1),cons_inner_loop(1,2)
			if(bc%bc_markers(i,j,k) == 0) then
				div_p_v = 0.0_dkind

				lame_coeffs		= 1.0_dkind				
				
				select case(coordinate_system)
					case ('cartesian')	
						lame_coeffs			= 1.0_dkind
					case ('cylindrical')
						! x -> z, y -> r
						lame_coeffs(1,1)	=  mesh%mesh(1,i,j,k) - 0.5_dkind*cell_size(1)			
						lame_coeffs(1,2)	=  mesh%mesh(1,i,j,k)
						lame_coeffs(1,3)	=  mesh%mesh(1,i,j,k) + 0.5_dkind*cell_size(1)	
					case ('spherical')
						! x -> r
						lame_coeffs(1,1)	=  (mesh%mesh(1,i,j,k) - 0.5_dkind*cell_size(1))**2
						lame_coeffs(1,2)	=  (mesh%mesh(1,i,j,k))**2
						lame_coeffs(1,3)	=  (mesh%mesh(1,i,j,k) + 0.5_dkind*cell_size(1))**2
				end select				
				
				do dim = 1,dimensions
		
                    av_velocity1	= 0.5_dkind * (v_int%pr(dim)%cells(i,j,k)	+ v_int%pr(dim)	%cells(i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3)))
                    av_pressure1	= 0.5_dkind * (p			%cells(i,j,k)	+ p				%cells(i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3)))
					
                    av_velocity2	= 0.5_dkind * (v_int%pr(dim)%cells(i,j,k)	+ v_int%pr(dim)	%cells(i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3)))
                    av_pressure2	= 0.5_dkind * (p			%cells(i,j,k)	+ p				%cells(i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3)))
					
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
		real(dkind)				,intent(in)		:: time_step

		real(dkind)	:: av_velocity, dif_velocity
		
		integer		,dimension(3,2)	:: flow_inner_loop
		
		character(len=20)	:: coordinate_system

		real(dkind)	,dimension(3)	:: cell_size	, cell_surface_area	
		integer	:: dimensions
		integer	:: sign
		integer :: i,j,k,plus,dim

		dimensions		= this%domain%get_domain_dimensions()
		
		cell_size			= this%mesh%mesh_ptr%get_cell_edges_length()
		coordinate_system	= this%domain%get_coordinate_system_name()
		
		flow_inner_loop	= this%domain%get_local_inner_faces_bounds()
		
		associate(  m_flux		=> this%m_flux%s_ptr	, &
        !            foam_marker => this%foam_marker%s_ptr, &
					rho			=> this%rho%s_ptr		, &
					v_int		=> this%v_int%v_ptr		, &
					mesh		=> this%mesh%mesh_ptr	, &
					bc			=> this%boundary%bc_ptr)

		call this%mpi_support%exchange_conservative_vector_field(v_int)
		call this%mpi_support%exchange_conservative_scalar_field(rho)

	!$omp parallel default(none)  private(i,j,k,dim,av_velocity,dif_velocity,cell_surface_area) , &
	!$omp& firstprivate(this)	,&
	!$omp& shared(bc,mesh,m_flux,rho,v_int,dimensions,cell_size,time_step,flow_inner_loop,coordinate_system)
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
	
                    av_velocity     = 0.5_dkind *(v_int%pr(dim)%cells(i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3)) + v_int%pr(dim)%cells(i,j,k))
                    dif_velocity    = time_step *(v_int%pr(dim)%cells(i,j,k) - v_int%pr(dim)%cells(i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3))) / cell_size(dim)
                    if( av_velocity >= 0 ) then
                        m_flux%cells(dim,i,j,k) = rho%cells(i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3))		* av_velocity * (cell_surface_area(dim) ) * time_step / (1.0_dkind + dif_velocity)
                    else
                        m_flux%cells(dim,i,j,k) = rho%cells(i,j,k)										* av_velocity * (cell_surface_area(dim) ) * time_step / (1.0_dkind + dif_velocity)
                    end if
                    
                  !  if(abs(foam_marker%cells(i,j,k) - foam_marker%cells(i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3))) == 1.0_dkind) then
                  !      m_flux%cells(dim,i,j,k) = 0.0_dkind
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
		real(dkind)				,intent(in)		:: time_step

		real(dkind)	:: D11, D12, D21, D22, rho_old, av_velocity1, av_velocity2
		
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

	call this%mpi_support%exchange_conservative_scalar_field(E_f_int)
	call this%mpi_support%exchange_conservative_vector_field(Y_int)

	call this%mpi_support%exchange_flow_scalar_field(m_flux)

	!$omp parallel default(none)  private(i,j,k,dim,dim1,dim2,spec,rho_old,av_velocity1,av_velocity2,cell_volume,D11,D21,D12,D22,i_ind1,i_ind2,j_ind1,j_ind2,k_ind1,k_ind2) , &
	!$omp& firstprivate(this)	,&
	!$omp& shared(bc,mesh,m_flux,rho,E_f_int,E_f,v,v_int,v_f,Y,Y_int,dimensions,species_number,cons_inner_loop,coordinate_system)
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

                    av_velocity1 = v_int%pr(dim2)%cells(i_ind1,j_ind1,k_ind1) + v_int%pr(dim2)%cells(i,j,k)
                    av_velocity2 = v_int%pr(dim2)%cells(i_ind2,j_ind2,k_ind2) + v_int%pr(dim2)%cells(i,j,k)

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
					v_f%pr(dim)%cells(dim,i,j,k) = 0.5_dkind * ( v%pr(dim)%cells(i,j,k) + v%pr(dim)%cells(i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3)))
				end do
			end if
		end do
		end do
		end do		
		
	!$omp end parallel

		end associate
	end subroutine

	subroutine calculate_p_dyn(this,time_step)

#ifdef mpi
	use MPI
#endif

		class(coarse_particles)	,intent(inout)	:: this
		real(dkind)				,intent(in)		:: time_step
		
        real(dkind) :: rho_min(1)
        real(dkind)	,dimension(:)	,allocatable	,save	:: rho_min_array
		
		real(dkind)	:: div_pres_flux, pres_flux1, pres_flux2
		
		real(dkind), dimension (3,3)	:: lame_coeffs	
		
		integer	                    :: dimensions, iterations
		integer	,dimension(3,2)	    :: cons_inner_loop, cons_utter_loop
		real(dkind)	,dimension(3)	:: cell_size        
        real(dkind)					:: time_step_adj
		
		integer						:: processor_rank, processor_number, mpi_communicator		
		character(len=20)			:: coordinate_system
		
		logical	:: converged
		
		integer	:: iteration
		integer	:: bound_number
		integer :: i,j,k,dim,specie_number
		
		integer	:: error

		dimensions		= this%domain%get_domain_dimensions()
		
		cons_inner_loop	= this%domain%get_local_inner_cells_bounds()
		cons_utter_loop = this%domain%get_local_utter_cells_bounds()
			
        cell_size		= this%mesh%mesh_ptr%get_cell_edges_length()

		time_step_adj	=  0.125_dkind*cell_size(1)**2	
		
		coordinate_system	= this%domain%get_coordinate_system_name()		
		
		processor_rank		= this%domain%get_processor_rank()
		mpi_communicator	= this%domain%get_mpi_communicator()

		if (.not.allocated(rho_min_array)) then
			processor_number = this%domain%get_mpi_communicator_size()
			allocate(rho_min_array(processor_number))
		end if		
		
        rho_min(1) = 100.0_dkind
        
		associate (	rho				=> this%rho%s_ptr			, &
					p				=> this%p%s_ptr				, &
                    p_dyn			=> this%p_dyn%s_ptr			, &
					p_stat			=> this%p_stat%s_ptr		, &
                    p_int			=> this%p_int%s_ptr			, &
					div_v			=> this%div_v%s_ptr         , &
					div_v_old       => this%div_v_old			, &
                    e_i             => this%e_i%s_ptr			, &
                    E_f             => this%E_f%s_ptr			, &
					mesh			=> this%mesh%mesh_ptr		, &
					bc				=> this%boundary%bc_ptr)
			
		!!$omp parallel default(none)  private(i,j,k) , &
		!!$omp& firstprivate(this)	,&
		!!$omp& shared(time_step,div_v_old,div_v,p_int,rho,e_i,E_f,rho_min,p,p_dyn,p_stat,bc,cons_inner_loop)
		!!$omp do collapse(3) schedule(guided) reduction(min:rho_min)
					
			do k = cons_inner_loop(3,1),cons_inner_loop(3,2)
			do j = cons_inner_loop(2,1),cons_inner_loop(2,2)
			do i = cons_inner_loop(1,1),cons_inner_loop(1,2)
				if(bc%bc_markers(i,j,k) == 0) then					
                    div_v_old(i,j,k) = div_v%cells(i,j,k)
                    div_v%cells(i,j,k) = rho%cells(i,j,k)/p_stat%cells(i,j,k)/time_step * (e_i%cells(i,j,k) - E_f%cells(i,j,k))
					p_int%cells(i,j,k) = p_dyn%cells(i,j,k)
                    if (rho%cells(i,j,k) < rho_min(1))  then
						rho_min(1) = rho%cells(i,j,k)
					end if
				end if
			end do
			end do
            end do

		!	p_int%cells = p_dyn%cells
			
		!!$omp end do nowait
		!!$omp end parallel			
		
		rho_min_array(processor_rank+1) = rho_min(1) 
		
#ifdef mpi
		call mpi_gather(rho_min,1,MPI_DOUBLE_PRECISION,rho_min_array,1,MPI_DOUBLE_PRECISION,0,mpi_communicator,error)
#endif
		
		if (processor_rank == 0) then
			do i = 0,size(rho_min_array) - 1 
				if (rho_min_array(i+1) < rho_min(1)) rho_min(1) = rho_min_array(i+1)
			end do
		end if

#ifdef mpi	
		call mpi_bcast(rho_min,1,MPI_DOUBLE_PRECISION,0,mpi_communicator,error)
#endif
	
		iteration = 0
		converged = .false.
		
		!$omp parallel default(none)  private(i,j,k,dim,pres_flux1,pres_flux2,div_pres_flux,lame_coeffs,iterations) , &
		!$omp& firstprivate(this) , &
		!$omp& shared(time_step,time_step_adj,coordinate_system,p_dyn,p_int,div_v,div_v_old,rho_min,bc,cell_size,dimensions,cons_inner_loop,mesh,converged,iteration)
		do while (.not.converged)!do iterations = 1, 100
		!do iterations = 1, 10
		!	call this%mpi_support%exchange_conservative_scalar_field(p_int)	
		
			!$omp barrier		
			!$omp atomic write
			converged = .true.
			!$omp end atomic			
				
			!$omp do collapse(3) schedule(guided)
			do k = cons_inner_loop(3,1),cons_inner_loop(3,2)
			do j = cons_inner_loop(2,1),cons_inner_loop(2,2)
			do i = cons_inner_loop(1,1),cons_inner_loop(1,2)
				if(bc%bc_markers(i,j,k) == 0) then
								
					lame_coeffs		= 1.0_dkind		
				 
					select case(coordinate_system)
						case ('cartesian')	
							lame_coeffs			= 1.0_dkind
						case ('cylindrical')
							lame_coeffs(1,1)	=  mesh%mesh(1,i,j,k) - 0.5_dkind*cell_size(1)
							lame_coeffs(1,2)	=  mesh%mesh(1,i,j,k)
							lame_coeffs(1,3)	=  mesh%mesh(1,i,j,k) + 0.5_dkind*cell_size(1)	
						case ('spherical')
							lame_coeffs(1,1)	=  (mesh%mesh(1,i,j,k) - 0.5_dkind*cell_size(1))**2
							lame_coeffs(1,2)	=  (mesh%mesh(1,i,j,k))**2
							lame_coeffs(1,3)	=  (mesh%mesh(1,i,j,k) + 0.5_dkind*cell_size(1))**2
					end select					
				 
					p_dyn%cells(i,j,k) = p_int%cells(i,j,k) - time_step_adj * (div_v%cells(i,j,k) - div_v_old(i,j,k))/time_step * rho_min(1) 
					
					div_pres_flux = 0.0_dkind
					do dim = 1,dimensions
						pres_flux1	= lame_coeffs(dim,1)* (p_int%cells(i,j,k) - p_int%cells(i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3))) / cell_size(dim)
     
						pres_flux2	= lame_coeffs(dim,3)* (p_int%cells(i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3)) - p_int%cells(i,j,k)) / cell_size(dim)
     
						div_pres_flux = div_pres_flux  + (pres_flux2 - pres_flux1) / cell_size(dim) / lame_coeffs(dim,2)
					end do			
					
					p_dyn%cells(i,j,k) = p_dyn%cells(i,j,k) + time_step_adj * div_pres_flux
				end if
			end do
			end do
			end do            
			!$omp end do
			
			!$omp do collapse(3) schedule(guided) reduction(.and.:converged)	
			do k = cons_inner_loop(3,1),cons_inner_loop(3,2)
			do j = cons_inner_loop(2,1),cons_inner_loop(2,2)
			do i = cons_inner_loop(1,1),cons_inner_loop(1,2)
			
				if ((abs(p_int%cells(i,j,1) - p_dyn%cells(i,j,1)) > 1.0e-02_dkind)) then
					converged = .false.
				!	print *, abs(p_int%cells(i,j,1) - p_dyn%cells(i,j,1))
				end if			
   
				p_int%cells(i,j,k) = p_dyn%cells(i,j,k)
			end do
			end do
			end do 
			!$omp end do		
			
			!$omp master
			iteration = iteration + 1
			!$omp end master
		
		end do
		!$omp end parallel
		
		print *, 'Poisson iteration:', iteration

		do k = cons_inner_loop(3,1),cons_inner_loop(3,2)
		do j = cons_inner_loop(2,1),cons_inner_loop(2,2)
		do i = cons_inner_loop(1,1),cons_inner_loop(1,2)		
			p%cells(i,j,k) = p_dyn%cells(i,j,k) + p_stat%cells(i,j,k)
		end do
		end do
		end do 	
		
		end associate

	end subroutine	
	
end module
