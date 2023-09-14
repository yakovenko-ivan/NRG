module viscosity_solver_class

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
	public	:: viscosity_solver, viscosity_solver_c
	
	type(field_scalar_cons)	,target	:: E_f_prod_visc, E_f_prod_visc_FDS, nu
	type(field_vector_cons)	,target	:: v_prod_visc	
	type(field_tensor)		,target	:: sigma
	
	
	
	type	:: viscosity_solver
		type(field_scalar_cons_pointer)				:: E_f_prod, E_f_prod_FDS, T, nu, rho, mol_mix_conc
		type(field_vector_cons_pointer)				:: v, v_prod, Y
		type(field_tensor_cons_pointer)				:: sigma
		type(computational_domain)					:: domain
		type(mpi_communications)					:: mpi_support
		type(boundary_conditions_pointer)			:: boundary
		type(computational_mesh_pointer)			:: mesh
		type(thermophysical_properties_pointer)		:: thermo
		type(chemical_properties_pointer)			:: chem
		
		real(dkind) ,dimension(:)   ,allocatable    :: viscosity_coeff_constant
		
		real(dkind) ,dimension(:,:,:)	,allocatable    :: sigma_theta_theta		
	contains
		procedure	,private	::	calculate_viscosity_coeff_constant
		procedure	,private	::	calculate_viscosity_coeff
		procedure	,private	::	calculate_sigma
		procedure	,private	::	apply_boundary_conditions
		procedure				::	solve_viscosity
	end type
	
	interface	viscosity_solver_c
		module procedure	constructor
	end interface
	
contains

	type(viscosity_solver)	function constructor(manager)
	
		type(data_manager)	, intent(inout)	:: manager
	
		type(field_scalar_cons_pointer)	:: scal_ptr
		type(field_vector_cons_pointer)	:: vect_ptr
		type(field_tensor_cons_pointer)	:: tens_ptr	
	
		character(len=20)	:: coordinate_system
		integer	,dimension(3,2)	:: cons_allocation_bounds, flow_allocation_bounds 
		
		integer	:: dimensions
		integer	:: dim
		
		call manager%get_cons_field_pointer_by_name(scal_ptr,vect_ptr,tens_ptr,'temperature')
		constructor%T%s_ptr				=> scal_ptr%s_ptr	
		call manager%get_cons_field_pointer_by_name(scal_ptr,vect_ptr,tens_ptr,'mixture_molar_concentration')
		constructor%mol_mix_conc%s_ptr		=> scal_ptr%s_ptr
		call manager%get_cons_field_pointer_by_name(scal_ptr,vect_ptr,tens_ptr,'density')
		constructor%rho%s_ptr					=> scal_ptr%s_ptr
		
		call manager%create_scalar_field(E_f_prod_visc		,'energy_production_viscosity'		,'E_f_prod_visc')
		constructor%E_f_prod%s_ptr		=> E_f_prod_visc	

		call manager%create_scalar_field(E_f_prod_visc_FDS	,'energy_production_viscosity_FDS'	,'E_f_prod_visc_FDS')
		constructor%E_f_prod_FDS%s_ptr	=> E_f_prod_visc_FDS	        
        
		call manager%create_scalar_field(nu				,'viscosity'					,'nu')
		constructor%nu%s_ptr			=> nu
		
		call manager%get_cons_field_pointer_by_name(scal_ptr,vect_ptr,tens_ptr,'velocity')
		constructor%v%v_ptr				=> vect_ptr%v_ptr
		call manager%get_cons_field_pointer_by_name(scal_ptr,vect_ptr,tens_ptr,'specie_molar_concentration')
		constructor%Y%v_ptr				=> vect_ptr%v_ptr
		
		call manager%create_vector_field(v_prod_visc	,'velocity_production_viscosity','v_prod_visc'	,'spatial')
		constructor%v_prod%v_ptr		=> v_prod_visc

		call manager%create_tensor_field(sigma			,'stress'						,'sigma')
		constructor%sigma%t_ptr			=> sigma

		constructor%mesh%mesh_ptr	=> manager%computational_mesh_pointer%mesh_ptr
		constructor%boundary%bc_ptr => manager%boundary_conditions_pointer%bc_ptr		
		
		constructor%domain				= manager%domain
		constructor%mpi_support			= manager%mpi_communications
		
		constructor%thermo%thermo_ptr	=> manager%thermophysics%thermo_ptr
		constructor%chem%chem_ptr		=> manager%chemistry%chem_ptr
		
		dimensions = constructor%domain%get_domain_dimensions()
		
		constructor%E_f_prod%s_ptr%cells(:,:,:)		= 0.0_dkind
        constructor%E_f_prod_FDS%s_ptr%cells(:,:,:)	= 0.0_dkind
        
		do dim = 1, dimensions
			constructor%v_prod%v_ptr%pr(dim)%cells(:,:,:) = 0.0_dkind
		end do
		
		cons_allocation_bounds		= manager%domain%get_local_utter_cells_bounds()
		coordinate_system	= constructor%domain%get_coordinate_system_name()
		select case(coordinate_system)
			case ('cylindrical','spherical')
				allocate(constructor%sigma_theta_theta(	cons_allocation_bounds(1,1):cons_allocation_bounds(1,2), &
														cons_allocation_bounds(2,1):cons_allocation_bounds(2,2), &
														cons_allocation_bounds(3,1):cons_allocation_bounds(3,2)))
				constructor%sigma_theta_theta			= 0.0_dkind
		end select			
		
		allocate(constructor%viscosity_coeff_constant(manager%chemistry%chem_ptr%species_number))
		
		call constructor%calculate_viscosity_coeff_constant()
	end function

	subroutine	solve_viscosity(this,time_step)
	
		class(viscosity_solver) ,intent(inout) :: this
		real(dkind)				,intent(in)		:: time_step
		
		real(dkind)	:: div_sigma, v_div_sigma, sigma_dv
        
		character(len=20)				:: coordinate_system
		real(dkind), dimension (3,3)	:: lame_coeffs
		
		integer						:: dimensions
		integer		,dimension(3,2)	:: cons_inner_loop
		real(dkind)	,dimension(3)	:: cell_size			
		
		integer	:: sign
		integer :: i,j,k,plus,dim1,dim2

		call this%calculate_sigma()
		call this%apply_boundary_conditions()

		dimensions		= this%domain%get_domain_dimensions()

		cons_inner_loop = this%domain%get_local_inner_cells_bounds()

		cell_size		= this%mesh%mesh_ptr%get_cell_edges_length()

		coordinate_system	= this%domain%get_coordinate_system_name()
		
		associate(  v_prod			=> this%v_prod%v_ptr		, &
					v				=> this%v%v_ptr				, &
					E_f_prod		=> this%E_f_prod%s_ptr		, &
					E_f_prod_FDS	=> this%E_f_prod_FDS%s_ptr	, &
					rho				=> this%rho%s_ptr			, &
					sigma			=> this%sigma%t_ptr			, &
					mesh			=> this%mesh%mesh_ptr		, &
					bc				=> this%boundary%bc_ptr)

		call this%mpi_support%exchange_conservative_tensor_field(sigma)
				
	!$omp parallel default(none)  private(i,j,k,dim1,dim2,plus,sign,v_div_sigma,div_sigma,sigma_dv,lame_coeffs) , &
	!$omp& firstprivate(this)	,&
	!$omp& shared(bc,mesh,v_prod,v,rho,E_f_prod,E_f_prod_FDS,sigma,time_step,cons_inner_loop,dimensions,cell_size,coordinate_system) 
	!$omp do collapse(3) schedule(static)
		do k = cons_inner_loop(3,1),cons_inner_loop(3,2)
		do j = cons_inner_loop(2,1),cons_inner_loop(2,2)
		do i = cons_inner_loop(1,1),cons_inner_loop(1,2)
		
			E_f_prod%cells(i,j,k)	= 0.0_dkind
            E_f_prod_FDS%cells(i,j,k)	= 0.0_dkind
		
			if(bc%bc_markers(i,j,k) == 0) then
                
				v_div_sigma		= 0.0_dkind
				sigma_dv		= 0.0_dkind
                
				lame_coeffs		= 1.0_dkind				
			
				select case(coordinate_system)
					case ('cartesian')	
						lame_coeffs			= 1.0_dkind
					case ('cylindrical')
						! x -> z, y -> r
						lame_coeffs(1,1)	=  mesh%mesh(1,i,j,k) - cell_size(1)			
						lame_coeffs(1,2)	=  mesh%mesh(1,i,j,k)
						lame_coeffs(1,3)	=  mesh%mesh(1,i,j,k) + cell_size(1)	
					
					case ('spherical')
						! x -> r
						lame_coeffs(1,1)	=  (mesh%mesh(1,i,j,k) - cell_size(1))**2
						lame_coeffs(1,2)	=  (mesh%mesh(1,i,j,k))**2
						lame_coeffs(1,3)	=  (mesh%mesh(1,i,j,k) + cell_size(1))**2
				end select							
				
				
				do dim1 = 1,dimensions
					div_sigma	= 0.0_dkind
					
					do dim2 = 1,dimensions
							
						div_sigma	= div_sigma + (	sigma%pr(dim1,dim2)%cells(i+I_m(dim2,1),j+I_m(dim2,2),k+I_m(dim2,3)) &
												-	sigma%pr(dim1,dim2)%cells(i-I_m(dim2,1),j-I_m(dim2,2),k-I_m(dim2,3)))/(2.0_dkind * cell_size(dim2))
                        
                        sigma_dv	= sigma_dv + sigma%pr(dim1,dim2)%cells(i,j,k) * (v%pr(dim1)%cells(i+I_m(dim2,1),j+I_m(dim2,2),k+I_m(dim2,3)) - v%pr(dim1)%cells(i-I_m(dim2,1),j-I_m(dim2,2),k-I_m(dim2,3))) / (2.0_dkind*cell_size(dim2))
                    end do
                    
					select case(coordinate_system)
					case ('cylindrical')
						if(dim1==1) div_sigma = div_sigma + (sigma%pr(1,1)%cells(i,j,k) - this%sigma_theta_theta(i,j,k))/mesh%mesh(1,i,j,k)
						if(dim1==2) div_sigma = div_sigma + (sigma%pr(1,2)%cells(i,j,k))/mesh%mesh(1,i,j,k)
					case ('spherical')
						if(dim1==1) div_sigma = div_sigma + 2.0_dkind*(sigma%pr(1,1)%cells(i,j,k) - this%sigma_theta_theta(i,j,k))/mesh%mesh(1,i,j,k)
                    end select			
				
                    v_div_sigma	= v_div_sigma + div_sigma * v%pr(dim1)%cells(i,j,k)    
					
					v_prod%pr(dim1)%cells(i,j,k)	=  div_sigma !* time_step

				end do
				
				E_f_prod%cells(i,j,k) = v_div_sigma + sigma_dv
                
                E_f_prod_FDS%cells(i,j,k) = sigma_dv

			end if
		end do
		end do
		end do
	!$omp end do nowait
	!$omp end parallel

		call this%mpi_support%exchange_conservative_scalar_field(E_f_prod)
        call this%mpi_support%exchange_conservative_scalar_field(E_f_prod_FDS)
		call this%mpi_support%exchange_conservative_vector_field(v_prod)		
		
		end associate
	end subroutine

	
	subroutine calculate_sigma(this)
		class(viscosity_solver)	,intent(inout)	::	this
		
		integer         :: i, j, k, dim1, dim2, dim3
        real(dkind)     :: div_v, nu_node, dv1_dx2, dv2_dx1, v_face_h2, v_face_l2

		real(dkind), dimension (3,3)	:: lame_coeffs
		character(len=20)				:: coordinate_system
		
		integer						:: dimensions
		integer		,dimension(3,2)	:: cons_inner_loop, flow_inner_loop
		real(dkind)	,dimension(3)	:: cell_size		
		
		call this%calculate_viscosity_coeff()
		
		dimensions		= this%domain%get_domain_dimensions()

		cons_inner_loop = this%domain%get_local_inner_cells_bounds()

		flow_inner_loop = this%domain%get_local_inner_faces_bounds()

		cell_size		= this%mesh%mesh_ptr%get_cell_edges_length()

		coordinate_system	= this%domain%get_coordinate_system_name()
		
		associate(	v		=> this%v%v_ptr					, &
					nu      => this%nu%s_ptr				, &
					sigma   => this%sigma%t_ptr				, &
					mesh    => this%mesh%mesh_ptr)

		call this%mpi_support%exchange_conservative_vector_field(v)					
					

		!$omp parallel default(none)  private(i,j,k,dim1,dim2,dim3,div_v,dv1_dx2,dv2_dx1,lame_coeffs) , &
		!$omp& firstprivate(this)	,&
		!$omp& shared(mesh,sigma,v,nu,cons_inner_loop,flow_inner_loop,cell_size,dimensions,coordinate_system) 

		!$omp do collapse(3) schedule(static)
		do k = cons_inner_loop(3,1),cons_inner_loop(3,2)
		do j = cons_inner_loop(2,1),cons_inner_loop(2,2)
		do i = cons_inner_loop(1,1),cons_inner_loop(1,2)
		
            do dim1 = 1,dimensions
            do dim2 = 1,dimensions
  
                div_v = 0.0_dkind
			
				lame_coeffs		= 1.0_dkind				
			
				select case(coordinate_system)
					case ('cartesian')	
						lame_coeffs			= 1.0_dkind
					case ('cylindrical')
						! x -> r, y -> z
						lame_coeffs(1,1)	=  mesh%mesh(1,i,j,k) - cell_size(1)			
						lame_coeffs(1,2)	=  mesh%mesh(1,i,j,k)
						lame_coeffs(1,3)	=  mesh%mesh(1,i,j,k) + cell_size(1)	

					case ('spherical')
						! x -> r
						lame_coeffs(1,1)	=  (mesh%mesh(1,i,j,k) - cell_size(1))**2
						lame_coeffs(1,2)	=  (mesh%mesh(1,i,j,k))**2
						lame_coeffs(1,3)	=  (mesh%mesh(1,i,j,k) + cell_size(1))**2
				end select					

                if (dim1 <= dim2) then
                    if (dim1 == dim2) then
                        do dim3 = 1,dimensions
                            div_v = div_v + (v%pr(dim3)%cells(i+I_m(dim3,1),j+I_m(dim3,2),k+I_m(dim3,3)) * lame_coeffs(dim3,3) - v%pr(dim3)%cells(i-I_m(dim3,1),j-I_m(dim3,2),k-I_m(dim3,3)) * lame_coeffs(dim3,1))/(2.0_dkind*cell_size(dim3)*lame_coeffs(dim3,2))
                        end do
                    else
                        div_v = 0.0_dkind
                    end if

                    dv1_dx2 = (v%pr(dim1)%cells(i+I_m(dim2,1),j+I_m(dim2,2),k+I_m(dim2,3)) - v%pr(dim1)%cells(i-I_m(dim2,1),j-I_m(dim2,2),k-I_m(dim2,3))) / (2.0_dkind*cell_size(dim2))
                    dv2_dx1 = (v%pr(dim2)%cells(i+I_m(dim1,1),j+I_m(dim1,2),k+I_m(dim1,3)) - v%pr(dim2)%cells(i-I_m(dim1,1),j-I_m(dim1,2),k-I_m(dim1,3))) / (2.0_dkind*cell_size(dim1))

                    sigma%pr(dim1,dim2)%cells(i,j,k) = nu%cells(i,j,k) * (dv1_dx2 + dv2_dx1 - 2.0_dkind/3.0_dkind*div_v)

                    select case(coordinate_system)
						case ('cylindrical')
							this%sigma_theta_theta(i,j,k) = -2.0_dkind/3.0_dkind*div_v
							this%sigma_theta_theta(i,j,k) = nu%cells(i,j,k) * (this%sigma_theta_theta(i,j,k) + 2.0_dkind*v%pr(1)%cells(i,j,k)/mesh%mesh(1,i,j,k))
						case ('spherical')
							this%sigma_theta_theta(i,j,k) = -2.0_dkind/3.0_dkind*div_v
							this%sigma_theta_theta(i,j,k) = nu%cells(i,j,k) * (this%sigma_theta_theta(i,j,k) + 2.0_dkind*v%pr(1)%cells(i,j,k)/mesh%mesh(1,i,j,k)) 
					end select						
                
					if (dim1 /= dim2) sigma%pr(dim2,dim1)%cells(i,j,k) = sigma%pr(dim1,dim2)%cells(i,j,k)
                    
                end if
			end do
            end do
        end do
        end do
        end do
		!$omp end do nowait

		!$omp end parallel
					
		end associate
					
	end subroutine
	
	subroutine calculate_viscosity_coeff_constant(this)
		class(viscosity_solver) ,intent(inout) :: this

		real(dkind)	:: reduced_collision_diameter
		real(dkind)	:: inv_reduced_molar_mass

		integer		:: species_number
		integer		:: specie_number1, specie_number2

		species_number = this%chem%chem_ptr%species_number

		associate(  potential_well_depth    => this%thermo%thermo_ptr%potential_well_depth              , &
					molar_masses            => this%thermo%thermo_ptr%molar_masses                      , &
					collision_diameter      => this%thermo%thermo_ptr%collision_diameter)

			do specie_number1 = 1,species_number
				this%viscosity_coeff_constant(specie_number1)    = 0.00001_dkind * 0.0844_dkind * sqrt(molar_masses(specie_number1))/collision_diameter(specie_number1)/collision_diameter(specie_number1)
			end do

		end associate

	end subroutine	
	
	subroutine calculate_viscosity_coeff(this)
		class(viscosity_solver) ,intent(inout) :: this

		real(dkind)                     :: mol_frac, smc, reduced_temperature, sum1, sum2
		real(dkind)                     :: specie_cp, specie_cv
		real(dkind)                     :: omega_2_2
		
		integer	:: dimensions
		integer	:: sign, bound_number
		integer	:: species_number

		integer	,dimension(3,2)	:: cons_inner_loop			
		
		integer :: specie_number
		integer :: i,j,k,plus,dim,dim1,dim2

		species_number	= this%chem%chem_ptr%species_number
		dimensions		= this%domain%get_domain_dimensions()
		
		cons_inner_loop = this%domain%get_local_inner_cells_bounds()

		associate(  T                       => this%T%s_ptr                            , &
					nu                      => this%nu%s_ptr                           , &
					mol_mix_conc            => this%mol_mix_conc%s_ptr                 , &
					Y						=> this%Y%v_ptr									, &
					potential_well_depth    => this%thermo%thermo_ptr%potential_well_depth	, &
					molar_masses            => this%thermo%thermo_ptr%molar_masses			, &
					collision_diameter      => this%thermo%thermo_ptr%collision_diameter	, &
					bc						=> this%boundary%bc_ptr)

	!$omp parallel default(none)  private(i,j,k,specie_number,sum1,sum2,mol_frac,reduced_temperature,omega_2_2,specie_cp,specie_cv,smc,sign,bound_number) , &
	!$omp& firstprivate(this)	,&
	!$omp& shared(T,nu,mol_mix_conc,Y,collision_diameter,molar_masses,potential_well_depth,species_number,cons_inner_loop,bc,dimensions) 
	!$omp do collapse(3) schedule(static)

		do k = cons_inner_loop(3,1),cons_inner_loop(3,2)
		do j = cons_inner_loop(2,1),cons_inner_loop(2,2)
		do i = cons_inner_loop(1,1),cons_inner_loop(1,2)

			if(bc%bc_markers(i,j,k) == 0) then

				sum1 = 0.0_dkind
				sum2 = 0.0_dkind

				do specie_number = 1,species_number
					if (molar_masses(specie_number) /= 0.0_dkind) then				
						mol_frac            =	Y%pr(specie_number)%cells(i,j,k)/molar_masses(specie_number) * mol_mix_conc%cells(i,j,k)
						if (mol_frac /= 0.0_dkind) then
							reduced_temperature =	T%cells(i,j,k) / potential_well_depth(specie_number)
							if (reduced_temperature < 70.0_dkind) then
								omega_2_2           =	1.16145_dkind / (reduced_temperature ** 0.14874_dkind)    +   &
														0.52487_dkind / exp(0.77320_dkind * reduced_temperature)  +   &
														2.16178_dkind / exp(2.43787_dkind * reduced_temperature)
							else
								omega_2_2           =	1.16145_dkind / (reduced_temperature ** 0.14874_dkind)
							end if
							specie_cp = this%thermo%thermo_ptr%calculate_specie_cp(T%cells(i,j,k), specie_number)
							specie_cv = specie_cp - r_gase_J

							smc = this%viscosity_coeff_constant(specie_number)   * sqrt(T%cells(i,j,k)) / omega_2_2

							sum1 = sum1 + smc * mol_frac
							sum2 = sum2 + mol_frac / smc
						end if
					end if
				end do

				if (sum2 <= 1.0E-10_dkind) then	
					nu%cells(i,j,k)      = 0.0_dkind
				else
					nu%cells(i,j,k)      = 0.5_dkind * (sum1 + 1.0_dkind / sum2)
				end if
				
				if(bc%bc_markers(i,j,k) == 0) then
					do dim = 1,dimensions
						do plus = 1,2
							sign			= (-1)**plus
							bound_number	= bc%bc_markers(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3))
							if( bound_number /= 0 ) then
								do dim1 = 1,dimensions
								do dim2 = 1,dimensions
									if (dim1 == dim2) then
										nu%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3))	= nu%cells(i,j,k)
									end if
								end do
								end do
							end if
						end do
					end do
				end if				
				
			end if
		end do
		end do
		end do

	!$omp end do nowait
	!$omp end parallel

		continue
		
		end associate

	end subroutine

	subroutine apply_boundary_conditions(this)
	
		class(viscosity_solver)		,intent(inout)		:: this

		integer					:: dimensions
		integer	,dimension(3,2)	:: cons_inner_loop		
		
		integer	:: sign, bound_number
		integer :: i,j,k,plus,dim,dim1,dim2
			
		dimensions		= this%domain%get_domain_dimensions()

		cons_inner_loop = this%domain%get_local_inner_cells_bounds()

		associate(  sigma		=> this%sigma%t_ptr			, &
					bc			=> this%boundary%bc_ptr	, &
					mesh		=> this%mesh%mesh_ptr)

		!$omp parallel default(none)  private(i,j,k,plus,dim,sign,bound_number) , &
		!$omp& firstprivate(this)	,&
		!$omp& shared(sigma,bc,dimensions,cons_inner_loop) 
		!$omp do collapse(3) schedule(static)

			do k = cons_inner_loop(3,1),cons_inner_loop(3,2)
			do j = cons_inner_loop(2,1),cons_inner_loop(2,2)
			do i = cons_inner_loop(1,1),cons_inner_loop(1,2)
				if(bc%bc_markers(i,j,k) == 0) then
					do dim = 1,dimensions
						do plus = 1,2
							sign			= (-1)**plus
							bound_number	= bc%bc_markers(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3))
							if( bound_number /= 0 ) then
								do dim1 = 1,dimensions
								do dim2 = 1,dimensions
								!	if (dim1 == dim2) then
										sigma%pr(dim1,dim2)%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3))	= 2*sigma%pr(dim1,dim2)%cells(i,j,k) - sigma%pr(dim1,dim2)%cells(i-sign*I_m(dim,1),j-sign*I_m(dim,2),k-sign*I_m(dim,3))
								!	end if
								end do
								end do
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
