module fickean_diffusion_solver_class

	use kind_parameters
	use global_data
	use field_pointers
	use boundary_conditions_class
	use data_manager_class
	use computational_mesh_class
	use computational_domain_class
	use thermophysical_properties_class
	use chemical_properties_class
	
    use benchmarking
	use mpi_communications_class

	implicit none
	
#ifdef OMP
	include "omp_lib.h"
#endif
	
	private
	public	:: diffusion_solver, diffusion_solver_c
	
	type(field_scalar_cons)	,target	:: E_f_prod_diff
	type(field_vector_cons)	,target	:: Y_prod_diff, D
	
	type	:: diffusion_solver
		type(field_scalar_cons_pointer)				:: T, p, mol_mix_conc, rho, E_f_prod
		type(field_vector_cons_pointer)				:: Y, Y_prod, D
		type(computational_domain)					:: domain
		type(mpi_communications)					:: mpi_support
		type(boundary_conditions_pointer)			:: boundary
		type(computational_mesh_pointer)			:: mesh
		type(thermophysical_properties_pointer)		:: thermo
		type(chemical_properties_pointer)			:: chem
		
		real(dp) ,dimension(:,:)	,allocatable	:: diffusivity_binary_constant

	contains
		procedure	,private	::	calculate_diffusivity_binary_constant
		procedure	,private	::	calculate_diffusivity_coeff
		procedure	,private	::	apply_boundary_conditions
		procedure				::	solve_diffusion
	end type
	
	interface	diffusion_solver_c
		module procedure	constructor
	end interface
	
contains

	type(diffusion_solver)	function constructor(manager)
	
		type(data_manager)	, intent(inout)	:: manager
	
		type(field_scalar_cons_pointer)	:: scal_ptr
		type(field_vector_cons_pointer)	:: vect_ptr
		type(field_tensor_cons_pointer)	:: tens_ptr	
	
		integer	:: species_number
		integer	:: spec
		
		call manager%get_cons_field_pointer_by_name(scal_ptr,vect_ptr,tens_ptr,'pressure')
		constructor%p%s_ptr					=> scal_ptr%s_ptr	
		call manager%get_cons_field_pointer_by_name(scal_ptr,vect_ptr,tens_ptr,'temperature')
		constructor%T%s_ptr					=> scal_ptr%s_ptr	
		call manager%get_cons_field_pointer_by_name(scal_ptr,vect_ptr,tens_ptr,'density')
		constructor%rho%s_ptr				=> scal_ptr%s_ptr	
		call manager%get_cons_field_pointer_by_name(scal_ptr,vect_ptr,tens_ptr,'mixture_molar_concentration')
		constructor%mol_mix_conc%s_ptr		=> scal_ptr%s_ptr
		
		call manager%get_cons_field_pointer_by_name(scal_ptr,vect_ptr,tens_ptr,'specie_molar_concentration')
		constructor%Y%v_ptr				=> vect_ptr%v_ptr

		call manager%create_scalar_field(E_f_prod_diff	,'energy_production_diffusion'	,'E_f_prod_diff')
		constructor%E_f_prod%s_ptr		=> E_f_prod_diff		
		call manager%create_vector_field(Y_prod_diff,'specie_production_diffusion'	,'Y_prod_diff'	,'chemical')
		constructor%Y_prod%v_ptr		=> Y_prod_diff		
		call manager%create_vector_field(D			,'diffusivity'					,'D'			,'chemical')
		constructor%D%v_ptr				=> D		
		
		constructor%mesh%mesh_ptr	=> manager%computational_mesh_pointer%mesh_ptr
		constructor%boundary%bc_ptr => manager%boundary_conditions_pointer%bc_ptr		
		
		constructor%domain				= manager%domain
		constructor%mpi_support			= manager%mpi_communications

		constructor%thermo%thermo_ptr	=> manager%thermophysics%thermo_ptr
		constructor%chem%chem_ptr		=> manager%chemistry%chem_ptr

		species_number = constructor%chem%chem_ptr%species_number
		
		constructor%E_f_prod%s_ptr%cells(:,:,:)		= 0.0_dp
		do spec = 1, species_number
		constructor%Y_prod%v_ptr%pr(spec)%cells(:,:,:)	= 0.0_dp		
		end do
		
		allocate(constructor%diffusivity_binary_constant(manager%chemistry%chem_ptr%species_number,manager%chemistry%chem_ptr%species_number))
		
		call constructor%calculate_diffusivity_binary_constant()

	end function

	subroutine	solve_diffusion(this,time_step)
	
		class(diffusion_solver) ,intent(inout) :: this
		real(dp)				,intent(in)		:: time_step

		real(dp)	:: average_molar_mass_left, average_molar_mass_right, average_molar_mass_middle
		real(dp)	:: div_dif_flux, diffusion_flux1, diffusion_flux2, diffusion_energy_flux1, diffusion_energy_flux2, div_dif_en_flux
		real(dp)	:: diff_velocity_corr1, diff_velocity_corr2
		
		real(dp) ,dimension(:)	,allocatable    :: diff_velocity1, diff_velocity2
		real(dp)					:: specie_enthalpy, specie_enthalpy1, specie_enthalpy2, h_s_Tref
		real(dp)					:: Y_prod_summ
		real(dp)					:: check_summ1, check_summ2
		
		real(dp), dimension (3,3)	:: lame_coeffs	
		
		integer						:: dimensions, species_number
		integer		,dimension(3,2)	:: cons_inner_loop
		real(dp)	,dimension(3)	:: cell_size			
		character(len=20)			:: coordinate_system
		
		character(len=20)			:: boundary_type_name
		integer	:: sign, bound_number
		integer :: i,j,k,plus,dim,specie_number,specie_number1
		
		call this%calculate_diffusivity_coeff()
		call this%apply_boundary_conditions()
		
		dimensions		= this%domain%get_domain_dimensions()

		species_number	= this%chem%chem_ptr%species_number

		cons_inner_loop = this%domain%get_local_inner_cells_bounds()

		cell_size		= this%mesh%mesh_ptr%get_cell_edges_length()

		coordinate_system	= this%domain%get_coordinate_system_name()
		
		allocate(diff_velocity1(species_number),diff_velocity2(species_number))
		diff_velocity1 = 0.0_dp
		diff_velocity2 = 0.0_dp
				
        associate(  rho				=> this%rho%s_ptr		, &
				    T				=> this%T%s_ptr			, &		
				    Y				=> this%Y%v_ptr			, &				
				    D				=> this%D%v_ptr			, &
				    mol_mix_conc    => this%mol_mix_conc%s_ptr)
        
		    call this%mpi_support%exchange_conservative_scalar_field(mol_mix_conc)
		    call this%mpi_support%exchange_conservative_scalar_field(rho)
		    call this%mpi_support%exchange_conservative_scalar_field(T)
		    call this%mpi_support%exchange_conservative_vector_field(D)
		    call this%mpi_support%exchange_conservative_vector_field(Y)
        end associate
    		
		associate(  rho				=> this%rho%s_ptr		, &
					T				=> this%T%s_ptr			, &		
					Y				=> this%Y%v_ptr			, &
					E_f_prod		=> this%E_f_prod%s_ptr	, &
					Y_prod			=> this%Y_prod%v_ptr	, &
					D				=> this%D%v_ptr			, &
					mol_mix_conc    => this%mol_mix_conc%s_ptr , &
					molar_masses    => this%thermo%thermo_ptr%molar_masses		, &
					mesh			=> this%mesh%mesh_ptr	, &
					bc				=> this%boundary%bc_ptr)

	!$omp parallel default(shared)  private(i,j,k,dim,specie_number,div_dif_flux,diff_velocity1, diff_velocity2, diff_velocity_corr1,diff_velocity_corr2,diffusion_flux1,diffusion_flux2, h_s_Tref, specie_enthalpy, specie_enthalpy1, specie_enthalpy2,lame_coeffs,sign,bound_number,boundary_type_name) !, &
	!!$omp& shared(this, time_step,cons_inner_loop,dimensions,species_number, cell_size, coordinate_system) 
	    
    !$omp do collapse(3) schedule(static)
		do k = cons_inner_loop(3,1),cons_inner_loop(3,2)
		do j = cons_inner_loop(2,1),cons_inner_loop(2,2)
		do i = cons_inner_loop(1,1),cons_inner_loop(1,2)
		
			E_f_prod%cells(i,j,k) = 0.0_dp
			do specie_number = 1,species_number
				Y_prod%pr(specie_number)%cells(i,j,k) = 0.0_dp
			end do		
		
			if((bc%bc_markers(i,j,k) == 0)) then

				do dim = 1,dimensions
					div_dif_flux		= 0.0_dp
					diff_velocity_corr1 = 0.0_dp
					diff_velocity_corr2 = 0.0_dp

					lame_coeffs		= 1.0_dp							
						
					select case(coordinate_system)
						case ('cartesian')	
							lame_coeffs			= 1.0_dp
						case ('cylindrical')
							lame_coeffs(1,1)	=  mesh%mesh(1,i,j,k) - 0.5_dp*cell_size(1)
							lame_coeffs(1,2)	=  mesh%mesh(1,i,j,k)
							lame_coeffs(1,3)	=  mesh%mesh(1,i,j,k) + 0.5_dp*cell_size(1)	
						case ('spherical')
							lame_coeffs(1,1)	=  (mesh%mesh(1,i,j,k) - 0.5_dp*cell_size(1))**2
							lame_coeffs(1,2)	=  (mesh%mesh(1,i,j,k))**2
							lame_coeffs(1,3)	=  (mesh%mesh(1,i,j,k) + 0.5_dp*cell_size(1))**2
					end select						
						
					do specie_number = 1,species_number
						if (molar_masses(specie_number) /= 0.0_dp) then
							diff_velocity1(specie_number) = 0.5_dp * (D%pr(specie_number)%cells(i,j,k) + D%pr(specie_number)%cells(i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3))) * lame_coeffs(dim,1) 	* &
															(	Y%pr(specie_number)%cells(i,j,k) - Y%pr(specie_number)%cells(i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3)))  / cell_size(dim) 
				
							diff_velocity2(specie_number) = 0.5_dp * (D%pr(specie_number)%cells(i,j,k) + D%pr(specie_number)%cells(i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3))) * lame_coeffs(dim,3) 	* &
															(	Y%pr(specie_number)%cells(i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3)) - Y%pr(specie_number)%cells(i,j,k))  / cell_size(dim)
																		
							diff_velocity_corr1	= diff_velocity_corr1 + diff_velocity1(specie_number)
					
							diff_velocity_corr2	= diff_velocity_corr2 + diff_velocity2(specie_number)	
						end if
					end do
						
					do specie_number = 1,species_number
						if (molar_masses(specie_number) /= 0.0_dp) then
							diffusion_flux1 =	diff_velocity1(specie_number) - diff_velocity_corr1 * 0.5_dp * (Y%pr(specie_number)%cells(i,j,k) + Y%pr(specie_number)%cells(i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3)))
					
							diffusion_flux1 =	diffusion_flux1 * 0.5_dp * (rho%cells(i,j,k) + rho%cells(i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3)))					
											
							diffusion_flux2 =	diff_velocity2(specie_number) - diff_velocity_corr2 * 0.5_dp * (Y%pr(specie_number)%cells(i,j,k) + Y%pr(specie_number)%cells(i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3)))
							
							diffusion_flux2 =	diffusion_flux2 * 0.5_dp * (rho%cells(i,j,k) + rho%cells(i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3)))
						
							div_dif_flux = (diffusion_flux2 - diffusion_flux1) / cell_size(dim)	/ lame_coeffs(dim,2)						
							
							Y_prod%pr(specie_number)%cells(i,j,k) = Y_prod%pr(specie_number)%cells(i,j,k) + div_dif_flux! * time_step
							
							!specie_enthalpy = (this%thermo%thermo_ptr%calculate_specie_cp(T%cells(i,j,k),specie_number))*T%cells(i,j,k)
							!specie_enthalpy1 = (this%thermo%thermo_ptr%calculate_specie_cp(T%cells(i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3)),specie_number))*T%cells(i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3))
							!specie_enthalpy2 = (this%thermo%thermo_ptr%calculate_specie_cp(T%cells(i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3)),specie_number))*T%cells(i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3))
       
							!div_dif_flux		=  (diffusion_flux2 * 0.5_dp * (specie_enthalpy + specie_enthalpy2) - diffusion_flux1 * 0.5_dp * (specie_enthalpy + specie_enthalpy1)) / cell_size(dim) / lame_coeffs(dim,2)	
                            !print *, div_dif_flux
                            
                            h_s_Tref			=	this%thermo%thermo_ptr%calculate_specie_enthalpy(T_ref, specie_number)
                            specie_enthalpy		=	this%thermo%thermo_ptr%calculate_specie_enthalpy(T%cells(i,j,k), specie_number)	- h_s_Tref
                            specie_enthalpy1	=	this%thermo%thermo_ptr%calculate_specie_enthalpy(T%cells(i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3)), specie_number) - h_s_Tref
                            specie_enthalpy2	=	this%thermo%thermo_ptr%calculate_specie_enthalpy(T%cells(i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3)), specie_number) - h_s_Tref
     
							div_dif_flux		=  (diffusion_flux2 * 0.5_dp * (specie_enthalpy + specie_enthalpy2) - diffusion_flux1 * 0.5_dp * (specie_enthalpy + specie_enthalpy1)) / cell_size(dim) / lame_coeffs(dim,2)	
!							print *, div_dif_flux
                            
							E_f_prod%cells(i,j,k) = E_f_prod%cells(i,j,k) +  div_dif_flux / molar_masses(specie_number)! * time_step 
						end if							
					end do
				end do
		
				do dim = 1,dimensions
					do plus = 1,2
						sign			= (-1)**plus
						bound_number	= bc%bc_markers(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3))
						if( bound_number /= 0 ) then
							boundary_type_name = bc%boundary_types(bound_number)%get_type_name()
							select case(boundary_type_name)
								case('inlet','outlet')
									do specie_number = 1,species_number
										Y_prod%pr(specie_number)%cells(i,j,k) = 0.0_dp
									end do
									E_f_prod%cells(i,j,k) = 0.0_dp
							end select
						end if
					end do
				end do
			
			end if			
		end do
		end do
		end do
	!$omp end do

	!$omp end parallel

        end associate

        associate( E_f_prod		=> this%E_f_prod%s_ptr	, &
                    Y_prod		=> this%Y_prod%v_ptr)
        
		call this%mpi_support%exchange_conservative_scalar_field(E_f_prod)
		call this%mpi_support%exchange_conservative_vector_field(Y_prod)
		
		end associate
		
        deallocate(diff_velocity1,diff_velocity2)
	end subroutine


	subroutine calculate_diffusivity_binary_constant(this)
		class(diffusion_solver) ,intent(inout) :: this

		real(dp)	:: reduced_collision_diameter
		real(dp)	:: inv_reduced_molar_mass
		
		integer		:: species_number
		integer		:: specie_number1, specie_number2

		species_number = this%chem%chem_ptr%species_number

		associate(  potential_well_depth    => this%thermo%thermo_ptr%potential_well_depth              , &
					molar_masses            => this%thermo%thermo_ptr%molar_masses                      , &
					collision_diameter      => this%thermo%thermo_ptr%collision_diameter)

		do specie_number1 = 1,species_number-1
		do specie_number2 = specie_number1+1,species_number
			if ((molar_masses(specie_number1) /= 0.0_dp).and.(molar_masses(specie_number2) /= 0.0_dp)) then
				reduced_collision_diameter  = 0.5_dp * (collision_diameter(specie_number1) + collision_diameter(specie_number2))
				reduced_collision_diameter  = reduced_collision_diameter * reduced_collision_diameter
				inv_reduced_molar_mass      = (molar_masses(specie_number1) + molar_masses(specie_number2))/(molar_masses(specie_number1) * molar_masses(specie_number2))
				this%diffusivity_binary_constant(specie_number1,specie_number2)   = 0.00001_dp * 0.595_dp * sqrt(inv_reduced_molar_mass) / reduced_collision_diameter
				continue
			! Constant is evaluated as 3/16 * sqrt ( 2 * Pi * kB ^ 3 * NA )/ ( Pi * 10 ^ -18 ). Ref. CHEMKIN 
			end if
		end do
		end do

		end associate

	end subroutine	
	
	subroutine calculate_diffusivity_coeff(this)
		class(diffusion_solver) ,intent(inout) :: this

		real(dp)                     :: mol_frac, stc, reduced_temperature, sum1, sum2
		real(dp)                     :: specie_cp, specie_cv
		real(dp)                     :: omega_2_2

		real(dp)						:: average_molar_mass
		real(dp)                     :: accumulation_value
		real(dp)                     :: omega_1_1
		real(dp)                     :: reduced_collision_diameter
		real(dp)                     :: inv_reduced_molar_mass
		real(dp)                     :: specie_mass_density,mixture_mass_density
		
		integer	:: species_number
		integer	,dimension(3,2)	:: cons_inner_loop		
		
		integer :: specie_number1, specie_number2
		integer :: i,j,k

		real(dp) ,dimension(this%chem%chem_ptr%species_number,this%chem%chem_ptr%species_number)   :: D_binary
		
		species_number = this%chem%chem_ptr%species_number

		cons_inner_loop = this%domain%get_local_inner_cells_bounds()

		associate(  T                       => this%T%s_ptr                            , &
					p						=> this%p%s_ptr							, &
					rho						=> this%rho%s_ptr							, &
					D						=> this%D%v_ptr							, &
					mol_mix_conc            => this%mol_mix_conc%s_ptr                 , &
					Y						=> this%Y%v_ptr									, &
					potential_well_depth    => this%thermo%thermo_ptr%potential_well_depth	, &
					molar_masses            => this%thermo%thermo_ptr%molar_masses			, &
                bc                      => this%boundary%bc_ptr		, &
					collision_diameter      => this%thermo%thermo_ptr%collision_diameter)

	!$omp parallel default(shared)  private(i,j,k,sum1,sum2,mol_frac,reduced_temperature,reduced_collision_diameter,inv_reduced_molar_mass,average_molar_mass,accumulation_value,specie_mass_density,mixture_mass_density,omega_1_1,D_binary,specie_number1,specie_number2) !, &
	!!$omp& shared(this,species_number,cons_inner_loop) 
                
	!$omp do collapse(3) schedule(static)

		do k = cons_inner_loop(3,1),cons_inner_loop(3,2)
		do j = cons_inner_loop(2,1),cons_inner_loop(2,2)
		do i = cons_inner_loop(1,1),cons_inner_loop(1,2)
			if(bc%bc_markers(i,j,k) == 0) then

				D_binary = 0.0_dp

				! Calculating upper triangle of matrix D_binary(i,j)	

				do specie_number1 = 1,species_number-1
				do specie_number2 = specie_number1+1,species_number
					if ((molar_masses(specie_number1) /= 0.0_dp).and.(molar_masses(specie_number2) /= 0.0_dp)) then
						reduced_temperature = T%cells(i,j,k) / sqrt(potential_well_depth(specie_number1)*potential_well_depth(specie_number2))
						if (reduced_temperature < 70.0_dp) then
							omega_1_1           = 	1.06036_dp / (reduced_temperature ** 0.15610_dp)      + &
													0.19300_dp / exp(0.47635_dp * reduced_temperature )   + &
													1.03587_dp / exp(1.52996_dp * reduced_temperature )   + &
													1.76474_dp / exp(3.89411_dp * reduced_temperature )
						else
							omega_1_1           = 	1.06036_dp / (reduced_temperature ** 0.15610_dp) 
						end if
					
						reduced_collision_diameter  = 0.5_dp * (collision_diameter(specie_number1) + collision_diameter(specie_number2))
						reduced_collision_diameter  = reduced_collision_diameter * reduced_collision_diameter
						inv_reduced_molar_mass      = (molar_masses(specie_number1) + molar_masses(specie_number2))/(molar_masses(specie_number1) * molar_masses(specie_number2))

						D_binary(specie_number1,specie_number2) = this%diffusivity_binary_constant(specie_number1,specie_number2) * (T%cells(i,j,k) ** 1.5_dp) / p%cells(i,j,k) / omega_1_1
					end if
				end do
				end do

				do specie_number1 = 1,species_number
					
					average_molar_mass = 0.0_dp
					do specie_number2 = 1,species_number
						if (molar_masses(specie_number2) /= 0.0_dp) then
							average_molar_mass = average_molar_mass + (Y%pr(specie_number2)%cells(i,j,k) /molar_masses(specie_number2))
						end if
					end do
					
					average_molar_mass =  1.0_dp / average_molar_mass 
					
					accumulation_value = 0.0_dp
					do specie_number2 = 1,species_number			
						if ((molar_masses(specie_number1) /= 0.0_dp).and.(molar_masses(specie_number2) /= 0.0_dp)) then
							if(specie_number2 > specie_number1) accumulation_value = accumulation_value + (1.0_dp - Y%pr(specie_number1)%cells(i,j,k)) * molar_masses(specie_number1) * Y%pr(specie_number2)%cells(i,j,k) / D_binary(specie_number1,specie_number2) /molar_masses(specie_number2)
							if(specie_number2 < specie_number1) accumulation_value = accumulation_value + (1.0_dp - Y%pr(specie_number1)%cells(i,j,k)) * molar_masses(specie_number1) * Y%pr(specie_number2)%cells(i,j,k) / D_binary(specie_number2,specie_number1) /molar_masses(specie_number2)
						
							if(specie_number2 > specie_number1) accumulation_value = accumulation_value + Y%pr(specie_number1)%cells(i,j,k) * Y%pr(specie_number2)%cells(i,j,k) / D_binary(specie_number1,specie_number2) 
							if(specie_number2 < specie_number1) accumulation_value = accumulation_value + Y%pr(specie_number1)%cells(i,j,k) * Y%pr(specie_number2)%cells(i,j,k) / D_binary(specie_number2,specie_number1) 
						end if
					!	if(specie_number2 > specie_number1) accumulation_value = accumulation_value + Y%pr(specie_number2)%cells(i,j,k)/molar_masses(specie_number2) / D_binary(specie_number1,specie_number2)
					!	if(specie_number2 < specie_number1) accumulation_value = accumulation_value + Y%pr(specie_number2)%cells(i,j,k)/molar_masses(specie_number2) / D_binary(specie_number2,specie_number1)
						
					!	if(specie_number2 > specie_number1) accumulation_value = accumulation_value + Y%pr(specie_number1)%cells(i,j,k)/molar_masses(specie_number1) / (1.0_dp - Y%pr(specie_number1)%cells(i,j,k)) * Y%pr(specie_number2)%cells(i,j,k) / D_binary(specie_number1,specie_number2)
					!	if(specie_number2 < specie_number1) accumulation_value = accumulation_value + Y%pr(specie_number1)%cells(i,j,k)/molar_masses(specie_number1) / (1.0_dp - Y%pr(specie_number1)%cells(i,j,k)) * Y%pr(specie_number2)%cells(i,j,k) / D_binary(specie_number2,specie_number1)

					end do
					
					if (accumulation_value <= 1.0E-10_dp) then
						D%pr(specie_number1)%cells(i,j,k) = 0.0_dp
					else
						if (molar_masses(specie_number1) /= 0.0_dp) then
							D%pr(specie_number1)%cells(i,j,k) = (1.0_dp - Y%pr(specie_number1)%cells(i,j,k)) * molar_masses(specie_number1) / accumulation_value / average_molar_mass
						end if
					!	D%pr(specie_number1)%cells(i,j,k) = 1.0_dp / accumulation_value / average_molar_mass
					end if

				end do
			end if
		end do
		end do
		end do

	!$omp end do 

	!$omp end parallel

    end associate

		continue		
		

	end subroutine

	subroutine apply_boundary_conditions(this)

		class(diffusion_solver) ,intent(inout) :: this
		
		integer					:: dimensions, species_number
		integer	,dimension(3,2)	:: cons_inner_loop
		
		
		integer :: specie_number
		integer	:: sign, bound_number
		integer :: i,j,k,plus,dim

		dimensions		= this%domain%get_domain_dimensions()

		species_number = this%chem%chem_ptr%species_number

		cons_inner_loop = this%domain%get_local_inner_cells_bounds()

		associate(  D			=> this%D%v_ptr				, &
					bc			=> this%boundary%bc_ptr		)
					
		!$omp parallel default(shared)  private(i,j,k,plus,dim,specie_number,sign,bound_number) !, &
		!!$omp& shared(this,cons_inner_loop,dimensions,species_number)

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
								do specie_number = 1,species_number
									D%pr(specie_number)%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3))	=  D%pr(specie_number)%cells(i,j,k)
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