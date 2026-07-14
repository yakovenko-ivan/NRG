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
	type(field_vector_cons)	,target	:: Y_prod_diff, D, D_T
	
	type	:: diffusion_solver
		type(field_scalar_cons_pointer)				:: T, p, mix_mol_mass, rho, E_f_prod
		type(field_vector_cons_pointer)				:: Y, Y_prod, D, D_T
		type(computational_domain)					:: domain
		type(mpi_communications)					:: mpi_support
		type(boundary_conditions_pointer)			:: boundary
		type(computational_mesh_pointer)			:: mesh
		type(thermophysical_properties_pointer)		:: thermo
		type(chemical_properties_pointer)			:: chem
		
        integer :: H_index = 0, H2_index = 0
        logical :: soret_enabled = .true.

	contains
		procedure	,private	::	calculate_diffusivity_coeff
        procedure   ,private    ::  calculate_thermal_diffusion_coeff
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
		call manager%get_cons_field_pointer_by_name(scal_ptr,vect_ptr,tens_ptr,'mixture_molar_mass')
		constructor%mix_mol_mass%s_ptr		=> scal_ptr%s_ptr
		
		call manager%get_cons_field_pointer_by_name(scal_ptr,vect_ptr,tens_ptr,'specie_mass_fraction')
		constructor%Y%v_ptr				=> vect_ptr%v_ptr

		call manager%create_scalar_field(E_f_prod_diff	,'energy_production_diffusion'	,'E_f_prod_diff')
		constructor%E_f_prod%s_ptr		=> E_f_prod_diff		
		call manager%create_vector_field(Y_prod_diff,'specie_production_diffusion'	,'Y_prod_diff'	,'chemical')
		constructor%Y_prod%v_ptr		=> Y_prod_diff		
		call manager%create_vector_field(D			,'diffusivity'					,'D'			,'chemical')
		constructor%D%v_ptr				=> D		
        call manager%create_vector_field(D_T         ,'thermal_diffusion_coefficient'  ,'D_T'         ,'chemical')
        constructor%D_T%v_ptr            => D_T
		
		constructor%mesh%mesh_ptr	=> manager%computational_mesh_pointer%mesh_ptr
		constructor%boundary%bc_ptr => manager%boundary_conditions_pointer%bc_ptr		
		
		constructor%domain				= manager%domain
		constructor%mpi_support			= manager%mpi_communications

		constructor%thermo%thermo_ptr	=> manager%thermophysics%thermo_ptr
		constructor%chem%chem_ptr		=> manager%chemistry%chem_ptr

		species_number = constructor%chem%chem_ptr%species_number
        constructor%H_index  = constructor%chem%chem_ptr%get_chemical_specie_index('H')
        constructor%H2_index = constructor%chem%chem_ptr%get_chemical_specie_index('H2')
		
		constructor%E_f_prod%s_ptr%cells(:,:,:)		= 0.0_dp
		do spec = 1, species_number
		constructor%Y_prod%v_ptr%pr(spec)%cells(:,:,:)	= 0.0_dp		
            constructor%D%v_ptr%pr(spec)%cells(:,:,:)      = 0.0_dp
            constructor%D_T%v_ptr%pr(spec)%cells(:,:,:)    = 0.0_dp
		end do
		

	end function

	subroutine	solve_diffusion(this,time_step)
	
		class(diffusion_solver) ,intent(inout) :: this
		real(dp)				,intent(in)		:: time_step

        real(dp) :: div_dif_flux, diffusion_flux1, diffusion_flux2
        real(dp) :: diff_velocity_corr1, diff_velocity_corr2
        real(dp) :: soret_flux_corr1, soret_flux_corr2
        real(dp) :: Y_face_sum1, Y_face_sum2, Y_face1, Y_face2
        real(dp) :: Y_corr1, Y_corr2, rho_face1, rho_face2
        real(dp) :: grad_log_T1, grad_log_T2
        real(dp) :: specie_enthalpy, specie_enthalpy1, specie_enthalpy2, h_s_Tref		

		real(dp) ,dimension(this%chem%chem_ptr%species_number) :: diff_velocity1, diff_velocity2
        real(dp), dimension(this%chem%chem_ptr%species_number) :: soret_raw_flux1, soret_raw_flux2
		
		real(dp), dimension (3,3)	:: lame_coeffs	
		
		integer						:: dimensions, species_number
		integer		,dimension(3,2)	:: cons_inner_loop
		real(dp)	,dimension(3)	:: cell_size			
		character(len=20)			:: coordinate_system
		
        integer :: bound_number1, bound_number2
        logical :: face1_is_boundary, face2_is_boundary
        integer :: i,j,k,dim,specie_number
		
		call this%calculate_diffusivity_coeff()
        call this%calculate_thermal_diffusion_coeff()
		call this%apply_boundary_conditions()
		
		dimensions		= this%domain%get_domain_dimensions()
		species_number	= this%chem%chem_ptr%species_number
		cons_inner_loop = this%domain%get_local_inner_cells_bounds()
		cell_size		= this%mesh%mesh_ptr%get_cell_edges_length()
		coordinate_system	= this%domain%get_coordinate_system_name()
		
        associate(  rho				=> this%rho%s_ptr		, &
				    T				=> this%T%s_ptr			, &		
				    Y				=> this%Y%v_ptr			, &				
				    D				=> this%D%v_ptr			, &
                    D_T          	=> this%D_T%v_ptr       , &
				    mix_mol_mass    => this%mix_mol_mass%s_ptr)
        
		    call this%mpi_support%exchange_conservative_scalar_field(mix_mol_mass)
		    call this%mpi_support%exchange_conservative_scalar_field(rho)
		    call this%mpi_support%exchange_conservative_scalar_field(T)
		    call this%mpi_support%exchange_conservative_vector_field(D)
            call this%mpi_support%exchange_conservative_vector_field(D_T)
		    call this%mpi_support%exchange_conservative_vector_field(Y)
        end associate
    		
		associate(  rho				=> this%rho%s_ptr		, &
					T				=> this%T%s_ptr			, &		
					Y				=> this%Y%v_ptr			, &
					E_f_prod		=> this%E_f_prod%s_ptr	, &
					Y_prod			=> this%Y_prod%v_ptr	, &
					D				=> this%D%v_ptr			, &
                    D_T          	=> this%D_T%v_ptr       , &
					molar_masses    => this%thermo%thermo_ptr%molar_masses		, &
					mesh			=> this%mesh%mesh_ptr	, &
					bc				=> this%boundary%bc_ptr)

    !$omp parallel default(shared) private(i,j,k,dim,specie_number,div_dif_flux,diff_velocity1,diff_velocity2, &
    !$omp& diff_velocity_corr1,diff_velocity_corr2,soret_raw_flux1,soret_raw_flux2,soret_flux_corr1,soret_flux_corr2, &
    !$omp& diffusion_flux1,diffusion_flux2,h_s_Tref,specie_enthalpy,specie_enthalpy1,specie_enthalpy2,lame_coeffs, &
    !$omp& bound_number1,bound_number2,face1_is_boundary,face2_is_boundary,Y_face_sum1,Y_face_sum2,Y_face1,Y_face2, &
    !$omp& Y_corr1,Y_corr2,rho_face1,rho_face2,grad_log_T1,grad_log_T2)
    	    
    !$omp do collapse(3) schedule(static)
		do k = cons_inner_loop(3,1),cons_inner_loop(3,2)
		do j = cons_inner_loop(2,1),cons_inner_loop(2,2)
		do i = cons_inner_loop(1,1),cons_inner_loop(1,2)
		
			E_f_prod%cells(i,j,k) = 0.0_dp
			do specie_number = 1,species_number
				Y_prod%pr(specie_number)%cells(i,j,k) = 0.0_dp
			end do		
		
            if (bc%bc_markers(i,j,k) == 0) then

				do dim = 1,dimensions
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
						
                    bound_number1 = bc%bc_markers(i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3))
                    bound_number2 = bc%bc_markers(i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3))
                    face1_is_boundary = (bound_number1 /= 0)
                    face2_is_boundary = (bound_number2 /= 0)

                    diff_velocity1 = 0.0_dp
                    diff_velocity2 = 0.0_dp
                    soret_raw_flux1   = 0.0_dp
                    soret_raw_flux2   = 0.0_dp
                    diff_velocity_corr1 = 0.0_dp
                    diff_velocity_corr2 = 0.0_dp
                    soret_flux_corr1    = 0.0_dp
                    soret_flux_corr2    = 0.0_dp
                    Y_face_sum1         = 0.0_dp
                    Y_face_sum2         = 0.0_dp
                    grad_log_T1         = 0.0_dp
                    grad_log_T2         = 0.0_dp

                    if (.not. face1_is_boundary) then
                        grad_log_T1 = lame_coeffs(dim,1) * &
                            (log(max(T%cells(i,j,k),tiny(1.0_dp))) - &
                             log(max(T%cells(i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3)),tiny(1.0_dp)))) / cell_size(dim)
                    end if
                    if (.not. face2_is_boundary) then
                        grad_log_T2 = lame_coeffs(dim,3) * &
                            (log(max(T%cells(i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3)),tiny(1.0_dp))) - &
                             log(max(T%cells(i,j,k),tiny(1.0_dp)))) / cell_size(dim)
                    end if

					do specie_number = 1,species_number
                        if (molar_masses(specie_number) <= 0.0_dp) cycle

                        Y_face1 = 0.5_dp * (Y%pr(specie_number)%cells(i,j,k) + &
                                           Y%pr(specie_number)%cells(i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3)))
                        Y_face2 = 0.5_dp * (Y%pr(specie_number)%cells(i,j,k) + &
                                           Y%pr(specie_number)%cells(i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3)))
                        Y_face_sum1 = Y_face_sum1 + Y_face1
                        Y_face_sum2 = Y_face_sum2 + Y_face2

                        if (.not. face1_is_boundary) then
                            diff_velocity1(specie_number) = 0.5_dp * &
                                (D%pr(specie_number)%cells(i,j,k) + &
                                 D%pr(specie_number)%cells(i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3))) * &
                                lame_coeffs(dim,1) * &
                                (Y%pr(specie_number)%cells(i,j,k) - &
                                 Y%pr(specie_number)%cells(i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3))) / cell_size(dim)

                            soret_raw_flux1(specie_number) = 0.5_dp * &
                                (D_T%pr(specie_number)%cells(i,j,k) + &
                                 D_T%pr(specie_number)%cells(i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3))) * grad_log_T1
																		
							diff_velocity_corr1	= diff_velocity_corr1 + diff_velocity1(specie_number)
                            soret_flux_corr1    = soret_flux_corr1    + soret_raw_flux1(specie_number)
                        end if

                        if (.not. face2_is_boundary) then
                            diff_velocity2(specie_number) = 0.5_dp * &
                                (D%pr(specie_number)%cells(i,j,k) + &
                                 D%pr(specie_number)%cells(i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3))) * &
                                lame_coeffs(dim,3) * &
                                (Y%pr(specie_number)%cells(i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3)) - &
                                 Y%pr(specie_number)%cells(i,j,k)) / cell_size(dim)

                            soret_raw_flux2(specie_number) = 0.5_dp * &
                                (D_T%pr(specie_number)%cells(i,j,k) + &
                                 D_T%pr(specie_number)%cells(i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3))) * grad_log_T2
					
							diff_velocity_corr2	= diff_velocity_corr2 + diff_velocity2(specie_number)	
                            soret_flux_corr2    = soret_flux_corr2    + soret_raw_flux2(specie_number)
						end if
					end do
						
                    rho_face1 = 0.5_dp * (rho%cells(i,j,k) + &
                                          rho%cells(i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3)))
                    rho_face2 = 0.5_dp * (rho%cells(i,j,k) + &
                                          rho%cells(i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3)))

                    do specie_number = 1,species_number
                        if (molar_masses(specie_number) <= 0.0_dp) cycle

                        Y_face1 = 0.5_dp * (Y%pr(specie_number)%cells(i,j,k) + &
                                           Y%pr(specie_number)%cells(i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3)))
                        Y_face2 = 0.5_dp * (Y%pr(specie_number)%cells(i,j,k) + &
                                           Y%pr(specie_number)%cells(i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3)))

                        Y_corr1 = 0.0_dp
                        Y_corr2 = 0.0_dp
                        if (Y_face_sum1 > 1.0e-30_dp) Y_corr1 = Y_face1/Y_face_sum1
                        if (Y_face_sum2 > 1.0e-30_dp) Y_corr2 = Y_face2/Y_face_sum2

                        diffusion_flux1 = 0.0_dp
                        diffusion_flux2 = 0.0_dp

                        if (.not. face1_is_boundary) then
                            diffusion_flux1 = rho_face1 * &
                                (diff_velocity1(specie_number) - diff_velocity_corr1*Y_corr1) + &
                                soret_raw_flux1(specie_number) - soret_flux_corr1*Y_corr1
                        end if

                        if (.not. face2_is_boundary) then
                            diffusion_flux2 = rho_face2 * &
                                (diff_velocity2(specie_number) - diff_velocity_corr2*Y_corr2) + &
                                soret_raw_flux2(specie_number) - soret_flux_corr2*Y_corr2
                        end if

                        div_dif_flux = (diffusion_flux2-diffusion_flux1) / cell_size(dim) / lame_coeffs(dim,2)
                        Y_prod%pr(specie_number)%cells(i,j,k) = &
                            Y_prod%pr(specie_number)%cells(i,j,k) + div_dif_flux

                        h_s_Tref        = this%thermo%thermo_ptr%specie_enthalpy_molar(T_ref,specie_number)
                        specie_enthalpy = this%thermo%thermo_ptr%specie_enthalpy_molar( &
                                               T%cells(i,j,k),specie_number) - h_s_Tref

                        specie_enthalpy1 = specie_enthalpy
                        specie_enthalpy2 = specie_enthalpy
                        if (.not. face1_is_boundary) then
                            specie_enthalpy1 = this%thermo%thermo_ptr%specie_enthalpy_molar( &
                                T%cells(i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3)),specie_number) - h_s_Tref
                        end if
                        if (.not. face2_is_boundary) then
                            specie_enthalpy2 = this%thermo%thermo_ptr%specie_enthalpy_molar( &
                                T%cells(i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3)),specie_number) - h_s_Tref
                        end if

                        div_dif_flux = (diffusion_flux2*0.5_dp*(specie_enthalpy+specie_enthalpy2) - &
                                        diffusion_flux1*0.5_dp*(specie_enthalpy+specie_enthalpy1)) / &
                                        cell_size(dim) / lame_coeffs(dim,2)

                        E_f_prod%cells(i,j,k) = E_f_prod%cells(i,j,k) + div_dif_flux/molar_masses(specie_number)
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
		
	end subroutine solve_diffusion

	subroutine calculate_diffusivity_coeff(this)
		class(diffusion_solver), intent(inout) :: this

        real(dp), dimension(this%chem%chem_ptr%species_number) :: Y_cell, D_cell
        integer, dimension(3,2) :: cons_inner_loop
        integer :: species_number
        integer :: i,j,k,specie_number

        species_number = this%chem%chem_ptr%species_number
        cons_inner_loop = this%domain%get_local_inner_cells_bounds()

        associate( T            => this%T%s_ptr, &
                   p            => this%p%s_ptr, &
                   Y            => this%Y%v_ptr, &
                   D            => this%D%v_ptr, &
                   mix_mol_mass => this%mix_mol_mass%s_ptr, &
                   bc           => this%boundary%bc_ptr)

    !$omp parallel default(shared) private(i,j,k,specie_number,Y_cell,D_cell)
    !$omp do collapse(3) schedule(static)
        do k = cons_inner_loop(3,1),cons_inner_loop(3,2)
        do j = cons_inner_loop(2,1),cons_inner_loop(2,2)
        do i = cons_inner_loop(1,1),cons_inner_loop(1,2)
            D_cell = 0.0_dp
            if (bc%bc_markers(i,j,k) == 0) then
                do specie_number = 1,species_number
                    Y_cell(specie_number) = Y%pr(specie_number)%cells(i,j,k)
                end do

                call this%thermo%thermo_ptr%mixture_averaged_diffusion_coefficients( &
                    T%cells(i,j,k),p%cells(i,j,k),Y_cell, &
                    mix_mol_mass%cells(i,j,k),D_cell)
            end if

            do specie_number = 1,species_number
                D%pr(specie_number)%cells(i,j,k) = D_cell(specie_number)
            end do
        end do
        end do
        end do
    !$omp end do
    !$omp end parallel

        end associate
	end subroutine calculate_diffusivity_coeff

    subroutine calculate_thermal_diffusion_coeff(this)
        class(diffusion_solver), intent(inout) :: this

        real(dp), parameter :: alpha_H  = 0.895_dp
        real(dp), parameter :: alpha_H2 = 0.910_dp

        real(dp), dimension(this%chem%chem_ptr%species_number) :: Y_cell, D_T_cell
        integer, dimension(2) :: target_indices
        real(dp), dimension(2) :: alpha_soret
        integer, dimension(3,2) :: cons_inner_loop
        integer :: species_number
        integer :: i,j,k,specie_number

        species_number = this%chem%chem_ptr%species_number
        cons_inner_loop = this%domain%get_local_inner_cells_bounds()

        target_indices = (/this%H_index,this%H2_index/)
        alpha_soret    = (/alpha_H,alpha_H2/)
        if (.not. this%soret_enabled) target_indices = 0

        associate( T            => this%T%s_ptr, &
                   Y            => this%Y%v_ptr, &
                   D_T          => this%D_T%v_ptr, &
                   mix_mol_mass => this%mix_mol_mass%s_ptr, &
                   bc           => this%boundary%bc_ptr)

    !$omp parallel default(shared) private(i,j,k,specie_number,Y_cell,D_T_cell)
        !$omp do collapse(3) schedule(static)
        do k = cons_inner_loop(3,1),cons_inner_loop(3,2)
        do j = cons_inner_loop(2,1),cons_inner_loop(2,2)
        do i = cons_inner_loop(1,1),cons_inner_loop(1,2)
            D_T_cell = 0.0_dp
            if (bc%bc_markers(i,j,k) == 0 .and. this%soret_enabled) then
                do specie_number = 1,species_number
                    Y_cell(specie_number) = Y%pr(specie_number)%cells(i,j,k)
                end do

                call this%thermo%thermo_ptr%reduced_thermal_diffusion_coefficients( &
                    T%cells(i,j,k),Y_cell,mix_mol_mass%cells(i,j,k), &
                    target_indices,alpha_soret,D_T_cell)
            end if

            do specie_number = 1,species_number
                D_T%pr(specie_number)%cells(i,j,k) = D_T_cell(specie_number)
            end do
        end do
        end do
        end do
    !$omp end do
    !$omp end parallel

        end associate
    end subroutine calculate_thermal_diffusion_coeff

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
                    D_T         => this%D_T%v_ptr            , &
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
                                    D_T%pr(specie_number)%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3)) = D_T%pr(specie_number)%cells(i,j,k)
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