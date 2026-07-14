module fourier_heat_transfer_solver_class

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
	public	:: heat_transfer_solver, heat_transfer_solver_c

	type(field_scalar_cons)	,target	:: E_f_prod_heat, kappa

	type	:: heat_transfer_solver
		type(field_scalar_cons_pointer)				:: E_f_prod, T, kappa, mix_mol_mass, rho
		type(field_vector_cons_pointer)				:: Y
		type(computational_domain)					:: domain
		type(mpi_communications)					:: mpi_support
		type(boundary_conditions_pointer)			:: boundary
		type(computational_mesh_pointer)			:: mesh
		type(thermophysical_properties_pointer)		:: thermo
		type(chemical_properties_pointer)			:: chem

	contains
		procedure	,private	::	calculate_thermal_c_coeff
		procedure	,private	::	apply_boundary_conditions
		procedure				::	solve_heat_transfer
	end type

	interface	heat_transfer_solver_c
		module procedure	constructor
	end interface

contains

	type(heat_transfer_solver)	function constructor(manager)

		type(data_manager)	, intent(inout)	:: manager

		type(field_scalar_cons_pointer)	:: scal_ptr
		type(field_vector_cons_pointer)	:: vect_ptr
		type(field_tensor_cons_pointer)	:: tens_ptr

		call manager%get_cons_field_pointer_by_name(scal_ptr,vect_ptr,tens_ptr,'temperature')
		constructor%T%s_ptr					=> scal_ptr%s_ptr
		call manager%get_cons_field_pointer_by_name(scal_ptr,vect_ptr,tens_ptr,'mixture_molar_mass')
		constructor%mix_mol_mass%s_ptr		=> scal_ptr%s_ptr
		call manager%get_cons_field_pointer_by_name(scal_ptr,vect_ptr,tens_ptr,'density')
		constructor%rho%s_ptr				=> scal_ptr%s_ptr

		call manager%create_scalar_field(E_f_prod_heat	,'energy_production_heat_transfer'	,'E_f_prod_heat')
		constructor%E_f_prod%s_ptr		=> E_f_prod_heat
		call manager%create_scalar_field(kappa			,'thermal_conductivity'				,'kappa')
		constructor%kappa%s_ptr			=> kappa

		call manager%get_cons_field_pointer_by_name(scal_ptr,vect_ptr,tens_ptr,'specie_mass_fraction')
		constructor%Y%v_ptr				=> vect_ptr%v_ptr

		constructor%mesh%mesh_ptr	=> manager%computational_mesh_pointer%mesh_ptr
		constructor%boundary%bc_ptr => manager%boundary_conditions_pointer%bc_ptr

		constructor%domain				= manager%domain
		constructor%mpi_support			= manager%mpi_communications

		constructor%thermo%thermo_ptr	=> manager%thermophysics%thermo_ptr
		constructor%chem%chem_ptr		=> manager%chemistry%chem_ptr

		constructor%E_f_prod%s_ptr%cells(:,:,:) = 0.0_dp
		
	end function

	subroutine	solve_heat_transfer(this,time_step)

		class(heat_transfer_solver) ,intent(inout) :: this
		real(dp)				,intent(in)		:: time_step

		real(dp)	:: div_thermo_flux, thermo_flux1, thermo_flux2
		real(dp)	:: energy_prod_summ

		real(dp), dimension (3,3)	:: lame_coeffs		
		
		integer						:: dimensions
		integer		,dimension(3,2)	:: cons_inner_loop
		real(dp)	,dimension(3)	:: cell_size
		character(len=20)	:: coordinate_system
		
		character(len=20)		:: boundary_type_name
		integer	:: sign, bound_number
		integer :: i,j,k,plus,dim

		call this%calculate_thermal_c_coeff()
		call this%apply_boundary_conditions()

		dimensions		= this%domain%get_domain_dimensions()

		cons_inner_loop = this%domain%get_local_inner_cells_bounds()

		cell_size		= this%mesh%mesh_ptr%get_cell_edges_length()

		coordinate_system	= this%domain%get_coordinate_system_name()
		
        associate(	T			=> this%T%s_ptr			, &
		        	kappa		=> this%kappa%s_ptr )
            
		call this%mpi_support%exchange_conservative_scalar_field(kappa)
		call this%mpi_support%exchange_conservative_scalar_field(T)					
				
        end associate
		
		associate(  E_f_prod	=> this%E_f_prod%s_ptr	, &
					rho			=> this%rho%s_ptr		, &
					T			=> this%T%s_ptr			, &
					kappa		=> this%kappa%s_ptr		, &
					bc			=> this%boundary%bc_ptr , &
					mesh		=> this%mesh%mesh_ptr)

	!$omp parallel default(shared)  private(i,j,k,dim,div_thermo_flux,thermo_flux1,thermo_flux2,lame_coeffs,sign,bound_number,boundary_type_name) !, &
	!!$omp& shared(this,time_step,cons_inner_loop,dimensions,cell_size,coordinate_system)
       
	!$omp do collapse(3) schedule(static)
		do k = cons_inner_loop(3,1),cons_inner_loop(3,2)
		do j = cons_inner_loop(2,1),cons_inner_loop(2,2)
		do i = cons_inner_loop(1,1),cons_inner_loop(1,2)

			E_f_prod%cells(i,j,k) = 0.0_dp
		
			if(bc%bc_markers(i,j,k) == 0) then
				
				div_thermo_flux = 0.0_dp
				
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
				
                do dim = 1,dimensions
					thermo_flux1	= 0.5_dp * (kappa%cells(i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3)) + kappa%cells(i,j,k)) * lame_coeffs(dim,1)  &
												* (T%cells(i,j,k) - T%cells(i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3))) / cell_size(dim)

					thermo_flux2	= 0.5_dp * (kappa%cells(i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3)) + kappa%cells(i,j,k)) * lame_coeffs(dim,3)  &
												* (T%cells(i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3)) - T%cells(i,j,k)) / cell_size(dim)

					div_thermo_flux = div_thermo_flux  + (thermo_flux2 - thermo_flux1) / cell_size(dim) / lame_coeffs(dim,2)
                end do
				E_f_prod%cells(i,j,k)	=  div_thermo_flux 
			
				do dim = 1,dimensions
					do plus = 1,2
						sign			= (-1)**plus
						bound_number	= bc%bc_markers(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3))
						if( bound_number /= 0 ) then
							boundary_type_name = bc%boundary_types(bound_number)%get_type_name()
							select case(boundary_type_name)
								case('inlet','outlet')
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
       
        associate(  E_f_prod	=> this%E_f_prod%s_ptr)
		call this%mpi_support%exchange_conservative_scalar_field(E_f_prod)	
        end associate
		
		
	end subroutine


	subroutine calculate_thermal_c_coeff(this)
		class(heat_transfer_solver) ,intent(inout) :: this

        real(dp), dimension(this%chem%chem_ptr%species_number) :: Y_cell

		integer	,dimension(3,2)	:: cons_inner_loop
		integer	:: species_number
		integer :: specie_number
		integer :: i,j,k

		species_number	= this%chem%chem_ptr%species_number
		cons_inner_loop = this%domain%get_local_inner_cells_bounds()

		associate(  T                       => this%T%s_ptr                            , &
					kappa					=> this%kappa%s_ptr							, &
					mix_mol_mass            => this%mix_mol_mass%s_ptr                 , &
					Y						=> this%Y%v_ptr									, &
					bc						=> this%boundary%bc_ptr)

	!$omp parallel default(shared) private(i,j,k,specie_number,Y_cell)
	!$omp do collapse(3) schedule(static)

		do k = cons_inner_loop(3,1),cons_inner_loop(3,2)
		do j = cons_inner_loop(2,1),cons_inner_loop(2,2)
		do i = cons_inner_loop(1,1),cons_inner_loop(1,2)

			if(bc%bc_markers(i,j,k) == 0) then
                do specie_number = 1,species_number
                    Y_cell(specie_number) = Y%pr(specie_number)%cells(i,j,k)
                end do

                kappa%cells(i,j,k) = this%thermo%thermo_ptr%mixture_thermal_conductivity( &
                                         T%cells(i,j,k),Y_cell,mix_mol_mass%cells(i,j,k))
			end if
		end do
		end do
		end do

	!$omp end do
	!$omp end parallel

		end associate

	end subroutine

	subroutine apply_boundary_conditions(this)

		class(heat_transfer_solver)		,intent(inout)		:: this

		integer					:: dimensions
		integer	,dimension(3,2)	:: cons_inner_loop
		character(len=20)		:: boundary_type_name
		real(dp)				:: wall_conductivity_ratio

		integer	:: sign, bound_number
		integer :: i,j,k,plus,dim

		dimensions		= this%domain%get_domain_dimensions()

		cons_inner_loop = this%domain%get_local_inner_cells_bounds()

		associate(  T			=> this%T%s_ptr			, &
					kappa		=> this%kappa%s_ptr		, &
					bc			=> this%boundary%bc_ptr	, &
					mesh		=> this%mesh%mesh_ptr)

		!$omp parallel default(shared)  private(i,j,k,plus,dim,sign,bound_number,wall_conductivity_ratio,boundary_type_name) !, &
		!!$omp& shared(this,cons_inner_loop,dimensions)
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
								kappa%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3))	= - kappa%cells(i,j,k)
								boundary_type_name = bc%boundary_types(bound_number)%get_type_name()
								select case(boundary_type_name)
									case('wall')
										if(bc%boundary_types(bound_number)%is_conductive()) then	
											wall_conductivity_ratio = bc%boundary_types(bound_number)%get_wall_conductivity_ratio()
											kappa%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3))	= wall_conductivity_ratio * kappa%cells(i,j,k)
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
