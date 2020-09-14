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

	use mpi_communications_class

	implicit none

#ifdef OMP
	include "omp_lib.h"
#endif

	private
	public	:: heat_transfer_solver, heat_transfer_solver_c

	type(field_scalar_cons)	,target	:: E_f_prod_heat, kappa

	type	:: heat_transfer_solver
		type(field_scalar_cons_pointer)				:: E_f_prod, T, kappa, mol_mix_conc, rho
		type(field_vector_cons_pointer)				:: Y
		type(computational_domain)					:: domain
		type(mpi_communications)					:: mpi_support
		type(boundary_conditions_pointer)			:: boundary
		type(computational_mesh_pointer)			:: mesh
		type(thermophysical_properties_pointer)		:: thermo
		type(chemical_properties_pointer)			:: chem

		real(dkind) ,dimension(:)   ,allocatable    :: thermal_c_coeff_constant
	contains
		procedure	,private	::	calculate_thermal_c_coeff_constant
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
		call manager%get_cons_field_pointer_by_name(scal_ptr,vect_ptr,tens_ptr,'mixture_molar_concentration')
		constructor%mol_mix_conc%s_ptr		=> scal_ptr%s_ptr
		call manager%get_cons_field_pointer_by_name(scal_ptr,vect_ptr,tens_ptr,'density')
		constructor%rho%s_ptr				=> scal_ptr%s_ptr

		call manager%create_scalar_field(E_f_prod_heat	,'energy_production_heat_transfer'	,'E_f_prod_heat')
		constructor%E_f_prod%s_ptr		=> E_f_prod_heat
		call manager%create_scalar_field(kappa			,'thermal_conductivity'				,'kappa')
		constructor%kappa%s_ptr			=> kappa

		call manager%get_cons_field_pointer_by_name(scal_ptr,vect_ptr,tens_ptr,'specie_molar_concentration')
		constructor%Y%v_ptr				=> vect_ptr%v_ptr

		constructor%mesh%mesh_ptr	=> manager%computational_mesh_pointer%mesh_ptr
		constructor%boundary%bc_ptr => manager%boundary_conditions_pointer%bc_ptr

		constructor%domain				= manager%domain
		constructor%mpi_support			= manager%mpi_communications

		constructor%thermo%thermo_ptr	=> manager%thermophysics%thermo_ptr
		constructor%chem%chem_ptr		=> manager%chemistry%chem_ptr

		constructor%E_f_prod%s_ptr%cells(:,:,:) = 0.0_dkind
		
		allocate(constructor%thermal_c_coeff_constant(manager%chemistry%chem_ptr%species_number))

		call constructor%calculate_thermal_c_coeff_constant()
	end function

	subroutine	solve_heat_transfer(this,time_step)

		class(heat_transfer_solver) ,intent(inout) :: this
		real(dkind)				,intent(in)		:: time_step

		real(dkind)	:: div_thermo_flux, thermo_flux1, thermo_flux2
		real(dkind)	:: energy_prod_summ

		real(dkind), dimension (3,3)	:: lame_coeffs		
		
		integer						:: dimensions
		integer		,dimension(3,2)	:: cons_inner_loop
		real(dkind)	,dimension(3)	:: cell_size
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
		
		
		
		associate(  E_f_prod	=> this%E_f_prod%s_ptr	, &
					rho			=> this%rho%s_ptr		, &
					T			=> this%T%s_ptr			, &
					kappa		=> this%kappa%s_ptr		, &
					bc			=> this%boundary%bc_ptr , &
					mesh		=> this%mesh%mesh_ptr)

		call this%mpi_support%exchange_conservative_scalar_field(kappa)
		call this%mpi_support%exchange_conservative_scalar_field(T)					
				
	!$omp parallel default(none)  private(i,j,k,dim,div_thermo_flux,thermo_flux1,thermo_flux2,lame_coeffs,sign,bound_number,boundary_type_name) , &
	!$omp& firstprivate(this)	,&
	!$omp& shared(E_f_prod,rho,kappa,T,time_step,cons_inner_loop,dimensions,cell_size,mesh,bc,coordinate_system)
	!$omp do collapse(3) schedule(static)
		do k = cons_inner_loop(3,1),cons_inner_loop(3,2)
		do j = cons_inner_loop(2,1),cons_inner_loop(2,2)
		do i = cons_inner_loop(1,1),cons_inner_loop(1,2)

			E_f_prod%cells(i,j,k) = 0.0_dkind
		
			if(bc%bc_markers(i,j,k) == 0) then
				
				div_thermo_flux = 0.0_dkind
				
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
				
                do dim = 1,dimensions
					thermo_flux1	= 0.5_dkind * (kappa%cells(i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3)) + kappa%cells(i,j,k)) * lame_coeffs(dim,1)  &
												* (T%cells(i,j,k) - T%cells(i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3))) / cell_size(dim)

					thermo_flux2	= 0.5_dkind * (kappa%cells(i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3)) + kappa%cells(i,j,k)) * lame_coeffs(dim,3)  &
												* (T%cells(i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3)) - T%cells(i,j,k)) / cell_size(dim)

					div_thermo_flux = div_thermo_flux  + (thermo_flux2 - thermo_flux1) / cell_size(dim) / lame_coeffs(dim,2)
					continue
                end do
				E_f_prod%cells(i,j,k)	=  div_thermo_flux !* time_step 
			
				do dim = 1,dimensions
					do plus = 1,2
						sign			= (-1)**plus
						bound_number	= bc%bc_markers(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3))
						if( bound_number /= 0 ) then
							boundary_type_name = bc%boundary_types(bound_number)%get_type_name()
							select case(boundary_type_name)
								case('inlet','outlet')
									E_f_prod%cells(i,j,k) = 0.0_dkind
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

		call this%mpi_support%exchange_conservative_scalar_field(E_f_prod)		
		
		end associate
	end subroutine


	subroutine calculate_thermal_c_coeff_constant(this)
		class(heat_transfer_solver) ,intent(inout) :: this

		real(dkind)	:: reduced_collision_diameter
		real(dkind)	:: inv_reduced_molar_mass

		integer		:: species_number
		integer		:: specie_number

		species_number = this%chem%chem_ptr%species_number

		associate(  potential_well_depth    => this%thermo%thermo_ptr%potential_well_depth              , &
					molar_masses            => this%thermo%thermo_ptr%molar_masses                      , &
					collision_diameter      => this%thermo%thermo_ptr%collision_diameter)

			do specie_number = 1,species_number
				if (molar_masses(specie_number) /= 0.0_dkind) then
					this%thermal_c_coeff_constant(specie_number)    = 0.0001_dkind  * 8.323_dkind  * sqrt(0.001_dkind / molar_masses(specie_number))/collision_diameter(specie_number)/collision_diameter(specie_number)
				end if
			end do

		end associate

	end subroutine

	subroutine calculate_thermal_c_coeff(this)
		class(heat_transfer_solver) ,intent(inout) :: this

		real(dkind)                     :: mol_frac, stc, reduced_temperature, sum1, sum2
		real(dkind)                     :: specie_cp, specie_cv
		real(dkind)                     :: omega_2_2

		integer	,dimension(3,2)	:: cons_inner_loop

		integer	:: species_number

		integer :: specie_number
		integer :: i,j,k

		species_number	= this%chem%chem_ptr%species_number

		cons_inner_loop = this%domain%get_local_inner_cells_bounds()

		associate(  T                       => this%T%s_ptr                            , &
					kappa					=> this%kappa%s_ptr							, &
					mol_mix_conc            => this%mol_mix_conc%s_ptr                 , &
					Y						=> this%Y%v_ptr									, &
					potential_well_depth    => this%thermo%thermo_ptr%potential_well_depth	, &
					molar_masses            => this%thermo%thermo_ptr%molar_masses			, &
					collision_diameter      => this%thermo%thermo_ptr%collision_diameter	, & 
					bc						=> this%boundary%bc_ptr)

	!$omp parallel default(none)  private(i,j,k,sum1,sum2,mol_frac,reduced_temperature,omega_2_2,specie_cp,specie_cv,stc,specie_number) , &
	!$omp& firstprivate(this)	,&
	!$omp& shared(T,kappa,mol_mix_conc,Y,collision_diameter,molar_masses,potential_well_depth,species_number,cons_inner_loop,bc)
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

							stc = this%thermal_c_coeff_constant(specie_number)   * sqrt(T%cells(i,j,k)) / omega_2_2 * (4.0_dkind * specie_cv / 15.0_dkind / r_gase_J + 3.0_dkind / 5.0_dkind)

							sum1 = sum1 + stc * mol_frac
							sum2 = sum2 + mol_frac / stc
						end if
					end if
				end do

				if (sum2 <= 1.0E-10_dkind) then
					kappa%cells(i,j,k)   = 0.0_dkind
				else
					kappa%cells(i,j,k)   = 0.5_dkind * (sum1 + 1.0_dkind / sum2)
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

		class(heat_transfer_solver)		,intent(inout)		:: this

		integer					:: dimensions
		integer	,dimension(3,2)	:: cons_inner_loop
		character(len=20)		:: boundary_type_name
		real(dkind)				:: wall_conductivity_ratio

		integer	:: sign, bound_number
		integer :: i,j,k,plus,dim

		dimensions		= this%domain%get_domain_dimensions()

		cons_inner_loop = this%domain%get_local_inner_cells_bounds()

		associate(  T			=> this%T%s_ptr			, &
					kappa		=> this%kappa%s_ptr		, &
					bc			=> this%boundary%bc_ptr	, &
					mesh		=> this%mesh%mesh_ptr)

		!$omp parallel default(none)  private(i,j,k,plus,dim,sign,bound_number,wall_conductivity_ratio,boundary_type_name) , &
		!$omp& firstprivate(this)	,&
		!$omp& shared(kappa,bc,cons_inner_loop,dimensions)
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
