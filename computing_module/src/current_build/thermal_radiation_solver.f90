module thermal_radiation_solver_class

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
	public	:: thermal_radiation_solver, thermal_radiation_solver_c

	type(field_scalar_cons)	,target	:: E_f_prod_radiation

	type	:: thermal_radiation_solver
		type(field_scalar_cons_pointer)				:: E_f_prod, T, p, mol_mix_conc
		type(field_vector_cons_pointer)				:: Y
		type(computational_domain)					:: domain
		type(boundary_conditions_pointer)			:: boundary
		type(thermophysical_properties_pointer)		:: thermo
		type(chemical_properties_pointer)			:: chem
        
        real(dp) ,dimension(:,:)   ,allocatable     :: abso_c
        real(dp)                                    :: T_b
	contains
		procedure	,private	::	calculate_absorption_coeff
		procedure				::	solve_radiation
	end type

	interface	thermal_radiation_solver_c
		module procedure	constructor
	end interface

contains

	type(thermal_radiation_solver)	function constructor(manager)

		type(data_manager)	, intent(inout)	:: manager

		type(field_scalar_cons_pointer)	:: scal_ptr
		type(field_vector_cons_pointer)	:: vect_ptr
		type(field_tensor_cons_pointer)	:: tens_ptr
        
        integer :: H2O_index

		call manager%get_cons_field_pointer_by_name(scal_ptr,vect_ptr,tens_ptr,'temperature')
		constructor%T%s_ptr					=> scal_ptr%s_ptr
		call manager%get_cons_field_pointer_by_name(scal_ptr,vect_ptr,tens_ptr,'mixture_molar_concentration')
		constructor%mol_mix_conc%s_ptr		=> scal_ptr%s_ptr
		call manager%get_cons_field_pointer_by_name(scal_ptr,vect_ptr,tens_ptr,'pressure')
		constructor%p%s_ptr				=> scal_ptr%s_ptr

		call manager%create_scalar_field(E_f_prod_radiation	,'energy_production_radiation'	,'E_f_prod_radiation')
		constructor%E_f_prod%s_ptr		=> E_f_prod_radiation

		call manager%get_cons_field_pointer_by_name(scal_ptr,vect_ptr,tens_ptr,'specie_molar_concentration')
		constructor%Y%v_ptr				=> vect_ptr%v_ptr

        constructor%boundary%bc_ptr => manager%boundary_conditions_pointer%bc_ptr

		constructor%domain				= manager%domain

		constructor%thermo%thermo_ptr	=> manager%thermophysics%thermo_ptr
		constructor%chem%chem_ptr		=> manager%chemistry%chem_ptr

		constructor%E_f_prod%s_ptr%cells(:,:,:) = 0.0_dp
		
		allocate(constructor%abso_c(manager%chemistry%chem_ptr%species_number,0:5))
        
        constructor%abso_c = 0.0_dp
        
        H2O_index = constructor%chem%chem_ptr%get_chemical_specie_index('H2O')
        
        constructor%T_b = 300.0_dp
        
        !# Mean absorption coefficients for H2O [Barlow et.al. // Comb. Flame 2001]
        if (H2O_index /= 0) then
            constructor%abso_c(H2O_index,0) = -0.23093_dp
            constructor%abso_c(H2O_index,1) = 1.12390_dp
            constructor%abso_c(H2O_index,2) = 9.41530_dp
            constructor%abso_c(H2O_index,3) = -2.99880_dp
            constructor%abso_c(H2O_index,4) = 0.51382_dp
            constructor%abso_c(H2O_index,5) = -1.86840E-05_dp
        end if 

	end function

	subroutine	solve_radiation(this,time_step)

		class(thermal_radiation_solver) ,intent(inout)  :: this
		real(dp)				        ,intent(in)		:: time_step

        real(dp)                    :: mol_frac, sum, ap
        
		integer						:: specie        
        
		integer						:: species_number
		integer		,dimension(3,2)	:: cons_inner_loop

		integer :: i,j,k

        species_number	= this%chem%chem_ptr%species_number
        
		cons_inner_loop = this%domain%get_local_inner_cells_bounds()
		
		associate(  E_f_prod	    => this%E_f_prod%s_ptr	    , &
					p			    => this%p%s_ptr		        , &
					T			    => this%T%s_ptr			    , &
                    mol_mix_conc    => this%mol_mix_conc%s_ptr  , &
					Y				=> this%Y%v_ptr			    , &
					bc			    => this%boundary%bc_ptr     , &
                    chem            => this%chem%chem_ptr       , &
                    molar_masses    => this%thermo%thermo_ptr%molar_masses)

	!$omp parallel default(shared)  private(i,j,k,sum,mol_frac,ap,specie) !, &
	!!$omp& shared(this,time_step,cons_inner_loop,dimensions,cell_size,coordinate_system)
       
	!$omp do collapse(3) schedule(static)
		do k = cons_inner_loop(3,1),cons_inner_loop(3,2)
		do j = cons_inner_loop(2,1),cons_inner_loop(2,2)
		do i = cons_inner_loop(1,1),cons_inner_loop(1,2)

			E_f_prod%cells(i,j,k) = 0.0_dp
		    sum = 0.0_dp
            
			if(bc%bc_markers(i,j,k) == 0) then
                E_f_prod%cells(i,j,k)	=  -4.0_dp * sigma_SB * (T%cells(i,j,k)**4 - this%T_b**4) 
                do specie = 1,species_number			
                    mol_frac    = Y%pr(specie)%cells(i,j,k)/molar_masses(specie) * mol_mix_conc%cells(i,j,k)
                    ap          = this%calculate_absorption_coeff(specie,T%cells(i,j,k))        
                    sum         = sum + mol_frac * p%cells(i,j,k) * ap
                end do
			end if				
		end do
		end do
		end do
	!$omp end do nowait
	!$omp end parallel

		end associate
       
	end subroutine


	recursive pure function calculate_absorption_coeff(this,specie,temperature)
		class(thermal_radiation_solver) ,intent(in) :: this
		real(dp)                                    :: calculate_absorption_coeff
		real(dp)                        ,intent(in) :: temperature
        integer                         ,intent(in) :: specie
		!!$OMP threadprivate (calculate_mixture_cp)

		real(dp)    :: temp
		integer     :: i

        temp = temperature
        
		calculate_absorption_coeff = 0.0_dp

        do i = 0, 5
		    calculate_absorption_coeff = calculate_absorption_coeff + this%abso_c(specie,i) * (1000.0_dp / temp) ** i
        end do
	end function

    
end module
