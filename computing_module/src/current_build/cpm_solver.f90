module cpm_solver_class

    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite

    use kind_parameters
    use global_data
    use data_manager_class
    use data_io_class
    use computational_domain_class
    use computational_mesh_class
    use boundary_conditions_class
    use field_pointers
    use table_approximated_real_gas_class
    use chemical_kinetics_solver_class
    use chemical_properties_class
    use viscosity_solver_class
    use fourier_heat_transfer_solver_class
    use fickean_diffusion_solver_class
    use continuous_particles_solver_class, only: continuous_particles_solver, &
        continuous_particles_solver_c
    use mpi_communications_class
    use benchmarking
    use solver_options_class

    implicit none

    private
    public :: cpm_solver, cpm_solver_c

    ! Data-manager fields owned by the merged coarse-particle solver.
    type(field_scalar_cons), target :: cpm_density_interm
    type(field_scalar_cons), target :: cpm_energy_interm
    type(field_scalar_cons), target :: cpm_pressure_work
    type(field_scalar_flow), target :: cpm_mass_flux
    type(field_vector_cons), target :: cpm_velocity_interm
    type(field_vector_cons), target :: cpm_species_interm
    type(field_vector_cons), target :: cpm_velocity_increment
    type(field_vector_flow), target :: cpm_face_velocity

    type(timer) :: cpm_timer
    type(timer) :: cpm_gas_dynamics_timer
    type(timer) :: cpm_eos_timer
    type(timer) :: cpm_chemistry_timer
    type(timer) :: cpm_diffusion_timer
    type(timer) :: cpm_heattransfer_timer
    type(timer) :: cpm_viscosity_timer

    type :: cpm_solver
        logical :: diffusion_flag = .false.
        logical :: viscosity_flag = .false.
        logical :: heat_trans_flag = .false.
        logical :: reactive_flag = .false.
        logical :: hydrodynamics_flag = .true.
        logical :: CFL_condition_flag = .false.

        real(dp) :: courant_fraction = 0.5_dp
        real(dp) :: time = 0.0_dp
        real(dp) :: time_step = 0.0_dp
        real(dp) :: initial_time_step = 0.0_dp
        real(dp), dimension(3) :: g = 0.0_dp

        integer :: additional_particles_phases_number = 0

        type(viscosity_solver) :: visc_solver
        type(heat_transfer_solver) :: heat_trans_solver
        type(diffusion_solver) :: diff_solver
        type(chemical_kinetics_solver) :: chem_kin_solver
        type(table_approximated_real_gas) :: state_eq
        type(continuous_particles_solver), dimension(:), allocatable :: particles_solver

        type(computational_domain) :: domain
        type(mpi_communications) :: mpi_support
        type(chemical_properties_pointer) :: chem
        type(computational_mesh_pointer) :: mesh
        type(boundary_conditions_pointer) :: boundary

        type(field_scalar_cons_pointer) :: rho
        type(field_scalar_cons_pointer) :: rho_int
        type(field_scalar_cons_pointer) :: T
        type(field_scalar_cons_pointer) :: p
        type(field_scalar_cons_pointer) :: v_s
        type(field_scalar_cons_pointer) :: mix_mol_mass
        type(field_scalar_cons_pointer) :: E_f
        type(field_scalar_cons_pointer) :: E_f_int
        type(field_scalar_cons_pointer) :: E_f_prod_gd
        type(field_scalar_cons_pointer) :: E_f_prod_chem
        type(field_scalar_cons_pointer) :: E_f_prod_heat
        type(field_scalar_cons_pointer) :: E_f_prod_visc
        type(field_scalar_cons_pointer) :: E_f_prod_diff

        type(field_vector_cons_pointer) :: v
        type(field_vector_cons_pointer) :: v_int
        type(field_vector_cons_pointer) :: v_prod_gd
        type(field_vector_cons_pointer) :: v_prod_visc
        type(field_vector_cons_pointer) :: Y
        type(field_vector_cons_pointer) :: Y_int
        type(field_vector_cons_pointer) :: Y_prod_diff
        type(field_vector_cons_pointer) :: Y_prod_chem

        type(field_scalar_flow_pointer) :: m_flux
        type(field_vector_flow_pointer) :: v_f

        type(field_scalar_cons_pointer), dimension(:), allocatable :: E_f_prod_particles
        type(field_vector_cons_pointer), dimension(:), allocatable :: Y_prod_particles
        type(field_vector_cons_pointer), dimension(:), allocatable :: v_prod_particles

    contains
        procedure, private :: apply_boundary_conditions_main
        procedure, private :: apply_boundary_conditions_interm_v
        procedure, private :: apply_boundary_conditions_interm_E_Y
        procedure, private :: calculate_pressure_velocity_increment
        procedure, private :: calculate_interm_v
        procedure, private :: calculate_pressure_work_increment
        procedure, private :: calculate_interm_E_Y
        procedure, private :: calculate_lagrangian_mass_flux
        procedure, private :: apply_conservative_transport
        procedure, private :: update_face_velocity
        procedure, private :: geometry_power
        procedure, private :: cell_metric
        procedure, private :: face_metric
        procedure :: solve_problem
        procedure :: calculate_time_step
        procedure :: set_CFL_coefficient
        procedure :: get_CFL_coefficient
        procedure :: get_time_step
        procedure :: get_time
    end type cpm_solver

    interface cpm_solver_c
        module procedure constructor
    end interface cpm_solver_c

contains

    type(cpm_solver) function constructor(manager, problem_data_io)
        type(data_manager), intent(inout) :: manager
        type(data_io), intent(inout) :: problem_data_io

        type(field_scalar_cons_pointer) :: scal_c_ptr
        type(field_vector_cons_pointer) :: vect_c_ptr
        type(field_tensor_cons_pointer) :: tens_c_ptr
        type(particles_phase) :: particles_params

        real(dp) :: calculation_time
        integer :: particles_phase_counter
        character(len=40) :: var_name

        constructor%diffusion_flag = manager%solver_options%get_molecular_diffusion_flag()
        constructor%viscosity_flag = manager%solver_options%get_viscosity_flag()
        constructor%heat_trans_flag = manager%solver_options%get_heat_transfer_flag()
        constructor%reactive_flag = manager%solver_options%get_chemical_reaction_flag()
        constructor%hydrodynamics_flag = manager%solver_options%get_hydrodynamics_flag()
        constructor%courant_fraction = manager%solver_options%get_CFL_condition_coefficient()
        constructor%CFL_condition_flag = manager%solver_options%get_CFL_condition_flag()
        constructor%g = manager%solver_options%get_grav_acc()
        constructor%additional_particles_phases_number = &
            manager%solver_options%get_additional_particles_phases_number()

        if (.not.ieee_is_finite(constructor%courant_fraction) .or. &
            constructor%courant_fraction <= 0.0_dp .or. constructor%courant_fraction > 1.0_dp) then
            error stop 'CPM constructor: CFL coefficient must be in (0,1]'
        end if

        constructor%domain = manager%domain
        constructor%mpi_support = manager%mpi_communications
        constructor%chem%chem_ptr => manager%chemistry%chem_ptr
        constructor%boundary%bc_ptr => manager%boundary_conditions_pointer%bc_ptr
        constructor%mesh%mesh_ptr => manager%computational_mesh_pointer%mesh_ptr

        call manager%get_cons_field_pointer_by_name(scal_c_ptr, vect_c_ptr, tens_c_ptr, 'density')
        constructor%rho%s_ptr => scal_c_ptr%s_ptr
        call manager%get_cons_field_pointer_by_name(scal_c_ptr, vect_c_ptr, tens_c_ptr, 'temperature')
        constructor%T%s_ptr => scal_c_ptr%s_ptr
        call manager%get_cons_field_pointer_by_name(scal_c_ptr, vect_c_ptr, tens_c_ptr, 'pressure')
        constructor%p%s_ptr => scal_c_ptr%s_ptr
        call manager%get_cons_field_pointer_by_name(scal_c_ptr, vect_c_ptr, tens_c_ptr, 'full_energy')
        constructor%E_f%s_ptr => scal_c_ptr%s_ptr
        call manager%get_cons_field_pointer_by_name(scal_c_ptr, vect_c_ptr, tens_c_ptr, 'velocity')
        constructor%v%v_ptr => vect_c_ptr%v_ptr
        call manager%get_cons_field_pointer_by_name(scal_c_ptr, vect_c_ptr, tens_c_ptr, &
            'specie_mass_fraction')
        constructor%Y%v_ptr => vect_c_ptr%v_ptr
        call manager%get_cons_field_pointer_by_name(scal_c_ptr, vect_c_ptr, tens_c_ptr, &
            'mixture_molar_mass')
        constructor%mix_mol_mass%s_ptr => scal_c_ptr%s_ptr

        ! The former coarse_particles object created these fields indirectly.
        ! They are now created once, before the EOS and source subsolvers are constructed.
        call manager%create_scalar_field(cpm_density_interm, 'density_interm_cpm', 'rho_int_cpm')
        constructor%rho_int%s_ptr => cpm_density_interm
        call manager%create_scalar_field(cpm_energy_interm, 'full_energy_interm', 'E_f_int')
        constructor%E_f_int%s_ptr => cpm_energy_interm
        call manager%create_scalar_field(cpm_pressure_work, &
            'energy_production_gas_dynamics', 'E_f_prod_gd')
        constructor%E_f_prod_gd%s_ptr => cpm_pressure_work
        call manager%create_scalar_field(cpm_mass_flux, 'mass_flux', 'm_flux')
        constructor%m_flux%s_ptr => cpm_mass_flux

        call manager%create_vector_field(cpm_velocity_interm, 'velocity_interm', 'v_int', 'spatial')
        constructor%v_int%v_ptr => cpm_velocity_interm
        call manager%create_vector_field(cpm_species_interm, &
            'specie_mass_fraction_interm', 'Y_int', 'chemical')
        constructor%Y_int%v_ptr => cpm_species_interm
        call manager%create_vector_field(cpm_velocity_increment, &
            'velocity_production_gas_dynamics', 'v_prod_gd', 'spatial')
        constructor%v_prod_gd%v_ptr => cpm_velocity_increment
        call manager%create_vector_field(cpm_face_velocity, 'velocity_flow', 'v_f', 'spatial')
        constructor%v_f%v_ptr => cpm_face_velocity

        constructor%state_eq = table_approximated_real_gas_c(manager)
        call manager%get_cons_field_pointer_by_name(scal_c_ptr, vect_c_ptr, tens_c_ptr, &
            'velocity_of_sound')
        constructor%v_s%s_ptr => scal_c_ptr%s_ptr

        if (constructor%viscosity_flag) then
            constructor%visc_solver = viscosity_solver_c(manager)
            call manager%get_cons_field_pointer_by_name(scal_c_ptr, vect_c_ptr, tens_c_ptr, &
                'energy_production_viscosity')
            constructor%E_f_prod_visc%s_ptr => scal_c_ptr%s_ptr
            call manager%get_cons_field_pointer_by_name(scal_c_ptr, vect_c_ptr, tens_c_ptr, &
                'velocity_production_viscosity')
            constructor%v_prod_visc%v_ptr => vect_c_ptr%v_ptr
        end if

        if (constructor%heat_trans_flag) then
            constructor%heat_trans_solver = heat_transfer_solver_c(manager)
            call manager%get_cons_field_pointer_by_name(scal_c_ptr, vect_c_ptr, tens_c_ptr, &
                'energy_production_heat_transfer')
            constructor%E_f_prod_heat%s_ptr => scal_c_ptr%s_ptr
        end if

        if (constructor%diffusion_flag) then
            constructor%diff_solver = diffusion_solver_c(manager)
            call manager%get_cons_field_pointer_by_name(scal_c_ptr, vect_c_ptr, tens_c_ptr, &
                'energy_production_diffusion')
            constructor%E_f_prod_diff%s_ptr => scal_c_ptr%s_ptr
            call manager%get_cons_field_pointer_by_name(scal_c_ptr, vect_c_ptr, tens_c_ptr, &
                'specie_production_diffusion')
            constructor%Y_prod_diff%v_ptr => vect_c_ptr%v_ptr
        end if

        if (constructor%reactive_flag) then
            constructor%chem_kin_solver = chemical_kinetics_solver_c(manager)
            call manager%get_cons_field_pointer_by_name(scal_c_ptr, vect_c_ptr, tens_c_ptr, &
                'energy_production_chemistry')
            constructor%E_f_prod_chem%s_ptr => scal_c_ptr%s_ptr
            call manager%get_cons_field_pointer_by_name(scal_c_ptr, vect_c_ptr, tens_c_ptr, &
                'specie_production_chemistry')
            constructor%Y_prod_chem%v_ptr => vect_c_ptr%v_ptr
        end if

        if (constructor%additional_particles_phases_number > 0) then
            allocate(constructor%particles_solver(constructor%additional_particles_phases_number))
            call constructor%particles_solver(1)%pre_constructor( &
                constructor%additional_particles_phases_number)
            allocate(constructor%E_f_prod_particles(constructor%additional_particles_phases_number))
            allocate(constructor%v_prod_particles(constructor%additional_particles_phases_number))
            allocate(constructor%Y_prod_particles(constructor%additional_particles_phases_number))

            do particles_phase_counter = 1, constructor%additional_particles_phases_number
                particles_params = manager%solver_options%get_particles_params(particles_phase_counter)
                constructor%particles_solver(particles_phase_counter) = continuous_particles_solver_c( &
                    manager, particles_params, particles_phase_counter)

                write(var_name, '(A,I2.2)') 'energy_production_particles', particles_phase_counter
                call manager%get_cons_field_pointer_by_name(scal_c_ptr, vect_c_ptr, tens_c_ptr, var_name)
                constructor%E_f_prod_particles(particles_phase_counter)%s_ptr => scal_c_ptr%s_ptr

                write(var_name, '(A,I2.2)') 'velocity_production_particles', particles_phase_counter
                call manager%get_cons_field_pointer_by_name(scal_c_ptr, vect_c_ptr, tens_c_ptr, var_name)
                constructor%v_prod_particles(particles_phase_counter)%v_ptr => vect_c_ptr%v_ptr

                write(var_name, '(A,I2.2)') 'concentration_production_particles', &
                    particles_phase_counter
                call manager%get_cons_field_pointer_by_name(scal_c_ptr, vect_c_ptr, tens_c_ptr, var_name)
                constructor%Y_prod_particles(particles_phase_counter)%v_ptr => vect_c_ptr%v_ptr
            end do
        end if

        problem_data_io = data_io_c(manager, calculation_time)
        if (problem_data_io%get_load_counter() /= 0) then
            call problem_data_io%add_io_scalar_cons_field(constructor%E_f)
        end if
        call problem_data_io%input_all_data()
        if (problem_data_io%get_load_counter() == 1) then
            call problem_data_io%add_io_scalar_cons_field(constructor%E_f)
        end if

        if (problem_data_io%get_load_counter() == 1) then
            call constructor%state_eq%apply_state_equation_for_initial_conditions()
            do particles_phase_counter = 1, constructor%additional_particles_phases_number
                call constructor%particles_solver(particles_phase_counter)%set_initial_distributions()
            end do
        else
            call constructor%state_eq%apply_state_equation()
            call constructor%state_eq%apply_boundary_conditions_for_initial_conditions()
        end if

        call constructor%apply_boundary_conditions_main()
        call constructor%mpi_support%exchange_conservative_scalar_field(constructor%p%s_ptr)
        call constructor%mpi_support%exchange_conservative_scalar_field(constructor%rho%s_ptr)
        call constructor%mpi_support%exchange_conservative_scalar_field(constructor%T%s_ptr)
        call constructor%mpi_support%exchange_conservative_scalar_field(constructor%E_f%s_ptr)
        call constructor%mpi_support%exchange_conservative_vector_field(constructor%Y%v_ptr)
        call constructor%mpi_support%exchange_conservative_vector_field(constructor%v%v_ptr)
        call constructor%mpi_support%exchange_boundary_conditions_markers(constructor%boundary%bc_ptr)
        call constructor%mpi_support%exchange_mesh(constructor%mesh%mesh_ptr)

        constructor%time = calculation_time
        constructor%initial_time_step = manager%solver_options%get_initial_time_step()
        constructor%time_step = constructor%initial_time_step
        constructor%rho_int%s_ptr%cells = constructor%rho%s_ptr%cells
        constructor%E_f_int%s_ptr%cells = constructor%E_f%s_ptr%cells
        call constructor%update_face_velocity()

        call manager%create_timer(cpm_timer, 'CPM solver time', 'sol_t')
        call manager%create_timer(cpm_gas_dynamics_timer, 'CPM gas dynamics time', 'gd_t')
        call manager%create_timer(cpm_eos_timer, 'CPM eos solver time', 'eos_t')
        call manager%create_timer(cpm_chemistry_timer, 'CPM chemistry solver time', 'chem_t')
        call manager%create_timer(cpm_diffusion_timer, 'CPM diffusion solver time', 'diff_t')
        call manager%create_timer(cpm_heattransfer_timer, 'CPM heattransfer solver time', 'ht_t')
        call manager%create_timer(cpm_viscosity_timer, 'CPM viscosity solver time', 'visc_t')
    end function constructor


    subroutine solve_problem(this)
        class(cpm_solver), intent(inout) :: this
        integer :: phase

        call cpm_timer%tic()
        this%time = this%time + this%time_step

        call this%apply_boundary_conditions_main()
        call this%mpi_support%exchange_conservative_scalar_field(this%p%s_ptr)
        call this%mpi_support%exchange_conservative_scalar_field(this%rho%s_ptr)
        call this%mpi_support%exchange_conservative_scalar_field(this%T%s_ptr)
        call this%mpi_support%exchange_conservative_scalar_field(this%E_f%s_ptr)
        call this%mpi_support%exchange_conservative_vector_field(this%Y%v_ptr)
        call this%mpi_support%exchange_conservative_vector_field(this%v%v_ptr)

        do phase = 1, this%additional_particles_phases_number
            call this%particles_solver(phase)%apply_boundary_conditions_main(this%time)
        end do

        call cpm_viscosity_timer%tic()
        if (this%viscosity_flag) call this%visc_solver%solve_viscosity(this%time_step)
        call cpm_viscosity_timer%toc(new_iter=.true.)

        call cpm_heattransfer_timer%tic()
        if (this%heat_trans_flag) call this%heat_trans_solver%solve_heat_transfer(this%time_step)
        call cpm_heattransfer_timer%toc(new_iter=.true.)

        call cpm_diffusion_timer%tic()
        if (this%diffusion_flag) call this%diff_solver%solve_diffusion(this%time_step)
        call cpm_diffusion_timer%toc(new_iter=.true.)

        call cpm_chemistry_timer%tic()
        if (this%reactive_flag) call this%chem_kin_solver%solve_chemical_kinetics(this%time_step)
        call cpm_chemistry_timer%toc(new_iter=.true.)

        do phase = 1, this%additional_particles_phases_number
            call this%particles_solver(phase)%particles_euler_step_v_E(this%time_step)
            call this%particles_solver(phase)%apply_boundary_conditions_interm_v_p()
            call this%particles_solver(phase)%particles_lagrange_step(this%time_step)
            call this%particles_solver(phase)%particles_final_step(this%time_step)
        end do

        call cpm_gas_dynamics_timer%tic()
        call this%calculate_pressure_velocity_increment(this%time_step)
        call this%calculate_interm_v(this%time_step)
        call this%apply_boundary_conditions_interm_v()
        call this%calculate_pressure_work_increment(this%time_step)
        call this%calculate_interm_E_Y(this%time_step)
        call this%apply_boundary_conditions_interm_E_Y()
        call this%calculate_lagrangian_mass_flux(this%time_step)
        call this%apply_conservative_transport()
        call cpm_gas_dynamics_timer%toc(new_iter=.true.)

        call cpm_eos_timer%tic()
        call this%state_eq%apply_state_equation()
        call cpm_eos_timer%toc(new_iter=.true.)

        call this%apply_boundary_conditions_main()
        call this%update_face_velocity()
        if (this%CFL_condition_flag) call this%calculate_time_step()
        call cpm_timer%toc(new_iter=.true.)
    end subroutine solve_problem


    subroutine calculate_pressure_velocity_increment(this, time_step)
        class(cpm_solver), intent(inout) :: this
        real(dp), intent(in) :: time_step

        integer :: dimensions
        integer :: i, j, k, dim
        integer, dimension(3,2) :: loop
        real(dp), dimension(3) :: dx
        real(dp) :: pressure_gradient

        dimensions = this%domain%get_domain_dimensions()
        loop = this%domain%get_local_inner_cells_bounds()
        dx = this%mesh%mesh_ptr%get_cell_edges_length()
        call this%mpi_support%exchange_conservative_scalar_field(this%p%s_ptr)

        associate(p => this%p%s_ptr, rho => this%rho%s_ptr, dv => this%v_prod_gd%v_ptr, &
                  bc => this%boundary%bc_ptr)
        !$omp parallel do collapse(3) schedule(static) private(i,j,k,dim,pressure_gradient)
        do k = loop(3,1), loop(3,2)
            do j = loop(2,1), loop(2,2)
                do i = loop(1,1), loop(1,2)
                    if (bc%bc_markers(i,j,k) /= 0) cycle
                    do dim = 1, dimensions
                        if (this%hydrodynamics_flag) then
                            pressure_gradient = 0.5_dp * ( &
                                p%cells(i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3)) - &
                                p%cells(i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3))) / dx(dim)
                            dv%pr(dim)%cells(i,j,k) = &
                                time_step * (this%g(dim) - pressure_gradient / rho%cells(i,j,k))
                        else
                            dv%pr(dim)%cells(i,j,k) = 0.0_dp
                        end if
                    end do
                end do
            end do
        end do
        !$omp end parallel do
        end associate
    end subroutine calculate_pressure_velocity_increment


    subroutine calculate_interm_v(this, time_step)
        class(cpm_solver), intent(inout) :: this
        real(dp), intent(in) :: time_step

        integer :: dimensions, species_number
        integer :: i, j, k, dim, spec, phase
        integer, dimension(3,2) :: loop
        real(dp), allocatable :: partial_density(:)
        real(dp) :: source, rho_new, momentum, negative_tolerance, scale

        dimensions = this%domain%get_domain_dimensions()
        species_number = this%chem%chem_ptr%species_number
        loop = this%domain%get_local_inner_cells_bounds()

        associate(rho => this%rho%s_ptr, rho_int => this%rho_int%s_ptr, &
                  v => this%v%v_ptr, v_int => this%v_int%v_ptr, &
                  dv => this%v_prod_gd%v_ptr, Y => this%Y%v_ptr, Y_int => this%Y_int%v_ptr, &
                  bc => this%boundary%bc_ptr)
        !$omp parallel default(shared) private(i,j,k,dim,spec,phase,partial_density,source, &
        !$omp& rho_new,momentum,negative_tolerance,scale)
        allocate(partial_density(species_number))
        !$omp do collapse(3) schedule(static)
        do k = loop(3,1), loop(3,2)
            do j = loop(2,1), loop(2,2)
                do i = loop(1,1), loop(1,2)
                    if (bc%bc_markers(i,j,k) /= 0) cycle

                    rho_new = 0.0_dp
                    do spec = 1, species_number
                        source = 0.0_dp
                        if (this%diffusion_flag) then
                            source = source + this%Y_prod_diff%v_ptr%pr(spec)%cells(i,j,k)
                        end if
                        if (this%reactive_flag) then
                            source = source + this%Y_prod_chem%v_ptr%pr(spec)%cells(i,j,k)
                        end if
                        do phase = 1, this%additional_particles_phases_number
                            source = source + &
                                this%Y_prod_particles(phase)%v_ptr%pr(spec)%cells(i,j,k)
                        end do

                        partial_density(spec) = rho%cells(i,j,k) * Y%pr(spec)%cells(i,j,k) + &
                            time_step * source
                        scale = max(abs(rho%cells(i,j,k) * Y%pr(spec)%cells(i,j,k)), &
                            abs(time_step * source), 1.0_dp)
                        negative_tolerance = 1000.0_dp * epsilon(1.0_dp) * scale
                        if (partial_density(spec) < -negative_tolerance) then
                            error stop 'CPM intermediate state: negative species partial density'
                        end if
                        partial_density(spec) = max(partial_density(spec), 0.0_dp)
                        rho_new = rho_new + partial_density(spec)
                    end do

                    if (.not.ieee_is_finite(rho_new) .or. rho_new <= tiny(1.0_dp)) then
                        error stop 'CPM intermediate state: non-positive density'
                    end if
                    rho_int%cells(i,j,k) = rho_new
                    do spec = 1, species_number
                        Y_int%pr(spec)%cells(i,j,k) = partial_density(spec) / rho_new
                    end do

                    do dim = 1, dimensions
                        momentum = rho%cells(i,j,k) * &
                            (v%pr(dim)%cells(i,j,k) + dv%pr(dim)%cells(i,j,k))
                        if (this%viscosity_flag) then
                            momentum = momentum + time_step * &
                                this%v_prod_visc%v_ptr%pr(dim)%cells(i,j,k)
                        end if
                        do phase = 1, this%additional_particles_phases_number
                            momentum = momentum + time_step * rho%cells(i,j,k) * &
                                this%v_prod_particles(phase)%v_ptr%pr(dim)%cells(i,j,k)
                        end do
                        v_int%pr(dim)%cells(i,j,k) = momentum / rho_new
                    end do
                end do
            end do
        end do
        !$omp end do
        deallocate(partial_density)
        !$omp end parallel
        end associate
    end subroutine calculate_interm_v


    subroutine calculate_pressure_work_increment(this, time_step)
        class(cpm_solver), intent(inout) :: this
        real(dp), intent(in) :: time_step

        integer :: dimensions, geom_power
        integer :: i, j, k, dim
        integer :: il, jl, kl, ir, jr, kr
        integer, dimension(3,2) :: loop
        real(dp), dimension(3) :: base_area
        real(dp) :: base_volume, volume, area_left, area_right
        real(dp) :: p_left, p_right, u_left, u_right, pressure_flux_divergence
        real(dp) :: gravity_work

        dimensions = this%domain%get_domain_dimensions()
        geom_power = this%geometry_power()
        loop = this%domain%get_local_inner_cells_bounds()
        base_area = this%mesh%mesh_ptr%get_cell_surface_area()
        base_volume = this%mesh%mesh_ptr%get_cell_volume()

        call this%mpi_support%exchange_conservative_scalar_field(this%p%s_ptr)
        call this%mpi_support%exchange_conservative_vector_field(this%v_int%v_ptr)

        associate(p => this%p%s_ptr, rho => this%rho%s_ptr, v => this%v%v_ptr, &
                  v_int => this%v_int%v_ptr, dE => this%E_f_prod_gd%s_ptr, &
                  bc => this%boundary%bc_ptr)
        !$omp parallel do collapse(3) schedule(static) private(i,j,k,dim,il,jl,kl,ir,jr,kr, &
        !$omp& volume,area_left,area_right,p_left,p_right,u_left,u_right, &
        !$omp& pressure_flux_divergence,gravity_work)
        do k = loop(3,1), loop(3,2)
            do j = loop(2,1), loop(2,2)
                do i = loop(1,1), loop(1,2)
                    if (bc%bc_markers(i,j,k) /= 0) cycle
                    if (.not.this%hydrodynamics_flag) then
                        dE%cells(i,j,k) = 0.0_dp
                        cycle
                    end if

                    volume = base_volume * this%cell_metric(i, j, k, geom_power)
                    pressure_flux_divergence = 0.0_dp
                    gravity_work = 0.0_dp

                    do dim = 1, dimensions
                        il = i - I_m(dim,1)
                        jl = j - I_m(dim,2)
                        kl = k - I_m(dim,3)
                        ir = i + I_m(dim,1)
                        jr = j + I_m(dim,2)
                        kr = k + I_m(dim,3)

                        area_left = base_area(dim) * &
                            this%face_metric(dim, i, j, k, geom_power)
                        area_right = base_area(dim) * &
                            this%face_metric(dim, ir, jr, kr, geom_power)
                        p_left = 0.5_dp * (p%cells(il,jl,kl) + p%cells(i,j,k))
                        p_right = 0.5_dp * (p%cells(i,j,k) + p%cells(ir,jr,kr))
                        u_left = 0.5_dp * (v_int%pr(dim)%cells(il,jl,kl) + &
                            v_int%pr(dim)%cells(i,j,k))
                        u_right = 0.5_dp * (v_int%pr(dim)%cells(i,j,k) + &
                            v_int%pr(dim)%cells(ir,jr,kr))
                        pressure_flux_divergence = pressure_flux_divergence + &
                            p_right * u_right * area_right - p_left * u_left * area_left
                        gravity_work = gravity_work + 0.5_dp * time_step * this%g(dim) * &
                            (v%pr(dim)%cells(i,j,k) + v_int%pr(dim)%cells(i,j,k))
                    end do

                    dE%cells(i,j,k) = -time_step * pressure_flux_divergence / &
                        (rho%cells(i,j,k) * volume) + gravity_work
                end do
            end do
        end do
        !$omp end parallel do
        end associate
    end subroutine calculate_pressure_work_increment


    subroutine calculate_interm_E_Y(this, time_step)
        class(cpm_solver), intent(inout) :: this
        real(dp), intent(in) :: time_step

        integer :: i, j, k, phase
        integer, dimension(3,2) :: loop
        real(dp) :: energy_density

        loop = this%domain%get_local_inner_cells_bounds()

        associate(rho => this%rho%s_ptr, rho_int => this%rho_int%s_ptr, &
                  E => this%E_f%s_ptr, E_int => this%E_f_int%s_ptr, &
                  dE => this%E_f_prod_gd%s_ptr, bc => this%boundary%bc_ptr)
        !$omp parallel do collapse(3) schedule(static) private(i,j,k,phase,energy_density)
        do k = loop(3,1), loop(3,2)
            do j = loop(2,1), loop(2,2)
                do i = loop(1,1), loop(1,2)
                    if (bc%bc_markers(i,j,k) /= 0) cycle
                    energy_density = rho%cells(i,j,k) * (E%cells(i,j,k) + dE%cells(i,j,k))
                    if (this%viscosity_flag) then
                        energy_density = energy_density + time_step * &
                            this%E_f_prod_visc%s_ptr%cells(i,j,k)
                    end if
                    if (this%heat_trans_flag) then
                        energy_density = energy_density + time_step * &
                            this%E_f_prod_heat%s_ptr%cells(i,j,k)
                    end if
                    if (this%diffusion_flag) then
                        energy_density = energy_density + time_step * &
                            this%E_f_prod_diff%s_ptr%cells(i,j,k)
                    end if
                    if (this%reactive_flag) then
                        energy_density = energy_density + time_step * &
                            this%E_f_prod_chem%s_ptr%cells(i,j,k)
                    end if
                    do phase = 1, this%additional_particles_phases_number
                        energy_density = energy_density + time_step * &
                            this%E_f_prod_particles(phase)%s_ptr%cells(i,j,k)
                    end do
                    E_int%cells(i,j,k) = energy_density / rho_int%cells(i,j,k)
                    if (.not.ieee_is_finite(E_int%cells(i,j,k))) then
                        error stop 'CPM intermediate state: non-finite total energy'
                    end if
                end do
            end do
        end do
        !$omp end parallel do
        end associate
    end subroutine calculate_interm_E_Y


    subroutine calculate_lagrangian_mass_flux(this, time_step)
        class(cpm_solver), intent(inout) :: this
        real(dp), intent(in) :: time_step

        integer :: dimensions, geom_power
        integer :: i, j, k, dim, il, jl, kl
        integer, dimension(3,2) :: cell_loop, face_loop
        real(dp), dimension(3) :: dx, base_area
        real(dp) :: u_face, jacobian, rho_upwind, area

        dimensions = this%domain%get_domain_dimensions()
        geom_power = this%geometry_power()
        cell_loop = this%domain%get_local_inner_cells_bounds()
        dx = this%mesh%mesh_ptr%get_cell_edges_length()
        base_area = this%mesh%mesh_ptr%get_cell_surface_area()

        call this%mpi_support%exchange_conservative_scalar_field(this%rho_int%s_ptr)
        call this%mpi_support%exchange_conservative_vector_field(this%v_int%v_ptr)

        associate(rho_int => this%rho_int%s_ptr, v_int => this%v_int%v_ptr, &
                  m_flux => this%m_flux%s_ptr, bc => this%boundary%bc_ptr)
        do dim = 1, dimensions
            face_loop = cell_loop
            face_loop(dim,2) = cell_loop(dim,2) + 1
            !$omp parallel do collapse(3) schedule(static) private(i,j,k,il,jl,kl,u_face, &
            !$omp& jacobian,rho_upwind,area)
            do k = face_loop(3,1), face_loop(3,2)
                do j = face_loop(2,1), face_loop(2,2)
                    do i = face_loop(1,1), face_loop(1,2)
                        il = i - I_m(dim,1)
                        jl = j - I_m(dim,2)
                        kl = k - I_m(dim,3)
                        if (bc%bc_markers(il,jl,kl) /= 0 .and. bc%bc_markers(i,j,k) /= 0) then
                            m_flux%cells(dim,i,j,k) = 0.0_dp
                            cycle
                        end if

                        u_face = 0.5_dp * (v_int%pr(dim)%cells(il,jl,kl) + &
                            v_int%pr(dim)%cells(i,j,k))
                        jacobian = 1.0_dp + time_step * &
                            (v_int%pr(dim)%cells(i,j,k) - v_int%pr(dim)%cells(il,jl,kl)) / dx(dim)
                        if (.not.ieee_is_finite(jacobian) .or. jacobian <= 0.0_dp) then
                            error stop 'CPM Lagrangian step: non-positive deformation Jacobian'
                        end if
                        if (u_face >= 0.0_dp) then
                            rho_upwind = rho_int%cells(il,jl,kl)
                        else
                            rho_upwind = rho_int%cells(i,j,k)
                        end if
                        area = base_area(dim) * this%face_metric(dim, i, j, k, geom_power)
                        m_flux%cells(dim,i,j,k) = &
                            rho_upwind * u_face * area * time_step / jacobian
                    end do
                end do
            end do
            !$omp end parallel do
        end do
        end associate
    end subroutine calculate_lagrangian_mass_flux


    subroutine apply_conservative_transport(this)
        class(cpm_solver), intent(inout) :: this

        integer :: dimensions, species_number, geom_power
        integer :: i, j, k, dim, spec
        integer :: il, jl, kl, ir, jr, kr
        integer, dimension(3,2) :: loop
        real(dp), allocatable :: rhoY_new(:)
        real(dp) :: base_volume, volume, mass_new, rhoE_new
        real(dp), dimension(3) :: momentum_new
        real(dp) :: m_left, m_right, q_left, q_right, species_sum, tolerance

        dimensions = this%domain%get_domain_dimensions()
        species_number = this%chem%chem_ptr%species_number
        geom_power = this%geometry_power()
        loop = this%domain%get_local_inner_cells_bounds()
        base_volume = this%mesh%mesh_ptr%get_cell_volume()

        call this%mpi_support%exchange_conservative_scalar_field(this%rho_int%s_ptr)
        call this%mpi_support%exchange_conservative_scalar_field(this%E_f_int%s_ptr)
        call this%mpi_support%exchange_conservative_vector_field(this%v_int%v_ptr)
        call this%mpi_support%exchange_conservative_vector_field(this%Y_int%v_ptr)
        call this%mpi_support%exchange_flow_scalar_field(this%m_flux%s_ptr)

        associate(rho => this%rho%s_ptr, rho_int => this%rho_int%s_ptr, &
                  E => this%E_f%s_ptr, E_int => this%E_f_int%s_ptr, &
                  v => this%v%v_ptr, v_int => this%v_int%v_ptr, &
                  Y => this%Y%v_ptr, Y_int => this%Y_int%v_ptr, &
                  m_flux => this%m_flux%s_ptr, bc => this%boundary%bc_ptr)
        !$omp parallel default(shared) private(i,j,k,dim,spec,il,jl,kl,ir,jr,kr,rhoY_new, &
        !$omp& volume,mass_new,rhoE_new,momentum_new,m_left,m_right,q_left,q_right, &
        !$omp& species_sum,tolerance)
        allocate(rhoY_new(species_number))
        !$omp do collapse(3) schedule(static)
        do k = loop(3,1), loop(3,2)
            do j = loop(2,1), loop(2,2)
                do i = loop(1,1), loop(1,2)
                    if (bc%bc_markers(i,j,k) /= 0) cycle
                    volume = base_volume * this%cell_metric(i, j, k, geom_power)
                    mass_new = rho_int%cells(i,j,k)
                    rhoE_new = rho_int%cells(i,j,k) * E_int%cells(i,j,k)
                    momentum_new = 0.0_dp
                    do dim = 1, dimensions
                        momentum_new(dim) = rho_int%cells(i,j,k) * v_int%pr(dim)%cells(i,j,k)
                    end do
                    do spec = 1, species_number
                        rhoY_new(spec) = rho_int%cells(i,j,k) * Y_int%pr(spec)%cells(i,j,k)
                    end do

                    do dim = 1, dimensions
                        il = i - I_m(dim,1)
                        jl = j - I_m(dim,2)
                        kl = k - I_m(dim,3)
                        ir = i + I_m(dim,1)
                        jr = j + I_m(dim,2)
                        kr = k + I_m(dim,3)
                        m_left = m_flux%cells(dim,i,j,k)
                        m_right = m_flux%cells(dim,ir,jr,kr)
                        mass_new = mass_new + (m_left - m_right) / volume

                        if (m_left >= 0.0_dp) then
                            q_left = E_int%cells(il,jl,kl)
                        else
                            q_left = E_int%cells(i,j,k)
                        end if
                        if (m_right >= 0.0_dp) then
                            q_right = E_int%cells(i,j,k)
                        else
                            q_right = E_int%cells(ir,jr,kr)
                        end if
                        rhoE_new = rhoE_new + (m_left * q_left - m_right * q_right) / volume

                        do spec = 1, species_number
                            if (m_left >= 0.0_dp) then
                                q_left = Y_int%pr(spec)%cells(il,jl,kl)
                            else
                                q_left = Y_int%pr(spec)%cells(i,j,k)
                            end if
                            if (m_right >= 0.0_dp) then
                                q_right = Y_int%pr(spec)%cells(i,j,k)
                            else
                                q_right = Y_int%pr(spec)%cells(ir,jr,kr)
                            end if
                            rhoY_new(spec) = rhoY_new(spec) + &
                                (m_left * q_left - m_right * q_right) / volume
                        end do

                        do spec = 1, dimensions
                            if (m_left >= 0.0_dp) then
                                q_left = v_int%pr(spec)%cells(il,jl,kl)
                            else
                                q_left = v_int%pr(spec)%cells(i,j,k)
                            end if
                            if (m_right >= 0.0_dp) then
                                q_right = v_int%pr(spec)%cells(i,j,k)
                            else
                                q_right = v_int%pr(spec)%cells(ir,jr,kr)
                            end if
                            momentum_new(spec) = momentum_new(spec) + &
                                (m_left * q_left - m_right * q_right) / volume
                        end do
                    end do

                    if (.not.ieee_is_finite(mass_new) .or. mass_new <= tiny(1.0_dp)) then
                        error stop 'CPM remap: non-positive density'
                    end if
                    tolerance = 1000.0_dp * epsilon(1.0_dp) * max(mass_new, 1.0_dp)
                    do spec = 1, species_number
                        if (rhoY_new(spec) < -tolerance) then
                            error stop 'CPM remap: negative species partial density'
                        end if
                        rhoY_new(spec) = max(rhoY_new(spec), 0.0_dp)
                    end do
                    species_sum = sum(rhoY_new)
                    if (species_sum <= tiny(1.0_dp)) then
                        error stop 'CPM remap: zero species density sum'
                    end if
                    if (abs(species_sum - mass_new) > 1.0e-9_dp * max(mass_new, 1.0_dp)) then
                        error stop 'CPM remap: mass and species conservation are inconsistent'
                    end if

                    rho%cells(i,j,k) = mass_new
                    E%cells(i,j,k) = rhoE_new / mass_new
                    do dim = 1, dimensions
                        v%pr(dim)%cells(i,j,k) = momentum_new(dim) / mass_new
                    end do
                    do spec = 1, species_number
                        Y%pr(spec)%cells(i,j,k) = rhoY_new(spec) / species_sum
                    end do
                end do
            end do
        end do
        !$omp end do
        deallocate(rhoY_new)
        !$omp end parallel
        end associate
    end subroutine apply_conservative_transport


    subroutine apply_boundary_conditions_main(this)
        class(cpm_solver), intent(inout) :: this

        integer :: dimensions, species_number
        integer :: i, j, k, dim, dim1, spec, plus, sign, bound_number
        integer :: ig, jg, kg
        integer, dimension(3,2) :: loop
        character(len=20) :: boundary_name
        real(dp) :: wall_temperature

        dimensions = this%domain%get_domain_dimensions()
        species_number = this%chem%chem_ptr%species_number
        loop = this%domain%get_local_inner_cells_bounds()

        associate(T => this%T%s_ptr, M => this%mix_mol_mass%s_ptr, p => this%p%s_ptr, &
                  rho => this%rho%s_ptr, v => this%v%v_ptr, Y => this%Y%v_ptr, &
                  bc => this%boundary%bc_ptr)
        !$omp parallel do collapse(3) schedule(static) private(i,j,k,dim,dim1,spec,plus,sign, &
        !$omp& bound_number,ig,jg,kg,boundary_name,wall_temperature)
        do k = loop(3,1), loop(3,2)
            do j = loop(2,1), loop(2,2)
                do i = loop(1,1), loop(1,2)
                    if (bc%bc_markers(i,j,k) /= 0) cycle
                    do dim = 1, dimensions
                        do plus = 1, 2
                            sign = (-1)**plus
                            ig = i + sign * I_m(dim,1)
                            jg = j + sign * I_m(dim,2)
                            kg = k + sign * I_m(dim,3)
                            bound_number = bc%bc_markers(ig,jg,kg)
                            if (bound_number == 0) cycle
                            boundary_name = bc%boundary_types(bound_number)%get_type_name()

                            select case (boundary_name)
                            case ('wall', 'symmetry_plane')
                                p%cells(ig,jg,kg) = p%cells(i,j,k)
                                M%cells(ig,jg,kg) = M%cells(i,j,k)
                                do spec = 1, species_number
                                    Y%pr(spec)%cells(ig,jg,kg) = Y%pr(spec)%cells(i,j,k)
                                end do
                                T%cells(ig,jg,kg) = T%cells(i,j,k)
                                if (boundary_name == 'wall') then
                                    if (bc%boundary_types(bound_number)%is_conductive()) then
                                        wall_temperature = &
                                            bc%boundary_types(bound_number)%get_wall_temperature()
                                        if (wall_temperature <= 0.0_dp) then
                                            error stop 'CPM wall boundary: non-positive wall temperature'
                                        end if
                                        T%cells(ig,jg,kg) = wall_temperature
                                    end if
                                end if
                                rho%cells(ig,jg,kg) = p%cells(ig,jg,kg) * M%cells(ig,jg,kg) / &
                                    (r_gase_J * T%cells(ig,jg,kg))

                                do dim1 = 1, dimensions
                                    if (boundary_name == 'wall' .and. &
                                        .not.bc%boundary_types(bound_number)%is_slip()) then
                                        v%pr(dim1)%cells(ig,jg,kg) = -v%pr(dim1)%cells(i,j,k)
                                    else if (dim1 == dim) then
                                        v%pr(dim1)%cells(ig,jg,kg) = -v%pr(dim1)%cells(i,j,k)
                                    else
                                        v%pr(dim1)%cells(ig,jg,kg) = v%pr(dim1)%cells(i,j,k)
                                    end if
                                end do

                            case ('outlet')
                                p%cells(ig,jg,kg) = p%cells(i,j,k)
                                rho%cells(ig,jg,kg) = rho%cells(i,j,k)
                                T%cells(ig,jg,kg) = T%cells(i,j,k)
                                M%cells(ig,jg,kg) = M%cells(i,j,k)
                                do spec = 1, species_number
                                    Y%pr(spec)%cells(ig,jg,kg) = Y%pr(spec)%cells(i,j,k)
                                end do
                                do dim1 = 1, dimensions
                                    v%pr(dim1)%cells(ig,jg,kg) = v%pr(dim1)%cells(i,j,k)
                                end do

                            case ('inlet')
                                ! Inlet ghost values are initialized by the boundary/EOS layer and remain fixed.
                                continue
                            case default
                                error stop 'CPM boundary: unsupported boundary type'
                            end select
                        end do
                    end do
                end do
            end do
        end do
        !$omp end parallel do
        end associate
    end subroutine apply_boundary_conditions_main


    subroutine apply_boundary_conditions_interm_v(this)
        class(cpm_solver), intent(inout) :: this

        integer :: dimensions, species_number
        integer :: i, j, k, dim, dim1, spec, plus, sign, bound_number
        integer :: ig, jg, kg
        integer, dimension(3,2) :: loop
        character(len=20) :: boundary_name

        dimensions = this%domain%get_domain_dimensions()
        species_number = this%chem%chem_ptr%species_number
        loop = this%domain%get_local_inner_cells_bounds()

        associate(rho => this%rho%s_ptr, rho_int => this%rho_int%s_ptr, &
                  v => this%v%v_ptr, v_int => this%v_int%v_ptr, &
                  Y => this%Y%v_ptr, Y_int => this%Y_int%v_ptr, bc => this%boundary%bc_ptr)
        !$omp parallel do collapse(3) schedule(static) private(i,j,k,dim,dim1,spec,plus,sign, &
        !$omp& bound_number,ig,jg,kg,boundary_name)
        do k = loop(3,1), loop(3,2)
            do j = loop(2,1), loop(2,2)
                do i = loop(1,1), loop(1,2)
                    if (bc%bc_markers(i,j,k) /= 0) cycle
                    do dim = 1, dimensions
                        do plus = 1, 2
                            sign = (-1)**plus
                            ig = i + sign * I_m(dim,1)
                            jg = j + sign * I_m(dim,2)
                            kg = k + sign * I_m(dim,3)
                            bound_number = bc%bc_markers(ig,jg,kg)
                            if (bound_number == 0) cycle
                            boundary_name = bc%boundary_types(bound_number)%get_type_name()

                            select case (boundary_name)
                            case ('wall', 'symmetry_plane')
                                rho_int%cells(ig,jg,kg) = rho_int%cells(i,j,k)
                                do spec = 1, species_number
                                    Y_int%pr(spec)%cells(ig,jg,kg) = Y_int%pr(spec)%cells(i,j,k)
                                end do
                                do dim1 = 1, dimensions
                                    if (boundary_name == 'wall' .and. &
                                        .not.bc%boundary_types(bound_number)%is_slip()) then
                                        v_int%pr(dim1)%cells(ig,jg,kg) = -v_int%pr(dim1)%cells(i,j,k)
                                    else if (dim1 == dim) then
                                        v_int%pr(dim1)%cells(ig,jg,kg) = -v_int%pr(dim1)%cells(i,j,k)
                                    else
                                        v_int%pr(dim1)%cells(ig,jg,kg) = v_int%pr(dim1)%cells(i,j,k)
                                    end if
                                end do
                            case ('outlet')
                                rho_int%cells(ig,jg,kg) = rho_int%cells(i,j,k)
                                do spec = 1, species_number
                                    Y_int%pr(spec)%cells(ig,jg,kg) = Y_int%pr(spec)%cells(i,j,k)
                                end do
                                do dim1 = 1, dimensions
                                    v_int%pr(dim1)%cells(ig,jg,kg) = v_int%pr(dim1)%cells(i,j,k)
                                end do
                            case ('inlet')
                                rho_int%cells(ig,jg,kg) = rho%cells(ig,jg,kg)
                                do spec = 1, species_number
                                    Y_int%pr(spec)%cells(ig,jg,kg) = Y%pr(spec)%cells(ig,jg,kg)
                                end do
                                do dim1 = 1, dimensions
                                    v_int%pr(dim1)%cells(ig,jg,kg) = v%pr(dim1)%cells(ig,jg,kg)
                                end do
                            case default
                                error stop 'CPM intermediate boundary: unsupported boundary type'
                            end select
                        end do
                    end do
                end do
            end do
        end do
        !$omp end parallel do
        end associate
    end subroutine apply_boundary_conditions_interm_v


    subroutine apply_boundary_conditions_interm_E_Y(this)
        class(cpm_solver), intent(inout) :: this

        integer :: dimensions, species_number
        integer :: i, j, k, dim, spec, plus, sign, bound_number
        integer :: ig, jg, kg
        integer, dimension(3,2) :: loop
        character(len=20) :: boundary_name

        dimensions = this%domain%get_domain_dimensions()
        species_number = this%chem%chem_ptr%species_number
        loop = this%domain%get_local_inner_cells_bounds()

        associate(E => this%E_f%s_ptr, E_int => this%E_f_int%s_ptr, &
                  Y => this%Y%v_ptr, Y_int => this%Y_int%v_ptr, bc => this%boundary%bc_ptr)
        !$omp parallel do collapse(3) schedule(static) private(i,j,k,dim,spec,plus,sign, &
        !$omp& bound_number,ig,jg,kg,boundary_name)
        do k = loop(3,1), loop(3,2)
            do j = loop(2,1), loop(2,2)
                do i = loop(1,1), loop(1,2)
                    if (bc%bc_markers(i,j,k) /= 0) cycle
                    do dim = 1, dimensions
                        do plus = 1, 2
                            sign = (-1)**plus
                            ig = i + sign * I_m(dim,1)
                            jg = j + sign * I_m(dim,2)
                            kg = k + sign * I_m(dim,3)
                            bound_number = bc%bc_markers(ig,jg,kg)
                            if (bound_number == 0) cycle
                            boundary_name = bc%boundary_types(bound_number)%get_type_name()
                            select case (boundary_name)
                            case ('wall', 'outlet', 'symmetry_plane')
                                E_int%cells(ig,jg,kg) = E_int%cells(i,j,k)
                                do spec = 1, species_number
                                    Y_int%pr(spec)%cells(ig,jg,kg) = Y_int%pr(spec)%cells(i,j,k)
                                end do
                            case ('inlet')
                                E_int%cells(ig,jg,kg) = E%cells(ig,jg,kg)
                                do spec = 1, species_number
                                    Y_int%pr(spec)%cells(ig,jg,kg) = Y%pr(spec)%cells(ig,jg,kg)
                                end do
                            case default
                                error stop 'CPM energy boundary: unsupported boundary type'
                            end select
                        end do
                    end do
                end do
            end do
        end do
        !$omp end parallel do
        end associate
    end subroutine apply_boundary_conditions_interm_E_Y


    subroutine update_face_velocity(this)
        class(cpm_solver), intent(inout) :: this

        integer :: dimensions
        integer :: i, j, k, dim, component, il, jl, kl
        integer, dimension(3,2) :: cell_loop, face_loop

        dimensions = this%domain%get_domain_dimensions()
        cell_loop = this%domain%get_local_inner_cells_bounds()
        call this%mpi_support%exchange_conservative_vector_field(this%v%v_ptr)

        associate(v => this%v%v_ptr, v_f => this%v_f%v_ptr, bc => this%boundary%bc_ptr)
        do dim = 1, dimensions
            face_loop = cell_loop
            face_loop(dim,2) = cell_loop(dim,2) + 1
            !$omp parallel do collapse(3) schedule(static) private(i,j,k,component,il,jl,kl)
            do k = face_loop(3,1), face_loop(3,2)
                do j = face_loop(2,1), face_loop(2,2)
                    do i = face_loop(1,1), face_loop(1,2)
                        il = i - I_m(dim,1)
                        jl = j - I_m(dim,2)
                        kl = k - I_m(dim,3)
                        if (bc%bc_markers(il,jl,kl) /= 0 .and. bc%bc_markers(i,j,k) /= 0) cycle
                        do component = 1, dimensions
                            v_f%pr(component)%cells(dim,i,j,k) = 0.5_dp * &
                                (v%pr(component)%cells(il,jl,kl) + &
                                 v%pr(component)%cells(i,j,k))
                        end do
                    end do
                end do
            end do
            !$omp end parallel do
        end do
        end associate
    end subroutine update_face_velocity


    subroutine calculate_time_step(this)
#ifdef mpi
        use MPI
#endif
        class(cpm_solver), intent(inout) :: this

        integer :: dimensions
        integer :: i, j, k, dim
#ifdef mpi
        integer :: error, communicator
#endif
        integer, dimension(3,2) :: loop
        real(dp), dimension(3) :: dx
        real(dp) :: local_limit, global_limit, acoustic_rate, compression_rate

        dimensions = this%domain%get_domain_dimensions()
        loop = this%domain%get_local_inner_cells_bounds()
        dx = this%mesh%mesh_ptr%get_cell_edges_length()
        local_limit = huge(1.0_dp)

        associate(v => this%v%v_ptr, c => this%v_s%s_ptr, bc => this%boundary%bc_ptr)
        !$omp parallel do collapse(3) schedule(static) reduction(min:local_limit) &
        !$omp& private(i,j,k,dim,acoustic_rate,compression_rate)
        do k = loop(3,1), loop(3,2)
            do j = loop(2,1), loop(2,2)
                do i = loop(1,1), loop(1,2)
                    if (bc%bc_markers(i,j,k) /= 0) cycle
                    acoustic_rate = 0.0_dp
                    do dim = 1, dimensions
                        acoustic_rate = acoustic_rate + &
                            (abs(v%pr(dim)%cells(i,j,k)) + c%cells(i,j,k)) / dx(dim)
                        compression_rate = &
                            (v%pr(dim)%cells(i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3)) - &
                             v%pr(dim)%cells(i,j,k)) / dx(dim)
                        if (compression_rate > 0.0_dp) then
                            local_limit = min(local_limit, 1.0_dp / compression_rate)
                        end if
                        compression_rate = &
                            (v%pr(dim)%cells(i,j,k) - &
                             v%pr(dim)%cells(i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3))) / dx(dim)
                        if (compression_rate > 0.0_dp) then
                            local_limit = min(local_limit, 1.0_dp / compression_rate)
                        end if
                    end do
                    if (acoustic_rate > tiny(1.0_dp)) then
                        local_limit = min(local_limit, 1.0_dp / acoustic_rate)
                    end if
                end do
            end do
        end do
        !$omp end parallel do
        end associate

        global_limit = local_limit
#ifdef mpi
        communicator = this%domain%get_mpi_communicator()
        call MPI_Allreduce(local_limit, global_limit, 1, MPI_DOUBLE_PRECISION, &
            MPI_MIN, communicator, error)
        if (error /= MPI_SUCCESS) error stop 'CPM time step: MPI_Allreduce failed'
#endif
        if (.not.ieee_is_finite(global_limit) .or. global_limit <= 0.0_dp) then
            error stop 'CPM time step: invalid stability limit'
        end if
        this%time_step = this%courant_fraction * global_limit
    end subroutine calculate_time_step




    integer function geometry_power(this)
        class(cpm_solver), intent(in) :: this
        character(len=20) :: coordinate_system

        coordinate_system = this%domain%get_coordinate_system_name()
        select case (coordinate_system)
        case ('cartesian')
            geometry_power = 0
        case ('cylindrical')
            geometry_power = 1
        case ('spherical')
            geometry_power = 2
        case default
            error stop 'CPM geometry: unsupported coordinate system'
        end select
    end function geometry_power


    real(dp) function cell_metric(this, i, j, k, power)
        class(cpm_solver), intent(in) :: this
        integer, intent(in) :: i, j, k, power
        real(dp) :: radius

        if (power == 0) then
            cell_metric = 1.0_dp
        else
            radius = max(this%mesh%mesh_ptr%mesh(1,i,j,k), 0.0_dp)
            cell_metric = radius**power
        end if
    end function cell_metric


    real(dp) function face_metric(this, dim, i, j, k, power)
        class(cpm_solver), intent(in) :: this
        integer, intent(in) :: dim, i, j, k, power
        integer :: il, jl, kl
        real(dp) :: radius

        if (power == 0) then
            face_metric = 1.0_dp
            return
        end if
        il = i - I_m(dim,1)
        jl = j - I_m(dim,2)
        kl = k - I_m(dim,3)
        radius = 0.5_dp * (this%mesh%mesh_ptr%mesh(1,il,jl,kl) + &
            this%mesh%mesh_ptr%mesh(1,i,j,k))
        face_metric = max(radius, 0.0_dp)**power
    end function face_metric


    subroutine set_CFL_coefficient(this, coefficient)
        class(cpm_solver), intent(inout) :: this
        real(dp), intent(in) :: coefficient

        if (.not.ieee_is_finite(coefficient) .or. &
            coefficient <= 0.0_dp .or. coefficient > 1.0_dp) then
            error stop 'CPM CFL coefficient must be in (0,1]'
        end if
        this%courant_fraction = coefficient
    end subroutine set_CFL_coefficient


    pure real(dp) function get_CFL_coefficient(this)
        class(cpm_solver), intent(in) :: this
        get_CFL_coefficient = this%courant_fraction
    end function get_CFL_coefficient


    pure real(dp) function get_time_step(this)
        class(cpm_solver), intent(in) :: this
        get_time_step = this%time_step
    end function get_time_step


    pure real(dp) function get_time(this)
        class(cpm_solver), intent(in) :: this
        get_time = this%time
    end function get_time

end module cpm_solver_class
