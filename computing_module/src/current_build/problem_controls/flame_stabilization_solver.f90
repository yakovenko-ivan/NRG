module flame_stabilization_solver_class

    use kind_parameters, only: dp
    use data_manager_class
    use field_pointers
    use computational_domain_class
    use computational_mesh_class
    use boundary_conditions_class
    use chemical_properties_class
    use problem_control_class, only: flame_stabilization_control
    use supplementary_routines

    implicit none

    private
    public :: flame_stabilization_solver, flame_stabilization_solver_c


    type :: flame_stabilization_runtime_state
        real(dp), allocatable :: time_hist(:), front_coord_hist(:)
        real(dp), allocatable :: diag_time_hist(:), diag_front_coord_hist(:)
        real(dp), allocatable :: diag_inlet_velocity_hist(:)
        real(dp), allocatable :: arm_metric_hist(:)
        real(dp), allocatable :: stabilized_product_mass_fractions(:)
        character(len=200) :: data_table_filename = ''
        character(len=200) :: chem_table_filename = ''
        integer :: track_counter = 0
        integer :: correction_counter = 0
        integer :: stabilization_counter = 0
        integer :: hist_count = 0
        integer :: diag_hist_count = 0
        integer :: arm_hist_count = 0
        integer :: domain_invalid_counter = 0
        integer :: hard_recovery_corrections = 0
        integer :: hard_recovery_no_progress_count = 0
        integer :: same_sign_error_counter = 0
        integer :: post_flamelet_hold_counter = 0
        integer :: structure_stationary_counter = 0
        integer :: measurement_attempt = 0
        integer :: flame_loc_unit = -1
        integer :: physics_output_unit = -1
        real(dp) :: previous_correction_time = -huge(1.0_dp)
        real(dp) :: filtered_velocity_save = 0.0_dp
        real(dp) :: diag_filtered_velocity_save = 0.0_dp
        real(dp) :: adaptive_gain = 5.0e-02_dp
        real(dp) :: inlet_velocity_target = 0.0_dp
        real(dp) :: inlet_velocity_applied = 0.0_dp
        real(dp) :: ramp_start_time = 0.0_dp
        real(dp) :: ramp_start_velocity = 0.0_dp
        real(dp) :: active_inlet_ramp_time = 2.0e-03_dp
        real(dp) :: active_response_settle_time = 5.0e-04_dp
        real(dp) :: previous_control_velocity = 0.0_dp
        real(dp) :: previous_control_inlet_velocity = 0.0_dp
        real(dp) :: front_reference_coord = 0.0_dp
        real(dp) :: heat_release_peak_save = 0.0_dp
        real(dp) :: inlet_temperature = 0.0_dp
        real(dp) :: inlet_pressure = 0.0_dp
        real(dp) :: arm_metric_rel_trend = 1.0e30_dp
        real(dp) :: final_thermal_thickness_reference = 0.0_dp
        real(dp) :: final_reaction_thickness_reference = 0.0_dp
        real(dp) :: hard_recovery_reference_deficit = huge(1.0_dp)
        real(dp) :: hard_recovery_reference_velocity = huge(1.0_dp)
        real(dp) :: hard_recovery_last_evaluation_time = -huge(1.0_dp)
        real(dp) :: bracket_u_a = 0.0_dp
        real(dp) :: bracket_v_a = 0.0_dp
        real(dp) :: bracket_u_b = 0.0_dp
        real(dp) :: bracket_v_b = 0.0_dp
        logical :: initialized = .false.
        logical :: inlet_velocity_initialized = .false.
        logical :: output_initialized = .false.
        logical :: flamelet_output_written = .false.
        logical :: sl_output_written = .false.
        logical :: have_previous_control_point = .false.
        logical :: has_bracket = .false.
        logical :: front_reference_initialized = .false.
        logical :: controller_armed = .false.
        logical :: flame_established = .false.
        logical :: hard_recovery_waiting_response = .false.
        logical :: domain_failure = .false.
        logical :: failure_report_written = .false.
        logical :: anchor_result_written = .false.
        integer :: control_stage = 0
        real(dp) :: stabilized_inlet_velocity = 0.0_dp
        real(dp) :: measurement_inlet_velocity = 0.0_dp
        real(dp) :: measurement_start_time = 0.0_dp
        real(dp) :: measurement_start_coord = 0.0_dp
        real(dp) :: sl_displacement_save = 0.0_dp
        real(dp) :: measurement_velocity_save = 0.0_dp
        real(dp) :: measurement_r2_save = 0.0_dp
        real(dp) :: measurement_rms_save = 0.0_dp
        real(dp) :: measurement_split_slope_diff_save = 0.0_dp
        real(dp) :: current_measurement_delta = 0.0_dp
        real(dp) :: measurement_delta_sign_save = 1.0_dp
        real(dp) :: scientific_time_stabilized = 0.0_dp
        real(dp) :: scientific_u_anchor = 0.0_dp
        real(dp) :: scientific_flame_thickness = 0.0_dp
        real(dp) :: scientific_preheat_zone = 0.0_dp
        real(dp) :: scientific_reaction_zone = 0.0_dp
        real(dp) :: scientific_qint = 0.0_dp
        real(dp) :: scientific_qmax = 0.0_dp
        real(dp) :: scientific_rho_fresh = 0.0_dp
        real(dp) :: scientific_rho_products = 0.0_dp
        real(dp) :: scientific_y_h2_fresh = 0.0_dp
        real(dp) :: scientific_y_h2_products = 0.0_dp
        real(dp) :: scientific_h2_consumption_flux = 0.0_dp
        real(dp) :: scientific_consumption_speed = 0.0_dp
        real(dp) :: scientific_t_products = 0.0_dp
        real(dp) :: scientific_h2_percent = 0.0_dp
        real(dp) :: scientific_d_inlet_preheat = 0.0_dp
        real(dp) :: scientific_d_outlet_reaction = 0.0_dp
        logical :: scientific_state_captured = .false.
        logical :: scientific_consumption_speed_valid = .false.
    end type flame_stabilization_runtime_state

    type :: flame_stabilization_solver
        private
        logical :: enabled = .false.
        type(flame_stabilization_control) :: flame_stabilization

        type(computational_domain) :: domain
        type(computational_mesh_pointer) :: mesh
        type(boundary_conditions_pointer) :: boundary
        type(chemical_properties_pointer) :: chem
        type(field_scalar_cons_pointer) :: T
        type(field_scalar_cons_pointer) :: rho
        type(field_scalar_cons_pointer) :: E_f_prod_chem
        type(field_vector_cons_pointer) :: Y
        type(field_vector_cons_pointer) :: Y_prod_chem

        integer :: load_counter = 0
        real(dp) :: inlet_velocity = 0.0_dp

        logical :: flamelet_output_pending = .false.
        character(len=200) :: requested_data_table_filename = ''
        character(len=200) :: requested_chem_table_filename = ''
        type(flame_stabilization_runtime_state) :: state
    contains
        procedure :: is_enabled => flame_stabilization_solver_is_enabled
        procedure :: solve
        procedure :: get_inlet_velocity
        procedure :: consume_flamelet_output_request
        procedure, private :: set_inlet_velocity
    end type flame_stabilization_solver

    interface flame_stabilization_solver_c
        module procedure constructor
    end interface flame_stabilization_solver_c

contains

    type(flame_stabilization_solver) function constructor(manager, load_counter, initial_inlet_velocity)
        type(data_manager), intent(inout) :: manager
        integer, intent(in) :: load_counter
        real(dp), intent(in) :: initial_inlet_velocity

        type(field_scalar_cons_pointer) :: scal_ptr
        type(field_vector_cons_pointer) :: vect_ptr
        type(field_tensor_cons_pointer) :: tens_ptr
        integer :: dimensions, number_of_boundary_types
        integer :: inlet_count, outlet_count, bound_number
        character(len=20) :: boundary_type_name

        constructor%flame_stabilization = &
            manager%problem_controls_config%get_flame_stabilization()
        constructor%enabled = constructor%flame_stabilization%is_enabled()
        constructor%load_counter = load_counter
        constructor%inlet_velocity = initial_inlet_velocity

        if (.not. constructor%enabled) return

        dimensions = manager%domain%get_domain_dimensions()
        inlet_count = 0
        outlet_count = 0
        number_of_boundary_types = &
            manager%boundary_conditions_pointer%bc_ptr%get_boundary_types()
        do bound_number = 1, number_of_boundary_types
            boundary_type_name = manager%boundary_conditions_pointer%bc_ptr% &
                boundary_types(bound_number)%get_type_name()
            select case (boundary_type_name)
                case ('inlet')
                    inlet_count = inlet_count + 1
                case ('outlet')
                    outlet_count = outlet_count + 1
            end select
        end do
        call constructor%flame_stabilization%validate_problem_compatibility( &
            manager%solver_options%get_chemical_reaction_flag(), &
            dimensions, inlet_count, outlet_count)

        constructor%domain = manager%domain
        constructor%mesh%mesh_ptr => manager%computational_mesh_pointer%mesh_ptr
        constructor%boundary%bc_ptr => manager%boundary_conditions_pointer%bc_ptr
        constructor%chem%chem_ptr => manager%chemistry%chem_ptr

        call manager%get_cons_field_pointer_by_name( &
            scal_ptr, vect_ptr, tens_ptr, 'temperature')
        constructor%T%s_ptr => scal_ptr%s_ptr

        call manager%get_cons_field_pointer_by_name( &
            scal_ptr, vect_ptr, tens_ptr, 'density')
        constructor%rho%s_ptr => scal_ptr%s_ptr

        call manager%get_cons_field_pointer_by_name( &
            scal_ptr, vect_ptr, tens_ptr, 'specie_mass_fraction')
        constructor%Y%v_ptr => vect_ptr%v_ptr

        call manager%get_cons_field_pointer_by_name( &
            scal_ptr, vect_ptr, tens_ptr, 'specie_production_chemistry')
        constructor%Y_prod_chem%v_ptr => vect_ptr%v_ptr

        call manager%get_cons_field_pointer_by_name( &
            scal_ptr, vect_ptr, tens_ptr, 'energy_production_chemistry')
        constructor%E_f_prod_chem%s_ptr => scal_ptr%s_ptr
    end function constructor

    logical function flame_stabilization_solver_is_enabled(this)
        class(flame_stabilization_solver), intent(in) :: this
        flame_stabilization_solver_is_enabled = this%enabled
    end function flame_stabilization_solver_is_enabled

    real(dp) function get_inlet_velocity(this)
        class(flame_stabilization_solver), intent(in) :: this
        get_inlet_velocity = this%inlet_velocity
    end function get_inlet_velocity

    subroutine consume_flamelet_output_request(this, requested, data_table_filename, chem_table_filename)
        class(flame_stabilization_solver), intent(inout) :: this
        logical, intent(out) :: requested
        character(len=*), intent(out) :: data_table_filename
        character(len=*), intent(out) :: chem_table_filename

        requested = this%flamelet_output_pending
        data_table_filename = ''
        chem_table_filename = ''

        if (requested) then
            data_table_filename = trim(this%requested_data_table_filename)
            chem_table_filename = trim(this%requested_chem_table_filename)
            this%flamelet_output_pending = .false.
        end if
    end subroutine consume_flamelet_output_request

    subroutine set_inlet_velocity(this, inlet_velocity)
        class(flame_stabilization_solver), intent(inout) :: this
        real(dp), intent(in) :: inlet_velocity

        integer :: bound_number, boundary_types

        this%inlet_velocity = inlet_velocity
        boundary_types = this%boundary%bc_ptr%get_boundary_types()
        do bound_number = 1, boundary_types
            if (this%boundary%bc_ptr%boundary_types(bound_number)%get_type_name() == 'inlet') then
                call this%boundary%bc_ptr%boundary_types(bound_number)%set_farfield_velocity(inlet_velocity)
            end if
        end do
    end subroutine set_inlet_velocity

    !--------------------------------------------------------------------------
    ! Flame anchoring / one-dimensional laminar-burning-velocity controller.
    ! LBV retains strict 1D travelling-wave semantics. Anchor mode also supports
    ! planar 2D flames and controls only bulk translation of the heat-release
    ! centroid while cellular morphology is allowed to remain unsteady.
    ! The feedback law is inherited from the former FDS-owned implementation.
    !--------------------------------------------------------------------------
    subroutine solve(this, time, stabilized, failed)
        class(flame_stabilization_solver), intent(inout) :: this
        real(dp), intent(in) :: time
        logical, intent(out) :: stabilized
        logical, intent(out), optional :: failed

        integer, parameter :: max_hist_size = 120
        integer, parameter :: max_diag_hist_size = 200
        integer, parameter :: min_hist_for_control = 12
        integer, parameter :: min_hist_for_capture = 6
        integer, parameter :: stable_required_count = 50
        integer, parameter :: flamelet_required_count = 50
        integer, parameter :: post_flamelet_hold_required_count = 100
        integer, parameter :: final_stationarity_required_count = 100
        real(dp), parameter :: final_stationarity_relative_trend_tolerance = 2.0e-03_dp
        real(dp), parameter :: final_thickness_tolerance_cells = 1.0_dp

        real(dp), parameter :: time_control_capture = 2.0e-04_dp
        real(dp), parameter :: response_settle_capture = 2.0e-04_dp
        real(dp), parameter :: inlet_ramp_time_fine = 2.0e-03_dp
        real(dp), parameter :: inlet_ramp_time_capture = 1.0e-03_dp
        real(dp), parameter :: controller_gain_initial = 2.5e-01_dp
        real(dp), parameter :: controller_gain_min = 2.0e-02_dp
        real(dp), parameter :: controller_gain_max = 5.0e-01_dp
        real(dp), parameter :: controller_gain_capture = 5.0e-01_dp
        real(dp), parameter :: controller_max_fraction = 5.0e-02_dp
        real(dp), parameter :: controller_max_fraction_capture = 1.0e-01_dp
        real(dp), parameter :: controller_max_fraction_emergency = 2.5e-01_dp
        real(dp), parameter :: min_abs_velocity_step = 2.0e-05_dp
        real(dp), parameter :: min_abs_velocity_step_capture = 1.0e-03_dp
        real(dp), parameter :: emergency_min_fraction = 1.0e-01_dp
        real(dp), parameter :: stabilization_displacement_off_cells = 5.0e-03_dp
        real(dp), parameter :: stabilization_displacement_on_cells = 1.0e-02_dp
        integer, parameter :: arming_min_window_samples = 50
        integer, parameter :: arming_max_window_samples = 200
        real(dp), parameter :: arming_grid_periods = 2.0_dp
        real(dp), parameter :: arming_relative_trend_tolerance = 2.0e-02_dp
        real(dp), parameter :: establishment_peak_fraction = 1.0e-03_dp
        real(dp), parameter :: preheat_progress_threshold = 1.0e-02_dp
        real(dp), parameter :: thermal_progress_low = 1.0e-01_dp
        real(dp), parameter :: thermal_progress_high = 9.0e-01_dp
        real(dp), parameter :: reaction_envelope_fraction = 1.0e-02_dp
        real(dp), parameter :: preferred_margin_thicknesses = 3.0_dp
        real(dp), parameter :: hard_margin_thicknesses = 2.0_dp
        integer, parameter :: preferred_margin_cells = 10
        integer, parameter :: hard_margin_cells = 5
        integer, parameter :: inlet_probe_cells = 5
        integer, parameter :: domain_invalid_required_count = 50
        integer, parameter :: hard_recovery_max_no_progress = 10
        real(dp), parameter :: hard_recovery_response_fraction = 1.0e-01_dp
        real(dp), parameter :: hard_recovery_response_min = 2.0e-03_dp
        real(dp), parameter :: hard_recovery_response_max = 2.0e-02_dp
        real(dp), parameter :: hard_recovery_velocity_progress_fraction = 1.0e-01_dp
        real(dp), parameter :: hard_recovery_deficit_progress_cells = 2.5e-01_dp
        real(dp), parameter :: filter_alpha = 1.0e-01_dp
        real(dp), parameter :: secant_relaxation = 2.0e-01_dp
        real(dp), parameter :: feedback_sign_default = -1.0_dp
        real(dp), parameter :: heat_release_cut_fraction = 1.0e-08_dp
        real(dp), parameter :: heat_release_valid_relative = 1.0e-07_dp
        real(dp), parameter :: heat_release_valid_absolute = 1.0e-20_dp
        real(dp), parameter :: position_relaxation_time = 5.0e-02_dp
        real(dp), parameter :: position_tolerance_cells = 4.0_dp
        real(dp), parameter :: capture_position_tolerance_cells = 12.0_dp
        real(dp), parameter :: capture_velocity_threshold = 2.0e-02_dp
        real(dp), parameter :: outlet_guard_cells = 20.0_dp
        real(dp), parameter :: outlet_guard_fraction = 1.5e-01_dp
        real(dp), parameter :: min_secant_du = 5.0e-05_dp
        real(dp), parameter :: max_response_slope_abs = 1.0e+03_dp
        real(dp), parameter :: persistent_error_factor = 5.0_dp
        real(dp), parameter :: tiny_weight = tiny(1.0_dp)
        logical, parameter :: use_secant_control = .false.

        integer, parameter :: multidim_anchor_arm_samples = 20
        integer, parameter :: multidim_anchor_success_count = 50
        real(dp), parameter :: multidim_anchor_position_fraction = 5.0e-02_dp
        integer, parameter :: multidim_anchor_position_cells = 8
        real(dp), parameter :: multidim_anchor_velocity_rel_tolerance = 2.0e-02_dp
        real(dp), parameter :: multidim_anchor_velocity_abs_tolerance = 2.0e-03_dp
        real(dp), parameter :: multidim_anchor_inlet_change_fraction = 5.0e-02_dp
        real(dp), parameter :: multidim_anchor_inlet_change_abs = 5.0e-03_dp

        integer, parameter :: workflow_flamelet_sl = 1
        integer, parameter :: workflow_anchor_observation = 2
        integer :: active_workflow
        logical :: enable_flamelet_output
        logical :: enable_drift_measurement
        logical :: pause_after_sl_measurement
        integer, parameter :: stage_anchor_control = 0
        integer, parameter :: stage_flamelet_ready = 1
        integer, parameter :: stage_measurement_ramp = 2
        integer, parameter :: stage_drift_measurement = 3
        integer, parameter :: stage_measurement_done = 4
        integer, parameter :: stage_measurement_failed = 5
        integer, parameter :: min_hist_for_sl_measurement = 40
        real(dp), parameter :: measurement_delta_fraction = 1.0e-01_dp
        real(dp), parameter :: measurement_delta_max_fraction = 2.0e-01_dp
        real(dp), parameter :: measurement_delta_growth = 1.5_dp
        real(dp), parameter :: measurement_duration_safety_factor = 1.5_dp
        real(dp), parameter :: measurement_target_time = 2.0e-02_dp
        real(dp) :: measurement_max_duration
        integer, parameter :: measurement_max_attempts = 3
        real(dp), parameter :: measurement_ramp_time = 1.0e-03_dp
        real(dp), parameter :: measurement_settle_time = 1.0e-03_dp
        real(dp) :: measurement_min_displacement_cells
        real(dp), parameter :: measurement_min_duration = 5.0e-03_dp
        real(dp), parameter :: measurement_r2_min = 9.95e-01_dp
        real(dp), parameter :: measurement_split_slope_rel_tol = 2.5e-01_dp
        real(dp), parameter :: measurement_split_slope_abs_tol = 2.0e-04_dp
        real(dp), parameter :: measurement_residual_cells = 5.0e-02_dp

        integer :: dimensions, species_number, boundary_types
        integer :: front_axis
        integer :: i, j, k, dim, bound_number, specie_number, specie_index
        integer :: H2_index, H_index
        integer :: cons_inner_loop(3,2)
        integer :: active_track_number
        integer :: arming_window_samples_active
        integer :: transverse_count, transverse_index, local_valid_count

        real(dp) :: cell_size(3), cell_volume
        real(dp) :: current_flame_location(3)
        real(dp) :: heat_release_centroid(3)
        real(dp) :: H_centroid(3)
        real(dp) :: Tgrad_centroid(3)
        real(dp) :: current_front_coord
        real(dp) :: flame_velocity_lsq, flame_velocity_filtered
        real(dp) :: diag_flame_velocity_lsq, diag_flame_velocity_filtered
        real(dp) :: control_velocity, position_error, position_velocity, position_control_error
        real(dp) :: front_spread, heat_release_integral
        real(dp) :: heat_release_max, heat_release_valid_limit, H_max, Tgrad_max
        real(dp) :: temperature_max, temperature_rise, thermal_progress
        real(dp) :: x_preheat, x_thermal_low, x_thermal_high, thermal_thickness
        real(dp) :: x_reaction_left, x_reaction_right, reaction_thickness
        real(dp) :: reaction_thickness_safety, x_reaction_right_safety
        real(dp) :: outlet_reaction_tip_distance
        real(dp) :: front_x05, front_x50, front_x95
        real(dp) :: local_reaction_thickness_p50, local_reaction_thickness_p95
        real(dp) :: domain_boundary_min, domain_boundary_max
        real(dp) :: inlet_preheat_distance, outlet_reaction_distance
        real(dp) :: inlet_theta_max, inlet_margin_ratio, outlet_margin_ratio
        real(dp) :: preferred_inlet_margin, hard_inlet_margin
        real(dp) :: preferred_outlet_margin, hard_outlet_margin
        real(dp) :: inlet_margin_deficit, outlet_margin_deficit
        real(dp) :: hard_inlet_deficit, hard_outlet_deficit
        real(dp) :: preferred_shift_min, preferred_shift_max
        real(dp) :: envelope_length, hard_required_length
        real(dp) :: envelope_position_error, centroid_position_error
        real(dp) :: correction_free_time
        real(dp) :: measured_inlet_velocity, proposed_inlet_velocity
        real(dp) :: du_raw, du_limited, max_velocity_step, X_H2
        real(dp) :: time_delay, time_track, time_control, response_settle_time, inlet_ramp_time
        real(dp) :: time_control_effective, diagnostic_window_time
        real(dp) :: velocity_tolerance_on, velocity_tolerance_off
        real(dp) :: coord(3), weight, qdot_cut, qdot, hval, tgrad
        real(dp) :: sum_weight, sum_s, sum_s2
        real(dp) :: sum_H_weight, sum_Tgrad_weight
        real(dp) :: s_coord, arm_metric
        real(dp) :: ramp_elapsed, ramp_fraction
        real(dp) :: target_step_for_log
        real(dp) :: sl_displacement, linear_r2, linear_rms, split_slope_diff
        real(dp) :: measurement_elapsed, measurement_displacement, measurement_delta
        real(dp) :: measurement_attempt_duration
        real(dp) :: measurement_target_displacement, measurement_delta_max
        real(dp) :: inlet_measurement_room, outlet_measurement_room
        logical :: measurement_direction_ok
        logical :: measurement_linear_ok
        real(dp) :: hard_recovery_deficit, hard_recovery_danger_velocity
        real(dp) :: hard_recovery_response_time
        real(dp) :: front_position_tolerance, capture_position_tolerance
        real(dp) :: domain_front_min, domain_front_max, domain_front_length
        real(dp) :: outlet_guard_distance, outlet_distance
        real(dp) :: front_safe_min, front_safe_max
        real(dp) :: anchor_position_band, anchor_mean_front_position
        real(dp) :: anchor_mean_inlet_velocity, anchor_inlet_velocity_slope
        real(dp) :: anchor_inlet_velocity_rms, anchor_velocity_tolerance
        real(dp) :: anchor_inlet_window_change, anchor_inlet_change_tolerance
        real(dp), allocatable :: local_qmax(:), local_qsum(:), local_xqsum(:)
        real(dp), allocatable :: local_reaction_left(:), local_reaction_right(:)
        real(dp), allocatable :: local_samples(:), local_width_samples(:)
        real(dp), allocatable :: farfield_concentrations(:), concs(:)
        character(len=10), allocatable :: farfield_species_names(:)
        character(len=5) :: axis_names(3)
        character(len=200) :: flame_debug_file
        character(len=200) :: flame_physics_file
        character(len=1000) :: av_header
        character(len=100) :: chemical_mechanism
        logical :: found_inlet_farfield, trace_success, flame_detected, control_performed
        logical :: measurement_enabled, ramp_settled, response_settled, reset_history_after_log
        logical :: capture_mode, emergency_mode, safety_control_active
        logical :: containment_mode, hard_recovery_hold, hard_recovery_progress
        logical :: thermal_envelope_found, reaction_envelope_found
        logical :: inlet_contaminated, hard_envelope_violation
        logical :: preferred_domain_fits, hard_domain_fits, domain_warning, domain_ok
        logical :: physical_envelope_available, physical_boundary_violation
        logical :: establishment_candidate
        logical :: structure_stationary_candidate
        logical :: final_structure_stationary, anchor_velocity_interior
        logical :: multidim_anchor, anchor_quasi_stationary_candidate


        stabilized = .false.
        if (present(failed)) failed = .false.
        if (.not. this%enabled) return
        if (this%state%domain_failure) then
            if (present(failed)) failed = .true.
            return
        end if

        if (this%flame_stabilization%is_anchor()) then
            active_workflow = workflow_anchor_observation
        else if (this%flame_stabilization%is_laminar_burning_velocity()) then
            active_workflow = workflow_flamelet_sl
        else
            error stop 'Flame stabilization called with disabled/unknown problem-control mode'
        end if

        enable_flamelet_output = (active_workflow == workflow_flamelet_sl)
        enable_drift_measurement = (active_workflow == workflow_flamelet_sl)
        pause_after_sl_measurement = (active_workflow == workflow_flamelet_sl)

        dimensions = this%domain%get_domain_dimensions()
        multidim_anchor = (active_workflow == workflow_anchor_observation .and. dimensions > 1)
        axis_names = this%domain%get_axis_names()
        species_number = this%chem%chem_ptr%species_number
        boundary_types = this%boundary%bc_ptr%get_boundary_types()
        cons_inner_loop = this%domain%get_local_inner_cells_bounds()
        cell_size = this%mesh%mesh_ptr%get_cell_edges_length()

        front_axis = 1
        if (front_axis > dimensions) front_axis = 1

        cell_volume = 1.0_dp
        do dim = 1, dimensions
            cell_volume = cell_volume * cell_size(dim)
        end do
        front_position_tolerance = max(position_tolerance_cells * cell_size(front_axis), 1.0e-8_dp)
        capture_position_tolerance = max(capture_position_tolerance_cells * cell_size(front_axis), &
            front_position_tolerance)
        domain_front_min = (real(cons_inner_loop(front_axis,1),dp) - 0.5_dp) * cell_size(front_axis)
        domain_front_max = (real(cons_inner_loop(front_axis,2),dp) - 0.5_dp) * cell_size(front_axis)
        domain_boundary_min = domain_front_min - 0.5_dp * cell_size(front_axis)
        domain_boundary_max = domain_front_max + 0.5_dp * cell_size(front_axis)
        domain_front_length = max(domain_boundary_max - domain_boundary_min, cell_size(front_axis))
        outlet_guard_distance = max(outlet_guard_cells * cell_size(front_axis), &
            outlet_guard_fraction * domain_front_length)

        H2_index = this%chem%chem_ptr%get_chemical_specie_index('H2')
        H_index = this%chem%chem_ptr%get_chemical_specie_index('H')

        if (.not. allocated(this%state%time_hist)) then
            allocate(this%state%time_hist(max_hist_size), this%state%front_coord_hist(max_hist_size))
            this%state%time_hist = 0.0_dp
            this%state%front_coord_hist = 0.0_dp
        end if

        if (.not. allocated(this%state%diag_time_hist)) then
            allocate(this%state%diag_time_hist(max_diag_hist_size), &
                this%state%diag_front_coord_hist(max_diag_hist_size), &
                this%state%diag_inlet_velocity_hist(max_diag_hist_size))
            this%state%diag_time_hist = 0.0_dp
            this%state%diag_front_coord_hist = 0.0_dp
            this%state%diag_inlet_velocity_hist = 0.0_dp
        end if

        if (.not. allocated(this%state%arm_metric_hist)) then
            allocate(this%state%arm_metric_hist(arming_max_window_samples))
            this%state%arm_metric_hist = 0.0_dp
        end if

        if (.not. allocated(this%state%stabilized_product_mass_fractions)) then
            allocate(this%state%stabilized_product_mass_fractions(species_number))
            this%state%stabilized_product_mass_fractions = 0.0_dp
        end if

        if (.not. this%state%initialized) then
            allocate(concs(species_number))
            concs = 0.0_dp

            found_inlet_farfield = .false.
            do bound_number = 1, boundary_types
                if (this%boundary%bc_ptr%boundary_types(bound_number)%get_type_name() == 'inlet') then
                    call this%boundary%bc_ptr%boundary_types(bound_number)%get_farfield_concentrations(farfield_concentrations)
                    call this%boundary%bc_ptr%boundary_types(bound_number)%get_farfield_species_names(farfield_species_names)
                    this%state%inlet_temperature = &
                        this%boundary%bc_ptr%boundary_types(bound_number)%get_farfield_temperature()
                    this%state%inlet_pressure = &
                        this%boundary%bc_ptr%boundary_types(bound_number)%get_farfield_pressure()
                    found_inlet_farfield = .true.
                    exit
                end if
            end do

            if (found_inlet_farfield .and. allocated(farfield_species_names)) then
                do specie_number = 1, size(farfield_species_names)
                    specie_index = this%chem%chem_ptr%get_chemical_specie_index(farfield_species_names(specie_number))
                    if (specie_index >= 1 .and. specie_index <= species_number) then
                        concs(specie_index) = farfield_concentrations(specie_number)
                    end if
                end do
            end if

            if (sum(concs) > 0.0_dp .and. H2_index >= 1 .and. H2_index <= species_number) then
                X_H2 = concs(H2_index) / sum(concs) * 100.0_dp
            else
                X_H2 = 0.0_dp
            end if

            this%state%scientific_h2_percent = X_H2
            chemical_mechanism = trim(this%chem%chem_ptr%get_chemical_mechanism())
            this%state%data_table_filename = 'H2-Air_flamelet_' // trim(chemical_mechanism) // '_' // &
                trim(str_r(X_H2)) // '_pcnt_' // trim(str_e(cell_size(1))) // '_dx.dat'
            this%state%chem_table_filename = 'H2-Air_chem_table_' // trim(chemical_mechanism) // '_' // &
                trim(str_r(X_H2)) // '_pcnt_' // trim(str_e(cell_size(1))) // '_dx.dat'

            this%state%previous_correction_time = time
            this%state%adaptive_gain = controller_gain_initial
            this%state%initialized = .true.
        end if

        if (.not. this%state%inlet_velocity_initialized) then
            this%state%inlet_velocity_target = this%inlet_velocity
            this%state%inlet_velocity_applied = this%inlet_velocity
            this%state%ramp_start_velocity = this%state%inlet_velocity_applied
            this%state%ramp_start_time = time
            this%state%inlet_velocity_initialized = .true.
        end if

        time_delay = this%flame_stabilization%get_time_delay()
        time_track = this%flame_stabilization%get_time_track()
        time_control = this%flame_stabilization%get_time_control()
        response_settle_time = this%flame_stabilization%get_response_settle_time()
        measurement_max_duration = this%flame_stabilization%get_measurement_max_duration()
        measurement_min_displacement_cells = &
            this%flame_stabilization%get_measurement_min_displacement_cells()
        inlet_ramp_time = this%state%active_inlet_ramp_time

        if (this%state%correction_counter == 0 .and. &
            this%state%control_stage == stage_anchor_control) then
            this%state%active_response_settle_time = response_settle_time
        end if

        diagnostic_window_time = max(real(max_diag_hist_size - 1, dp) * time_track, time_track)
        velocity_tolerance_off = max(stabilization_displacement_off_cells * &
            cell_size(front_axis) / diagnostic_window_time, 1.0e-12_dp)
        velocity_tolerance_on = max(stabilization_displacement_on_cells * &
            cell_size(front_axis) / diagnostic_window_time, 2.0_dp * velocity_tolerance_off)

        if (inlet_ramp_time > 0.0_dp) then
            ramp_elapsed = max(time - this%state%ramp_start_time, 0.0_dp)
            ramp_fraction = min(ramp_elapsed / inlet_ramp_time, 1.0_dp)
        else
            ramp_fraction = 1.0_dp
        end if
        this%state%inlet_velocity_applied = this%state%ramp_start_velocity + &
            ramp_fraction * (this%state%inlet_velocity_target - this%state%ramp_start_velocity)
        call this%set_inlet_velocity(this%state%inlet_velocity_applied)

        ramp_settled = (ramp_fraction >= 1.0_dp)

        if ((time - time_delay) / time_track <= real(this%state%track_counter + 1, dp)) return

        if (multidim_anchor) then
            transverse_count = cons_inner_loop(2,2) - cons_inner_loop(2,1) + 1
            allocate(local_qmax(transverse_count), local_qsum(transverse_count), &
                local_xqsum(transverse_count), local_reaction_left(transverse_count), &
                local_reaction_right(transverse_count), local_samples(transverse_count), &
                local_width_samples(transverse_count))
        else
            transverse_count = 0
        end if

        associate (T => this%T%s_ptr, &
                   Y => this%Y%v_ptr, &
                   E_f_prod_chem => this%E_f_prod_chem%s_ptr, &
                   bc => this%boundary%bc_ptr)

            if (.not. this%state%output_initialized) call initialize_output_file()

            heat_release_max = 0.0_dp
            H_max = 0.0_dp
            Tgrad_max = 0.0_dp
            temperature_max = this%state%inlet_temperature
            if (multidim_anchor) then
                local_qmax = 0.0_dp
                local_qsum = 0.0_dp
                local_xqsum = 0.0_dp
                local_reaction_left = huge(1.0_dp)
                local_reaction_right = -huge(1.0_dp)
                local_samples = 0.0_dp
                local_width_samples = 0.0_dp
            end if
            do k = cons_inner_loop(3,1), cons_inner_loop(3,2)
            do j = cons_inner_loop(2,1), cons_inner_loop(2,2)
            do i = cons_inner_loop(1,1), cons_inner_loop(1,2)
                if (bc%bc_markers(i,j,k) /= 0) cycle

                qdot = max(E_f_prod_chem%cells(i,j,k), 0.0_dp)
                heat_release_max = max(heat_release_max, qdot)
                if (multidim_anchor) then
                    transverse_index = j - cons_inner_loop(2,1) + 1
                    local_qmax(transverse_index) = max(local_qmax(transverse_index), qdot)
                end if
                temperature_max = max(temperature_max, T%cells(i,j,k))

                if (H_index >= 1 .and. H_index <= species_number) then
                    hval = abs(Y%pr(H_index)%cells(i,j,k))
                    H_max = max(H_max, hval)
                end if

                tgrad = temperature_gradient_norm(i,j,k)
                Tgrad_max = max(Tgrad_max, tgrad)
            end do
            end do
            end do

            this%state%heat_release_peak_save = max(this%state%heat_release_peak_save, heat_release_max)
            heat_release_valid_limit = max(heat_release_valid_absolute, &
                heat_release_valid_relative * this%state%heat_release_peak_save)

            heat_release_centroid = 0.0_dp
            H_centroid = 0.0_dp
            Tgrad_centroid = 0.0_dp
            sum_weight = 0.0_dp
            sum_H_weight = 0.0_dp
            sum_Tgrad_weight = 0.0_dp
            sum_s = 0.0_dp
            sum_s2 = 0.0_dp
            heat_release_integral = 0.0_dp

            temperature_rise = max(temperature_max - this%state%inlet_temperature, 1.0e-12_dp)
            x_preheat = huge(1.0_dp)
            x_thermal_low = huge(1.0_dp)
            x_thermal_high = huge(1.0_dp)
            x_reaction_left = huge(1.0_dp)
            x_reaction_right = -huge(1.0_dp)
            inlet_theta_max = 0.0_dp

            qdot_cut = heat_release_cut_fraction * heat_release_max

            do k = cons_inner_loop(3,1), cons_inner_loop(3,2)
            do j = cons_inner_loop(2,1), cons_inner_loop(2,2)
            do i = cons_inner_loop(1,1), cons_inner_loop(1,2)
                if (bc%bc_markers(i,j,k) /= 0) cycle

                coord = cell_center_coordinates(i,j,k)
                s_coord = coord(front_axis)

                qdot = max(E_f_prod_chem%cells(i,j,k), 0.0_dp)
                thermal_progress = min(max((T%cells(i,j,k) - this%state%inlet_temperature) / &
                    temperature_rise, 0.0_dp), 1.0_dp)

                if (thermal_progress >= preheat_progress_threshold) &
                    x_preheat = min(x_preheat, s_coord)
                if (thermal_progress >= thermal_progress_low) &
                    x_thermal_low = min(x_thermal_low, s_coord)
                if (thermal_progress >= thermal_progress_high) &
                    x_thermal_high = min(x_thermal_high, s_coord)

                if (i <= cons_inner_loop(1,1) + inlet_probe_cells - 1) then
                    inlet_theta_max = max(inlet_theta_max, thermal_progress)
                end if

                if (heat_release_max > 0.0_dp .and. &
                    qdot >= reaction_envelope_fraction * heat_release_max) then
                    x_reaction_left = min(x_reaction_left, s_coord)
                    x_reaction_right = max(x_reaction_right, s_coord)
                end if

                if (multidim_anchor) then
                    transverse_index = j - cons_inner_loop(2,1) + 1
                    if (local_qmax(transverse_index) >= heat_release_valid_limit .and. &
                        qdot >= reaction_envelope_fraction * local_qmax(transverse_index)) then
                        local_reaction_left(transverse_index) = &
                            min(local_reaction_left(transverse_index), s_coord)
                        local_reaction_right(transverse_index) = &
                            max(local_reaction_right(transverse_index), s_coord)
                    end if
                    if (qdot > qdot_cut) then
                        local_qsum(transverse_index) = local_qsum(transverse_index) + qdot
                        local_xqsum(transverse_index) = &
                            local_xqsum(transverse_index) + qdot * s_coord
                    end if
                end if

                if (qdot > qdot_cut) then
                    weight = qdot * cell_volume
                    heat_release_centroid = heat_release_centroid + weight * coord
                    sum_weight = sum_weight + weight
                    sum_s = sum_s + weight * s_coord
                    sum_s2 = sum_s2 + weight * s_coord * s_coord
                    heat_release_integral = heat_release_integral + weight
                end if

                if (H_index >= 1 .and. H_index <= species_number) then
                    hval = max(Y%pr(H_index)%cells(i,j,k), 0.0_dp)
                    if (hval > 0.0_dp) then
                        weight = hval * cell_volume
                        H_centroid = H_centroid + weight * coord
                        sum_H_weight = sum_H_weight + weight
                    end if
                end if

                tgrad = temperature_gradient_norm(i,j,k)
                if (tgrad > 0.0_dp) then
                    weight = tgrad * cell_volume
                    Tgrad_centroid = Tgrad_centroid + weight * coord
                    sum_Tgrad_weight = sum_Tgrad_weight + weight
                end if
            end do
            end do
            end do

            thermal_envelope_found = (x_preheat < 0.5_dp * huge(1.0_dp)) .and. &
                (x_thermal_low < 0.5_dp * huge(1.0_dp)) .and. &
                (x_thermal_high < 0.5_dp * huge(1.0_dp))
            reaction_envelope_found = (x_reaction_left < 0.5_dp * huge(1.0_dp)) .and. &
                (x_reaction_right > -0.5_dp * huge(1.0_dp))

            front_x05 = 0.0_dp
            front_x50 = 0.0_dp
            front_x95 = 0.0_dp
            local_reaction_thickness_p50 = 0.0_dp
            local_reaction_thickness_p95 = 0.0_dp
            reaction_thickness_safety = 0.0_dp
            x_reaction_right_safety = x_reaction_right

            if (multidim_anchor) then
                local_valid_count = 0
                do transverse_index = 1, transverse_count
                    if (local_qsum(transverse_index) > tiny_weight .and. &
                        local_qmax(transverse_index) >= heat_release_valid_limit .and. &
                        local_reaction_left(transverse_index) < 0.5_dp * huge(1.0_dp) .and. &
                        local_reaction_right(transverse_index) > -0.5_dp * huge(1.0_dp)) then
                        local_valid_count = local_valid_count + 1
                        local_samples(local_valid_count) = &
                            local_xqsum(transverse_index) / local_qsum(transverse_index)
                        local_width_samples(local_valid_count) = max( &
                            local_reaction_right(transverse_index) - &
                            local_reaction_left(transverse_index), cell_size(front_axis))
                    end if
                end do

                if (local_valid_count > 0) then
                    front_x05 = percentile_value(local_samples, local_valid_count, 0.05_dp)
                    front_x50 = percentile_value(local_samples, local_valid_count, 0.50_dp)
                    front_x95 = percentile_value(local_samples, local_valid_count, 0.95_dp)
                    local_reaction_thickness_p50 = &
                        percentile_value(local_width_samples, local_valid_count, 0.50_dp)
                    local_reaction_thickness_p95 = &
                        percentile_value(local_width_samples, local_valid_count, 0.95_dp)
                    reaction_thickness_safety = local_reaction_thickness_p95

                    local_valid_count = 0
                    do transverse_index = 1, transverse_count
                        if (local_qmax(transverse_index) >= heat_release_valid_limit .and. &
                            local_reaction_right(transverse_index) > -0.5_dp * huge(1.0_dp)) then
                            local_valid_count = local_valid_count + 1
                            local_samples(local_valid_count) = &
                                local_reaction_right(transverse_index)
                        end if
                    end do
                    if (local_valid_count > 0) then
                        x_reaction_right_safety = &
                            percentile_value(local_samples, local_valid_count, 0.95_dp)
                    end if
                end if
            end if

            if (thermal_envelope_found) then
                thermal_thickness = max(x_thermal_high - x_thermal_low, cell_size(front_axis))
                inlet_preheat_distance = x_preheat - domain_boundary_min
            else
                thermal_thickness = 0.0_dp
                inlet_preheat_distance = 0.0_dp
            end if

            if (reaction_envelope_found) then
                reaction_thickness = max(x_reaction_right - x_reaction_left, cell_size(front_axis))
                outlet_reaction_tip_distance = domain_boundary_max - x_reaction_right
                if (multidim_anchor .and. reaction_thickness_safety > 0.0_dp) then
                    outlet_reaction_distance = domain_boundary_max - x_reaction_right_safety
                else
                    reaction_thickness_safety = reaction_thickness
                    outlet_reaction_distance = outlet_reaction_tip_distance
                end if
            else
                reaction_thickness = 0.0_dp
                reaction_thickness_safety = 0.0_dp
                outlet_reaction_distance = 0.0_dp
                outlet_reaction_tip_distance = 0.0_dp
            end if

            preferred_inlet_margin = max(preferred_margin_thicknesses * thermal_thickness, &
                real(preferred_margin_cells, dp) * cell_size(front_axis))
            hard_inlet_margin = max(hard_margin_thicknesses * thermal_thickness, &
                real(hard_margin_cells, dp) * cell_size(front_axis))
            preferred_outlet_margin = max(preferred_margin_thicknesses * reaction_thickness_safety, &
                real(preferred_margin_cells, dp) * cell_size(front_axis))
            hard_outlet_margin = max(hard_margin_thicknesses * reaction_thickness_safety, &
                real(hard_margin_cells, dp) * cell_size(front_axis))

            if (thermal_thickness > 0.0_dp) then
                inlet_margin_ratio = inlet_preheat_distance / thermal_thickness
            else
                inlet_margin_ratio = 0.0_dp
            end if
            if (reaction_thickness_safety > 0.0_dp) then
                outlet_margin_ratio = outlet_reaction_distance / reaction_thickness_safety
            else
                outlet_margin_ratio = 0.0_dp
            end if

            inlet_contaminated = inlet_theta_max >= preheat_progress_threshold
            inlet_margin_deficit = max(preferred_inlet_margin - inlet_preheat_distance, 0.0_dp)
            outlet_margin_deficit = max(preferred_outlet_margin - outlet_reaction_distance, 0.0_dp)
            hard_inlet_deficit = max(hard_inlet_margin - inlet_preheat_distance, 0.0_dp)
            hard_outlet_deficit = max(hard_outlet_margin - outlet_reaction_distance, 0.0_dp)

            if (thermal_envelope_found .and. reaction_envelope_found) then
                if (multidim_anchor) then
                    envelope_length = max(x_reaction_right_safety - x_preheat, 0.0_dp)
                else
                    envelope_length = max(x_reaction_right - x_preheat, 0.0_dp)
                end if
                hard_required_length = hard_inlet_margin + envelope_length + hard_outlet_margin
                hard_domain_fits = hard_required_length <= domain_front_length

                preferred_shift_min = preferred_inlet_margin - inlet_preheat_distance
                preferred_shift_max = outlet_reaction_distance - preferred_outlet_margin
                preferred_domain_fits = preferred_shift_min <= preferred_shift_max
            else
                envelope_length = 0.0_dp
                hard_required_length = huge(1.0_dp)
                preferred_shift_min = 0.0_dp
                preferred_shift_max = 0.0_dp
                preferred_domain_fits = .false.
                hard_domain_fits = .false.
            end if

            physical_envelope_available = thermal_envelope_found .and. reaction_envelope_found

            hard_envelope_violation = physical_envelope_available .and. &
                (inlet_contaminated .or. (.not. hard_domain_fits) .or. &
                 hard_inlet_deficit > 0.0_dp .or. hard_outlet_deficit > 0.0_dp)

            domain_ok = physical_envelope_available .and. &
                (.not. inlet_contaminated) .and. hard_domain_fits .and. &
                (inlet_preheat_distance >= hard_inlet_margin) .and. &
                (outlet_reaction_distance >= hard_outlet_margin)

            domain_warning = thermal_envelope_found .and. reaction_envelope_found .and. &
                preferred_domain_fits .and. &
                (preferred_shift_min > 0.0_dp .or. preferred_shift_max < 0.0_dp)

            physical_boundary_violation = physical_envelope_available .and. &
                (inlet_preheat_distance <= 0.5_dp * cell_size(front_axis) + 1.0e-12_dp .or. &
                 outlet_reaction_tip_distance <= 0.5_dp * cell_size(front_axis) + 1.0e-12_dp)

            envelope_position_error = 0.0_dp
            if (preferred_domain_fits) then
                if (preferred_shift_min > 0.0_dp) then
                    envelope_position_error = -preferred_shift_min
                else if (preferred_shift_max < 0.0_dp) then
                    envelope_position_error = -preferred_shift_max
                end if
            end if

            trace_success = .false.
            flame_detected = .false.
            current_flame_location = 0.0_dp
            front_spread = 0.0_dp

            if (sum_weight > tiny_weight) then
                current_flame_location = heat_release_centroid / sum_weight
                front_spread = max(sum_s2 / sum_weight - (sum_s / sum_weight)**2, 0.0_dp)
                front_spread = sqrt(front_spread)
                trace_success = .true.
                flame_detected = (heat_release_max >= heat_release_valid_limit)
            else if (sum_H_weight > tiny_weight) then
                current_flame_location = H_centroid / sum_H_weight
                trace_success = .true.
            else if (sum_Tgrad_weight > tiny_weight) then
                current_flame_location = Tgrad_centroid / sum_Tgrad_weight
                trace_success = .true.
            end if

            if (.not. trace_success) then
                this%state%track_counter = this%state%track_counter + 1
                return
            end if

            current_front_coord = current_flame_location(front_axis)
            if (flame_detected .and. .not. this%state%front_reference_initialized) then
                this%state%front_reference_coord = current_front_coord
                this%state%front_reference_initialized = .true.
            end if

            front_safe_min = domain_front_min + outlet_guard_distance
            front_safe_max = domain_front_max - outlet_guard_distance
            if (front_safe_min >= front_safe_max) then
                front_safe_min = domain_front_min
                front_safe_max = domain_front_max
            end if

            diag_flame_velocity_lsq = 0.0_dp
            diag_flame_velocity_filtered = 0.0_dp
            if (flame_detected) then
                call append_diagnostic_history(time, current_front_coord)
                diag_flame_velocity_lsq = diagnostic_least_squares_velocity()
                if (this%state%diag_hist_count <= 2) then
                    diag_flame_velocity_filtered = diag_flame_velocity_lsq
                    this%state%diag_filtered_velocity_save = diag_flame_velocity_filtered
                else
                    diag_flame_velocity_filtered = (1.0_dp - filter_alpha) * this%state%diag_filtered_velocity_save + &
                        filter_alpha * diag_flame_velocity_lsq
                    this%state%diag_filtered_velocity_save = diag_flame_velocity_filtered
                end if
            else
                this%state%diag_hist_count = 0
                this%state%diag_filtered_velocity_save = 0.0_dp
            end if

            ! A translating discrete flame can carry a periodic grid-phase
            ! modulation in Qint.  Average over approximately two cell-passage
            ! periods, but bound the history so slowly translating flames do not
            ! impose an excessive fixed delay.
            arming_window_samples_active = arming_max_window_samples
            if (flame_detected .and. this%state%diag_hist_count >= 2) then
                arming_window_samples_active = ceiling( &
                    arming_grid_periods * cell_size(front_axis) / &
                    max(abs(diag_flame_velocity_filtered), velocity_tolerance_on) / &
                    max(time_track, tiny_weight))
                arming_window_samples_active = max(arming_min_window_samples, &
                    min(arming_max_window_samples, arming_window_samples_active))
            end if

            if (flame_detected .and. physical_envelope_available .and. &
                heat_release_max >= establishment_peak_fraction * &
                    max(this%state%heat_release_peak_save, heat_release_valid_absolute)) then
                arm_metric = heat_release_integral
                call append_arming_history(arm_metric, arming_window_samples_active)
            else
                call clear_arming_history()
            end if

            control_performed = .false.
            reset_history_after_log = .false.
            du_raw = 0.0_dp
            du_limited = 0.0_dp
            target_step_for_log = 0.0_dp
            flame_velocity_lsq = 0.0_dp
            flame_velocity_filtered = 0.0_dp
            position_error = 0.0_dp
            position_control_error = 0.0_dp
            position_velocity = 0.0_dp
            control_velocity = 0.0_dp
            sl_displacement = this%state%sl_displacement_save
            linear_r2 = this%state%measurement_r2_save
            linear_rms = this%state%measurement_rms_save
            split_slope_diff = this%state%measurement_split_slope_diff_save
            measurement_elapsed = 0.0_dp
            measurement_displacement = 0.0_dp
            measurement_target_displacement = &
                measurement_min_displacement_cells * cell_size(front_axis)
            measurement_delta = this%state%current_measurement_delta
            measurement_delta_max = measurement_delta_max_fraction * &
                max(abs(this%state%stabilized_inlet_velocity), min_abs_velocity_step)
            measurement_direction_ok = .false.
            measurement_linear_ok = .false.
            capture_mode = .false.
            emergency_mode = .false.
            containment_mode = .false.
            hard_recovery_hold = .false.
            hard_recovery_progress = .false.
            hard_recovery_deficit = 0.0_dp
            hard_recovery_danger_velocity = 0.0_dp
            hard_recovery_response_time = hard_recovery_response_min

            if (this%state%front_reference_initialized) then
                position_error = current_front_coord - this%state%front_reference_coord
                if (multidim_anchor) then
                    anchor_position_band = max( &
                        real(multidim_anchor_position_cells, dp) * cell_size(front_axis), &
                        multidim_anchor_position_fraction * domain_front_length)
                    centroid_position_error = safe_window_error( &
                        current_front_coord, &
                        this%state%front_reference_coord - anchor_position_band, &
                        this%state%front_reference_coord + anchor_position_band)
                    position_control_error = combined_position_error( &
                        centroid_position_error, envelope_position_error)
                else if (physical_envelope_available) then
                    centroid_position_error = 0.0_dp
                    position_control_error = envelope_position_error
                else
                    centroid_position_error = safe_window_error( &
                        current_front_coord, front_safe_min, front_safe_max)
                    position_control_error = centroid_position_error
                end if
                position_velocity = position_control_error / max(position_relaxation_time, time_track)
                outlet_distance = domain_front_max - current_front_coord
                capture_mode = hard_envelope_violation
                if (.not. physical_envelope_available) then
                    capture_mode = capture_mode .or. &
                        (centroid_position_error /= 0.0_dp) .or. &
                        (outlet_distance < outlet_guard_distance)
                end if
                emergency_mode = hard_envelope_violation
                if (.not. physical_envelope_available) then
                    emergency_mode = emergency_mode .or. &
                        (outlet_distance < 0.5_dp * outlet_guard_distance) .or. &
                        (abs(centroid_position_error) > capture_position_tolerance)
                end if
            end if

            if (multidim_anchor) then
                establishment_candidate = flame_detected .and. physical_envelope_available .and. &
                    (heat_release_max >= establishment_peak_fraction * &
                        max(this%state%heat_release_peak_save, heat_release_valid_absolute)) .and. &
                    (this%state%diag_hist_count >= multidim_anchor_arm_samples)
            else
                establishment_candidate = flame_detected .and. physical_envelope_available .and. &
                    (heat_release_max >= establishment_peak_fraction * &
                        max(this%state%heat_release_peak_save, heat_release_valid_absolute)) .and. &
                    (this%state%arm_hist_count >= arming_window_samples_active) .and. &
                    (this%state%arm_metric_rel_trend <= arming_relative_trend_tolerance)
            end if

            if (.not. this%state%flame_established .and. establishment_candidate) then
                this%state%flame_established = .true.
                this%state%controller_armed = .true.
                this%state%previous_correction_time = time
                this%state%has_bracket = .false.
                this%state%have_previous_control_point = .false.
                call clear_control_history()
            end if

            ! A developing flame may be repositioned before it is scientifically
            ! established.  Preferred-envelope containment is deliberately
            ! position-only and uses the fine controller limits.
            containment_mode = (.not. this%state%flame_established) .and. &
                physical_envelope_available .and. domain_warning .and. domain_ok

            ! For a recoverable hard-margin violation, wait long enough after
            ! each correction for the broad flame structure to respond before
            ! spending another correction.
            hard_recovery_deficit = max(hard_inlet_deficit, 0.0_dp) + &
                max(hard_outlet_deficit, 0.0_dp)
            if (hard_outlet_deficit >= hard_inlet_deficit .and. &
                hard_outlet_deficit > 0.0_dp) then
                hard_recovery_danger_velocity = max(diag_flame_velocity_filtered, 0.0_dp)
            else if (hard_inlet_deficit > 0.0_dp) then
                hard_recovery_danger_velocity = max(-diag_flame_velocity_filtered, 0.0_dp)
            else
                hard_recovery_danger_velocity = 0.0_dp
            end if

            if (multidim_anchor) then
                hard_recovery_response_time = hard_recovery_response_fraction * &
                    max(reaction_thickness_safety, cell_size(front_axis)) / &
                    max(abs(diag_flame_velocity_filtered), velocity_tolerance_on)
                hard_recovery_response_time = max(hard_recovery_response_min, &
                    min(hard_recovery_response_max, hard_recovery_response_time))
            else if (thermal_envelope_found) then
                hard_recovery_response_time = hard_recovery_response_fraction * &
                    thermal_thickness / max(abs(diag_flame_velocity_filtered), &
                        velocity_tolerance_on)
                hard_recovery_response_time = max(hard_recovery_response_min, &
                    min(hard_recovery_response_max, hard_recovery_response_time))
            else
                hard_recovery_response_time = hard_recovery_response_min
            end if

            hard_recovery_hold = .false.
            if (physical_envelope_available .and. hard_domain_fits .and. &
                (.not. domain_ok) .and. this%state%hard_recovery_waiting_response) then
                if ((time - this%state%hard_recovery_last_evaluation_time) < &
                    hard_recovery_response_time) then
                    hard_recovery_hold = .true.
                else
                    hard_recovery_progress = &
                        (this%state%hard_recovery_reference_deficit - &
                            hard_recovery_deficit >= &
                            hard_recovery_deficit_progress_cells * cell_size(front_axis))
                    if (this%state%hard_recovery_reference_velocity > &
                        velocity_tolerance_on) then
                        hard_recovery_progress = hard_recovery_progress .or. &
                            (hard_recovery_danger_velocity <= &
                                (1.0_dp - hard_recovery_velocity_progress_fraction) * &
                                this%state%hard_recovery_reference_velocity)
                    end if

                    if (hard_recovery_progress) then
                        this%state%hard_recovery_no_progress_count = 0
                        this%state%hard_recovery_reference_deficit = &
                            hard_recovery_deficit
                        this%state%hard_recovery_reference_velocity = &
                            hard_recovery_danger_velocity
                        this%state%hard_recovery_last_evaluation_time = time
                        hard_recovery_hold = .true.
                    else
                        this%state%hard_recovery_no_progress_count = &
                            this%state%hard_recovery_no_progress_count + 1
                        this%state%hard_recovery_waiting_response = .false.
                    end if
                end if
            end if

            ! Safety intervention is transient and does not itself declare the
            ! ignition/development transient complete.
            safety_control_active = capture_mode .or. emergency_mode .or. containment_mode

            if (physical_boundary_violation) then
                this%state%domain_failure = .true.
                this%state%domain_invalid_counter = 0
                call reset_hard_recovery_state()
            else if (physical_envelope_available .and. (.not. hard_domain_fits)) then
                ! During ignition/development the envelope can transiently be
                ! broader than the final flame.  Do not classify geometric
                ! impossibility until the flame has actually been established.
                call reset_hard_recovery_state()
                if (this%state%flame_established) then
                    this%state%domain_invalid_counter = &
                        this%state%domain_invalid_counter + 1
                    if (this%state%domain_invalid_counter >= &
                        domain_invalid_required_count) then
                        this%state%domain_failure = .true.
                    end if
                else
                    this%state%domain_invalid_counter = 0
                end if
            else if (physical_envelope_available .and. (.not. domain_ok)) then
                ! Recoverable hard-margin violation.  Failure is based on
                ! repeated corrections that fail to produce a measurable
                ! response, not on the total number of useful corrections.
                this%state%domain_invalid_counter = 0
                if (this%state%hard_recovery_no_progress_count >= &
                    hard_recovery_max_no_progress) then
                    this%state%domain_failure = .true.
                end if
            else
                this%state%domain_invalid_counter = 0
                call reset_hard_recovery_state()
            end if

            if (capture_mode) then
                time_control_effective = time_control_capture
            else
                time_control_effective = time_control
            end if

            response_settled = ramp_settled .and. &
                (time >= this%state%ramp_start_time + inlet_ramp_time + &
                    this%state%active_response_settle_time)

            measurement_enabled = (this%state%controller_armed .or. safety_control_active) .and. &
                flame_detected .and. this%state%front_reference_initialized .and. &
                response_settled .and. &
                ((this%state%control_stage == stage_anchor_control) .or. &
                 (this%state%control_stage == stage_flamelet_ready))

            if (measurement_enabled) then
                call append_front_history(time, current_front_coord)
                flame_velocity_lsq = least_squares_velocity()

                if (this%state%hist_count <= 2) then
                    flame_velocity_filtered = flame_velocity_lsq
                    this%state%filtered_velocity_save = flame_velocity_filtered
                else
                    flame_velocity_filtered = (1.0_dp - filter_alpha) * this%state%filtered_velocity_save + &
                        filter_alpha * flame_velocity_lsq
                    this%state%filtered_velocity_save = flame_velocity_filtered
                end if

                position_error = current_front_coord - this%state%front_reference_coord
                if (multidim_anchor) then
                    anchor_position_band = max( &
                        real(multidim_anchor_position_cells, dp) * cell_size(front_axis), &
                        multidim_anchor_position_fraction * domain_front_length)
                    centroid_position_error = safe_window_error( &
                        current_front_coord, &
                        this%state%front_reference_coord - anchor_position_band, &
                        this%state%front_reference_coord + anchor_position_band)
                    position_control_error = combined_position_error( &
                        centroid_position_error, envelope_position_error)
                else if (physical_envelope_available) then
                    centroid_position_error = 0.0_dp
                    position_control_error = envelope_position_error
                else
                    centroid_position_error = safe_window_error( &
                        current_front_coord, front_safe_min, front_safe_max)
                    position_control_error = centroid_position_error
                end if
                position_velocity = position_control_error / max(position_relaxation_time, time_track)
                if (containment_mode) then
                    control_velocity = position_velocity
                else
                    control_velocity = flame_velocity_filtered + position_velocity
                end if
                outlet_distance = domain_front_max - current_front_coord
                capture_mode = hard_envelope_violation .or. &
                    (this%state%flame_established .and. &
                     abs(flame_velocity_filtered) > capture_velocity_threshold)
                if (.not. physical_envelope_available) then
                    capture_mode = capture_mode .or. &
                        (centroid_position_error /= 0.0_dp) .or. &
                        (outlet_distance < outlet_guard_distance)
                end if
                emergency_mode = hard_envelope_violation
                if (.not. physical_envelope_available) then
                    emergency_mode = emergency_mode .or. &
                        (outlet_distance < 0.5_dp * outlet_guard_distance) .or. &
                        (abs(centroid_position_error) > capture_position_tolerance)
                end if
                containment_mode = (.not. this%state%flame_established) .and. &
                    physical_envelope_available .and. domain_warning .and. domain_ok
                if (capture_mode) then
                    time_control_effective = time_control_capture
                else
                    time_control_effective = time_control
                end if
            else if (this%state%control_stage == stage_anchor_control) then
                call clear_control_history()
                if (.not. flame_detected) this%state%has_bracket = .false.
            else if (.not. flame_detected) then
                call clear_control_history()
                this%state%has_bracket = .false.
            end if

            measurement_attempt_duration = measurement_max_duration
            if (abs(this%state%current_measurement_delta) > tiny_weight) then
                measurement_attempt_duration = measurement_duration_safety_factor * &
                    measurement_min_displacement_cells * cell_size(front_axis) / &
                    abs(this%state%current_measurement_delta)
                measurement_attempt_duration = max(measurement_min_duration, &
                    min(measurement_max_duration, measurement_attempt_duration))
            end if

            if (this%state%control_stage == stage_measurement_ramp) then
                if (flame_detected .and. ramp_settled .and. &
                    time >= this%state%ramp_start_time + inlet_ramp_time + measurement_settle_time) then
                    call clear_control_history()
                    this%state%measurement_start_time = time
                    this%state%measurement_start_coord = current_front_coord
                    this%state%control_stage = stage_drift_measurement
                end if
            else if (this%state%control_stage == stage_drift_measurement) then
                if (flame_detected .and. ramp_settled) then
                    call append_front_history(time, current_front_coord)
                    flame_velocity_lsq = least_squares_velocity()
                    flame_velocity_filtered = flame_velocity_lsq
                    this%state%measurement_velocity_save = flame_velocity_lsq
                    control_velocity = flame_velocity_lsq
                    measurement_elapsed = time - this%state%measurement_start_time
                    measurement_displacement = current_front_coord - this%state%measurement_start_coord
                    call drift_linearity_diagnostics(linear_r2, linear_rms, split_slope_diff)
                    this%state%measurement_r2_save = linear_r2
                    this%state%measurement_rms_save = linear_rms
                    this%state%measurement_split_slope_diff_save = split_slope_diff
                    sl_displacement = this%state%inlet_velocity_target - flame_velocity_lsq
                    this%state%sl_displacement_save = sl_displacement
                    measurement_direction_ok = &
                        flame_velocity_lsq * this%state%measurement_delta_sign_save > 0.0_dp
                    measurement_linear_ok = domain_ok .and. &
                        measurement_direction_ok .and. &
                        (this%state%hist_count >= min_hist_for_sl_measurement) .and. &
                        (measurement_elapsed >= measurement_min_duration) .and. &
                        (abs(measurement_displacement) >= measurement_target_displacement) .and. &
                        (linear_r2 >= measurement_r2_min) .and. &
                        (linear_rms <= measurement_residual_cells * cell_size(front_axis)) .and. &
                        (split_slope_diff <= max(measurement_split_slope_abs_tol, &
                            measurement_split_slope_rel_tol * &
                            max(abs(flame_velocity_lsq), velocity_tolerance_on)))

                    if (measurement_linear_ok) then
                        this%state%measurement_velocity_save = flame_velocity_lsq
                        this%state%control_stage = stage_measurement_done
                        this%state%stabilization_counter = stable_required_count
                        call write_laminar_velocity_once()

                    else if (measurement_elapsed >= measurement_attempt_duration) then
                        measurement_delta_max = measurement_delta_max_fraction * &
                            max(abs(this%state%stabilized_inlet_velocity), min_abs_velocity_step)

                        if (abs(measurement_displacement) < measurement_target_displacement .and. &
                            this%state%measurement_attempt < measurement_max_attempts .and. &
                            this%state%current_measurement_delta < &
                                (1.0_dp - 1.0e-06_dp) * measurement_delta_max) then

                            this%state%measurement_attempt = this%state%measurement_attempt + 1
                            this%state%current_measurement_delta = min( &
                                measurement_delta_growth * this%state%current_measurement_delta, &
                                measurement_delta_max)

                            this%state%measurement_inlet_velocity = max( &
                                this%state%stabilized_inlet_velocity + &
                                this%state%measurement_delta_sign_save * &
                                    this%state%current_measurement_delta, 0.0_dp)
                            this%state%inlet_velocity_target = &
                                this%state%measurement_inlet_velocity
                            this%state%ramp_start_velocity = this%state%inlet_velocity_applied
                            this%state%ramp_start_time = time
                            this%state%active_inlet_ramp_time = measurement_ramp_time
                            this%state%control_stage = stage_measurement_ramp
                            this%state%sl_displacement_save = 0.0_dp
                            this%state%measurement_velocity_save = 0.0_dp
                            this%state%measurement_r2_save = 0.0_dp
                            this%state%measurement_rms_save = 0.0_dp
                            this%state%measurement_split_slope_diff_save = 0.0_dp
                            call clear_control_history()
                        else
                            this%state%control_stage = stage_measurement_failed
                        end if
                    end if
                else
                    call clear_control_history()
                end if
            end if

            if (measurement_enabled .and. (.not. hard_recovery_hold) .and. &
                this%state%hist_count >= merge(min_hist_for_capture, min_hist_for_control, capture_mode)) then
                if ((time - this%state%previous_correction_time) >= time_control_effective) then
                    if (control_action_needed()) then
                        measured_inlet_velocity = this%state%inlet_velocity_target

                        call update_adaptive_gain(control_velocity)
                        call update_bracket(measured_inlet_velocity, control_velocity)
                        call choose_new_inlet_target(measured_inlet_velocity, control_velocity, &
                            proposed_inlet_velocity, du_raw, du_limited)

                        if (abs(proposed_inlet_velocity - measured_inlet_velocity) > 0.0_dp) then
                            this%state%inlet_velocity_target = proposed_inlet_velocity
                            target_step_for_log = this%state%inlet_velocity_target - measured_inlet_velocity
                            if (capture_mode) then
                                this%state%active_inlet_ramp_time = inlet_ramp_time_capture
                                this%state%active_response_settle_time = response_settle_capture
                            else
                                this%state%active_inlet_ramp_time = inlet_ramp_time_fine
                                this%state%active_response_settle_time = response_settle_time
                            end if
                            this%state%ramp_start_velocity = this%state%inlet_velocity_applied
                            this%state%ramp_start_time = time
                            this%state%previous_correction_time = time
                            this%state%correction_counter = this%state%correction_counter + 1
                            if (hard_envelope_violation .and. hard_domain_fits) then
                                this%state%hard_recovery_corrections = &
                                    this%state%hard_recovery_corrections + 1
                                this%state%hard_recovery_waiting_response = .true.
                                this%state%hard_recovery_reference_deficit = &
                                    hard_recovery_deficit
                                this%state%hard_recovery_reference_velocity = &
                                    hard_recovery_danger_velocity
                                this%state%hard_recovery_last_evaluation_time = time
                            end if
                            control_performed = .true.
                            reset_history_after_log = .true.
                        end if

                        this%state%previous_control_inlet_velocity = measured_inlet_velocity
                        this%state%previous_control_velocity = control_velocity
                        this%state%have_previous_control_point = .true.
                    end if
                end if
            end if

            correction_free_time = max(time - this%state%previous_correction_time, 0.0_dp)

            ! Scientific convergence is stricter than flame establishment.
            ! A slowly evolving lean flame must not advance merely because Vfl
            ! passes transiently through zero.
            structure_stationary_candidate = this%state%flame_established .and. &
                flame_detected .and. physical_envelope_available .and. domain_ok .and. &
                (.not. domain_warning) .and. response_settled .and. &
                (.not. control_performed) .and. &
                (this%state%arm_metric_rel_trend <= &
                    final_stationarity_relative_trend_tolerance)

            if (structure_stationary_candidate) then
                if (this%state%structure_stationary_counter == 0) then
                    this%state%final_thermal_thickness_reference = thermal_thickness
                    this%state%final_reaction_thickness_reference = reaction_thickness
                    this%state%structure_stationary_counter = 1
                else if (abs(thermal_thickness - &
                    this%state%final_thermal_thickness_reference) <= &
                        final_thickness_tolerance_cells * cell_size(front_axis) .and. &
                    abs(reaction_thickness - &
                    this%state%final_reaction_thickness_reference) <= &
                        final_thickness_tolerance_cells * cell_size(front_axis)) then
                    this%state%structure_stationary_counter = &
                        this%state%structure_stationary_counter + 1
                else
                    ! Start a new structural-stationarity interval from the
                    ! current envelope instead of accumulating slow drift.
                    this%state%final_thermal_thickness_reference = thermal_thickness
                    this%state%final_reaction_thickness_reference = reaction_thickness
                    this%state%structure_stationary_counter = 1
                end if
            else
                this%state%structure_stationary_counter = 0
                this%state%final_thermal_thickness_reference = thermal_thickness
                this%state%final_reaction_thickness_reference = reaction_thickness
            end if

            final_structure_stationary = &
                (this%state%structure_stationary_counter >= &
                    final_stationarity_required_count)

            ! U=0 is the lower actuator limit, not a scientifically validated
            ! anchor.  If Vfl merely crosses zero there, continue observing;
            ! a subsequently negative Vfl will drive U positive again.
            anchor_velocity_interior = &
                (this%state%inlet_velocity_target > min_abs_velocity_step)

            if (this%state%control_stage == stage_anchor_control) then
                if (multidim_anchor) then
                    anchor_mean_front_position = diagnostic_mean_front_position()
                    anchor_mean_inlet_velocity = diagnostic_mean_inlet_velocity()
                    anchor_inlet_velocity_slope = diagnostic_inlet_velocity_slope()
                    anchor_inlet_velocity_rms = diagnostic_inlet_velocity_rms()
                    anchor_velocity_tolerance = max( &
                        multidim_anchor_velocity_abs_tolerance, &
                        multidim_anchor_velocity_rel_tolerance * &
                            max(abs(anchor_mean_inlet_velocity), min_abs_velocity_step))
                    anchor_inlet_window_change = abs(anchor_inlet_velocity_slope) * &
                        diagnostic_window_time
                    anchor_inlet_change_tolerance = max( &
                        multidim_anchor_inlet_change_abs, &
                        multidim_anchor_inlet_change_fraction * &
                            max(abs(anchor_mean_inlet_velocity), min_abs_velocity_step))
                    anchor_position_band = max( &
                        real(multidim_anchor_position_cells, dp) * cell_size(front_axis), &
                        multidim_anchor_position_fraction * domain_front_length)

                    anchor_quasi_stationary_candidate = &
                        this%state%flame_established .and. this%state%controller_armed .and. &
                        flame_detected .and. anchor_velocity_interior .and. &
                        response_settled .and. domain_ok .and. &
                        this%state%diag_hist_count >= max_diag_hist_size .and. &
                        abs(diag_flame_velocity_lsq) <= anchor_velocity_tolerance .and. &
                        abs(anchor_mean_front_position - this%state%front_reference_coord) <= &
                            anchor_position_band .and. &
                        anchor_inlet_window_change <= anchor_inlet_change_tolerance

                    if (anchor_quasi_stationary_candidate) then
                        this%state%stabilization_counter = &
                            this%state%stabilization_counter + 1
                    else if (.not. flame_detected .or. this%state%domain_failure .or. &
                        abs(diag_flame_velocity_lsq) > 2.0_dp * anchor_velocity_tolerance .or. &
                        abs(anchor_mean_front_position - this%state%front_reference_coord) > &
                            1.5_dp * anchor_position_band .or. &
                        anchor_inlet_window_change > 2.0_dp * anchor_inlet_change_tolerance) then
                        this%state%stabilization_counter = 0
                    end if

                    if (this%state%stabilization_counter >= multidim_anchor_success_count) &
                        call write_anchor_result_once()
                else if (this%state%flame_established .and. final_structure_stationary .and. &
                    anchor_velocity_interior .and. measurement_enabled .and. &
                    this%state%hist_count >= min_hist_for_control .and. &
                    this%state%diag_hist_count >= max_diag_hist_size .and. &
                    correction_free_time >= diagnostic_window_time .and. &
                    domain_ok .and. (.not. domain_warning) .and. &
                    abs(diag_flame_velocity_filtered) < velocity_tolerance_off .and. &
                    position_control_error == 0.0_dp) then
                    this%state%stabilization_counter = this%state%stabilization_counter + 1
                else if (control_performed .or. .not. flame_detected .or. &
                    abs(diag_flame_velocity_filtered) > velocity_tolerance_on .or. &
                    position_control_error /= 0.0_dp) then
                    this%state%stabilization_counter = 0
                end if
            else if (this%state%control_stage == stage_flamelet_ready) then
                if (this%state%flame_established .and. final_structure_stationary .and. &
                    anchor_velocity_interior .and. measurement_enabled .and. &
                    this%state%hist_count >= min_hist_for_control .and. &
                    this%state%diag_hist_count >= max_diag_hist_size .and. &
                    correction_free_time >= diagnostic_window_time .and. &
                    domain_ok .and. (.not. domain_warning) .and. &
                    abs(diag_flame_velocity_filtered) < velocity_tolerance_off .and. &
                    position_control_error == 0.0_dp) then
                    this%state%post_flamelet_hold_counter = this%state%post_flamelet_hold_counter + 1
                    this%state%stabilization_counter = this%state%post_flamelet_hold_counter
                else if (control_performed .or. .not. flame_detected .or. &
                    abs(diag_flame_velocity_filtered) > velocity_tolerance_on .or. &
                    position_control_error /= 0.0_dp) then
                    this%state%post_flamelet_hold_counter = 0
                    this%state%stabilization_counter = 0
                end if
            else if (this%state%control_stage == stage_measurement_done) then
                this%state%stabilization_counter = stable_required_count
            end if

            if (this%state%control_stage == stage_measurement_done) then
                flame_velocity_lsq = this%state%measurement_velocity_save
                flame_velocity_filtered = this%state%measurement_velocity_save
                control_velocity = this%state%measurement_velocity_save
                sl_displacement = this%state%sl_displacement_save
                linear_r2 = this%state%measurement_r2_save
                linear_rms = this%state%measurement_rms_save
                split_slope_diff = this%state%measurement_split_slope_diff_save
            end if

            active_track_number = this%state%track_counter
            call write_tracking_line()
            call write_physics_line()

            if (this%state%domain_failure) then
                call write_domain_failure_once()
                this%state%track_counter = this%state%track_counter + 1
                if (present(failed)) failed = .true.
                return
            end if

            if (this%state%control_stage == stage_measurement_failed) then
                call write_measurement_failure_once()
                this%state%track_counter = this%state%track_counter + 1
                if (present(failed)) failed = .true.
                return
            end if

            if (reset_history_after_log) call clear_control_history()

            this%state%track_counter = this%state%track_counter + 1

            if (this%state%control_stage == stage_anchor_control .and. &
                this%state%stabilization_counter >= flamelet_required_count) then
                if (active_workflow == workflow_flamelet_sl) then
                    this%state%control_stage = stage_flamelet_ready
                    this%state%stabilization_counter = 0
                    this%state%post_flamelet_hold_counter = 0
                    call clear_control_history()
                else
                    stabilized = .false.
                end if

            else if (this%state%control_stage == stage_flamelet_ready .and. &
                this%state%post_flamelet_hold_counter >= post_flamelet_hold_required_count) then
                if (enable_flamelet_output) call request_flamelet_tables_once()
                if (enable_drift_measurement .and. dimensions == 1) then
                    call capture_stabilized_scientific_state()
                    call start_drift_measurement_ramp(time)
                else if (pause_after_sl_measurement) then
                    stabilized = .true.
                end if

            else if (this%state%control_stage == stage_measurement_done) then
                if (pause_after_sl_measurement) stabilized = .true.
            end if

        end associate

    contains

                subroutine initialize_output_file()
            integer :: local_dim

            if (this%load_counter == 1) then
                flame_debug_file = 'flame_stabilization_debug.dat'
                flame_physics_file = 'flame_physics.dat'
            else
                write(flame_debug_file,'(A,I0,A)') &
                    'flame_stabilization_debug_', this%load_counter, '.dat'
                write(flame_physics_file,'(A,I0,A)') &
                    'flame_physics_', this%load_counter, '.dat'
            end if

            ! Internal stabilization/debug stream.
            open(newunit = this%state%flame_loc_unit, &
                file = trim(flame_debug_file), status = 'replace', &
                form = 'formatted')

            write(this%state%flame_loc_unit,'(A)') &
                'TITLE = "NRG flame stabilization debug history"'

            av_header = 'VARIABLES = "time"'
            do local_dim = 1, dimensions
                av_header = trim(av_header) // &
                    ' "xf_' // trim(axis_names(local_dim)) // '"'
            end do
            av_header = trim(av_header) // &
                ' "Vfl_lsq" "Vfl_filtered" "Vfl_diag_lsq" "Vfl_diag_filtered"' // &
                ' "Vfl_measurement_lsq" "Vcontrol" "pos_error" "x_ref"'
            av_header = trim(av_header) // ' "U_in_applied" "U_in_target"'
            av_header = trim(av_header) // &
                ' "dU_target" "adaptive_gain" "hist_count"' // &
                ' "measurement_on" "bracket_on" "capture_on" "emergency_on"' // &
                ' "flame_detected"'
            av_header = trim(av_header) // &
                ' "front_spread" "Qint" "Qmax" "Qvalid" "Hmax" "Tgradmax"'
            av_header = trim(av_header) // &
                ' "stage" "SL_disp" "lin_R2" "lin_RMS" "split_dV"'
            av_header = trim(av_header) // &
                ' "corr_count" "stab_count" "track_count"' // &
                ' "controller_armed" "ramp_settled" "response_settled"' // &
                ' "vel_tol_on" "vel_tol_off" "Qint_rel_trend"' // &
                ' "flame_established" "establish_window_samples"' // &
                ' "hard_recovery_corrections"' // &
                ' "containment_on" "hard_recovery_hold"' // &
                ' "hard_recovery_no_progress"' // &
                ' "final_stationary" "structure_stab_count" "anchor_interior"' // &
                ' "x_preheat" "x_T10" "x_T90" "delta_T"' // &
                ' "x_reaction_left" "x_reaction_right" "delta_reaction"' // &
                ' "d_inlet_preheat" "d_outlet_reaction"' // &
                ' "inlet_theta_max" "inlet_margin_ratio" "outlet_margin_ratio"' // &
                ' "domain_warning" "domain_ok" "domain_failure"' // &
                ' "correction_free_time"'

            write(this%state%flame_loc_unit,'(A)') trim(av_header)

            ! Stable physics / interpretation stream.
            open(newunit = this%state%physics_output_unit, &
                file = trim(flame_physics_file), status = 'replace', &
                form = 'formatted')

            write(this%state%physics_output_unit,'(A)') &
                'TITLE = "NRG flame physics history"'
            if (multidim_anchor) then
                write(this%state%physics_output_unit,'(A)') &
                    'VARIABLES = ' // &
                    '"time_s" "xf_mean_m" "Vf_mean_m_s" "U_in_m_s" ' // &
                    '"Qint_W_m" "Qmax_W_m3" "Tmax_K" "front_spread_m" ' // &
                    '"front_x05_m" "front_x50_m" "front_x95_m" ' // &
                    '"x_reaction_min_m" "x_reaction_max_m" ' // &
                    '"axial_reaction_span_m" "local_reaction_p50_m" ' // &
                    '"local_reaction_p95_m" "d_inlet_preheat_m" ' // &
                    '"d_outlet_reaction_p95_m"'
            else
                write(this%state%physics_output_unit,'(A)') &
                    'VARIABLES = ' // &
                    '"time_s" "xf_m" "Vfl_m_s" "U_in_m_s" "S_kinematic_m_s" ' // &
                    '"Qint_W_m2" "Qmax_W_m3" "Sc_m_s" "Sc_valid" ' // &
                    '"H2_consumption_flux_kg_m2_s" "H2_convective_flux_kg_m2_s" ' // &
                    '"H2_inventory_kg_m2" "total_mass_kg_m2" ' // &
                    '"rho_fresh_kg_m3" "YH2_fresh" "YH2_products" "Tmax_K" ' // &
                    '"x_preheat_m" "x_T10_m" "x_T90_m" "delta_T_m" ' // &
                    '"x_reaction_left_m" "x_reaction_right_m" ' // &
                    '"delta_reaction_m" "front_spread_m" "Hmax" ' // &
                    '"d_inlet_preheat_m" "d_outlet_reaction_m"'
            end if

            this%state%output_initialized = .true.
        end subroutine initialize_output_file


        function cell_center_coordinates(ii,jj,kk) result(xc)
            integer, intent(in) :: ii, jj, kk
            real(dp) :: xc(3)

            xc = 0.0_dp
            xc(1) = (real(ii,dp) - 0.5_dp) * cell_size(1)
            if (dimensions >= 2) xc(2) = (real(jj,dp) - 0.5_dp) * cell_size(2)
            if (dimensions >= 3) xc(3) = (real(kk,dp) - 0.5_dp) * cell_size(3)
        end function cell_center_coordinates

        function temperature_gradient_norm(ii,jj,kk) result(grad_norm)
            integer, intent(in) :: ii, jj, kk
            real(dp) :: grad_norm
            real(dp) :: g2, gd

            g2 = 0.0_dp
            if (dimensions >= 1) then
                if (ii > cons_inner_loop(1,1) .and. ii < cons_inner_loop(1,2)) then
                    gd = (this%T%s_ptr%cells(ii+1,jj,kk) - this%T%s_ptr%cells(ii-1,jj,kk)) / (2.0_dp * cell_size(1))
                    g2 = g2 + gd * gd
                end if
            end if
            if (dimensions >= 2) then
                if (jj > cons_inner_loop(2,1) .and. jj < cons_inner_loop(2,2)) then
                    gd = (this%T%s_ptr%cells(ii,jj+1,kk) - this%T%s_ptr%cells(ii,jj-1,kk)) / (2.0_dp * cell_size(2))
                    g2 = g2 + gd * gd
                end if
            end if
            if (dimensions >= 3) then
                if (kk > cons_inner_loop(3,1) .and. kk < cons_inner_loop(3,2)) then
                    gd = (this%T%s_ptr%cells(ii,jj,kk+1) - this%T%s_ptr%cells(ii,jj,kk-1)) / (2.0_dp * cell_size(3))
                    g2 = g2 + gd * gd
                end if
            end if
            grad_norm = sqrt(g2)
        end function temperature_gradient_norm

        subroutine append_front_history(t_new, s_new)
            real(dp), intent(in) :: t_new, s_new

            if (this%state%hist_count < max_hist_size) then
                this%state%hist_count = this%state%hist_count + 1
                this%state%time_hist(this%state%hist_count) = t_new
                this%state%front_coord_hist(this%state%hist_count) = s_new
            else
                this%state%time_hist(1:max_hist_size-1) = this%state%time_hist(2:max_hist_size)
                this%state%front_coord_hist(1:max_hist_size-1) = this%state%front_coord_hist(2:max_hist_size)
                this%state%time_hist(max_hist_size) = t_new
                this%state%front_coord_hist(max_hist_size) = s_new
            end if
        end subroutine append_front_history

        subroutine append_diagnostic_history(t_new, s_new)
            real(dp), intent(in) :: t_new, s_new

            if (this%state%diag_hist_count < max_diag_hist_size) then
                this%state%diag_hist_count = this%state%diag_hist_count + 1
                this%state%diag_time_hist(this%state%diag_hist_count) = t_new
                this%state%diag_front_coord_hist(this%state%diag_hist_count) = s_new
                this%state%diag_inlet_velocity_hist(this%state%diag_hist_count) = &
                    this%state%inlet_velocity_applied
            else
                this%state%diag_time_hist(1:max_diag_hist_size-1) = &
                    this%state%diag_time_hist(2:max_diag_hist_size)
                this%state%diag_front_coord_hist(1:max_diag_hist_size-1) = &
                    this%state%diag_front_coord_hist(2:max_diag_hist_size)
                this%state%diag_inlet_velocity_hist(1:max_diag_hist_size-1) = &
                    this%state%diag_inlet_velocity_hist(2:max_diag_hist_size)
                this%state%diag_time_hist(max_diag_hist_size) = t_new
                this%state%diag_front_coord_hist(max_diag_hist_size) = s_new
                this%state%diag_inlet_velocity_hist(max_diag_hist_size) = &
                    this%state%inlet_velocity_applied
            end if
        end subroutine append_diagnostic_history

        subroutine append_arming_history(metric_new, window_samples)
            real(dp), intent(in) :: metric_new
            integer, intent(in) :: window_samples

            integer :: window_use, half_samples
            integer :: first_old, first_new
            real(dp) :: mean_old, mean_new, metric_scale

            if (this%state%arm_hist_count < arming_max_window_samples) then
                this%state%arm_hist_count = this%state%arm_hist_count + 1
                this%state%arm_metric_hist(this%state%arm_hist_count) = metric_new
            else
                this%state%arm_metric_hist(1:arming_max_window_samples-1) = &
                    this%state%arm_metric_hist(2:arming_max_window_samples)
                this%state%arm_metric_hist(arming_max_window_samples) = metric_new
            end if

            window_use = max(arming_min_window_samples, &
                min(arming_max_window_samples, window_samples))

            if (this%state%arm_hist_count >= window_use) then
                half_samples = max(window_use / 2, 1)
                first_old = this%state%arm_hist_count - window_use + 1
                first_new = this%state%arm_hist_count - half_samples + 1

                mean_old = sum(this%state%arm_metric_hist( &
                    first_old:first_old + half_samples - 1)) / real(half_samples, dp)
                mean_new = sum(this%state%arm_metric_hist( &
                    first_new:this%state%arm_hist_count)) / real(half_samples, dp)

                metric_scale = max(abs(mean_old), abs(mean_new), tiny_weight)
                this%state%arm_metric_rel_trend = abs(mean_new - mean_old) / metric_scale
            else
                this%state%arm_metric_rel_trend = 1.0e30_dp
            end if
        end subroutine append_arming_history

        subroutine clear_arming_history()
            this%state%arm_hist_count = 0
            this%state%arm_metric_hist = 0.0_dp
            this%state%arm_metric_rel_trend = 1.0e30_dp
        end subroutine clear_arming_history

        subroutine reset_hard_recovery_state()
            this%state%hard_recovery_corrections = 0
            this%state%hard_recovery_no_progress_count = 0
            this%state%hard_recovery_reference_deficit = huge(1.0_dp)
            this%state%hard_recovery_reference_velocity = huge(1.0_dp)
            this%state%hard_recovery_last_evaluation_time = -huge(1.0_dp)
            this%state%hard_recovery_waiting_response = .false.
        end subroutine reset_hard_recovery_state

        subroutine clear_control_history()
            this%state%hist_count = 0
            this%state%filtered_velocity_save = 0.0_dp
            this%state%time_hist = 0.0_dp
            this%state%front_coord_hist = 0.0_dp
        end subroutine clear_control_history

        function least_squares_velocity() result(vfit)
            real(dp) :: vfit
            integer :: n
            real(dp) :: t_av, s_av, numerator, denominator

            if (this%state%hist_count < 2) then
                vfit = 0.0_dp
                return
            end if

            t_av = sum(this%state%time_hist(1:this%state%hist_count)) / real(this%state%hist_count, dp)
            s_av = sum(this%state%front_coord_hist(1:this%state%hist_count)) / real(this%state%hist_count, dp)
            numerator = 0.0_dp
            denominator = 0.0_dp
            do n = 1, this%state%hist_count
                numerator = numerator + (this%state%time_hist(n) - t_av) * (this%state%front_coord_hist(n) - s_av)
                denominator = denominator + (this%state%time_hist(n) - t_av)**2
            end do

            if (denominator > tiny(denominator)) then
                vfit = numerator / denominator
            else
                vfit = 0.0_dp
            end if
        end function least_squares_velocity

        function diagnostic_least_squares_velocity() result(vfit)
            real(dp) :: vfit
            integer :: n
            real(dp) :: t_av, s_av, numerator, denominator

            if (this%state%diag_hist_count < 2) then
                vfit = 0.0_dp
                return
            end if

            t_av = sum(this%state%diag_time_hist(1:this%state%diag_hist_count)) / real(this%state%diag_hist_count, dp)
            s_av = sum(this%state%diag_front_coord_hist(1:this%state%diag_hist_count)) / real(this%state%diag_hist_count, dp)
            numerator = 0.0_dp
            denominator = 0.0_dp
            do n = 1, this%state%diag_hist_count
                numerator = numerator + (this%state%diag_time_hist(n) - t_av) * &
                    (this%state%diag_front_coord_hist(n) - s_av)
                denominator = denominator + (this%state%diag_time_hist(n) - t_av)**2
            end do

            if (denominator > tiny(denominator)) then
                vfit = numerator / denominator
            else
                vfit = 0.0_dp
            end if
        end function diagnostic_least_squares_velocity

        function diagnostic_mean_front_position() result(x_mean)
            real(dp) :: x_mean
            if (this%state%diag_hist_count < 1) then
                x_mean = 0.0_dp
            else
                x_mean = sum(this%state%diag_front_coord_hist( &
                    1:this%state%diag_hist_count)) / real(this%state%diag_hist_count, dp)
            end if
        end function diagnostic_mean_front_position

        function diagnostic_mean_inlet_velocity() result(u_mean)
            real(dp) :: u_mean
            if (this%state%diag_hist_count < 1) then
                u_mean = this%state%inlet_velocity_applied
            else
                u_mean = sum(this%state%diag_inlet_velocity_hist( &
                    1:this%state%diag_hist_count)) / real(this%state%diag_hist_count, dp)
            end if
        end function diagnostic_mean_inlet_velocity

        function diagnostic_inlet_velocity_slope() result(uslope)
            real(dp) :: uslope
            integer :: n
            real(dp) :: t_av, u_av, numerator, denominator
            if (this%state%diag_hist_count < 2) then
                uslope = 0.0_dp
                return
            end if
            t_av = sum(this%state%diag_time_hist(1:this%state%diag_hist_count)) / &
                real(this%state%diag_hist_count, dp)
            u_av = diagnostic_mean_inlet_velocity()
            numerator = 0.0_dp
            denominator = 0.0_dp
            do n = 1, this%state%diag_hist_count
                numerator = numerator + (this%state%diag_time_hist(n) - t_av) * &
                    (this%state%diag_inlet_velocity_hist(n) - u_av)
                denominator = denominator + (this%state%diag_time_hist(n) - t_av)**2
            end do
            if (denominator > tiny(denominator)) then
                uslope = numerator / denominator
            else
                uslope = 0.0_dp
            end if
        end function diagnostic_inlet_velocity_slope

        function diagnostic_inlet_velocity_rms() result(u_rms)
            real(dp) :: u_rms, u_mean
            if (this%state%diag_hist_count < 1) then
                u_rms = 0.0_dp
                return
            end if
            u_mean = diagnostic_mean_inlet_velocity()
            u_rms = sqrt(sum((this%state%diag_inlet_velocity_hist( &
                1:this%state%diag_hist_count) - u_mean)**2) / &
                real(this%state%diag_hist_count, dp))
        end function diagnostic_inlet_velocity_rms

        function percentile_value(values, n_values, fraction) result(value_out)
            real(dp), intent(in) :: values(:)
            integer, intent(in) :: n_values
            real(dp), intent(in) :: fraction
            real(dp) :: value_out
            real(dp), allocatable :: work(:)
            real(dp) :: key, rank, blend
            integer :: ii, jj, i0, i1
            if (n_values <= 0) then
                value_out = 0.0_dp
                return
            end if
            allocate(work(n_values))
            work = values(1:n_values)
            do ii = 2, n_values
                key = work(ii)
                jj = ii - 1
                do while (jj >= 1)
                    if (work(jj) <= key) exit
                    work(jj + 1) = work(jj)
                    jj = jj - 1
                end do
                work(jj + 1) = key
            end do
            rank = 1.0_dp + min(max(fraction, 0.0_dp), 1.0_dp) * real(n_values - 1, dp)
            i0 = int(floor(rank))
            i1 = min(i0 + 1, n_values)
            blend = rank - real(i0, dp)
            value_out = (1.0_dp - blend) * work(i0) + blend * work(i1)
        end function percentile_value

        subroutine drift_linearity_diagnostics(r2_out, rms_out, split_slope_diff_out)
            real(dp), intent(out) :: r2_out, rms_out, split_slope_diff_out
            integer :: n, mid
            real(dp) :: v_all, t_av, s_av, intercept, ss_tot, ss_res, residual
            real(dp) :: v_first, v_second

            r2_out = 0.0_dp
            rms_out = huge(1.0_dp)
            split_slope_diff_out = huge(1.0_dp)
            if (this%state%hist_count < 4) return

            v_all = least_squares_velocity()
            t_av = sum(this%state%time_hist(1:this%state%hist_count)) / real(this%state%hist_count, dp)
            s_av = sum(this%state%front_coord_hist(1:this%state%hist_count)) / real(this%state%hist_count, dp)
            intercept = s_av - v_all * t_av
            ss_tot = 0.0_dp
            ss_res = 0.0_dp
            do n = 1, this%state%hist_count
                residual = this%state%front_coord_hist(n) - (intercept + v_all * this%state%time_hist(n))
                ss_res = ss_res + residual * residual
                ss_tot = ss_tot + (this%state%front_coord_hist(n) - s_av)**2
            end do
            rms_out = sqrt(ss_res / real(this%state%hist_count, dp))
            if (ss_tot > tiny(ss_tot)) then
                r2_out = max(0.0_dp, 1.0_dp - ss_res / ss_tot)
            else
                r2_out = 0.0_dp
            end if

            mid = this%state%hist_count / 2
            v_first = least_squares_velocity_range(1, mid)
            v_second = least_squares_velocity_range(mid + 1, this%state%hist_count)
            split_slope_diff_out = abs(v_second - v_first)
        end subroutine drift_linearity_diagnostics

        function least_squares_velocity_range(n_first, n_last) result(vfit)
            integer, intent(in) :: n_first, n_last
            real(dp) :: vfit
            integer :: n, n_local
            real(dp) :: t_av, s_av, numerator, denominator

            n_local = n_last - n_first + 1
            if (n_local < 2) then
                vfit = 0.0_dp
                return
            end if

            t_av = 0.0_dp
            s_av = 0.0_dp
            do n = n_first, n_last
                t_av = t_av + this%state%time_hist(n)
                s_av = s_av + this%state%front_coord_hist(n)
            end do
            t_av = t_av / real(n_local, dp)
            s_av = s_av / real(n_local, dp)

            numerator = 0.0_dp
            denominator = 0.0_dp
            do n = n_first, n_last
                numerator = numerator + (this%state%time_hist(n) - t_av) * (this%state%front_coord_hist(n) - s_av)
                denominator = denominator + (this%state%time_hist(n) - t_av)**2
            end do
            if (denominator > tiny(denominator)) then
                vfit = numerator / denominator
            else
                vfit = 0.0_dp
            end if
        end function least_squares_velocity_range

        pure function safe_window_error(value, lower_bound, upper_bound) result(error_out)
            real(dp), intent(in) :: value, lower_bound, upper_bound
            real(dp) :: error_out

            if (value < lower_bound) then
                error_out = value - lower_bound
            else if (value > upper_bound) then
                error_out = value - upper_bound
            else
                error_out = 0.0_dp
            end if
        end function safe_window_error

        pure function combined_position_error(centroid_error, envelope_error) result(error_out)
            real(dp), intent(in) :: centroid_error, envelope_error
            real(dp) :: error_out

            if (envelope_error /= 0.0_dp) then
                error_out = envelope_error
            else
                error_out = centroid_error
            end if
        end function combined_position_error

        logical function control_action_needed()
            control_action_needed = (abs(flame_velocity_filtered) > velocity_tolerance_on) .or. &
                (abs(position_control_error) > 0.0_dp .and. abs(control_velocity) > velocity_tolerance_on)
        end function control_action_needed

        subroutine update_adaptive_gain(v_current)
            real(dp), intent(in) :: v_current

            if (this%state%have_previous_control_point) then
                if (v_current * this%state%previous_control_velocity < 0.0_dp) then
                    this%state%same_sign_error_counter = 0
                    this%state%adaptive_gain = max(0.7_dp * this%state%adaptive_gain, controller_gain_min)
                else
                    this%state%same_sign_error_counter = this%state%same_sign_error_counter + 1
                    if (this%state%same_sign_error_counter >= 2) then
                        this%state%adaptive_gain = min(1.10_dp * this%state%adaptive_gain, controller_gain_max)
                    else
                        this%state%adaptive_gain = min(1.03_dp * this%state%adaptive_gain, controller_gain_max)
                    end if
                end if
            else
                this%state%same_sign_error_counter = 0
            end if
        end subroutine update_adaptive_gain

        subroutine update_bracket(u_current, v_current)
            real(dp), intent(in) :: u_current, v_current

            if (abs(v_current) <= velocity_tolerance_off) return

            if (this%state%have_previous_control_point) then
                if (abs(u_current - this%state%previous_control_inlet_velocity) >= min_secant_du .and. &
                    v_current * this%state%previous_control_velocity < 0.0_dp) then
                    this%state%bracket_u_a = this%state%previous_control_inlet_velocity
                    this%state%bracket_v_a = this%state%previous_control_velocity
                    this%state%bracket_u_b = u_current
                    this%state%bracket_v_b = v_current
                    this%state%has_bracket = .true.
                end if
            end if

            if (this%state%has_bracket) then
                if (v_current * this%state%bracket_v_a > 0.0_dp) then
                    this%state%bracket_u_a = u_current
                    this%state%bracket_v_a = v_current
                else if (v_current * this%state%bracket_v_b > 0.0_dp) then
                    this%state%bracket_u_b = u_current
                    this%state%bracket_v_b = v_current
                end if
                if (this%state%bracket_v_a * this%state%bracket_v_b > 0.0_dp) this%state%has_bracket = .false.
                if (abs(this%state%bracket_u_a - this%state%bracket_u_b) < min_secant_du) this%state%has_bracket = .false.
            end if
        end subroutine update_bracket

        subroutine choose_new_inlet_target(u_current, v_current, u_new, du_unlimited, du_final)
            real(dp), intent(in) :: u_current, v_current
            real(dp), intent(out) :: u_new, du_unlimited, du_final
            real(dp) :: dU, dV, response_slope, secant_step, bracket_target
            real(dp) :: gain_effective, max_fraction_effective, min_step_effective

            if (emergency_mode) then
                gain_effective = controller_gain_capture
                max_fraction_effective = controller_max_fraction_emergency
                min_step_effective = max(min_abs_velocity_step_capture, &
                    emergency_min_fraction * max(abs(u_current), min_abs_velocity_step))
            else if (capture_mode) then
                gain_effective = controller_gain_capture
                max_fraction_effective = controller_max_fraction_capture
                min_step_effective = min_abs_velocity_step_capture
            else
                gain_effective = this%state%adaptive_gain
                max_fraction_effective = controller_max_fraction
                min_step_effective = min_abs_velocity_step
            end if

            max_velocity_step = max(max_fraction_effective * &
                max(abs(u_current), min_abs_velocity_step), min_step_effective)

            du_unlimited = feedback_sign_default * gain_effective * v_current
            if (emergency_mode) then
                if (abs(du_unlimited) < min_step_effective) then
                    if (du_unlimited /= 0.0_dp) then
                        du_unlimited = sign(min_step_effective, du_unlimited)
                    else
                        du_unlimited = feedback_sign_default * &
                            sign(min_step_effective, v_current)
                    end if
                end if
            else if (this%state%has_bracket .and. .not. capture_mode) then
                bracket_target = 0.5_dp * (this%state%bracket_u_a + this%state%bracket_u_b)
                du_unlimited = bracket_target - u_current
            else if (use_secant_control .and. this%state%have_previous_control_point .and. .not. capture_mode) then
                dU = u_current - this%state%previous_control_inlet_velocity
                dV = v_current - this%state%previous_control_velocity
                if (abs(dU) >= min_secant_du .and. abs(dV) > velocity_tolerance_off) then
                    response_slope = dV / dU
                    if (response_slope * feedback_sign_default < 0.0_dp .and. &
                        abs(response_slope) > 1.0e-12_dp .and. &
                        abs(response_slope) < max_response_slope_abs) then
                        secant_step = -secant_relaxation * v_current / response_slope
                        if (abs(secant_step) <= 5.0_dp * max_velocity_step) then
                            du_unlimited = secant_step
                        end if
                    end if
                end if
            end if

            du_final = min(max(du_unlimited, -max_velocity_step), max_velocity_step)

            if (abs(control_velocity) > persistent_error_factor * velocity_tolerance_on .or. &
                abs(position_control_error) > 0.0_dp) then
                if (abs(du_final) < min_step_effective) then
                    if (du_unlimited /= 0.0_dp) then
                        du_final = sign(min_step_effective, du_unlimited)
                    else
                        du_final = feedback_sign_default * sign(min_step_effective, v_current)
                    end if
                end if
            end if

            u_new = max(u_current + du_final, 0.0_dp)
        end subroutine choose_new_inlet_target

        subroutine request_flamelet_tables_once()
            if (.not. this%state%flamelet_output_written) then
                this%requested_chem_table_filename = trim(this%state%chem_table_filename)
                this%requested_data_table_filename = trim(this%state%data_table_filename)
                this%flamelet_output_pending = .true.
                this%state%flamelet_output_written = .true.
            end if
        end subroutine request_flamelet_tables_once

        subroutine capture_stabilized_scientific_state()
            integer, parameter :: plateau_half_width = 2
            integer :: ii, i_fresh, i_products, i_first, i_last
            integer :: sample_count, spec_local
            real(dp) :: fresh_coord, products_coord
            real(dp) :: fuel_mass_fraction_drop, consumption_denominator

            if (.not. domain_ok .or. domain_warning) return
            if (.not. thermal_envelope_found .or. .not. reaction_envelope_found) return

            fresh_coord = 0.5_dp * (domain_boundary_min + x_preheat)
            products_coord = 0.5_dp * (x_reaction_right + domain_boundary_max)

            i_fresh = nint(fresh_coord / cell_size(front_axis) + 0.5_dp)
            i_products = nint(products_coord / cell_size(front_axis) + 0.5_dp)
            i_fresh = min(max(i_fresh, cons_inner_loop(1,1)), cons_inner_loop(1,2))
            i_products = min(max(i_products, cons_inner_loop(1,1)), cons_inner_loop(1,2))

            this%state%scientific_rho_fresh = 0.0_dp
            this%state%scientific_y_h2_fresh = 0.0_dp
            i_first = max(i_fresh - plateau_half_width, cons_inner_loop(1,1))
            i_last = min(i_fresh + plateau_half_width, cons_inner_loop(1,2))
            sample_count = i_last - i_first + 1
            do ii = i_first, i_last
                this%state%scientific_rho_fresh = this%state%scientific_rho_fresh + &
                    this%rho%s_ptr%cells(ii,1,1)
                if (H2_index >= 1 .and. H2_index <= species_number) then
                    this%state%scientific_y_h2_fresh = &
                        this%state%scientific_y_h2_fresh + &
                        this%Y%v_ptr%pr(H2_index)%cells(ii,1,1)
                end if
            end do
            this%state%scientific_rho_fresh = this%state%scientific_rho_fresh / real(sample_count,dp)
            if (H2_index >= 1 .and. H2_index <= species_number) then
                this%state%scientific_y_h2_fresh = &
                    this%state%scientific_y_h2_fresh / real(sample_count,dp)
            end if

            this%state%scientific_rho_products = 0.0_dp
            this%state%scientific_t_products = 0.0_dp
            this%state%stabilized_product_mass_fractions = 0.0_dp
            i_first = max(i_products - plateau_half_width, cons_inner_loop(1,1))
            i_last = min(i_products + plateau_half_width, cons_inner_loop(1,2))
            sample_count = i_last - i_first + 1
            do ii = i_first, i_last
                this%state%scientific_rho_products = this%state%scientific_rho_products + &
                    this%rho%s_ptr%cells(ii,1,1)
                this%state%scientific_t_products = this%state%scientific_t_products + &
                    this%T%s_ptr%cells(ii,1,1)
                do spec_local = 1, species_number
                    this%state%stabilized_product_mass_fractions(spec_local) = &
                        this%state%stabilized_product_mass_fractions(spec_local) + &
                        this%Y%v_ptr%pr(spec_local)%cells(ii,1,1)
                end do
            end do

            this%state%scientific_rho_products = &
                this%state%scientific_rho_products / real(sample_count,dp)
            this%state%scientific_t_products = &
                this%state%scientific_t_products / real(sample_count,dp)
            this%state%stabilized_product_mass_fractions = &
                this%state%stabilized_product_mass_fractions / real(sample_count,dp)

            this%state%scientific_h2_consumption_flux = 0.0_dp
            this%state%scientific_consumption_speed = 0.0_dp
            this%state%scientific_consumption_speed_valid = .false.
            this%state%scientific_y_h2_products = 0.0_dp

            if (H2_index >= 1 .and. H2_index <= species_number) then
                this%state%scientific_y_h2_products = &
                    this%state%stabilized_product_mass_fractions(H2_index)

                ! The current diagnostic is deliberately restricted to the
                ! planar 1D formulation. Curved flames require area weighting
                ! and a flame-surface normalization.
                if (dimensions == 1 .and. &
                    trim(this%domain%get_coordinate_system_name()) == 'cartesian') then
                    do ii = cons_inner_loop(1,1), cons_inner_loop(1,2)
                        this%state%scientific_h2_consumption_flux = &
                            this%state%scientific_h2_consumption_flux - &
                            this%Y_prod_chem%v_ptr%pr(H2_index)%cells(ii,1,1) * &
                            cell_size(front_axis)
                    end do

                    fuel_mass_fraction_drop = &
                        this%state%scientific_y_h2_fresh - &
                        this%state%scientific_y_h2_products
                    consumption_denominator = &
                        this%state%scientific_rho_fresh * fuel_mass_fraction_drop

                    if (this%state%scientific_h2_consumption_flux > 0.0_dp .and. &
                        consumption_denominator > tiny_weight) then
                        this%state%scientific_consumption_speed = &
                            this%state%scientific_h2_consumption_flux / &
                            consumption_denominator
                        this%state%scientific_consumption_speed_valid = .true.
                    end if
                end if
            end if

            this%state%scientific_time_stabilized = time
            this%state%scientific_u_anchor = this%state%inlet_velocity_target
            this%state%scientific_flame_thickness = max(x_thermal_high - x_thermal_low, 0.0_dp)
            this%state%scientific_preheat_zone = max(x_thermal_low - x_preheat, 0.0_dp)
            this%state%scientific_reaction_zone = reaction_thickness
            this%state%scientific_qint = heat_release_integral
            this%state%scientific_qmax = heat_release_max
            this%state%scientific_d_inlet_preheat = inlet_preheat_distance
            this%state%scientific_d_outlet_reaction = outlet_reaction_distance
            this%state%scientific_state_captured = .true.
        end subroutine capture_stabilized_scientific_state

        subroutine write_anchor_result_once()
            integer :: result_unit
            character(len=240) :: result_title
            if (this%state%anchor_result_written) return
            result_title = trim(this%flame_stabilization%get_scientific_result_title())
            if (len_trim(result_title) == 0) result_title = 'NRG flame anchor result'
            open(newunit = result_unit, file = 'flame_anchor_result.json', &
                status = 'replace', form = 'formatted')
            write(result_unit,'(A)') '{'
            write(result_unit,'(A)') '  "schema": "nrg.flame_anchor_result.v1",'
            write(result_unit,'(A)') '  "status": "anchored_quasi_stationary",'
            write(result_unit,'(A)') '  "title": "' // trim(result_title) // '",'
            write(result_unit,'(A)') '  "provenance": {'
            write(result_unit,'(A)') '    "solver": "fds_low_mach",'
            write(result_unit,'(A)') '    "control_mode": "anchor",'
            write(result_unit,'(A)') '    "setup": "' // &
                trim(this%flame_stabilization%get_case_setup()) // '",'
            write(result_unit,'(A)') '    "coordinate_system": "' // &
                trim(this%domain%get_coordinate_system_name()) // '",'
            write(result_unit,'(A)') '    "chemical_mechanism": "' // &
                trim(this%chem%chem_ptr%get_chemical_mechanism()) // '",'
            write(result_unit,'(A,ES24.16,A)') '    "hydrogen_percent": ', &
                this%state%scientific_h2_percent, ','
            write(result_unit,'(A,ES24.16,A)') '    "pressure_Pa": ', this%state%inlet_pressure, ','
            write(result_unit,'(A,ES24.16)') '    "grid_spacing_m": ', cell_size(front_axis)
            write(result_unit,'(A)') '  },'
            write(result_unit,'(A)') '  "result": {'
            write(result_unit,'(A,ES24.16,A)') '    "time_anchored_s": ', time, ','
            write(result_unit,'(A,ES24.16,A)') '    "mean_front_position_m": ', &
                anchor_mean_front_position, ','
            write(result_unit,'(A,ES24.16,A)') '    "reference_front_position_m": ', &
                this%state%front_reference_coord, ','
            write(result_unit,'(A,ES24.16,A)') '    "anchor_position_band_m": ', anchor_position_band, ','
            write(result_unit,'(A,ES24.16,A)') '    "mean_front_velocity_m_s": ', &
                diag_flame_velocity_lsq, ','
            write(result_unit,'(A,ES24.16,A)') '    "mean_inlet_velocity_m_s": ', &
                anchor_mean_inlet_velocity, ','
            write(result_unit,'(A,ES24.16,A)') '    "rms_inlet_velocity_m_s": ', &
                anchor_inlet_velocity_rms, ','
            write(result_unit,'(A,ES24.16,A)') '    "inlet_velocity_trend_m_s2": ', &
                anchor_inlet_velocity_slope, ','
            write(result_unit,'(A,ES24.16,A)') '    "observation_window_s": ', diagnostic_window_time, ','
            write(result_unit,'(A,ES24.16,A)') '    "front_x05_m": ', front_x05, ','
            write(result_unit,'(A,ES24.16,A)') '    "front_x50_m": ', front_x50, ','
            write(result_unit,'(A,ES24.16,A)') '    "front_x95_m": ', front_x95, ','
            write(result_unit,'(A,ES24.16,A)') '    "local_reaction_thickness_p50_m": ', &
                local_reaction_thickness_p50, ','
            write(result_unit,'(A,ES24.16,A)') '    "local_reaction_thickness_p95_m": ', &
                local_reaction_thickness_p95, ','
            write(result_unit,'(A,ES24.16,A)') '    "global_reaction_tip_m": ', x_reaction_right, ','
            write(result_unit,'(A,ES24.16,A)') '    "integrated_heat_release_W_m": ', &
                heat_release_integral, ','
            write(result_unit,'(A,ES24.16)') '    "peak_heat_release_W_m3": ', heat_release_max
            write(result_unit,'(A)') '  }'
            write(result_unit,'(A)') '}'
            close(result_unit)
            this%state%anchor_result_written = .true.
        end subroutine write_anchor_result_once


        subroutine write_measurement_failure_once()
            integer :: result_unit
            character(len=240) :: result_title

            if (this%state%failure_report_written) return
            result_title = trim(this%flame_stabilization%get_scientific_result_title())
            if (len_trim(result_title) == 0) result_title = 'NRG laminar flame result'

            open(newunit = result_unit, file = 'laminar_flame_result.json', &
                status = 'replace', form = 'formatted')
            write(result_unit,'(A)') '{'
            write(result_unit,'(A)') '  "schema": "nrg.laminar_flame_result.v1",'
            write(result_unit,'(A)') '  "status": "measurement_failed",'
            write(result_unit,'(A)') '  "title": "' // trim(result_title) // '",'
            write(result_unit,'(A)') '  "provenance": {'
            write(result_unit,'(A)') '    "solver": "fds_low_mach",'
            write(result_unit,'(A)') '    "control_mode": "laminar_burning_velocity",'
            write(result_unit,'(A)') '    "setup": "' // &
                trim(this%flame_stabilization%get_case_setup()) // '",'
            write(result_unit,'(A)') '    "coordinate_system": "' // &
                trim(this%domain%get_coordinate_system_name()) // '",'
            write(result_unit,'(A)') '    "chemical_mechanism": "' // &
                trim(this%chem%chem_ptr%get_chemical_mechanism()) // '",'
            write(result_unit,'(A,ES24.16,A)') '    "hydrogen_percent": ', &
                this%state%scientific_h2_percent, ','
            write(result_unit,'(A,ES24.16,A)') '    "fresh_temperature_K": ', &
                this%state%inlet_temperature, ','
            write(result_unit,'(A,ES24.16,A)') '    "pressure_Pa": ', &
                this%state%inlet_pressure, ','
            write(result_unit,'(A,ES24.16)') '    "grid_spacing_m": ', cell_size(front_axis)
            write(result_unit,'(A)') '  },'
            write(result_unit,'(A)') '  "failure": {'
            write(result_unit,'(A,ES24.16,A)') '    "time_s": ', time, ','
            write(result_unit,'(A,I0,A)') '    "measurement_attempt": ', this%state%measurement_attempt, ','
            write(result_unit,'(A,ES24.16,A)') '    "measurement_delta_m_s": ', &
                this%state%current_measurement_delta, ','
            write(result_unit,'(A,ES24.16,A)') '    "measurement_delta_sign": ', &
                this%state%measurement_delta_sign_save, ','
            write(result_unit,'(A,ES24.16,A)') '    "measurement_elapsed_s": ', measurement_elapsed, ','
            write(result_unit,'(A,ES24.16,A)') '    "measurement_displacement_m": ', measurement_displacement, ','
            write(result_unit,'(A,ES24.16,A)') '    "measurement_target_displacement_m": ', &
                measurement_target_displacement, ','
            write(result_unit,'(A,ES24.16,A)') '    "measurement_velocity_m_s": ', &
                this%state%measurement_velocity_save, ','
            write(result_unit,'(A,ES24.16,A)') '    "r2": ', this%state%measurement_r2_save, ','
            write(result_unit,'(A,ES24.16,A)') '    "rms_m": ', this%state%measurement_rms_save, ','
            write(result_unit,'(A,ES24.16)') '    "split_velocity_difference_m_s": ', &
                this%state%measurement_split_slope_diff_save
            write(result_unit,'(A)') '  },'
            write(result_unit,'(A)') '  "diagnostics": {'
            write(result_unit,'(A,ES24.16,A)') '    "anchor_velocity_m_s": ', &
                this%state%scientific_u_anchor, ','
            if (this%state%scientific_consumption_speed_valid) then
                write(result_unit,'(A)') '    "consumption_speed_available": true,'
                write(result_unit,'(A,ES24.16,A)') '    "consumption_speed_m_s": ', &
                    this%state%scientific_consumption_speed, ','
                write(result_unit,'(A,ES24.16)') &
                    '    "integrated_h2_consumption_kg_m2_s": ', &
                    this%state%scientific_h2_consumption_flux
            else
                write(result_unit,'(A)') '    "consumption_speed_available": false,'
                write(result_unit,'(A)') '    "consumption_speed_m_s": null,'
                write(result_unit,'(A)') &
                    '    "integrated_h2_consumption_kg_m2_s": null'
            end if
            write(result_unit,'(A)') '  }'
            write(result_unit,'(A)') '}'
            close(result_unit)
            this%state%failure_report_written = .true.
        end subroutine write_measurement_failure_once

        subroutine write_domain_failure_once()
            integer :: result_unit
            character(len=240) :: result_title
            character(len=40) :: failure_status
            character(len=40) :: result_file, result_schema, control_mode

            if (this%state%failure_report_written) return
            result_title = trim(this%flame_stabilization%get_scientific_result_title())
            if (len_trim(result_title) == 0) result_title = 'NRG laminar flame result'

            if (physical_boundary_violation) then
                failure_status = 'physical_boundary_reached'
            else if (.not. hard_domain_fits) then
                failure_status = 'hard_envelope_does_not_fit'
            else
                failure_status = 'hard_domain_recovery_failed'
            end if

            if (this%flame_stabilization%is_anchor()) then
                result_file = 'flame_anchor_result.json'
                result_schema = 'nrg.flame_anchor_result.v1'
                control_mode = 'anchor'
            else
                result_file = 'laminar_flame_result.json'
                result_schema = 'nrg.laminar_flame_result.v1'
                control_mode = 'laminar_burning_velocity'
            end if

            open(newunit = result_unit, file = trim(result_file), &
                status = 'replace', form = 'formatted')
            write(result_unit,'(A)') '{'
            write(result_unit,'(A)') '  "schema": "' // trim(result_schema) // '",'
            write(result_unit,'(A)') '  "status": "' // trim(failure_status) // '",'
            write(result_unit,'(A)') '  "title": "' // trim(result_title) // '",'
            write(result_unit,'(A)') '  "provenance": {'
            write(result_unit,'(A)') '    "solver": "fds_low_mach",'
            write(result_unit,'(A)') '    "control_mode": "' // trim(control_mode) // '",'
            write(result_unit,'(A)') '    "setup": "' // &
                trim(this%flame_stabilization%get_case_setup()) // '",'
            write(result_unit,'(A)') '    "coordinate_system": "' // &
                trim(this%domain%get_coordinate_system_name()) // '",'
            write(result_unit,'(A)') '    "chemical_mechanism": "' // &
                trim(this%chem%chem_ptr%get_chemical_mechanism()) // '",'
            write(result_unit,'(A,ES24.16,A)') '    "hydrogen_percent": ', &
                this%state%scientific_h2_percent, ','
            write(result_unit,'(A,ES24.16,A)') '    "fresh_temperature_K": ', &
                this%state%inlet_temperature, ','
            write(result_unit,'(A,ES24.16,A)') '    "pressure_Pa": ', &
                this%state%inlet_pressure, ','
            write(result_unit,'(A,ES24.16)') '    "grid_spacing_m": ', cell_size(front_axis)
            write(result_unit,'(A)') '  },'
            write(result_unit,'(A)') '  "failure": {'
            write(result_unit,'(A,ES24.16,A)') '    "time_s": ', time, ','
            write(result_unit,'(A,I0,A)') '    "hard_recovery_corrections": ', &
                this%state%hard_recovery_corrections, ','
            write(result_unit,'(A,I0,A)') '    "hard_recovery_no_progress": ', &
                this%state%hard_recovery_no_progress_count, ','
            write(result_unit,'(A,ES24.16,A)') '    "establishment_relative_trend": ', &
                this%state%arm_metric_rel_trend, ','
            write(result_unit,'(A,ES24.16,A)') '    "x_preheat_m": ', x_preheat, ','
            write(result_unit,'(A,ES24.16,A)') '    "thermal_thickness_m": ', &
                thermal_thickness, ','
            write(result_unit,'(A,ES24.16,A)') '    "x_reaction_left_m": ', &
                x_reaction_left, ','
            write(result_unit,'(A,ES24.16,A)') '    "x_reaction_right_m": ', &
                x_reaction_right, ','
            write(result_unit,'(A,ES24.16,A)') '    "reaction_thickness_m": ', &
                reaction_thickness, ','
            write(result_unit,'(A,ES24.16,A)') '    "inlet_preheat_distance_m": ', &
                inlet_preheat_distance, ','
            write(result_unit,'(A,ES24.16,A)') '    "outlet_reaction_distance_m": ', &
                outlet_reaction_distance, ','
            write(result_unit,'(A,ES24.16,A)') '    "inlet_theta_max": ', &
                inlet_theta_max, ','
            write(result_unit,'(A,ES24.16)') '    "hard_required_length_m": ', &
                hard_required_length
            write(result_unit,'(A)') '  }'
            write(result_unit,'(A)') '}'
            close(result_unit)

            this%state%failure_report_written = .true.
        end subroutine write_domain_failure_once

        subroutine write_laminar_velocity_once()
            integer :: result_unit, spec_local
            real(dp) :: expansion_ratio
            real(dp) :: consumption_to_lbv_ratio
            real(dp) :: consumption_lbv_relative_difference
            logical :: consumption_lbv_comparison_valid
            character(len=240) :: result_title
            character(len=20) :: specie_name
            character(len=20) :: mechanism_name
            character(len=20) :: coordinate_name
            character(len=80) :: case_setup_name

            if (this%state%sl_output_written) return
            if (.not. this%state%scientific_state_captured) then
                error stop 'flame stabilization: scientific state was not captured before LBV measurement'
            end if

            result_title = trim(this%flame_stabilization%get_scientific_result_title())
            if (len_trim(result_title) == 0) result_title = 'NRG laminar flame result'
            mechanism_name = trim(this%chem%chem_ptr%get_chemical_mechanism())
            coordinate_name = trim(this%domain%get_coordinate_system_name())
            case_setup_name = trim(this%flame_stabilization%get_case_setup())

            expansion_ratio = this%state%scientific_rho_fresh / &
                max(this%state%scientific_rho_products, tiny_weight)

            consumption_to_lbv_ratio = 0.0_dp
            consumption_lbv_relative_difference = 0.0_dp
            consumption_lbv_comparison_valid = &
                this%state%scientific_consumption_speed_valid .and. &
                abs(this%state%sl_displacement_save) > tiny_weight
            if (consumption_lbv_comparison_valid) then
                consumption_to_lbv_ratio = &
                    this%state%scientific_consumption_speed / &
                    this%state%sl_displacement_save
                consumption_lbv_relative_difference = &
                    abs(this%state%scientific_consumption_speed - &
                        this%state%sl_displacement_save) / &
                    abs(this%state%sl_displacement_save)
            end if

            open(newunit = result_unit, file = 'laminar_flame_result.json', &
                status = 'replace', form = 'formatted')

            write(result_unit,'(A)') '{'
            write(result_unit,'(A)') '  "schema": "nrg.laminar_flame_result.v1",'
            write(result_unit,'(A)') '  "status": "success",'
            write(result_unit,'(A)') '  "title": "' // trim(result_title) // '",'
            write(result_unit,'(A)') '  "provenance": {'
            write(result_unit,'(A)') '    "solver": "fds_low_mach",'
            write(result_unit,'(A)') '    "control_mode": "laminar_burning_velocity",'
            write(result_unit,'(A)') '    "setup": "' // trim(case_setup_name) // '",'
            write(result_unit,'(A)') '    "coordinate_system": "' // trim(coordinate_name) // '",'
            write(result_unit,'(A)') '    "chemical_mechanism": "' // trim(mechanism_name) // '",'
            write(result_unit,'(A,ES24.16,A)') '    "hydrogen_percent": ', &
                this%state%scientific_h2_percent, ','
            write(result_unit,'(A,ES24.16,A)') '    "fresh_temperature_K": ', &
                this%state%inlet_temperature, ','
            write(result_unit,'(A,ES24.16,A)') '    "pressure_Pa": ', &
                this%state%inlet_pressure, ','
            write(result_unit,'(A,ES24.16)') '    "grid_spacing_m": ', cell_size(front_axis)
            write(result_unit,'(A)') '  },'

            write(result_unit,'(A)') '  "result": {'
            write(result_unit,'(A,ES24.16,A)') '    "laminar_burning_velocity_m_s": ', &
                this%state%sl_displacement_save, ','
            if (this%state%scientific_consumption_speed_valid) then
                write(result_unit,'(A,ES24.16,A)') '    "consumption_speed_m_s": ', &
                    this%state%scientific_consumption_speed, ','
            else
                write(result_unit,'(A)') '    "consumption_speed_m_s": null,'
            end if
            write(result_unit,'(A,ES24.16,A)') '    "anchor_velocity_m_s": ', &
                this%state%scientific_u_anchor, ','
            write(result_unit,'(A,ES24.16,A)') '    "flame_thickness_mm": ', &
                1.0e3_dp * this%state%scientific_flame_thickness, ','
            write(result_unit,'(A,ES24.16,A)') '    "preheat_zone_mm": ', &
                1.0e3_dp * this%state%scientific_preheat_zone, ','
            write(result_unit,'(A,ES24.16,A)') '    "energy_release_zone_mm": ', &
                1.0e3_dp * this%state%scientific_reaction_zone, ','
            write(result_unit,'(A,ES24.16,A)') '    "integrated_heat_release_W_m2": ', &
                this%state%scientific_qint, ','
            write(result_unit,'(A,ES24.16,A)') '    "peak_heat_release_W_m3": ', &
                this%state%scientific_qmax, ','
            write(result_unit,'(A,ES24.16)') '    "product_temperature_K": ', &
                this%state%scientific_t_products
            write(result_unit,'(A)') '  },'

            write(result_unit,'(A)') '  "states": {'
            write(result_unit,'(A,ES24.16,A)') '    "fresh_density_kg_m3": ', &
                this%state%scientific_rho_fresh, ','
            write(result_unit,'(A,ES24.16,A)') '    "product_density_kg_m3": ', &
                this%state%scientific_rho_products, ','
            write(result_unit,'(A,ES24.16,A)') '    "expansion_ratio": ', expansion_ratio, ','
            write(result_unit,'(A)') '    "product_mass_fractions": {'
            do spec_local = 1, species_number
                specie_name = this%chem%chem_ptr%get_chemical_specie_name(spec_local)
                if (spec_local < species_number) then
                    write(result_unit,'(A,A,A,ES24.16,A)') '      "', &
                        trim(specie_name), '": ', &
                        this%state%stabilized_product_mass_fractions(spec_local), ','
                else
                    write(result_unit,'(A,A,A,ES24.16)') '      "', &
                        trim(specie_name), '": ', &
                        this%state%stabilized_product_mass_fractions(spec_local)
                end if
            end do
            write(result_unit,'(A)') '    }'
            write(result_unit,'(A)') '  },'

            write(result_unit,'(A)') '  "stabilization": {'
            write(result_unit,'(A,ES24.16,A)') '    "time_stabilized_s": ', &
                this%state%scientific_time_stabilized, ','
            write(result_unit,'(A,ES24.16,A)') '    "time_result_s": ', time, ','
            write(result_unit,'(A,ES24.16,A)') '    "inlet_preheat_distance_mm": ', &
                1.0e3_dp * this%state%scientific_d_inlet_preheat, ','
            write(result_unit,'(A,ES24.16)') '    "outlet_reaction_distance_mm": ', &
                1.0e3_dp * this%state%scientific_d_outlet_reaction
            write(result_unit,'(A)') '  },'

            write(result_unit,'(A)') '  "validation": {'
            if (this%state%scientific_consumption_speed_valid) then
                write(result_unit,'(A)') '    "consumption_speed_available": true,'
                write(result_unit,'(A,ES24.16,A)') &
                    '    "integrated_h2_consumption_kg_m2_s": ', &
                    this%state%scientific_h2_consumption_flux, ','
                write(result_unit,'(A,ES24.16,A)') &
                    '    "fresh_h2_mass_fraction": ', &
                    this%state%scientific_y_h2_fresh, ','
                write(result_unit,'(A,ES24.16,A)') &
                    '    "product_h2_mass_fraction": ', &
                    this%state%scientific_y_h2_products, ','
            else
                write(result_unit,'(A)') '    "consumption_speed_available": false,'
                write(result_unit,'(A)') &
                    '    "integrated_h2_consumption_kg_m2_s": null,'
                write(result_unit,'(A)') '    "fresh_h2_mass_fraction": null,'
                write(result_unit,'(A)') '    "product_h2_mass_fraction": null,'
            end if

            if (consumption_lbv_comparison_valid) then
                write(result_unit,'(A,ES24.16,A)') &
                    '    "consumption_to_lbv_ratio": ', &
                    consumption_to_lbv_ratio, ','
                write(result_unit,'(A,ES24.16)') &
                    '    "consumption_lbv_relative_difference": ', &
                    consumption_lbv_relative_difference
            else
                write(result_unit,'(A)') '    "consumption_to_lbv_ratio": null,'
                write(result_unit,'(A)') &
                    '    "consumption_lbv_relative_difference": null'
            end if
            write(result_unit,'(A)') '  },'

            write(result_unit,'(A)') '  "measurement": {'
            write(result_unit,'(A,I0,A)') '    "attempt": ', this%state%measurement_attempt, ','
            write(result_unit,'(A,ES24.16,A)') '    "inlet_velocity_m_s": ', &
                this%state%measurement_inlet_velocity, ','
            write(result_unit,'(A,ES24.16,A)') '    "flame_velocity_m_s": ', &
                this%state%measurement_velocity_save, ','
            write(result_unit,'(A,ES24.16,A)') '    "perturbation_m_s": ', &
                this%state%measurement_delta_sign_save * this%state%current_measurement_delta, ','
            write(result_unit,'(A,ES24.16,A)') '    "r2": ', this%state%measurement_r2_save, ','
            write(result_unit,'(A,ES24.16,A)') '    "rms_mm": ', &
                1.0e3_dp * this%state%measurement_rms_save, ','
            write(result_unit,'(A,ES24.16)') '    "split_velocity_difference_m_s": ', &
                this%state%measurement_split_slope_diff_save
            write(result_unit,'(A)') '  }'
            write(result_unit,'(A)') '}'

            close(result_unit)
            this%state%sl_output_written = .true.
        end subroutine write_laminar_velocity_once

        subroutine start_drift_measurement_ramp(t_now)
            real(dp), intent(in) :: t_now

            this%state%stabilized_inlet_velocity = this%state%inlet_velocity_target
            this%state%measurement_attempt = 1

            measurement_target_displacement = &
                measurement_min_displacement_cells * cell_size(front_axis)
            measurement_delta_max = measurement_delta_max_fraction * &
                max(abs(this%state%stabilized_inlet_velocity), min_abs_velocity_step)

            this%state%current_measurement_delta = measurement_delta_fraction * &
                max(abs(this%state%stabilized_inlet_velocity), min_abs_velocity_step)
            this%state%current_measurement_delta = min( &
                this%state%current_measurement_delta, measurement_delta_max)

            inlet_measurement_room = max( &
                inlet_preheat_distance - hard_inlet_margin, 0.0_dp)
            outlet_measurement_room = max( &
                outlet_reaction_distance - hard_outlet_margin, 0.0_dp)

            if (outlet_measurement_room >= inlet_measurement_room) then
                this%state%measurement_delta_sign_save = 1.0_dp
            else
                this%state%measurement_delta_sign_save = -1.0_dp
            end if

            this%state%measurement_inlet_velocity = max( &
                this%state%stabilized_inlet_velocity + &
                this%state%measurement_delta_sign_save * &
                    this%state%current_measurement_delta, 0.0_dp)
            this%state%inlet_velocity_target = this%state%measurement_inlet_velocity
            this%state%ramp_start_velocity = this%state%inlet_velocity_applied
            this%state%ramp_start_time = t_now
            this%state%active_inlet_ramp_time = measurement_ramp_time
            this%state%sl_displacement_save = 0.0_dp
            this%state%measurement_velocity_save = 0.0_dp
            this%state%measurement_r2_save = 0.0_dp
            this%state%measurement_rms_save = 0.0_dp
            this%state%measurement_split_slope_diff_save = 0.0_dp
            this%state%control_stage = stage_measurement_ramp
            this%state%stabilization_counter = 0
            this%state%post_flamelet_hold_counter = 0
            call clear_control_history()
        end subroutine start_drift_measurement_ramp

        subroutine write_physics_line()
            integer, parameter :: plateau_half_width = 2
            integer :: ii, i_fresh, i_products, i_first, i_last
            integer :: sample_count
            real(dp) :: fresh_coord, products_coord
            real(dp) :: physics_vfl, physics_kinematic_speed
            real(dp) :: physics_h2_consumption_flux
            real(dp) :: physics_h2_convective_flux
            real(dp) :: physics_h2_inventory, physics_total_mass
            real(dp) :: physics_rho_fresh
            real(dp) :: physics_y_h2_fresh, physics_y_h2_products
            real(dp) :: physics_consumption_speed
            real(dp) :: fuel_mass_fraction_drop, consumption_denominator
            real(dp) :: consumption_valid_flag
            logical :: consumption_valid

            physics_vfl = diag_flame_velocity_lsq
            physics_kinematic_speed = &
                this%state%inlet_velocity_applied - physics_vfl

            if (multidim_anchor) then
                write(this%state%physics_output_unit,'(100E20.12)') &
                    time, current_front_coord, physics_vfl, &
                    this%state%inlet_velocity_applied, heat_release_integral, &
                    heat_release_max, temperature_max, front_spread, &
                    front_x05, front_x50, front_x95, x_reaction_left, &
                    x_reaction_right, reaction_thickness, &
                    local_reaction_thickness_p50, local_reaction_thickness_p95, &
                    inlet_preheat_distance, outlet_reaction_distance
                return
            end if

            physics_h2_consumption_flux = 0.0_dp
            physics_h2_convective_flux = 0.0_dp
            physics_h2_inventory = 0.0_dp
            physics_total_mass = 0.0_dp
            physics_rho_fresh = 0.0_dp
            physics_y_h2_fresh = 0.0_dp
            physics_y_h2_products = 0.0_dp
            physics_consumption_speed = 0.0_dp
            consumption_valid = .false.

            do ii = cons_inner_loop(1,1), cons_inner_loop(1,2)
                if (this%boundary%bc_ptr%bc_markers(ii,1,1) /= 0) cycle

                physics_total_mass = physics_total_mass + &
                    this%rho%s_ptr%cells(ii,1,1) * cell_size(front_axis)

                if (H2_index >= 1 .and. H2_index <= species_number) then
                    physics_h2_inventory = physics_h2_inventory + &
                        this%rho%s_ptr%cells(ii,1,1) * &
                        this%Y%v_ptr%pr(H2_index)%cells(ii,1,1) * &
                        cell_size(front_axis)

                    physics_h2_consumption_flux = &
                        physics_h2_consumption_flux - &
                        this%Y_prod_chem%v_ptr%pr(H2_index)%cells(ii,1,1) * &
                        cell_size(front_axis)
                end if
            end do

            if (physical_envelope_available .and. &
                H2_index >= 1 .and. H2_index <= species_number) then

                fresh_coord = 0.5_dp * (domain_boundary_min + x_preheat)
                products_coord = &
                    0.5_dp * (x_reaction_right + domain_boundary_max)

                i_fresh = nint( &
                    fresh_coord / cell_size(front_axis) + 0.5_dp)
                i_products = nint( &
                    products_coord / cell_size(front_axis) + 0.5_dp)

                i_fresh = min(max(i_fresh, cons_inner_loop(1,1)), &
                    cons_inner_loop(1,2))
                i_products = min(max(i_products, cons_inner_loop(1,1)), &
                    cons_inner_loop(1,2))

                i_first = max(i_fresh - plateau_half_width, &
                    cons_inner_loop(1,1))
                i_last = min(i_fresh + plateau_half_width, &
                    cons_inner_loop(1,2))
                sample_count = i_last - i_first + 1

                do ii = i_first, i_last
                    physics_rho_fresh = physics_rho_fresh + &
                        this%rho%s_ptr%cells(ii,1,1)
                    physics_y_h2_fresh = physics_y_h2_fresh + &
                        this%Y%v_ptr%pr(H2_index)%cells(ii,1,1)
                end do

                physics_rho_fresh = physics_rho_fresh / &
                    real(sample_count, dp)
                physics_y_h2_fresh = physics_y_h2_fresh / &
                    real(sample_count, dp)

                i_first = max(i_products - plateau_half_width, &
                    cons_inner_loop(1,1))
                i_last = min(i_products + plateau_half_width, &
                    cons_inner_loop(1,2))
                sample_count = i_last - i_first + 1

                do ii = i_first, i_last
                    physics_y_h2_products = physics_y_h2_products + &
                        this%Y%v_ptr%pr(H2_index)%cells(ii,1,1)
                end do
                physics_y_h2_products = physics_y_h2_products / &
                    real(sample_count, dp)

                physics_h2_convective_flux = physics_rho_fresh * &
                    this%state%inlet_velocity_applied * &
                    physics_y_h2_fresh

                fuel_mass_fraction_drop = &
                    physics_y_h2_fresh - physics_y_h2_products
                consumption_denominator = &
                    physics_rho_fresh * fuel_mass_fraction_drop

                if (dimensions == 1 .and. &
                    trim(this%domain%get_coordinate_system_name()) == &
                        'cartesian' .and. &
                    physics_h2_consumption_flux > 0.0_dp .and. &
                    consumption_denominator > tiny_weight) then
                    physics_consumption_speed = &
                        physics_h2_consumption_flux / &
                        consumption_denominator
                    consumption_valid = .true.
                end if
            end if

            consumption_valid_flag = merge(1.0_dp, 0.0_dp, consumption_valid)

            write(this%state%physics_output_unit,'(100E20.12)') &
                time, current_front_coord, &
                physics_vfl, this%state%inlet_velocity_applied, &
                physics_kinematic_speed, &
                heat_release_integral, heat_release_max, &
                physics_consumption_speed, consumption_valid_flag, &
                physics_h2_consumption_flux, physics_h2_convective_flux, &
                physics_h2_inventory, physics_total_mass, &
                physics_rho_fresh, physics_y_h2_fresh, &
                physics_y_h2_products, temperature_max, &
                x_preheat, x_thermal_low, x_thermal_high, thermal_thickness, &
                x_reaction_left, x_reaction_right, reaction_thickness, &
                front_spread, H_max, &
                inlet_preheat_distance, outlet_reaction_distance
        end subroutine write_physics_line


        subroutine write_tracking_line()
            real(dp) :: measurement_flag, bracket_flag, capture_flag, emergency_flag, flame_detected_flag
            real(dp) :: controller_armed_flag, ramp_settled_flag, response_settled_flag
            real(dp) :: flame_established_flag
            real(dp) :: containment_flag, hard_recovery_hold_flag
            real(dp) :: final_stationary_flag, anchor_interior_flag
            real(dp) :: domain_warning_flag, domain_ok_flag, domain_failure_flag

            if (measurement_enabled) then
                measurement_flag = 1.0_dp
            else
                measurement_flag = 0.0_dp
            end if

            if (this%state%has_bracket) then
                bracket_flag = 1.0_dp
            else
                bracket_flag = 0.0_dp
            end if

            if (capture_mode) then
                capture_flag = 1.0_dp
            else
                capture_flag = 0.0_dp
            end if

            if (emergency_mode) then
                emergency_flag = 1.0_dp
            else
                emergency_flag = 0.0_dp
            end if

            if (flame_detected) then
                flame_detected_flag = 1.0_dp
            else
                flame_detected_flag = 0.0_dp
            end if

            controller_armed_flag = merge(1.0_dp, 0.0_dp, this%state%controller_armed)
            ramp_settled_flag = merge(1.0_dp, 0.0_dp, ramp_settled)
            response_settled_flag = merge(1.0_dp, 0.0_dp, response_settled)
            flame_established_flag = merge(1.0_dp, 0.0_dp, this%state%flame_established)
            containment_flag = merge(1.0_dp, 0.0_dp, containment_mode)
            hard_recovery_hold_flag = merge(1.0_dp, 0.0_dp, hard_recovery_hold)
            final_stationary_flag = merge(1.0_dp, 0.0_dp, final_structure_stationary)
            anchor_interior_flag = merge(1.0_dp, 0.0_dp, anchor_velocity_interior)
            domain_warning_flag = merge(1.0_dp, 0.0_dp, domain_warning)
            domain_ok_flag = merge(1.0_dp, 0.0_dp, domain_ok)
            domain_failure_flag = merge(1.0_dp, 0.0_dp, this%state%domain_failure)

            write(this%state%flame_loc_unit,'(100E20.12)') &
                time, current_flame_location(1:dimensions), &
                flame_velocity_lsq, flame_velocity_filtered, diag_flame_velocity_lsq, &
                diag_flame_velocity_filtered, this%state%measurement_velocity_save, &
                control_velocity, position_error, this%state%front_reference_coord, &
                this%state%inlet_velocity_applied, this%state%inlet_velocity_target, &
                target_step_for_log, this%state%adaptive_gain, real(this%state%hist_count,dp), measurement_flag, bracket_flag, &
                capture_flag, emergency_flag, flame_detected_flag, &
                front_spread, heat_release_integral, heat_release_max, heat_release_valid_limit, H_max, Tgrad_max, &
                real(this%state%control_stage,dp), sl_displacement, linear_r2, linear_rms, split_slope_diff, &
                real(this%state%correction_counter,dp), real(this%state%stabilization_counter,dp), real(active_track_number,dp), &
                controller_armed_flag, ramp_settled_flag, response_settled_flag, &
                velocity_tolerance_on, velocity_tolerance_off, this%state%arm_metric_rel_trend, &
                flame_established_flag, real(arming_window_samples_active, dp), &
                real(this%state%hard_recovery_corrections, dp), &
                containment_flag, hard_recovery_hold_flag, &
                real(this%state%hard_recovery_no_progress_count, dp), &
                final_stationary_flag, real(this%state%structure_stationary_counter, dp), &
                anchor_interior_flag, &
                x_preheat, x_thermal_low, x_thermal_high, thermal_thickness, &
                x_reaction_left, x_reaction_right, reaction_thickness, &
                inlet_preheat_distance, outlet_reaction_distance, inlet_theta_max, &
                inlet_margin_ratio, outlet_margin_ratio, domain_warning_flag, domain_ok_flag, &
                domain_failure_flag, correction_free_time
        end subroutine write_tracking_line

    end subroutine solve

end module flame_stabilization_solver_class
