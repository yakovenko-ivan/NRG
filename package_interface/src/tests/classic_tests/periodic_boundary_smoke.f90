program periodic_boundary_smoke
    use kind_parameters
    use computational_domain_class
    use boundary_conditions_class

    implicit none

    type(computational_domain) :: domain
    type(boundary_conditions) :: boundaries
    integer, dimension(3,2) :: utter

    domain = computational_domain_c( &
        dimensions        = 2, &
        cells_number      = (/8,6,1/), &
        coordinate_system = 'cartesian', &
        lengths           = reshape((/0.0_dp,0.0_dp,0.0_dp,1.0_dp,1.0_dp,1.0_dp/),(/3,2/)), &
        axis_names        = (/'x','y','z'/), &
        periodic          = (/.true.,.false.,.false./))

    boundaries = boundary_conditions_c(domain, number_of_boundary_types=1, default_boundary=1)
    call boundaries%create_boundary_type( &
        type_name='wall', slip=.true., conductive=.false., &
        wall_temperature=0.0_dp, wall_conductivity_ratio=0.0_dp, priority=1)

    utter = domain%get_local_utter_cells_bounds()

    ! x is periodic: away from y-wall intersections its ghost markers are fluid.
    if (any(boundaries%bc_markers(utter(1,1),1:6,1) /= 0)) &
        error stop 'periodic smoke: x-low periodic markers are not zero'
    if (any(boundaries%bc_markers(utter(1,2),1:6,1) /= 0)) &
        error stop 'periodic smoke: x-high periodic markers are not zero'

    ! y is a real physical boundary and therefore retains the wall marker.
    if (any(boundaries%bc_markers(:,utter(2,1),1) /= 1)) &
        error stop 'periodic smoke: y-low wall markers are incorrect'
    if (any(boundaries%bc_markers(:,utter(2,2),1) /= 1)) &
        error stop 'periodic smoke: y-high wall markers are incorrect'

    print *, 'periodic boundary smoke test: PASS'
end program periodic_boundary_smoke
