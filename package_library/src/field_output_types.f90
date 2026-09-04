module field_output_types

    use kind_parameters
    use computational_domain_class
    use computational_mesh_class
    use boundary_conditions_class
    use field_pointers

    implicit none

    private
    public :: output_field_descriptor, field_snapshot
    public :: FIELD_LOCATION_CELL_CENTER

    integer, parameter :: FIELD_LOCATION_CELL_CENTER = 1

    type :: output_field_descriptor
        type(field_scalar_cons_pointer) :: values
        character(len=80) :: nrg_name = ''
        character(len=40) :: short_name = ''
        integer :: grid_location = FIELD_LOCATION_CELL_CENTER
    end type output_field_descriptor

    ! A format-independent view of one NRG field snapshot.  The arrays remain
    ! owned by the solver/data manager; writers only receive pointer wrappers.
    type :: field_snapshot
        type(computational_domain) :: domain
        type(computational_mesh_pointer) :: mesh
        type(boundary_conditions_pointer) :: boundaries
        type(output_field_descriptor), allocatable :: fields(:)

        real(dp) :: time_seconds = 0.0_dp
        real(dp) :: display_time = 0.0_dp
        integer :: output_index = 0
        logical :: include_ghost_cells = .false.
        character(len=10) :: time_units_abbreviation = 's'

        ! Common serialization metadata. Writers must use these values
        ! verbatim so Tecplot and CGNS expose the same dataset/zone identity.
        character(len=32) :: dataset_name = 'NRG'
        character(len=32) :: zone_name = ''
    end type field_snapshot

end module field_output_types
