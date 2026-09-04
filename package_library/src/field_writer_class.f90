module field_writer_class

    use field_output_types

    implicit none

    private
    public :: field_writer

    type, abstract :: field_writer
    contains
        procedure(write_snapshot_interface), deferred :: write_snapshot
    end type field_writer

    abstract interface
        subroutine write_snapshot_interface(this, snapshot, output_folder, file_stem)
            import :: field_writer, field_snapshot
            class(field_writer), intent(inout) :: this
            type(field_snapshot), intent(in) :: snapshot
            character(len=*), intent(in) :: output_folder
            character(len=*), intent(in) :: file_stem
        end subroutine write_snapshot_interface
    end interface

end module field_writer_class
