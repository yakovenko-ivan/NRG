#!/usr/bin/env python3
"""Apply NRG periodic-boundary infrastructure, stage 1.

Target base: feature/flame-anchoring-2d @ 82f470052860e7d0985ce9640a23f4c81fad0fd4

Changes:
  * computational_domain gains periodic(3), constructor/file I/O/getter/logging.
  * periodicity is passed to MPI_CART_CREATE.
  * multi-rank periodic execution is rejected until seam exchange is implemented.
  * boundary_conditions leaves periodic outer ghost planes as fluid marker 0.
  * existing conservative MPI exchange entry points also perform local periodic
    halo wrapping when a periodic direction has a single rank (including a
    non-MPI/OpenMP build).
  * adds a small package-interface smoke test.

The script is intentionally defensive: every source edit must match exactly once.
"""
from __future__ import annotations

from pathlib import Path
import re
import sys

ROOT = Path.cwd()


def read(rel: str) -> str:
    p = ROOT / rel
    if not p.exists():
        raise SystemExit(f"missing expected file: {p}")
    return p.read_text(encoding="utf-8")


def write(rel: str, text: str) -> None:
    p = ROOT / rel
    p.write_text(text, encoding="utf-8")
    print(f"updated {rel}")


def replace_once(text: str, old: str, new: str, label: str) -> str:
    n = text.count(old)
    if n != 1:
        raise SystemExit(f"{label}: expected exactly one match, found {n}")
    return text.replace(old, new, 1)


def regex_once(text: str, pattern: str, repl: str, label: str, flags: int = 0) -> str:
    out, n = re.subn(pattern, repl, text, count=1, flags=flags)
    if n != 1:
        raise SystemExit(f"{label}: expected exactly one regex match, found {n}")
    return out


def insert_before_subroutine_end(text: str, name: str, insertion: str) -> str:
    pat = rf"(\n\s*subroutine\s+{re.escape(name)}\b.*?)(\n\s*end\s+subroutine(?:\s+{re.escape(name)})?\s*\n)"
    m = re.search(pat, text, flags=re.IGNORECASE | re.DOTALL)
    if not m:
        raise SystemExit(f"could not locate subroutine {name}")
    body = m.group(1)
    if insertion.strip() in body:
        raise SystemExit(f"subroutine {name}: insertion already present")
    return text[:m.start()] + body + insertion + m.group(2) + text[m.end():]


# ---------------------------------------------------------------------------
# computational_domain_class.f90
# ---------------------------------------------------------------------------
rel = "package_library/src/computational_domain_class.f90"
s = read(rel)

s = replace_once(
    s,
    "\t\tcharacter(len=20)\t\t\t\t\t\t\t\t\t:: coordinate_system\t! Coordinate system (cartesian/cylindrical/spherical)\n",
    "\t\tcharacter(len=20)\t\t\t\t\t\t\t\t\t:: coordinate_system\t! Coordinate system (cartesian/cylindrical/spherical)\n"
    "\t\tlogical\t\t\t\t,dimension(3)\t\t\t\t\t:: periodic = .false.\t! Periodic topology by spatial axis\n",
    "domain: periodic member",
)

s = replace_once(
    s,
    "\t\tprocedure\t\t\t\t:: get_coordinate_system_name\n",
    "\t\tprocedure\t\t\t\t:: get_coordinate_system_name\n"
    "\t\tprocedure\t\t\t\t:: get_periodic_directions\n",
    "domain: periodic getter binding",
)

s = replace_once(
    s,
    "\ttype(computational_domain)\tfunction constructor(dimensions,cells_number,coordinate_system,lengths,axis_names)\n",
    "\ttype(computational_domain)\tfunction constructor(dimensions,cells_number,coordinate_system,lengths,axis_names,periodic)\n",
    "domain: constructor signature",
)

s = replace_once(
    s,
    "        character(len=*)    ,dimension(3)   ,intent(in) :: axis_names\n\n\t\tinteger\t:: io_unit\n",
    "        character(len=*)    ,dimension(3)   ,intent(in) :: axis_names\n"
    "        logical             ,dimension(3)   ,intent(in), optional :: periodic\n\n"
    "\t\tinteger\t:: io_unit\n",
    "domain: constructor periodic argument",
)

s = replace_once(
    s,
    "\t\tcall constructor%set_properties(dimensions,cells_number,lengths,coordinate_system,axis_names)\n",
    "\t\tcall constructor%set_properties(dimensions,cells_number,lengths,coordinate_system,axis_names,periodic)\n",
    "domain: constructor set_properties",
)

# read_properties block
s = replace_once(
    s,
    "        character(len=5)    ,dimension(3)   :: axis_names\n\n\t\tnamelist /domain_properties_1/ dimensions, cells_number, coordinate_system\n\t\tnamelist /domain_properties_2/ lengths, axis_names\n\t\t\n\t\tread(unit = domain_data_unit, nml = domain_properties_1)\n",
    "        character(len=5)    ,dimension(3)   :: axis_names\n"
    "        logical             ,dimension(3)   :: periodic\n\n"
    "\t\tnamelist /domain_properties_1/ dimensions, cells_number, coordinate_system, periodic\n"
    "\t\tnamelist /domain_properties_2/ lengths, axis_names\n\t\t\n"
    "\t\t! Old domain_data.inf files do not contain PERIODIC; retain non-periodic behavior.\n"
    "\t\tperiodic = .false.\n\t\t\n"
    "\t\tread(unit = domain_data_unit, nml = domain_properties_1)\n",
    "domain: read periodic namelist",
)

s = replace_once(
    s,
    "\t\tcall this%set_properties(dimensions,cells_number,lengths,coordinate_system,axis_names)\n",
    "\t\tcall this%set_properties(dimensions,cells_number,lengths,coordinate_system,axis_names,periodic)\n",
    "domain: read set_properties",
)

# write_properties block (remaining matching axis_names/namelist sequence)
s = replace_once(
    s,
    "        character(len=5)    ,dimension(3)   :: axis_names\n\n\t\tnamelist /domain_properties_1/ dimensions, cells_number, coordinate_system\n\t\tnamelist /domain_properties_2/ lengths, axis_names\n",
    "        character(len=5)    ,dimension(3)   :: axis_names\n"
    "        logical             ,dimension(3)   :: periodic\n\n"
    "\t\tnamelist /domain_properties_1/ dimensions, cells_number, coordinate_system, periodic\n"
    "\t\tnamelist /domain_properties_2/ lengths, axis_names\n",
    "domain: write periodic namelist",
)

s = replace_once(
    s,
    "\t\tcoordinate_system\t= this%coordinate_system\n\t\t\n\t\taxis_names\t= this%axis_names\n",
    "\t\tcoordinate_system\t= this%coordinate_system\n"
    "\t\tperiodic            = this%periodic\n\t\t\n"
    "\t\taxis_names\t= this%axis_names\n",
    "domain: write periodic value",
)

s = replace_once(
    s,
    "\tsubroutine set_properties(this,dimensions,cells_number,lengths,coordinate_system,axis_names)\n",
    "\tsubroutine set_properties(this,dimensions,cells_number,lengths,coordinate_system,axis_names,periodic)\n",
    "domain: set_properties signature",
)

s = replace_once(
    s,
    "        character(len=*)    ,dimension(3)   ,intent(in)     :: axis_names\n\t\t\n\t\tinteger\t:: dim\n",
    "        character(len=*)    ,dimension(3)   ,intent(in)     :: axis_names\n"
    "        logical             ,dimension(3)   ,intent(in), optional :: periodic\n\t\t\n"
    "\t\tinteger\t:: dim\n",
    "domain: set_properties periodic argument",
)

s = replace_once(
    s,
    "\t\tthis%coordinate_system\t\t= coordinate_system\n\t\tthis%axis_names\t\t\t\t= axis_names\n\n\t\tthis%cells_number(dimensions+1:3) = 1\n",
    "\t\tthis%coordinate_system\t\t= coordinate_system\n"
    "\t\tthis%axis_names\t\t\t\t= axis_names\n"
    "\t\tthis%periodic               = .false.\n"
    "\t\tif (present(periodic)) this%periodic(1:dimensions) = periodic(1:dimensions)\n\n"
    "\t\tthis%cells_number(dimensions+1:3) = 1\n",
    "domain: assign periodic",
)

s = replace_once(
    s,
    "\t\tthis%faces_number(dimensions+1:3) = 1\n\n\t\tthis%mpi_communicator_size = 1\n",
    "\t\tthis%faces_number(dimensions+1:3) = 1\n\n"
    "\t\tif (any(this%periodic) .and. trim(this%coordinate_system) /= 'cartesian') then\n"
    "\t\t\terror stop 'Periodic boundaries are currently supported only for Cartesian domains.'\n"
    "\t\tend if\n\n"
    "\t\tthis%mpi_communicator_size = 1\n",
    "domain: cartesian periodic validation",
)

s = replace_once(
    s,
    "\t\twrite(log_unit,'(A,3A)')\t' Domain axis names         : ',\tthis%axis_names\t\n",
    "\t\twrite(log_unit,'(A,3A)')\t' Domain axis names         : ',\tthis%axis_names\t\n"
    "\t\twrite(log_unit,'(A,3L2)')\t' Domain periodic axes      : ',\tthis%periodic\n",
    "domain: periodic log",
)

s = replace_once(
    s,
    "#ifdef mpi\n    \tcall MPI_COMM_SIZE(MPI_COMM_WORLD, this%mpi_communicator_size, error)\n\t\tcall MPI_COMM_RANK(MPI_COMM_WORLD, this%processor_rank, error)\n#endif\n",
    "#ifdef mpi\n"
    "    \tcall MPI_COMM_SIZE(MPI_COMM_WORLD, this%mpi_communicator_size, error)\n"
    "\t\tcall MPI_COMM_RANK(MPI_COMM_WORLD, this%processor_rank, error)\n"
    "\t\tif (any(this%periodic) .and. this%mpi_communicator_size > 1) then\n"
    "\t\t\terror stop 'Multi-rank periodic halo exchange is not implemented yet; use one MPI rank/OpenMP.'\n"
    "\t\tend if\n"
    "#endif\n",
    "domain: multi-rank periodic guard",
)

s = replace_once(
    s,
    "\t\tis_periodic = .false.\n",
    "\t\tis_periodic = this%periodic\n",
    "domain: MPI Cartesian periodicity",
)

getter_anchor = (
    "\tpure function get_coordinate_system_name(this)\n"
    "\t\tclass(computational_domain)\t\t\t,intent(in)\t\t:: this\n"
    "\t\tcharacter(len=20)\t:: get_coordinate_system_name\n\t\t\n"
    "\t\tget_coordinate_system_name = this%coordinate_system\n"
    "\tend function\t\n"
)
getter_new = getter_anchor + (
    "\n\tpure function get_periodic_directions(this)\n"
    "\t\tclass(computational_domain)\t,intent(in)\t:: this\n"
    "\t\tlogical, dimension(3)\t\t\t\t\t:: get_periodic_directions\n\n"
    "\t\tget_periodic_directions = this%periodic\n"
    "\tend function\n"
)
s = replace_once(s, getter_anchor, getter_new, "domain: periodic getter implementation")
write(rel, s)


# ---------------------------------------------------------------------------
# boundary_conditions_class.f90
# ---------------------------------------------------------------------------
rel = "package_library/src/boundary_conditions_class.f90"
s = read(rel)

s = replace_once(
    s,
    "\t\tinteger\t,dimension(3)\t:: processor_grid_coord \n\n\t\tinteger\t:: bound_number ,bound\n",
    "\t\tinteger\t,dimension(3)\t:: processor_grid_coord \n"
    "\t\tlogical\t,dimension(3)\t:: periodic\n\n"
    "\t\tinteger\t:: bound_number ,bound\n",
    "bc: periodic local",
)

s = replace_once(
    s,
    "\t\tprocessor_number\t\t= domain%get_processor_number()\n"
    "\t\tprocessor_grid_coord\t= domain%get_processor_grid_coord()\n",
    "\t\tprocessor_number\t\t= domain%get_processor_number()\n"
    "\t\tprocessor_grid_coord\t= domain%get_processor_grid_coord()\n"
    "\t\tperiodic                = domain%get_periodic_directions()\n",
    "bc: read periodic topology",
)

old = """\t\tif (processor_grid_coord(1) == 0) \t\t\t\t\t\tthis%bc_markers(allocation_bounds(1,1),:,:)\t= default_boundary
\t\tif (processor_grid_coord(1) == processor_number(1)-1) \tthis%bc_markers(allocation_bounds(1,2),:,:)\t= default_boundary

\t\tif(dimensions >= 2) then
\t\t\tif (processor_grid_coord(2) == 0) \t\t\t\t\t\tthis%bc_markers(:,allocation_bounds(2,1),:)\t= default_boundary
\t\t\tif (processor_grid_coord(2) == processor_number(2)-1) \tthis%bc_markers(:,allocation_bounds(2,2),:)\t= default_boundary

\t\t\tif(dimensions == 3) then
\t\t\t\tif (processor_grid_coord(3) == 0) \t\t\t\t\t\tthis%bc_markers(:,:,allocation_bounds(3,1))\t= default_boundary
\t\t\t\tif (processor_grid_coord(3) == processor_number(3)-1) \tthis%bc_markers(:,:,allocation_bounds(3,2))\t= default_boundary\t\t\t\t
\t\t\tend if
\t\tend if
"""
new = """\t\t! Periodic ghost planes remain marker 0: they are fluid connectivity,
\t\t! not physical boundary conditions.  Intersections with a non-periodic
\t\t! boundary retain that physical boundary's marker.
\t\tif (.not. periodic(1)) then
\t\t\tif (processor_grid_coord(1) == 0) \t\t\t\t\t\tthis%bc_markers(allocation_bounds(1,1),:,:)\t= default_boundary
\t\t\tif (processor_grid_coord(1) == processor_number(1)-1) \tthis%bc_markers(allocation_bounds(1,2),:,:)\t= default_boundary
\t\tend if

\t\tif(dimensions >= 2) then
\t\t\tif (.not. periodic(2)) then
\t\t\t\tif (processor_grid_coord(2) == 0) \t\t\t\t\t\tthis%bc_markers(:,allocation_bounds(2,1),:)\t= default_boundary
\t\t\t\tif (processor_grid_coord(2) == processor_number(2)-1) \tthis%bc_markers(:,allocation_bounds(2,2),:)\t= default_boundary
\t\t\tend if

\t\t\tif(dimensions == 3) then
\t\t\t\tif (.not. periodic(3)) then
\t\t\t\t\tif (processor_grid_coord(3) == 0) \t\t\t\t\t\tthis%bc_markers(:,:,allocation_bounds(3,1))\t= default_boundary
\t\t\t\t\tif (processor_grid_coord(3) == processor_number(3)-1) \tthis%bc_markers(:,:,allocation_bounds(3,2))\t= default_boundary\t\t\t\t
\t\t\t\tend if
\t\t\tend if
\t\tend if
"""
s = replace_once(s, old, new, "bc: periodic marker initialization")
write(rel, s)


# ---------------------------------------------------------------------------
# mpi_communications_class.f90 -- local, one-rank periodic conservative halos
# ---------------------------------------------------------------------------
rel = "package_library/src/mpi_communications_class.f90"
s = read(rel)

s = insert_before_subroutine_end(
    s,
    "exchange_conservative_scalar_field",
    "\n        ! Also covers serial/OpenMP and one-rank periodic axes.\n"
    "        call apply_local_periodic_cons_halo(this%domain, scal_ptr%cells)\n",
)

s = insert_before_subroutine_end(
    s,
    "exchange_boundary_conditions_markers",
    "\n        call apply_local_periodic_marker_halo(this%domain, bc_ptr%bc_markers)\n",
)

s = insert_before_subroutine_end(
    s,
    "exchange_conservative_vector_field",
    "\n        vector_projections_number = vect_ptr%get_projections_number()\n"
    "        do dim = 1, vector_projections_number\n"
    "            call apply_local_periodic_cons_halo(this%domain, vect_ptr%pr(dim)%cells)\n"
    "        end do\n",
)

s = insert_before_subroutine_end(
    s,
    "exchange_conservative_tensor_field",
    "\n        tensor_projections_number = tens_ptr%get_projections_number()\n"
    "        do dim4 = 1, tensor_projections_number(1)\n"
    "        do dim5 = 1, tensor_projections_number(2)\n"
    "            call apply_local_periodic_cons_halo(this%domain, tens_ptr%pr(dim4,dim5)%cells)\n"
    "        end do\n"
    "        end do\n",
)

helpers = r'''

    !--------------------------------------------------------------------------
    ! One-rank periodic halo helpers.  Assumed-shape dummy arrays are indexed
    ! relative to their local extent, so this remains correct for NRG fields
    ! whose allocated lower bound is zero in active dimensions.
    !--------------------------------------------------------------------------
    subroutine apply_local_periodic_cons_halo(domain, field)
        type(computational_domain), intent(in) :: domain
        real(dp), dimension(:,:,:), intent(inout) :: field
        logical, dimension(3) :: periodic
        integer, dimension(3) :: processor_number

        periodic = domain%get_periodic_directions()
        processor_number = domain%get_processor_number()

        if (periodic(1) .and. processor_number(1) == 1 .and. size(field,1) > 2) then
            field(1,:,:) = field(size(field,1)-1,:,:)
            field(size(field,1),:,:) = field(2,:,:)
        end if
        if (periodic(2) .and. processor_number(2) == 1 .and. size(field,2) > 2) then
            field(:,1,:) = field(:,size(field,2)-1,:)
            field(:,size(field,2),:) = field(:,2,:)
        end if
        if (periodic(3) .and. processor_number(3) == 1 .and. size(field,3) > 2) then
            field(:,:,1) = field(:,:,size(field,3)-1)
            field(:,:,size(field,3)) = field(:,:,2)
        end if
    end subroutine apply_local_periodic_cons_halo

    subroutine apply_local_periodic_marker_halo(domain, markers)
        type(computational_domain), intent(in) :: domain
        integer(i1), dimension(:,:,:), intent(inout) :: markers
        logical, dimension(3) :: periodic
        integer, dimension(3) :: processor_number

        periodic = domain%get_periodic_directions()
        processor_number = domain%get_processor_number()

        if (periodic(1) .and. processor_number(1) == 1 .and. size(markers,1) > 2) then
            markers(1,:,:) = markers(size(markers,1)-1,:,:)
            markers(size(markers,1),:,:) = markers(2,:,:)
        end if
        if (periodic(2) .and. processor_number(2) == 1 .and. size(markers,2) > 2) then
            markers(:,1,:) = markers(:,size(markers,2)-1,:)
            markers(:,size(markers,2),:) = markers(:,2,:)
        end if
        if (periodic(3) .and. processor_number(3) == 1 .and. size(markers,3) > 2) then
            markers(:,:,1) = markers(:,:,size(markers,3)-1)
            markers(:,:,size(markers,3)) = markers(:,:,2)
        end if
    end subroutine apply_local_periodic_marker_halo
'''

s = replace_once(
    s,
    "\nend module mpi_communications_class",
    helpers + "\nend module mpi_communications_class",
    "mpi: periodic halo helpers",
)
write(rel, s)


# ---------------------------------------------------------------------------
# Smoke-test package interface
# ---------------------------------------------------------------------------
smoke_rel = "package_interface/src/tests/classic_tests/periodic_boundary_smoke.f90"
smoke = r'''program periodic_boundary_smoke
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
'''
smoke_path = ROOT / smoke_rel
if smoke_path.exists():
    raise SystemExit(f"refusing to overwrite existing {smoke_rel}")
smoke_path.write_text(smoke, encoding="utf-8")
print(f"created {smoke_rel}")

print("\nStage 1 periodic-boundary infrastructure applied successfully.")
