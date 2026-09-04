#include <stdint.h>
#include <stdio.h>
#include <string.h>

#include <cgnslib.h>

static int nrg_cgns_report(int status, const char *operation)
{
    if (status != CG_OK) {
        const char *message = cg_get_error();
        fprintf(stderr, "NRG CGNS error in %s: %s\n",
                operation, message != NULL ? message : "unknown CGNS error");
    }
    return status;
}

int nrg_cgns_begin(const char *filename,
                   const char *dataset_name,
                   const char *zone_name,
                   int cell_dim,
                   int phys_dim,
                   const int64_t cell_counts[3],
                   double time_seconds,
                   int output_index,
                   const char *coordinate_system,
                   int *fn,
                   int *base,
                   int *zone,
                   int *solution)
{
    int status;
    int d;
    int file_type = CG_FILE_NONE;
    cgsize_t zone_size[9] = {0};
    cgsize_t one = 1;
    double time_value = time_seconds;
    char index_text[64];

    if (cell_dim < 1 || cell_dim > 3 || phys_dim < 1 || phys_dim > 3) {
        fprintf(stderr, "NRG CGNS error: unsupported dimensions cell=%d physical=%d\n",
                cell_dim, phys_dim);
        return CG_ERROR;
    }

    if (dataset_name == NULL || dataset_name[0] == '\0' ||
        zone_name == NULL || zone_name[0] == '\0') {
        fprintf(stderr, "NRG CGNS error: dataset and zone names must not be empty\n");
        return CG_ERROR;
    }
    if (strlen(dataset_name) > 32 || strlen(zone_name) > 32) {
        fprintf(stderr, "NRG CGNS error: dataset/zone names exceed the 32-character CGNS limit\n");
        return CG_ERROR;
    }

    status = cg_set_file_type(CG_FILE_HDF5);
    if (nrg_cgns_report(status, "cg_set_file_type(CG_FILE_HDF5)") != CG_OK)
        return status;

    status = cg_open(filename, CG_MODE_WRITE, fn);
    if (nrg_cgns_report(status, "cg_open") != CG_OK)
        return status;

    status = cg_get_file_type(*fn, &file_type);
    if (nrg_cgns_report(status, "cg_get_file_type") != CG_OK)
        goto fail;
    if (file_type != CG_FILE_HDF5) {
        fprintf(stderr, "NRG CGNS error: the opened CGNS file is not HDF5-backed. "
                        "Check the CGNS build and CGNS_FILETYPE environment variable.\n");
        status = CG_ERROR;
        goto fail;
    }

    status = cg_base_write(*fn, dataset_name, cell_dim, phys_dim, base);
    if (nrg_cgns_report(status, "cg_base_write") != CG_OK)
        goto fail;

    status = cg_simulation_type_write(*fn, *base, TimeAccurate);
    if (nrg_cgns_report(status, "cg_simulation_type_write") != CG_OK)
        goto fail;

    for (d = 0; d < cell_dim; ++d) {
        if (cell_counts[d] < 1) {
            fprintf(stderr, "NRG CGNS error: non-positive cell count in direction %d\n", d + 1);
            status = CG_ERROR;
            goto fail;
        }
        zone_size[d] = (cgsize_t)(cell_counts[d] + 1);
        zone_size[cell_dim + d] = (cgsize_t)cell_counts[d];
        zone_size[2 * cell_dim + d] = 0;
    }

    status = cg_zone_write(*fn, *base, zone_name, zone_size, Structured, zone);
    if (nrg_cgns_report(status, "cg_zone_write") != CG_OK)
        goto fail;

    status = cg_sol_write(*fn, *base, *zone, "FlowSolutionCellCenter", CellCenter, solution);
    if (nrg_cgns_report(status, "cg_sol_write") != CG_OK)
        goto fail;

    status = cg_biter_write(*fn, *base, "TimeIterValues", 1);
    if (nrg_cgns_report(status, "cg_biter_write") != CG_OK)
        goto fail;

    status = cg_goto(*fn, *base, "BaseIterativeData_t", 1, "end");
    if (nrg_cgns_report(status, "cg_goto(BaseIterativeData_t)") != CG_OK)
        goto fail;

    status = cg_array_write("TimeValues", RealDouble, 1, &one, &time_value);
    if (nrg_cgns_report(status, "cg_array_write(TimeValues)") != CG_OK)
        goto fail;

    status = cg_goto(*fn, *base, "end");
    if (nrg_cgns_report(status, "cg_goto(CGNSBase_t)") != CG_OK)
        goto fail;

    status = cg_descriptor_write("NRGOutputSchema", "NRG-CGNS-1.0");
    if (nrg_cgns_report(status, "cg_descriptor_write(NRGOutputSchema)") != CG_OK)
        goto fail;

    status = cg_descriptor_write("NRGCoordinateSystem", coordinate_system);
    if (nrg_cgns_report(status, "cg_descriptor_write(NRGCoordinateSystem)") != CG_OK)
        goto fail;

    status = cg_descriptor_write("NRGUnitSystem", "SI");
    if (nrg_cgns_report(status, "cg_descriptor_write(NRGUnitSystem)") != CG_OK)
        goto fail;

    snprintf(index_text, sizeof(index_text), "%d", output_index);
    status = cg_descriptor_write("NRGOutputIndex", index_text);
    if (nrg_cgns_report(status, "cg_descriptor_write(NRGOutputIndex)") != CG_OK)
        goto fail;

    return CG_OK;

fail:
    cg_close(*fn);
    return status;
}

int nrg_cgns_write_coordinate(int fn, int base, int zone,
                              const char *name, const double *data)
{
    int coord = 0;
    int status = cg_coord_write(fn, base, zone, RealDouble, name, data, &coord);
    return nrg_cgns_report(status, "cg_coord_write");
}

int nrg_cgns_set_solution_rind(int fn, int base, int zone, int solution,
                               int index_dim, const int rind_data[6])
{
    int status;
    int rind[6] = {0, 0, 0, 0, 0, 0};
    int i;

    for (i = 0; i < 2 * index_dim; ++i)
        rind[i] = rind_data[i];

    status = cg_goto(fn, base, "Zone_t", zone,
                     "FlowSolution_t", solution, "end");
    if (nrg_cgns_report(status, "cg_goto(FlowSolution_t)") != CG_OK)
        return status;

    status = cg_rind_write(rind);
    return nrg_cgns_report(status, "cg_rind_write");
}

int nrg_cgns_write_field_double(int fn, int base, int zone, int solution,
                                const char *name, const double *data,
                                const char *nrg_long_name,
                                const char *nrg_short_name,
                                const char *units)
{
    int field = 0;
    int status = cg_field_write(fn, base, zone, solution,
                                RealDouble, name, data, &field);
    if (nrg_cgns_report(status, "cg_field_write(RealDouble)") != CG_OK)
        return status;

    status = cg_goto(fn, base, "Zone_t", zone,
                     "FlowSolution_t", solution,
                     "DataArray_t", field, "end");
    if (nrg_cgns_report(status, "cg_goto(DataArray_t)") != CG_OK)
        return status;

    status = cg_descriptor_write("NRGLongName", nrg_long_name);
    if (nrg_cgns_report(status, "cg_descriptor_write(NRGLongName)") != CG_OK)
        return status;

    status = cg_descriptor_write("NRGShortName", nrg_short_name);
    if (nrg_cgns_report(status, "cg_descriptor_write(NRGShortName)") != CG_OK)
        return status;

    if (units != NULL && units[0] != '\0') {
        status = cg_descriptor_write("NRGUnits", units);
        if (nrg_cgns_report(status, "cg_descriptor_write(NRGUnits)") != CG_OK)
            return status;
    }

    return CG_OK;
}

int nrg_cgns_write_field_int(int fn, int base, int zone, int solution,
                             const char *name, const int *data)
{
    int field = 0;
    int status = cg_field_write(fn, base, zone, solution,
                                Integer, name, data, &field);
    return nrg_cgns_report(status, "cg_field_write(Integer)");
}

int nrg_cgns_close(int fn)
{
    int status = cg_close(fn);
    return nrg_cgns_report(status, "cg_close");
}
