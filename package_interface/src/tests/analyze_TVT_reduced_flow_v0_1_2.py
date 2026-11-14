#!/usr/bin/env python3
"""Analyze NRG Tecplot cell-centred flows from reduced active-matter experiments.

Reads ZONE/POINT files with arbitrary additional saved variables. Uses variable
names rather than fixed column indices. Grid geometry is reconstructed from the
NRG reduced interface; only active fluid cells contribute to diagnostics.

This tool cannot measure staggered wall-face fluxes or discrete Poisson residuals.
Those diagnostics must be instrumented in the Fortran solver.
"""
from __future__ import annotations

import argparse
import csv
import re
from pathlib import Path

import numpy as np

VARIABLE = re.compile(r'"([^"]+)"')
DIMENSIONS = re.compile(r'\bI\s*=\s*(\d+)\s*,\s*J\s*=\s*(\d+)\s*,\s*K\s*=\s*(\d+)')
TIME = re.compile(r'SOLUTIONTIME\s*=\s*([-+0-9.eEdD]+)')


def parse_tecplot(file_path: Path, input_time_unit: str = 'microseconds'):
    """Yield (physical time in milliseconds, shape, variables, numeric_data).

    NRG Tecplot SOLUTIONTIME is expressed in the data_save_c save_time_units,
    not necessarily seconds or milliseconds. The reduced TVT interface uses
    microseconds, so 500 in SOLUTIONTIME means 0.5 ms.
    """
    unit_to_ms = {
        'seconds': 1000.0,
        'milliseconds': 1.0,
        'microseconds': 1.0e-3,
        'nanoseconds': 1.0e-6,
    }
    time_factor = unit_to_ms[input_time_unit]
    with file_path.open('r', encoding='utf-8', errors='replace') as stream:
        vars_section = []
        names = None
        pending = None
        for line in stream:
            if line.lstrip().startswith('ZONE '):
                pending = line
                break
            if 'VARIABLES' in line.upper() or vars_section:
                vars_section.append(line)
        if pending is None:
            raise ValueError('No Tecplot ZONE header found')
        names = VARIABLE.findall(''.join(vars_section))
        if len(names) < 10:
            raise ValueError('Missing/invalid Tecplot VARIABLES declaration')

        while pending:
            if not pending.lstrip().startswith('ZONE '):
                raise ValueError(f'Unexpected header: {pending[:120]!r}')
            header = [pending]
            while True:
                line = stream.readline()
                if not line:
                    raise ValueError('Unexpected EOF in ZONE header')
                header.append(line)
                if line.lstrip().startswith('DT='):
                    break
            hdr = ''.join(header)
            mt = TIME.search(hdr)
            md = DIMENSIONS.search(hdr)
            if mt is None or md is None:
                raise ValueError(f'Unable to parse time and dimensions in: {hdr[:250]!r}')
            time_ms = float(mt.group(1).replace('D', 'E')) * time_factor
            shape = tuple(int(v) for v in md.groups())
            nrows = int(np.prod(shape))
            data = np.loadtxt((stream.readline() for _ in range(nrows)), dtype=float)
            if data.ndim != 2 or data.shape != (nrows, len(names)):
                raise ValueError(f'Expected {nrows}x{len(names)} columns, got {data.shape}')
            yield time_ms, shape, names, data
            pending = None
            for line in stream:
                if line.lstrip().startswith('ZONE '):
                    pending = line
                    break


def mask_for_grid(n: int, box_um: float, spherical_mask: bool):
    centers = (np.arange(n) + 0.5) * box_um/n
    xx, yy, zz = np.meshgrid(centers, centers, centers, indexing='ij')
    if spherical_mask:
        fluid = ((xx-box_um/2)**2 + (yy-box_um/2)**2 + (zz-box_um/2)**2) <= (box_um/2.4)**2
    else:
        fluid = np.ones((n,n,n), dtype=bool)
    # The NRG interface marks all outer cell planes as stationary walls.
    for axis in range(3):
        left = [slice(None)]*3
        right = [slice(None)]*3
        left[axis] = 0
        right[axis] = -1
        fluid[tuple(left)] = False
        fluid[tuple(right)] = False
    wall_adjacent = np.zeros_like(fluid)
    for axis in range(3):
        lower = [slice(None)]*3
        upper = [slice(None)]*3
        lower[axis] = slice(None,-1)
        upper[axis] = slice(1,None)
        lo, hi = tuple(lower), tuple(upper)
        wall_adjacent[lo] |= fluid[lo] & ~fluid[hi]
        wall_adjacent[hi] |= fluid[hi] & ~fluid[lo]
    return fluid, wall_adjacent


def field_index(names, field):
    try:
        return names.index(field)
    except ValueError as err:
        raise ValueError(f'Missing field {field!r}; available names: {names}') from err


def evaluate(file_path, n, box_um, spherical_mask, input_time_unit="microseconds"):
    fluid, wall = mask_for_grid(n, box_um, spherical_mask)
    dx_um = box_um/n
    volume = (dx_um * 1.0e-6)**3
    rows = []
    mass0 = None
    for t, shape, names, data in parse_tecplot(file_path, input_time_unit):
        if shape != (n,n,n):
            raise ValueError(f'File dimensions {shape} differ from requested {n}^3')
        indices = np.rint(data[:,0:3]*1e6/dx_um-0.5).astype(int)
        if (indices < 0).any() or (indices >= n).any():
            raise ValueError('Unexpected mesh coordinates (check --n/--box-um)')
        ids = tuple(indices[:,j] for j in range(3))
        selected = fluid[ids]
        next_wall = wall[ids]
        if not selected.any():
            raise ValueError('No active fluid cells found')
        density = data[selected,field_index(names,'rho')]
        temperature = data[selected,field_index(names,'T')]
        mu = data[selected,field_index(names,'nu')]
        velocity = data[:,[field_index(names, f'v({axis})') for axis in 'xyz']]
        speed = np.linalg.norm(velocity, axis=1)
        fspeed = speed[selected]
        max_idx = np.where(selected)[0][np.argmax(fspeed)]
        mass = float(density.sum()*volume)
        if mass0 is None:
            mass0 = mass
        energy = float((0.5*density*fspeed**2).sum()*volume)
        out = dict(time_ms=round(t,9),
                   fluid_cells=int(selected.sum()),
                   max_velocity_um_s=float(fspeed.max()*1e6),
                   rms_velocity_um_s=float(np.sqrt(np.mean(fspeed**2))*1e6),
                   max_velocity_is_wall_adjacent=(bool(next_wall[max_idx]) if fspeed.max() > 0.0 else None),
                   temperature_min_K=float(temperature.min()),
                   temperature_max_K=float(temperature.max()),
                   density_min_kg_m3=float(density.min()),
                   density_max_kg_m3=float(density.max()),
                   viscosity_min_Pa_s=float(mu.min()),
                   viscosity_max_Pa_s=float(mu.max()),
                   total_fluid_mass_kg=mass,
                   relative_mass_change=float((mass-mass0)/mass0),
                   kinetic_energy_J=energy)
        prefix = 'v_prod_p01'
        if all(f'{prefix}({axis})' in names for axis in 'xyz'):
            source = data[selected,:][:,[field_index(names, f'{prefix}({axis})') for axis in 'xyz']]
            out['particle_source_max_m_s2'] = float(np.linalg.norm(source,axis=1).max())
            out['particle_source_total_force_N'] = float(np.linalg.norm((density[:,None]*source).sum(axis=0)*volume))
        rows.append(out)
    if not rows:
        raise ValueError('No snapshots found')
    return rows


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('flow',type=Path)
    p.add_argument('--n',type=int,default=16)
    p.add_argument('--box-um',type=float,default=None)
    p.add_argument('--geometry',choices=('sphere','box'),default='sphere')
    p.add_argument('--csv',type=Path,default=None)
    p.add_argument('--input-time-unit', choices=('seconds','milliseconds','microseconds','nanoseconds'),
                   default='microseconds', help='Units used by Tecplot SOLUTIONTIME; TVT reduced default: microseconds')
    opts = p.parse_args()
    box = opts.box_um if opts.box_um is not None else opts.n*3.75
    rows = evaluate(opts.flow,opts.n,box,opts.geometry=='sphere',opts.input_time_unit)
    for row in rows:
        print(f"t={row['time_ms']:8.4f} ms |u|max={row['max_velocity_um_s']:11.3f} um/s "
              f"T=[{row['temperature_min_K']:.2f},{row['temperature_max_K']:.2f}] K "
              f"rho=[{row['density_min_kg_m3']:.2f},{row['density_max_kg_m3']:.2f}] "
              f"dM/M0={row['relative_mass_change']:+.3e} "
              f"K={row['kinetic_energy_J']:.3e} J")
    if opts.csv:
        fields=list(dict.fromkeys(k for row in rows for k in row))
        with opts.csv.open('w',newline='',encoding='utf-8') as fout:
            writer=csv.DictWriter(fout,fieldnames=fields)
            writer.writeheader();writer.writerows(rows)
        print('Wrote',opts.csv)

if __name__ == '__main__':
    main()
