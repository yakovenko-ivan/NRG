#!/usr/bin/env python3
"""Audit Tecplot ASCII POINT-format NRG TVT 3D active-matter flow snapshots.

Usage: python tvt_flow_audit.py flow.dat [--particles lagrangian_particles_phase_01.dat]

Outputs CSV and two PNG figures alongside flow.dat. Requires numpy, matplotlib.
Optional particle interpolation additionally needs scipy.
Assumes this TVT interface's centered R=50 um fluid sphere inside 120 um box.
"""
from __future__ import annotations
import argparse
import csv
import re
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt


def load_tecplot_zones(path):
    """Read a compact multi-ZONE NRG Tecplot file using POINT records."""
    with path.open('r', encoding='ascii', errors='replace') as f:
        while (line := f.readline()):
            if not line.strip().startswith('ZONE '):
                continue
            time_line = f.readline()
            mesh_line = f.readline()
            data_pack = f.readline()
            dtype_line = f.readline()
            assert 'DATAPACKING=POINT' in data_pack, data_pack
            tm = re.search(r'SOLUTIONTIME\s*=\s*([\d.eE+\-]+)', time_line)
            dims = re.search(r'I\s*=\s*(\d+)\s*,\s*J\s*=\s*(\d+)\s*,\s*K\s*=\s*(\d+)', mesh_line)
            if not (tm and dims):
                raise ValueError(f'Unexpected Tecplot ZONE header: {line!r}, {time_line!r}, {mesh_line!r}')
            nx,ny,nz = [int(d) for d in dims.groups()]
            count=nx*ny*nz
            data = np.loadtxt((f.readline() for _ in range(count)), dtype=float)
            if data.shape!=(count,14):
                raise ValueError(f'Expected 14 NRG fields, received {data.shape}')
            yield float(tm.group(1)),data,(nx,ny,nz)


def audit(filename, particles_file=None):
    filename=Path(filename).resolve()
    records=[]
    profile=[]
    last=None
    csv_path=filename.parent/'tvt_flow_audit.csv'
    hist_edges=np.arange(0,52,2)
    particles = np.loadtxt(particles_file,skiprows=1) if particles_file else None
    for time_ms, data, (nx,ny,nz) in load_tecplot_zones(filename):
        xyz=data[:,:3]
        r_um=1e6*np.linalg.norm(xyz-60e-6,axis=1)
        fluid = r_um<50
        core = r_um<40
        shell = (r_um>=48)&fluid
        vel = np.linalg.norm(data[:,7:10],axis=1)*1e6
        sources=np.linalg.norm(data[:,11:14],axis=1)
        row={'time_ms':time_ms,
             'fluid_cells':int(fluid.sum()),
             'core_mean_fluid_speed_um_s':float(vel[core].mean()),
             'core_max_fluid_speed_um_s':float(vel[core].max()),
             'boundary_48_50_mean_fluid_speed_um_s':float(vel[shell].mean()),
             'boundary_48_50_max_fluid_speed_um_s':float(vel[shell].max()),
             'fluid_max_speed_um_s':float(vel[fluid].max()),
             'fluid_temperature_min_K':float(data[fluid,5].min()),
             'fluid_temperature_max_K':float(data[fluid,5].max()),
             'fluid_density_min_kg_m3':float(data[fluid,6].min()),
             'fluid_density_max_kg_m3':float(data[fluid,6].max()),
             'core_viscosity_mean_Pa_s':float(data[core,10].mean()),
             'fluid_viscosity_min_Pa_s':float(data[fluid,10].min()),
             'fluid_viscosity_max_Pa_s':float(data[fluid,10].max()),
             'fluid_static_pressure_MPa':float(np.median(data[fluid,3])/1e6),
             'fluid_max_source_acceleration_m_s2':float(sources[fluid].max()),
             'nonzero_source_cells':int(np.count_nonzero(sources>1e-6))}
        if particles is not None:
            selected=particles[np.abs(particles[:,0]-time_ms/1000.)<5e-5]
            if len(selected):
                row['mean_particle_speed_um_s']=float((np.linalg.norm(selected[:,5:8],axis=1)*1e6).mean())
                row['max_particle_speed_um_s']=float((np.linalg.norm(selected[:,5:8],axis=1)*1e6).max())
                try:
                    from scipy.ndimage import map_coordinates
                    fields=data.reshape(nz,ny,nx,14)
                    dx=120e-6/nx
                    coord=(selected[:,2:5]/dx-.5)[:,::-1].T
                    u=np.column_stack([map_coordinates(fields[:,:,:,7+j],coord,order=1,mode='nearest') for j in range(3)])
                    row['mean_interpolated_fluid_at_particles_um_s']=float((np.linalg.norm(u,axis=1)*1e6).mean())
                    row['mean_particle_vs_centered_velocity_slip_um_s']=float((np.linalg.norm(selected[:,5:8]-u,axis=1)*1e6).mean())
                except ImportError:
                    pass
        records.append(row)
        for lo,hi in zip(hist_edges[:-1],hist_edges[1:]):
            sel=(r_um>=lo)&(r_um<hi)
            profile.append((time_ms,lo,hi,float(vel[sel].mean()) if sel.any() else np.nan))
        last=(time_ms,data,(nx,ny,nz))
    if not records:raise ValueError('No Tecplot zones detected')
    with csv_path.open('w',newline='',encoding='utf-8') as f:
        fieldnames=list(dict.fromkeys(k for x in records for k in x))
        writer=csv.DictWriter(f,fieldnames=fieldnames);writer.writeheader();writer.writerows(records)

    fig,ax=plt.subplots(figsize=(8.1,4.9),layout='constrained')
    for ms in sorted({r[0] for r in profile}):
        subset=np.array([r for r in profile if r[0]==ms])
        ax.plot((subset[:,1]+subset[:,2])/2,subset[:,3],marker='.',label=f'{ms:g} ms')
    ax.set_yscale('log')
    ax.set(xlabel='Distance from sphere center (µm)',ylabel='Shell-averaged |u| (µm/s)',title='Fluid acceleration is concentrated near the immersed boundary')
    ax.grid(True,which='both',alpha=0.3);ax.legend(ncol=2)
    radial_path=filename.parent/'tvt_flow_radial_velocity.png';fig.savefig(radial_path,dpi=180);plt.close(fig)

    ms,data,(nx,ny,nz)=last
    cube=data.reshape(nz,ny,nx,14)
    j=12 if ny==32 else ny//2
    slice_xz=cube[:,j,:,:]
    xx=slice_xz[:,:,0]*1e6
    zz=slice_xz[:,:,2]*1e6
    yy=slice_xz[:,:,1]*1e6
    rad=np.sqrt((xx-60)**2+(yy-60)**2+(zz-60)**2)
    vel=np.linalg.norm(slice_xz[:,:,7:10],axis=-1)*1e6
    masked=np.ma.array(vel,mask=rad>=50)
    fig,ax=plt.subplots(figsize=(7,5.6),layout='constrained')
    pcm=ax.pcolormesh(xx,zz,masked,shading='auto',norm=__import__('matplotlib').colors.LogNorm(vmin=0.2,vmax=max(vel[rad<50].max(),1)))
    ax.contour(xx,zz,rad,levels=[50],linewidths=1.4)
    fig.colorbar(pcm,ax=ax,label='Fluid |u| (µm/s), logarithmic')
    ax.set(xlabel='x (µm)',ylabel='z (µm)',title=f'{ms:g} ms: x-z section at y={yy[0,0]:.2f} µm, solid cells omitted')
    ax.set_aspect('equal')
    slice_path=filename.parent/'tvt_flow_slice_last.png';fig.savefig(slice_path,dpi=180);plt.close(fig)
    print('Processed',len(records),'zones; wrote:')
    for p in (csv_path,radial_path,slice_path):print(' ',p)
    print('Last snapshot:',records[-1])


if __name__=='__main__':
    cli=argparse.ArgumentParser()
    cli.add_argument('flow',type=Path)
    cli.add_argument('--particles',type=Path)
    args=cli.parse_args()
    audit(args.flow,args.particles)
