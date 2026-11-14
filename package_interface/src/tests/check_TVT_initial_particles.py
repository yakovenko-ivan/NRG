#!/usr/bin/env python3
"""Check the first NRG 3D active-particle snapshot for initialization defects.

Example:
  python check_TVT_initial_particles.py tvt_active_3d/lagrangian_particles_phase_01.dat

Requires only Python's standard library. Assumes v0.1 active-particle
snapshot rows: time ID x y z ux uy uz T mass orient(3) active_force(3) unwrapped(3).
"""
from __future__ import annotations
import argparse
import math
from pathlib import Path


def read_first_snapshot(path: Path):
    with path.open(encoding='utf-8') as stream:
        header=stream.readline()
        if 'particle' not in header.lower():
            raise ValueError('Unrecognized particle output header')
        initial_time=None
        particles=[]
        for line in stream:
            if not line.strip() or line.lstrip().startswith(('#','ZONE','TITLE','VARIABLES')):
                continue
            columns=line.split()
            if len(columns)<13:
                continue
            time=float(columns[0].replace('D','E'))
            if initial_time is None:
                initial_time=time
            if abs(time-initial_time)>1e-12*max(1.0,abs(time)):
                break
            id=int(columns[1])
            xyz=tuple(float(v.replace('D','E')) for v in columns[2:5])
            velocity=tuple(float(v.replace('D','E')) for v in columns[5:8])
            orient=tuple(float(v.replace('D','E')) for v in columns[10:13])
            particles.append((id,xyz,velocity,orient))
    return initial_time, particles


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('snapshot',type=Path)
    parser.add_argument('--diameter-um',type=float,default=2.0)
    parser.add_argument('--box-um',type=float,default=120.0)
    parser.add_argument('--cloud-radius-um',type=float,default=40.0)
    parser.add_argument('--expected-count',type=int,default=600)
    args=parser.parse_args()
    t,particles=read_first_snapshot(args.snapshot)
    if not particles:
        raise SystemExit('No particles were found in the first snapshot.')
    center=0.5*args.box_um*1e-6
    diameter=args.diameter_um*1e-6
    cloud_radius=args.cloud_radius_um*1e-6
    radii=[math.dist(xyz,(center,center,center)) for _,xyz,_,_ in particles]
    orn=[math.sqrt(sum(v*v for v in orientation)) for *_,orientation in particles]
    speeds=[math.sqrt(sum(v*v for v in velocity)) for _,_,velocity,_ in particles]
    overlaps=[]
    minsep=float('inf')
    for i in range(len(particles)):
        for j in range(i+1,len(particles)):
            sep=math.dist(particles[i][1],particles[j][1])
            minsep=min(minsep,sep)
            if sep<diameter:
                overlaps.append((particles[i][0],particles[j][0],sep))
    print(f'Time of first snapshot: {t:g} s')
    print(f'Particles: {len(particles)} (expected {args.expected_count})')
    print(f'Max center distance: {max(radii)*1e6:.4f} um (initial-cloud radius {args.cloud_radius_um} um)')
    print(f'Max initial speed: {max(speeds)*1e6:.5f} um/s')
    print(f'Orientation norm range: {min(orn):.6f} ... {max(orn):.6f}')
    print(f'Closest center separation: {minsep*1e6:.4f} um')
    print(f'Overlapping center pairs (r < dp): {len(overlaps)}')
    for i,j,sep in overlaps[:10]:
        print(f'  overlap particle {i} with {j}: {sep*1e6:.4f} um')
    failed=(len(particles)!=args.expected_count or max(radii)>cloud_radius*(1+1e-6)
            or max(abs(m-1.0) for m in orn)>1e-3)
    if overlaps:
        print('WARNING: random_sphere is NOT hard-sphere excluded-volume initialization.')
    if failed:
        raise SystemExit('FAIL: count, position, or orientation sanity check did not pass.')
    print('PASS: count, location, and orientation checks (overlaps reported separately).')


if __name__=='__main__':
    main()
