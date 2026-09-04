# NRG — Numerical Reactive Gas-dynamics

![Fortran](https://img.shields.io/badge/Fortran-2003%2F2008-734f96?logo=fortran)
![CMake](https://img.shields.io/badge/CMake-3.15%2B-064F8C?logo=cmake)

**NRG** is an open-source research CFD package for numerical simulation of reactive gas flows, combustion, heat and mass transfer, and coupled dispersed phases. The code is written primarily in modern Fortran and is organized as a reusable physics/data library, a set of gas-dynamic solvers, problem-generation interfaces, and post-processing utilities.

This documentation reflects the current `refactor/core-solvers-validation` branch.

## Current capabilities

The current solver executable (`computing_module`) contains the following gas-dynamic backends:

- compressible **CABARET** solver;
- **CABARET low-Mach** solver;
- **FDS-style low-Mach** solver;
- **coarse-particle method (CPM)** solver.

The shared package currently includes infrastructure for:

- multicomponent reactive mixtures and detailed finite-rate chemical kinetics;
- temperature-dependent thermodynamic and transport properties;
- molecular diffusion, including optional Soret diffusion;
- viscosity and thermal conduction;
- thermal radiation;
- Cartesian, cylindrical, and spherical computational domains;
- Eulerian and Lagrangian dispersed-particle models;
- configurable post-processing and checkpoint/restart data;
- run-level termination by physical time, wall time, or external stop/pause requests;
- OpenMP shared-memory execution;
- Tecplot field output and an optional CGNS/HDF5 field-output backend.

NRG is a research code under active development. Individual problem interfaces and advanced backends may have narrower platform support than the core CMake build; see [Installation](INSTALLATION.md#current-branch-limitations).

## Repository layout

```text
NRG/
├── CMakeLists.txt                 # top-level build configuration
├── cmake/                         # optional dependency configuration
├── computing_module/              # production solver executable
│   └── src/current_build/         # active solver implementations
├── package_library/               # shared data structures, physics and I/O
├── package_interface/             # problem/campaign generators
│   ├── src/tests/classic_tests/   # reference problem interfaces
│   └── task_setup/                # chemistry/thermophysical input database
├── package_utilities/             # standalone post-processing utilities
├── INSTALLATION.md
└── TUTORIAL.md
```

The intended workflow is deliberately split into two stages:

1. **Generate a case** with a selected `package_interface` program. The interface writes a self-contained problem directory and its `task_setup/*.inf` files.
2. **Run the case** by starting `computing_module` with that generated problem directory as the current working directory.

## Requirements

At the CMake level, NRG currently accepts:

- **CMake 3.15 or newer** for the dependency-free core build;
- Intel oneAPI Fortran **`ifx`** (`IntelLLVM` compiler ID);
- GNU **`gfortran`**;
- Intel classic **`ifort`** is retained as a compatibility path;
- OpenMP is enabled by default and may be disabled with `-DNRG_ENABLE_OPENMP=OFF`.

Additional version requirements apply to some configurations:

- **CMake 3.29+ is recommended for Visual Studio 2022 + `ifx`**, because Visual Studio generators can then select `ifx` explicitly with `-T fortran=ifx`.
- The intended bundled CGNS/HDF5 path requires **CMake 3.20+**.

The current reference `1D_laminar_velocity.f90` interface itself still contains Intel/Windows-specific code (`use ifport` and `xcopy`), so the tutorial uses **Windows + Intel `ifx`** as the reference platform for that case.

## Build quick start

### Windows 11, Visual Studio 2022, Intel oneAPI `ifx`

Open an Intel oneAPI command prompt configured for Visual Studio 2022 and run:

```cmd
git clone https://github.com/yakovenko-ivan/NRG.git
cd NRG
git switch refactor/core-solvers-validation

cmake -S . -B build ^
  -G "Visual Studio 17 2022" ^
  -A x64 ^
  -T fortran=ifx

cmake --build build --target computing_module --config Release --parallel
```

For this generator the executable is written to:

```text
build/bin/Release/computing_module.exe
```

Do **not** normally pass `-DCMAKE_Fortran_COMPILER=ifx` to a Visual Studio generator. With CMake 3.29+, `-T fortran=ifx` is the supported compiler-selection mechanism and avoids several Windows/Visual-Studio compiler-discovery problems.

### Linux core build

For Intel oneAPI:

```bash
source /opt/intel/oneapi/setvars.sh
cmake -S . -B build \
  -DCMAKE_Fortran_COMPILER=ifx \
  -DCMAKE_BUILD_TYPE=Release
cmake --build build --target computing_module --parallel
```

For GNU Fortran:

```bash
cmake -S . -B build \
  -DCMAKE_Fortran_COMPILER=gfortran \
  -DCMAKE_BUILD_TYPE=Release
cmake --build build --target computing_module --parallel
```

NRG places configuration-specific products under `build/bin/<Config>/` and `build/lib/<Config>/` for both single- and multi-configuration generators. A single-config Release build therefore uses `build/bin/Release/`, not `build/bin/`.

See [INSTALLATION.md](INSTALLATION.md) for the complete platform and CMake guide.

## Selecting a problem interface

The CMake cache variable `PACKAGE_INTERFACE_SOURCE` selects one Fortran problem generator, with the path relative to `package_interface/`.

The current default is:

```text
src/tests/classic_tests/1D_laminar_velocity.f90
```

Example:

```bash
cmake -S . -B build \
  -DCMAKE_BUILD_TYPE=Release \
  -DPACKAGE_INTERFACE_SOURCE=src/tests/classic_tests/3D_droplet_evaporation.f90

cmake --build build --target package_interface --parallel
```

The CMake target remains `package_interface`, while the produced executable is named from the selected source, for example:

```text
package_interface_1D_laminar_velocity
package_interface_3D_droplet_evaporation
```

Interfaces split across several Fortran files may additionally use the semicolon-separated cache variable `PACKAGE_INTERFACE_EXTRA_SOURCES`.

## Output formats

NRG uses a backend-neutral full-field output layer.

| Format | CMake requirement | File extension | Current notes |
|---|---|---|---|
| Tecplot TDV112 | none | `.plt` | default field-output backend used by the reference tests |
| CGNS/HDF5 | `-DNRG_ENABLE_CGNS=ON` | `.cgns` | optional; currently restricted to one MPI rank |

A problem interface selects the format through `data_save_c(..., save_format=...)`. The standard/reference path uses `save_format='tecplot'`, so no CGNS/HDF5 dependency is required for a normal first build. CGNS/HDF5 is opt-in; if a generated case requests it while NRG was built with `NRG_ENABLE_CGNS=OFF`, `computing_module` terminates with a clear error. See [INSTALLATION.md](INSTALLATION.md#9-cgnshdf5-output) for optional CGNS configuration and current backend limitations.

## Running `computing_module`

Run the solver **from inside a generated problem directory**:

```cmd
C:\path\to\NRG\build\bin\Release\computing_module.exe --num_threads=8
```

OpenMP builds default to one requested thread unless `--num_threads=<N>` is supplied by the user.

Useful command-line options currently implemented by `computing_module` include:

```text
-v, --version
-h, --help
--num_threads=<N>
--benchmark=<iterations>
```

The built-in help text currently lists only `--version` and `--help`; the thread and benchmark options are nevertheless parsed by the executable.

See [TUTORIAL.md](TUTORIAL.md) for a complete case-generation and execution example.

## Development status and portability

The root CMake configuration is being made portable, but the branch is still in a transition period:

- core compilation options support GNU and Intel compiler families;
- the reference 1D package interface is currently Intel/Windows-specific;
- OpenMP is integrated in the root build;
- MPI code paths exist in the sources, but this branch does not yet expose a root CMake option that configures MPI as a user-facing build mode;
- CGNS/HDF5 is new and optional; Tecplot remains the dependency-free output path;
- there is currently no `cmake --install`/installation target; NRG is normally run from the build tree;
- version metadata is not yet fully synchronized: `computing_module --version` reports `1.1.0`, while the top-level CMake project currently declares `VERSION 1.0`.

These limitations should be treated as current implementation status rather than permanent design constraints.

## Documentation

- [Installation and build guide](INSTALLATION.md)
- [First-simulation tutorial](TUTORIAL.md)

## Authors

NRG is developed by researchers of the **Computational Physics Laboratory, Joint Institute for High Temperatures of the Russian Academy of Sciences (JIHT RAS)**.

Leading developer: **I. S. Yakovenko, PhD**.
