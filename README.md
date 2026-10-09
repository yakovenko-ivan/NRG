# NRG — Numerical Reactive Gas-dynamics

![Fortran](https://img.shields.io/badge/Fortran-2003%2F2008-734f96?logo=fortran)
![CMake](https://img.shields.io/badge/CMake-3.15%2B-064F8C?logo=cmake)

**NRG** is an open-source research CFD package for reactive gas flows, combustion, heat and mass transfer, and coupled dispersed phases. Its numerical solvers and physical submodels are implemented primarily in modern Fortran. Case generation, time integration, output, and scientific diagnostics are separate components.

**Documentation target:** [`refactor/core-solvers-validation`](https://github.com/yakovenko-ivan/NRG/tree/refactor/core-solvers-validation), inspected at commit `ae878d213346eb02c8986dcbcc11bfe43abfa03e` (7 October 2026). The branch is under active development; verify the interface source and CMake options when using a later revision.

## Capabilities

### Flow solvers

The `computing_module` executable contains four selectable gas-dynamic backends:

| Solver selector | Method |
|---|---|
| `CABARET` | Compressible CABARET |
| `CABARET_low_mach` | Low-Mach CABARET |
| `fds_low_mach` | FDS-style low-Mach formulation |
| `cpm` | Coarse-particle method (CPM) |

The shared library supplies multicomponent thermodynamics (including JANAF/NASA polynomial properties), finite-rate kinetics, molecular diffusion and optional Soret diffusion, viscosity, Fourier heat transfer, thermal radiation, and Eulerian/Lagrangian dispersed-phase models. Computational domains can use Cartesian, cylindrical, or spherical coordinates. Available functionality depends on the chosen solver and problem interface; not every combination of physics, solver, and geometry is validated.

### Chemistry and validation

Chemistry integration is selected **per generated case**, independently of the flow solver:

| `chemistry_backend` | Availability | Notes |
|---|---|---|
| `slatec` | Included; default | Existing stiff ODE backend |
| `qss1` | Included | First-order quasi-steady-state integration with conservation correction |
| `qss2` | Included | Adaptive QSS predictor–corrector backend with user-configurable error tolerance |
| `cvode` | Optional build | SUNDIALS CVODE, requiring `NRG_ENABLE_CVODE=ON` and Fortran-enabled SUNDIALS |

The **mesh-free** `0D_chemical_equilibrium_validation.f90` interface uses the same thermophysical and reaction-rate infrastructure as the CFD solvers. It supports fixed-`T,p` Gibbs equilibrium (`tp`), adiabatic constant-pressure equilibrium (`hp`), and the H₂/O₂ branching-versus-HO₂-stabilization crossover temperature (`crossover`). It writes CSV states and residual diagnostics; this is a standalone validation executable, **not** a `computing_module` CFD case generator.

Additional included interfaces provide 0D ignition-delay campaign generation and a mixture-diffusion transport microbenchmark. The current 1D flame interface includes active flame-stabilization controls for its counter-flow configurations.

### Execution, provenance and output

- OpenMP CPU execution; optional iteration benchmarking and restart-at-benchmark-end.
- Case-level run control by simulated time, elapsed wall time, or external `run_control.stop` / `run_control.pause` requests.
- Git source revision and `package_library` tree/dirty-state provenance in generated case and execution logs.
- Full-field **Tecplot TDV112** (`.plt`) output by default; optional **CGNS/HDF5** (`.cgns`) output, currently restricted to one MPI rank.
- Human-readable Tecplot-style `TITLE=` / `VARIABLES=` headers in time-series post-processor `.dat` files, with operation names and coordinate units.
- Optional CMake-enabled chemistry profiling and FDS pressure-solver profiling (`fds_pressure_profile.csv`).

## Repository layout

```text
NRG/
├── CMakeLists.txt                    Root build options and targets
├── cmake/                            Optional CGNS/HDF5 dependency setup
├── computing_module/src/current_build/
│   ├── main.f90                      Solver entry point
│   ├── *solver.f90                    Flow and physics solvers
│   └── problem_controls/             Flame, ignition, forcing controls
├── package_library/src/             Shared fields, chemistry, mesh, I/O
├── package_interface/
│   ├── src/tests/classic_tests/      0D/1D/3D interfaces
│   ├── src/tests/performance/        Transport microbenchmark
│   └── task_setup/                   Chemistry, thermo, example namelist
├── package_utilities/src/           Output and trajectory utilities
├── INSTALLATION.md
└── TUTORIAL.md
```

## Build quick start

**Compiler support declared by CMake:** GNU `gfortran`, Intel oneAPI `ifx`, and legacy Intel `ifort`. OpenMP is enabled by default. The dependency-free core requires CMake **3.15+**. **CMake 3.29+** is recommended for Visual Studio 2022 with Intel `ifx`.

### Linux: GNU Fortran

```bash
git clone https://github.com/yakovenko-ivan/NRG.git
cd NRG
git switch refactor/core-solvers-validation

cmake -S . -B build -DCMAKE_Fortran_COMPILER=gfortran -DCMAKE_BUILD_TYPE=Release
cmake --build build --target computing_module package_interface --parallel
```

The executables reside in `build/bin/Release/`, including on single-configuration generators. A compiler is sufficient for the core build; an optional CVODE or CGNS build has additional requirements.

### Windows: Visual Studio 2022 + Intel `ifx`

Open an Intel oneAPI / Visual Studio developer prompt:

```cmd
git clone https://github.com/yakovenko-ivan/NRG.git
cd NRG
git switch refactor/core-solvers-validation

cmake -S . -B build -G "Visual Studio 17 2022" -A x64 -T fortran=ifx
cmake --build build --target computing_module package_interface --config Release --parallel
```

The solver is `build\bin\Release\computing_module.exe`. For Visual Studio, use `-T fortran=ifx`, **not** `-DCMAKE_Fortran_COMPILER=ifx`; change toolchains in a fresh build directory.

See [INSTALLATION.md](INSTALLATION.md) for GNU/Intel instructions, feature switches, optional libraries, and troubleshooting.

## How to run NRG

Most CFD problems follow **two distinct stages**:

1. Select and build a `package_interface` program. Run it from a working directory containing a copy of the repository's `package_interface/task_setup/` folder; it generates one or more self-contained case directories.
2. Change to a **generated leaf case directory** (containing `problem_setup.log` and `task_setup/`) and run `computing_module` there. Results are written relative to that working directory.

Use `PACKAGE_INTERFACE_SOURCE` (relative to `package_interface/`) to select an interface, for example:

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release \
  -DPACKAGE_INTERFACE_SOURCE=src/tests/classic_tests/1D_laminar_velocity.f90
cmake --build build --target package_interface --parallel
```

The CMake **target** remains `package_interface`, while the executable is named `package_interface_<source-basename>`.

**Current 1D default:** the source presently generates **three** Cartesian, counter-flow, FDS/KEROMNES cases with **8 vol.% H₂ in air** and `dx = 2.0e-4`, `1.0e-4`, and `5.0e-5 m`. These values are hard-coded in the campaign loop and should not be confused with earlier near-wall/17% defaults. The generator now uses portable `ISO_C_BINDING` working-directory routines and `cmake -E` filesystem commands; CMake must remain on `PATH` when running it. See [TUTORIAL.md](TUTORIAL.md).

`computing_module` accepts:

```text
-v, --version
-h, --help
--num_threads=N
--benchmark=N
--benchmark_checkpoint        (requires --benchmark=N)
```

Use `--num_threads=N` for an OpenMP run; the requested thread count defaults to one. `--benchmark_checkpoint` asks for a restart checkpoint at the benchmark endpoint.

## Other interfaces and workflows

- `src/tests/classic_tests/0D_chemical_equilibrium_validation.f90`: mesh-free equilibrium and crossover validation; works with the example [`equilibrium_validation.nml`](package_interface/task_setup/equilibrium_validation.nml).
- `src/tests/classic_tests/0D_ignition_delay_agent.f90`: self-contained constant-volume reactor case generated from `setup_input.nml`.
- `src/tests/classic_tests/0D_ignition_delay_campaign.f90`: campaign-oriented constant-volume reactor case with `slatec` / `qss1` / `qss2` / `cvode` controls, run limits and case identification.
- `src/tests/classic_tests/3D_droplet_evaporation.f90`: dispersed-phase example.
- `src/tests/performance/transport_mixture_diffusion_benchmark.f90`: standalone transport timing/consistency output to `transport_microbenchmark.csv`.

Each interface has its **own input contract**. The equilibrium and transport benchmarks execute directly and do not require a generated CFD case or `computing_module`.

## Development status

NRG is a research code, not a prepackaged CFD distribution. The root build does not currently define `cmake --install` or a user-facing MPI build switch. The CGNS/HDF5 dependency path is opt-in (the previously missing `hdf5-config.cmake.in` compatibility template is now committed). CPU support for GNU and Intel compilers is configured in CMake, but individual research interfaces and physics combinations should be built and tested on each target platform. Version reporting in the CMake project and `computing_module --version` is **1.1.0**.

## Documentation

- [Installation and build configuration](INSTALLATION.md)
- [Tutorial: equilibrium validation and 1D CFD case](TUTORIAL.md)

## Authors

NRG is developed by researchers of the **Computational Physics Laboratory, Joint Institute for High Temperatures of the Russian Academy of Sciences (JIHT RAS)**.

Leading developer: **I. S. Yakovenko, PhD**.
