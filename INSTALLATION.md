# NRG Installation and Build Guide

This guide documents the `refactor/core-solvers-validation` branch at commit `ae878d213346eb02c8986dcbcc11bfe43abfa03e` (7 October 2026). It describes out-of-source builds of the `computing_module`, selected interfaces and utilities, optional CVODE and CGNS dependencies, and profiling switches. For an end-to-end example, see [TUTORIAL.md](TUTORIAL.md).

## 1. Toolchains and minimum requirements

The root `CMakeLists.txt` declares **CMake 3.15+**, project version **1.1.0**, and supports these Fortran compiler families:

| Platform/toolchain | Core build | Notes |
|---|---|---|
| Windows + Visual Studio 2022 + Intel oneAPI `ifx` | Configured | **CMake 3.29+** recommended for `-T fortran=ifx`; install Intel VS integration |
| Linux + Intel oneAPI `ifx` | Configured | Initialize oneAPI in the shell |
| Linux + GNU `gfortran` | Configured | Does not require the Intel compiler |
| macOS + GNU `gfortran` | Configured | Individual research interfaces and physics paths require platform testing |
| Classic Intel `ifort` | Compatibility code path | Not the preferred new installation |

A default **core** build needs a Fortran compiler, CMake, and a build tool (Ninja, Make or Visual Studio). It does **not** require SUNDIALS, CGNS, or HDF5. CMake finds OpenMP by default; set `NRG_ENABLE_OPENMP=OFF` if your toolchain does not provide it.

CMake **3.20+** is required for the *bundled* CGNS/HDF5 option. Its dependency setup also requires a C compiler and downloads pinned third-party source archives on first configuration. The optional CVODE backend requires an existing SUNDIALS build with the **Fortran 2003 API** and the modules listed in Section 9.

## 2. Clone the development branch

```bash
git clone https://github.com/yakovenko-ivan/NRG.git
cd NRG
git switch refactor/core-solvers-validation
git branch --show-current
```

For a reproducible research run, record the full commit SHA:

```bash
git rev-parse HEAD
```

The build embeds the current Git HEAD, branch, `package_library` tree SHA, and library dirty-state information. Reconfigure/rebuild after changing revisions to keep this metadata in step with the executable.

## 3. Windows: Visual Studio 2022 + Intel `ifx`

Install Visual Studio 2022 with **Desktop development with C++**, the Intel oneAPI HPC Toolkit with Fortran/Visual Studio integration, Git, and CMake **3.29+**. Open a oneAPI / Visual Studio developer prompt.

```cmd
where ifx
ifx --version
cmake --version

cmake -S . -B build -G "Visual Studio 17 2022" -A x64 -T fortran=ifx
cmake --build build --target computing_module package_interface --config Release --parallel
```

Executables:

```text
build\bin\Release\computing_module.exe
build\bin\Release\package_interface_1D_laminar_velocity.exe
```

**Do not** normally set `-DCMAKE_Fortran_COMPILER=ifx` with a Visual Studio generator. CMake's `-T fortran=ifx` selects the correct Intel Fortran toolset; using a cached alternative compiler setting can produce the misleading error `ifx is not a full path and was not found in the PATH`.

On Windows, `NRG_STATIC_RUNTIME=ON` by default selects `/MT` (`/MTd` in Debug). Use `-DNRG_STATIC_RUNTIME=OFF` for `/MD`/`/MDd`. The Intel OpenMP runtime (`libiomp5md.dll`) can still be dynamically required when OpenMP is enabled.

## 4. Linux: GNU Fortran or Intel oneAPI

Install GNU tools on Ubuntu/Debian:

```bash
sudo apt update
sudo apt install git gfortran cmake make
```

Then build:

```bash
cmake -S . -B build -DCMAKE_Fortran_COMPILER=gfortran -DCMAKE_BUILD_TYPE=Release
cmake --build build --target computing_module package_interface --parallel
```

The executable paths are `build/bin/Release/computing_module` and `build/bin/Release/package_interface_1D_laminar_velocity`. Unlike older NRG versions, the 1D interface no longer imports Intel `IFPORT` or invokes Windows `xcopy`. It uses `ISO_C_BINDING` for working directories and `cmake -E` for portable copying and directory creation. **Keep CMake on `PATH` at interface execution time**.

For Intel oneAPI instead:

```bash
source /opt/intel/oneapi/setvars.sh
cmake -S . -B build-ifx -DCMAKE_Fortran_COMPILER=ifx -DCMAKE_BUILD_TYPE=Release
cmake --build build-ifx --target computing_module package_interface --parallel
```

A configured toolchain is not a guarantee that every physics/interface combination passes regression tests; validate representative cases on the target machine.

## 5. macOS: GNU Fortran

```bash
brew install git cmake gcc
xcode-select --install
cmake -S . -B build -DCMAKE_Fortran_COMPILER=gfortran -DCMAKE_BUILD_TYPE=Release
cmake --build build --target computing_module --parallel
```

If Homebrew installs a versioned `gfortran-<version>`, select that path. Intel oneAPI `ifx` is not available on macOS. Core compiler support is configured; individual interfaces and physics backends have not all been demonstrated on macOS.

## 6. Build modes, paths and targets

**Single-configuration** generators (Unix Makefiles, Ninja) use `-DCMAKE_BUILD_TYPE=Release` at configure time. If omitted, NRG defaults to `Debug`.

**Multi-configuration** generators (Visual Studio) use `--config Release` at build time. Products are always stored in configuration-specific subdirectories:

```text
build/bin/<Config>/
build/lib/<Config>/
build/mod_files/
```

Common targets:

```bash
cmake --build build --target computing_module --parallel
cmake --build build --target package_interface --parallel
cmake --build build --target package_utilities --parallel
```

On Visual Studio append `--config Release` to each command. The `package_interface` and `package_utilities` **target names do not change**, but their executable names are derived from the selected source files.

The current CMake configuration does not expose a `cmake --install` distribution or a documented root MPI enablement switch. MPI code sections exist but are not an out-of-the-box MPI build mode in this branch.

## 7. Select a problem interface or utility

`PACKAGE_INTERFACE_SOURCE` is **relative to `package_interface/`**. The default is `src/tests/classic_tests/1D_laminar_velocity.f90`.

| Interface source under `package_interface/` | Purpose |
|---|---|
| `src/tests/classic_tests/1D_laminar_velocity.f90` | 1D counter-flow / near-wall flame campaign generator |
| `src/tests/classic_tests/0D_chemical_equilibrium_validation.f90` | Mesh-free TP/HP equilibrium and H₂/O₂ crossover validation |
| `src/tests/classic_tests/0D_ignition_delay_agent.f90` | Namelist-driven constant-volume ignition case generator |
| `src/tests/classic_tests/0D_ignition_delay_campaign.f90` | 0D chemistry backend/campaign validation case generator |
| `src/tests/classic_tests/3D_droplet_evaporation.f90` | Dispersed-droplet test case |
| `src/tests/performance/transport_mixture_diffusion_benchmark.f90` | Standalone transport microbenchmark |

Select the mesh-free equilibrium executable, for example:

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release \
  -DPACKAGE_INTERFACE_SOURCE=src/tests/classic_tests/0D_chemical_equilibrium_validation.f90
cmake --build build --target package_interface --parallel
```

The generated executable is named `package_interface_0D_chemical_equilibrium_validation` (with `.exe` on Windows). The source selection applies to both Linux and Windows; on Visual Studio, keep the original `-G`, `-A` and `-T` options when reconfiguring, and build with `--config Release`.

Interfaces with companion sources can set `PACKAGE_INTERFACE_EXTRA_SOURCES` to a semicolon-separated list of paths relative to `package_interface/`. For example:

```bash
-DPACKAGE_INTERFACE_EXTRA_SOURCES="src/helper_a.f90;src/helper_b.f90"
```

This is an optional mechanism, not required by the included standalone examples.

`PACKAGE_UTILITY_SOURCE` selects a source relative to `package_utilities/`; the default is `src/leading_point_velocity.f90`. Other available programs are `src/output_merger.f90` and `src/tecplot_merger.f90`. Example:

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release \
  -DPACKAGE_UTILITY_SOURCE=src/output_merger.f90
cmake --build build --target package_utilities --parallel
```

Changing an interface or utility source generally only requires CMake reconfiguration, not deleting the build directory. When switching compilers, generators or runtime policies, start a **new** build directory.

## 8. OpenMP and runtime options

Options declared by the root build:

| Cache option | Default | Effect |
|---|---|---|
| `NRG_ENABLE_OPENMP` | `ON` | OpenMP support; root CMake requires `OpenMP::OpenMP_Fortran` |
| `NRG_STATIC_RUNTIME` | `ON` | Static Intel/MSVC runtime on Windows |
| `NRG_ENABLE_CVODE` | `OFF` | SUNDIALS CVODE chemistry integration |
| `NRG_ENABLE_CHEMISTRY_PROFILE` | `OFF` | Chemistry performance counters/timing (`CHEMISTRY_PROFILE` preprocessor symbol) |
| `NRG_ENABLE_FDS_PRESSURE_PROFILE` | `OFF` | FDS pressure/predictor profiling; writes `fds_pressure_profile.csv` |
| `NRG_ENABLE_CGNS` | `OFF` | CGNS/HDF5 field output |
| `NRG_USE_SYSTEM_CGNS` | `OFF` | Use an installed CGNS rather than FetchContent dependencies (only meaningful when CGNS is enabled) |

Disable OpenMP, for instance:

```bash
cmake -S . -B build-serial -DCMAKE_Fortran_COMPILER=gfortran \
  -DCMAKE_BUILD_TYPE=Release -DNRG_ENABLE_OPENMP=OFF
cmake --build build-serial --target computing_module --parallel
```

Build parallelism (`cmake --build ... --parallel`) and runtime OpenMP parallelism are different. At runtime, the solver uses **one requested thread** unless invoked with `--num_threads=N`. `--benchmark=N` stops after an iteration count; add `--benchmark_checkpoint` if a restart state is needed at the end of the timed run. See [TUTORIAL.md](TUTORIAL.md).

## 9. Optional SUNDIALS CVODE chemistry backend

The default `slatec`, `qss1`, and `qss2` CPU chemistry backends are **built in** and do not require SUNDIALS. Only the `cvode` option requires an external library.

NRG uses `find_package(SUNDIALS CONFIG REQUIRED)` and verifies all of these imported targets:

```text
SUNDIALS::fcvode_mod
SUNDIALS::fnvecserial_mod
SUNDIALS::fsunmatrixdense_mod
SUNDIALS::fsunlinsoldense_mod
```

It additionally links the corresponding native CVODE/serial/dense components (`cvode`, `nvecserial`). Install a compatible SUNDIALS built with the **Fortran 2003 interface** and those serial/dense modules; an installation with only C/C++ CVODE is insufficient.

Example, with a separately installed SUNDIALS prefix:

```bash
cmake -S . -B build-cvode -DCMAKE_Fortran_COMPILER=gfortran \
  -DCMAKE_BUILD_TYPE=Release -DNRG_ENABLE_CVODE=ON \
  -DCMAKE_PREFIX_PATH=/path/to/sundials/install
cmake --build build-cvode --target computing_module --parallel
```

A generated case selects its integrator in `solver_options` (`chemistry_backend='cvode'`), or via the `0D_ignition_delay_campaign.f90` namelist. **Building CVODE support does not automatically switch cases to CVODE.** Conversely, requesting `cvode` from a binary built with `NRG_ENABLE_CVODE=OFF` is unsupported.

The relevant campaign namelist group is `/chemistry_backend_config/`, with `chemistry_backend` values `slatec`, `qss1`, `qss2`, or `cvode`. Configurable CVODE tolerances/max steps and the QSS1/QSS2 adaptivity controls are defined in that interface and in `solver_options_class.f90`.

## 10. Optional CGNS/HDF5 output

Tecplot output needs no external library. To write `.cgns`, configure `NRG_ENABLE_CGNS=ON` **and** generate cases requesting `save_format='cgns'` in `data_save_c(...)`. These settings are independent; changing CMake alone does not change an existing case's save format.

### Bundled dependencies (CMake 3.20+)

```bash
cmake -S . -B build-cgns -DCMAKE_BUILD_TYPE=Release -DNRG_ENABLE_CGNS=ON
cmake --build build-cgns --target computing_module --parallel
```

The branch's `cmake/NRGCGNSDependencies.cmake` pins **HDF5 1.14.6** and **CGNS 4.5.2** via `FetchContent` and uses `cmake/hdf5-buildtree/hdf5-config.cmake.in`. **That compatibility template is present in the inspected branch**; the missing-template warning in earlier installation documentation is obsolete. Initial configuration needs network access to retrieve the archives. The optional dependency stack has its own compiler/platform constraints; availability of the files is not proof of successful builds on every target machine.

### Installed CGNS

```bash
cmake -S . -B build-cgns -DCMAKE_BUILD_TYPE=Release \
  -DNRG_ENABLE_CGNS=ON -DNRG_USE_SYSTEM_CGNS=ON \
  -DCMAKE_PREFIX_PATH=/path/to/cgns/install
```

As needed, set `cgns_DIR`, `NRG_CGNS_INCLUDE_DIR`, and/or `NRG_CGNS_LIBRARY`. The system fallback may also need an HDF5 C installation. The current writer is **single-MPI-rank only**. For a dependency-free first run, keep the reference interface's `save_format='tecplot'`.

If a case selects CGNS with the backend disabled, `computing_module` reports `Data save: CGNS requested but NRG was built with NRG_ENABLE_CGNS=OFF`. Reconfigure/rebuild with CGNS enabled or regenerate the case with Tecplot output.

## 11. Performance instrumentation

Instrumentation is **compile-time opt-in** and therefore normally belongs in a separate build tree to avoid perturbing normal benchmarks:

```bash
cmake -S . -B build-profile -DCMAKE_Fortran_COMPILER=gfortran \
  -DCMAKE_BUILD_TYPE=Release \
  -DNRG_ENABLE_CHEMISTRY_PROFILE=ON \
  -DNRG_ENABLE_FDS_PRESSURE_PROFILE=ON
cmake --build build-profile --target computing_module --parallel
```

The FDS backend writes `fds_pressure_profile.csv` in the **case working directory** when pressure profiling is compiled in and the relevant solver path executes. Chemistry profiling enables additional counters/timing output within `chemical_kinetics_solver`; it is not the same as `--benchmark=N`.

A separate standalone transport microbenchmark is selectable as `PACKAGE_INTERFACE_SOURCE=src/tests/performance/transport_mixture_diffusion_benchmark.f90`. From a directory containing `task_setup/`, invoke:

```bash
./package_interface_transport_mixture_diffusion_benchmark --calls=100000 --repeats=5 --warmup=10000
```

The output file is `transport_microbenchmark.csv`. The benchmark is not a CFD case and does not use `computing_module`.

## 12. Metadata and scientific reproducibility

CMake creates `nrg_build_info` from `package_library/src/nrg_build_info.f90.in`. Generated cases and solver logs report:

```text
Git branch
Git HEAD
Package library tree
Package library modified
```

The tree SHA tracks the committed shared library source. A dirty working tree means it is **not** a complete source fingerprint; retain the local changes as well. After branch switches, fresh commits, or modified build options, re-run configuration and rebuild before generating/running provenance-sensitive cases.

Post-processor `.dat` files now start with Tecplot-style `TITLE=` and `VARIABLES=` metadata, including time units, operation names, and coordinates. Scripts that previously read them as headerless numeric tables must skip/parse these header records. Full-field `.plt` and `.cgns` output is independent of these time-series files.

## 13. Build inspection and troubleshooting

Inspect cached options:

```bash
cmake -S . -B build -LA
```

Useful compiler flags are set by NRG: GNU uses `-cpp` and `-ffree-line-length-512`; Intel uses preprocessing/traceback options. Debug enables runtime checking and traps; Release enables optimization. Prefer **Release** for performance campaigns and **Debug/RelWithDebInfo** for diagnosing correctness issues.

**`ifx is not a full path` / wrong Intel compiler.** On Visual Studio 2022 use CMake 3.29+ and configure a clean directory with `-T fortran=ifx`. Verify `where ifx`, oneAPI setup, and Intel Visual Studio integration. Do not reuse an old cache with `CMAKE_Fortran_COMPILER=ifx`.

**`libiomp5md.dll` not found.** Run under an initialized oneAPI environment or provide Intel OpenMP runtime libraries. `NRG_STATIC_RUNTIME=ON` does not make Intel OpenMP entirely static.

**OpenMP not found.** Install/configure compiler OpenMP support or set `-DNRG_ENABLE_OPENMP=OFF`.

**Cannot locate SUNDIALS Fortran imported targets.** Set `CMAKE_PREFIX_PATH` or `SUNDIALS_DIR` to a Fortran-enabled SUNDIALS installation, not just a C-library prefix. Alternatively disable `NRG_ENABLE_CVODE` and select `slatec`, `qss1`, or `qss2`.

**Missing chemistry/thermophysical files at runtime.** Run an interface from a directory containing a copy of `package_interface/task_setup/`. NRG resolves data under relative paths such as `task_setup/chemical_mechanisms/KEROMNES.txt`. For a generated CFD case, change into the leaf case directory before launching the solver.

**`cmake` command not found when executing the 1D interface.** The current portable generator invokes `cmake -E make_directory` and `cmake -E copy_directory` at runtime. Install CMake and keep it on `PATH`.

**CGNS requested but disabled.** The case requests CGNS, while the solver was compiled with `NRG_ENABLE_CGNS=OFF`. Rebuild or select Tecplot. The old missing-`hdf5-config.cmake.in` diagnosis does **not** apply to this commit.

**Cannot find executable at `build/bin/`.** NRG writes to `build/bin/Release/`, `build/bin/Debug/`, etc., even for GNU single-configuration builds.

**Unexpectedly slow executable.** For Makefiles/Ninja configure with `-DCMAKE_BUILD_TYPE=Release`; on Visual Studio specify `--config Release` when building.

**Stale compiler/toolset cache.** Use a new build directory, or remove the old one (`rm -rf build` on Linux/macOS; `rmdir /S /Q build` in Windows CMD), then configure again. Reconfiguration is normally enough for just changing `PACKAGE_INTERFACE_SOURCE`.

## 14. Continue with the tutorial

[TUTORIAL.md](TUTORIAL.md) covers a quick **mesh-free equilibrium validation** and the full **1D CFD case generation → integration** sequence. It also identifies the three cases generated by the current 1D loop bounds and explains the new post-processor headers.
