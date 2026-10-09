# NRG Tutorial — Equilibrium Validation and a 1D CFD Case

This tutorial matches the [`refactor/core-solvers-validation`](https://github.com/yakovenko-ivan/NRG/tree/refactor/core-solvers-validation) branch at `ae878d213346eb02c8986dcbcc11bfe43abfa03e` (7 October 2026). It introduces **two different execution models**:

1. **Mesh-free validation** — build `0D_chemical_equilibrium_validation.f90` and run the resulting executable directly. No computational mesh, generated CFD case, or `computing_module` is required.
2. **Flow simulation** — build a problem-interface generator, generate case directories containing `task_setup/*.inf`, then start `computing_module` from one generated case directory.

The examples use the built-in KEROMNES input data. The 1D generator supports Windows and Unix-like file operations via `ISO_C_BINDING` and `cmake -E`; it is **no longer tied to `IFPORT` or Windows `xcopy`**. This describes supported source paths, not a claim that every compiler and physical case has been regression-tested.

## 1. Build prerequisites

Install Git, CMake, a supported Fortran compiler, and an associated native build tool:

- **Linux:** `gfortran` + Make/Ninja and CMake 3.15+.
- **Windows:** Visual Studio 2022 with Intel oneAPI `ifx` and Fortran integration; CMake 3.29+ recommended.

Fetch the source:

```bash
git clone https://github.com/yakovenko-ivan/NRG.git
cd NRG
git switch refactor/core-solvers-validation
```

See [INSTALLATION.md](INSTALLATION.md) for detailed prerequisites, the optional CVODE backend, and CGNS configuration. Neither optional dependency is required here.

## 2. Tutorial A — mesh-free chemical equilibrium

This is the shortest way to validate that NRG reads its thermodynamic and reaction data without launching a CFD calculation.

### 2.1 Build the validation interface

**Linux / GNU:**

```bash
cmake -S . -B build -DCMAKE_Fortran_COMPILER=gfortran \
  -DCMAKE_BUILD_TYPE=Release \
  -DPACKAGE_INTERFACE_SOURCE=src/tests/classic_tests/0D_chemical_equilibrium_validation.f90
cmake --build build --target package_interface --parallel
```

**Windows CMD / Intel `ifx`:**

```cmd
cmake -S . -B build -G "Visual Studio 17 2022" -A x64 -T fortran=ifx ^
  -DPACKAGE_INTERFACE_SOURCE=src/tests/classic_tests/0D_chemical_equilibrium_validation.f90
cmake --build build --target package_interface --config Release --parallel
```

The output executable is:

```text
build/bin/Release/package_interface_0D_chemical_equilibrium_validation
```

Add `.exe` on Windows.

### 2.2 Prepare the working directory

All NRG chemistry/thermophysical data paths are resolved relative to the process **working directory**. From the repository root, copy the template input tree:

```bash
cmake -E copy_directory package_interface/task_setup build/bin/Release/task_setup
```

This command also works in Windows CMD. You should then have:

```text
build/bin/Release/
├── package_interface_0D_chemical_equilibrium_validation[.exe]
└── task_setup/
    ├── chemical_mechanisms/KEROMNES.txt
    ├── thermophysical_data/KEROMNES_THERMO.txt
    ├── thermophysical_data/KEROMNES_TRANSDATA.txt
    └── equilibrium_validation.nml
```

### 2.3 Run TP, HP and crossover calculations

On Linux:

```bash
cd build/bin/Release
./package_interface_0D_chemical_equilibrium_validation task_setup/equilibrium_validation.nml
```

On Windows, change into `build\bin\Release` and run:

```cmd
package_interface_0D_chemical_equilibrium_validation.exe task_setup\equilibrium_validation.nml
```

The bundled namelist starts with the stoichiometric H₂/O₂/N₂ mixture `2 : 1 : 3.762` (moles), pressure `101325 Pa`, fixed-temperature TP test at `1500 K`, HP test with initial temperature `300 K`, and crossover search over `400–3000 K`. `mode='all'` runs all three calculations. You may choose `tp`, `hp` or `crossover` instead.

This program writes CSV files directly to its working directory, typically:

```text
equilibrium_validation_input.csv
equilibrium_validation_initial_composition.csv
equilibrium_validation_tp_state.csv
equilibrium_validation_tp_diagnostics.csv
equilibrium_validation_hp_state.csv
equilibrium_validation_hp_diagnostics.csv
equilibrium_validation_crossover_diagnostics.csv
equilibrium_validation_crossover_reactions.csv
```

Inspect `converged`, elemental/mass/enthalpy residuals, and the crossover balance residual; do **not** treat a printed temperature or composition as validated solely because the program completed. For quantitative validation, compare the equilibrium states against a reference implementation such as Cantera using matching mechanism, species data, pressure, and initial composition.

**Key distinction:** no `problem_setup.log`, `solver_setup.inf`, CFD grid, or `computing_module` is created by this executable.

## 3. Tutorial B — a generated 1D reactive-flow case

The 1D interface is a source-driven campaign generator. Build its executable, run it with `task_setup/` available, and then select a generated case to integrate.

### 3.1 Select and build the 1D interface

From the **repository root**, return from the validation directory first if necessary.

Linux:

```bash
cmake -S . -B build -DCMAKE_Fortran_COMPILER=gfortran \
  -DCMAKE_BUILD_TYPE=Release \
  -DPACKAGE_INTERFACE_SOURCE=src/tests/classic_tests/1D_laminar_velocity.f90
cmake --build build --target package_interface computing_module --parallel
```

Windows CMD:

```cmd
cmake -S . -B build -G "Visual Studio 17 2022" -A x64 -T fortran=ifx ^
  -DPACKAGE_INTERFACE_SOURCE=src/tests/classic_tests/1D_laminar_velocity.f90
cmake --build build --target package_interface computing_module --config Release --parallel
```

Executables in `build/bin/Release/`:

```text
package_interface_1D_laminar_velocity[.exe]
computing_module[.exe]
```

`PACKAGE_INTERFACE_SOURCE` controls only which generator is built. The solver still contains all four gas-dynamic backends. Changing *problem parameters* within the 1D interface requires rebuilding only `package_interface`.

### 3.2 Copy the master `task_setup/` files

If you did not do so for Tutorial A, execute from the repository root:

```bash
cmake -E copy_directory package_interface/task_setup build/bin/Release/task_setup
```

The 1D generator expects `task_setup` in **its current working directory**, not at an arbitrary path relative to its executable. Keep `cmake` executable on `PATH`, because the source uses `cmake -E make_directory` and `cmake -E copy_directory` while generating cases.

### 3.3 Generate the current 1D campaign

Linux:

```bash
cd build/bin/Release
./package_interface_1D_laminar_velocity
```

Windows CMD (from `build\bin\Release`):

```cmd
package_interface_1D_laminar_velocity.exe
```

**The current hard-coded loop bounds generate THREE cases, not one:**

| Parameter | Current values |
|---|---|
| Physical setup | `counter_flow` (`cf`) |
| Coordinate system | Cartesian |
| Gas-dynamic solver | `fds_low_mach` (`FDS`) |
| Mechanism | `KEROMNES` |
| Fresh-mixture hydrogen mole fraction | **8 vol.%** |
| Grid spacings | **`2.0e-4`, `1.0e-4`, `5.0e-5 m`** |

The case tree is:

```text
1D_LBV_test/
└── cf/
    └── cartesian/
        └── FDS/
            └── KEROMNES/
                └── <8%-H2-and-phi-label>/
                    ├── dx_2.0e-04/
                    ├── dx_1.0e-04/
                    └── dx_5.0e-05/
```

The concentration/equivalence-ratio directory is generated by `str_r(...)` in the source, so inspect it rather than guessing its exact spelling. This is **not** the earlier `near_wall` (`nw`), 17% H₂, `dx=1.25e-5 m` setup described in old NRG tutorials.

To enumerate the generated leaves on Linux:

```bash
find 1D_LBV_test -name problem_setup.log -print
```

Windows CMD:

```cmd
dir /S /B 1D_LBV_test\problem_setup.log
```

Case generation creates substantial data and three different mesh sizes. Select the coarsest case first for a basic functional check; numerical fidelity still requires a convergence study.

### 3.4 Understand the generated CFD case

Each generated **leaf case directory** contains something like:

```text
<case>/
├── problem_setup.log
├── task_setup/
│   ├── domain_data.inf
│   ├── field_data.inf
│   ├── chemical_data.inf
│   ├── thermophysical_setup.inf
│   ├── solver_setup.inf
│   ├── problem_controls.inf
│   ├── boundary_conditions_setup.inf
│   ├── data_save_setup.inf
│   ├── data_io_setup.inf
│   ├── post_processors_manager_setup.inf
│   ├── post_processor*.inf
│   ├── chemical_mechanisms/...
│   └── thermophysical_data/...
├── data_save/
└── data_output/
```

The source writes initial fields and a human-readable `problem_setup.log`. Its active flame-stabilization control is separate from the gas-dynamic solver. **Counter-flow cases** select `laminar_burning_velocity` stabilization; **near-wall cases** disable active inlet stabilization (`mode='none'`).

Optional `task_setup/run_control.inf` specifies explicit simulated-time/wall-time termination. A legacy-generated case without this file remains valid and retains legacy termination behavior.

### 3.5 Execute one generated case

Change to the directory containing that case's `problem_setup.log` and `task_setup/`. On Linux, start the solver with an **absolute path** to the executable:

```bash
/path/to/NRG/build/bin/Release/computing_module --num_threads=4
```

On Windows CMD:

```cmd
C:\path\to\NRG\build\bin\Release\computing_module.exe --num_threads=4
```

These are path templates: replace `/path/to/NRG` or `C:\path\to\NRG` with your actual repository location. You must launch the solver from **inside the generated leaf case directory**. Merely specifying a file path to the solver while remaining in `build/bin/Release` is not equivalent.

The solver uses one requested OpenMP thread unless `--num_threads=N` is specified. For an iteration-limited smoke test, you can use:

```bash
/path/to/NRG/build/bin/Release/computing_module --num_threads=4 --benchmark=20
```

For a restart checkpoint at the benchmark endpoint, add `--benchmark_checkpoint`. Do not interpret the first 20 iterations as a physically converged flame simulation.

### 3.6 Check output and provenance

The case's `problem_setup.log` records configuration and source revision data (Git branch/HEAD, `package_library` tree/dirty status). Source provenance is useful only when the interface and solver were reconfigured/rebuilt for the revision they actually use.

Normal output locations:

```text
<case>/data_save/            Full-field snapshots (*.plt by default)
<case>/data_output/          Restart/checkpoint data
<case>/proc1.dat             Flame-front time-series diagnostics
<case>/problem_setup.log     Setup and solver information
```

The reference interface requests Tecplot `.plt` field output at simulated millisecond intervals. **Post-processor `.dat` files now have headers**, for example a `TITLE="..."` line followed by a `VARIABLES=` line describing time, monitored quantities, and coordinates. The output columns are not anonymous. If a Python script previously used `numpy.loadtxt('proc1.dat')` or `pandas.read_csv` as a plain headerless table, update it to account for the first two records and whitespace-separated numeric data.

Variable names depend on the generated post-processor operations. For `proc1`, they include a temperature-gradient operation and pressure/density/temperature transducers with signed offsets.

The default `.plt` snapshots and `proc1.dat` are **different output streams**. The optional `.cgns` full-field backend is selected through both the build flag and the interface's `save_format`.

## 4. Adjust the physical/numerical sweep

Edit `package_interface/src/tests/classic_tests/1D_laminar_velocity.f90` and rebuild the interface. The six source-loop selectors are:

| Selector | Meaning | Current default |
|---|---|---|
| `task1` | Physical setup: `counter_flow`, precomputed flamelet, or `near_wall` | `1` |
| `task2` | Coordinate system | `1` (Cartesian) |
| `task3` | Gas-dynamic method | `1` (FDS low-Mach) |
| `task4` | Chemical mechanism | `1` (KEROMNES) |
| `task5` | H₂ mole fraction, percent | `8` |
| `task6` | Spatial resolution index | `1:3` |

To switch to a **single** near-wall case, for example, select `task1=3` and a single `task6` index. The supported `task6` spacing map in the source is `0 → 4.0e-4`, `1 → 2.0e-4`, `2 → 1.0e-4`, `3 → 5.0e-5`, `4 → 2.5e-5`, `5 → 1.25e-5`, `6 → 6.25e-6 m`.

```bash
cmake --build build --target package_interface --parallel
```

For Windows Visual Studio add `--config Release`. Re-run the generator from a directory containing `task_setup/`. This modifies/creates case directories, but does **not** require recompiling `computing_module` unless you changed the underlying solver/library code or its compile-time options.

## 5. Chemistry integrators and 0D campaigns

`solver_options` separates the chemical ODE backend from the CFD solver. The case-configurable values are:

| Backend | Extra CMake dependency | Important control names |
|---|---|---|
| `slatec` | None | `chemistry_slatec_accuracy`, `chemistry_slatec_max_steps`, etc. |
| `qss1` | None | `chemistry_qss1_relative_change_limit`, `chemistry_qss1_minimum_internal_step`, etc. |
| `qss2` | None | `chemistry_qss2_error_tolerance`, `chemistry_qss2_minimum_internal_step`, `chemistry_qss2_max_steps`, etc. |
| `cvode` | SUNDIALS with Fortran 2003 interfaces | `chemistry_cvode_relative_tolerance`, `chemistry_cvode_absolute_tolerance`, `chemistry_cvode_max_steps` |

The `0D_ignition_delay_campaign.f90` source reads a **user-supplied** `setup_input.nml` (or another namelist path given as its first argument). It requires distinct groups `/case_config/`, `/reactor_config/`, `/mixture_config/`, `/physics_config/`, `/chemistry_backend_config/`, `/run_control_config/`, and `/output_config/`. It creates an individual constant-volume reactor case for `computing_module`, preserving the configuration for repeatable backend comparisons. It is **not** the same executable as the mesh-free equilibrium validation tool.

For instance, its `/chemistry_backend_config/` can contain:

```fortran
&chemistry_backend_config
  chemistry_backend = 'qss2'
  chemistry_qss2_error_tolerance = 1.0e-3
  chemistry_qss2_minimum_internal_step = 1.0e-10
  chemistry_qss2_active_concentration_fraction = 1.0e-7
  chemistry_qss2_max_steps = 100000
/
```

This is **one group**, not a complete campaign input file. Supply and validate all seven groups before running the generator. For precise input names and constraints, consult the current source. QSS1/QSS2 are CPU backends included in the ordinary build; CVODE requires `-DNRG_ENABLE_CVODE=ON` at compile time.

## 6. Run control, benchmarking and profiling

The solver understands:

```text
--version
--help
--num_threads=N
--benchmark=N
--benchmark_checkpoint    requires --benchmark=N
```

A generated case can optionally contain `task_setup/run_control.inf` with termination modes `none`, `simulation_time`, `wall_time`, or `either`. The solver also recognizes `run_control.stop` and `run_control.pause` request files in the case root; time-budget/pause termination can require a checkpoint. No `run_control.inf` means no new termination policy is imposed.

To profile instead of just timing an iteration-limited run, configure a **separate** instrumented build:

```bash
cmake -S . -B build-profile -DCMAKE_Fortran_COMPILER=gfortran \
  -DCMAKE_BUILD_TYPE=Release \
  -DNRG_ENABLE_CHEMISTRY_PROFILE=ON \
  -DNRG_ENABLE_FDS_PRESSURE_PROFILE=ON
cmake --build build-profile --target computing_module --parallel
```

If a profiled FDS case executes the corresponding pressure-solver path, `fds_pressure_profile.csv` is written in the current **case directory**. Optional chemistry counters/timing are enabled by `NRG_ENABLE_CHEMISTRY_PROFILE`; these are instrumentation options, not chemistry backend selectors.

For an isolated mixture-diffusion cost experiment, select `src/tests/performance/transport_mixture_diffusion_benchmark.f90` as `PACKAGE_INTERFACE_SOURCE`, build `package_interface`, run its executable from a directory containing `task_setup/`, and pass (for example) `--calls=100000 --repeats=5 --warmup=10000`. The resulting CSV is `transport_microbenchmark.csv`.

## 7. Optional CGNS/HDF5 output

The reference 1D interface uses `save_format='tecplot'`. To request CGNS instead, edit the interface's `data_save_c(...)` call to use `save_format='cgns'`, regenerate the case, and use a solver built with:

```bash
-DNRG_ENABLE_CGNS=ON
```

A system-CGNS configuration may also use `-DNRG_USE_SYSTEM_CGNS=ON`. The bundled dependency route pins HDF5 1.14.6 and CGNS 4.5.2 and requires CMake 3.20+, a C compiler, and network access on first configuration. The earlier warning that `cmake/hdf5-buildtree/hdf5-config.cmake.in` was missing is **outdated**; that file is present in the inspected branch. The CGNS writer currently supports a single MPI rank.

See [INSTALLATION.md](INSTALLATION.md#10-optional-cgnshdf5-output) for detailed configuration and troubleshooting.

## 8. Troubleshooting

**Mechanism or thermo file missing:** your current working directory probably lacks `task_setup/` and its `chemical_mechanisms/` and `thermophysical_data/` subdirectories. Check where you launch the process from, not only where its executable resides.

**`cmake -E` fails during 1D generation:** CMake needs to be on `PATH` at **runtime**, including when the generator is invoked manually outside a developer prompt.

**`computing_module` cannot open `problem_setup.log` or `task_setup/*.inf`:** run from the generated **leaf** directory, not the executable directory or the parent of multiple cases.

**Only one CPU core is used:** supply `--num_threads=N`; verify the binary was built with `NRG_ENABLE_OPENMP=ON`.

**No CGNS output / unsupported format error:** verify both `NRG_ENABLE_CGNS=ON` in the solver build and `save_format='cgns'` in the generated case configuration.

**QSS2/CVODE not selected:** check the actual generated `solver_setup.inf` rather than assuming that a different CMake build automatically selects a chemistry backend.

**Post-processing parsing fails:** `proc*.dat` files now include Tecplot-style `TITLE=` and `VARIABLES=` lines. Skip or parse those records before interpreting numeric rows.

**Old documentation claims a Windows-only 1D interface or a single near-wall 17% case:** those claims describe an earlier revision. The current generator uses portable filesystem operations and defaults to three 8% H₂ counter-flow cases.

---

For a broader architecture overview, see [README.md](README.md). For optional dependency flags, compiler setups, and build issues, use [INSTALLATION.md](INSTALLATION.md).
