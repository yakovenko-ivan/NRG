# NRG Installation and Build Guide

This guide describes the build system in the current `refactor/core-solvers-validation` branch of NRG. It covers the core solver, selectable problem interfaces and utilities, OpenMP, and the optional CGNS/HDF5 output backend.

The commands below are written for an **out-of-source build**. Generated modules, libraries and executables should remain inside the CMake build directory.

## 1. Supported build configuration

The root `CMakeLists.txt` currently requires **CMake 3.15 or newer** and accepts these Fortran compiler IDs:

- `IntelLLVM` — Intel oneAPI **`ifx`**;
- `GNU` — **`gfortran`**;
- `Intel` — classic **`ifort`**, retained as a compatibility path.

OpenMP is enabled by default.

### Recommended versions by configuration

| Configuration | Recommended CMake | Notes |
|---|---:|---|
| Windows + Visual Studio 2022 + `ifx` | **3.29+** | enables `-T fortran=ifx` for Visual Studio generators |
| Linux core build, `ifx`/`gfortran` | 3.15+ | newer CMake is still recommended |
| macOS core build, `gfortran` | 3.15+ | current classic interfaces may not be portable |
| Bundled CGNS/HDF5 | **3.20+** | additionally affected by the current branch issue described below |

## 2. Clone the branch

```bash
git clone https://github.com/yakovenko-ivan/NRG.git
cd NRG
git switch refactor/core-solvers-validation
```

To confirm the branch:

```bash
git branch --show-current
```

## 3. Important CMake behavior

### 3.1 Build configurations

For **single-configuration** generators such as Unix Makefiles or Ninja, choose the configuration at configure time:

```bash
-DCMAKE_BUILD_TYPE=Release
```

If `CMAKE_BUILD_TYPE` is omitted, NRG currently defaults to **Debug**.

For **multi-configuration** generators such as Visual Studio, choose the configuration at build time:

```cmd
cmake --build build --config Release
```

### 3.2 Output directories

The current root CMake file deliberately places products in configuration-specific subdirectories for **both** generator families:

```text
build/bin/<Config>/
build/lib/<Config>/
build/mod_files/
```

Examples:

```text
build/bin/Release/computing_module
build/bin/Release/computing_module.exe
build/bin/Debug/computing_module.exe
```

This differs from older NRG documentation that used `build/bin/` directly for single-config generators.

### 3.3 OpenMP

OpenMP is enabled by default:

```text
NRG_ENABLE_OPENMP=ON
```

Disable it with:

```bash
-DNRG_ENABLE_OPENMP=OFF
```

When enabled, CMake requires `OpenMP::OpenMP_Fortran` and defines the `OMP` preprocessing symbol.

### 3.4 Windows runtime linkage

On Windows, NRG defaults to:

```text
NRG_STATIC_RUNTIME=ON
```

For Intel/MSVC ABI builds this selects the static multithreaded compiler runtime (`/MT`, or `/MTd` in Debug). Set:

```cmd
-DNRG_STATIC_RUNTIME=OFF
```

for `/MD`/`/MDd` linkage.

**Important:** when OpenMP is enabled with Intel Fortran, Intel's OpenMP runtime remains dynamic (`libiomp5md.dll`) even if `NRG_STATIC_RUNTIME=ON`.

## 4. Windows 11: Visual Studio 2022 + Intel oneAPI `ifx`

This is the recommended configuration for the current reference problem interfaces.

### 4.1 Install prerequisites

Install:

1. **Visual Studio 2022** with the **Desktop development with C++** workload;
2. Intel **oneAPI HPC Toolkit**, including Intel Fortran Compiler and Visual Studio integration;
3. CMake **3.29 or newer**;
4. Git.

Open an **Intel oneAPI command prompt for Visual Studio 2022** (or another terminal where both environments are initialized).

Verify:

```cmd
where ifx
ifx --version
cmake --version
```

### 4.2 Configure with the Visual Studio Fortran toolset

Use CMake's Visual Studio toolset selection:

```cmd
cmake -S . -B build ^
  -G "Visual Studio 17 2022" ^
  -A x64 ^
  -T fortran=ifx
```

Then build the solver:

```cmd
cmake --build build --target computing_module --config Release --parallel
```

Build the currently selected problem interface:

```cmd
cmake --build build --target package_interface --config Release --parallel
```

Build the currently selected utility:

```cmd
cmake --build build --target package_utilities --config Release --parallel
```

Build every configured target:

```cmd
cmake --build build --config Release --parallel
```

### 4.3 Why `-T fortran=ifx` is preferred

With Visual Studio generators, do not normally force:

```text
-DCMAKE_Fortran_COMPILER=ifx
```

CMake 3.29 introduced the Visual Studio generator-toolset field `fortran=ifx`. The reliable syntax is:

```cmd
-T fortran=ifx
```

This makes the generated Visual Studio Fortran projects use Intel `ifx` and avoids the common state where CMake identifies IntelLLVM during probing but later reports that `ifx` is not a full path or cannot be found.

When changing from an older compiler-selection method, delete the old build directory first because CMake caches the compiler and generator toolset:

```cmd
rmdir /S /Q build
```

then configure again with `-T fortran=ifx`.

## 5. Linux

### 5.1 Intel oneAPI `ifx`

Initialize oneAPI:

```bash
source /opt/intel/oneapi/setvars.sh
```

Configure and build the core solver:

```bash
cmake -S . -B build \
  -DCMAKE_Fortran_COMPILER=ifx \
  -DCMAKE_BUILD_TYPE=Release

cmake --build build --target computing_module --parallel
```

The executable is expected at:

```text
build/bin/Release/computing_module
```

### 5.2 GNU `gfortran`

Ubuntu/Debian example:

```bash
sudo apt update
sudo apt install git gfortran cmake make
```

Configure and build:

```bash
cmake -S . -B build \
  -DCMAKE_Fortran_COMPILER=gfortran \
  -DCMAKE_BUILD_TYPE=Release

cmake --build build --target computing_module --parallel
```

### Current interface portability warning

The root build supports GNU Fortran, but the current default interface

```text
package_interface/src/tests/classic_tests/1D_laminar_velocity.f90
```

contains unconditional Intel/Windows-specific code, including:

```fortran
use ifport
```

and an `xcopy` command used to create the case directory. Therefore **the current 1D reference interface should not be presented as a Linux/gfortran quick start**. The core solver can still be built independently with the `computing_module` target.

## 6. macOS

Intel oneAPI `ifx` is not available for macOS. The intended core compiler path is Homebrew GNU Fortran:

```bash
brew install git cmake gcc
xcode-select --install
```

Configure:

```bash
cmake -S . -B build \
  -DCMAKE_Fortran_COMPILER=gfortran \
  -DCMAKE_BUILD_TYPE=Release

cmake --build build --target computing_module --parallel
```

If Homebrew installs a version-suffixed compiler (for example `gfortran-15`), pass that executable name instead.

As on Linux, the current `1D_laminar_velocity.f90` problem interface is not portable to macOS without source changes.

## 7. Selecting a package interface

The cache variable:

```text
PACKAGE_INTERFACE_SOURCE
```

contains a Fortran source path **relative to `package_interface/`**.

Current default:

```text
src/tests/classic_tests/1D_laminar_velocity.f90
```

Other interfaces currently present in the branch include:

```text
src/tests/classic_tests/0D_ignition_delay_agent.f90
src/tests/classic_tests/0D_ignition_delay_campaign.f90
src/tests/classic_tests/3D_droplet_evaporation.f90
```

Example:

```bash
cmake -S . -B build \
  -DCMAKE_BUILD_TYPE=Release \
  -DPACKAGE_INTERFACE_SOURCE=src/tests/classic_tests/3D_droplet_evaporation.f90

cmake --build build --target package_interface --parallel
```

The **target name** is always:

```text
package_interface
```

The **output executable name** is derived from the selected source basename. For the default source it is:

```text
package_interface_1D_laminar_velocity
```

### Interfaces with companion source files

The package-interface CMake file also defines:

```text
PACKAGE_INTERFACE_EXTRA_SOURCES
```

This is a semicolon-separated list of extra Fortran sources, each relative to `package_interface/`.

Example syntax:

```bash
-DPACKAGE_INTERFACE_SOURCE=src/my_case.f90 \
-DPACKAGE_INTERFACE_EXTRA_SOURCES="src/helper_a.f90;src/helper_b.f90"
```

## 8. Selecting a package utility

The corresponding cache variable is:

```text
PACKAGE_UTILITY_SOURCE
```

Current default:

```text
src/leading_point_velocity.f90
```

Other utilities in this branch include:

```text
src/output_merger.f90
src/tecplot_merger.f90
```

Example:

```bash
cmake -S . -B build \
  -DCMAKE_BUILD_TYPE=Release \
  -DPACKAGE_UTILITY_SOURCE=src/output_merger.f90

cmake --build build --target package_utilities --parallel
```

The target name is `package_utilities`; the executable name is derived from the selected source, for example `package_utilities_output_merger`.

## 9. CGNS/HDF5 output

### 9.1 Default state

The optional backend is **disabled by default**:

```text
NRG_ENABLE_CGNS=OFF
```

A case that requests:

```fortran
save_format = 'cgns'
```

cannot be run with a solver built this way. `computing_module` will stop with:

```text
Data save: CGNS requested but NRG was built with NRG_ENABLE_CGNS=OFF
```

### 9.2 Intended bundled dependency path

The branch contains a `FetchContent` configuration intended to download and build:

- HDF5 **1.14.6**;
- CGNS **4.5.2**.

The intended configuration is:

```bash
cmake -S . -B build \
  -DCMAKE_BUILD_TYPE=Release \
  -DNRG_ENABLE_CGNS=ON
```

This path requires CMake 3.20+ and network access during first configuration.

### 9.3 Current branch blocker in the bundled path

At commit:

```text
0270e5947d3ca2b80d005856c577188952170d99
```

`cmake/NRGCGNSDependencies.cmake` tries to configure:

```text
cmake/hdf5-buildtree/hdf5-config.cmake.in
```

but that template file is not present in the branch. Therefore the bundled CGNS/HDF5 path is **not currently self-contained** and is expected to fail during configuration until the missing file is restored.

This should be fixed in the repository rather than worked around in production documentation.

### 9.4 Using a system CGNS installation

If a compatible CGNS installation built with HDF5 is available, select it with:

```bash
-DNRG_ENABLE_CGNS=ON \
-DNRG_USE_SYSTEM_CGNS=ON
```

If necessary also provide one of:

```text
CMAKE_PREFIX_PATH
cgns_DIR
NRG_CGNS_INCLUDE_DIR
NRG_CGNS_LIBRARY
```

The current CGNS writer supports **one MPI rank only**.

### 9.5 Default Tecplot path

The Tecplot writer requires no external CGNS/HDF5 dependency and is the default field-output path for the reference tutorial. Problem interfaces select it with:

```fortran
save_format = 'tecplot'
```

A normal NRG build therefore does **not** require HDF5 or CGNS. Enable `NRG_ENABLE_CGNS` only for cases that explicitly request CGNS output.

## 10. Building without OpenMP

```bash
cmake -S . -B build-serial \
  -DCMAKE_Fortran_COMPILER=gfortran \
  -DCMAKE_BUILD_TYPE=Release \
  -DNRG_ENABLE_OPENMP=OFF

cmake --build build-serial --target computing_module --parallel
```

Use a separate build directory when changing compilers, generators, architectures, or OpenMP/runtime settings.

## 11. Compiler flags configured by NRG

### GNU Fortran

Common options include:

```text
-cpp
-ffree-line-length-512
```

Debug adds runtime checks, traceback/debug information, and floating-point traps.

### Intel Fortran

Windows builds use `/fpp` and `/traceback`; Unix-like Intel builds use `-fpp` and `-traceback`. Debug configurations enable stronger runtime checking. Release uses optimization with floating-point exception handling configured by the root build.

## 12. Inspecting a configuration

Useful commands:

```bash
cmake -S . -B build -LA
```

During configuration NRG prints, among other things:

- generator and configuration mode;
- operating system;
- Fortran compiler path/ID/version;
- OpenMP status;
- Windows runtime-linkage mode;
- selected package-interface source;
- selected utility source;
- runtime output directory;
- CGNS/HDF5 status.

## 13. Clean reconfiguration

Delete the build directory when changing:

- Fortran compiler;
- Visual Studio `fortran=` toolset;
- generator;
- architecture;
- OpenMP toolchain;
- static/dynamic runtime policy;
- major dependency mode.

Linux/macOS:

```bash
rm -rf build
```

Windows:

```cmd
rmdir /S /Q build
```

Changing only `PACKAGE_INTERFACE_SOURCE` or `PACKAGE_UTILITY_SOURCE` normally requires re-running CMake but does not require deleting the build directory.

## 14. Troubleshooting

### `ifx is not a full path and was not found in the PATH`

For Visual Studio 2022, use CMake 3.29+ and configure a **fresh build directory** with:

```cmd
cmake -S . -B build ^
  -G "Visual Studio 17 2022" ^
  -A x64 ^
  -T fortran=ifx
```

Do not carry forward a cache created with `-DCMAKE_Fortran_COMPILER=ifx`.

Also verify that Intel Fortran's Visual Studio integration is installed and that `where ifx` works in the active oneAPI/VS terminal.

### CMake selects `ifort` instead of `ifx`

Older Visual Studio generators may allow Intel integration to choose the classic compiler by default. Specify:

```cmd
-T fortran=ifx
```

with CMake 3.29+.

### OpenMP is not found

Either install/fix the compiler's OpenMP environment or configure:

```bash
-DNRG_ENABLE_OPENMP=OFF
```

### `libiomp5md.dll` is missing on Windows

Intel OpenMP remains dynamically linked even when `NRG_STATIC_RUNTIME=ON`. Run in an initialized oneAPI environment or install/provide the Intel OpenMP runtime.

### GNU build fails in `1D_laminar_velocity.f90` at `use ifport`

This is a current source portability limitation of that interface, not a root CMake compiler-detection problem. Build `computing_module` alone or port/select another interface.

### The 1D interface reports that `xcopy` is missing

The current 1D interface uses a Windows `xcopy` command unconditionally. It is therefore not currently a portable Linux/macOS case generator.

### `CGNS requested but NRG was built with NRG_ENABLE_CGNS=OFF`

The generated case explicitly requests CGNS output, but the executable was built without the optional backend. Reconfigure and rebuild with `-DNRG_ENABLE_CGNS=ON`, or use a Tecplot-configured case.

### Bundled CGNS configuration fails around `hdf5-config.cmake.in`

This is a known issue at the current branch tip: the referenced compatibility template is missing from the repository. Use system CGNS or Tecplot until the repository is corrected.

### Executable is not in `build/bin/`

The current build always adds a configuration subdirectory. Check:

```text
build/bin/Release/
build/bin/Debug/
build/bin/RelWithDebInfo/
```

### Release build is unexpectedly slow

For Unix Makefiles/Ninja, ensure the configure command contained:

```bash
-DCMAKE_BUILD_TYPE=Release
```

For Visual Studio, build with:

```cmd
--config Release
```

## 15. No installation target yet

The current branch does not define a `cmake --install` packaging/install layout. NRG is normally used directly from the CMake build tree, with problem directories stored separately or generated under a chosen working directory.

After a successful build, continue with [TUTORIAL.md](TUTORIAL.md).
