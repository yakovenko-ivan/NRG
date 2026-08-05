# NRG Installation and Build Guide

This guide describes how to configure and build NRG from source on **Windows**, **Linux**, and **macOS**. It matches the options and directory layout defined by the current root `CMakeLists.txt`.

## Supported build configuration

NRG currently supports:

- **CMake 3.13 or newer**
- **GNU Fortran** (`gfortran`)
- **Intel oneAPI Fortran** (`ifx`)
- **Intel classic Fortran** (`ifort`) for compatibility
- Optional **OpenMP**, enabled by default
- Single-configuration generators such as Unix Makefiles and Ninja
- Multi-configuration generators such as Visual Studio

Other Fortran compiler families are rejected by the root CMake configuration.

## Prerequisites

All platforms require:

- [Git](https://git-scm.com/downloads)
- [CMake](https://cmake.org/download/) 3.13 or newer
- A supported Fortran compiler
- A build backend appropriate for the selected CMake generator, such as Visual Studio/MSBuild, Make, or Ninja

Check the installed tools before configuring:

```bash
git --version
cmake --version
gfortran --version
```

For Intel oneAPI:

```bash
ifx --version
```

## Important CMake behavior

### Build type

For a **single-configuration generator**, select the build type during configuration:

```bash
-DCMAKE_BUILD_TYPE=Release
```

Supported values are:

- `Debug`
- `Release`
- `RelWithDebInfo`
- `MinSizeRel`

When `CMAKE_BUILD_TYPE` is omitted for a single-configuration generator, NRG defaults to **Debug**.

For a **multi-configuration generator**, such as Visual Studio, select the configuration during the build:

```bash
cmake --build build --config Release
```

`CMAKE_BUILD_TYPE` is not used by multi-configuration generators.

### OpenMP

OpenMP is enabled by default:

```bash
-DNRG_ENABLE_OPENMP=ON
```

When enabled, CMake must find the Fortran OpenMP package. To build without OpenMP:

```bash
-DNRG_ENABLE_OPENMP=OFF
```

The build system defines the `OMP` preprocessing symbol only when OpenMP is enabled.

## Installation on Windows

### Option A: Visual Studio 2022 and Intel oneAPI `ifx`

This is the recommended native Windows configuration for production calculations.

Install:

1. **Visual Studio 2022** with the **Desktop development with C++** workload.
2. **Intel oneAPI Base Toolkit**.
3. **Intel oneAPI HPC Toolkit**, which provides `ifx`.
4. CMake and Git.

Open a terminal in which the Visual Studio and oneAPI environments are initialized, then run:

```cmd
git clone https://github.com/yakovenko-ivan/NRG.git
cd NRG

cmake -S . -B build ^
  -G "Visual Studio 17 2022" ^
  -A x64 ^
  -DCMAKE_Fortran_COMPILER=ifx

cmake --build build --target computing_module --config Release --parallel
```

The Visual Studio generator is multi-configuration, so `Release` is selected with `--config Release`.

To build all configured targets:

```cmd
cmake --build build --config Release --parallel
```

### Option B: GNU `gfortran`

Install a 64-bit GNU Fortran distribution together with a compatible build backend such as Make or Ninja.

Example using a single-configuration generator:

```cmd
git clone https://github.com/yakovenko-ivan/NRG.git
cd NRG

cmake -S . -B build ^
  -G "Unix Makefiles" ^
  -DCMAKE_Fortran_COMPILER=gfortran ^
  -DCMAKE_BUILD_TYPE=Release

cmake --build build --target computing_module --parallel
```

Use a generator that matches the build tool installed on the system. For example, use `Ninja` only when the `ninja` executable is available.

## Installation on Linux

### Option A: Intel oneAPI `ifx`

Install the Intel oneAPI HPC Toolkit and initialize its environment:

```bash
source /opt/intel/oneapi/setvars.sh
```

Install CMake, Git, and a build backend if they are not already present.

Ubuntu/Debian:

```bash
sudo apt update
sudo apt install git cmake make
```

Fedora/RHEL:

```bash
sudo dnf install git cmake make
```

Configure and build:

```bash
git clone https://github.com/yakovenko-ivan/NRG.git
cd NRG

cmake -S . -B build \
  -DCMAKE_Fortran_COMPILER=ifx \
  -DCMAKE_BUILD_TYPE=Release

cmake --build build --target computing_module --parallel
```

### Option B: GNU `gfortran`

Ubuntu/Debian:

```bash
sudo apt update
sudo apt install git gfortran cmake make
```

Fedora/RHEL:

```bash
sudo dnf install git gcc-gfortran cmake make
```

Configure and build:

```bash
git clone https://github.com/yakovenko-ivan/NRG.git
cd NRG

cmake -S . -B build \
  -DCMAKE_Fortran_COMPILER=gfortran \
  -DCMAKE_BUILD_TYPE=Release

cmake --build build --target computing_module --parallel
```

## Installation on macOS

The supported compiler path on macOS is GNU `gfortran`. Intel `ifx` is not available for macOS.

Install the required tools with Homebrew:

```bash
brew install git cmake gcc
xcode-select --install
```

Configure and build:

```bash
git clone https://github.com/yakovenko-ivan/NRG.git
cd NRG

cmake -S . -B build \
  -DCMAKE_Fortran_COMPILER=gfortran \
  -DCMAKE_BUILD_TYPE=Release

cmake --build build --target computing_module --parallel
```

Homebrew may install a version-suffixed compiler executable. When `gfortran` is not found, use the name reported by:

```bash
brew list gcc | grep '/gfortran'
```

and pass that executable to `CMAKE_Fortran_COMPILER`.

## Building without OpenMP

OpenMP is enabled by default and is required during configuration unless explicitly disabled.

Example:

```bash
cmake -S . -B build-serial \
  -DCMAKE_Fortran_COMPILER=gfortran \
  -DCMAKE_BUILD_TYPE=Release \
  -DNRG_ENABLE_OPENMP=OFF

cmake --build build-serial --target computing_module --parallel
```

Using separate build directories such as `build` and `build-serial` prevents incompatible configurations from sharing one CMake cache.

## Available targets

The root project adds the following NRG components:

- `package_library`
- `package_interface`
- `package_utilities`
- `computing_module`

Build a specific target with:

```bash
cmake --build build --target computing_module
```

For a multi-configuration generator:

```bash
cmake --build build --target computing_module --config Release
```

Build all targets by omitting `--target`:

```bash
cmake --build build --parallel
```

## Selecting a package-interface problem

The root CMake configuration exposes the cache variable:

```text
PACKAGE_INTERFACE_SOURCE
```

Its value is a source path **relative to `package_interface/`**.

The default is:

```text
src/tests/classic_tests/1D_laminar_velocity.f90
```

To select another problem interface:

```bash
cmake -S . -B build \
  -DCMAKE_BUILD_TYPE=Release \
  -DPACKAGE_INTERFACE_SOURCE=src/tests/classic_tests/1D_laminar_velocity.f90

cmake --build build --target package_interface --parallel
```

With Visual Studio:

```cmd
cmake -S . -B build ^
  -G "Visual Studio 17 2022" ^
  -A x64 ^
  -DCMAKE_Fortran_COMPILER=ifx ^
  -DPACKAGE_INTERFACE_SOURCE=src/tests/classic_tests/1D_laminar_velocity.f90

cmake --build build --target package_interface --config Release --parallel
```

Changing `PACKAGE_INTERFACE_SOURCE` requires re-running the CMake configuration step, but normally does not require deleting the build directory.

## Selecting a package utility

The corresponding cache variable for the utilities target is:

```text
PACKAGE_UTILITY_SOURCE
```

Its value is a source path **relative to `package_utilities/`**.

The default is:

```text
src/leading_point_velocity.f90
```

Example:

```bash
cmake -S . -B build \
  -DCMAKE_BUILD_TYPE=Release \
  -DPACKAGE_UTILITY_SOURCE=src/leading_point_velocity.f90

cmake --build build --target package_utilities --parallel
```

Both selectable sources can be configured at once:

```bash
cmake -S . -B build \
  -DCMAKE_BUILD_TYPE=Release \
  -DPACKAGE_INTERFACE_SOURCE=src/tests/classic_tests/1D_laminar_velocity.f90 \
  -DPACKAGE_UTILITY_SOURCE=src/leading_point_velocity.f90
```

## Build output directories

The root CMake project keeps generated files outside the source tree:

| Output | Single-configuration generator | Multi-configuration generator |
|---|---|---|
| Executables | `build/bin/` | normally `build/bin/<Config>/` |
| Libraries | `build/lib/` | normally `build/lib/<Config>/` |
| Fortran module files | `build/mod_files/` | `build/mod_files/` |

Examples:

```text
build/bin/computing_module
build/bin/computing_module.exe
build/bin/Release/computing_module.exe
```

The exact executable names are defined by the component-level CMake files.

## Inspecting the configuration

During configuration, CMake reports:

- generator type;
- build type for single-configuration generators;
- operating system;
- Fortran compiler path;
- compiler family and version;
- OpenMP status;
- selected package-interface source;
- selected package-utility source.

To display cached options:

```bash
cmake -S . -B build -LA
```

## Reconfiguration and clean builds

CMake caches the compiler, generator, build type, and selected source files.

A new build directory is strongly recommended when changing:

- Fortran compiler;
- compiler family;
- CMake generator;
- 32-bit/64-bit architecture;
- OpenMP toolchain.

Example:

```bash
rm -rf build
cmake -S . -B build \
  -DCMAKE_Fortran_COMPILER=gfortran \
  -DCMAKE_BUILD_TYPE=Release
```

On Windows:

```cmd
rmdir /S /Q build
```

For a normal rebuild with the same configuration:

```bash
cmake --build build --target clean
cmake --build build --parallel
```

## Compiler behavior configured by NRG

The root CMake file supplies the required preprocessing and diagnostic options.

For GNU Fortran it enables:

- preprocessing with `-cpp`;
- free-form lines up to 512 characters;
- runtime checks, traceback, and floating-point traps in `Debug`.

For Intel Fortran it enables:

- preprocessing with `/fpp` on Windows or `-fpp` on Unix-like systems;
- traceback information;
- runtime checks in `Debug`;
- floating-point exception handling in all configurations.

These options are propagated through the `package_library` target to NRG executables.

## Troubleshooting

### CMake cannot find a Fortran compiler

Verify that the compiler is available in the active terminal:

```bash
gfortran --version
```

or:

```bash
ifx --version
```

For Intel oneAPI on Linux, initialize the environment:

```bash
source /opt/intel/oneapi/setvars.sh
```

On Windows, use a terminal initialized for both Visual Studio and Intel oneAPI.

### OpenMP is not found

OpenMP is enabled by default. Either install/configure the compiler's OpenMP runtime or configure a serial build:

```bash
cmake -S . -B build -DNRG_ENABLE_OPENMP=OFF
```

### CMake reports an unsupported compiler

The root project currently accepts only GNU, IntelLLVM, and Intel compiler IDs. Use `gfortran`, `ifx`, or the compatibility `ifort` path.

### The build unexpectedly runs slowly

A single-configuration build defaults to `Debug` when `CMAKE_BUILD_TYPE` is omitted. Reconfigure explicitly:

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
```

### `--config Release` has no effect

`--config` applies only to multi-configuration generators. For Unix Makefiles or Ninja, set:

```bash
-DCMAKE_BUILD_TYPE=Release
```

during configuration.

### `CMAKE_BUILD_TYPE` has no effect with Visual Studio

Visual Studio is multi-configuration. Select the build configuration with:

```cmd
cmake --build build --config Release
```

### CMake keeps using the previous compiler

The compiler is stored in the CMake cache. Delete the build directory and configure again.

### The requested build program is missing

Select a generator whose backend is installed. For example, Unix Makefiles require `make`, Ninja requires `ninja`, and Visual Studio generators require the corresponding Visual Studio/MSBuild installation.

## Next steps

After a successful build, continue with the project tutorial for generating a problem setup with `package_interface` and running the resulting case with `computing_module`.

For build problems, search the [NRG issue tracker](https://github.com/yakovenko-ivan/NRG/issues) before opening a new issue.
