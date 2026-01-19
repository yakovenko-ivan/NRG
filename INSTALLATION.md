# NRG: Complete Installation Guide

This guide provides step-by-step instructions to build the NRG software from source on **Windows**, **Linux**, and **macOS**. Multiple compiler and toolchain options are provided for each platform.

**Choose your operating system:**
- [Installation on Windows](#installation-on-windows)
- [Installation on Linux](#installation-on-linux)
- [Installation on macOS](#installation-on-macos)

---

## Prerequisites

All installations require:
*   **Git**: To clone the repository. [Download Git](https://git-scm.com/downloads)
*   **CMake 3.12 or higher**: The build system generator. [Download CMake](https://cmake.org/download/)

Platform-specific requirements are detailed below.

---

## Installation on Windows

### **Option A: Native Windows with Visual Studio & Intel oneAPI (Recommended for Performance)**

This method uses the native Microsoft toolchain (`MSBuild`) with the Intel `ifx` compiler from the command line.

1.  **Install Prerequisites:**
    *   **Visual Studio 2022**: Install the free **Community Edition**. During setup, select the **"Desktop development with C++"** workload. This provides the necessary build tools (`MSBuild`).
    *   **Intel oneAPI**: Download and install both the **Base Toolkit** and **HPC Toolkit** (free). This provides the `ifx` Fortran compiler.

2.  **Build from the Developer Command Prompt:**
    Open the **"x64 Native Tools Command Prompt for VS 2022"** from the Start Menu. This is crucial as it sets the correct environment variables.
    ```cmd
    git clone https://github.com/yakovenko-ivan/NRG.git
    cd NRG
    mkdir build && cd build
    cmake .. -G "Visual Studio 17 2022" -A x64 -DCMAKE_Fortran_COMPILER=ifx
    cmake --build . --target computing_module --config Release
    ```
    The executable `computing_module.exe` will be in `build\Release\`.

### **Option B: GNU Toolchain with gfortran (Recommended for Simplicity)**

This method uses GNU tools (`gfortran`, `make`) via the Equation.com installer and can be used from a standard terminal or Visual Studio Code.

1.  **Install Prerequisites:**
    *   **gfortran for Windows**: Download and run the latest 64-bit installer from [Equation.com](https://www.equation.com/servlet/equation.cmd?fa=fortran).
    *   **(Optional) Visual Studio Code**: A lightweight editor with excellent terminal integration.

2.  **Build from the Command Prompt or VS Code Terminal:**
    Open a standard **Command Prompt** or the terminal in VS Code.
    ```cmd
    git clone https://github.com/yakovenko-ivan/NRG.git
    cd NRG
    mkdir build && cd build
    cmake .. -G "MinGW Makefiles"
    cmake --build . --target computing_module
    ```
    The executable `computing_module.exe` will be in the `build\` directory.

---

## Installation on Linux

### **Option A: Using Intel oneAPI `ifx` (Recommended for Performance)**

1.  **Install Prerequisites:**
    *   **Intel oneAPI HPC Toolkit**: Follow the [official instructions](https://www.intel.com/content/www/us/en/developer/tools/oneapi/hpc-toolkit.html) for your distribution (e.g., using the APT or YUM repository).
    *   **CMake & `make`**:
        ```bash
        # Ubuntu/Debian
        sudo apt install cmake make
        # Fedora/RHEL
        sudo dnf install cmake make
        ```

2.  **Set Up the Compiler Environment & Build:**
    Before building, source the oneAPI environment variables. You can add this line to your `~/.bashrc`.
    ```bash
    source /opt/intel/oneapi/setvars.sh
    ```
    Then, clone and build NRG:
    ```bash
    git clone https://github.com/yakovenko-ivan/NRG.git
    cd NRG
    mkdir build && cd build
    cmake .. -DCMAKE_Fortran_COMPILER=ifx
    cmake --build . --target computing_module
    ```
    The `computing_module` executable will be in the `build/` directory.

### **Option B: Using GNU `gfortran` (Common Default)**

1.  **Install Prerequisites:**
    ```bash
    # Ubuntu/Debian
    sudo apt install gfortran cmake make
    # Fedora/RHER
    sudo dnf install gcc-gfortran cmake make
    ```

2.  **Build NRG:**
    ```bash
    git clone https://github.com/yakovenko-ivan/NRG.git
    cd NRG
    mkdir build && cd build
    cmake .. -DCMAKE_BUILD_TYPE=Release
    cmake --build . --target computing_module
    ```
    The `computing_module` executable will be in the `build/` directory.

**Visual Studio Code** can be installed via Snap or your distribution's package manager for a convenient development environment.

---

## Installation on macOS

**Primary Recommended Path: Using GNU `gfortran`**

The Intel `ifx` compiler is **not available** for macOS. The classic `ifort` compiler is deprecated and not recommended for new projects.

1.  **Install Prerequisites:**
    The easiest method is via **Homebrew**. Install [Homebrew](https://brew.sh/) if you don't have it, then run:
    ```bash
    brew install gcc cmake
    # Ensure Xcode Command Line Tools are installed for 'make':
    xcode-select --install
    ```
    This installs `gfortran` (as part of the `gcc` formula) and CMake.

2.  **Build NRG:**
    ```bash
    git clone https://github.com/yakovenko-ivan/NRG.git
    cd NRG
    mkdir build && cd build
    cmake .. -DCMAKE_Fortran_COMPILER=gfortran -DCMAKE_BUILD_TYPE=Release
    cmake --build . --target computing_module
    ```
    The `computing_module` executable will be in the `build/` directory.

**Visual Studio Code** can be installed via Homebrew (`brew install --cask visual-studio-code`) or downloaded from its website.

---

## Troubleshooting Common Issues

*   **CMake Error: "No CMAKE_Fortran_COMPILER could be found."**
    *   **Windows (Option A)**: Ensure you are using the **"x64 Native Tools Command Prompt"**, not a regular one.
    *   **Linux (Option A)**: Ensure you have sourced the oneAPI environment (`source /opt/intel/oneapi/setvars.sh`).
    *   **General**: Verify your compiler is installed and in your system's PATH (`ifx --version` or `gfortran --version`).

*   **Build fails with compilation errors**
    *   Try building in `Debug` mode first, which may have fewer aggressive optimizations: `cmake .. -DCMAKE_BUILD_TYPE=Debug`

*   **"Permission denied" error on Linux/macOS**
    ```bash
    chmod +x build/computing_module
    ```

---

## Next Steps

After a successful build, proceed to the [Tutorial](https://github.com/yakovenko-ivan/NRG/wiki/Tutorial-1:-First-Simulation) to run your first simulation.

**Need Help?**
Please search for your error message in the [GitHub Issues](https://github.com/yakovenko-ivan/NRG/issues) or open a new one if your problem isn't documented.