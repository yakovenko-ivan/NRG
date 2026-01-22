![Fortran](https://img.shields.io/badge/Fortran-2003-brightgreen?logo=fortran&labelColor=734f96)
![CMake](https://img.shields.io/badge/CMake-4.0.0-brightgreen?logo=cmake&labelColor=royalblue)
![GitHub repo size](https://img.shields.io/github/repo-size/yakovenko-ivan/NRG)
![Static Badge](https://img.shields.io/badge/stability-stable-brightgreen?style=flat)
![Static Badge](https://img.shields.io/badge/version-1.1.0-royalblue?style=flat)

# NRG: Numerical Reactive Gas-dynamics

**NRG** is a high-performance, open-source Computational Fluid Dynamics (CFD) software package designed for high-fidelity simulations of **reactive flows** and **combustion phenomena**. Built with modern Fortran, it leverages efficient numerical algorithms for scalable performance on high-performance computing (HPC) systems.

This software is developed at the **Joint Institute for High Temperatures (JIHT) of the Russian Academy of Sciences (RAS)**.

## ✨ Key Features

*   **Reactive Flow Solver**: Simulates compressible flows with detailed chemical kinetics, species transport, and heat transfer.
*   **Modular & Extensible Architecture**: Clear separation between the core solver (`computing_module`), physics libraries (`package_library`), and problem setup (`package_interface`).
*   **Modern Fortran Codebase**: Written using Fortran 2003/2008 standards for clarity, modularity, and performance.
*   **Portable Build System**: Uses CMake for straightforward configuration on Linux, macOS, and Windows.
*   **Targeted for HPC**: Designed with scalability in mind for parallel computing environments.

## 🚀 Getting Started

### Prerequisites
*   A Fortran compiler (`gfortran` ≥ 9.0 or Intel `ifx`)
*   `CMake` ≥ 3.12
*   `make` (or your system's build tool)
