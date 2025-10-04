![Fortran](https://img.shields.io/badge/Fortran-2003-brightgreen?logo=fortran&labelColor=734f96)
![CMake](https://img.shields.io/badge/CMake-4.0.0-brightgreen?logo=cmake&labelColor=royalblue)
![GitHub repo size](https://img.shields.io/github/repo-size/yakovenko-ivan/NRG)
![Static Badge](https://img.shields.io/badge/stability-stable-brightgreen?style=flat)
![Static Badge](https://img.shields.io/badge/version-1.1.0-royalblue?style=flat)


# NRG project

Welcome to Numerical Reactive Gas-dynamics software package. 

## Getting Started

These instructions will get you a copy of the project up and running on your local machine for development and testing purposes.

### Prerequisites

Fortran compiler with 2003 standard features support.
### Installing

Installation guide in pdf format can be found in docs folder. 

**CMakeBuild tool installation (gfortran, makefile generator):**

`mkdir build`

`cd build`

`cmake .. -G"Unix Makefiles" -DCMAKE_Fortran_COMPILER=gfortran -DCMAKE_BUILD_TYPE=DEBUG`

`cmake --build . --target computing_module`


**CMakeBuild tool installation (ifx, visual studio generator):**

`mkdir build`

`cd build`

`cmake .. -G "Visual Studio 17 2022" -A x64 -DCMAKE_Fortran_COMPILER=ifx -DCMAKE_BUILD_TYPE=DEBUG`

`cmake --build . --target computing_module`

