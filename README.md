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

