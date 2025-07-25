# NRG project

Welcome to Numerical Reactive Gas-dynamics software package. 

## Getting Started

These instructions will get you a copy of the project up and running on your local machine for development and testing purposes.

### Prerequisites

Fortran compiler with 2003 standard features support.
### Installing

Installation guide in pdf format can be found in docs folder. 

CMakeBuild tool installation:
mkdir build
cmake .. -G"Unix Makefiles" -DCMAKE_Fortran_COMPILER=gfortran -DCMAKE_BUILD_TYPE=DEBUG
cmake --build . --target computing_module