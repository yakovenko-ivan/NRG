# NRG project

Welcome to Numerical Reactive Gas-dynamics software package. 

## Getting Started

These instructions will get you a copy of the project up and running on your local machine for development and testing purposes.

### Prerequisites

Fortran compiler with 2003 standard features support.
### Installing

Installation guide in pdf format can be found in docs folder. 

### Building with **CMake**

#### Environment variables

|   Name    |                 Role                 |                               Values                                    |
| :-------- | :----------------------------------: | :---------------------------------------------------------------------: |
|    ENV    |      Specify build configuration     | DEBUG_OMP, DEBUG_MPI, RELEASE_SEQUENTIAL, RELEASE_OMP, DEBUG_SEQUENTIAL |
|    WIN    | Use windows or Unix like file system |                                 -                                       |

#### Build targets

- Console application `NRG-computing-module`

- Static library `NRG-package-library`

- Console application `NRG-package-utilities`

- Examples: `NRG-examples-<FILE_NAME_WITHOUT_EXT>`