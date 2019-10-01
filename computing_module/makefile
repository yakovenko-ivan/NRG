FCOMP	= mpif90

FLIB	= release_MPI

#FLAGS 	= -O0 -I$(FLIB) -g -traceback -fpp 		-check all -fpe:0 	# Debug Sequential
#FLAGS 	= -O0 -I$(FLIB) -g -traceback -fpp -DOMP	-check all -fpe:0 	# Debug OMP
#FLAGS 	= -O0 -I$(FLIB) -g -traceback -fpp -Dmpi	-check all -fpe:0 	# Debug MPI

#FLAGS	= -O3 -I$(FLIB)	-fpp 			-fpe:0				# Release Sequential
#FLAGS 	= -O3 -I$(FLIB)	-fpp -DOMP -qopenmp	-fpe:0				# Release OMP
FLAGS	= -O3 -I$(FLIB)	-fpp -Dmpi		-fpe:0				# Release MPI

TARGET	= computing_module

SLATEC 	= $(wildcard src/third_party/slatec/*.f90)
SLATEC_OBJ = $(addprefix src/third_party/slatec/,$(notdir $(SLATEC:.f90=.o)))

vpath %.f90 src/current_build
vpath %.f90 src/test_versions

computing_module: all
	$(FCOMP) $(FLAGS) -o $(TARGET) *.o $(FLIB)/package_library.a

all: $(SLATEC_OBJ)
all: particles_solver.o
all: droplets_solver.o
all: lagrangian_droplets_solver.o
all: viscosity_solver.o
all: fickean_diffusion_solver.o
all: chemical_kinetics_solver.o
all: table_approximated_real_gas.o 
all: fourier_heat_transfer_solver.o
all: coarse_particles_method.o
all: cpm_solver.o
all: cpm_low_mach_solver.o
all: fds_low_mach_solver.o
all: cabaret_solver.o
all: main.o


particles_solver.o				: particles_solver.f90
droplets_solver.o				: droplets_solver.f90
lagrangian_droplets_solver.o			: lagrangian_droplets_solver.f90
viscosity_solver.o				: viscosity_solver.f90
fickean_diffusion_solver.o			: fickean_diffusion_solver.f90
chemical_kinetics_solver.o 			: chemical_kinetics_solver.f90 
table_approximated_real_gas.o			: table_approximated_real_gas.f90
fourier_heat_transfer_solver.o			: fourier_heat_transfer_solver.f90
coarse_particles_method.o			: coarse_particles_method.f90
cpm_solver.o					: cpm_solver.f90
cpm_low_mach_solver.o				: cpm_low_mach_solver.f90
fds_low_mach_solver.o				: fds_low_mach_solver.f90
cabaret_solver.o				: cabaret_solver.f90
main.o						: main.f90

%.o: %.f90
	$(FCOMP) $(FLAGS) -c $<

src/third_party/slatec/%.o: src/third_party/slatec/%.f90
	$(FCOMP) $(FLAGS) -c $<

clean:
	rm -f *.o *.mod
