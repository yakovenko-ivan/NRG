FCOMP	= ifort

FLIB	= release_OMP

#FLAGS 	= -O0 -I$(FLIB) -g -traceback -fpp 		-check all -fpe:0 -warn all	# Debug Sequential
#FLAGS 	= -O0 -I$(FLIB) -g -traceback -fpp -DOMP	-check all -fpe:0 -warn all	# Debug OMP
#FLAGS 	= -O0 -I$(FLIB) -g -traceback -fpp -Dmpi	-check all -fpe:0 -warn all	# Debug MPI

#FLAGS	= -O3 -I$(FLIB)	-fpp 			-fpe:0				# Release Sequential
FLAGS 	= -O3 -I$(FLIB)	-fpp -DOMP -qopenmp	-fpe:0				# Release OMP
#FLAGS	= -O3 -I$(FLIB)	-fpp -Dmpi		-fpe:0				# Release MPI

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
all: fickean_diffusion_solver_new_multicoord.o
all: viscosity_solver_multicoord.o
all: fourier_heat_transfer_solver_multicoord.o
all: chemical_kinetics_solver_class.o
all: table_approximated_real_gas_class.o 
all: coarse_particles_method_multicoord.o
all: cpm_solver_class.o
all: cpm_low_mach_solver_class.o
all: fds_low_mach_solver_multidim_parallel.o
all: cabaret_solver_class_alternative_multicoord.o
all: main.o


particles_solver.o				: particles_solver.f90
droplets_solver.o				: droplets_solver.f90
viscosity_solver_multicoord.o			: viscosity_solver_multicoord.f90
fickean_diffusion_solver_new_multicoord.o	: fickean_diffusion_solver_new_multicoord.f90
chemical_kinetics_solver_class.o 		: chemical_kinetics_solver_class.f90 
table_approximated_real_gas_class.o		: table_approximated_real_gas_class.f90
fourier_heat_transfer_solver_multicoord.o	: fourier_heat_transfer_solver_multicoord.f90
coarse_particles_method_multicoord.o		: coarse_particles_method_multicoord.f90
cpm_solver_class.o				: cpm_solver_class.f90
cpm_low_mach_solver_class.o			: cpm_low_mach_solver_class.f90
fds_low_mach_solver_multidim_parallel.o		: fds_low_mach_solver_multidim_parallel.f90
cabaret_solver_class_alternative_multicoord.o	: cabaret_solver_class_alternative_multicoord.f90
main.o						: main.f90

%.o: %.f90
	$(FCOMP) $(FLAGS) -c $<

src/third_party/slatec/%.o: src/third_party/slatec/%.f90
	$(FCOMP) $(FLAGS) -c $<

clean:
	rm -f *.o *.mod
