FCOMP	= mpif90

#FLAGS 	= -O0 -g -traceback -fpp 	-check all -fpe:0 -warn all	# Debug Sequential
#FLAGS 	= -O0 -g -traceback -fpp -DOMP	-check all -fpe:0 -warn all	# Debug OMP
FLAGS 	= -O0 -g -traceback -fpp -Dmpi	-check all -fpe:0 -warn all	# Debug MPI

#FLAGS	= -O3 -g -traceback -fpp 		   -fpe:0		# Release Sequential
#FLAGS 	= -O3 -g -traceback -fpp -DOMP	-qopenmp   -fpe:0		# Release OMP
#FLAGS	= -O3 -g -traceback -fpp -Dmpi		   -fpe:0		# Release MPI

TARGET	= package_library.a

vpath %.f90 src

library: all
	ar rc $(TARGET) *.o

all: kind_parameters.o global_data.o computational_domain_class.o computational_mesh_class.o boundary_type_class.o boundary_conditions_class.o 
all: chemical_properties_class.o field_scalar_class.o field_tensor_class.o field_vector_class.o field_pointers.o thermophysical_properties_class.o 
all: mpi_communications_class.o data_manager_class.o data_save_class.o data_io_class.o post_processor_operation_class.o post_processor_class.o post_processor_manager_class.o
all: solver_options_class.o

kind_parameters.o			: kind_parameters.f90
global_data.o				: global_data.f90
computational_domain_class.o		: computational_domain_class.f90
computational_mesh_class.o		: computational_mesh_class.f90
boundary_type_class.o			: boundary_type_class.f90	
boundary_conditions_class.o		: boundary_conditions_class.f90
chemical_properties_class.o		: chemical_properties_class.f90
field_scalar_class.o			: field_scalar_class.f90
field_tensor_class.o			: field_tensor_class.f90
field_vector_class.o			: field_vector_class.f90
field_pointers.o			: field_pointers.f90
thermophysical_properties_class.o	: thermophysical_properties_class.f90
mpi_communications_class.o		: mpi_communications_class.f90
data_manager_class.o			: data_manager_class.f90
data_save_class.o			: data_save_class.f90
data_io_class.o				: data_io_class.f90
post_processor_operation_class.o	: post_processor_operation_class.f90
post_processor_class.o			: post_processor_class.f90
post_processor_manager_class.o		: post_processor_manager_class.f90
solver_options_class.o			: solver_options_class.f90

%.o: %.f90
	$(FCOMP) $(FLAGS) -c $<

clean:
	rm -f *.o 