FCOMP	= ifort
FLIB	= debug_seq
FLAGS	= -g -traceback -O0 -I$(FLIB) -check all
TARGET	= tecplot_merger

vpath %.f90 tecplot_merger_program

tmerger: all
	$(FCOMP) $(FLAGS) -o $(TARGET) *.o $(FLIB)/package_library.a

all: tecplot_merger.o

tecplot_merger.o	: tecplot_merger.f90

%.o: %.f90
	$(FCOMP) $(FLAGS) -c $<

clean:
	rm -f *.o *.mod