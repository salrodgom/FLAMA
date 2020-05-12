FFLAGS = -ffree-form -march=native -O3 -fopenmp
PROGRAM = flama.exe
install:
	gfortran flama.f90 -o $(PROGRAM) $(FFLAGS)
all:
	gfortran flama.f90 -o $(PROGRAM) $(FFLAGS)
	cp peros_initial.lib peros.lib ; ./$(PROGRAM) < input
clean:
	rm -rf $(PROGRAM) *.o *.mod
