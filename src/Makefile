install:
	gfortran flama.f90 -o flama -ffree-form -march=native
all:
	gfortran flama.f90 -o flama -ffree-form -march=native
	cp peros_initial.lib peros.lib ; ./flama < input
clean:
	rm -rf flama *.o *.mod
