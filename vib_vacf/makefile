all: vib_vacf.f90
	ifort vib_vacf.f90 -o exe.x -O3 -lfftw3f_threads -lfftw3f -lm -lpthread -openmp -stand f95
clean:
	rm *.x
