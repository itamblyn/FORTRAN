all: psd.f90
	ifort psd.f90 -o exe.x -O3 -lfftw3f_threads -lfftw3f -lm -lpthread -openmp -stand f95
clean:
	rm *.x
