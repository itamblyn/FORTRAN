all: center_of_mass.f90
	ifort center_of_mass.f90 -o exe.x -O3 -lm -static -stand f95
test: test.f90
	ifort test.f90 -o exe.x -O3 -lm -static -stand f95
clean:
	rm exe.x
