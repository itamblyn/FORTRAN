cmolang: cmolang.f90
	ifort vectors.f90 cmolang.f90 -o cmol.x -O3 -lm -stand f95
molang: molang.f90
	ifort vectors.f90 molang.f90 -o mol.x -O3 -lm -stand f95
clean:
	rm *.x *.mod
