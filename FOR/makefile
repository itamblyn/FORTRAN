cbn: cbn_from_xyz_cpmd.f90
	ifort cbn_from_xyz_cpmd.f90 -o cbn_from_xyz_cpmd.x -O3 -stand f95

msd: msd_com.f90
	ifort msd_com.f90 -o msd.x -O3 -stand f95

unwrap: unwrap_PBC.f90
	ifort unwrap_PBC.f90 -o unwrap_PBC.x -O3 -stand f95

small: small_displacement_check.f90 
	ifort small_displacement_check.f90 -o small_displacement_check.x -O3 -stand f95

stepsize: stepsize_check.f90
	ifort stepsize_check.f90 -o stepsize_check.x -O3 -stand f95

clean:
	rm *.x
