#!/bin/bash

/usr/local/analysis/fortran/SDF/python/sample.py all_ions_and_molecules.hist all_ions_and_molecules.rix all_ions_and_molecules.crix &
/usr/local/analysis/fortran/SDF/python/sample.py all_ions.hist all_ions.rix all_ions.crix &
/usr/local/analysis/fortran/SDF/python/sample.py all_molecules.hist all_molecules.rix all_molecules.crix &
/usr/local/analysis/fortran/SDF/python/sample.py nearest_ions_and_molecules.hist nearest_ions_and_molecules.rix nearest_ions_and_molecules.crix &
/usr/local/analysis/fortran/SDF/python/sample.py nearest_ions.hist nearest_ions.rix nearest_ions.crix &
/usr/local/analysis/fortran/SDF/python/sample.py nearest_molecules.hist nearest_molecules.rix nearest_molecules.crix &

wait

/usr/local/analysis/fortran/SDF/python/pypng.py all_ions_and_molecules.crix all_ions_and_molecules.png &
/usr/local/analysis/fortran/SDF/python/pypng.py all_ions.crix all_ions.png &
/usr/local/analysis/fortran/SDF/python/pypng.py all_molecules.crix all_molecules.png &
/usr/local/analysis/fortran/SDF/python/pypng.py nearest_ions_and_molecules.crix nearest_ions_and_molecules.png &
/usr/local/analysis/fortran/SDF/python/pypng.py nearest_ions.crix nearest_ions.png &
/usr/local/analysis/fortran/SDF/python/pypng.py nearest_molecules.crix nearest_molecules.png &

wait
