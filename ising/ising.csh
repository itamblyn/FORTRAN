#!/bin/csh
#
F90 -c -g ising.f90 >& compiler.txt
if ( $status != 0 ) then
  echo "Errors compiling ising.f90"
  exit
endif
rm compiler.txt
#
F90 ising.o
if ( $status != 0 ) then
  echo "Errors linking and loading ising.o"
  exit
endif
rm ising.o
#
chmod ugo+x a.out
mv a.out ~/bin/$ARCH/ising
#
echo "The program has been installed as ~/bin/$ARCH/ising."
