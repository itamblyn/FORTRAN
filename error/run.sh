#!/bin/bash

make
./generator.py

rm -f result.dat

x=10
while [ $x -le 100 ]
do
  echo "nsteps $x "
  echo -n "$x " >> result.dat
  echo "data.dat" > error.in
  echo "$x " >> error.in
  echo "$x " >> error.in
  x=$(( $x + 20 ))
  ./exe.x < error.in
  rm error.in
  R --no-save -q < fit.R > fit.out
  grep -A 1 Intercept fit.out | awk '/I /{print $2,$3}' | head -1 | awk '{printf $0}' >> result.dat
  grep -A 1 Intercept fit.out | awk '/I /{print " ",$2,$3}' | tail -1 >> result.dat
  rm fit.out
done
#gnuplot plot.plt
