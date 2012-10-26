#! /usr/local/bin/python

import sys, numpy

if len(sys.argv) == 1:

  print "usage: " + sys.argv[0] + " file.rix"

inputFile = open(sys.argv[1], 'r')

nx = 100   # r
ny = 30    # theta

array = numpy.zeros((nx,ny),dtype=numpy.float)

i = 0   # will go from 0-> nx - 1


for line in inputFile.readlines():
    j = 0  # will go from 0-> ny - 1
    for element in line.split():
        value = float(element)
        array[i][j] = value
        j+=1
    i+=1 

inputFile.close()

outputFile = open('line.dat','w')
outputFile.write('# r, 0 deg, 45 deg, 90 deg\n')

theta0 = 0
theta1 = ny/4 
theta2 = ny/2

for i in range(nx - 1):

   outputFile.write(str(i) + ' ' + str(array[i][theta0]) + ' ' + str(array[i][theta1]) + ' ' + str(array[i][theta2]) + '\n')

outputFile.close()
