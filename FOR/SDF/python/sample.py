#! /usr/bin/env python

# this file now requires that hist files ions+molecules.hist are normalised to begin with. 
# note that ions.hist and molecules.hist are NOT normalised, but sum to ions+molecules.hist 

import sys

if len(sys.argv) == 1: 
     print 'usage: ' + sys.argv[0] + ' input.hist output.rix output.crix'

import numpy

r_min = 1E6
r_max = float(0)
theta_min = 1E6
theta_max = float(0)

HIST_sum = float(0)

#number_of_r_bins = 100 
number_of_r_bins = 99 
number_of_theta_bins = 30

inputFile_HIST = open (sys.argv[1],'r')
line = inputFile_HIST.readline()

r_step = float(line.split()[5])
theta_step = float(line.split()[7])*(numpy.pi/180.0)

HIST_array = []
#read line into array
for line in inputFile_HIST.readlines():
    j = 0
    for i in line.split():
        if j == 0:
           if float(i) < r_min: r_min = float(i)
           if float(i) > r_max: r_max = float(i)
        if j == 1:
           if float(i) < theta_min: theta_min = float(i)
           if float(i) > theta_max: theta_max = float(i)
        if j == 2: 
           HIST_array.append(float(i))
           HIST_sum += float(i)
        j += 1
inputFile_HIST.close()

#print r_min*0.529177, r_max*0.529177
#print  180.0*theta_min/numpy.pi, 180.0*theta_max/numpy.pi

matrix = numpy.zeros((number_of_r_bins,number_of_theta_bins), dtype=numpy.float)

for r_count in range(number_of_r_bins):
    for theta_count in range(number_of_theta_bins):
         matrix[r_count][theta_count] = HIST_array[r_count*number_of_theta_bins + theta_count]

#if matrix_sum != 0:
#     matrix /= HIST_sum

#scale = float(sys.argv[3])

number_of_x_bins = 2048 
number_of_y_bins = number_of_x_bins/2

x_step = (2*r_max)/number_of_x_bins  
y_step = (r_max)/number_of_y_bins    # so x_step = y_step :)

cart_matrix = numpy.zeros((number_of_x_bins,number_of_y_bins), dtype=numpy.float)

for x_count in range(number_of_x_bins):

    for y_count in range(number_of_y_bins):

        x = (x_count*x_step + 0.5*x_step) - r_max
        y = y_count*y_step + 0.5*y_step

        r = (x**2 + y**2)**0.5
        theta = numpy.arctan2(y,x)

        r_bin = int(r/r_step)
        theta_bin = int(theta/theta_step)

        if r_bin < number_of_r_bins: occ = matrix[r_bin][theta_bin]
        else: occ = 0

        cart_matrix[x_count][y_count] = occ
        

matrix_sum = float(0)

matrix_max = 0
matrix_min = 1E6

matrixFile = open(sys.argv[2], 'w')

for row in matrix:
     for column in row:
          if column > matrix_max: matrix_max = column
          if column < matrix_min: matrix_min = column
          matrixFile.write(str(column) + ' ')
          matrix_sum += column
     matrixFile.write('\n')

matrixFile.close()

cart_matrixFile = open(sys.argv[3], 'w')

for row in cart_matrix:
     for column in row:
          cart_matrixFile.write(str(column) + ' ')
     cart_matrixFile.write('\n')

matrixFile.close()


if HIST_sum != matrix_sum: 
     print '\nERROR\n'
     print HIST_sum, matrix_sum

print matrix_min, matrix_max
