#!/usr/bin/env python

import sys, numpy

inputFile = open('output.rix', 'r')

array = []

for line in inputFile.readlines():

   array.append([])

   for element in line.split():

      array[-1].append(float(element))



inputFile.close()

array = numpy.transpose(array)/numpy.sum(array)

outputFile = open('data.dat','w')

for i in range(numpy.shape(array)[0]):

#    outputFile.write(str(i) + ' ')

    for j in range(numpy.shape(array)[1]):

        outputFile.write( str(array[i][j]) + ' ')

    outputFile.write('\n')

outputFile.close()
    

