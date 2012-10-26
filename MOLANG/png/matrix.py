#!/usr/bin/python

import sys

print 'usage: ' + sys.argv[0] + ' input.hist output.rix'

inputFile = open(sys.argv[1],'r')
outputFile = open(sys.argv[2],'w')

counter = 0

for line in inputFile.readlines():

   value = line.split()[2]
   outputFile.write(value + ' ')
   counter +=1

   if (counter%36 == 0): outputFile.write('\n') # this magic number is based on the outer loop of .hist files

inputFile.close()
outputFile.close()

