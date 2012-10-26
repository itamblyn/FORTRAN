#! /usr/bin/env python

import numpy, matplotlib, pylab, sys

inputFile = open ('debug','r')

pitch_index = 1 
yaw_index = 3
bearing_index = 5

pitch_array = []
yaw_array = []
bearing_array = []

for line in inputFile.readlines():
    pitch_array.append(float(line.split()[pitch_index]))
    yaw_array.append(float(line.split()[yaw_index]))
    bearing_array.append(float(line.split()[bearing_index]))

inputFile.close()

print 'input read complete'

#pylab.hist(pitch_array)
#pylab.savefig('pitch.png')
#pylab.hist(yaw_array)
#pylab.savefig('yaw.png')
pylab.hist(bearing_array)
pylab.savefig('bearing.png')
