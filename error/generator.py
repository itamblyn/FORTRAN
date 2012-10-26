#!/usr/local/bin/python

import numpy
import Gnuplot, Gnuplot.funcutils

#numpy.random.seed(100)

dV = 0.5
dI = 0.005

R = 100

measured_data = []

for V in range(20):

    I = V / float(R)

    measured_V = V + (numpy.random.random() - 0.5)*dV
    measured_I = I + (numpy.random.random() - 0.5)*dI
    measured_data.append([measured_I, measured_V, dI, dV])

gnuplot = Gnuplot.Gnuplot()

gnuplot('set data style xyerrorbars')
gnuplot('set yrange [0:*]')
gnuplot('set xrange [0:*]')
gnuplot.plot(measured_data)
raw_input('Press any key to continue\n')

outputFile = open('data.dat','w')

outputFile.write('I V dI dV\n')

for row in measured_data:

    outputFile.write(str(row[0]) + ' ' + str(row[1]) + ' ' + str(dI) + ' ' + str(dV) + '\n')

outputFile.close()
