#! /usr/bin/env python

import numpy, matplotlib, pylab, sys

if len(sys.argv) == 1:
     print 'usage: ' + sys.argv[0] + ' circle.rix [vmin, vmax] circle.png'

inputFile = open (sys.argv[1],'r')

line_counter = 0

PROBABILITY_array = []
for line in inputFile.readlines():
    PROBABILITY_array.append([])

    row_counter = 0
    for i in line.split():
        PROBABILITY_array[-1].append(float(i))
        row_counter += 1

    line_counter += 1

inputFile.close()

print line_counter, row_counter

PROBABILITY_array = numpy.reshape(PROBABILITY_array,(line_counter,row_counter))


if len(sys.argv) == 5:
     vmin_value = float(sys.argv[2])
     vmax_value = float(sys.argv[3])
     savefile = sys.argv[4]

else:
     vmin_value = PROBABILITY_array.min()
     vmax_value = PROBABILITY_array.max()
     savefile = sys.argv[2]

#pylab.figure(num=1,figsize=(3,3),facecolor='w',edgecolor='k')

im = pylab.imshow(PROBABILITY_array,vmin=vmin_value,vmax=vmax_value)
im.set_interpolation('nearest')
pylab.xticks([1],' ')
pylab.yticks([1],' ')
pylab.colorbar(cax=pylab.axes([0.85,0.1,0.05,0.8]))
pylab.savefig(savefile)
#pylab.show()

#pylab.pcolor(PROBABILITY_array,vmin=vmin_value,vmax=vmax_value)
#pylab.xticks([0,10,20,30,40,50,60],["0","60","120","180","240","300","360"])
#pylab.yticks([0,5,10,15,20,25,30],["-90","-60","-30","0","30","60","90"])
#pylab.savefig(savefile)
#pylab.show()
