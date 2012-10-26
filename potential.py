#!/usr/local/bin/python

import numpy

for r in numpy.arange(0.0,2.0,0.025):
    v = numpy.sin(min(r - 1,numpy.pi/2))**2
    dv = 2.*numpy.sin(min(r - 1,numpy.pi/2))*numpy.cos(min(r,numpy.pi/2))
    print r,v,dv
