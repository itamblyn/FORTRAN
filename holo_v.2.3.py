#!/usr/local/bin/python

import numpy
import random
import array
import time

start_time = time.time()

input_file_source_name = raw_input('Input filename: ')

# coordinates from LEEPS/OPS

xmin = -0.5
xmax = 0.5
ymin = -0.5
ymax = 0.5
zmin = 0.0
zmax = 0.0

xmin *= 100
xmax *= 100
ymin *= 100
ymax *= 100
zmin *= 100
zmax *= 100

# decides how many times to "enhance" the images
power = 1 

threshold_intensity = input('Set threshold intensity [0-255]: ')  ## all intensities less this value will be set to zero
threshold_intensity = numpy.power(threshold_intensity, power) ## to enhance the contrast

threshold_neighbours = input('Set erosion parameter [0-5]: ')  ## all points with less than this number of neighbours will be discarded

printNumbers = False 
simulate_file = False

number_of_source_files = input('Number of reconstructed slices? ')

min_cluster_size = input('Set minimum cluster size to print: ')

output_file = 'coordinates'

## import file

if (simulate_file == False):

     nx = 1024
     ny = 1024
     nz  = number_of_source_files

     intensity_array = numpy.zeros((nx,ny,nz), dtype=numpy.Float)
     
#     input_file_source_name = input('Input filename, without extension: ')
     
     i = 0
     j = 0
     k = 0
     while k < number_of_source_files:
          input_file = input_file_source_name + ' (' + repr(k) + ').pgm'
          fileobj = open(input_file, mode='rb')
          fileobj.readline()
          fileobj.readline()
          fileobj.readline()
          binvalues = array.array('B')
          binvalues.read(fileobj, nx * ny)

          import_array = numpy.array(binvalues, dtype=numpy.Float)
          import_array = numpy.reshape(import_array, (nx, ny))
          
          print 'input file ' + input_file_source_name + ' (' + repr(k) + ').pgm sucessfully read'
          fileobj.close()
          
          while i < nx:
               while j < ny:
                    intensity_array[i][j][k] = numpy.power(import_array[i][j], power)
                    j += 1
               j = 0
               i +=1
          i = 0
          j = 0
          k += 1
     
     import_array = 0

## end import file

## simulate file ## 

if (simulate_file == True):

     nx = 10
     ny = 10
     nz = 1

     intensity_array = numpy.zeros((nx,ny,nz))

     i = 0
     j = 0
     k = 0

     while i < nx:
          while j < ny:
               while k < nz:
                    intensity_array[i][j][k] += random.randint(0,255)
                    k+=1
               k = 0
               j += 1
          k= 0
          j= 0
          i+= 1


## end simulate file ##

## get rid of all voxels below threshold value ##

i = 0
j = 0
k = 0

while i < nx:
     while j < ny:
          while k < nz:
               if intensity_array[i][j][k] < threshold_intensity: intensity_array[i][j][k] = 0
               k+=1
          k = 0
          j += 1
     k = 0
     j = 0
     i+= 1

print 'voxels below ' + repr(threshold_intensity) + ' sucessfully removed'


## end get rid of below threshold voxels ##


## erode (gets rid of single voxels and voxels that do not have enough nearest neighbours) ##

eroded_array = numpy.zeros((nx,ny,nz), dtype=numpy.Float)

neighbour_count = 0

i = 0
j = 0
k = 0

while i < nx:
     while j < ny:
          while k < nz:
          
               ## note, this is based on 6 neighbours...another loop would have to be used for more ##
               if (k+1 < nz):
                    if (intensity_array[i][j][k+1] != 0): neighbour_count +=1 ## north
               if (k-1 >= 0):
                    if (intensity_array[i][j][k-1] != 0): neighbour_count +=1 ## south
               if (i+1 < nx):
                    if (intensity_array[i+1][j][k] !=0): neighbour_count +=1 ## east 
               if (i-1 >= 0):
                    if(intensity_array[i-1][j][k] !=0): neighbour_count +=1 ## west 
               if (j+1 < ny):
                    if(intensity_array[i][j+1][k] !=0): neighbour_count +=1 ## up 
               if (j-1 >= 0):
                    if(intensity_array[i][j-1][k] !=0): neighbour_count +=1 ## down
                    
               if (neighbour_count >= threshold_neighbours): eroded_array[i][j][k] = intensity_array[i][j][k]
               neighbour_count = 0
               k+=1
          k = 0
          j += 1
     k = 0
     j = 0
     i+= 1

print 'cluster erode complete'

## end erode ##

## delete intensity array ##

intensity_array = 0

## end delete of intensity array ##


## start labelling clusters

print 'labelling clusters . . .'

label_array = numpy.ones((nx,ny,nz))
label_array *= nx*ny*nz

i = 0
j = 0
k = 0

cluster_indice = 0 

# based on the number of nearest neighbours / 2

cluster_coordinates_array = []

while i < nx:
     while j < ny:
          while k < nz:
          
               ## note, this is based on 6 neighbours . . . another loop would have to be used for more ##

               if (eroded_array[i][j][k] != 0):
                    
                    neighbour_label_array = numpy.ones(3)
                    neighbour_label_array *= nx*ny*nz
               
                    if (i-1 >= 0):
                         neighbour_label_array[0] = label_array[i - 1][j][k]
                    if (j-1 >= 0):
                         neighbour_label_array[1] = label_array[i][j - 1][k]
                    if (k-1 >= 0):
                         neighbour_label_array[2] = label_array[i][j][k - 1]

                         
                    if (min(neighbour_label_array) == nx*ny*nz):

                         label_array[i][j][k] = cluster_indice
                         cluster_coordinates_array.append([])
                         cluster_coordinates_array[len(cluster_coordinates_array) -1].append([i,j,k])
                         cluster_indice +=1                         
                         
                    else: 
                         label_array[i][j][k] = min(neighbour_label_array)
                         cluster_coordinates_array[min(neighbour_label_array)].append([i,j,k])                    
                          
                         if (neighbour_label_array[0] != min(neighbour_label_array) and neighbour_label_array[0] != nx*ny*nz):
                              a = 0
                              while a < len(cluster_coordinates_array[neighbour_label_array[0]]):
#                                   print ' a:' + repr(a) + ' len(cluster_coordinates_array[neighbour_label_array[0]]): ' + repr(len(cluster_coordinates_array[neighbour_label_array[0]]))
#                                  print cluster_coordinates_array[neighbour_label_array[0]][a][0]
                                   label_array[cluster_coordinates_array[neighbour_label_array[0]][a][0]][cluster_coordinates_array[neighbour_label_array[0]][a][1]][cluster_coordinates_array[neighbour_label_array[0]][a][2]] = min(neighbour_label_array)
                                   a += 1
                                   
                              a = 0
                              while a < len(cluster_coordinates_array[neighbour_label_array[0]]):
                                   cluster_coordinates_array[min(neighbour_label_array)].append(cluster_coordinates_array[neighbour_label_array[0]][a])
                                   a += 1
                                                                                 
                         if (neighbour_label_array[1] != min(neighbour_label_array) and neighbour_label_array[1] != nx*ny*nz):
                              a = 0
                              while a < len(cluster_coordinates_array[neighbour_label_array[1]]):
#                                   print ' a:' + repr(a) + ' len(cluster_coordinates_array[neighbour_label_array[1]]): ' + repr(len(cluster_coordinates_array[neighbour_label_array[1]]))
#                                   print cluster_coordinates_array[neighbour_label_array[1]][a][0]
                                   label_array[cluster_coordinates_array[neighbour_label_array[1]][a][0]][cluster_coordinates_array[neighbour_label_array[1]][a][1]][cluster_coordinates_array[neighbour_label_array[1]][a][2]] = min(neighbour_label_array)
                                   a += 1
                              
                              a = 0     
                              while a < len(cluster_coordinates_array[neighbour_label_array[1]]):
                                   cluster_coordinates_array[min(neighbour_label_array)].append(cluster_coordinates_array[neighbour_label_array[1]][a])
                                   a += 1      


                         if (neighbour_label_array[2] != min(neighbour_label_array) and neighbour_label_array[2] != nx*ny*nz):
                              a = 0
                              while a < len(cluster_coordinates_array[neighbour_label_array[2]]):
#                                   print ' a:' + repr(a) + ' len(cluster_coordinates_array[neighbour_label_array[2]]): ' + repr(len(cluster_coordinates_array[neighbour_label_array[2]]))
#                                   print cluster_coordinates_array[neighbour_label_array[2]][a][0]
                                   label_array[cluster_coordinates_array[neighbour_label_array[2]][a][0]][cluster_coordinates_array[neighbour_label_array[2]][a][1]][cluster_coordinates_array[neighbour_label_array[2]][a][2]] = min(neighbour_label_array)
                                   a += 1
                                   
                              a = 0     
                              while a < len(cluster_coordinates_array[neighbour_label_array[2]]):
                                   cluster_coordinates_array[min(neighbour_label_array)].append(cluster_coordinates_array[neighbour_label_array[2]][a])
                                   a += 1
                                   
                           
               k+=1
          k = 0
          j += 1
     k = 0
     j = 0
     i+= 1


## end labelling clusters


## print results

if (printNumbers == True):

     i = 0
     j = 0
     k = 0

     line = ' '

     while i < nx:
          while j < ny:
               while k < nz:
                    eroded_value = eroded_array[i][j][k]
                    line += ' %(#)03d' % {"#": eroded_value}
                    k+=1
               k = 0
               j += 1
          print line
          line = ' '
          k = 0
          j = 0
          i+= 1


     print ''

     i = 0
     j = 0
     k = 0

     line = ' '

     while i < nx:
          while j < ny:
               while k < nz:
                    label_value = label_array[i][j][k]
                    if (label_value != nx*ny*nz): line += ' %(#)01d' % {"#": label_value}
                    else: line += '  '
                    k+=1
               k = 0
               j += 1
          print line
          line = ' '
          k = 0
          j = 0
          i+= 1
     
## end prints

## centre of mass calculation

xmass_running_array = numpy.zeros(cluster_indice, dtype=numpy.Float)
ymass_running_array = numpy.zeros(cluster_indice, dtype=numpy.Float)
zmass_running_array = numpy.zeros(cluster_indice, dtype=numpy.Float)
mass_running_array  = numpy.zeros(cluster_indice, dtype=numpy.Float)

xlength = xmax - xmin
ylength = ymax - ymin
zlength = zmax - zmin

xstep = float(xlength)/float(nx)
ystep = float(ylength)/float(ny)
zstep = float(zlength)/float(nz)

for i in range(nx):
     for j in range(ny):
          for k in range(nz):
               if (label_array[i][j][k] != nx*ny*nz):
                    xmass_running_array[label_array[i][j][k]] += (i + 0.5)*xstep*eroded_array[i][j][k]
                    ymass_running_array[label_array[i][j][k]] += (j + 0.5)*ystep*eroded_array[i][j][k]
                    zmass_running_array[label_array[i][j][k]] += (k + 0.5)*zstep*eroded_array[i][j][k]
                    mass_running_array[label_array[i][j][k]]  += eroded_array[i][j][k]                                          

cluster_counter = 0

for i in range(cluster_indice):
     if (mass_running_array[i] != 0 and len(cluster_coordinates_array[i]) >= min_cluster_size): cluster_counter+=1

output_file += '.xyz'

fileobj = open(output_file, mode='w')
fileobj.write(str(cluster_counter) + '\n')
fileobj.write('cluster_size, xcoor, ycoor, zcoor\n')

for i in range(cluster_indice):
     if (mass_running_array[i] != 0 and len(cluster_coordinates_array[i]) >= min_cluster_size):
          fileobj.write(str(len(cluster_coordinates_array[i])) + '   ' + repr(float(xmass_running_array[i])/float(mass_running_array[i])) + '   ' + repr(float(ymass_running_array[i])/float(mass_running_array[i])) + '   ' + repr(float(zmass_running_array[i])/float(mass_running_array[i])) + '\n')
fileobj.close()

print 'number_of_source_files: ' + repr(number_of_source_files)
print 'threshold_intensity: ' + repr(threshold_intensity)
print 'threshold_neighbours: ' + repr(threshold_neighbours)
print 'min_cluster_size: ' + repr(min_cluster_size)

end_time = time.time()
print 'Time taken: ' + repr(end_time - start_time) + ' s'
