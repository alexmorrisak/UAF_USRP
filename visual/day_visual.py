#!/usr/bin/env python

# Reads a single hdf5 file.  Plots the radar's
# return power intensity as a function of range (y-axis) and 
# time (x-axis).

import numpy as np
import matplotlib.pyplot as plt
import glob
import h5py

ifile = h5py.File("soundings_20131126.h5", 'r')
nrecords = len(ifile.keys())
print "nrecords:", nrecords

#Specify parameters for things that should be plotted
mintime = 1000
maxtime = 2000
ntimes = 0

#Determine the size of the image matrix
for key,value in ifile.items():
	time = int(ifile[key].attrs['Time'])
	if time < mintime or time > maxtime:
		continue
	ntimes += 1
	nranges = int(ifile[key].attrs['N_Ranges'])
	prf = int(ifile[key].attrs['PRF (Hz)'])
print ntimes
print nranges

#Gather data for the image matrix
image_matrix = np.empty([nranges, ntimes], dtype = 'float')
i = 0
for key,value in ifile.items():
	time = int(ifile[key].attrs['Time'])
	if time < mintime or time > maxtime:
		continue
	data_arr = ifile[key]
	image_matrix[:,i] = data_arr
	i += 1

#Arrange the image matrix for plotting
image_matrix = np.flipud(image_matrix)
plt.xlabel('Time (AKST)'); plt.ylabel('Virtual Range (km)')
ext = [mintime, maxtime, 0, 1.5e5 / prf]
plt.imshow(image_matrix,extent = ext, interpolation='None'); #plt.colorbar();
##plt.imshow(data[100:600,1::], aspect = 0.1); plt.colorbar(); plt.show()
plt.show()
	
