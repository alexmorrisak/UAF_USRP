#!/usr/bin/env python

# The USRP process saves raw binary files.  The process.py program
# processes them and yields a numpy array file.  This program takes 
# many numpy array files from a single day and packs them into a single
# hdf5 file.

import numpy as np
import matplotlib.pyplot as plt
import h5py
import glob

flist = glob.glob("/home/alex/sounder/20131126/rx_*.npy")
flist = np.sort(flist)
filename = "soundings_20131126.h5"
outfile = h5py.File(filename, 'w')
for f in flist:
	data = np.load(f)
	data = 10*np.log10(data[:,0])
	temp = f.split('_')[1]
	time = temp.split('.')[0]
	freq = "4495"
	recordname = time + '_' + freq
	outfile[recordname] = data
	record = outfile[recordname]
	record.attrs['Units'] = 'dB'
	record.attrs['Time'] = time
	record.attrs['Frequency (kHz)'] = freq
	record.attrs['P_Code'] = 'Barker_13'
	record.attrs['N_Pulses'] = 8192
	record.attrs['PRF (Hz)'] = 400
	record.attrs['N_Ranges'] = np.size(data)
outfile.close()
