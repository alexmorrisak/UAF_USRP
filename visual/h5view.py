#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
import h5py
from optparse import OptionParser
import glob

def get_options():
    usage="%prog: [options]"

    parser = OptionParser(usage=usage)
    parser.add_option("-d", "--date", type="string", default="20140205",
                      help="Select date (YYYYMMDD) [default=%default]")
    parser.add_option("-t", "--time", type="string", default="1200",
                      help="Select time (hhmm) [default = %default]"
                      )
    parser.add_option("-f", "--prf", type=int, default=200,
                      help="Select pulse-repetition frequency (Hz) [default = %default]"
                      )
    (options, args) = parser.parse_args()

    return options

options = get_options()
fstring = "ionogram." + options.date + "." + options.time + ".h5"
#print fstring
f = h5py.File("/home/radar/UAF_USRP/uhd/client/"+fstring,'r')
#nfreqs = len(f.keys())/2;


dset = f['/Omode']
print dset.shape[0]
print dset.shape[1]
image_o = np.empty([dset.shape[0], dset.shape[1]], float)
dset.read_direct(image_o);

dset = f['/Xmode']
print dset.shape[0]
print dset.shape[1]
image_x = np.empty([dset.shape[0], dset.shape[1]], float)
dset.read_direct(image_x);
print "Number of range bins:", dset.shape[1]
print "Number of frequencies:", dset.shape[0]

freqs = f.attrs['Frequencies(kHz)']
minfreq = freqs[0]
print minfreq
maxfreq = freqs[freqs.shape[0]-1]
print "Number of frequencies:", freqs.shape[0]

ranges = f.attrs['Ranges(km)']
minrange = ranges[0]
print minrange
maxrange = ranges[ranges.shape[0]-1]
print "Number of range bins:", ranges.shape[0]

##ave_o = np.average(image_o);
##ave_x = np.average(image_x);
##print ave_o
##print ave_x
##for i in range(0,nfreqs):
#	#image_o[:,i] /= np.median(image_o[:,i])
#	#image_x[:,i] /= np.median(image_x[:,i])
#	#image[:,i] /= np.average(image[:,i])
#
ext = [minfreq/1e3, maxfreq/1e3, minrange, maxrange]
asp = 2. * (maxfreq/1e3 - minfreq/1e3) / (maxrange - minrange);
#asp = 1
print asp
#
image_o = np.rot90(image_o)
image_x = np.rot90(image_x)
##image[np.where(image < 0)] = 0
##image[np.where(image > 20)] = 20
plt.subplot(121)
plt.title('O-Mode (dB)')
plt.imshow(image_o,extent=ext,aspect=asp,interpolation = "none")
#plt.imshow(image_o,aspect=asp,interpolation = "none")
plt.colorbar()
plt.subplot(122)
plt.title('X-Mode (dB)')
plt.imshow(image_x,extent=ext,aspect=asp,interpolation = "none")
#plt.imshow(image_x,aspect=asp,interpolation = "none")
plt.colorbar()
plt.show()
    


