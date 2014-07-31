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
f = h5py.File("/home/radar/UAF_USRP/uhd/client/"+fstring,'r')
nfreqs = len(f.keys());

dset = f[f.keys()[0]]
image = np.empty([dset.shape[0], nfreqs], float)
#minfreq_str=f.keys()[0];
minfreq = float((f.keys()[0]).replace("omode_","")) / 1000
print minfreq
maxfreq = float((f.keys()[nfreqs-1]).replace("omode_","")) / 1000
#maxfreq = float(f.keys()[nfreqs-1]) / 1000
print maxfreq
for i in range(0,nfreqs):
    dset = f[f.keys()[i]]
    arr = np.empty(dset.shape[0],float)
    dset.read_direct(arr)
    image[:,i] = arr
ave = np.average(image);
print ave
for i in range(0,nfreqs):
	image[:,i] /= np.median(image[:,i])
	#image[:,i] /= np.average(image[:,i])

ext = [minfreq, maxfreq, -100, 1.5e5 /options.prf]

image = np.flipud(10*np.log10(image))
#image[np.where(image < 0)] = 0
#image[np.where(image > 20)] = 20
plt.imshow(image,extent=ext,aspect=0.01,interpolation = "none")
plt.colorbar()
plt.show()
    


