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
nfreqs = len(f.keys())/2;
print "nfreqs: " +  str(nfreqs)

dset = f[f.keys()[0]]
image_o = np.empty([dset.shape[0], nfreqs], float)
image_x = np.empty([dset.shape[0], nfreqs], float)
#minfreq_str=f.keys()[0];
minfreq = float((f.keys()[0]).replace("omode_","")) / 1000
print minfreq
maxfreq = float((f.keys()[nfreqs-1]).replace("omode_","")) / 1000
#maxfreq = float(f.keys()[nfreqs-1]) / 1000
print maxfreq
print "Doing omode data"
row=0
for i in range(0,2*nfreqs):
    if "xmode" in f.keys()[i]:
        continue
    print f.keys()[i]
    dset = f[f.keys()[i]]
    arr = np.empty(dset.shape[0],float)
    dset.read_direct(arr)
    image_o[:,row] = arr
    row += 1

print "Doing xmode data"
row=0
for i in range(0,2*nfreqs):
    print i
    if "omode" in f.keys()[i]:
        continue
    print f.keys()[i]
    dset = f[f.keys()[i]]
    arr = np.empty(dset.shape[0],float)
    dset.read_direct(arr)
    print f.keys()[i]
    print row
    image_x[:,row] = arr
    row += 1

#ave_o = np.average(image_o);
#ave_x = np.average(image_x);
#print ave_o
#print ave_x
#for i in range(0,nfreqs):
	#image_o[:,i] /= np.median(image_o[:,i])
	#image_x[:,i] /= np.median(image_x[:,i])
	#image[:,i] /= np.average(image[:,i])

ext = [minfreq, maxfreq, 0, 1.5e5 /options.prf]

image_o = np.flipud(10*np.log10(image_o))
image_x = np.flipud(10*np.log10(image_x))
#image[np.where(image < 0)] = 0
#image[np.where(image > 20)] = 20
plt.subplot(121)
plt.imshow(image_o,extent=ext,aspect=0.01,interpolation = "none")
plt.colorbar()
plt.subplot(122)
plt.imshow(image_x,extent=ext,aspect=0.01,interpolation = "none")
plt.colorbar()
plt.show()
    


