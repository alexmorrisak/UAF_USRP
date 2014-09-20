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
f = h5py.File("/home/alex/UAF_USRP/control_programs/swept_freq/"+fstring,'r')
#nfreqs = len(f.keys())/2;


dset = f['/OPower']
opower = np.empty([dset.shape[0], dset.shape[1]], float)
dset.read_direct(opower);

dset = f['/OVelocity']
ovelocity = np.empty([dset.shape[0], dset.shape[1]], float)
dset.read_direct(ovelocity);

dset = f['/XPower']
xpower = np.empty([dset.shape[0], dset.shape[1]], float)
dset.read_direct(xpower);

dset = f['/XVelocity']
xvelocity = np.empty([dset.shape[0], dset.shape[1]], float)
dset.read_direct(xvelocity);

freqs = f.attrs['Frequencies (kHz)']
minfreq = freqs[0]
maxfreq = freqs[freqs.shape[0]-1]
print "Start frequency:", minfreq
print "End frequency:", maxfreq
print "Number of frequencies:", freqs.shape[0]

ranges = f.attrs['Ranges (km)']
minrange = ranges[0]
print minrange
maxrange = ranges[ranges.shape[0]-1]
print "Start range:", minrange
print "End range:", maxrange
print "Number of range bins:", ranges.shape[0]

##ave_o = np.average(opower);
##ave_x = np.average(xpower);
##print ave_o
##print ave_x
##for i in range(0,nfreqs):
#	#opower[:,i] /= np.median(opower[:,i])
#	#xpower[:,i] /= np.median(xpower[:,i])
#	#image[:,i] /= np.average(image[:,i])
#
ext = [minfreq/1e3, maxfreq/1e3, minrange, maxrange]
asp = 2. * (maxfreq/1e3 - minfreq/1e3) / (maxrange - minrange);
#asp = 1
#print asp
#
opower = np.rot90(opower)
xpower = np.rot90(xpower)
ovelocity = np.rot90(ovelocity)
xvelocity = np.rot90(xvelocity)
##image[np.where(image < 0)] = 0
##image[np.where(image > 20)] = 20
plt.subplot(221)
plt.title('O-Power (dB)')
plt.imshow(opower,extent=ext,aspect=asp,interpolation = "none")
plt.colorbar()

plt.subplot(222)
plt.title('O-Velocity (m/s)')
plt.imshow(ovelocity,extent=ext,aspect=asp,interpolation = "none")
plt.colorbar()


plt.subplot(223)
plt.title('X-Power (dB)')
plt.imshow(xpower,extent=ext,aspect=asp,interpolation = "none")
plt.colorbar()
    
plt.subplot(224)
plt.title('X-Velocity (m/s)')
plt.imshow(xvelocity,extent=ext,aspect=asp,interpolation = "none")
plt.colorbar()

plt.show()

