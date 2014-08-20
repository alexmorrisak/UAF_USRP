#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
import glob
from optparse import OptionParser

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
rootdir = "/home/alex/ionosonde_data"
fdir = rootdir + "/" + options.date + "/" + options.time + "/rx.*.npy"
flist = glob.glob(fdir)
flist = np.sort(flist)
power = np.load(flist[0])
power = power[:,0] / np.median(power)
power = np.resize(power,[np.size(power),1])
doppler = np.load(flist[0])
doppler = doppler[:,1]
doppler = np.resize(doppler,[np.size(doppler),1])
doppler[np.where(power < np.mean(power[400:2000]) + 1*np.std(power[400:2000]))] = 0
for f in flist:
	d = np.load(f)
	p = d[:,0]
	v = d[:,1]
	p = np.resize(p,[np.size(p),1])
	p = p/np.median(p)
	v = np.resize(v,[np.size(v),1])
	doppler = np.append(doppler, v, axis = 1)
	power = np.append(power, p, axis = 1)
doppler = np.flipud(doppler)
power = np.flipud(20*np.log10(power))
power[np.where(power < -10)] = -10
doppler = doppler - 64
#doppler[np.where(doppler < -20)] = -10
#doppler[np.where(doppler > 10)] = 10
doppler[np.where(power < np.median(power) + 3*np.std(power))] = np.nan 

minfleaf = flist[0].split('/')[-1]
minfreq = float(minfleaf.split('.')[-2]) / 1e6
maxfleaf = flist[-1].split('/')[-1]
maxfreq = float(maxfleaf.split('.')[-2]) / 1e6
print "minfreq:", minfreq
print "maxfreq:", maxfreq

a = 0.005
plt.subplot(211)
ext = [minfreq, maxfreq, 0, 1.5e5 / options.prf]
plt.imshow(power,extent = ext, aspect = a, interpolation='None'); #plt.colorbar();
titlestr = "Fairbanks AK, " + options.date + ", " + options.time + " AST"
plt.xlabel("Frequency (MHz)")
plt.ylabel("Virtual Height (km)")
plt.colorbar();

plt.subplot(212)
ext = [minfreq, maxfreq, 0, 1.5e5 / options.prf]
plt.imshow(doppler,extent = ext, aspect = a,interpolation='None'); #plt.colorbar();
aspect = 0.3
plt.xlabel("Frequency (MHz)")
plt.ylabel("Virtual Height (km)")
plt.colorbar();
plt.show()

