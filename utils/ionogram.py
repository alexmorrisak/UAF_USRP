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
    parser.add_option("-f", "--prf", type=int, default=100,
                      help="Select pulse-repetition frequency (Hz) [default = %default]"
                      )
    (options, args) = parser.parse_args()

    return options

options = get_options()
rootdir = "/home/alex/ionosonde_data"
fdir = rootdir + "/" + options.date + "/" + options.time + "/rx.*.npy"
print fdir
flist = glob.glob(fdir)
flist = np.sort(flist)
data = np.load(flist[0])
data = data[:,0]
data = np.resize(data,[np.size(data),1])
for f in flist:
	d = np.load(f)
	d = d[:,0]
	d = np.resize(d,[np.size(d),1])
	#d = d / np.median(d)
	data = np.append(data, d, axis = 1)
data = np.flipud(data)
#temp = data[0:67,:]
#data[0:2500-67,:] = data[67:2500,:]
#data[2500-67::,:] = temp
#print np.size(data,0)
#print np.size(data,1)
for column in range(0,np.size(data,1)):
	data[:,column] = data[:,column] / np.median(data[:,column])
data = 20*np.log10(data[:,1::])
data[np.where(data > 20)] = 20
data[np.where(data < -0)] = -0
print np.shape(data)
print np.size(data,0)
print np.size(data,1)

minfleaf = flist[0].split('/')[-1]
minfreq = float(minfleaf.split('.')[2]) / 1e6
maxfleaf = flist[-1].split('/')[-1]
maxfreq = float(maxfleaf.split('.')[2]) / 1e6
ext = [minfreq, maxfreq, 0, 1.5e5 / options.prf]
#plt.imshow(image_matrix,extent = ext, interpolation='None'); #plt.colorbar();

a = (float(np.size(data,1)) / np.size(data,0)) / 4
plt.imshow(data, extent = ext, aspect=a, interpolation='None'); 
aspect = 0.3
titlestr = "Fairbanks AK, " + options.date + ", " + options.time + " AST"
plt.title(titlestr)
plt.xlabel("Frequency (MHz)")
plt.ylabel("Virtual Height (km)")
plt.colorbar();
plt.show()

