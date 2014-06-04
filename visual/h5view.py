#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
import h5py

f = h5py.File("../uhd/client/ionogram.20140523.1114.h5",'r')
nfreqs = len(f.keys())

dset = f[f.keys()[0]]
image = np.empty([dset.shape[0], nfreqs], float)
for i in range(0,nfreqs):
    #print f.keys()[i]
    dset = f[f.keys()[i]]
    arr = np.empty(dset.shape[0],float)
    dset.read_direct(arr)
    image[:,i] = arr
image = np.flipud(10*np.log10(image))
plt.imshow(image,aspect=0.1)
plt.colorbar()
plt.show()
    


