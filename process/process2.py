#!/usr/bin/env python

from optparse import OptionParser
import sys
import math
import numpy as np
import matplotlib.pyplot as plt
import h5py

#
# Program for radar data post-processing.
# Input is baseband received data at modest sample rate.
# Outputs are decimated, filtered rx signal, and range/doppler bin grid.
# This program performs matched filtering, decimation, cross-correlation with transmitted phase code, and spectral
# estimation.
#


p_code = {
'RECT' : [1.],
'BARKER_7' : [1.,1.,1.,-1.,-1.,1.,-1.],
'BARKER_11' : [1.,1.,1.,-1.,-1.,-1.,1.,-1.,-1.,1.,-1],
'BARKER_13' : [1.,1.,1.,1.,1.,-1.,-1.,1.,1.,-1.,1.,-1.,1.],
'MLS_31' : [1,1,1,1,-1,-1,-1,1,1,-1,1,1,1,-1,1,-1,1,-1,-1,-1,-1,1,-1,-1,1,-1,1 ,1,-1,-1, 1]
}

def rrc_filter(alpha, osr, ntaps):
        n = ntaps
        freq_response = np.empty([n], dtype = 'float')
        dBW = float(n) / osr # Digital bandwidth

        freq_response[0 : (1-alpha)*dBW] = 1/(2*dBW)
        freq_response[(1-alpha)*dBW : (1+alpha)*dBW] = \
                1/(4*dBW) * (1+np.cos(np.linspace(0, np.pi, \
                np.size(freq_response[(1-alpha)*dBW : (1+alpha)*dBW]))))
        freq_response[(1+alpha)*dBW : n] = 0

        freq_response = np.sqrt(freq_response)
        freq_response = np.concatenate((freq_response, freq_response[::-1]))

        impulse_response = np.real(np.fft.fftshift(np.fft.ifft(freq_response)))
        return impulse_response[n-int(float(n)/5): n+int(float(n)/5)]


# Class that reads the input file of type complex short and converts to complex float

class rx_data():
    def __init__(self, filename, length, inx):

        in_file = open(filename, 'rb')
        in_arr = np.fromfile(in_file, dtype='i2')
        in_file.close()

	self.i = in_arr[inx:length+inx:2]
	self.q = in_arr[inx+1:length+inx:2]

    def get_rx_complex(self):
	return (self.i + 1j * self.q).astype('complex64') / 2**15

    def get_rx_i(self):
	return self.i.astype('float32') / 2**15

    def get_rx_q(self):
	return self.q.astype('float32') / 2**15

class matched_filter():
    def __init__(self, samp_rate, sym_rate, pc_key):
        lp_filter = self.set_matched_filter(	 1.0,
	               				 samp_rate,
       	        				 sym_rate,
       	        				 0.2,
       	        				 5 * int(samp_rate / sym_rate),
						 pc_key)
	self.taps = lp_filter

    def set_matched_filter(self, gain, samp_rate, sym_rate, alpha, ntaps, pc_key):
	rrc = rrc_filter(alpha, samp_rate/sym_rate, 513)
	pc = p_code[pc_key]
	pcode = []
	for i,val in enumerate(pc):
		pcode.append(val)
		for j in range(0,samp_rate/sym_rate-1):
			pcode.append(0)
	pcode = np.array(pcode)

	filter_taps = np.convolve(rrc, pcode)
	return filter_taps

def combine_pulses(idata, npulses, nranges, decim_rate):
	odata = np.empty([nranges, npulses / decim_rate], dtype = 'complex')
	for i in range(0,npulses/decim_rate):
		datacombine = np.reshape(idata[i*nranges*decim_rate: (i+1)*nranges*decim_rate], [nranges, decim_rate], order='F')
		ndatacombine = np.sum(datacombine, axis=1)
		odata[:,i] = ndatacombine
	return odata

class doppler_process():
    def __init__(self, in_data, sym_rate, pulse_rep_freq, n_pulses, samp_rate):
	data = in_data 
	self.data_fft = self.spectrum(data, int(samp_rate / pulse_rep_freq), n_pulses)
	self.datarr = self.find_reflections(self.data_fft)

    def spectrum(self, data, n_ranges, n_pulses):
	data_mat = np.reshape(data, (n_pulses, n_ranges), order='C')
	nspectrum = 4*n_pulses
	data_fft = np.empty([nspectrum, n_ranges], dtype = 'float') 
	window = np.hamming(n_pulses)
	window[0::2] = -1 * window[0::2]
	for row in range(0, n_ranges):
		data_row = np.empty(nspectrum, dtype='complex')
		data_row[0:n_pulses] = np.multiply(window, data_mat[:,row])
		data_row[n_pulses::] = np.zeros([nspectrum - n_pulses])
		data_fft[:,row] = np.abs((np.fft.fft(data_row)))
	#plt.imshow(20*np.log10(abs(data_fft))); plt.colorbar(); plt.show()
	return data_fft

    def find_reflections(self, data):
	nbins = np.size(data, 0)
	nranges = np.size(data,1)
	peak_power = np.empty([nranges,2], dtype = 'float')
	for row in range(0, nranges):
		peak_power[row,0] = max(data[:, row])
		peak_power[row,1] = np.argmax(data[:, row])
	return peak_power
 	

def get_options():
    usage="%prog: [options]"

    parser = OptionParser(usage=usage)
    parser.add_option("-t", "--time", type="int", default=1200,
                      help="Declare time of sounding (AKST) [default=%default]")
    parser.add_option("-c", "--freq", type="int", default=4500,
                      help="Declare RF center frequency (kHz) [default=%default]")
    parser.add_option("-r", "--samp-rate", type="int", default=250000,
                      help="set sample rate (bandwidth) [default=%default]")
    parser.add_option("-p", "--pulse-code", type="string", default="BARKER_13",
                      help="Select transmit pulse-code [default = BARKER_13]"
                      )
    parser.add_option("-f", "--prf", type="int", default=200,
                      help="Select pulse repetition frequency [default = %default Hz]"
                      )
    parser.add_option("-n", "--n-pulses", type="int", default=512,
                      help="Select number of pulses in transmitted signal [default = %default]"
                      )
    parser.add_option("-s", "--sym-rate", type="int", default=50000,
                      help="Select symbol (chip) rate of transmission [default = %default]"
                      )
    parser.add_option("-b", "--buf-time", type="int", default=60,
                      help="Select T/R buffer time in units of us [default = %default]"
                      )
    parser.add_option("-R", "--rx-file", type="string", default="rx.dat",
                      help="Select rx file to process [default = %default]"
                      )
    parser.add_option("-O", "--out-file", type="string", default="datarr",
                      help="Select file to output [default = %default]"
                      )
    parser.add_option("-v", "--verbose", action="store_true", default=False,
                      help="Use verbose console output [default=%default]")

    (options, args) = parser.parse_args()

    return options

def main():
	#EXAMPLE: process.py --rx-file rx.dat --n-pulses 
    try:
	# Grab the program options
	options = get_options()
	print "Got program options"

	# Read the input file and convert to complex float format
	length = 2*options.n_pulses*options.samp_rate/options.prf
	print length
    	rx = rx_data(options.rx_file, length, 0)
    	raw_data = rx.get_rx_complex()
	print np.size(raw_data)
	SPECT_SIZE = 64
	combined_data = combine_pulses(raw_data, options.n_pulses, options.samp_rate/options.prf, options.n_pulses/SPECT_SIZE)
	raw_data = np.reshape(combined_data, [np.size(combined_data)], order='F')
	print "Got the data"

	# Setup the matched filter parameters and filter the rx signal
    	mfilter = matched_filter(options.samp_rate, options.sym_rate, options.pulse_code) 
	print np.shape(raw_data)
	rx_filt = np.convolve(raw_data, mfilter.taps[::-1], mode='same')
	print "Done filtering"
#
#	# Determine the start index of the transmit sequence and Doppler-process the signal
    	ts = doppler_process(rx_filt, options.sym_rate, options.prf, SPECT_SIZE, options.samp_rate)
	print "Done doppler-processing"
	
    	#rx_filt.tofile("trunc.dat", format = 'complex64')
    	#np.save("fft", ts.data_fft)
	nranges = ts.datarr[:,0]
	if not os.path.isfile(out_file):
		f = h5py.File(options.out_file, 'w')
		rawpower = f.create_dataset('raw_power', (1,nranges), maxshape=(None,nranges))
		rawdoppler = f.create_dataset('raw_doppler', (1,nranges), maxshape=(None,nranges))
		ranges = f.create_dataset('Ranges', (nranges))
		times = f.create_dataset('Time', (1), maxshape=(None))

		f[rawpower] = ts.datarr[:,0]
		f[rawdoppler] = ts.datarr[:,1]
		f[ranges] = np.linspace(0,3e5/(2*options.prf),nranges)
		f[times] = options.time

		rawpower.attrs['Units'] = 'dB'
        	rawpower.attrs['Frequency (kHz)'] = options.freq
        	rawpower.attrs['P_Code'] = options.pulse_code
        	rawpower.attrs['N_Pulses'] = options.n_pulses
        	rawpower.attrs['PRF (Hz)'] = options.prf
        	rawpower.attrs['N_Ranges'] = np.size(data)
		rawdoppler.attrs['Units'] = 'Hz'
        	rawdoppler.attrs['Frequency (kHz)'] = options.freq
        	rawdoppler.attrs['P_Code'] = options.pulse_code
        	rawdoppler.attrs['N_Pulses'] = options.n_pulses
        	rawdoppler.attrs['PRF (Hz)'] = options.prf
        	rawdoppler.attrs['N_Ranges'] = np.size(data)

		f.close()
	else:
		f = h5py.File
		f[rawpower] = np.concatenate((f[rawpower],ts.datarr[:,0]))
		f[rawdoppler] = np.concatenate((f[rawdoppler],ts.datarr[:,1]))
		f[times] = np.concatenate((f[times],options.time))
	
    	np.save(options.out_file, ts.datarr)
	print "Done saving to file"
	
    except RuntimeError, e:
	print e
	sys.exit(1)

if __name__ == "__main__":
    main()


