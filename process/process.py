#!/usr/bin/env python

from optparse import OptionParser
#from pulse_gen import pulse_seq
import sys
import os
import math
import numpy as np
import matplotlib.pyplot as plt
#import scipy.signal as sig
#from array import array

#n2s = eng_notation.num_to_str

#
# Program for radar data post-processing.
# Input is baseband received data at modest sample rate.
# Outputs are decimated, filtered rx signal, and range/doppler bin grid.
# This program performs matched filtering, decimation, cross-correlation with transmitted phase code, and spectral
# estimation.
#


p_code = {
'RECT' : [1.],
'BARKER_3' : [1.,1.,-1.],
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
        rtaps = impulse_response[n-int(float(n)/5): n+int(float(n)/5)]
        rtaps = rtaps*np.hamming(np.size(rtaps))
        return rtaps


# Class that reads the input file of type complex short and converts to complex float

class rx_data():
    def __init__(self, filename, length, inx):

	inx *= 2
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
    def __init__(self, data, samp_rate, sym_rate, pc_key):
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

# Combine consecutive pulses in the received data series.
# Effectively this is a low-pass filter and it is done before the
# matched-filter convolution so it has a net improvement in efficiency.
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
	self.data_fft = self.spectrum(data, int(1.0*samp_rate / (1.0*pulse_rep_freq)), n_pulses)
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
		data_fft[:,row] = np.abs(np.fft.fftshift(np.fft.fft(data_row)))
	#plt.imshow(20*np.log10(abs(data_fft))); plt.colorbar(); plt.show()
	#for column in range(0, 4*n_pulses):
		#data_fft[column,:] = data_fft[column,:] / np.average(data_fft[column,:])
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
    parser.add_option("-r", "--samp-rate", type="int", default=250000,
                      help="set sample rate (bandwidth) [default=%default]")
    parser.add_option("-p", "--pulse-code", type="string", default="BARKER_13",
                      help="Select transmit pulse-code [default = BARKER_13]"
                      )
    parser.add_option("-f", "--prf", type="float", default=100,
                      help="Select pulse repetition frequency [default = %default Hz]"
                      )
    parser.add_option("-n", "--n-pulses", type="int", default=128,
                      help="Select number of pulses in transmitted signal [default = %default]"
                      )
    parser.add_option("-s", "--sym-rate", type="int", default=25000,
                      help="Select symbol (chip) rate of transmission [default = %default]"
                      )
    parser.add_option("-b", "--buf-time", type="int", default=80,
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
	while True:
		line=sys.stdin.readline()
		if (line):
			if line.split(":")[0] == "END":
				print "Ending process.py"
				return 0;
			if line.split(":")[0] == "RXFILE":
				rx_file = line.split(":")[1]
				rx_file = rx_file.split("\n")[0]
				sys.stdin.flush()
				print "rx_file: ", rx_file

				#rx_file = options.rx_file
    				rx = rx_data(rx_file, length, 22) # 22 corresponds to hardware delay in the USRP
    				raw_data = rx.get_rx_complex()
				print np.size(raw_data)
				SPECT_SIZE = 32
				combined_data = combine_pulses(raw_data, \
					options.n_pulses, options.samp_rate/options.prf, \
					options.n_pulses/SPECT_SIZE)
				raw_data = np.reshape(combined_data, [np.size(combined_data)], order='F')
				print "Got the data"

				# Setup the matched filter parameters and filter the rx signal
    				mfilter = matched_filter(raw_data, \
					options.samp_rate, options.sym_rate, \
					options.pulse_code) 
				print np.shape(raw_data)
				rx_filt = np.convolve(raw_data, mfilter.taps[::-1], mode='full')
				print "Done filtering"
				startinx = (np.size(rx_filt)-np.size(raw_data))/2 + 1e-6*options.samp_rate*options.buf_time + options.samp_rate/options.sym_rate*np.size(p_code[options.pulse_code])/2
				print 1e-6*options.samp_rate*options.buf_time
				print "startinx:", startinx
				rx_filt = rx_filt[startinx:startinx+np.size(raw_data)]
				if (np.size(rx_filt) < np.size(raw_data)):
					rx_filt = np.append(rx_filt, np.zeros(np.size(raw_data)-np.size(rx_filt)))
#
#				# Determine the start index of the transmit sequence and Doppler-process the signal
    				ts = doppler_process(rx_filt, options.sym_rate, options.prf, SPECT_SIZE, options.samp_rate)
				print "Done doppler-processing"
				
    				rx_filt.tofile("trunc.dat", format = 'complex64')
    				#ts.data_fft.tofile("fft.dat", format = 'complex64')
    				#ts.datarr.tofile("datarr.dat", format = 'complex64')
    				np.save("fft", ts.data_fft)
				out_file = rx_file.replace(".dat","")
				print "outfile: ", out_file
				sys_command = "rm "+rx_file
				os.system(sys_command)
    				np.save(out_file, ts.datarr)
				print "Done saving to file"
			else:
				sys.stdin.flush()
	else:
		time.sleep(0.1)	
	
    except RuntimeError, e:
	print e
	sys.exit(1)

if __name__ == "__main__":
    main()


