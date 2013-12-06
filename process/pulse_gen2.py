#!/usr/bin/env python

from optparse import OptionParser
import sys
import math
import numpy as np
from array import array
import matplotlib.pyplot as plt

SMEAR_FACTOR = 20

p_code = {
'RECT' : [1],
'BARKER_7' : [1.,1.,1.,-1.,-1.,1.,-1.],
'BARKER_11' : [1.,1.,1.,-1.,-1.,-1.,1.,-1.,-1.,1.,-1],
'BARKER_13' : [1.,1.,1.,1.,1.,-1.,-1.,1.,1.,-1.,1.,-1.,1.],
'MLS_31' : [1,1,1,1,-1,-1,-1,1,1,-1,1,1,1,-1,1,-1,1,-1,-1,-1,-1,1,-1,-1,1,-1,1 ,1,-1,-1, 1]
}

def rrc_filter(alpha, osr):
        n = 50*osr + 1
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

class pulse_seq():
    def __init__(self, pc_key, sym_rate, pulse_rep_freq, n_pulses, samp_rate, tr_buff):

	# The pulse code to start with...
	c = np.array(p_code[pc_key])

	# Adjust parameters so that they are integer multiples of each other
	self.sym_rate = self.refine_sym_rate(sym_rate, samp_rate)
	self.prf = self.refine_prf(c, self.sym_rate, pulse_rep_freq)
	self.tr_buff = self.refine_tr_buff(self.sym_rate, tr_buff)
	self.samp_rate = samp_rate

	# Create transmit signal
	pulse = self.pad_pulse(c, self.sym_rate, self.prf)
	pulse_sequence = self.replicate(pulse, n_pulses)
	delayed_pulse_sequence = self.delay(pulse_sequence, self.tr_buff, self.sym_rate)
	oversampled_sequence = self.expand_sequence(delayed_pulse_sequence, self.sym_rate, samp_rate, 0)
	filtered_sequence = self.lp_filter(oversampled_sequence, 0.5, int(self.samp_rate / self.sym_rate))
	tx_i = self.amplify(filtered_sequence, 1.0)

	tx_i = np.multiply(tx_i, (2**15-1))
	self.tx_i = tx_i.astype('i2')
	
	self.tx_q = np.zeros(len(self.tx_i))
	self.tx_q = self.tx_q.astype('i2')

	# Create TR signal bits and slip into the LSB of the in-phase tx stream
	on = self.buffer([abs(A) for A in c] + (SMEAR_FACTOR + 1)* [1], self.sym_rate, self.tr_buff)
	onoff = self.pad_pulse(on, self.sym_rate, self.prf)
	onoff_sequence = self.replicate(onoff, n_pulses)
	onoff_sequence_full = self.append(onoff_sequence, self.sym_rate, self.tr_buff)
	self.tr = self.expand_sequence(onoff_sequence_full, self.sym_rate, samp_rate, 1)
	inx = np.nonzero(np.array(self.tr) == 1)

	tr_mask = np.empty(len(self.tx_i)).astype('i2')
	tr_mask.fill(0xfffe)
	np.bitwise_and(tr_mask, self.tx_i, self.tx_i) # Set LSB of tx_i signal to zero
	self.tr = np.zeros(len(self.tx_i), dtype = 'i2')
	self.tr[inx] = 0x0001
	np.bitwise_or(self.tr, self.tx_i, self.tx_i) # Set LSB of tx_i signal to tr value

	# Create sync signal bits to slip into the LSB of the quad-phase tx stream
	inx = np.nonzero(np.array(self.tr) != 1)[0][0]

	self.sync_mask = np.empty(len(self.tx_i)).astype('i2')
	self.sync_mask.fill(0xfffe)
	np.bitwise_and(self.sync_mask, self.tx_q, self.tx_q)
	self.sync = np.zeros(len(self.tx_i), dtype = 'i2')
	self.sync[inx] = 0x0001
	np.bitwise_or(self.sync, self.tx_q, self.tx_q)

	# Interleave 2 Tx dimensions into IQ party pack
	tx_tr = np.empty([2*np.size(self.tx_i)], dtype='i2')
	tx_tr[0::2] = self.tx_i
	tx_tr[1::2] = self.tx_q
	self.tx_tr = tx_tr

    def refine_tr_buff(self, sym_rate, tr_buff):
	new_tr_buff = math.ceil(tr_buff * sym_rate) / sym_rate
	if new_tr_buff != tr_buff:
		print "Warning: Requested T/R switch buffer is", \
			"not an integer number of samples! ",\
			"\nChanging buffer time to nearest integer number..",\
			"\nDesired time: ", 1e6 * tr_buff, "usec",\
			"\nActual time: ", 1e6 * new_tr_buff, "usec"
	return new_tr_buff

    def refine_sym_rate(self, sym_rate, samp_rate):
	new_sym_rate = samp_rate / round(samp_rate / sym_rate)
	if samp_rate % sym_rate != 0:
		print "Warning: Requested symbol rate is", \
			"not an integer number of samples! ",\
			"\nChanging symbol rate to nearest integer number..",\
			"\nDesired symbol rate: ", sym_rate, "Hz",\
			"\nActual symbol rate: ", int(new_sym_rate), "Hz"
	return new_sym_rate

    def refine_prf(self, code, sym_rate, pulse_rep_freq):
	new_prf = sym_rate / round(sym_rate / pulse_rep_freq)
	if sym_rate % pulse_rep_freq != 0:
		print "Warning: Requested pulse period is", \
			"not an integer number of pulse chips! ",\
			"\nRounding to nearest integer number..\
			\nDesired pulse period: ", 1 / float(pulse_rep_freq), "sec",\
			"\nActual pulse period: ", 1 / new_prf, "sec"
    	return new_prf

    # Append zeros to the T/R signal to compensate for the forward shift of the pulse train
    # due to the prepended buffer

    def append(self, c, sym_rate, buff):
    	return np.concatenate((np.zeros([int(buff * sym_rate)]), c))

    # Buffer the T/R signal by the desired amount of time to allow for slew time
    # in the external T/R switches

    def buffer(self, c, sym_rate, buff):
	n_ones = int(buff * sym_rate)
    	ones = [1.] * n_ones
    	return c + ones

    # For the time when the radar is receiving, fill the transmit signal
    # with zeros

    def pad_pulse(self, c, sym_rate, pulse_rep_freq):
	t_pulse = len(c) / float(sym_rate)
	t_period = 1 / float(pulse_rep_freq)
	t_rx = t_period - t_pulse
    	zeros = np.zeros([int(t_rx * sym_rate)])
    	return np.concatenate((c, zeros))

    # Delay the transmit sequence by time t_delay.  This may not be necessary, but is provided
    # to give the receiver a head start so that no samples are missed.

    def delay(self, c, t_delay, sym_rate):
	n_zeros = math.ceil(t_delay * sym_rate)
	if t_delay != n_zeros / float(sym_rate):
		print "Warning: Requested time delay is",\
			"not an integer number of pulse chips."\
			"\nRounding to nearest integer number.."\
			"\nDesired time delay: ", 1e6 * t_delay, "sec",\
			"\nActual time delay: ", 1e6 * n_zeros / sym_rate, "sec"
    	return np.concatenate((np.zeros([int(t_delay * sym_rate)]), c))
	
    # Create a sequence of pulses of arbitrary number.

    def replicate(self, tx_v, n_pulses):
	return np.tile(tx_v, n_pulses)

    # Oversample the signal to the desired sample rate
	
    def expand_sequence(self, tx_v, sym_rate, samp_rate, appval):
	samp_per_sym = int(samp_rate / sym_rate)
    	_tx_v = []
    	for i, val in enumerate(tx_v):
		_tx_v.append(val)
		if appval == 1:
    			for j in range(0, samp_per_sym-1):
    				_tx_v.append(val)
		else:
    			for j in range(0, samp_per_sym-1):
    				_tx_v.append(0)
			
    	return np.repeat(tx_v, int(samp_rate / sym_rate))

    # Choose amplification factor
	
    def amplify(self, tx_v, amplitude):
	if amplitude > 1.0 or amplitude < -1.0:
		amplitude = 1.0
		print "Amplitude must be less than or equal to 1.0."
	factor = 1 / max([abs(A) for A in tx_v])
	gain = factor * amplitude
	return np.multiply(gain, tx_v)
	#return [A * gain for A in tx_v]

    # Function used to create the filter for band-limiting the baseband signal.
    # For now the Root-Raised Cosine filter is the only option.
    # Possible additions for the future include Gaussian, etc. 

    def lp_filter(self, in_signal, alpha, osr):
	filter_taps = rrc_filter(alpha, osr)
	return np.convolve(in_signal, filter_taps[::-1], mode='valid')

def get_options():
    usage="%prog: [options]"

    parser = OptionParser(usage=usage)
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
    parser.add_option("-b", "--buff-time", type="float", default=60,
                      help="Select T/R buffer time in units of us [default = %default]"
                      )
    parser.add_option("-T", "--tx-file", type="string", default="tx.dat",
                      help="Select tx file to write [default = %default]"
                      )
    parser.add_option("-v", "--verbose", action="store_true", default=False,
                      help="Use verbose console output [default=%default]")

    (options, args) = parser.parse_args()

    return options

	
def write_to_file():
    options = get_options()
    pulse = pulse_seq(options.pulse_code, options.sym_rate, options.prf, options.n_pulses, options.samp_rate, options.buff_time / 1e6)
    out_file = open(options.tx_file, 'wb')
    out_array = array('i', pulse.tx_i)
    out_array.tofile(out_file)
    out_file.close()

if __name__ == "__main__":
    write_to_file()
