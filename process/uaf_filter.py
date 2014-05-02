#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt

# Define dictionary of pulse codes
# So far there are just binary phase codes, but quad-phase could work too!
p_code = {
'RECT' : [1],
'BARKER_7' : [1.,1.,1.,-1.,-1.,1.,-1.],
'BARKER_11' : [1.,1.,1.,-1.,-1.,-1.,1.,-1.,-1.,1.,-1],
'BARKER_13' : [1.,1.,1.,1.,1.,-1.,-1.,1.,1.,-1.,1.,-1.,1.],
'MLS_31' : [1,1,1,1,-1,-1,-1,1,1,-1,1,1,1,-1,1,-1,1,-1,-1,-1,-1,1,-1,-1,1,-1,1 ,1,-1,-1, 1]
}

# Root-raised cosine filter
# Define rolloff rate alpha and over-sample rate osr (number of samples per baud)

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

class filter():
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

    def matched(p_code, lp_filter, osr)

    def lowpass(lp_filter, 
