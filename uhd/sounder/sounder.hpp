void capture_spectrum(
    uhd::usrp::multi_usrp::sptr usrp,
    uhd::rx_streamer::sptr rx_stream,
    size_t start_freq_khz,
    size_t stop_freq_khz,
    size_t bandwidth_khz,
    float* periodogram
);
int lp_filter(
    std::vector<std::complex<float> *> indata,
    std::vector<std::complex<float> *> outdata,
    int slowdim,
    int fastdim,
    float samprate,
    float bw,
    int decimrate
);
int matched_filter(
    std::complex<float>*** indata,
    std::vector<std::complex<float> *> outdata,
    float** pcode,
    int pcode_length,
    int slowdim,
    int fastdim,
    int osr
);

int doppler_process(
    std::vector<std::complex<float> *> indata,
    float* outpow,
    float* outvel,
    int slowdim,
    int fastdim,
    int osr,
    int first_range_inx
);
void transceive(
    uhd::usrp::multi_usrp::sptr usrp,
    uhd::tx_streamer::sptr tx_stream,
    uhd::rx_streamer::sptr rx_stream,
    unsigned int npulses,
    float pulse_time,
    std::vector<std::complex<int16_t> >* txbuff0,
    std::vector<std::complex<int16_t> >* txbuff1,
    float tx_ontime,
    std::complex<int16_t>** outdata,
    size_t samps_per_pulse
);
