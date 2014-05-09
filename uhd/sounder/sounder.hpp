void recv_clr_freq(
    uhd::usrp::multi_usrp::sptr usrp,
    uhd::rx_streamer::sptr rx_stream,
    double *center_freq,
    double bandwidth
);
int recv_to_file(
    uhd::usrp::multi_usrp::sptr usrp,
    uhd::rx_streamer::sptr rx_stream,
    std::vector<std::complex<float> *> outdata,
    size_t samps_per_buff,
    uhd::time_spec_t start_time,
    int num_requested_samples
);
int rx_worker(
    uhd::usrp::multi_usrp::sptr usrp,
    uhd::rx_streamer::sptr rx_stream,
    std::vector<std::complex<float> *> outdata,
    size_t samps_per_buff,
    unsigned int npulses,
    unsigned int pulse_time,
    unsigned int nave,
    uhd::time_spec_t start_time
);
void tx_worker(
    uhd::usrp::multi_usrp::sptr usrp,
    uhd::tx_streamer::sptr tx_stream,
    std::complex<int16_t>* tx_buff,
    size_t bufflen,
    unsigned int npulses,
    unsigned int pulse_time,
    uhd::time_spec_t start_time
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
    std::vector<std::complex<float> *> indata,
    std::vector<std::complex<float> *> outdata,
    float *pcode,
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
    int fastdim
);
void transceive (
    uhd::usrp::multi_usrp::sptr usrp,
    uhd::tx_streamer::sptr tx_stream,
    uhd::time_spec_t start_time,
    unsigned int npulses,
    float pulse_time,
    std::complex<int16_t>* txbuff,
    size_t bufflen,
    uhd::rx_streamer::sptr rx_stream,
    std::vector<std::complex<float> *> outdata,
    size_t samps_per_pulse,
    size_t nave
);
