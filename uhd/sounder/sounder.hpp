void recv_clr_freq(
    uhd::usrp::multi_usrp::sptr usrp,
    uhd::rx_streamer::sptr rx_stream,
    double *center_freq,
    double bandwidth
);
int recv_to_file(
    uhd::usrp::multi_usrp::sptr usrp,
    uhd::rx_streamer::sptr rx_stream,
    const std::string &file,
    size_t samps_per_buff,
    uhd::time_spec_t start_time,
    int num_requested_samples
);
int rx_worker(
    uhd::usrp::multi_usrp::sptr usrp,
    uhd::rx_streamer::sptr rx_stream,
    const std::string &file,
    size_t samps_per_buff,
    unsigned int npulses,
    unsigned int pulse_time,
    unsigned int nave,
    uhd::time_spec_t start_time
);
void tx_worker(
    uhd::usrp::multi_usrp::sptr usrp,
    uhd::tx_streamer::sptr tx_stream,
    size_t bufflen,
    unsigned int npulses,
    unsigned int pulse_time,
    uhd::time_spec_t start_time
);

