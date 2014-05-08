// tx_worker.cpp
// 
#include <uhd/utils/thread_priority.hpp>
#include <uhd/utils/safe_main.hpp>
#include <uhd/utils/static.hpp>
#include <uhd/usrp/multi_usrp.hpp>
#include <uhd/exception.hpp>
#include <boost/thread/thread.hpp>
#include <boost/program_options.hpp>
#include <boost/math/special_functions/round.hpp>
#include <boost/foreach.hpp>
#include <boost/format.hpp>
#include <boost/algorithm/string/split.hpp>
#include <boost/algorithm/string/classification.hpp>
#include <iostream>
#include <string>
#include <fstream>
#include <complex>
#include <csignal>
#include <cmath>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <fftw3.h>


typedef std::complex<int16_t>  sc16;

/***********************************************************************
 * tx_worker function
 * A function to be used as a boost::thread_group thread for transmitting
 **********************************************************************/
void transceive (
    uhd::usrp::multi_usrp::sptr usrp,
    uhd::tx_streamer::sptr tx_stream,
    uhd::time_spec_t start_time,
    unsigned int npulses,
    unsigned int pulse_time,
    std::complex<int16_t>* txbuff,
    size_t bufflen,
    uhd::rx_streamer::sptr rx_stream,
    std::vector<std::complex<float> *> outdata,
    size_t samps_per_pulse,
    size_t nave
){
    //create metadeta tags for transmit streams
    uhd::tx_metadata_t md;
    md.start_of_burst = true;
    md.end_of_burst = true;
    md.has_time_spec = true;
    md.time_spec = start_time;
    float txon = (float) bufflen / usrp->get_tx_rate();

    //create metadata tags for receive stream
    uhd::rx_metadata_t rxmd;
    std::vector<std::complex<int16_t> > buff(samps_per_pulse,0);
    std::vector<std::complex<float> > fbuff(samps_per_pulse,0);

    uhd::stream_cmd_t stream_cmd = uhd::stream_cmd_t::STREAM_MODE_NUM_SAMPS_AND_DONE;
    stream_cmd.num_samps = samps_per_pulse;
    stream_cmd.stream_now = false;
    stream_cmd.time_spec = start_time;

    //std::vector<std::complex<int16_t> > buff(bufflen, std::complex<int16_t>(0x0ffe,0x0000));
    //buff[0] = std::complex<int16_t>(0x0fff,0x0001);
   
   usrp->set_gpio_attr("TXA","CTRL",0x0, 0x4);
   usrp->set_gpio_attr("TXA","DDR",0x4, 0x4);
   //loop for every pulse
   for (int i=0; i<npulses/nave; i++){
    for (int j=0; j<nave; j++){
        float timeout = 0.1;
        usrp->set_command_time(start_time-50e-6,0);
        usrp->set_gpio_attr("TXA","OUT",0x4, 0x4);

        size_t nsamples = tx_stream->send(txbuff, bufflen, md);

        usrp->issue_stream_cmd(stream_cmd);

        usrp->set_command_time(start_time+txon+50e-6,0);
        usrp->set_gpio_attr("TXA","OUT",0x0, 0x4);

        size_t nrx_samples = rx_stream->recv(&buff.front(), buff.size(), rxmd, timeout);
        for (int k=0; k<samps_per_pulse; k++)
            outdata[i][k] += buff[k];


        start_time += 1e-3*float(pulse_time);
        md.time_spec = start_time;
        stream_cmd.time_spec = start_time;
    }
   }
}

