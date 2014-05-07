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
void tx_worker (
    uhd::usrp::multi_usrp::sptr usrp,
    uhd::tx_streamer::sptr tx_stream,
    std::complex<int16_t>* txbuff,
    size_t bufflen,
    unsigned int npulses,
    unsigned int pulse_time,
    uhd::time_spec_t start_time
){
    //create metadeta tags for transmit streams
    uhd::tx_metadata_t md;
    md.start_of_burst = true;
    md.end_of_burst = true;
    md.has_time_spec = true;
    md.time_spec = start_time;
    //std::cout << "bufflen: " << bufflen << std::endl;

    std::vector<std::complex<int16_t> > buff(bufflen, std::complex<int16_t>(0x0ffe,0x0000));
    buff[0] = std::complex<int16_t>(0x0fff,0x0001);
   
   //loop for every pulse
   for (int i=0; i<npulses; i++){
    size_t nsamples = tx_stream->send(txbuff, bufflen, md);
    //size_t nsamples = tx_stream->send(&buff.front(), bufflen, md);
    md.time_spec += 1e-3*float(pulse_time);
   }
}

