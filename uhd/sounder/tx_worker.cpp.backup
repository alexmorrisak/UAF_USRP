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
 * transmit_worker function
 * A function to be used as a boost::thread_group thread for transmitting
 **********************************************************************/
void transmit_worker(
    uhd::usrp::multi_usrp::sptr usrp,
    uhd::tx_streamer::sptr tx_stream,
    const std::string &file,
    size_t *ntx_samps,
    uhd::time_spec_t start_time
){
    //create metadeta tags for transmit streams
    uhd::tx_metadata_t md;
    md.start_of_burst = true;
    md.end_of_burst = false;
    md.has_time_spec = true;
    md.time_spec = start_time;

    //open binary file for baseband tx stream
    size_t spb = tx_stream->get_max_num_samps();
    std::vector<sc16> buff(spb);
    std::ifstream infile(file.c_str(), std::ifstream::binary);

    //*ntx_samps = infile.gcount()/sizeof(samp_type);
    //std::cout << *ntx_samps << std::endl;
    //loop until the entire file has been read
    while(not md.end_of_burst){
        infile.read((char*)&buff.front(), buff.size()*sizeof(sc16));
        md.end_of_burst = infile.eof();
        tx_stream->send(&buff.front(), buff.size(), md);
	md.has_time_spec = false;
	md.start_of_burst = false;
    }
    md.end_of_burst=true;
    tx_stream->send("",0,md);

    infile.close();
}

