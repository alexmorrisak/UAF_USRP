// recv_to_file.cpp

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

/***********************************************************************
 * recv_to_file function
 **********************************************************************/
int recv_to_file(
    uhd::usrp::multi_usrp::sptr usrp,
    uhd::rx_streamer::sptr rx_stream,
    const std::string &file,
    size_t samps_per_buff,
    uhd::time_spec_t start_time,
    int num_requested_samples
){
    int num_total_samps = 0;


    uhd::rx_metadata_t md;
    std::vector<std::complex<short> > buff(1000);
    std::ofstream outfile(file.c_str(), std::ofstream::binary);
    bool overflow_message = true;
    float timeout = 0.1;

    //setup streaming
    uhd::stream_cmd_t stream_cmd = uhd::stream_cmd_t::STREAM_MODE_NUM_SAMPS_AND_DONE;
    stream_cmd.num_samps = num_requested_samples;
    stream_cmd.stream_now = false;
    stream_cmd.time_spec = start_time;
    std::cout << boost::format("Start time: %d") % start_time.get_real_secs() << std::endl;

    usrp->issue_stream_cmd(stream_cmd);
    md.error_code = uhd::rx_metadata_t::ERROR_CODE_NONE;
    uhd::time_spec_t t0 = usrp->get_time_now();
    while(num_total_samps < num_requested_samples){
	timeout = 0.2;
        size_t num_rx_samps = rx_stream->recv(&buff.front(), buff.size(), md, timeout);
        if (md.error_code == uhd::rx_metadata_t::ERROR_CODE_TIMEOUT) {
            std::cerr << boost::format("Timeout while streaming") << std::endl;
            return -1;
        }
        if (md.error_code == uhd::rx_metadata_t::ERROR_CODE_OVERFLOW){
            if (overflow_message){
                overflow_message = false;
                std::cerr << boost::format(
                    "Got an overflow indication. Please consider the following:\n"
                    "  Your write medium must sustain a rate of %fMB/s.\n"
                    "  Dropped samples will not be written to the file.\n"
                    "  Please modify this example for your purposes.\n"
                    "  This message will not appear again.\n"
                ) % (usrp->get_rx_rate()*sizeof(std::complex<short>)/1e6);
            }
            return -1;
        }
        if (md.error_code != uhd::rx_metadata_t::ERROR_CODE_NONE){
            std::cerr << "Unexpected error code 0x" << 
		md.error_code << std::endl;
	   return -1;
        }

        num_total_samps += num_rx_samps;

        outfile.write((const char*)&buff.front(), num_rx_samps*sizeof(std::complex<short>));
    }

    outfile.close();
    std::cout <<  "RXFILE:" << file.c_str() << std::endl;
    return 0;
}

