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
 * rx_worker function
 **********************************************************************/
int rx_worker(
    uhd::usrp::multi_usrp::sptr usrp,
    uhd::rx_streamer::sptr rx_stream,
    std::vector<std::complex<float> *> outdata,
    size_t samps_per_pulse,
    unsigned int npulses,
    unsigned int pulse_time,
    unsigned int nave,
    uhd::time_spec_t start_time
){
    int num_total_samps = 0;

    uhd::rx_metadata_t md;
    std::vector<std::complex<int16_t> > buff(samps_per_pulse);
    std::vector<std::complex<float> > fbuff(samps_per_pulse,0);
    //std::ofstream outfile(file.c_str(), std::ofstream::binary);
    bool overflow_message = true;

    //setup streaming
    uhd::stream_cmd_t stream_cmd = uhd::stream_cmd_t::STREAM_MODE_NUM_SAMPS_AND_DONE;
    stream_cmd.num_samps = samps_per_pulse;
    stream_cmd.stream_now = false;
    stream_cmd.time_spec = start_time;

    //uhd::time_spec_t t0 = usrp->get_time_now();
    for (int i=0; i<npulses/nave; i++){
        for (int j=0; j<nave; j++){
            float timeout = 0.2;
            usrp->issue_stream_cmd(stream_cmd);
            md.error_code = uhd::rx_metadata_t::ERROR_CODE_NONE;

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

            for (int k=0; k<samps_per_pulse; k++) // add recent samples to the running total
                outdata[i][k] += buff[k];
            stream_cmd.time_spec += 1e-3*float(pulse_time);
        }
        //outfile.write((const char*)&fbuff.front(), fbuff.size()*sizeof(std::complex<int16_t>));
    }
    //for (int i=0; i<samps_per_pulse; i+=10)
    //    std::cout << i << " " << outdata[0][i] << std::endl;

    //outfile.close();
    return 0;
}

