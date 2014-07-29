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

#include <thread>


typedef std::complex<int16_t>  sc16;
extern int verbose;

void tx_worker(
    unsigned int bufflen,
    uhd::tx_streamer::sptr tx_stream,
    uhd::time_spec_t start_time,
    std::complex<int16_t>* vec_ptr,
    int end
);

void rx_worker(
    uhd::rx_streamer::sptr rx_stream,
    unsigned int samps_per_pulse,
    std::vector<std::complex<int16_t>* >& recv_ptr
);

void transceive (
    uhd::usrp::multi_usrp::sptr usrp,
    uhd::tx_streamer::sptr tx_stream,
    uhd::rx_streamer::sptr rx_stream,
    unsigned int npulses,
    float pulse_time,
    //std::complex<int16_t>* txbuff,
    std::vector<std::complex<int16_t> >* txbuff,
    float tx_ontime,
    std::complex<int16_t>** outdata,
    size_t samps_per_pulse
){
    //create metadeta tags for transmit streams
    uhd::time_spec_t start_time = usrp->get_time_now() + 0.05;
    uhd::tx_metadata_t md;
    md.start_of_burst = true;
    md.end_of_burst = false;
    md.has_time_spec = true;
    md.time_spec = start_time;
    std::vector<std::complex<int16_t> *> vec_ptr;
    vec_ptr.resize(1);
    vec_ptr[0] = &txbuff->front();

    usrp->set_gpio_attr("TXA","CTRL",0x0, 0x60);
    usrp->set_gpio_attr("TXA","DDR",0x60, 0x60);

    //for (int i=0; i<bufflen; i++){
    //    printf("tx %i: %i,%i\n",i,txbuff[i].real(), txbuff[i].imag());
    //}

    //create metadata tags for receive stream
    uhd::rx_metadata_t rxmd;
    std::vector<std::complex<int16_t> > buff(samps_per_pulse,0);
    if (verbose) std::cout << "buff size: " << buff.size() << std::endl;
    uhd::stream_cmd_t stream_cmd = uhd::stream_cmd_t::STREAM_MODE_NUM_SAMPS_AND_DONE;
    stream_cmd.num_samps = npulses*samps_per_pulse;
    stream_cmd.stream_now = false;
    stream_cmd.time_spec = start_time;//+txon;
    if (verbose) std::cout << "time spec: " << stream_cmd.time_spec.get_real_secs() << std::endl;
   
   //loop for every pulse in the sequence
   size_t spb;
   std::vector<std::complex<int16_t>* > rx_dptr;
   rx_dptr.resize(usrp->get_rx_num_channels());
   spb = tx_stream->get_max_num_samps();
   if (verbose) std::cout << "npulses: " << npulses << std::endl;
   boost::thread_group rx_threads;
   boost::thread_group tx_threads;
   for (int ipulse=0; ipulse<npulses; ipulse++){
       for (size_t ichan=0; ichan<usrp->get_rx_num_channels(); ichan++){
        rx_dptr[ichan] = ipulse*samps_per_pulse + outdata[ichan];
       }

       float timeout = 1.1;

       usrp->set_command_time(start_time-50e-6,0);
       usrp->set_gpio_attr("TXA","OUT",0x60, 0x60);

       if (ipulse==0){
         if (verbose) std::cout << "time spec: " << stream_cmd.time_spec.get_real_secs() << std::endl;
         if (verbose) std::cout << "Issuing stream command to start collecting samples\n";
         usrp->issue_stream_cmd(stream_cmd);
       }

       usrp->set_command_time(start_time+tx_ontime+20e-6,0);
       usrp->set_gpio_attr("TXA","OUT",0x0, 0x60);

       size_t acc_samps=0;
       vec_ptr[0] = &txbuff->front();
       
       tx_threads.join_all();
       if (ipulse != npulses-1) {
            tx_threads.create_thread(boost::bind(tx_worker,
                txbuff->size(), tx_stream, start_time, vec_ptr[0], 0));
       }
       if (ipulse == npulses-1) {
            tx_threads.create_thread(boost::bind(tx_worker,
                txbuff->size(), tx_stream, start_time-2e-3, vec_ptr[0], 1));
       }

       rx_threads.join_all();
       rx_threads.create_thread(boost::bind(rx_worker,
        rx_stream, samps_per_pulse, rx_dptr));

       //for (int k=0; k<10; k++){
       // //std::cout << "raw data: " << outdata[0][i][k] << "\t" << outdata[1][i][k] << std::endl;
       // std::cout << "raw data: " << rx_dptr[0][k] << " " << rx_dptr[1][k] << std::endl;
       //}
       //for (int k=0; k<samps_per_pulse; k++)
       //    outdata[i][k] += buff[k];


       start_time += float(pulse_time);
   }
}

/***********************************************************************
 * tx_worker function
 * A function to be used in its own thread for transmitting.  Push all
 * tx values into the USRP buffer as USRP buffer space is available,
 * but allow other actions to occur concurrently.
 **********************************************************************/
void tx_worker(
    unsigned int bufflen,
    uhd::tx_streamer::sptr tx_stream,
    uhd::time_spec_t start_time,
    std::complex<int16_t>* vec_ptr,
    int end
){
    unsigned int acc_samps = 0;

    uhd::tx_metadata_t md;
    md.start_of_burst = true;
    md.has_time_spec = true;
    md.time_spec = start_time;

    size_t spb = tx_stream->get_max_num_samps();
    if (spb > bufflen) spb = bufflen;

	while(acc_samps < bufflen-spb){
            size_t nsamples = tx_stream->send(vec_ptr, spb, md);
            vec_ptr += spb;
            acc_samps += nsamples;
            //std::cout << acc_samps <<std::endl;
            md.start_of_burst = false;
            md.has_time_spec = false;
    }
    // Now on the last packet
    if (end) md.end_of_burst = true;
    spb = bufflen - acc_samps;
    size_t nsamples = tx_stream->send(vec_ptr, spb, md);
}

void rx_worker(
    uhd::rx_streamer::sptr rx_stream,
    unsigned int samps_per_pulse,
    std::vector<std::complex<int16_t>* >& recv_ptr
){
    uhd::rx_metadata_t rxmd;
    float timeout = 1.1;
    rxmd.error_code = uhd::rx_metadata_t::ERROR_CODE_NONE;
    size_t nrx_samples = rx_stream->recv(recv_ptr, samps_per_pulse, rxmd, timeout);
    if (rxmd.error_code != uhd::rx_metadata_t::ERROR_CODE_NONE){
        std::cerr << "Error!\t";
        std::cerr << rxmd.error_code << std::endl;
    }
}
