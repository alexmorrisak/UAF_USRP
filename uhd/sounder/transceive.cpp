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

/***********************************************************************
 * tx_worker function
 * A function to be used as a thread for transmitting.  Push all
 * tx values into the USRP buffer as USRP buffer space is available,
 * but allow other actions to occur concurrently.
 **********************************************************************/
void tx_worker(
    unsigned int bufflen,
    unsigned int spb,
    uhd::tx_streamer::sptr tx_stream,
    uhd::time_spec_t start_time,
    std::complex<int16_t>* vec_ptr
){
    unsigned int acc_samps = 0;

    uhd::tx_metadata_t md;
    md.start_of_burst = true;
    md.has_time_spec = true;
    md.time_spec = start_time;

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
        md.end_of_burst = true;
        spb = bufflen - acc_samps;
        size_t nsamples = tx_stream->send(vec_ptr, spb, md);
}

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
){
    //create metadeta tags for transmit streams
    uhd::tx_metadata_t md;
    md.start_of_burst = true;
    md.end_of_burst = false;
    md.has_time_spec = true;
    md.time_spec = start_time;
    float txon = (float) bufflen / usrp->get_tx_rate();
    std::vector<std::complex<int16_t> *> vec_ptr;
    vec_ptr.resize(1);
    vec_ptr[0] = txbuff;

    usrp->set_gpio_attr("TXA","CTRL",0x0, 0x40);
    usrp->set_gpio_attr("TXA","DDR",0x40, 0x40);

    //for (int i=0; i<bufflen; i++){
    //    printf("tx %i: %i,%i\n",i,txbuff[i].real(), txbuff[i].imag());
    //}

    //create metadata tags for receive stream
    uhd::rx_metadata_t rxmd;
    std::vector<std::complex<int16_t> > buff(samps_per_pulse,0);
    printf("buff size: %i\n", buff.size());
    std::vector<std::complex<int16_t> > fbuff(bufflen,1000);

    uhd::stream_cmd_t stream_cmd = uhd::stream_cmd_t::STREAM_MODE_NUM_SAMPS_AND_DONE;
    stream_cmd.num_samps = samps_per_pulse;
    stream_cmd.stream_now = false;
    //stream_cmd.time_spec = start_time;//+txon;
    std::cout << "time spec: " << stream_cmd.time_spec.get_real_secs() << std::endl;
   
   //loop for every pulse in the sequence
   for (int i=0; i<npulses/nave; i++){
    rxmd.error_code = uhd::rx_metadata_t::ERROR_CODE_NONE;
    for (int j=0; j<nave; j++){
        float timeout = 1.1;

        usrp->set_command_time(start_time-50e-6,0);
        usrp->set_gpio_attr("TXA","OUT",0x40, 0x40);

        stream_cmd.time_spec = start_time;
	    usrp->issue_stream_cmd(stream_cmd);

        usrp->set_command_time(start_time+txon+10e-6,0);
        //usrp->set_command_time(start_time+1e-9,0);
        usrp->set_gpio_attr("TXA","OUT",0x0, 0x40);

        size_t acc_samps=0;
        size_t spb;
        spb = tx_stream->get_max_num_samps();
        vec_ptr[0] = txbuff;
        //std::thread tx(tx_worker,
            //bufflen, 100, tx_stream, start_time, vec_ptr[0]);
        tx_worker(bufflen, spb, tx_stream, start_time, vec_ptr[0]);

        //tx.join();

        size_t nrx_samples = rx_stream->recv(&buff.front(), samps_per_pulse, rxmd, timeout);
        if (rxmd.error_code != uhd::rx_metadata_t::ERROR_CODE_NONE){
            printf("error!\n");
            std::cout << rxmd.error_code << std::endl;
            std::cout << uhd::rx_metadata_t::ERROR_CODE_OVERFLOW << std::endl;
        }
        for (int k=0; k<samps_per_pulse; k++)
            outdata[i][k] += buff[k];


        start_time += float(pulse_time);
        //md.time_spec = start_time;
    }
   //for (int k=0; k<samps_per_pulse; k+=10)
   //    printf("outdata[%i][%i]: %.1f\n", i, k, std::abs(outdata[i][k]));
   }
}

