// main.cpp
// 
// Program for USRP radar operations
// Adapted from example txrx_loopback_to_file.cpp
// Alex Morris
// 06 Dec 2013
#include <sys/socket.h>

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

#include <sounder.hpp>

#include "c_utils.h"
#include "../client/global_variables.h"

#define HOST_PORT 45001

namespace po = boost::program_options;
typedef std::complex<int16_t>  sc16;

/***********************************************************************
 * Signal interrupt handlers
 **********************************************************************/
static bool stop_signal_called = false;
void sig_int_handler(int){stop_signal_called = true;}

int verbose;


/***********************************************************************
 * Main function
 **********************************************************************/
int UHD_SAFE_MAIN(int argc, char *argv[]){
    uhd::set_thread_priority_safe();

    int return_status = 0;

    //universal variables to be set by po
    std::string ref, otw, type, args;
    double freq;

    //variables for handling arguments from the FIFO
    std::vector<std::string> prog_args;
    std::string arg_str;

    //status flags
    int new_seq_flag = 1;

    //transmit variables
    std::string tx_subdev, tx_file;
    double tx_rate; 
    unsigned int bufflen;
    boost::thread_group transmit_thread;

    //receive variables
    std::string rx_subdev, rx_file;
    size_t spb;
    double rx_rate; 
    unsigned int nave = 1;
    std::vector<std::vector<std::complex<float> > > outvecs;
    std::vector<std::complex<float> *> outvec_ptrs;

    //socket-related variables
    int sock, msgsock, rval, rfds, efds;
    int msg;
    struct soundingParms parms;
    char usrpmsg;

    //create a usrp device
    args="addr=192.168.10.2";
    uhd::usrp::multi_usrp::sptr usrp = uhd::usrp::multi_usrp::make(args);

    //Lock device to the motherboard clock and 
    //then initialize to an arbitrary time
    usrp->set_clock_source("internal");
    usrp->set_time_now(uhd::time_spec_t(0.0)); 

    //create a transmit streamer
    uhd::stream_args_t tx_stream_args("sc16", "sc16");
    uhd::tx_streamer::sptr tx_stream = 
    	usrp->get_tx_stream(tx_stream_args);
    
    //create a receive streamer
    uhd::stream_args_t rx_stream_args("sc16", "sc16");
    uhd::rx_streamer::sptr rx_stream = 
    	usrp->get_rx_stream(rx_stream_args);

    std::cout << boost::format("Using Device: %s") % usrp->get_pp_string() << std::endl;
    
    //USRP is initialized;
    //now execute the tx/rx operations per arguments passed by the tcp socket
    sock = tcpsocket(HOST_PORT);
    printf("socket: %i\n",sock);
    listen(sock,1);
    rval = 1;
    printf("Waiting for client connection..\n");
    msgsock = accept(sock, 0, 0);
    printf("Listening on socket %i\n", msgsock);
    while(true){
        rval = recv_data(msgsock, &usrpmsg, sizeof(usrpmsg));
        if (usrpmsg == EXIT) {
            printf("Done sounding; quiting program\n");
            return 0;
        }
        printf("Starting sounding.\n");

        rval = recv_data(msgsock, &parms, sizeof(parms));
        printf("msg values\n");
        printf("txfile: %s\n", parms.txfile);
        printf("rxfile: %s\n", parms.rxfile);
        printf("freq: %f\n", parms.freq);
        printf("txrate: %i\n", parms.txrate);
        printf("rxrate: %i\n", parms.rxrate);

	    tx_file = parms.txfile;
	    rx_file = parms.rxfile;
	    freq = parms.freq;
	    tx_rate = parms.txrate;
	    rx_rate = parms.rxrate;
	    std::cout << "Using options: \n" << 
	    	tx_file << "\n" << rx_file << "\n" <<
	    	freq << "\n" << tx_rate << "\n" << rx_rate << std::endl;


	    //configure the USRP according to the arguments from the FIFO
	    freq /= 1e3; //convert to kHz
	    recv_clr_freq(
	        usrp,
	        rx_stream,
	        &freq,
	        100);
	    freq *= 1000;
        usrp->set_rx_freq(freq);
        usrp->set_tx_freq(freq);
        usrp->set_rx_rate(rx_rate);
        usrp->set_tx_rate(tx_rate);
	    
        std::cout << boost::format("Actual RX Freq: %f MHz...") % (usrp->get_rx_freq()/1e6) << "\n" << std::endl;

	    if (new_seq_flag == 1){
	    	std::cout << boost::format("Requesting Tx sample rate: %f Ksps...") % (tx_rate/1e3) << std::endl;
	    	usrp->set_tx_rate(tx_rate);
	    	std::cout << boost::format("Actual TX Rate: %f Ksps...") % (usrp->get_tx_rate()/1e3) << std::endl << std::endl;
	    	
	    	std::cout << boost::format("Requesting Rx sample rate: %f Ksps...") % (rx_rate/1e3) << std::endl;
	    	usrp->set_rx_rate(rx_rate);
	    	std::cout << boost::format("Actual RX Rate: %f Ksps...") % (usrp->get_rx_rate()/1e3) << std::endl << std::endl;
	    }
	    new_seq_flag = 0;


        stop_signal_called = false;

        //prepare tx information
        bufflen = (unsigned int) (1e-6*(float)parms.symboltime * tx_rate);

        //prepare rx information
        int ptime_eff = floor(1e3*3e8 / (2*MAX_VELOCITY*freq));
        printf("ptime_eff: %i\n", ptime_eff);
        nave = 1;
        ptime_eff /= 2;
        while(ptime_eff > parms.pulsetime){
            nave *= 2;
            ptime_eff /= 2;
        }
        printf("nave: %i\n", nave);
        ptime_eff = nave*parms.pulsetime;
        printf("ptime_eff: %i\n", ptime_eff);

        int slowdim = parms.npulses/nave;
        printf("slowdim: %i\n", slowdim);
        outvecs.resize(slowdim);
        outvec_ptrs.resize(slowdim);
        for (int i=0; i<slowdim; i++){
            outvecs[i].resize(parms.nsamps_per_pulse,0);
            outvec_ptrs[i] = &outvecs[i].front();
        }

        //Set a start_time to be used for both transmit and receive usrp devices
        uhd::time_spec_t start_time =  usrp->get_time_now() + .1;
        
        std::cout << "starting tx\n";
        //call function to spawn thread and transmit from file
        transmit_thread.create_thread(boost::bind(
	    tx_worker, 
	    usrp, 
	    tx_stream,
        bufflen,
        parms.npulses,
        parms.pulsetime,
	    start_time));

        std::cout << "starting rx\n";
        //call function to receive to memory buffer
        return_status = rx_worker(
	    usrp,
	    rx_stream,
	    outvec_ptrs,
	    parms.nsamps_per_pulse,
	    parms.npulses,
	    parms.pulsetime,
        nave,
	    start_time);

        std::cout << "starting lp prep\n";
        int dmrate = 1;
        std::vector<std::vector<std::complex<float> > > filtvecs;
        std::vector<std::complex<float> *> filtvec_ptrs;
        filtvec_ptrs.resize(slowdim);
        filtvecs.resize(slowdim);
        for (int i=0; i<slowdim; i++){
            filtvecs[i].resize(parms.nsamps_per_pulse/dmrate);
            filtvec_ptrs[i] = &filtvecs[i].front();
        }

        std::cout << "starting lp\n";
        int rval = lp_filter(
            outvec_ptrs,
            filtvec_ptrs,
            slowdim,
            parms.nsamps_per_pulse,
            (float) parms.rxrate,
            10e3,
            dmrate);

        std::vector<float> fdata;
        fdata.resize(parms.nsamps_per_pulse/dmrate);
        rval = doppler_process(
            filtvec_ptrs,
            &fdata.front(),
            slowdim,
            parms.nsamps_per_pulse/dmrate);



        //for (int i=0; i<parms.nsamps_per_pulse/dmrate; i+=10){
        //    printf("%i: ",i);
        //    for (int j=0; j<slowdim; j++){
        //        printf("%.1f ",std::abs(filtvec_ptrs[j][i]));
        //    }
        //    printf("\n");
        //}

        for (int i=0; i<parms.nsamps_per_pulse/dmrate; i+=10){
            printf("%i: %.1f \n",i,fdata[i]);
        }

	    if (return_status){
            std::cerr << "This is a bad record..\n";
            //Do something..?
        } else
        {
	        std::cout << "Done receiving, waiting for transmit thread.." << std::endl;
            transmit_thread.join_all();
            //
            //finished
            std::cout << "Done! File written to: " << rx_file << std::endl;
	    }
	
    }
    return EXIT_SUCCESS;
}
