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

    //
    //float p_code[13] = {1,1,1,1,1,-1,-1,1,1,-1,1,-1,1};
    //float p_code[1] = {1};
    float *p_code;//[3] = {1,1,-1};
    //float p_code[10] = {0,0,0,1,1,1,1,1,1,1};
    //int pcodelen = 3;
    std::vector<float> pcode;

    //status flags
    int new_seq_flag = 1;

    //transmit variables
    std::string tx_subdev="A:A";
    unsigned int bufflen;
    unsigned int samps_per_sym;
    boost::thread_group transmit_thread;
    std::vector<std::complex<float> > tx_raw_buff;
    std::vector<std::complex<int16_t> > tx_filt_buff;
    std::vector<std::complex<float> > filter_taps;
    int ntaps;
    float txsamprate, txbw;

    //receive variables
    std::string rx_subdev="A:A";
    size_t spb;
    double rx_rate; 
    unsigned int nave = 1;
    std::vector<std::vector<std::complex<float> > > outvecs;
    std::vector<std::complex<float> *> outvec_ptrs;
    float ptime_eff;

    //data processing variables
    int dmrate, slowdim, fastdim;
    float bandwidth;
    std::vector<std::vector<std::complex<float> > > filtvecs;
    std::vector<std::complex<float> *> filtvec_ptrs;
    std::vector<std::vector<std::complex<float> > > ffvecs;
    std::vector<std::complex<float> *> ffvec_ptrs;
    std::vector<float> fpow;
    std::vector<float> fvel;

    //socket-related variables
    int sock, msgsock, rval, rfds, efds;
    int msg;
    struct soundingParms parms;
    char usrpmsg;

    //usrp-related variables
    uhd::time_spec_t start_time;

    //create a usrp device
    args="addr=192.168.10.2";
    uhd::usrp::multi_usrp::sptr usrp = uhd::usrp::multi_usrp::make(args);

    //Lock device to the motherboard clock and 
    //then initialize to an arbitrary time
    usrp->set_clock_source("internal");
    usrp->set_time_now(uhd::time_spec_t(0.0)); 

    //create a transmit streamer
    uhd::stream_args_t tx_stream_args("sc16", "sc16");
    usrp->set_tx_subdev_spec(tx_subdev);
    uhd::tx_streamer::sptr tx_stream = 
    	usrp->get_tx_stream(tx_stream_args);
    
    //create a receive streamer
    uhd::stream_args_t rx_stream_args("sc16", "sc16");
    usrp->set_rx_subdev_spec(rx_subdev);
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
    std::vector<std::string> banks;
    while(true){
        rval = recv_data(msgsock, &usrpmsg, sizeof(usrpmsg));
        switch (usrpmsg){
            case EXIT:
                printf("Done sounding; quiting program\n");
                //banks = usrp->get_gpio_banks(0);
                //for (size_t i=0; i<banks.size(); i++){
                //    std::cout << banks[i] << std::endl;
                //}
                return 0;
            case SEND:
                printf("Starting sounding.\n");

                rval = recv_data(msgsock, &parms, sizeof(parms));
                printf("msg values\n");
                printf("freq: %f\n", parms.freq);
                printf("txrate: %f\n", parms.txrate);
                printf("rxrate: %f\n", parms.rxrate);
                std::cout << parms.pc_str << std::endl;
                if (strcmp(parms.pc_str,"barker13")==0){
                    std::cout << "using barker 13 pcode\n";
                   static const int arr[] = BARKER_13;
                   pcode.assign(arr, arr+sizeof(arr)/sizeof(arr[0]));
                } else if (strcmp(parms.pc_str, "rect")!=0){
                    std::cout << "using no pcode\n";
                   static const int arr[] = RECT;
                   pcode.assign(arr, arr+sizeof(arr)/sizeof(arr[0]));
                } else {
                   std::cerr << "invalid pulse code requested\n";
                   return 1;
                }
                for (int i=0; i<pcode.size(); i++){
                    std::cout << pcode[i] <<std::endl;
                }
                        
                        
                //pcode.assign(parms.pc_vec.begin(), parms.pc_vec.end());

	            freq = parms.freq;
	            //tx_rate = parms.txrate;
	            //rx_rate = parms.rxrate;
	            std::cout << "Using options: \n" << 
	            	freq << "\n" << parms.txrate << "\n" << rx_rate << std::endl;


	            //configure the USRP according to the arguments from the FIFO
	            freq /= 1e3; //convert to kHz
	            //recv_clr_freq(
	            //    usrp,
	            //    rx_stream,
	            //    &freq,
	            //    100);
	            freq *= 1000;
                usrp->set_rx_freq(freq);
                usrp->set_tx_freq(freq);
                usrp->set_rx_rate(parms.rxrate);
                usrp->set_tx_rate(parms.txrate);
	            
                std::cout << boost::format("Actual RX Freq: %f MHz...") % (usrp->get_rx_freq()/1e6) << "\n" << std::endl;

	            if (new_seq_flag == 1){
	            	std::cout << boost::format("Requesting Tx sample rate: %f Ksps...") % (parms.txrate/1e3) << std::endl;
	            	usrp->set_tx_rate(parms.txrate);
	            	std::cout << boost::format("Actual TX Rate: %f Ksps...") % (usrp->get_tx_rate()/1e3) << std::endl << std::endl;
	            	
	            	std::cout << boost::format("Requesting Rx sample rate: %f Ksps...") % (parms.rxrate/1e3) << std::endl;
	            	usrp->set_rx_rate(parms.rxrate);
	            	std::cout << boost::format("Actual RX Rate: %f Ksps...") % (usrp->get_rx_rate()/1e3) << std::endl << std::endl;
	            }
	            new_seq_flag = 0;


                stop_signal_called = false;

                //Prepare filter taps
                samps_per_sym = (unsigned int) (parms.symboltime * parms.txrate);
                ntaps = 4*samps_per_sym;
                filter_taps.resize(ntaps,0);
                txbw = 1/(2*parms.symboltime);
                txsamprate = usrp->get_tx_rate();
                for (int i=0; i<ntaps; i++){
                    double x=2*(2*M_PI*((float)i/ntaps)-M_PI);
                    filter_taps[i] = std::complex<float>(
                        //txbw*(0.54-0.46*cos((2*M_PI*((float)(i)+0.5))/ntaps))*sin(x)/(x)/txsamprate,
                        //0);
                        1*(0.54-0.46*cos((2*M_PI*((float)(i)+0.5))/ntaps))*sin(x)/(x),
                        0);
                }
                filter_taps[ntaps/2] = std::complex<float>(1,0);

                //for (int i=0; i<ntaps; i++){
                //    printf("filter_taps %i: %f, %f\n", i, filter_taps[i].real(), filter_taps[i].imag());
                //}
                
                //prepare raw tx information
                bufflen = pcode.size()* samps_per_sym;
                //printf("pcodelen: %i\tsamps_per_sym: %i\n", pcodelen,samps_per_sym);
                tx_raw_buff.resize(bufflen+ntaps,0);
                //for (int isym=ntaps/2+samps_per_sym/2; isym<pcodelen+ntaps/2+samps_per_sym/2; isym++){
                for (size_t isym=0; isym<pcode.size(); isym++){
                    //printf("isym: %i %i\n",isym, isym*samps_per_sym + ntaps/2 + samps_per_sym/2);
                    tx_raw_buff[isym*samps_per_sym+ntaps/2 + samps_per_sym/2] = std::complex<float>(
                        pcode[isym]*15000, 0x0000);
                }
                for (int i=0; i<tx_raw_buff.size(); i++){
                    //printf("tx_raw_buff %i: %f, %f\n", i, tx_raw_buff[i].real(), tx_raw_buff[i].imag());
                }

                //filter the raw tx vector
                tx_filt_buff.resize(bufflen,0);
                for (int i=0; i<bufflen; i++){
                    std::complex<float> temp(0,0);
                    for (int j=0; j<ntaps; j++){
                        temp += filter_taps[j] * tx_raw_buff[i+j];
                    }
                    tx_filt_buff[i] = std::complex<int16_t>((int16_t)temp.real(), (int16_t)temp.imag());
                    //printf("tx_filt_buff %i: %i, %i\n", i, tx_filt_buff[i].real(), tx_filt_buff[i].imag());
                }

                //prepare rx information
                ptime_eff = 3e8 / (2*MAX_VELOCITY*freq);
                printf("ptime_eff: %f\n", ptime_eff);
                nave = 1;
                ptime_eff /= 2;
                while(ptime_eff > parms.pulsetime && parms.npulses/nave > 1){
                    nave *= 2;
                    ptime_eff /= 2;
                }
                printf("nave: %i\n", nave);
                ptime_eff = nave*parms.pulsetime;
                printf("ptime_eff: %f\n", ptime_eff);

                slowdim = parms.npulses/nave;
                printf("slowdim: %i\n", slowdim);
                outvecs.resize(slowdim);
                outvec_ptrs.resize(slowdim);
                for (int i=0; i<slowdim; i++){
                    outvecs[i].resize(parms.nsamps_per_pulse,0);
                    outvecs[i].assign(parms.nsamps_per_pulse,0);
                    outvec_ptrs[i] = &outvecs[i].front();
                }

                //Set a start_time to be used for both transmit and receive usrp devices
                start_time =  usrp->get_time_now() + .1;
                
                //call function to spawn thread and transmit from file
                transceive(
                    usrp,
                    tx_stream,
                    start_time,
                    parms.npulses,
                    parms.pulsetime,
                    &tx_filt_buff.front(),
                    tx_filt_buff.size(),
                    rx_stream,
                    outvec_ptrs,
                    parms.nsamps_per_pulse,
                    nave);
                //transmit_thread.create_thread(boost::bind(tx_worker, 
	            //usrp, 
	            //tx_stream,
                //&tx_raw_buff.front(),
                //tx_raw_buff.size(),
                //parms.npulses,
                //parms.pulsetime,
	            //start_time));

                ////call function to receive to memory buffer
                //return_status = rx_worker(
	            //usrp,
	            //rx_stream,
	            //outvec_ptrs,
	            //parms.nsamps_per_pulse,
	            //parms.npulses,
	            //parms.pulsetime,
                //nave,
	            //start_time);

	            std::cout << "Done receiving, waiting for transmit thread.." << std::endl;
                transmit_thread.join_all();

	            if (return_status){
                    std::cerr << "This is a bad record..\n";
                    //Do something..?
                }
                send(msgsock, &return_status, sizeof(return_status),0);
                std::cout << "Done rxing\n";
                break;

            case PROCESS:
                std::cout << "Starting processing\n";
                dmrate = parms.symboltime * parms.rxrate/ 2; 
                fastdim = parms.nsamps_per_pulse/dmrate;
                bandwidth = 1/(2*parms.symboltime);
                printf("fastdim: %i\n",fastdim);
                printf("dmrate: %i\n", dmrate);
                printf("bandwidth: %f\n", bandwidth);
                filtvec_ptrs.resize(slowdim);
                filtvecs.resize(slowdim);
                for (int i=0; i<slowdim; i++){
                    filtvecs[i].resize(parms.nsamps_per_pulse/dmrate);
                    filtvec_ptrs[i] = &filtvecs[i].front();
                }


                rval = lp_filter(
                    outvec_ptrs,
                    filtvec_ptrs,
                    slowdim,
                    parms.nsamps_per_pulse,
                    (float) parms.rxrate,
                    bandwidth,
                    dmrate);

                ffvec_ptrs.resize(slowdim);
                ffvecs.resize(slowdim);
                for (int i=0; i<slowdim; i++){
                    ffvecs[i].resize(parms.nsamps_per_pulse/dmrate);
                    ffvec_ptrs[i] = &ffvecs[i].front();
                }
                rval = matched_filter(
                    filtvec_ptrs,
                    ffvec_ptrs,
                    &pcode.front(),
                    pcode.size(),
                    slowdim,
                    parms.nsamps_per_pulse/dmrate,
                    2);

                //typedef float rrec[2];
                fpow.resize(parms.nsamps_per_pulse/dmrate,0);
                fvel.resize(parms.nsamps_per_pulse/dmrate,0);
                rval = doppler_process(
                    ffvec_ptrs,
                    &fpow.front(),
                    &fvel.front(),
                    slowdim,
                    parms.nsamps_per_pulse/dmrate);

                send(msgsock, &return_status, sizeof(return_status),0);
        for (int i=0; i<parms.nsamps_per_pulse/dmrate; i++){
            fpow[i] /= (float(parms.npulses*parms.npulses)*std::pow(2,30)*std::pow(50,2))/16;
            //printf("%i: %.1f @ %.1f\n",i,30+10*log10(fpow[i]),fvel[i]);
            //printf("%i: %e @ %.1f\n",i,fpow[i],fvel[i]);
            //filtvec_ptrs[0][i] /= ((float)nave*std::pow(2,15));
            //ffvec_ptrs[0][i] /= ((float)nave*std::pow(2,15));
            //outvec_ptrs[0][i] /= ((float)nave*std::pow(2,15))/4;
            //printf("%i: %e @ %.1f\n",i,std::abs(filtvec_ptrs[0][i]),fvel[i]);
            //printf("%i: %e @ %.1f\n",i,std::abs(ffvec_ptrs[0][i]),fvel[i]);
            //printf("%i: %.3f @ %.1f\n",i,std::abs(outvec_ptrs[0][i]),fvel[i]);
        }
                break;

            case GET_DATA:
                send(msgsock, &fastdim, sizeof(fastdim),0);
                send(msgsock, &fpow.front(), fpow.size()*sizeof(float),0);
                break;

            default:
                std::cout << "not a valid usrpmsg.  Exiting.\n";
                return 1;
        }
	
    }
    return EXIT_SUCCESS;
}
