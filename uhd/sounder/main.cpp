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

//#include "c_utils.h"
#include "utils.hpp"
#include "../client/global_variables.h"

#define HOST_PORT 45001

namespace po = boost::program_options;
typedef std::complex<int16_t>  sc16;

/***********************************************************************
 * Signal interrupt handlers
 **********************************************************************/
static bool stop_signal_called = false;
void sig_int_handler(int){stop_signal_called = true;}

int verbose = 1;
int debug = 0;


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
    std::vector<float> pcode0;
    std::vector<float> pcode1;
    float* pcode_ptrs[2];

    //status flags
    int new_seq_flag = 1;

    //periodogram variables
    size_t center_freq_khz;
    size_t span_khz;
    std::vector<float> periodogram;

    //transmit variables
    std::string tx_subdev="A:A";
    unsigned int bufflen;
    unsigned int samps_per_sym;
    boost::thread_group transmit_thread;
    std::vector<std::complex<float> > tx_raw_buff0;
    std::vector<std::complex<float> > tx_raw_buff1;
    std::vector<std::complex<int16_t> > tx_filt_buff0;
    std::vector<std::complex<int16_t> > tx_filt_buff1;
    std::vector<std::complex<float> > filter_taps;
    int ntaps;
    float txsamprate, txbw, tx_ontime;
    size_t tx_ontime_usec;

    //receive variables
    std::string rx_subdev="A:A A:B";
    size_t spb;
    double rx_rate; 
    unsigned int nave = 1;
    std::vector<std::vector<std::complex<int16_t> > > rawvecs;
    std::vector<std::complex<int16_t> *> rawvec_ptrs;
    float ptime_eff;

    //data processing variables
    int dmrate, slowdim, fastdim;
    float bandwidth;

    std::vector<std::vector<std::complex<float> > > outvecs[2][2];
    std::vector<std::complex<float> *> outvec_ptrs[2][2];

    std::vector<std::vector<std::complex<float> > > filtvecs[2][2];
    std::vector<std::complex<float> *> filtvec_ptrs[2][2];
    std::complex<float>** filtvec_dptr[2][2];

    std::vector<std::vector<std::complex<float> > > ffvecs[2];
    std::vector<std::complex<float> *> ffvec_ptrs[2];
    std::vector<float> fpow[2];
    std::vector<float> fvel[2];

    //socket-related variables
    int sock, msgsock, rval, rfds, efds;
    int msg;
    struct soundingParms parms;
    struct periodogramParms lparms;
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
    usrp->set_tx_subdev_spec(tx_subdev);
    uhd::tx_streamer::sptr tx_stream = 
    	usrp->get_tx_stream(tx_stream_args);
    
    usrp->set_rx_subdev_spec(rx_subdev);

    //create a receive streamer
    if (verbose) std::cout << "number of channels: " << usrp->get_rx_num_channels() << std::endl;
    uhd::stream_args_t rx_stream_args("sc16", "sc16");
    for (size_t rx_chan = 0; rx_chan < usrp->get_rx_num_channels(); rx_chan++){
        rx_stream_args.channels.push_back(rx_chan); //linear mapping
    }

    uhd::rx_streamer::sptr rx_stream = 
    	usrp->get_rx_stream(rx_stream_args);

    if (verbose) std::cout << boost::format("Using Device: %s") % usrp->get_pp_string() << std::endl;
    
    //USRP is initialized;
    //now execute the tx/rx operations per arguments passed by the tcp socket
    sock = tcpsocket(HOST_PORT);
    if (verbose) printf("socket: %i\n",sock);
    while(true){
        listen(sock,1);
        //rval = 1;
        if (verbose) printf("Waiting for client connection..\n");
        msgsock = accept(sock, 0, 0);
        if (verbose) printf("Listening on socket %i\n", msgsock);
        //std::vector<std::string> banks;
        while(true){
            rval = recv_data(msgsock, &usrpmsg, sizeof(usrpmsg));
            if (rval <= 0) std::cout << "breaking..\n";
            if (rval <= 0) break;
            switch (usrpmsg){
                case EXIT:
                    if (verbose) printf("Client done\n");
                    //banks = usrp->get_gpio_banks(0);
                    //for (size_t i=0; i<banks.size(); i++){
                    //    std::cout << banks[i] << std::endl;
                    //}
                    //return 0;
                    break;
                case LISTEN:
                    if (verbose) printf("Starting Listening.\n");
                    rval = recv_data(msgsock, &lparms, sizeof(parms));
                    center_freq_khz = (lparms.end_freq_khz + lparms.start_freq_khz) / 2;
                    span_khz = lparms.end_freq_khz - lparms.start_freq_khz;
                    //if (span_khz > 400){
                    //    span_khz = 400;
                    //    lparms.start_freq_khz = 
                    if (verbose){
                        std::cout << "start freq: " << lparms.start_freq_khz << std::endl;
                        std::cout << "end freq: " << lparms.end_freq_khz << std::endl;
                        std::cout << "bandwidth: " << lparms.bandwidth_khz << std::endl;
                        std::cout << "center_freq_khz: " << center_freq_khz << std::endl;
                        std::cout << "span_khz: " << span_khz << std::endl;
                    }
                    periodogram.resize(span_khz);
                    capture_spectrum(
                        usrp,
                        rx_stream,
                        lparms.start_freq_khz,
                        lparms.end_freq_khz, 
                        lparms.bandwidth_khz,
                        &periodogram.front());

                    send(msgsock, &periodogram.front(), periodogram.size()*sizeof(float),0);
                    return_status=0;
                    send(msgsock, &return_status, sizeof(return_status),0);
                    break;

                case SEND:
                    if (verbose) printf("Starting sounding.\n");

                    rval = recv_data(msgsock, &parms, sizeof(parms));
                    if (verbose) printf("msg values\n");
                    if (verbose) printf("freq: %f\n", parms.freq);
                    //if (verbose) printf("txrate: %i\n", parms.txrate_khz);
                    //if (verbose) printf("rxrate: %i\n", parms.rxrate_khz);
                    if (verbose) std::cout << parms.pc_str << std::endl;
                    if (strcmp(parms.pc_str,"barker13") == 0){
                        if (verbose) std::cout << "using barker 13 pcode\n";
                        pcode0 = BARKER_13;
                        pcode1 = BARKER_13;
                    } 
                    else if (strcmp(parms.pc_str, "golay16") == 0){
                        if (verbose) std::cout << "using golay 16 pcode\n";
                        pcode0 = GOLAY_16_0;
                        pcode1 = GOLAY_16_1;
                    }
                    else if (strcmp(parms.pc_str, "golay10") == 0){
                        if (verbose) std::cout << "using golay 10 pcode\n";
                        pcode0 = GOLAY_10_0;
                        pcode1 = GOLAY_10_1;
                    }
                    else if (strcmp(parms.pc_str, "golay8") == 0){
                        if (verbose) std::cout << "using golay 8 pcode\n";
                        pcode0 = GOLAY_8_0;
                        pcode1 = GOLAY_8_1;
                    }
                    else if (strcmp(parms.pc_str, "golay4") == 0){
                        if (verbose) std::cout << "using golay 4 pcode\n";
                        pcode0 = GOLAY_4_0;
                        pcode1 = GOLAY_4_1;
                    }
                    else if (strcmp(parms.pc_str, "rect") == 0){
                        if (verbose) std::cout << "using no pcode\n";
                        pcode0 = RECT;
                        pcode1 = RECT;
                    }
                    else {
                        std::cerr << "invalid pulse code requested\n";
                        return 1;
                    }
                    for (int i=0; i<pcode0.size(); i++){
                        if (verbose) std::cout << pcode0[i] << " " << pcode1[i] << std::endl;
                    }
                            
	                freq = 1e3*parms.freq;
                    if (freq < 7e6){
                        if (verbose) std::cout << "Using antenna A\n";
                        usrp->set_tx_subdev_spec(std::string("A:A"));
                    }
                    else {
                        if (verbose) std::cout << "Using antenna B\n";
                        usrp->set_tx_subdev_spec(std::string("A:B"));
                    }

	                //configure the USRP according to the arguments from the client
                    for (size_t i=0; i<usrp->get_rx_num_channels(); i++){
                        usrp->set_rx_freq(freq,i);
                    }
                    usrp->set_tx_freq(freq);
                    usrp->set_rx_rate(RX_RATE);
                    usrp->set_tx_rate(TX_RATE);
	                
                    if (verbose) std::cout << boost::format("Actual RX Freq: %f MHz...") % (usrp->get_rx_freq()/1e6) << "\n" << std::endl;

	                if (new_seq_flag == 1){
	                	usrp->set_tx_rate(TX_RATE);
	                	if (verbose) std::cout << boost::format("Actual TX Rate: %f Ksps...") % (usrp->get_tx_rate()/1e3) << std::endl << std::endl;
	                	usrp->set_rx_rate(RX_RATE);
	                	if (verbose) std::cout << boost::format("Actual RX Rate: %f Ksps...") % (usrp->get_rx_rate()/1e3) << std::endl << std::endl;
	                }
	                new_seq_flag = 0;


                    stop_signal_called = false;

                    //Prepare lp filter taps for filtering the tx samples
                    //samps_per_sym = (unsigned int) (parms.symboltime * parms.txrate);
                    //samps_per_sym = parms.symboltime_usec * parms.txrate_khz / 1000;
                    samps_per_sym = parms.symboltime_usec * TX_RATE / 1e6;
                    ntaps = 4*samps_per_sym;
                    filter_taps.resize(ntaps,0);
                    txbw = 1/(2.e-6*parms.symboltime_usec);
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

                    if (debug){
                        for (int i=0; i<ntaps; i++){
                            printf("filter_taps %i: %f, %f\n", i, filter_taps[i].real(), filter_taps[i].imag());
                        }
                    }
                    
                    //prepare raw tx information
                    pcode_ptrs[0] = &pcode0.front();
                    pcode_ptrs[1] = &pcode1.front();


                    bufflen = pcode0.size()* samps_per_sym;
                    tx_raw_buff0.resize(bufflen+ntaps,0);
                    tx_raw_buff1.resize(bufflen+ntaps,0);
                    for (size_t isym=0; isym<pcode0.size(); isym++){
                        tx_raw_buff0[isym*samps_per_sym+ntaps/2 + samps_per_sym/2] = std::complex<float>(
                            pcode0[isym]*15000, 0x0000);
                        tx_raw_buff1[isym*samps_per_sym+ntaps/2 + samps_per_sym/2] = std::complex<float>(
                            pcode1[isym]*15000, 0x0000);
                    }
                    //for (int i=0; i<tx_raw_buff.size(); i++){
                    //    printf("tx_raw_buff %i: %f, %f\n", i, tx_raw_buff[i].real(), tx_raw_buff[i].imag());
                    //}

                    //filter the raw tx vector
                    tx_filt_buff0.resize(bufflen,0);
                    tx_filt_buff1.resize(bufflen,0);
                    for (int i=0; i<bufflen; i++){
                        std::complex<float> temp0(0,0);
                        std::complex<float> temp1(0,0);
                        for (int j=0; j<ntaps; j++){
                            temp0 += filter_taps[j] * tx_raw_buff0[i+j];
                            temp1 += filter_taps[j] * tx_raw_buff1[i+j];
                        }
                        tx_filt_buff0[i] = std::complex<int16_t>((int16_t)temp0.real(), (int16_t)temp0.imag());
                        tx_filt_buff1[i] = std::complex<int16_t>((int16_t)temp1.real(), (int16_t)temp1.imag());
                        //printf("tx_filt_buff %i: %i, %i\n", i, tx_filt_buff[i].real(), tx_filt_buff[i].imag());
                    }

                    //tx_ontime = (float) tx_filt_buff0.size() / (1.e3*parms.txrate_khz) + 100e-6;
                    tx_ontime = (float) tx_filt_buff0.size() / (TX_RATE) + 100e-6;
                    if (debug) std::cout << "tx_ontime: " << tx_ontime << std::endl;
                    //tx_ontime_usec = tx_filt_buff0.size() / TX_RATE / 1000 + 50;
                    tx_filt_buff0.resize((parms.pulsetime_usec * 1e-6*TX_RATE), 0);
                    tx_filt_buff1.resize((parms.pulsetime_usec * 1e-6*TX_RATE), 0);
                    if (verbose) std::cout << "pulsetime_usec: " << parms.pulsetime_usec << std::endl;
                    if (verbose) std::cout << "TX_RATE: " << TX_RATE << std::endl;
                    if (verbose) std::cout << "tx_filt_buffx size: " << parms.pulsetime_usec * TX_RATE << std::endl;
                    if (verbose) std::cout << "tx_filt_buffx size: " << (size_t) (1.e-6*parms.pulsetime_usec * TX_RATE) << std::endl;

                    //prepare rx information
                    ptime_eff = 3e8 / (2*MAX_VELOCITY*freq);
                    if (verbose) printf("ptime_eff: %f\n", ptime_eff);
                    nave = 2;
                    ptime_eff /= 2;
                    while(ptime_eff > 1.e-6*parms.pulsetime_usec && parms.npulses/nave > 1){
                        nave *= 2;
                        ptime_eff /= 2;
                    }
                    if (verbose) printf("nave: %i\n", nave);
                    ptime_eff = nave*parms.pulsetime_usec;
                    if (verbose) printf("ptime_eff: %f usec\n", ptime_eff);

                    slowdim = parms.npulses/nave;
                    if (verbose) printf("slowdim: %i\n", slowdim);

                    rawvecs.resize(usrp->get_rx_num_channels());
                    rawvec_ptrs.resize(usrp->get_rx_num_channels());
                    for (size_t i=0; i<rawvecs.size(); i++){
                        rawvecs[i].resize(parms.nsamps_per_pulse*parms.npulses);
                        rawvec_ptrs[i] = &rawvecs[i].front();
                    }

                    transceive(
                        usrp,
                        tx_stream,
                        rx_stream,
                        parms.npulses,
                        1.e-6*parms.pulsetime_usec,
                        &tx_filt_buff0,
                        &tx_filt_buff1,
                        tx_ontime,
                        &rawvec_ptrs.front(),
                        parms.nsamps_per_pulse
                        );


	                if (verbose) std::cout << "Done receiving, waiting for transmit thread.." << std::endl;
                    transmit_thread.join_all();

	                if (return_status){
                        std::cerr << "This is a bad record..\n";
                        //Do something..?
                    }
                    if (verbose) std::cout << "Done rxing 0\n";
                    send(msgsock, &return_status, sizeof(return_status),0);
                    if (verbose) std::cout << "Done rxing\n";
                    break;

                case PROCESS:
                    if (verbose) std::cout << "Starting processing\n";
                    //dmrate = parms.symboltime_usec * parms.rxrate_khz / (OSR*1000);
                    dmrate = (size_t) (1.e-6*parms.symboltime_usec * RX_RATE / OSR);
                    fastdim = parms.nsamps_per_pulse/dmrate;
                    bandwidth = 1/(2.e-6*parms.symboltime_usec);
	    	        if (verbose) printf("symbol time: %i usec\n", parms.symboltime_usec);
                    if (verbose) printf("fastdim: %i\n",fastdim);
                    if (verbose) printf("dmrate: %i\n", dmrate);
                    if (verbose) printf("bandwidth: %f\n", bandwidth);

                    outvecs[0][0].resize(slowdim);
                    outvecs[0][1].resize(slowdim);
                    outvecs[1][0].resize(slowdim);
                    outvecs[1][1].resize(slowdim);
                    outvec_ptrs[0][0].resize(slowdim);
                    outvec_ptrs[0][1].resize(slowdim);
                    outvec_ptrs[1][0].resize(slowdim);
                    outvec_ptrs[1][1].resize(slowdim);
                    for (int i=0; i<slowdim; i++){
                        outvecs[0][0][i].clear();
                        outvecs[0][1][i].clear();
                        outvecs[1][0][i].clear();
                        outvecs[1][1][i].clear();
                        outvecs[0][0][i].resize(parms.nsamps_per_pulse,0);
                        outvecs[0][1][i].resize(parms.nsamps_per_pulse,0);
                        outvecs[1][0][i].resize(parms.nsamps_per_pulse,0);
                        outvecs[1][1][i].resize(parms.nsamps_per_pulse,0);
                        outvec_ptrs[0][0][i] = &outvecs[0][0][i].front();
                        outvec_ptrs[0][1][i] = &outvecs[0][1][i].front();
                        outvec_ptrs[1][0][i] = &outvecs[1][0][i].front();
                        outvec_ptrs[1][1][i] = &outvecs[1][1][i].front();
                    }
                    //for (int i=0; i<rawvecs[1].size(); i++){
                    //    std::cout << i << " " << rawvecs[1][i] << std::endl;
                    //}
                    for (int i=0; i<slowdim; i++){
                        for (int j=0; j<nave; j++){
                            for (int k=0; k<parms.nsamps_per_pulse; k++){
                                if (j%2 == 0){
                                    outvecs[0][0][i][k] += 
                                        std::complex<int16_t>(1,0) * 
                                        (rawvecs[0][i*nave*parms.nsamps_per_pulse+j*parms.nsamps_per_pulse+k] + 
                                            std::complex<int16_t>(0,1) * 
                                            rawvecs[1][i*nave*parms.nsamps_per_pulse+j*parms.nsamps_per_pulse+k]);
                                    outvecs[1][0][i][k] += 
                                        std::complex<int16_t>(1,0) * 
                                        (rawvecs[0][i*nave*parms.nsamps_per_pulse+j*parms.nsamps_per_pulse+k] - 
                                            std::complex<int16_t>(0,1) * 
                                            rawvecs[1][i*nave*parms.nsamps_per_pulse+j*parms.nsamps_per_pulse+k]);
                                    //std::cout << j << " " <<
                                    //    std::complex<int16_t>(1,0) * 
                                    //    (rawvecs[0][i*nave*parms.nsamps_per_pulse+j*parms.nsamps_per_pulse+k] + 
                                    //        std::complex<int16_t>(0,-1) * 
                                    //        rawvecs[1][i*nave*parms.nsamps_per_pulse+j*parms.nsamps_per_pulse+k])
                                    //        << std::endl;
                                }
                                if (j%2 == 1){
                                    outvecs[0][1][i][k] += 
                                        std::complex<int16_t>(1,0) * 
                                        (rawvecs[0][i*nave*parms.nsamps_per_pulse+j*parms.nsamps_per_pulse+k] + 
                                            std::complex<int16_t>(0,1) * 
                                            rawvecs[1][i*nave*parms.nsamps_per_pulse+j*parms.nsamps_per_pulse+k]);
                                    outvecs[1][1][i][k] += 
                                        std::complex<int16_t>(1,0) * 
                                        (rawvecs[0][i*nave*parms.nsamps_per_pulse+j*parms.nsamps_per_pulse+k] - 
                                            std::complex<int16_t>(0,1) * 
                                            rawvecs[1][i*nave*parms.nsamps_per_pulse+j*parms.nsamps_per_pulse+k]);
                                    //std::cout <<  j << " " << 
                                    //    std::complex<int16_t>(1,0) * 
                                    //    (rawvecs[0][i*nave*parms.nsamps_per_pulse+j*parms.nsamps_per_pulse+k] + 
                                    //        std::complex<int16_t>(0,-1) * 
                                    //        rawvecs[1][i*nave*parms.nsamps_per_pulse+j*parms.nsamps_per_pulse+k])
                                    //        << std::endl;
                                }
                            }
                            if (debug){
                                //if (j == nave-1 || j == nave-2) {
                                    for (int k=0; k<parms.nsamps_per_pulse; k++){
                                        std::cout << j << " " << k << " ";
                                        std::cout << outvecs[0][0][i][k] << "\t";
                                        std::cout << outvecs[0][1][i][k] << "\t";
                                        std::cout << outvecs[1][0][i][k] << "\t";
                                        std::cout << outvecs[1][1][i][k] << std::endl;
                                    }
                                //}
                            }
                        }
                    }

                    filtvecs[0][0].resize(slowdim);
                    filtvecs[0][1].resize(slowdim);
                    filtvecs[1][0].resize(slowdim);
                    filtvecs[1][1].resize(slowdim);
                    filtvec_ptrs[0][0].resize(slowdim);
                    filtvec_ptrs[0][1].resize(slowdim);
                    filtvec_ptrs[1][0].resize(slowdim);
                    filtvec_ptrs[1][1].resize(slowdim);
                    for (int i=0; i<slowdim; i++){
                        filtvecs[0][0][i].resize(parms.nsamps_per_pulse/dmrate);
                        filtvecs[0][1][i].resize(parms.nsamps_per_pulse/dmrate);
                        filtvecs[1][0][i].resize(parms.nsamps_per_pulse/dmrate);
                        filtvecs[1][1][i].resize(parms.nsamps_per_pulse/dmrate);
                        filtvec_ptrs[0][0][i] = &filtvecs[0][0][i].front();
                        filtvec_ptrs[0][1][i] = &filtvecs[0][1][i].front();
                        filtvec_ptrs[1][0][i] = &filtvecs[1][0][i].front();
                        filtvec_ptrs[1][1][i] = &filtvecs[1][1][i].front();
                    }
                    for (int mode=0; mode<2; mode++){
                        for (int pcode=0; pcode<2; pcode++){
                            rval = lp_filter(
                                outvec_ptrs[mode][pcode],
                                filtvec_ptrs[mode][pcode],
                                slowdim,
                                parms.nsamps_per_pulse,
                                RX_RATE,
                                bandwidth,
                                dmrate);
                        }
                    }
                    filtvec_dptr[0][0] = &filtvec_ptrs[0][0].front();
                    filtvec_dptr[0][1] = &filtvec_ptrs[0][1].front();
                    filtvec_dptr[1][0] = &filtvec_ptrs[1][0].front();
                    filtvec_dptr[1][1] = &filtvec_ptrs[1][1].front();

                    //for (int a=0; a<2; a++){
                    //    for (int i=0; i<slowdim; i++){
                    //        for (int j=0; j<parms.nsamps_per_pulse/dmrate; j++){
                    //            std::cout << a << " " << j << " " << filtvec_dptr[a][i][j] << std::endl;
                    //        }
                    //    }
                    //}
                    ffvec_ptrs[0].resize(slowdim);
                    ffvec_ptrs[1].resize(slowdim);
                    ffvecs[0].resize(slowdim);
                    ffvecs[1].resize(slowdim);
                    for (int i=0; i<slowdim; i++){
                        ffvecs[0][i].resize(parms.nsamps_per_pulse/dmrate);
                        ffvecs[1][i].resize(parms.nsamps_per_pulse/dmrate);
                        ffvec_ptrs[0][i] = &ffvecs[0][i].front();
                        ffvec_ptrs[1][i] = &ffvecs[1][i].front();
                    }

                    /* Performed matched filtering for both o- and x-mode samples*/
                    for (int mode=0; mode<2; mode++){
                        rval = matched_filter(
                            filtvec_dptr[mode],
                            ffvec_ptrs[mode],
                            pcode_ptrs,
                            pcode0.size(),
                            slowdim,
                            parms.nsamps_per_pulse/dmrate,
                            OSR);
                    }

                    //typedef float rrec[2];
                    fpow[0].resize(parms.nsamps_per_pulse/dmrate,0);
                    fpow[1].resize(parms.nsamps_per_pulse/dmrate,0);
                    fvel[0].resize(parms.nsamps_per_pulse/dmrate,0);
                    fvel[1].resize(parms.nsamps_per_pulse/dmrate,0);
                    for (int mode=0; mode<2; mode++){
                        rval = doppler_process(
                            ffvec_ptrs[mode],
                            &fpow[mode].front(),
                            &fvel[mode].front(),
                            slowdim,
                            parms.nsamps_per_pulse/dmrate);
                    }
                    send(msgsock, &return_status, sizeof(return_status),0);
                    for (int i=0; i<parms.nsamps_per_pulse/dmrate; i++){
                        fpow[0][i] /= (float(parms.npulses*parms.npulses)*std::pow(2,30)*std::pow(50,2))/16;
                        fpow[1][i] /= (float(parms.npulses*parms.npulses)*std::pow(2,30)*std::pow(50,2))/16;
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
                    send(msgsock, &fpow[0].front(), fpow[0].size()*sizeof(float),0);
                    send(msgsock, &fpow[1].front(), fpow[1].size()*sizeof(float),0);
                    break;

                default:
                    if (verbose) std::cout << "not a valid usrpmsg.  Exiting.\n";
                    return 1;
            }
	    
        }
    }
    //return EXIT_SUCCESS;
}
