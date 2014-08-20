// sounder.cpp
// 
// Program for USRP radar operations
// Adapted from example txrx_loopback_to_file.cpp
// Alex Morris
// 06 Dec 2013
// 
// Copyright 2010-2012 Ettus Research LLC
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.
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

namespace po = boost::program_options;
typedef std::complex<int16_t>  sc16;

/***********************************************************************
 * File to get the length of the transmit sequence file
 **********************************************************************/
long get_file_size(std::string &filename)
{
  struct stat filestatus;
  stat(filename.c_str(), &filestatus);
  return filestatus.st_size;
}

/***********************************************************************
 * Signal interrupt handlers
 **********************************************************************/
static bool stop_signal_called = false;
void sig_int_handler(int){stop_signal_called = true;}

/***********************************************************************
 * transmit_worker function
 * A function to be used as a boost::thread_group thread for transmitting
 **********************************************************************/
template<typename samp_type> void transmit_worker(
    uhd::usrp::multi_usrp::sptr usrp,
    uhd::tx_streamer::sptr tx_stream,
    const std::string &file,
    size_t samps_per_buff,
    uhd::time_spec_t start_time
){
    //create metadeta tags for transmit streams
    uhd::tx_metadata_t md;
    md.start_of_burst = false;
    md.end_of_burst = false;
    md.has_time_spec = true;
    md.time_spec = start_time;

    //open binary file for baseband tx stream
    std::vector<samp_type> buff(samps_per_buff);
    std::ifstream infile(file.c_str(), std::ifstream::binary);

    //loop until the entire file has been read
    while(not md.end_of_burst and not stop_signal_called){
        infile.read((char*)&buff.front(), buff.size()*sizeof(samp_type));
        size_t num_tx_samps = infile.gcount()/sizeof(samp_type);
        md.end_of_burst = infile.eof();
        tx_stream->send(&buff.front(), num_tx_samps, md);
    }

    infile.close();
}

/***********************************************************************
 * clear_freq function
 **********************************************************************/
void recv_clr_freq(
    uhd::usrp::multi_usrp::sptr usrp,
    uhd::rx_streamer::sptr rx_stream,
    double *center_freq,
    double bandwidth,
    size_t samps_per_buff
){
    int num_total_samps = 0;
    int nave=10;
    int nsamps = nave*(int)bandwidth;
    fftw_complex in=NULL,out=NULL;
    fftw_plan plan;
    
  
    uhd::rx_metadata_t md;
    std::vector<std::complex<double> > buff(nsamps/nave);
    bool overflow_message = true;
    float timeout = 0.2;

    //setup streaming
    uhd::stream_cmd_t stream_cmd = uhd::stream_cmd_t::STREAM_MODE_NUM_SAMPS_AND_DONE;
    stream_cmd.num_samps = nsamps;
    stream_cmd.stream_now = true;

    usrp->issue_stream_cmd(stream_cmd);
    md.error_code = uhd::rx_metadata_t::ERROR_CODE_NONE;
    while(not stop_signal_called and num_total_samps != nsamps){
	timeout = 0.2;
        size_t num_rx_samps = rx_stream->recv(&buff.front(), buff.size(), md, timeout);
        if (md.error_code == uhd::rx_metadata_t::ERROR_CODE_TIMEOUT) {
            std::cout << boost::format("Timeout while streaming") << std::endl;
            break;
        }
        if (md.error_code == uhd::rx_metadata_t::ERROR_CODE_OVERFLOW){
            continue;
        }
        if (md.error_code != uhd::rx_metadata_t::ERROR_CODE_NONE){
            throw std::runtime_error(str(boost::format(
                "Unexpected error code 0x%x"
            ) % md.error_code));
        }
        num_total_samps += num_rx_samps;
    
    if(in!=NULL){free(in);in=NULL;}
    in =(fftw_complex*)fftw_malloc(sizeof(fftw_complex)*bandwidth);
    if(out!=NULL){free(out);out=NULL;}
    out =(fftw_complex*)fftw_malloc(sizeof(fftw_complex)*bandwidth);

    plan = fftw_plan_dft_1d(bandwidth, in, out, FFTW_FORWARD, FFTW_ESTIMATE);

    for (int i=0;i<bandwidth;i++){
       //in[i][0] = (hann_window[i]*rx_short_vecs[0][i].real());
       //in[i][1] = (hann_window[i]*rx_short_vecs[0][i].imag());
       in[i][0] = (int)buff[i].real();
       in[i][1] = (int)buff[i].imag();
    }
	
    fftw_execute(plan);

     //Add the current spectrum to the running total
     for (int i=0;i<bandwidth;i++){
             tmp_pwr[i] += out[i][0]*out[i][0] + out[i][1]*out[i][1] / 
                     double(bandwidth*bandwidth*naverages);
     }
     //Done adding spectrum to running total

    fftw_destroy_plan(plan);
    if (in!=NULL) {free(in); in=NULL;}
    if (out!=NULL) {free(out); out=NULL;}

   }

  //Center the fft (fftshift)
  for(int i=0;i<(bandwidth/2);i++){
      pwr[bandwidth/2+i]=tmp_pwr[i];
      pwr[i]=tmp_pwr[bandwidth/2+i];
  }
  //Done centering

  if (tmp_pwr!=NULL) free(tmp_pwr);

}

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
    std::vector<sc16> buff(samps_per_buff);
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
    while(not stop_signal_called and num_total_samps < num_requested_samples){
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
                ) % (usrp->get_rx_rate()*sizeof(sc16)/1e6);
            }
            return -1;
        }
        if (md.error_code != uhd::rx_metadata_t::ERROR_CODE_NONE){
            std::cerr << "Unexpected error code 0x" << 
		md.error_code << std::endl;
	   return -1;
        }

        num_total_samps += num_rx_samps;

        outfile.write((const char*)&buff.front(), num_rx_samps*sizeof(sc16));
    }

    outfile.close();
    std::cout <<  "RXFILE:" << file.c_str() << std::endl;
    return 0;
}

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
    boost::thread_group transmit_thread;

    //receive variables
    std::string rx_subdev, rx_file;
    size_t nsamps, spb;
    double rx_rate; 

    //setup the program options
    po::options_description desc("Allowed options");
    desc.add_options()
        ("help", "help message")
        ("args", po::value<std::string>(&args)->default_value(""), "uhd device ip address args")
        ("tx-subdev", po::value<std::string>(&tx_subdev)->default_value("A:A"), "uhd transmit device address")
        ("rx-subdev", po::value<std::string>(&rx_subdev)->default_value("A:A"), "uhd receive device address")
        ("tx-file", po::value<std::string>(&tx_file)->default_value("tx.dat"), "name of the file to read binary samples from")
        ("rx-file", po::value<std::string>(&rx_file)->default_value("rx.dat"), "name of the file to write binary samples to")
        ("type", po::value<std::string>(&type)->default_value("short"), "host-side sample type: double, float, or short")
        ("nsamps", po::value<size_t>(&nsamps)->default_value(150000), "total number of samples to receive")
        ("spb", po::value<size_t>(&spb)->default_value(1000), "host-side memory buffer size")
        ("tx-rate", po::value<double>(&tx_rate)->default_value(250000), "rate of transmit outgoing samples")
        ("rx-rate", po::value<double>(&rx_rate)->default_value(250000), "rate of receive incoming samples")
        ("freq", po::value<double>(&freq), "Transmit and receive RF center frequency in Hz")
        ("ref", po::value<std::string>(&ref)->default_value("internal"), "clock reference (internal, external, mimo)")
        ("otw", po::value<std::string>(&otw)->default_value("sc16"), "specify the over-the-wire sample mode")
    ;
    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, desc), vm);
    po::notify(vm);

    //print the help message
    if (vm.count("help")){
        std::cout << boost::format("Sounder: %s") % desc << std::endl;
        return ~0;
    }

    //create a usrp device
    uhd::usrp::multi_usrp::sptr usrp = uhd::usrp::multi_usrp::make(args);

    //Lock device to the motherboard clock and 
    //then initialize to an arbitrary time
    usrp->set_clock_source("internal");
    usrp->set_time_now(uhd::time_spec_t(0.0)); 

    if (vm.count("tx-subdev")) usrp->set_tx_subdev_spec(tx_subdev);
    if (vm.count("rx-subdev")) usrp->set_rx_subdev_spec(rx_subdev);

    std::cout << boost::format("Using Device: %s") % usrp->get_pp_string() << std::endl;
    //USRP is initialized;
    //now execute the tx/rx operations per arguments passed by the FIFO
    while(true){
	//Grab one line of commands from the FIFO
	std::cout << "Waiting for program options " <<
		"(tx_file, rx_file, nsamples, freq, tx_rate, rx_rate)...:" <<
		std::endl;
	std::getline(std::cin, arg_str);
	if (arg_str == "EXIT"){
		std::cout << "Freq string is EXIT.  Quitting program" << std::endl;
		std::cout << "END: End the process.py why don't you.." << std::endl;
		return 0;
	}
	boost::split(prog_args, arg_str, boost::is_any_of("\t "));
	tx_file = prog_args[0].c_str();
	rx_file = prog_args[1].c_str();
	nsamps = atoi(prog_args[2].c_str());
	freq = atoi(prog_args[3].c_str());
	tx_rate = atoi(prog_args[4].c_str());
	rx_rate = atoi(prog_args[5].c_str());
	std::cout << "Using options: \n" << 
		tx_file << "\n" << rx_file << "\n" <<
		freq << "\n" << tx_rate << "\n" << rx_rate << std::endl;

	std::cout << "nsamps:" << nsamps << std::endl;
	//configure the USRP according to the arguments from the FIFO
        std::cout << boost::format("Requesting Tx center frequency: %f MHz...") % (freq/1e6) << std::endl;
        usrp->set_tx_freq(freq);
        std::cout << boost::format("Actual Tx freq: %f MHz...") % (usrp->get_tx_freq()/1e6) << "\n" << std::endl;

        std::cout << boost::format("Requesting Rx center frequency: %f MHz...") % (freq/1e6) << std::endl;
        usrp->set_rx_freq(freq);
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

        //Set a start_time to be used for both transmit and receive usrp devices
        uhd::time_spec_t start_time =  usrp->get_time_now() + .1;

        stop_signal_called = false;

    	//create a transmit streamer
    	uhd::stream_args_t tx_stream_args("sc16", "sc16");
    	uhd::tx_streamer::sptr tx_stream = 
		usrp->get_tx_stream(tx_stream_args);

    	//create a receive streamer
    	uhd::stream_args_t rx_stream_args("sc16", "sc16");
    	uhd::rx_streamer::sptr rx_stream = 
		usrp->get_rx_stream(rx_stream_args);

        //call function to spawn thread and transmit from file
        transmit_thread.create_thread(boost::bind(
		&transmit_worker<std::complex<short> >, 
		usrp, 
		tx_stream,
		tx_file,
		spb,
		start_time));

        //call function to receive to file
        return_status = recv_to_file(
		usrp,
		rx_stream,
		rx_file,
		spb,
		start_time,
		nsamps);

	if (return_status) std::cerr << "This is a bad record..\n";
	if (tx_file == "NONE") {
		std::cout << "done receiving spectrum file\n";
		std::cout << "RXFILE_SPECT:" << rx_file << std::endl;
	}
	else{
	std::cout << "Done receiving, waiting for transmit thread.." << std::endl;

        //Wait for the transmit thread
        stop_signal_called = true;
        transmit_thread.join_all();

        //finished
        std::cout << "File written to: " << rx_file << std::endl;
        std::cout << std::endl << "Done!" << std::endl << std::endl;
	}
	
  }
        return EXIT_SUCCESS;
}
