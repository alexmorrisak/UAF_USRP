// main.cpp
// 
// Program for USRP radar operations
// Adapted from example txrx_loopback_to_file.cpp
// Alex Morris
// 06 Dec 2013

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

namespace po = boost::program_options;
typedef std::complex<int16_t>  sc16;

/***********************************************************************
 * Signal interrupt handlers
 **********************************************************************/
static bool stop_signal_called = false;
void sig_int_handler(int){stop_signal_called = true;}


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

        //Set a start_time to be used for both transmit and receive usrp devices
        uhd::time_spec_t start_time =  usrp->get_time_now() + .1;

        //call function to spawn thread and transmit from file
        transmit_thread.create_thread(boost::bind(
		transmit_worker, 
		usrp, 
		tx_stream,
		tx_file,
		&spb,
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
