// sounder.cpp
// Program for radar operations
// For transmitting and receiving data streams to/from the USRP

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

/***********************************************************************
 * Signal handlers
 **********************************************************************/
static bool stop_signal_called = false;
void sig_int_handler(int){stop_signal_called = true;}

long GetFileSize(std::string &filename)  {
  struct stat filestatus;
  stat(filename.c_str(), &filestatus);
  return filestatus.st_size;
  }

/***********************************************************************
 * transmit_worker function
 * A function to be used as a boost::thread_group thread for transmitting
 **********************************************************************/
template<typename samp_type> void transmit_worker(
    uhd::usrp::multi_usrp::sptr usrp,
    const std::string &cpu_format,
    const std::string &wire_format,
    const std::string &file,
    size_t samps_per_buff,
    uhd::time_spec_t start_time,
    uhd::time_spec_t buffer_time,
    uhd::time_spec_t tx_time,
    uhd::time_spec_t time_step,
    double tx_rate,
    int npulses
){
    uhd::stream_args_t stream_args(cpu_format, wire_format);
    uhd::tx_streamer::sptr tx_stream = usrp->get_tx_stream(stream_args);
    uhd::tx_metadata_t md;
    std::vector<samp_type> buff(samps_per_buff);
    std::vector<samp_type> null_values(time_step.get_real_secs() / tx_rate);
    for (size_t i = 0; i < (time_step.get_real_secs() / tx_rate); i++){
	null_values[i] = 0;
    };

    uhd::usrp::dboard_iface::sptr dboard = usrp->get_tx_dboard_iface();
    dboard->set_gpio_ddr(dboard->UNIT_TX, 0x0001, 0x0001);
    md.time_spec = start_time;
    for(int i = 0; i < npulses; i++){

	usrp->set_command_time(md.time_spec - buffer_time);
	dboard->set_gpio_out(dboard->UNIT_TX, 0x0001, 0x0001);

    	md.start_of_burst = false;
    	md.end_of_burst = false;
    	md.has_time_spec = true;
    	md.time_spec += time_step;
    	std::ifstream infile(file.c_str(), std::ifstream::binary);

    	//loop until the entire file has been read
    	while(not md.end_of_burst and not stop_signal_called){
    	    infile.read((char*)&buff.front(), buff.size()*sizeof(samp_type));
    	    size_t num_tx_samps = infile.gcount()/sizeof(samp_type);
    	    md.end_of_burst = infile.eof();
    	    tx_stream->send(&buff.front(), num_tx_samps, md);
    	}
	usrp->set_command_time(md.time_spec + tx_time);
	dboard->set_gpio_out(dboard->UNIT_TX, 0x0000, 0x0001);
	tx_stream->send(&null_values, null_values.size(), md);
    	infile.close();
    }
}


/***********************************************************************
 * recv_to_file function
 **********************************************************************/
template<typename samp_type> void recv_to_file(
    uhd::usrp::multi_usrp::sptr usrp,
    const std::string &cpu_format,
    const std::string &wire_format,
    const std::string &file,
    size_t samps_per_buff,
    int num_requested_samples,
    uhd::time_spec_t start_time,
    uhd::time_spec_t time_step,
    uhd::time_spec_t pulse_time
){
    int num_total_samps = 0;
    //create a receive streamer
    uhd::stream_args_t stream_args(cpu_format,wire_format);
    uhd::rx_streamer::sptr rx_stream = usrp->get_rx_stream(stream_args);

    uhd::rx_metadata_t md;
    std::vector<samp_type> buff(samps_per_buff);
    std::ofstream outfile(file.c_str(), std::ofstream::binary);
    bool overflow_message = true;
    float timeout = start_time.get_real_secs() + 0.5; //expected settling time + padding for first recv

    //setup streaming
    uhd::stream_cmd_t stream_cmd((num_requested_samples == 0)?
        uhd::stream_cmd_t::STREAM_MODE_NUM_SAMPS_AND_DONE:
        uhd::stream_cmd_t::STREAM_MODE_NUM_SAMPS_AND_DONE
    );
    stream_cmd.num_samps = num_requested_samples;
    stream_cmd.stream_now = false;
    stream_cmd.time_spec = start_time;
    std::cout << boost::format("Start time: %d") % start_time.get_real_secs() << std::endl;

    usrp->issue_stream_cmd(stream_cmd);
    md.error_code = uhd::rx_metadata_t::ERROR_CODE_NONE;
    uhd::time_spec_t t0 = usrp->get_time_now();
    while(num_requested_samples != num_total_samps){
        //std::cout << boost::format("Rx Time spec: %d") % t0.get_real_secs() << std::endl;
	//std::cout << "yep: " << num_total_samps << std::endl;
        size_t num_rx_samps = rx_stream->recv(&buff.front(), buff.size(), md, timeout);

        if (md.error_code == uhd::rx_metadata_t::ERROR_CODE_TIMEOUT) {
            std::cout << boost::format("Timeout while streaming") << std::endl;
            break;
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
                ) % (usrp->get_rx_rate()*sizeof(samp_type)/1e6);
            }
            continue;
        }
        if (md.error_code != uhd::rx_metadata_t::ERROR_CODE_NONE){
            throw std::runtime_error(str(boost::format(
                "Unexpected error code 0x%x"
            ) % md.error_code));
        }

        num_total_samps += num_rx_samps;
	//std::cout << num_total_samps << std::endl;

        outfile.write((const char*)&buff.front(), num_rx_samps*sizeof(samp_type));
    }

    outfile.close();
}


/***********************************************************************
 * Main function
 **********************************************************************/
int UHD_SAFE_MAIN(int argc, char *argv[]){
    uhd::set_thread_priority_safe();

    //universal variables to be set by po
    std::string ref, otw, type;
    double freq;

    //transmit variables to be set by po
    std::string tx_args, tx_subdev, tx_file;
    double tx_rate; 
    float prf;
    int npulses;

    //receive variables to be set by po
    std::string rx_args, rx_subdev, rx_file;
    size_t nsamps, spb;
    double rx_rate; 

    //setup the program options
    po::options_description desc("Allowed options");
    desc.add_options()
        ("help", "help message")
        ("tx-args", po::value<std::string>(&tx_args)->default_value(""), "uhd transmit device address args")
        ("rx-args", po::value<std::string>(&rx_args)->default_value(""), "uhd receive device address args")
        ("tx-subdev", po::value<std::string>(&tx_subdev)->default_value("A:A"), "uhd transmit device address")
        ("rx-subdev", po::value<std::string>(&rx_subdev)->default_value("A:A"), "uhd receive device address")
        ("tx-file", po::value<std::string>(&tx_file)->default_value("tx.dat"), "name of the file to read binary samples from")
        ("rx-file", po::value<std::string>(&rx_file)->default_value("rx.dat"), "name of the file to write binary samples to")
        ("type", po::value<std::string>(&type)->default_value("short"), "sample type in file: double, float, or short")
        ("nsamps", po::value<size_t>(&nsamps)->default_value(1500000), "total number of samples to receive")
        ("spb", po::value<size_t>(&spb)->default_value(1000), "samples per buffer, 0 for default")
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

    //create 2 usrp devices
    uhd::usrp::multi_usrp::sptr tx_usrp = uhd::usrp::multi_usrp::make(tx_args);
    uhd::usrp::multi_usrp::sptr rx_usrp = uhd::usrp::multi_usrp::make(rx_args);

    //Lock mboard clocks
    tx_usrp->set_clock_source(ref);
    rx_usrp->set_clock_source(ref);

    std::cout << boost::format("Setting device timestamp to 0...") << std::endl;
    tx_usrp->set_time_now(uhd::time_spec_t(0.0));
    rx_usrp->set_time_now(uhd::time_spec_t(0.0));

    //always select the subdevice first, the channel mapping affects the other settings
    if (vm.count("tx-subdev")) tx_usrp->set_tx_subdev_spec(tx_subdev);
    if (vm.count("rx-subdev")) rx_usrp->set_rx_subdev_spec(rx_subdev);

    std::cout << boost::format("Using Device: %s") % tx_usrp->get_pp_string() << std::endl;
    std::cout << boost::format("Using Device: %s") % rx_usrp->get_pp_string() << std::endl;

    std::cout << boost::format("Setting TX Rate: %f Msps...") % (tx_rate/1e6) << std::endl;
    tx_usrp->set_tx_rate(tx_rate);
    std::cout << boost::format("Actual TX Rate: %f Msps...") % (tx_usrp->get_tx_rate()/1e6) << std::endl << std::endl;

    std::cout << boost::format("Setting RX Rate: %f Msps...") % (rx_rate/1e6) << std::endl;
    rx_usrp->set_rx_rate(rx_rate);
    std::cout << boost::format("Actual RX Rate: %f Msps...") % (rx_usrp->get_rx_rate()/1e6) << std::endl << std::endl;

    std::vector<std::string> prog_args;
    std::string arg_str;
    boost::thread_group transmit_thread;
    while(true){
	std::cout << "Waiting for program options (tx_file rx_file samp_rate freq prf npulses)...:" << std::endl;
	std::getline (std::cin, arg_str);
	if (arg_str == "EXIT")
		{
			std::cout << "freq string is EXIT" << std::endl;
			return 0;
		}
	boost::split(prog_args, arg_str, boost::is_any_of("\t "));
	tx_file = prog_args[0].c_str();
	rx_file = prog_args[1].c_str();
	tx_rate = atoi(prog_args[2].c_str());
	rx_rate = atoi(prog_args[2].c_str());
	freq = atoi(prog_args[3].c_str());
	prf = atoi(prog_args[4].c_str());
	npulses = atoi(prog_args[5].c_str());
	std::cout << tx_file << std::endl;
	std::cout << rx_file << std::endl;
	std::cout << tx_rate << std::endl;
	std::cout << rx_rate << std::endl;
	std::cout << freq << std::endl;
        std::cout << boost::format("Setting TX Freq: %f MHz...") % (freq/1e6) << std::endl;
        tx_usrp->set_tx_freq(freq);
        std::cout << boost::format("Actual TX Freq: %f MHz...") % (tx_usrp->get_tx_freq()/1e6) << std::endl << std::endl;

        std::cout << boost::format("Setting RX Freq: %f MHz...") % (freq/1e6) << std::endl;
        rx_usrp->set_rx_freq(freq);
        std::cout << boost::format("Actual RX Freq: %f MHz...") % (rx_usrp->get_rx_freq()/1e6) << std::endl << std::endl;

        //Determine host buffer size for transmit and receive streams
        if (spb == 0) {
            std::cout << "Must specify samples per buffer!" << std::endl;
            return ~0;
        }

        //Create a start_time that will be used for both transmit and receive usrp devices
        uhd::time_spec_t start_time =  tx_usrp->get_time_now() + .1;
        uhd::time_spec_t tx_time =  tx_usrp->get_time_now() + .1;
        uhd::time_spec_t buffer_time =  tx_usrp->get_time_now() + .1;
        //uhd::time_spec_t time_step =  tx_usrp->get_time_now() + .1;

        stop_signal_called = false;
        uhd::time_spec_t time_step =  1. / prf;
        //transmit from file
        transmit_thread.create_thread(boost::bind(&transmit_worker<std::complex<short> >,
		tx_usrp, "sc16", otw, tx_file, spb,
		start_time, buffer_time, tx_time, time_step,
		tx_rate, npulses));

        //recv to file
	nsamps = npulses / prf;
        recv_to_file<std::complex<short> >(rx_usrp, "sc16", otw, rx_file, spb, nsamps,
		 start_time, time_step, tx_time);

        //clean up transmit worker
        stop_signal_called = true;
        transmit_thread.join_all();

        //long filesize;
        //filesize = GetFileSize(tx_file);
        std::cout << "rx_file: " << rx_file << std::endl;


        //finished
        std::cout << std::endl << "Done!" << std::endl << std::endl;
  }
        return EXIT_SUCCESS;
}
