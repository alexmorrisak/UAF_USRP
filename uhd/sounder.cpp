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

/***********************************************************************
 * Signal handlers
 **********************************************************************/
static bool stop_signal_called = false;
void sig_int_handler(int){stop_signal_called = true;}

long GetFileSize(std::string &filename)  {
  struct stat filestatus;
  stat( filename.c_str(), &filestatus );
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
    uhd::time_spec_t time_step
	
){

    //std::cout << "TXWORKER: Creating transmit worker thread..." << std::endl;
    //std::cout << boost::format("Parameters: %d, %d") % samps_per_buff % start_time.get_real_secs() << std::endl;
    //std::cout << cpu_format << " " <<  wire_format << " " << file << std::endl;

    //create a transmit streamer
    uhd::stream_args_t stream_args(cpu_format, wire_format);
    uhd::tx_streamer::sptr tx_stream = usrp->get_tx_stream(stream_args);

    uhd::tx_metadata_t md;
    md.start_of_burst = false;
    md.end_of_burst = false;
    md.has_time_spec = true;
    md.time_spec = start_time;
    std::vector<samp_type> buff(samps_per_buff);

    std::ifstream infile(file.c_str(), std::ifstream::binary);

    //loop until the entire file has been read
    uhd::time_spec_t t0 = usrp->get_time_now();

    while(not md.end_of_burst and not stop_signal_called){

        infile.read((char*)&buff.front(), buff.size()*sizeof(samp_type));

        size_t num_tx_samps = infile.gcount()/sizeof(samp_type);

        md.end_of_burst = infile.eof();

        //std::cout << boost::format("TX WORKER: Timestamp on a single buffer: %d") % t0.get_real_secs() << std::endl;

        tx_stream->send(&buff.front(), num_tx_samps, md);
    }

    infile.close();
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
    uhd::time_spec_t start_time
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
    //std::cout << boost::format("Timeout in seconds: %f") % timeout << std::endl;

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
    while(not stop_signal_called and num_total_samps != num_requested_samples){
        //std::cout << boost::format("Rx Time spec: %d") % t0.get_real_secs() << std::endl;
	//std::cout << "yep: " << num_total_samps << std::endl;
        size_t num_rx_samps = rx_stream->recv(&buff.front(), buff.size(), md, timeout);
	t0 = usrp->get_time_now();
        //timeout = t0.get_real_secs() + 0.1; //small timeout for subsequent recv

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

    //std::cout << boost::format("Done receiving: %d") % t0.get_real_secs() << std::endl;
    outfile.close();
    //std::cout << boost::format("Done writing to file: %d") % t0.get_real_secs() << std::endl;
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

    //syncronize the clock for the two usrp devices by polling the motherboard clock
    // Uncomment the lines below to use an external clock, the *correct* way to sychronize 
    // the transmit and receive streams
    std::cout << boost::format("Setting device timestamp to 0...") << std::endl;
    //const uhd::time_spec_t last_pps_time = tx_usrp->get_time_last_pps();
    //while (last_pps_time == tx_usrp->get_time_last_pps()){
    //    //current_time = tx_usrp->get_time_now();
    //    current_time = tx_usrp->get_time_last_pps();
    //    std::cout << boost::format("last PPS: %d, current time: %d ") % last_pps_time.get_real_secs()  % current_time.get_real_secs() << std::endl;
    //    usleep(100000);
    //}
    //tx_usrp->set_time_next_pps(uhd::time_spec_t(0.0));
    //rx_usrp->set_time_next_pps(uhd::time_spec_t(0.0));
    tx_usrp->set_time_now(uhd::time_spec_t(0.0));
    rx_usrp->set_time_now(uhd::time_spec_t(0.0));

    //always select the subdevice first, the channel mapping affects the other settings
    if (vm.count("tx-subdev")) tx_usrp->set_tx_subdev_spec(tx_subdev);
    if (vm.count("rx-subdev")) rx_usrp->set_rx_subdev_spec(rx_subdev);

    std::cout << boost::format("Using Device: %s") % tx_usrp->get_pp_string() << std::endl;
    std::cout << boost::format("Using Device: %s") % rx_usrp->get_pp_string() << std::endl;

    //set the transmit sample rate
    //if (not vm.count("tx-rate")){
    //    std::cerr << "Please specify the transmit sample rate with --tx-rate" << std::endl;
    //    return ~0;
    //}
    std::cout << boost::format("Setting TX Rate: %f Msps...") % (tx_rate/1e6) << std::endl;
    tx_usrp->set_tx_rate(tx_rate);
    std::cout << boost::format("Actual TX Rate: %f Msps...") % (tx_usrp->get_tx_rate()/1e6) << std::endl << std::endl;

    //set the receive sample rate
    //if (not vm.count("rx-rate")){
    //    std::cerr << "Please specify the sample rate with --rx-rate" << std::endl;
    //    return ~0;
    //}
    std::cout << boost::format("Setting RX Rate: %f Msps...") % (rx_rate/1e6) << std::endl;
    rx_usrp->set_rx_rate(rx_rate);
    std::cout << boost::format("Actual RX Rate: %f Msps...") % (rx_usrp->get_rx_rate()/1e6) << std::endl << std::endl;

    //set the transmit center frequency
    //if (not vm.count("freq")){
    //    std::cerr << "Please specify the transmit center frequency with --freq" << std::endl;
    //    return ~0;
    //}

    //for(size_t chan = 0; chan < tx_usrp->get_tx_num_channels(); chan++) {
    std::vector<std::string> prog_args;
    std::string arg_str;
    boost::thread_group transmit_thread;
    while(true){
	std::cout << "Waiting for program options (tx_file, rx_file, samp_rate, freq)...:" << std::endl;
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

        //Check Ref and LO Lock detect
        //std::vector<std::string> tx_sensor_names, rx_sensor_names;
        //tx_sensor_names = tx_usrp->get_tx_sensor_names(0);
        //if (std::find(tx_sensor_names.begin(), tx_sensor_names.end(), "lo_locked") != tx_sensor_names.end()) {
        //    uhd::sensor_value_t lo_locked = tx_usrp->get_tx_sensor("lo_locked",0);
        //    std::cout << boost::format("Checking TX: %s ...") % lo_locked.to_pp_string() << std::endl;
        //    UHD_ASSERT_THROW(lo_locked.to_bool());
        //}
        //rx_sensor_names = rx_usrp->get_rx_sensor_names(0);
        //if (std::find(rx_sensor_names.begin(), rx_sensor_names.end(), "lo_locked") != rx_sensor_names.end()) {
        //    uhd::sensor_value_t lo_locked = rx_usrp->get_rx_sensor("lo_locked",0);
        //    std::cout << boost::format("Checking RX: %s ...") % lo_locked.to_pp_string() << std::endl;
        //    UHD_ASSERT_THROW(lo_locked.to_bool());
        //}

        //tx_sensor_names = tx_usrp->get_mboard_sensor_names(0);
        //if ((ref == "mimo") and (std::find(tx_sensor_names.begin(), tx_sensor_names.end(), "mimo_locked") != tx_sensor_names.end())) {
        //    uhd::sensor_value_t mimo_locked = tx_usrp->get_mboard_sensor("mimo_locked",0);
        //    std::cout << boost::format("Checking TX: %s ...") % mimo_locked.to_pp_string() << std::endl;
        //    UHD_ASSERT_THROW(mimo_locked.to_bool());
        //}
        //if ((ref == "external") and (std::find(tx_sensor_names.begin(), tx_sensor_names.end(), "ref_locked") != tx_sensor_names.end())) {
        //    uhd::sensor_value_t ref_locked = tx_usrp->get_mboard_sensor("ref_locked",0);
        //    std::cout << boost::format("Checking TX: %s ...") % ref_locked.to_pp_string() << std::endl;
        //    UHD_ASSERT_THROW(ref_locked.to_bool());
        //}

        //rx_sensor_names = rx_usrp->get_mboard_sensor_names(0);
        //if ((ref == "mimo") and (std::find(rx_sensor_names.begin(), rx_sensor_names.end(), "mimo_locked") != rx_sensor_names.end())) {
        //    uhd::sensor_value_t mimo_locked = rx_usrp->get_mboard_sensor("mimo_locked",0);
        //    std::cout << boost::format("Checking RX: %s ...") % mimo_locked.to_pp_string() << std::endl;
        //    UHD_ASSERT_THROW(mimo_locked.to_bool());
        //}
        //if ((ref == "external") and (std::find(rx_sensor_names.begin(), rx_sensor_names.end(), "ref_locked") != rx_sensor_names.end())) {
        //    uhd::sensor_value_t ref_locked = rx_usrp->get_mboard_sensor("ref_locked",0);
        //    std::cout << boost::format("Checking RX: %s ...") % ref_locked.to_pp_string() << std::endl;
        //    UHD_ASSERT_THROW(ref_locked.to_bool());
        //}

        //Set a start_time to be used for both transmit and receive usrp devices
        uhd::time_spec_t start_time =  tx_usrp->get_time_now() + .1;
        uhd::time_spec_t buffer_time =  tx_usrp->get_time_now() + .1;
        uhd::time_spec_t tx_time =  tx_usrp->get_time_now() + .1;
        uhd::time_spec_t time_step =  tx_usrp->get_time_now() + .1;

        uhd::time_spec_t t0 = tx_usrp->get_time_now();
        uhd::time_spec_t t1 = rx_usrp->get_time_now();
        std::cout << "tx, rx times: " << t0.get_real_secs() << t1.get_real_secs() << std::endl;

        stop_signal_called = false;

        //transmit from file
        transmit_thread.create_thread(boost::bind(&transmit_worker<std::complex<short> >, 
		tx_usrp, "sc16", otw, tx_file, spb, 
		start_time, buffer_time, tx_time, time_step));

        //recv to file
        recv_to_file<std::complex<short> >(rx_usrp, "sc16", otw, rx_file, spb, nsamps, start_time);

        //clean up transmit worker
        stop_signal_called = true;
        transmit_thread.join_all();

        long filesize;
        filesize = GetFileSize(tx_file);
        std::cout << "rx_file: " << rx_file << std::endl;


        //finished
        std::cout << std::endl << "Done!" << std::endl << std::endl;
  }
        return EXIT_SUCCESS;
}
