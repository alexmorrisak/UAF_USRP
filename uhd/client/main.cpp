#include <sys/socket.h>
#include <sys/types.h>
#include <iostream>
#include <string>
#include <fstream>
#include <string.h>
#include <vector>
#include <complex>
#include <string>
#include <iostream>

#include <netinet/in.h> //for sockaddr_in
#include <arpa/inet.h> //for inet_addr()
#include <unistd.h> //for close()

#include <errno.h>

#include <boost/program_options.hpp>
#include <boost/format.hpp>

#include "global_variables.h"
#include "hdf5.h"

namespace po = boost::program_options;

int main(int argc, char *argv[]){
    int verbose =1;
    /* Socket and status variables*/
    int sockfd=0, n=0;
    int rval =0, status=0;
    struct sockaddr_in serv_addr;
    struct soundingParms2 parms;
    struct periodogramParms spect_parms;
    char usrpmsg;

    /* hdf5 and file-writing variables */
    char dset[80];
    hid_t file_id, dataspace_id, dataset_id, attribute_id;
    herr_t eval = 0;
    hsize_t dims[2];
    time_t rawtime;
    struct tm * timeinfo;
    time(&rawtime);
    timeinfo = localtime(&rawtime);
    //std::cout << rawtime << std::endl;
    char fname[80];

    //Data storage variables
    std::vector<float> rxpow[2];
    std::vector<float> rxvel[2];
    std::vector<std::vector<float> > omode;
    std::vector<std::vector<float> > xmode;
    std::vector<float*> omode_ptrs;
    std::vector<float*> xmode_ptrs;
    std::vector<uint32_t> ranges;
    std::vector<uint32_t> frequencies;

    //Variables for testing
    std::vector<float> test;

    //variables for clear frequency search
    std::vector<float> spectrum;
    size_t min_inx;
    size_t nominal_freq;

    //variables for command line arguments
    float start_freq;
    float step_freq;
    float stop_freq;
    unsigned int nsteps;
    unsigned int npulses;
    float resolution;
    size_t osr;
    size_t first_range,last_range;
    int ifreq;

    po::options_description desc("Allowed options");
    desc.add_options()
        ("help", "help message")
        ("start", po::value<float>(&start_freq)->default_value(2e3),
            "Start frequency of swept-frequency ionogram")
        ("stop", po::value<float>(&stop_freq)->default_value(2e3),
            "Stop frequency of swept-frequency ionogram")
        ("nsteps", po::value<unsigned int>(&nsteps)->default_value(1),
            "Number of frequencies to step through between start and stop")
        ("npulses", po::value<unsigned int>(&npulses)->default_value(128),
            "Number of pulses in each integration period")
        ("resolution", po::value<float>(&resolution)->default_value(5),
            "Desired range resolution in km")
        ("osr", po::value<size_t>(&osr)->default_value(1),
            "Desired range oversampling for super-resolution")
        ("last-range", po::value<size_t>(&last_range)->default_value(750),
            "Maximum unambiguous range in km")
        ("first-range", po::value<size_t>(&first_range)->default_value(50),
            "First receivable range in km")
        ("write","write to file")
    ;
    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, desc), vm);
    po::notify(vm);

    if (vm.count("help")){
        std::cout << desc << "\n";
        return 1;
    }

    parms.num_pulses = npulses;
    parms.first_range_km = first_range;
    parms.last_range_km = last_range;
    parms.over_sample_rate = osr;
    parms.range_res_km = resolution;

    //size_t temp_symboltime_usec = resolution / 1.5e-1;
    //size_t temp_dmrate = (size_t) ceil(temp_symboltime_usec * RX_RATE / (OSR*1000));
    ////if (temp_dmrate%2 == 1) temp_dmrate -= 1;
    //while (temp_dmrate%OSR != 0) temp_dmrate -= 1;
    ////parms.symboltime_usec = OSR*1000*temp_dmrate/RX_RATE;


    //size_t temp_pfactor = (size_t) ceil(2*last_range / 3.0e-1 / parms.symboltime_usec);
    //parms.pulsetime_usec = temp_pfactor * parms.symboltime_usec;

    //std::cout << "Using pulse time: " << parms.pulsetime_usec << std::endl;
    //std::cout << "symbol time: " << parms.symboltime_usec << std::endl;
    //std::cout << "Using range resolution: " << 1.5e-1*parms.symboltime_usec << std::endl;


    //size_t max_txtime_usec = (size_t) floor(2*first_range / 3.0e-1);
    //size_t max_code_length = (size_t) floor(max_txtime_usec / parms.symboltime_usec);
    //if (max_code_length < 4){
    //    sprintf(parms.pc_str,"rect");
    //}
    //else if (max_code_length < 8){
    //    sprintf(parms.pc_str,"golay4");
    //}
    //else if (max_code_length < 10){
    //    sprintf(parms.pc_str,"golay8");
    //}
    //else if (max_code_length < 16){
    //    sprintf(parms.pc_str,"golay10");
    //}
    //else{
    //    sprintf(parms.pc_str,"golay16");
    //}
    //std::cout << "Using pcode: " << parms.pc_str << std::endl;
        

    //parms.pulsetime_usec = 5000;

    //sprintf(parms.pc_str,"golay8");
    //parms.nsamps_per_pulse = (size_t) (1e-6*parms.pulsetime_usec*RX_RATE);
    size_t datalen;

    printf("\nmsg values\n");
    //printf("freq: %i kHz\n", nominal_freq);
    //printf("nsamps per pulse: %i\n", parms.nsamps_per_pulse);

    if ( (sockfd = socket(AF_INET, SOCK_STREAM, 0)) == -1){
        printf("Error in creating socket\n");
        return 1;
    }

    serv_addr.sin_family = AF_INET;
    serv_addr.sin_addr.s_addr = inet_addr("127.0.0.1");
    serv_addr.sin_port = htons(HOST_PORT);

    if ( (connect(sockfd, (struct sockaddr *)&serv_addr, sizeof(serv_addr))) == -1){
        printf("bind failed\n");
        std::cout << strerror(errno) << std::endl;
        return 1;
    }
    printf ("connected!\n");

    strftime(fname, 80, "ionogram.%Y%m%d.%H%M.h5",timeinfo);
    if (vm.count("write"))
        file_id = H5Fcreate(fname, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

    step_freq = (stop_freq - start_freq) / nsteps;
    nominal_freq = start_freq;
    ifreq = 0;
    while(nominal_freq < stop_freq){
        /* Get the spectrum and then select the quietest frequency*/
        if (verbose) std::cout << "Grabbing spectrum\n";
        usrpmsg = LISTEN;
        send(sockfd, &usrpmsg, sizeof(usrpmsg), 0);
        size_t search_range = (int) step_freq / 2;
        if (search_range%2 !=0) search_range-=1;
        spectrum.resize(search_range);
        spect_parms.start_freq_khz = nominal_freq-search_range/2;
        spect_parms.end_freq_khz = nominal_freq+search_range/2;
        spect_parms.bandwidth_khz = 10;
        send(sockfd, &spect_parms, sizeof(spect_parms), 0);
        recv(sockfd, &spectrum.front(), search_range*sizeof(float), 0);
        min_inx = 0;
        for (int i=1; i<spectrum.size(); i++){
            if (spectrum[i] < spectrum[min_inx]) min_inx = i;
        }
        parms.freq_khz = min_inx + spect_parms.start_freq_khz;
        frequencies.push_back(parms.freq_khz);
        rval = recv(sockfd, &status, sizeof(int),0);

        /* Perform the sounding */
        if (verbose) std::cout << "Performing sounding\n";
        usrpmsg = SEND;
        send(sockfd, &usrpmsg, sizeof(usrpmsg), 0);
        send(sockfd, &parms, sizeof(parms), 0);
        rval = recv(sockfd, &status, sizeof(int),0);

        /* Process the data */
        if (verbose) std::cout << "Processing data\n";
        usrpmsg = PROCESS;
        send(sockfd, &usrpmsg, sizeof(usrpmsg), 0);
        rval = recv(sockfd, &status, sizeof(int),0);

        /* Get the data and append it to running data structure*/
        if (verbose) std::cout << "Getting data\n";
        usrpmsg = GET_DATA;
        send(sockfd, &usrpmsg, sizeof(usrpmsg), 0);
        recv(sockfd, &datalen, sizeof(datalen), 0);
        if (verbose) std::cout << "data length: " << datalen << std::endl;
        recv(sockfd, &first_range, sizeof(first_range), 0);
        if (verbose) std::cout << "first range: " << first_range << std::endl;
        recv(sockfd, &last_range, sizeof(last_range), 0);
        if (verbose) std::cout << "last range: " << last_range << std::endl;
        rxpow[0].resize(rxpow[0].size()+datalen);
        rxpow[1].resize(rxpow[1].size()+datalen);
        rxvel[0].resize(rxvel[0].size()+datalen);
        rxvel[1].resize(rxvel[1].size()+datalen);
        recv(sockfd, &rxpow[0].front() + (ifreq*datalen), datalen*sizeof(float), 0);
        recv(sockfd, &rxpow[1].front() + (ifreq*datalen), datalen*sizeof(float), 0);
        recv(sockfd, &rxvel[0].front() + (ifreq*datalen), datalen*sizeof(float), 0);
        recv(sockfd, &rxvel[1].front() + (ifreq*datalen), datalen*sizeof(float), 0);

        printf("Sounding %i of %i complete (%i kHz)\n", ifreq+1, nsteps, parms.freq_khz);
        nominal_freq += step_freq;
	    ifreq++;
    }

    ranges.resize(datalen);
    for (size_t i=0; i<ranges.size(); i++){
        ranges[i] = first_range + i*(last_range-first_range)/ranges.size();
    }

    /* Write the data to hdf5 file */
    if (vm.count("write")){
        /* Create 2-D dataset for omode and xmode powers.*/
        dims[0] = nsteps;
        dims[1] = datalen;
        dataspace_id = H5Screate_simple(2, dims, NULL);

        /* Write omode power to file */
        printf("OPower");
        sprintf(dset, "OPower");
        dataset_id = H5Dcreate2(file_id, dset, H5T_IEEE_F32BE, dataspace_id,
            H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

        eval =H5Dwrite(dataset_id, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
            &rxpow[0].front());
        if (eval) std::cerr << "Error writing to dataset: " << dset << std::endl;

        eval =H5Dclose(dataset_id);
        if (eval) std::cerr << "Error closing dataset: " << dset << std::endl;

        eval =H5Sclose(dataspace_id);
        if (eval) std::cerr << "Error closing dataspace: " << fname << std::endl;

        /* Write omode velocity to file */
        dataspace_id = H5Screate_simple(2, dims, NULL);

        printf("OVelocity");
        sprintf(dset, "OVelocity");

        dataset_id = H5Dcreate2(file_id, dset, H5T_IEEE_F32BE, dataspace_id,
            H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

        eval =H5Dwrite(dataset_id, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
            &rxvel[0].front());
        if (eval) std::cerr << "Error writing to dataset: " << dset << std::endl;

        eval =H5Dclose(dataset_id);
        if (eval) std::cerr << "Error closing dataset: " << dset << std::endl;

        eval =H5Sclose(dataspace_id);
        if (eval) std::cerr << "Error closing dataspace: " << fname << std::endl;

        /* Write xmode power to file */
        dataspace_id = H5Screate_simple(2, dims, NULL);

        printf("XPower");
        sprintf(dset, "XPower");
        dataset_id = H5Dcreate2(file_id, dset, H5T_IEEE_F32BE, dataspace_id,
            H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

        eval =H5Dwrite(dataset_id, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
            &rxpow[1].front());
        if (eval) std::cerr << "Error writing to dataset: " << dset << std::endl;

        eval =H5Dclose(dataset_id);
        if (eval) std::cerr << "Error closing dataset: " << dset << std::endl;

        eval =H5Sclose(dataspace_id);
        if (eval) std::cerr << "Error closing dataspace: " << fname << std::endl;

        /* Write xmode velocity to file */
        dataspace_id = H5Screate_simple(2, dims, NULL);

        sprintf(dset, "XVelocity");
        dataset_id = H5Dcreate2(file_id, dset, H5T_IEEE_F32BE, dataspace_id,
            H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

        eval =H5Dwrite(dataset_id, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
            &rxvel[1].front());
        if (eval) std::cerr << "Error writing to dataset: " << dset << std::endl;

        eval =H5Dclose(dataset_id);
        if (eval) std::cerr << "Error closing dataset: " << dset << std::endl;

        eval =H5Sclose(dataspace_id);
        if (eval) std::cerr << "Error closing dataspace: " << fname << std::endl;

        /* Create 1-D dataspace for frequency values */
        dims[0] = nsteps;
        dataspace_id = H5Screate_simple(1, dims, NULL);

        attribute_id = H5Acreate2(file_id, "Frequencies(kHz)", H5T_STD_U32BE, dataspace_id,
            H5P_DEFAULT, H5P_DEFAULT);

        eval =H5Awrite(attribute_id, H5T_NATIVE_UINT, &frequencies.front());
        if (eval) std::cerr << "Error writing to dataset: " << dset << std::endl;

        eval =H5Aclose(attribute_id);
        if (eval) std::cerr << "Error closing dataset: " << dset << std::endl;

        eval =H5Sclose(dataspace_id);
        if (eval) std::cerr << "Error closing dataspace: " << fname << std::endl;

        /* Create attributes for first range, last range */
        dims[0]=datalen;
        dataspace_id = H5Screate_simple(1, dims, NULL);

        attribute_id = H5Acreate2(file_id, "Ranges(km)", H5T_STD_U32BE, dataspace_id,
            H5P_DEFAULT, H5P_DEFAULT);

        H5Awrite(attribute_id, H5T_NATIVE_UINT, &ranges.front());

        H5Aclose(attribute_id);

        eval =H5Sclose(dataspace_id);
        if (eval) std::cerr << "Error closing dataspace: " << fname << std::endl;

        eval =H5Fclose(file_id);
        if (eval) std::cerr << "Error closing file: " << fname << std::endl;
    }

    /* Disconnect from sounding server */
    usrpmsg = 'x';
    send(sockfd, &usrpmsg, sizeof(usrpmsg), 0);
    close(sockfd);

    printf("done!\n");

    return 0;
}
