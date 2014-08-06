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

#include <boost/program_options.hpp>
#include <boost/format.hpp>

#include "global_variables.h"
#include "hdf5.h"

namespace po = boost::program_options;

int main(int argc, char *argv[]){
    /* Socket and status variables*/
    int sockfd=0, n=0;
    int rval =0, status=0;
    struct sockaddr_in serv_addr;
    struct soundingParms parms;
    struct periodogramParms spect_parms;
    char usrpmsg;

    /* hdf5 and file-writing variables */
    char dset[80];
    hid_t file_id, dataspace_id, dataset_id;
    herr_t eval = 0;
    hsize_t dims[2];
    time_t rawtime;
    struct tm * timeinfo;
    time(&rawtime);
    timeinfo = localtime(&rawtime);
    //std::cout << rawtime << std::endl;
    char fname[80];

    //Variables for testing
    std::vector<float> rxdata[2];
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
    unsigned int first_range,last_range;
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
        ("last-range", po::value<unsigned int>(&last_range)->default_value(750),
            "Maximum unambiguous range in km")
        ("first-range", po::value<unsigned int>(&first_range)->default_value(50),
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

    parms.txrate_khz = 500;
    parms.rxrate_khz = 500;
    parms.npulses = npulses;

    size_t temp_symboltime_usec = resolution / 1.5e-1;
    size_t temp_dmrate = (size_t) ceil(temp_symboltime_usec * parms.rxrate_khz / (OSR*1000));
    //if (temp_dmrate%2 == 1) temp_dmrate -= 1;
    while (temp_dmrate%OSR != 0) temp_dmrate -= 1;
    parms.symboltime_usec = OSR*1000*temp_dmrate/parms.rxrate_khz;
    //parms.symboltime_usec = 30;

    size_t temp_pfactor = (size_t) ceil(2*last_range / 3.0e-1 / parms.symboltime_usec);
    parms.pulsetime_usec = temp_pfactor * parms.symboltime_usec;

    std::cout << "Using pulse time: " << parms.pulsetime_usec << std::endl;
    std::cout << "symbol time: " << parms.symboltime_usec << std::endl;
    std::cout << "Using range resolution: " << 1.5e-1*parms.symboltime_usec << std::endl;

    size_t max_txtime_usec = (size_t) floor(2*first_range / 3.0e-1);
    size_t max_code_length = (size_t) floor(max_txtime_usec / parms.symboltime_usec);
    if (max_code_length < 4){
        sprintf(parms.pc_str,"rect");
    }
    else if (max_code_length < 8){
        sprintf(parms.pc_str,"golay4");
    }
    else if (max_code_length < 10){
        sprintf(parms.pc_str,"golay8");
    }
    //else if (max_code_length < 16){
    //    sprintf(parms.pc_str,"golay10");
    //}
    else{
        sprintf(parms.pc_str,"golay16");
    }
    std::cout << "Using pcode: " << parms.pc_str << std::endl;
        

    //parms.pulsetime_usec = 5000;

    //sprintf(parms.pc_str,"golay8");
    parms.nsamps_per_pulse = (size_t) (1e-6*parms.pulsetime_usec*RX_RATE);
    int datalen;

    printf("\nmsg values\n");
    printf("freq: %04.f kHz\n", nominal_freq);
    printf("txrate: %03.f kHz\n", parms.txrate_khz);
    printf("rxrate: %03.f kHz\n", parms.rxrate_khz);
    printf("nsamps per pulse: %i\n", parms.nsamps_per_pulse);

    if ((sockfd = socket(AF_INET, SOCK_STREAM, 0)) < 0){
        printf("Error in creating socket\n");
        return 1;
    }

    serv_addr.sin_family = AF_INET;
    serv_addr.sin_addr.s_addr = inet_addr("127.0.0.1");
    serv_addr.sin_port = htons(45001);

    if (connect(sockfd, (struct sockaddr *)&serv_addr, sizeof(serv_addr)) < 0){
        printf("bind failed\n");
        return 1;
    }
    printf ("connected!\n");

    strftime(fname, 80, "ionogram.%Y%m%d.%H%M.h5",timeinfo);
    if (vm.count("write"))
        file_id = H5Fcreate(fname, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

    step_freq = (stop_freq+1 - start_freq) / nsteps;
    nominal_freq = start_freq;
    ifreq = 0;
    while(nominal_freq < stop_freq){
        //parms.freq = 1e6 + i*100e3;

        /* Get the spectrum so that we can get the quietest frequency*/
        usrpmsg = 'l';
        send(sockfd, &usrpmsg, sizeof(usrpmsg), 0);
        spectrum.resize(200);
        spect_parms.start_freq_khz = nominal_freq-100;
        spect_parms.end_freq_khz = nominal_freq+100;
        spect_parms.bandwidth_khz = 10;
        send(sockfd, &spect_parms, sizeof(parms), 0);
        recv(sockfd, &spectrum.front(), 200*sizeof(float), 0);
        min_inx = 0;
        for (int i=1; i<spectrum.size(); i++){
            if (spectrum[i] < spectrum[min_inx]) min_inx = i;
        }
        parms.freq = min_inx + spect_parms.start_freq_khz;
        rval = recv(sockfd, &status, sizeof(int),0);
        //printf("rx status: %i\n", status);

        /* Perform the sounding */
        usrpmsg = 's';
        send(sockfd, &usrpmsg, sizeof(usrpmsg), 0);
        send(sockfd, &parms, sizeof(parms), 0);
        rval = recv(sockfd, &status, sizeof(int),0);
        //printf("rx status: %i\n", status);

        /* Process the data */
        usrpmsg = 'p';
        send(sockfd, &usrpmsg, sizeof(usrpmsg), 0);
        rval = recv(sockfd, &status, sizeof(int),0);
        //printf("process status: %i\n", status);

        /* Get the data */
        usrpmsg = 'd';
        send(sockfd, &usrpmsg, sizeof(usrpmsg), 0);
        recv(sockfd, &datalen, sizeof(datalen), 0);
        rxdata[0].resize(datalen);
        rxdata[1].resize(datalen);
        test.resize(datalen,100);
        recv(sockfd, &rxdata[0].front(), rxdata[0].size()*sizeof(float), 0);
        recv(sockfd, &rxdata[1].front(), rxdata[1].size()*sizeof(float), 0);

        /* Write the data to hdf5 file */
        dims[0] = datalen;
        dims[1] = 1;
        if (vm.count("write")){
            dataspace_id = H5Screate_simple(1, dims, NULL);
        }

        //sprintf(dset, "%05.f",parms.freq);
        sprintf(dset, "omode_%05i",nominal_freq);
        if (vm.count("write")){
            dataset_id = H5Dcreate2(file_id, dset, H5T_IEEE_F32BE, dataspace_id,
                H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        }

        if (vm.count("write")){
            eval =H5Dwrite(dataset_id, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                &rxdata[0].front());
            if (eval) std::cerr << "Error writing to dataset: " << dset << std::endl;
        }

        if (vm.count("write")){
            eval =H5Dclose(dataset_id);
            if (eval) std::cerr << "Error closing dataset: " << dset << std::endl;
        }

        sprintf(dset, "xmode_%05i",nominal_freq);
        if (vm.count("write")){
            dataset_id = H5Dcreate2(file_id, dset, H5T_IEEE_F32BE, dataspace_id,
                H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        }

        if (vm.count("write")){
            eval =H5Dwrite(dataset_id, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                &rxdata[1].front());
            if (eval) std::cerr << "Error writing to dataset: " << dset << std::endl;
        }

        if (vm.count("write")){
            eval =H5Dclose(dataset_id);
            if (eval) std::cerr << "Error closing dataset: " << dset << std::endl;
        }

        printf("Sounding %i of %i complete (%i kHz)\n", ifreq+1, nsteps, parms.freq);
        //for (size_t i=0; i<rxdata.size(); i++){
        //    printf("%lu: %.1f\n",i,30+10*log10(rxdata[i]));
        //}
        nominal_freq += step_freq;
	ifreq++;
    }

    if (vm.count("write")){
        eval =H5Sclose(dataspace_id);
        if (eval) std::cerr << "Error closing dataspace: " << fname << std::endl;
    }

    if (vm.count("write")){
        eval =H5Fclose(file_id);
        if (eval) std::cerr << "Error closing file: " << fname << std::endl;
    }

    /* Close the sounding server */
    usrpmsg = 'x';
    send(sockfd, &usrpmsg, sizeof(usrpmsg), 0);
    close(sockfd);

    printf("done!\n");

    return 0;
}
