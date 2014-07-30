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
    std::cout << rawtime << std::endl;
    char fname[80];

    //Variables for testing
    std::vector<float> rxdata;
    std::vector<float> test;

    //variables for command line arguments
    float start_freq;
    float step_freq;
    float stop_freq;
    unsigned int nsteps;
    unsigned int npulses;
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
            "Number of pulses in a single frequency sounding")
        ("write","write to file")
    ;
    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, desc), vm);
    po::notify(vm);

    if (vm.count("help")){
        std::cout << desc << "\n";
        return 1;
    }
            
    parms.freq = 1e3*start_freq;
    parms.txrate = 500e3;
    parms.rxrate = 200e3;
    parms.npulses = npulses;
    parms.symboltime = 40e-6;
    parms.pulsetime = 5e-3;
    sprintf(parms.pc_str,"barker13");
    //parms.nsamps_per_pulse = (parms.pulsetime-parms.symboltime)*parms.rxrate;
    //parms.nsamps_per_pulse = round(3*parms.pulsetime*parms.rxrate)/50;
    //parms.nsamps_per_pulse *= 10;
    parms.nsamps_per_pulse = parms.pulsetime*parms.rxrate;
    int datalen;

    printf("\nmsg values\n");
    printf("freq: %04.f kHz\n", parms.freq/1e3);
    printf("txrate: %03.f kHz\n", parms.txrate/1e3);
    printf("rxrate: %03.f kHz\n", parms.rxrate/1e3);
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
    parms.freq = 1e3*start_freq;
    ifreq = 0;
    while(parms.freq < 1e3*stop_freq+1){
        //parms.freq = 1e6 + i*100e3;

        /* Perform the sounding */
        usrpmsg = 's';
        send(sockfd, &usrpmsg, sizeof(usrpmsg), 0);
        send(sockfd, &parms, sizeof(parms), 0);
        rval = recv(sockfd, &status, sizeof(int),0);
        printf("rx status: %i\n", status);

        /* Process the data */
        usrpmsg = 'p';
        send(sockfd, &usrpmsg, sizeof(usrpmsg), 0);
        rval = recv(sockfd, &status, sizeof(int),0);
        printf("process status: %i\n", status);

        /* Get the data */
        usrpmsg = 'd';
        send(sockfd, &usrpmsg, sizeof(usrpmsg), 0);
        recv(sockfd, &datalen, sizeof(datalen), 0);
        rxdata.resize(datalen);
        test.resize(datalen,100);
        recv(sockfd, &rxdata.front(), rxdata.size()*sizeof(float), 0);

        /* Write the data to hdf5 file */
        dims[0] = datalen;
        dims[1] = 1;
        if (vm.count("write")){
            dataspace_id = H5Screate_simple(1, dims, NULL);
        }

        sprintf(dset, "%05.f",parms.freq/1e3);
        if (vm.count("write")){
            dataset_id = H5Dcreate2(file_id, dset, H5T_IEEE_F32BE, dataspace_id,
                H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        }

        if (vm.count("write")){
            eval =H5Dwrite(dataset_id, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                &rxdata.front());
            if (eval) std::cerr << "Error writing to dataset: " << dset << std::endl;
        }

        if (vm.count("write")){
            eval =H5Dclose(dataset_id);
            if (eval) std::cerr << "Error closing dataset: " << dset << std::endl;
        }

        printf("freq %i of %i (%.0f kHz)\n", ifreq+1, nsteps, parms.freq/1e3);
        //for (size_t i=0; i<rxdata.size(); i++){
        //    printf("%lu: %.1f\n",i,30+10*log10(rxdata[i]));
        //}
        parms.freq += 1e3*step_freq;
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

    printf("done!\n");

    return 0;
}
