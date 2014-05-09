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

#include "global_variables.h"
#include "hdf5.h"

int main(){
    /* Socket and status variables*/
    int sockfd=0, n=0;
    int rval =0, status=0;
    struct sockaddr_in serv_addr;
    struct soundingParms parms;
    char usrpmsg;
    //std::string dset;
    char dset[80];

    /* hdf5 variables */
    hid_t file_id, dataspace_id, dataset_id;
    herr_t eval;
    hsize_t dims[2];
    std::string filename = "/tmp/ionosonde_data.h5";
    file_id = H5Fcreate(filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    //file_id = H5Fopen(filename, H5F_ACC_RDWR, H5P_DEFAULT);


    //parms.rxfile = "/home/radar/UAF_USRP/uhd/sounder/rx.dat";
    parms.freq = 12e6;
    parms.txrate = 2000e3;
    parms.rxrate = 200e3;
    parms.npulses = 256;
    parms.symboltime = 200e-6;
    parms.pulsetime = 10e-3;
    //parms.nsamps_per_pulse = (parms.pulsetime-parms.symboltime)*parms.rxrate;
    parms.nsamps_per_pulse = 4*parms.pulsetime*parms.rxrate/5;
    int datalen;
    std::vector<float> rxdata;
    std::vector<float> test;

    printf("\nmsg values\n");
    printf("freq: %f\n", parms.freq);
    printf("txrate: %f\n", parms.txrate);
    printf("rxrate: %f\n", parms.rxrate);
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

    for (int i=0; i<10; i++){
        parms.freq = 1e6 + i*100e3;

        usrpmsg = 's';
        send(sockfd, &usrpmsg, sizeof(usrpmsg), 0);
        send(sockfd, &parms, sizeof(parms), 0);
        rval = recv(sockfd, &status, sizeof(int),0);
        printf("rx status: %i\n", status);

        usrpmsg = 'p';
        send(sockfd, &usrpmsg, sizeof(usrpmsg), 0);
        rval = recv(sockfd, &status, sizeof(int),0);
        printf("process status: %i\n", status);

        usrpmsg = 'd';
        send(sockfd, &usrpmsg, sizeof(usrpmsg), 0);
        recv(sockfd, &datalen, sizeof(datalen), 0);
        rxdata.resize(datalen);
        test.resize(datalen,100);

            

        recv(sockfd, &rxdata.front(), rxdata.size()*sizeof(float), 0);

        dims[0] = datalen;
        dims[1] = 1;
        dataspace_id = H5Screate_simple(1, dims, NULL);
        //printf("dataspace_id: %i\n",dataspace_id);
        //dset = std::to_string(parms.freq);
        sprintf(dset, "%05.f",parms.freq/1e3);
        std::cout << "dset string: " << dset << std::endl;

        dataset_id = H5Dcreate2(file_id, dset, H5T_IEEE_F32BE, dataspace_id,
            H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        //printf("dataset_id: %i\n",dataset_id);
        eval = -1;
        eval =H5Dwrite(dataset_id, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
            &rxdata.front());
        //eval =H5Dread(dataset_id, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
        //    &rxdata.front());
        printf("eval: %i\n",eval);
        eval =H5Dclose(dataset_id);
        printf("eval: %i\n",eval);

        //for (size_t i=0; i<rxdata.size(); i++){
        //    printf("%lu: %.1f\n",i,30+10*log10(rxdata[i]));
        //}


    }
    eval =H5Sclose(dataspace_id);
    printf("eval: %i\n",eval);
    eval =H5Fclose(file_id);
    printf("eval: %i\n",eval);
    usrpmsg = 'x';
    send(sockfd, &usrpmsg, sizeof(usrpmsg), 0);

    printf("done!\n");

    return 0;
}
