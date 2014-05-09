#include <sys/socket.h>
#include <sys/types.h>
#include <iostream>
#include <string>
#include <fstream>
#include <string.h>
#include <vector>
#include <complex>

#include <netinet/in.h> //for sockaddr_in
#include <arpa/inet.h> //for inet_addr()
#include <unistd.h> //for close()

#include "global_variables.h"

int main(){
    int sockfd=0, n=0;
    int rval =0, status=0;
    struct sockaddr_in serv_addr;
    struct soundingParms parms;
    char usrpmsg;
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

    for (int i=0; i<1; i++){
        parms.freq = 10e6 + i*100e3;

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
        recv(sockfd, &rxdata.front(), rxdata.size()*sizeof(float), 0);
        for (size_t i=0; i<rxdata.size(); i++){
            printf("%i: %.1f\n",i,30+10*log10(rxdata[i]));
        }


    }
    usrpmsg = 'x';
    send(sockfd, &usrpmsg, sizeof(usrpmsg), 0);

    printf("done!\n");

    return 0;
}
