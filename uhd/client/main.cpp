#include <sys/socket.h>
#include <sys/types.h>
#include <iostream>
#include <string>
#include <fstream>
#include <string.h>

#include <netinet/in.h> //for sockaddr_in
#include <arpa/inet.h> //for inet_addr()
#include <unistd.h> //for close()

#include "global_variables.h"

int main(){
    int sockfd=0, n=0;
    int sbuf = 69;
    struct sockaddr_in serv_addr;
    struct soundingParms parms;
    strncpy(parms.txfile,"/home/radar/UAF_USRP/process/tx.dat",80);
    strncpy(parms.rxfile,"/home/radar/UAF_USRP/uhd/sounder/rx.dat",80);
    //parms.rxfile = "/home/radar/UAF_USRP/uhd/sounder/rx.dat";
    parms.freq = 12e6;
    parms.txrate = 200e3;
    parms.rxrate = 200e3;
    parms.npulses = 1024;
    parms.symboltime = 1000;
    parms.pulsetime = 5;
    parms.nsamps_per_pulse = (1e-3*parms.pulsetime-1e-6*parms.symboltime)*parms.rxrate;
    parms.nsamps_per_pulse = (0.90*1e-3*parms.pulsetime)*parms.rxrate;

    printf("\nmsg values\n");
    printf("txfile: %s\n", parms.txfile);
    printf("rxfile: %s\n", parms.rxfile);
    printf("freq: %f\n", parms.freq);
    printf("txrate: %i\n", parms.txrate);
    printf("rxrate: %i\n", parms.rxrate);
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

    char usrpmsg = 's';
    send(sockfd, &usrpmsg, sizeof(usrpmsg), 0);
    send(sockfd, &parms, sizeof(parms), 0);
    usrpmsg = 'x';
    send(sockfd, &usrpmsg, sizeof(usrpmsg), 0);

    printf("done!\n");

    return 0;
}
