#include <pthread.h>
#include <stdio.h>
#include <sys/types.h>
#include <sys/socket.h>
#include <sys/un.h>
#include <netinet/in.h>
#include <netinet/tcp.h>
#include <netdb.h>
#include <unistd.h>
#include <stdlib.h>
#include <string.h>

#ifdef __cplusplus
	extern "C" {
#endif

//#define MSG_NOSIGNAL 0 
int opentcpsock(char *hostip, int port){
	//DECLARE VARIABLES FOR IP CONNECTIONS
	char	datacode,ipaddr[16];
	int		sock,data,temp;
	struct	sockaddr_in	server;
	struct	hostent		*hp;//, *gethostbyname();
	float	ftime;
	int		buffer2,channel2,sample2i,option;//,optionlen;
	socklen_t	optionlen;
	int		do_scan_rx[4];
	struct 	protent*	protocol_info;

	//hostip=HOST;
	//port=DEFAULT_PORT;

	//SET UP IP CONNECTION
	sock=socket(AF_INET, SOCK_STREAM, 0);
	if( (sock < 0) ) {
		perror("opening stream socket");
		exit(1);
	}
	server.sin_family=AF_INET;
	hp=gethostbyname(hostip);
	if( hp == 0 ){
		fprintf(stderr, "unknown host");
		exit(2);
	}
	memcpy(&server.sin_addr, hp->h_addr, hp->h_length);
	server.sin_port=htons(port);
	temp=connect(sock, (struct sockaddr *)&server, (socklen_t) sizeof(server));
	if( temp < 0){
		perror("connecting stream socket");
		sock=-1;
	}


	//protocol_info=getprotobyname("tcp");
	//printf("protocol name = %s\n",protocol_info->p_name);
	option=TCP_NODELAY;
	optionlen=4;
	temp=setsockopt(sock,6,TCP_NODELAY,&option,optionlen);
	temp=getsockopt(sock,6,TCP_NODELAY,&option,&optionlen);
//	printf("temp=%d  optionlen=%d option=%d\n",temp,optionlen,option);
	optionlen=4;
	option=32768;
	temp=setsockopt(sock,SOL_SOCKET,SO_SNDBUF,&option,optionlen);
	temp=getsockopt(sock,SOL_SOCKET,SO_SNDBUF,&option,&optionlen);
//	printf("temp=%d  optionlen=%d option=%d\n",temp,optionlen,option);
	optionlen=4;
	option=32768;
	temp=setsockopt(sock,SOL_SOCKET,SO_RCVBUF,&option,optionlen);
	temp=getsockopt(sock,SOL_SOCKET,SO_RCVBUF,&option,&optionlen);
//	printf("temp=%d  optionlen=%d option=%d\n",temp,optionlen,option);

    return sock;
}

int server_unixsocket(char *hostip,int port){

	char	path[256];
	int	sock,len,temp;
	struct	sockaddr_un	saun;
	int	option;//,optionlen;
	socklen_t	optionlen;


    if ((sock = socket(AF_UNIX, SOCK_STREAM, 0)) < 0) {
        perror("server: socket");
        exit(1);
    }

    saun.sun_family = AF_UNIX;
    sprintf(path,"%s",hostip);

    strcpy(saun.sun_path, path);
    fprintf(stdout,"Unix Path: %s\n",path);
    unlink(path);
    len = sizeof(saun.sun_family) + strlen(saun.sun_path);

    if (bind(sock, (struct sockaddr *)&saun, len) < 0) {
        perror("server: bind");
        exit(1);
    }
	optionlen=4;
	option=32768;
	temp=setsockopt(sock,SOL_SOCKET,SO_SNDBUF,&option,optionlen);
	temp=getsockopt(sock,SOL_SOCKET,SO_SNDBUF,&option,&optionlen);
	optionlen=4;
	option=32768;
	temp=setsockopt(sock,SOL_SOCKET,SO_RCVBUF,&option,optionlen);
	temp=getsockopt(sock,SOL_SOCKET,SO_RCVBUF,&option,&optionlen);

	/* return the socket */
	return sock;
}

int openunixsock(char *hostip, int port){
	//DECLARE VARIABLES FOR IP CONNECTIONS
	char	path[256];
	int	sock,len,temp;
	struct	sockaddr_un	saun;
	int	option;//,optionlen;
	socklen_t	optionlen;

        if ((sock = socket(AF_UNIX, SOCK_STREAM, 0)) < 0) {
	  perror("opening unix domain stream socket");
          exit(1);
        }
	saun.sun_family=AF_UNIX;
	sprintf(path,"%s",hostip);
        fprintf(stderr,"Sock Path: %s\n",path);
        strcpy(saun.sun_path, path);
        //unlink(path);
        len = sizeof(saun.sun_family) + strlen(saun.sun_path);
        if (connect(sock,(struct sockaddr *) &saun,(socklen_t) len) < 0) {
          perror("client: connect");
          sock=-1;
        }

	optionlen=4;
	option=32768;
	temp=setsockopt(sock,SOL_SOCKET,SO_SNDBUF,&option,optionlen);
	temp=getsockopt(sock,SOL_SOCKET,SO_SNDBUF,&option,&optionlen);
	optionlen=4;
	option=32768;
	temp=setsockopt(sock,SOL_SOCKET,SO_RCVBUF,&option,optionlen);
	temp=getsockopt(sock,SOL_SOCKET,SO_RCVBUF,&option,&optionlen);

    return sock;
}
int send_data(int fd,void  *buf,size_t buflen)
{
     int cc=0,total=0;
     while (buflen > 0) {
          cc = send(fd, buf, buflen, MSG_NOSIGNAL);
          if (cc == -1) {
               return cc;
          }
          if (cc == 0) {
            return -1;
          }

          buf += cc;
          total += cc;
          buflen -= cc;
     }
     return total;
}

int recv_data(int fd,void *buf,size_t buflen)
{
     int cc=0,total=0;
     while (buflen > 0) {
          cc = recv(fd, buf, buflen, MSG_NOSIGNAL);
          if (cc == -1) {
               return cc;
          }
          if (cc == 0) {
            return -1;
          }
          buf += cc;
          total += cc;
          buflen -= cc;
     }
     return total;

}


#ifdef __cplusplus
	}
#endif
