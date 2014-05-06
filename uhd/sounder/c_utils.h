#ifdef __cplusplus 
    extern "C" {
#endif

int tcpsocket(int port);
int recv_data(int fd, void *buf, size_t buflen);

#ifdef __cplusplus
    }
#endif
