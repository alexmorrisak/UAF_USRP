#define SEND 's'
#define LISTEN 'l'
#define EXIT 'x'
#define GET_DATA 'd'
#define PROCESS 'p'
#define MAX_VELOCITY 100 //The maximum unambiguous velocity

#define BARKER_13 {0.,1.,1.,1.,1.,1.,-1.,-1.,1.,1.,-1.,1.,-1.,1.,0.}
#define GOLAY_8_0 {1.,1.,1.,-1.,1.,1.,-1.,1.};
#define GOLAY_8_1 {1.,1.,1.,-1.,-1.,-1.,1.,-1.};
#define GOLAY_4_0 {1.,1.,1.,-1.};
#define GOLAY_4_1 {1.,1.,-1.,1.};
#define RECT {1.}

struct soundingParms{
    size_t freq;
    size_t txrate_khz;
    size_t rxrate_khz;
    size_t npulses;
    size_t nsamps_per_pulse;
    size_t symboltime_usec;
    size_t pulsetime_usec;
    char pc_str[80];
};

struct periodogramParms{
    size_t start_freq_khz;
    size_t end_freq_khz; //End frequency, NOT inclusive
    size_t bandwidth_khz; //Bandwidth of each spectral bin. There's a bin every 1 khz, but there is spectral filtering!!
};
