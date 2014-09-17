/*Port number that client and server use*/
#define HOST_PORT 45001
/*Definitions of commands that client gives to server*/
#define SEND 's'
#define LISTEN 'l'
#define EXIT 'x'
#define GET_DATA 'd'
#define PROCESS 'p'

#define MAX_VELOCITY 100 //The maximum unambiguous Doppler velocity. Low number improves computational efficiency
#define OSR 1 //Factor by which each range gate is oversampled.
#define TX_RATE 1000e3 //Sample rate "over-the-wire" between host computer and USRP
#define RX_RATE 1000e3 //Sample rate "over-the-wire" between host computer and USRP
#define XOVER_FREQ 30e6 //Frequency at which USRP switches from antenna port A to B
#define MAX_CLRSEARCH_SPAN 1000
#define MIN_CLRSEARCH_SPAN 1000

/* Pulse codes known to both client and server */
#define BARKER_13 {0.,1.,1.,1.,1.,1.,-1.,-1.,1.,1.,-1.,1.,-1.,1.,0.}
#define GOLAY_16_0 {1.,1.,1.,-1.,1.,1.,-1.,1., 1.,1.,1.,-1.,-1.,-1.,1.,-1.};
#define GOLAY_16_1 {1.,1.,1.,-1.,1.,1.,-1.,1., -1.,-1.,-1.,1.,1.,1.,-1.,1.};
#define GOLAY_10_0 {1.,1.,-1.,1.,-1.,1.,-1.,-1.,1.,1.};
#define GOLAY_10_1 {1.,1.,-1.,1.,1.,1.,1.,1.,-1.,-1.};
#define GOLAY_8_0 {1.,1.,1.,-1.,1.,1.,-1.,1.};
#define GOLAY_8_1 {1.,1.,1.,-1.,-1.,-1.,1.,-1.};
#define GOLAY_4_0 {1.,1.,1.,-1.};
#define GOLAY_4_1 {1.,1.,-1.,1.};
#define RECT {1.}

/*Data structures shared between client and server*/
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

struct soundingParms2{
    size_t freq_khz;
    size_t num_pulses;
    int first_range_km;
    int last_range_km;
    size_t num_rx_samples;
    size_t over_sample_rate;
    float range_res_km;
};

struct periodogramParms{
    size_t start_freq_khz;
    size_t end_freq_khz; //End frequency, NOT inclusive
    size_t bandwidth_khz; //Bandwidth of each spectral bin. There's a bin every 1 khz, but there is spectral filtering!!
};
