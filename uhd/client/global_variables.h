#define SEND 's'
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
    float freq;
    size_t txrate_khz;
    size_t rxrate_khz;
    size_t npulses;
    size_t nsamps_per_pulse;
    size_t symboltime_usec;
    size_t pulsetime_usec;
    char pc_str[80];
};
