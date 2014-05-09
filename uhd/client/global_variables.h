#define SEND 's'
#define EXIT 'x'
#define GET_DATA 'd'
#define PROCESS 'p'
#define MAX_VELOCITY 100 //The maximum unambiguous velocity

#define BARKER_13 {1.,1.,1.,1.,1.,-1.,-1.,1.,1.,-1.,1.,-1.,1.}
#define RECT {1.}

struct soundingParms{
    float freq;
    float txrate;
    float rxrate;
    unsigned int npulses;
    unsigned int nsamps_per_pulse;
    float symboltime;
    float pulsetime;
    char pc_str[80];
};
