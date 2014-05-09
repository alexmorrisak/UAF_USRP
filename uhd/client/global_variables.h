#define SEND 's'
#define EXIT 'x'
#define GET_DATA 'd'
#define PROCESS 'p'
#define MAX_VELOCITY 100 //The maximum unambiguous velocity

struct soundingParms{
    float freq;
    float txrate;
    float rxrate;
    unsigned int npulses;
    unsigned int nsamps_per_pulse;
    float symboltime;
    float pulsetime;
};
