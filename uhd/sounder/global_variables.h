struct soundingParms{
    char txfile[80];
    char rxfile[80];
    unsigned int nrxsamples;
    float freq;
    unsigned int txrate;
    unsigned int rxrate;
    unsigned int npulses;
    unsigned int nsamps_per_pulse;
};
