// recv_to_file.cpp

#include <iostream>
#include <string>
#include <fstream>
#include <complex>
#include <csignal>
#include <cmath>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <fftw3.h>
#include <vector>
#include <stdlib.h>

/***********************************************************************
 * doppler_process function
 **********************************************************************/
int doppler_process(
    std::vector<std::complex<float> *> indata,
    float* outpow,
    float* outvel,
    int slowdim,
    int fastdim
){
    fftw_complex *in =NULL, *out=NULL; 
    fftw_plan plan;
    std::vector<double> tmp_pwr(slowdim,0);
    std::vector<double> pwr(slowdim,0);

    if(in!=NULL){free(in);in=NULL;}
    in = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*slowdim);
    if(out!=NULL){free(out);out=NULL;}
    out = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*slowdim);

    plan = fftw_plan_dft_1d(slowdim, in, out, FFTW_FORWARD, FFTW_ESTIMATE);

    for (int irange=0; irange<fastdim; irange++){
        for (int i=0; i<slowdim;i++){
            in[i][0] = indata[i][irange].real();
            in[i][1] = indata[i][irange].imag();
        }
        fftw_execute(plan);

        int mininx = 0;
        for (int i=1; i<slowdim; i++){
            if ((std::abs(out[i][0])+std::abs(out[i][1])) > 
                (std::abs(out[mininx][0]) + std::abs(out[mininx][1])))
                    mininx = i;
            //printf("Mininx: %i / %i\n", mininx, slowdim);
        }
        outpow[irange] = 
            (out[mininx][0]*out[mininx][0] + out[mininx][1]*out[mininx][1]);
        if (mininx > slowdim/2)
            mininx = mininx - slowdim;
        outvel[irange] = (float) mininx;
    }

    fftw_destroy_plan(plan);
    if(in!=NULL){free(in);in=NULL;}
    if(out!=NULL){free(out);out=NULL;}

    return 0;
}

