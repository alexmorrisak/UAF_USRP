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
    std::vector<std::vector<std::complex<double> > > doppler_grid(slowdim);
    for (int i=0; i<slowdim; i++)
        doppler_grid[i].resize(fastdim,0);
    std::vector<double> tmp_pwr(slowdim,0);
    std::vector<double> pwr(slowdim,0);

    if(in!=NULL){free(in);in=NULL;}
    in = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*slowdim);
    if(out!=NULL){free(out);out=NULL;}
    out = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*slowdim);

    plan = fftw_plan_dft_1d(slowdim, in, out, FFTW_FORWARD, FFTW_ESTIMATE);

    // Calculate the FFT for each range bin
    for (int irange=0; irange<fastdim; irange++){
        for (int i=0; i<slowdim;i++){
            in[i][0] = indata[i][irange].real();
            in[i][1] = indata[i][irange].imag();
        }
        fftw_execute(plan);
        for (int i=0; i<slowdim; i++){
            doppler_grid[i][irange] = std::complex<double>(out[i][0], out[i][1]);
        }
    }
    // Find and remove any narrow-band interference
    for (int i=0; i<slowdim; i++){
        std::complex<double> DC(0,0);
        for (int j=20; j<fastdim; j++)
            DC += doppler_grid[i][j];
        DC /= (fastdim-20);
        //printf("slowdim %i, DC: %.1f, %.1f\n", i, std::abs(DC), 180/M_PI*std::arg(DC));
        for (int j=0; j<fastdim; j++)
             doppler_grid[i][j] -= DC;
    }

    // Find the max spectral component in each range bin
    for (int irange=0; irange<fastdim; irange++){
        int maxinx = 0;
        for (int i=1; i<slowdim; i++){
            if (std::abs(doppler_grid[i][irange]) >
                std::abs(doppler_grid[maxinx][irange]))
                    maxinx = i;
            //printf("Mininx: %i / %i\n", maxinx, slowdim);
        }
        outpow[irange] = 
            std::abs(doppler_grid[maxinx][irange])*std::abs(doppler_grid[maxinx][irange]);
            //(out[maxinx][0]*out[maxinx][0] + out[maxinx][1]*out[maxinx][1]);
        for (int j=0; j<slowdim; j++){
            //printf("outpow %i: %.1f @ %.1f\t",irange, std::abs(doppler_grid[j][irange]), 180/M_PI*std::arg(doppler_grid[j][irange]));
        }
        //printf("\n");
        if (maxinx > slowdim/2)
            maxinx = maxinx - slowdim;
        outvel[irange] = (float) maxinx;
    }

    fftw_destroy_plan(plan);
    if(in!=NULL){free(in);in=NULL;}
    if(out!=NULL){free(out);out=NULL;}

    return 0;
}

