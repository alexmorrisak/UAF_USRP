// lp_filter.cpp

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

extern int verbose;
//int verbose=2;

/***********************************************************************
 * lp_filter function
 **********************************************************************/
int lp_filter(
    std::vector<std::complex<float> *> indata,
    std::vector<std::complex<float> *> outdata,
    int slowdim,
    int fastdim,
    float samprate,
    float bw,
    int decimrate
){
    int debug = 0;

    //Create a vector of filter taps for low-pass filtering.  
    //This is a simple hamming-windowed sinc function
    float ntaps = 4*samprate/bw;
    std::vector<std::complex<float> > filter_taps(ntaps);
    for (int i=0;i<ntaps;i++){
        double x=4*(2*M_PI*((float)i/ntaps)-M_PI);
        filter_taps[i] = std::complex<float>(
            bw/samprate*(0.54+0.46*cos((2*M_PI*((float)(i-ntaps/2)+0.5))/ntaps))*sin(x)/x,
            0);
    }
    filter_taps[ntaps/2] = std::complex<float>(bw/samprate,0);

    if(true){
        for (int i=0; i<ntaps; i++){
            printf("%i: (%.4f,%.4f)\n", i, filter_taps[i].real(), filter_taps[i].imag());
        }
    }

    //populate the zero-padded temporary vector with input samples
    //This copies the input samples to a new memory location.  For greater efficiency, the samples
    //should be placed in the center of a "large" zero-valued array right off the bat.
    for (int ipulse=0;ipulse<slowdim; ipulse++){
        std::complex<float> DC(0,0);
        std::vector<std::complex<float> > tempvec(fastdim+ntaps,0);
        for (int i=0; i<fastdim; i++){
            tempvec[ntaps/2+i] = indata[ipulse][i];
        }
        //Print the zero-padded vector of input samples
        if (debug){
            for (int i=0; i<(fastdim+ntaps); i++){
                printf("in: %i (%.2f, %.2f)\n", i, tempvec[i].real(), tempvec[i].imag());
            }
        }
        //perform the convolution
        std::complex<float> temp(0,0);
        for (int isamp =0; isamp<fastdim; isamp+=decimrate){
            temp = std::complex<float>(0,0);
            for (int i=0; i<ntaps; i++){
                temp += filter_taps[i]*tempvec[isamp+i];
            }
            outdata[ipulse][isamp/decimrate] = temp;
        }
    }
    return 0;
}

//int main(){
//    std::vector<std::vector<std::complex<float> > > invecs(10);
//    std::vector<std::vector<std::complex<float> > > outvecs(10);
//    std::vector<std::complex<float> *> in;
//    std::vector<std::complex<float> *> out;
//    int nsamps = 16;
//    float sr = 100;
//    float bw = 20;
//    int dmrate = 1;
//
//    in.resize(10);
//    for (int i=0; i<10; i++){
//        invecs[i].resize(nsamps,10.);
//    	in[i] = &invecs[i].front();
//    }
//
//    out.resize(10);
//    for (int i=0; i<10; i++){
//        outvecs[i].resize(nsamps,10.);
//    	out[i] = &outvecs[i].front();
//    }
//    //outvecs[0].resize(nsamps);
//    //out.resize(nsamps);
//    //out[0] = &outvecs[0].front();
//
//    int rval = lp_filter(
//        in,
//        out,
//        1,
//        nsamps,
//        sr,
//        bw,
//        dmrate);
//
//    printf("ouputs\n");
//    for (int i=0; i<nsamps; i++){
//        printf("%i (%.2f, %.2f)\n", i, out[0][i].real(), out[0][i].imag());
//    }
//}
