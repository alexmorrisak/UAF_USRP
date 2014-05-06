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
    std::cout << "entering lp filter\n";
    float ntaps = 4*samprate/bw;
    std::cout << "ntaps: " << ntaps << std::endl;
    std::vector<std::complex<float> > filter_taps(ntaps);
    for (int i=0;i<ntaps;i++){
        double x=2*M_PI*((float)i/ntaps)-M_PI;
        filter_taps[i] = std::complex<float>(
            (0.54-0.46*cos((2*M_PI*((float)(i)+0.5))/ntaps))*sin(2*x)/(2*x)/ntaps,
            0);
        //filter_taps[i] = std::complex<float>(1./ntaps,0);
    }
    filter_taps[ntaps/2] = std::complex<float>(1./ntaps,0);

    //for (int i=0; i<ntaps; i++){
    //    printf("%i: (%.2f,%.2f)\n", i, 1e2*filter_taps[i].real(), 1e2*filter_taps[i].imag());
    //}

    std::cout << "performing convolution\n";
    //populate the zero-padded temporary vector
    for (int ipulse=0;ipulse<slowdim; ipulse++){
        std::vector<std::complex<float> > tempvec(fastdim+ntaps,0);
        for (int i=0; i<fastdim; i++){
            tempvec[ntaps/2+i] = indata[ipulse][i];
            //printf("in: %i (%.2f, %.2f)\n", i, tempvec[i].real(), tempvec[i].imag());
        }
        //for (int i=0; i<fastdim+ntaps; i++){
        //    printf("in: %i (%.2f, %.2f)\n", i, tempvec[i].real(), tempvec[i].imag());
        //}

        //std::cout << "ipulse: " << ipulse << std::endl;
        //perform the convolution
        for (int isamp =0; isamp<fastdim; isamp+=decimrate){
            //std::cout << "isamp: " << isamp << std::endl;
            std::complex<float> temp(0,0);
            for (int i=0; i<ntaps; i++){
                temp += filter_taps[i]*tempvec[isamp+i];
                //printf("%i <-> %i\n",i,isamp+i);
            }
            outdata[ipulse][isamp/decimrate] = temp;
            //printf("out %i,%i: %.2f\n",ipulse,isamp,std::abs(outdata[ipulse][isamp]));
        }
    }




    return 0;
}

//int main(){
//    std::vector<std::vector<std::complex<float> > > invecs(1);
//    std::vector<std::vector<std::complex<float> > > outvecs(1);
//    std::vector<std::complex<float> *> in;
//    std::vector<std::complex<float> *> out;
//    int nsamps = 16;
//    float sr = 100;
//    float bw = 20;
//    int dmrate = 1;
//
//    invecs[0].resize(nsamps,10.);
//    in.resize(nsamps);
//    in[0] = &invecs[0].front();
//
//    outvecs[0].resize(nsamps);
//    out.resize(nsamps);
//    out[0] = &outvecs[0].front();
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
