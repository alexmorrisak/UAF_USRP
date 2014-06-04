// matched_filter.cpp

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
 * matched_filter function
 **********************************************************************/
int matched_filter(
    std::vector<std::complex<float> *> indata,
    std::vector<std::complex<float> *> outdata,
    float *pcode,
    int pcode_length,
    int slowdim,
    int fastdim,
    int osr
){
    int ntaps = osr*pcode_length;
    //std::cout << "entering matched filter. ntaps: " << ntaps <<std::endl;
    //std::cout << slowdim <<std::endl;
    //std::cout << fastdim <<std::endl;
    std::vector<std::complex<float> > filter_taps(ntaps);
    for (int isym=0;isym<pcode_length;isym++){
        for (int j=0; j<osr; j++){
            filter_taps[isym*osr+j] = std::complex<float>(
                pcode[isym] / std::pow(ntaps,0.5),0);
            //std::cout << filter_taps[isym*osr+j] << std::endl;
        }
    }

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

        //perform the convolution
        for (int isamp =0; isamp<fastdim; isamp++){
            std::complex<float> temp(0,0);
            for (int i=0; i<ntaps; i++){
                temp += filter_taps[i]*tempvec[isamp+i];
            }
            outdata[ipulse][isamp] = temp;
            //printf("out %i,%i: %.2f\n",ipulse,isamp,10*log10(std::abs(outdata[ipulse][isamp])));
            //printf("out %i,%i: (%.1f, %.1f)\n",ipulse,isamp,outdata[ipulse][isamp].real(),
            //    outdata[ipulse][isamp].imag());
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
//    invecs[0].resize(nsamps,0.);
//    //for (int i=0; i<invecs[0].size(); i++){
//    //    invecs[0][i] = std::complex<float>(3*i,0);
//    //}
//    invecs[0][10] = 10;
//    invecs[0][11] = 10;
//    invecs[0][12] = -10;
//    in.resize(nsamps);
//    in[0] = &invecs[0].front();
//
//    outvecs[0].resize(nsamps);
//    out.resize(nsamps);
//    out[0] = &outvecs[0].front();
//
//    int pcode_len = 1;
//    //float pcode[13] = {1,1,1,1,1,-1,-1,1,1,-1,1,-1,1};
//    //float pcode[3] = {1,1,-1};
//    float pcode[1] = {1};
//
//    //int matched_filter(
//    //    std::vector<std::complex<float> *> indata,
//    //    std::vector<std::complex<float> *> outdata,
//    //    std::vector<float *> pcode,
//    //    int pcode_length,
//    //    int slowdim,
//    //    int fastdim,
//    //    float osr
//    //){
//    int rval = matched_filter(
//        in,
//        out,
//        &pcode[0],
//        pcode_len,
//        1,
//        nsamps,
//        1);
//
//    printf("ouputs\n");
//    for (int i=0; i<nsamps; i++){
//        printf("%i (%.2f, %.2f)\n", i, out[0][i].real(), out[0][i].imag());
//    }
//}
