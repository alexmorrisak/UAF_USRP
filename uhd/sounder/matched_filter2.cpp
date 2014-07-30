// matched_filter2.cpp

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
 * matched_filter2 function
 * Designed to work with complementary pulse codes
 **********************************************************************/
int matched_filter2(
    std::complex<float>*** indata,
    std::vector<std::complex<float> *> outdata,
    float** pcode,
    int pcode_length,
    int slowdim,
    int fastdim,
    int osr
){
    int ntaps = osr*pcode_length;
    //std::cout << "entering matched filter. ntaps: " << ntaps <<std::endl;
    //std::cout << slowdim <<std::endl;
    //std::cout << fastdim <<std::endl;
    std::vector<std::complex<float> > filter_taps[2];
    for (int i=0; i<2; i++){
        filter_taps[i].resize(ntaps);
    }
    for (int i=0; i<2; i++){
        for (int isym=0;isym<pcode_length;isym++){ //Create the over-sampled matched filter waveform
            for (int j=0; j<osr; j++){
                filter_taps[i][isym*osr+j] = std::complex<float>(
                    pcode[i][isym] / std::pow(ntaps,0.5),0);
                //std::cout << filter_taps[isym*osr+j] << std::endl;
            }
        }
    }

    for (int ipulse=0;ipulse<slowdim; ipulse++){ // Zero the output values
        for (int isamp=0; isamp<fastdim; isamp++){
            outdata[ipulse][isamp] = 0;
        }
    }

    for (int ipulse=0;ipulse<slowdim; ipulse++){
        for (int icode=0; icode<2; icode++){
            std::vector<std::complex<float> > tempvec(fastdim+ntaps,0);
            //populate the zero-padded version of the input samples
            for (int i=0; i<fastdim; i++){
                tempvec[ntaps/2+i] = indata[icode][ipulse][i];
            }
            for (int i=0; i<fastdim+ntaps; i++){
                //printf("in: %i (%.2f, %.2f)\n", i, tempvec[i].real(), tempvec[i].imag());
            }

            //perform the convolution
            for (int isamp=0; isamp<fastdim; isamp++){
                std::complex<float> temp(0,0);
                for (int i=0; i<ntaps; i++){
                    temp += filter_taps[icode][i]*tempvec[isamp+i];
                }
                outdata[ipulse][isamp] += temp;
                //printf("out %i,%i: %.2f\n",ipulse,isamp,10*log10(std::abs(outdata[ipulse][isamp])));
                //printf("temp %i,%i: %.2f\n",ipulse,isamp,10*log10(std::abs(temp)));
                //printf("temp %i,%i: %.2f @ %.0f\n",ipulse,isamp,std::abs(temp),180*std::arg(temp)/M_PI);
                //printf("out %i,%i: %.2f\n",ipulse,isamp,std::abs(outdata[ipulse][isamp]));
                //printf("out %i,%i: (%.1f, %.1f)\n",ipulse,isamp,outdata[ipulse][isamp].real(),
                //    outdata[ipulse][isamp].imag());
                //printf("temp %i,%i: (%.1f, %.1f)\n",ipulse,isamp,temp.real(),temp.imag());
            }
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
