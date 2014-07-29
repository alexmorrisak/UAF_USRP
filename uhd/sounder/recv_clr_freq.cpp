#include <uhd/utils/thread_priority.hpp>
#include <uhd/utils/safe_main.hpp>
#include <uhd/utils/static.hpp>
#include <uhd/usrp/multi_usrp.hpp>
#include <uhd/exception.hpp>
#include <boost/thread/thread.hpp>
#include <boost/program_options.hpp>
#include <boost/math/special_functions/round.hpp>
#include <boost/foreach.hpp>
#include <boost/format.hpp>
#include <boost/algorithm/string/split.hpp>
#include <boost/algorithm/string/classification.hpp>
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

typedef std::complex<int16_t>  sc16;
extern int verbose;

/***********************************************************************
 * clear_freq function
 **********************************************************************/
void recv_clr_freq(
    uhd::usrp::multi_usrp::sptr usrp,
    uhd::rx_streamer::sptr rx_stream,
    double *center_freq,
    double bandwidth
){
    if (verbose) std::cout << "Entering recv_clr_freq\n";
    int num_total_samps = 0;
    int nave=100;
    int nsamps = 100;
    double min_inx;
    fftw_complex *in=NULL,*out=NULL;
    fftw_plan plan;
    std::vector<double> tmp_pwr(nsamps,0);
    std::vector<double> pwr;
    pwr.reserve(nsamps);

    usrp->set_rx_rate(1e6);
    usrp->set_rx_freq(1e3*(*center_freq));
    //printf("%f\n",*center_freq);
  
    uhd::rx_metadata_t md;
    std::vector<std::complex<int16_t> *> buff_ptrs;
    std::vector<std::vector<std::complex<int16_t> > > buff;
    buff.resize(2); // Two receive channels
    buff_ptrs.resize(2); // Two receive channels
    for (int i=0; i<2; i++){
        buff[i].resize(nsamps);
        buff_ptrs[i] = &buff[i].front();
    }
    bool overflow_message = true;
    float timeout = usrp->get_time_now().get_real_secs() + 1.1;

    //setup streaming
    uhd::stream_cmd_t stream_cmd = uhd::stream_cmd_t::STREAM_MODE_NUM_SAMPS_AND_DONE;
    stream_cmd.num_samps = nave*nsamps;
    stream_cmd.stream_now = false;
    stream_cmd.time_spec = usrp->get_time_now() + 0.05;
    if (verbose) std::cout << "num_samps: " << nave*nsamps << std::endl;

    usrp->issue_stream_cmd(stream_cmd);
    md.error_code = uhd::rx_metadata_t::ERROR_CODE_NONE;
    while(num_total_samps != nave*nsamps){
	timeout = 0.1;
        size_t num_rx_samps = rx_stream->recv(buff_ptrs, buff[0].size(), md, timeout);
        if (md.error_code == uhd::rx_metadata_t::ERROR_CODE_TIMEOUT) {
            std::cerr << boost::format("Timeout while streaming") << std::endl;
            break;
        }
        if (md.error_code == uhd::rx_metadata_t::ERROR_CODE_OVERFLOW){
            continue;
        }
        if (md.error_code != uhd::rx_metadata_t::ERROR_CODE_NONE){
            throw std::runtime_error(str(boost::format(
                "Unexpected error code 0x%x"
            ) % md.error_code));
        }
    
    if(in!=NULL){free(in);in=NULL;}
    in =(fftw_complex*)fftw_malloc(sizeof(fftw_complex)*nsamps);
    if(out!=NULL){free(out);out=NULL;}
    out =(fftw_complex*)fftw_malloc(sizeof(fftw_complex)*nsamps);

    plan = fftw_plan_dft_1d(nsamps, in, out, FFTW_FORWARD, FFTW_ESTIMATE);

    for (int i=0;i<nsamps;i++){
       //in[i][0] = (hann_window[i]*rx_short_vecs[0][i].real());
       //in[i][1] = (hann_window[i]*rx_short_vecs[0][i].imag());
       in[i][0] = (double)buff[0][i].real();
       in[i][1] = (double)buff[0][i].imag();
    }
	
    fftw_execute(plan);

     //Add the current spectrum to the running total
    //if (num_total_samps != 0){
     for (int i=0;i<nsamps;i++){
             tmp_pwr[i] += (out[i][0]*out[i][0] + out[i][1]*out[i][1]) / 
                     (double)(nsamps*nsamps*nave);
     }
    //}
     //Done adding spectrum to running total
    num_total_samps += num_rx_samps;

    fftw_destroy_plan(plan);
    if (in!=NULL) {free(in); in=NULL;}
    if (out!=NULL) {free(out); out=NULL;}

   }

  //for (int i=0;i<nsamps;i++){
  //  printf("%i, %f\n",i, 10*log(tmp_pwr[i]));
  //}
  //Center the fft (fftshift)
  for(int i=0;i<(nsamps/2);i++){
      pwr[nsamps/2+i]=tmp_pwr[i];
      pwr[i]=tmp_pwr[nsamps/2+i];
  }
  //Done centering

  // Now find quietest frequency
  min_inx = (int) (nsamps/2 - bandwidth/2e1);
  int i;
  for (i = min_inx; i<(nsamps/2+bandwidth/2e1); i++){
    //printf("%f: %f\n",*center_freq + 10*i-500,log(pwr[i]));
    if (pwr[i] < pwr[min_inx])
    	min_inx = i;
  }
  *center_freq = *center_freq + (10*min_inx-500);

}

