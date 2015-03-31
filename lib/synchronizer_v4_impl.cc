/* -*- c++ -*- */
/* 
 * Copyright 2015 Free Software Foundation, Inc
 * 
 * This is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3, or (at your option)
 * any later version.
 * 
 * This software is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this software; see the file COPYING.  If not, write to
 * the Free Software Foundation, Inc., 51 Franklin Street,
 * Boston, MA 02110-1301, USA.
 */

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "qa_helpers.h"

#include <gnuradio/fft/fft.h>
#include <gnuradio/io_signature.h>
#include "synchronizer_v4_impl.h"
#include <math.h>
#include <stdio.h>

#include <volk/volk_typedefs.h>
#include <volk/volk.h>

namespace gr {
  namespace burst {

    const gr_complex synchronizer_v4_impl::preSyms_fliplr_conj[48] =
        {0.70710677-0.70710677j,  0.70710677+0.70710677j,
        -0.70710677-0.70710677j, -0.70710677-0.70710677j,
        -0.70710677-0.70710677j, -0.70710677+0.70710677j,
        -0.70710677+0.70710677j, -0.70710677+0.70710677j,
         0.70710677-0.70710677j,  0.70710677-0.70710677j,
         0.70710677-0.70710677j, -0.70710677+0.70710677j,
         0.70710677-0.70710677j, -0.70710677+0.70710677j,
        -0.70710677+0.70710677j,  0.70710677-0.70710677j,
        -0.70710677-0.70710677j, -0.70710677+0.70710677j,
        -0.70710677+0.70710677j, -0.70710677-0.70710677j,
        -0.70710677-0.70710677j,  0.70710677+0.70710677j,
         0.70710677-0.70710677j,  0.70710677-0.70710677j,
         0.70710677+0.70710677j,  0.70710677-0.70710677j,
        -0.70710677+0.70710677j,  0.70710677+0.70710677j,
        -0.70710677-0.70710677j,  0.70710677+0.70710677j,
        -0.70710677-0.70710677j,  0.70710677+0.70710677j,
        -0.70710677-0.70710677j,  0.70710677+0.70710677j,
         0.70710677-0.70710677j,  0.70710677-0.70710677j,
         0.70710677+0.70710677j, -0.70710677+0.70710677j,
        -0.70710677+0.70710677j, -0.70710677-0.70710677j,
         0.70710677+0.70710677j, -0.70710677-0.70710677j,
        -0.70710677-0.70710677j, -0.70710677+0.70710677j,
        -0.70710677+0.70710677j,  0.70710677+0.70710677j,
         0.70710677-0.70710677j, -0.70710677-0.70710677j};

    const gr_complex synchronizer_v4_impl::preSyms_x2_fliplr_conj[96] =
        {-0.05751960-1.06877923j,  0.70710665-0.70710665j,
			0.96787697+0.31263849j,  0.70710683+0.70710689j,
			0.04994498+0.02322075j, -0.70710677-0.70710677j,
		   -1.00477576-0.77428943j, -0.70710683-0.70710677j,
		   -0.45783547-0.87726957j, -0.70710671-0.70710683j,
		   -0.95984519+0.06683378j, -0.70710677+0.70710671j,
		   -0.42076343+0.74770236j, -0.70710671+0.70710665j,
		   -1.08900166+0.88655835j, -0.70710683+0.70710677j,
			0.21010289-0.073215j  ,  0.70710683-0.70710671j,
			0.64430445-0.73423678j,  0.70710671-0.70710665j,
			0.97398925-0.92401201j,  0.70710683-0.70710677j,
		   -0.20128489+0.18980072j, -0.70710671+0.70710677j,
		   -0.11851995+0.08858704j,  0.70710665-0.70710677j,
			0.46924281-0.38897231j, -0.70710671+0.70710677j,
		   -1.37390292+1.22124052j, -0.70710677+0.70710677j,
			0.43081510-0.13785198j,  0.70710677-0.70710671j,
			0.00463670-0.91634619j, -0.70710677-0.70710683j,
		   -0.83524835-0.04208908j, -0.70710671+0.70710653j,
		   -0.70443785+1.00534058j, -0.70710683+0.70710683j,
		   -0.62667698+0.03305677j, -0.70710677-0.70710677j,
		   -0.90893525-1.08450902j, -0.70710671-0.70710683j,
			0.04223688+0.18974672j,  0.70710665+0.70710671j,
			0.81873888+0.25460532j,  0.70710665-0.70710671j,
			0.73506474-1.18435335j,  0.70710683-0.70710671j,
			0.56178945+0.24885255j,  0.70710671+0.70710671j,
			1.03661072+0.15496387j,  0.70710683-0.70710671j,
		   -0.30967343-0.56743401j, -0.70710683+0.70710677j,
			0.08827468+1.53255141j,  0.70710671+0.70710671j,
			0.05508756-0.6689443j , -0.70710671-0.70710671j,
		   -0.17450830+0.38476533j,  0.70710677+0.70710677j,
			0.28946680-0.20185675j, -0.70710677-0.70710683j,
		   -0.41400599+0.04977918j,  0.70710689+0.70710671j,
			0.56852669+0.10137521j, -0.70710677-0.70710677j,
		   -0.80945057-0.28088036j,  0.70710671+0.70710671j,
			1.56957042+0.55650836j,  0.70710683-0.70710671j,
		   -0.03159594-1.3975693j ,  0.70710683-0.70710671j,
			1.46415174+0.35750556j,  0.70710665+0.70710683j,
		   -0.53795636+0.5759654j , -0.70710671+0.70710671j,
		   -0.34947091+0.99419743j, -0.70710677+0.70710677j,
		   -1.24164104-0.18950649j, -0.70710683-0.70710689j,
			0.43727133-0.15147361j,  0.70710665+0.70710689j,
		   -0.04705258+0.51324481j, -0.70710671-0.70710677j,
		   -0.77846509-1.41305459j, -0.70710677-0.70710677j,
		   -0.73197603+0.41725174j, -0.70710683+0.70710677j,
		   -0.69809616+0.48355645j, -0.70710671+0.70710677j,
		   -0.27108425+1.10157359j,  0.70710665+0.70710677j,
			1.29566348-0.25711846j,  0.70710683-0.70710683j,
		   -0.38806984-0.5718742j , -0.70710683-0.70710677j};

    synchronizer_v4::sptr
    synchronizer_v4::make(double Fs, int sps)
    {
      return gnuradio::get_initial_sptr
        (new synchronizer_v4_impl(Fs, sps));
    }

    /*
     * The private constructor
     */
    synchronizer_v4_impl::synchronizer_v4_impl(double Fs, int sps)
      : gr::sync_block("synchronizer_v4",
              gr::io_signature::make(0, 0, 0),
              gr::io_signature::make(0, 0, 0)), preFFTEngine(128), preIFFTEngine(128), wOpt_gr(48)
    {
    	message_port_register_in(pmt::mp("cpdus"));
		message_port_register_out(pmt::mp("cpdus"));
		message_port_register_out(pmt::mp("debug_post_cfo"));
		message_port_register_out(pmt::mp("debug_pre_xcorr"));
		set_msg_handler(pmt::mp("cpdus"), boost::bind(&synchronizer_v4_impl::handler, this, _1));

        std::cout << "v4 sync\n";

    	d_Fs = Fs;
    	d_sps = sps;

    	preSymsSize = 48;
    	preSymsX2Size = 96;
    	preFFTEngineFFTSize = 128;

    	wOpt = gsl_vector_complex_alloc(48);

    	debugMode = false;
    }

    /*
     * Our virtual destructor.
     */
    synchronizer_v4_impl::~synchronizer_v4_impl()
    {
    	gsl_vector_complex_free(wOpt);
    }

    void synchronizer_v4_impl::enableDebugMode() {
    	debugMode = true;
    }

    void synchronizer_v4_impl::qpskBurstCFOCorrect(gr_complex* x_in, int burstSize) {
		// TODO: maybe we can have a fixed burst size, only calculate cfo over a certain burst size
		// then we can statically create the fft engine in the constructor, will probably be a lot faster
		if(burstSize%2!=0) {
			burstSize = burstSize-1;
		}
		fft::fft_complex fftEngine = fft::fft_complex(burstSize);
		gr_complex* fftInBuf = fftEngine.get_inbuf();
		gr_complex* fftOutBuf = fftEngine.get_outbuf();
		std::vector<float> fftOutAbs(burstSize, 0.0);
		// fill the FFT buffer
		volk_32fc_s32f_power_32fc(fftInBuf, x_in, 4, burstSize);

		if(debugMode) {
			std::string filename = "/tmp/gr_cfoFFTInput.txt";
			std::vector<gr_complex> b(burstSize);
			b.assign(fftInBuf, fftInBuf+burstSize);
			qa_helpers::writeComplexFile(filename, b);
		}

		// compute the fft
		fftEngine.execute();
		// take absolute value of fft output and find argmax simultaneously
		int maxIdx = 0;
		float maxVal = -100000;
		for(int ii=0; ii<burstSize; ii++) {
			fftOutAbs[ii] = pow(fftOutBuf[ii].real(), 2) + pow(fftOutBuf[ii].imag(), 2);
			if(fftOutAbs[ii]>maxVal) {
				maxVal = fftOutAbs[ii];
				maxIdx = ii;
			}
		}

		if(debugMode) {
			std::string filename = "/tmp/gr_cfoFFTOut_abs.txt";
			qa_helpers::writeFloatFile(filename, fftOutAbs);
		}

		// account for fftshift
		if(maxIdx<burstSize/2) {
			maxIdx = maxIdx + burstSize/2;
		}
		else {
			maxIdx = maxIdx - burstSize/2;
		}
		double cfoEstimate = (-d_Fs/2.0+maxIdx*(d_Fs)/(burstSize-1)) / 4.0;

		if(debugMode) {
			std::string filename = "/tmp/gr_cfoEstimate.txt";
			std::vector<float> cfoEstimateVec(1);
			cfoEstimateVec[0] = cfoEstimate;
			qa_helpers::writeFloatFile(filename, cfoEstimateVec);
		}

		// correct for cfo
		shiftFreq(&x_in[0], burstSize, d_Fs, cfoEstimate, 0);
	}

    void synchronizer_v4_impl::shiftFreq(gr_complex* buf, int bufLen, double Fs, double freq, double tStart) {
		std::complex<float> phInc(cos(2*M_PI*-freq*1.0/Fs),sin(2*M_PI*-freq*1.0/Fs));
		std::complex<float> phStart(cos(2*M_PI*-freq*tStart),sin(2*M_PI*-freq*tStart));
		volk_32fc_s32fc_x2_rotator_32fc(buf, buf, phInc, &phStart, bufLen);
	}

    void synchronizer_v4_impl::conv(gr_complex* a, int aLen, const gr_complex* b, int bLen, std::vector<gr_complex> &result)
	{
    	int fftSize = aLen+bLen-1;
		result.resize(fftSize);

		fft::fft_complex fftEngine1 = fft::fft_complex(fftSize);
		gr_complex* fftInBuf = fftEngine1.get_inbuf();
		memcpy(fftInBuf, a, sizeof(gr_complex)*aLen);
		memset(fftInBuf+aLen, 0, sizeof(gr_complex)*(fftSize-aLen));
		fftEngine1.execute();
		gr_complex* fft1OutBuf = fftEngine1.get_outbuf();

		fft::fft_complex fftEngine2 = fft::fft_complex(fftSize);
		fftInBuf = fftEngine2.get_inbuf();
		memcpy(fftInBuf, b, sizeof(gr_complex)*bLen);
		memset(fftInBuf+bLen, 0, sizeof(gr_complex)*(fftSize-bLen));
		fftEngine2.execute();
		gr_complex* fft2OutBuf = fftEngine2.get_outbuf();

		// multiply the fft outputs
		fft::fft_complex ifftEngine = fft::fft_complex(fftSize, false);
		gr_complex* ifftInBuf = ifftEngine.get_inbuf();
		volk_32fc_x2_multiply_32fc(ifftInBuf, fft1OutBuf, fft2OutBuf, fftSize);

		// ifft
		ifftEngine.execute();

		// scale back
		volk_32fc_s32fc_multiply_32fc(&result[0], ifftEngine.get_outbuf(), 1.0/fftSize, fftSize);
	}

    void synchronizer_v4_impl::determineOptimalFilter(gsl_vector_complex* w, gr_complex* x, int xLen) {
    	gr_complex* fftInBuf = preFFTEngine.get_inbuf();
    	// fill the fft input buffer
    	memcpy(&fftInBuf[0], &x[0], 48*sizeof(gr_complex));
    	memset(&fftInBuf[48], 0, 80*sizeof(gr_complex));

    	preFFTEngine.execute();
    	gr_complex* fftOutBuf = preFFTEngine.get_outbuf();

    	if(debugMode) {
			std::string filename = "/tmp/gr_dofFFTOutput.txt";
			std::vector<gr_complex> b(128);
			b.assign(fftOutBuf, fftOutBuf+128);
			qa_helpers::writeComplexFile(filename, b);
		}

    	// take the output, and store the absolute value in an array
    	gr_complex* ifftInBuf = preIFFTEngine.get_inbuf();
    	for(int ii=0; ii<preFFTEngineFFTSize; ii++) {
    		ifftInBuf[ii].real() = pow( std::abs(fftOutBuf[ii]), 2);
    		ifftInBuf[ii].imag() = 0;
    	}

    	if(debugMode) {
			std::string filename = "/tmp/gr_dofIFFTInput.txt";
			std::vector<float> b(128);
			for(int ii=0; ii<128; ii++) {
				b[ii] = ifftInBuf[ii].real();
			}
			qa_helpers::writeFloatFile(filename, b);
		}

    	preIFFTEngine.execute();
    	gr_complex* ifftOutBuf = preIFFTEngine.get_outbuf();

    	if(debugMode) {
			std::string filename = "/tmp/gr_dofIFFTOutput.txt";
			std::vector<gr_complex> b(128);
			b.assign(ifftOutBuf, ifftOutBuf+128);
			qa_helpers::writeComplexFile(filename, b);
		}

    	// generate the row and col vectors for toeplitz matrix creation and
    	// scale as necessary
		std::vector<gr_complex> row(preSymsSize);
		std::vector<gr_complex> col(preSymsSize);
    	for(int ii=0; ii<preSymsSize; ii++) {
    		// downscale by m and the ifft size
    		col[ii] = std::conj(ifftOutBuf[ii])*1.0f/(48.0f*128.0f);
    		row[ii] = ifftOutBuf[ii]*1.0f/(48.0f*128.0f);
    	}

    	if(debugMode) {
			std::string filename = "/tmp/gr_dofToeplitzCol.txt";
			qa_helpers::writeComplexFile(filename, col);
			std::string filename2 = "/tmp/gr_dofToeplitzRow.txt";
			qa_helpers::writeComplexFile(filename2, row);
		}

    	gsl_matrix_complex* R = gsl_matrix_complex_alloc(preSymsSize, preSymsSize);
    	toeplitz(&col[0], preSymsSize, &row[0], preSymsSize, R);

    	std::vector<gr_complex> xc(95);
    	// compute correlation between preSyms and x[1:48]
    	std::vector<gr_complex> x_in(144, gr_complex(0,0));
    	memcpy(&x_in[0], &x[0], sizeof(gr_complex)*48);
    	// the difference in the cross-correlations between the zero-pad and the non-zeropad are small,
    	// not sure if we can get away w/ doing no zeropad?? investigate w/ perofrmacne
    	conv(&x[0], 48, preSyms_fliplr_conj, 48, xc);

    	if(debugMode) {
			std::string filename = "/tmp/gr_dofxc.txt";
			qa_helpers::writeComplexFile(filename, xc);

			// do this to make sure we load P properly as debugging measure
			std::vector<gr_complex> P(48);
			int jj = xc.size()-1;
			for(int ii=0; ii<48; ii++) {
				P[ii] = gr_complex(xc[jj].real(), -xc[jj].imag() );
				jj--;
			}
			std::string filename2 = "/tmp/gr_P.txt";
			qa_helpers::writeComplexFile(filename2, P);
		}

    	// make P vector
        gsl_vector_complex* P = gsl_vector_complex_alloc(48);
    	gsl_complex cval;
    	int jj = xc.size()-1;
    	for(int ii=0; ii<48; ii++) {
    		GSL_SET_COMPLEX(&cval, xc[jj].real(), -xc[jj].imag());
    		gsl_vector_complex_set(P, ii, cval);
    		jj--;
    	}

    	// solve for R
    	int s;
    	gsl_permutation * p = gsl_permutation_alloc(48);
    
        // sizes ...
        // R = [ 48x48 ]
        // p = [ 48x1 ]
        // s = [ 1x1 ]
        // P = [ 48x1 ]
        // w (wOpt) = [48x1]

    	gsl_linalg_complex_LU_decomp (R, p, &s);
    	gsl_linalg_complex_LU_solve (R, p, P, w);

    	// rescale w happens when we copy w back to a gr vector

    	// free gsl memory
    	gsl_matrix_complex_free(R);
		gsl_vector_complex_free(P);
		gsl_permutation_free(p);
    }


    void synchronizer_v4_impl::toeplitz(gr_complex* col, int M, gr_complex* row, int N, gsl_matrix_complex* T) {
    	// main diagonal and below the main diagonal
    	gsl_complex colVal, rowVal;
		for( int d=0; d<M; ++d ) {
			for( int i=0; i<M-d; ++i ) {
				GSL_SET_COMPLEX(&colVal, col[d].real(), col[d].imag());
				gsl_matrix_complex_set(T,i+d,i,colVal);
			}
		}

		// above the main diagonal
		for( int d=1; d<N; ++d ) {
			for( int i=0; i<N-d; ++i ) {
				GSL_SET_COMPLEX(&rowVal, row[d].real(), row[d].imag());
				gsl_matrix_complex_set(T,i,i+d,rowVal);
			}
		}
    }

    void synchronizer_v4_impl::qpskFirstOrderPLL(gr_complex* x, int size, float alpha, gr_complex* y) {
    	gr_complex phiHat = gr_complex(1,0);
    	gr_complex xHat, er, phiHatT;
	    for(int ii=0; ii<size; ii++) {
			// correct w/ estimated phase
            y[ii] = x[ii]*phiHat;

            // demodulating circuit
			if(y[ii].real()>=0 && y[ii].imag()>=0) {
				xHat.real() = M_SQRT1_2;
				xHat.imag() = M_SQRT1_2;
			}
			else if(y[ii].real()>=0 && y[ii].imag()<0) {
				xHat.real() = M_SQRT1_2;
				xHat.imag() = -M_SQRT1_2;
			}
			else if(y[ii].real()<0 && y[ii].imag()<0) {
				xHat.real() = -M_SQRT1_2;
				xHat.imag() = -M_SQRT1_2;
			}
			else {
				xHat.real() = -M_SQRT1_2;
				xHat.imag() = M_SQRT1_2;
			}

			// loop filter to update phase estimate
			er = std::conj(xHat)*y[ii];
			phiHatT = er/std::abs(er);
			phiHat = std::conj( std::pow( phiHatT, alpha)) * phiHat;
	    }
    }

    void synchronizer_v4_impl::handler(pmt::pmt_t msg)
	{
		// get parameters from the pdu
		std::vector<gr_complex> eqBurst(pmt::c32vector_elements(pmt::cdr(msg)));
		pmt::pmt_t meta = pmt::car(msg);

		if(debugMode) {
			std::string filename = "/tmp/gr_inputData.txt";
			qa_helpers::writeComplexFile(filename, eqBurst);
		}

		// perform cfo correction
		qpskBurstCFOCorrect(&eqBurst[0], eqBurst.size());

		if(debugMode) {
			std::string filename = "/tmp/gr_burstCFOCorrected.txt";
			qa_helpers::writeComplexFile(filename, eqBurst);
		}


        // publish debug port #1
        if(1){
            pmt::pmt_t cfo_vec = pmt::init_c32vector(eqBurst.size(), &eqBurst[0]);
            message_port_pub(pmt::mp("debug_post_cfo"), pmt::cons( meta, cfo_vec ) );
        }

        // search for start index/timing
		std::vector<gr_complex> preCrossCorr_cmplx(96+eqBurst.size()-1);
		conv(&eqBurst[0], eqBurst.size(), preSyms_x2_fliplr_conj, 96, preCrossCorr_cmplx);

		std::vector<float> preCrossCorr(preCrossCorr_cmplx.size());
		int maxIdx = 0;
		float maxVal = -99999;
		for(int ii=0; ii<preCrossCorr.size(); ii++) {
			preCrossCorr[ii] = std::abs(preCrossCorr_cmplx[ii]);
			if(preCrossCorr[ii]>maxVal) {
				maxIdx = ii;
				maxVal = preCrossCorr[ii];
			}
		}

		int preambleIdxStart = maxIdx - 96 + 1;
		if(debugMode) {
			std::string f1 = "/tmp/gr_preCrossCorr.txt";
			qa_helpers::writeFloatFile(f1, preCrossCorr);

			std::string filename = "/tmp/gr_preCrossCorrIdx.txt";
			std::vector<float> preambleIdxStartVec(1);
			preambleIdxStartVec[0] = preambleIdxStart;
			qa_helpers::writeFloatFile(filename, preambleIdxStartVec);
		}

        // publish debug port #2
        if(1){
            pmt::pmt_t xcorr_vec = pmt::init_f32vector(preCrossCorr.size(), &preCrossCorr[0]);
            message_port_pub(pmt::mp("debug_pre_xcorr"), pmt::cons( meta, xcorr_vec ) );
        }

		if(preambleIdxStart<0) {
			// means we didn't find the preamble, quit
			return;
		}

		int eqIn_decimated_len = (eqBurst.size() - preambleIdxStart + 1)/2;		// zero_padding
		std::vector<gr_complex> eqIn_decimated(eqIn_decimated_len, gr_complex(0,0));
		// decimate the signal
		for(int ii=0; ii<eqIn_decimated.size(); ii++) {
			eqIn_decimated[ii] = eqBurst[preambleIdxStart];
			preambleIdxStart+=2;
		}
		if(debugMode) {
			std::string filename = "/tmp/gr_eqInDecimated.txt";
			qa_helpers::writeComplexFile(filename, eqIn_decimated);
		}

//		gsl_vector_complex* wOpt = gsl_vector_complex_alloc(48);
		determineOptimalFilter(wOpt, &eqIn_decimated[0], eqIn_decimated.size());

		// filter
		std::vector<gr_complex> whFilt(eqIn_decimated.size()+48-1);
		float maxWopt = -999999;
		float tmp;
		for(int ii=0; ii<48; ii++) {
			gsl_complex c = gsl_vector_complex_get(wOpt, ii);
			wOpt_gr[ii] = gr_complex(GSL_REAL(c), GSL_IMAG(c));
			tmp = std::abs(wOpt_gr[ii]);
			if(tmp>maxWopt) {
				maxWopt = tmp;
			}
		}
//		gsl_vector_complex_free(wOpt);

		if(debugMode) {
			std::string filename = "/tmp/gr_wOpt.txt";
			qa_helpers::writeComplexFile(filename, wOpt_gr);
		}
		std::complex<float> maxWoptCmplx(1.0/maxWopt,0);
		volk_32fc_s32fc_multiply_32fc(&wOpt_gr[0], &wOpt_gr[0], maxWoptCmplx, wOpt_gr.size());

		conv(&eqIn_decimated[0], eqIn_decimated.size(), &wOpt_gr[0], wOpt_gr.size(), whFilt);

		if(debugMode) {
			std::string filename2 = "/tmp/gr_wOpt_scaled.txt";
			qa_helpers::writeComplexFile(filename2, wOpt_gr);

			std::string filename3 = "/tmp/gr_whFilt.txt";
			qa_helpers::writeComplexFile(filename3, whFilt);
		}


		// normalize
		std::vector<float> whFiltAbs(whFilt.size());
	    volk_32fc_magnitude_32f(&whFiltAbs[0], &whFilt[0], whFilt.size());
	    float normFactor;
	    volk_32f_accumulator_s32f(&normFactor, &whFiltAbs[0], whFilt.size());
	    normFactor /= whFilt.size();

	    std::complex<float> normFactorCmplx(1.0/normFactor,0);
	    volk_32fc_s32fc_multiply_32fc(&whFilt[0], &whFilt[0], normFactorCmplx, whFilt.size());

		if(debugMode) {
			std::string filename = "/tmp/gr_whFilt_scaled.txt";
			qa_helpers::writeComplexFile(filename, whFilt);
		}

		// apply pll
		std::vector<gr_complex> phRecoveredSyms(whFilt.size());
		float alpha = 0.002;
		qpskFirstOrderPLL(&whFilt[0], whFilt.size(), alpha, &phRecoveredSyms[0]);

		if(debugMode) {
			std::string filename = "/tmp/gr_phRecoveredSyms.txt";
			qa_helpers::writeComplexFile(filename, phRecoveredSyms);
		}

		// put into new pdu and send
		pmt::pmt_t newvec = pmt::init_c32vector(phRecoveredSyms.size(), &phRecoveredSyms[0]);
		msg = pmt::cons( meta, newvec );
		message_port_pub(pmt::mp("cpdus"), msg);
	}

    int
    synchronizer_v4_impl::work(int noutput_items,
			  gr_vector_const_void_star &input_items,
			  gr_vector_void_star &output_items)
    {
        // Do <+signal processing+>

        // Tell runtime system how many output items we produced.
        return noutput_items;
    }

  } /* namespace burst */
} /* namespace gr */

