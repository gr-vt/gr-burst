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

#ifndef INCLUDED_BURST_SYNCHRONIZER_V4_IMPL_H
#define INCLUDED_BURST_SYNCHRONIZER_V4_IMPL_H

#include <gnuradio/fft/fft.h>
#include <burst/synchronizer_v4.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_matrix_complex_double.h>
#include <gsl/gsl_vector_complex_double.h>
#include <gsl/gsl_linalg.h>
#include <gnuradio/filter/fft_filter.h>
#include <mapper/constellation.h>

namespace gr {
  namespace burst {

    class synchronizer_v4_impl : public synchronizer_v4
    {
     private:
      float qpskBurstCFOCorrect(gr_complex* x, int burstSize);
      void shiftFreq(gr_complex* buf, int bufLen, double Fs, double freq, double tStart);
      void determineOptimalFilter(gsl_vector_complex* w, gr_complex* x, int xLen);
      void qpskFirstOrderPLL(gr_complex* x, int size, float alpha, gr_complex* y);
      void toeplitz(gr_complex* col, int M, gr_complex* row, int N, gsl_matrix_complex* T);
      void conv(gr_complex* a, int aLen, const gr_complex* b, int bLen, std::vector<gr_complex> &result);
      void handler(pmt::pmt_t msg);

      void enableDebugMode();

     public:
      synchronizer_v4_impl(double Fs, int sps, std::vector<unsigned char> preamble_bits, std::vector<int> sym_mapping);
      ~synchronizer_v4_impl();

      // Where all the action really happens
      int work(int noutput_items,
	       gr_vector_const_void_star &input_items,
	       gr_vector_void_star &output_items);

      float d_Fs;
      int d_sps;

      int preFFTEngineFFTSize;
      fft::fft_complex preFFTEngine;

      std::vector<gr_complex> preSyms_fliplr_conj;			// preamble bits size / 2
      std::vector<gr_complex> preSyms_xR_fliplr_conj;		// ( preamble bits size / 2 )* sps

      int optimalFilterSize;
      std::vector<gr_complex> wOpt_gr;
      gsl_vector_complex* wOpt;

      int preSymsSize;
      int preSymsRateMatchedSize;

      mapper::constellation d_const;

      bool debugMode;

    };

  } // namespace burst
} // namespace gr

#endif /* INCLUDED_BURST_SYNCHRONIZER_V4_IMPL_H */

