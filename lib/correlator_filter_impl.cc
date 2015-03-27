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

#include <gnuradio/io_signature.h>
#include "correlator_filter_impl.h"
#include <stdio.h>

namespace gr {
  namespace burst {

    correlator_filter::sptr
    correlator_filter::make(int l1, float i1, float thresh)
    {
      return gnuradio::get_initial_sptr
        (new correlator_filter_impl(l1, i1,  thresh));
    }

    /*
     * The private constructor
     */
    correlator_filter_impl::correlator_filter_impl(int l1, float i1, float thresh)
      : gr::sync_block("correlator_filter",
              gr::io_signature::make(1,1, sizeof(float)),
              gr::io_signature::make(1,1, sizeof(float))),
    d_l1(l1), d_i1(i1),
    d_ma_mu(0.0), d_ma_sigma(0.0), d_thresh(thresh),
    d_last_event(1000)
    {
        set_history(d_l1+1);
    }

    correlator_filter_impl::~correlator_filter_impl()
    {
    }

    int
    correlator_filter_impl::work(int noutput_items,
			  gr_vector_const_void_star &input_items,
			  gr_vector_void_star &output_items)
    {
        const float *in = (const float *) input_items[0];
        float *out = (float*) output_items[0];


        for(size_t i=0; i<noutput_items; i++){
            int now = i + history()-1;

            d_ma_mu += in[now];
            d_ma_mu -= in[now-d_l1];
            
            float mu = d_ma_mu / d_l1;


            //float iir_update = mu / in[now]; // ( compute geometric deviation from mean?)
            float iir_update = abs(mu - in[now]); // ( compute geometric deviation from mean?)
            //float d_ma_sigma = (iir_update*d_i1) + (d_ma_sigma)*(1.0f-d_i1); // normal IIR update routine
            float d_ma_sigma = (iir_update*1e-4) + (d_ma_sigma)*(1.0f-1e-4); // normal IIR update routine
            float ref_level = d_ma_sigma + mu;  // reference level ... (Switch from geometric to linear??)
            float stat = in[now] / ref_level + d_thresh;
//            printf("iir_update = %f, d_ma_sigma = %f, ref_level = %f, stat = %f\n", iir_update, d_ma_sigma, ref_level, stat);
            out[i] = stat;
            }    
        //throw std::runtime_error("done");

        return noutput_items;
    }

  } /* namespace burst */
} /* namespace gr */

