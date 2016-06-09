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
#include <volk/volk.h>
#include "length_detect_c_impl.h"
#include <stdio.h>

#include <gnuradio/filter/fft_filter.h>

namespace gr {
  namespace burst {

    void conv(std::vector<gr_complex> &a, std::vector<gr_complex> &b, std::vector<gr_complex> &result)
    {
        // set up vectors
        result.resize(a.size()+b.size()-1);
        std::vector<gr_complex> taps;
        taps.resize(b.size());

        // set conj rev taps
        for(int i=0; i<b.size(); i++){
            taps[i] = std::conj(b[b.size()-i-1]);
            }

        // set up filter
        gr::filter::kernel::fft_filter_ccc filt(1, taps, 1); // decim, taps, nthreads
        filt.set_taps(taps);  // we can update taps if we want to keep the filter as a local var

        // perform convolution
        filt.filter( result.size(), &a[0], &result[0]);

        gr_complex q(1,1);
        q = q/1.0f;
    }

    length_detect_c::sptr
    length_detect_c::make()
    {
      return gnuradio::get_initial_sptr
        (new length_detect_c_impl());
    }

    /*
     * The private constructor
     */
    length_detect_c_impl::length_detect_c_impl()
      : gr::sync_block("length_detect_c",
              gr::io_signature::make(0,0,0),
              gr::io_signature::make(0,0,0)),
        d_cnt_good(0),
        d_cnt_bad(0)
    {
        message_port_register_in(pmt::mp("cpdus"));
        message_port_register_out(pmt::mp("cpdus"));
        set_msg_handler(pmt::mp("cpdus"), boost::bind(&length_detect_c_impl::handler, this, _1));
    }

    /*
     * Our virtual destructor.
     */
    length_detect_c_impl::~length_detect_c_impl()
    {
    }

    void length_detect_c_impl::handler(pmt::pmt_t msg){
        pmt::pmt_t meta = pmt::car(msg);
        std::vector<gr_complex> input(pmt::c32vector_elements(pmt::cdr(msg)));
        std::vector<float> mag(input.size());
        std::vector<float> magavg(input.size()-50);
        volk_32fc_magnitude_squared_32f(&mag[0], &input[0], input.size());

        // compute moving average
        float sum = 0;
        for(int i=0; i<50; i++){ sum += mag[i]; }
        for(int i=50; i<mag.size(); i++){
            sum += mag[i];
            sum -= mag[i-50];
            magavg[i-50] = sum;
        }

        // find max power / threshold power 6dB down
        uint16_t max_idx;
        volk_32f_index_max_16u(&max_idx, &magavg[0], magavg.size());
        float max_val = magavg[max_idx];
        float thresh_val = max_val / 16;
       
        // find when we go under this thresh after the burst
        int start_offset = 1000;
        int end_offset = 0;
        for(int i=start_offset; i<magavg.size(); i++){
            //printf("%f %f\n", magavg[i], thresh_val);
            if(magavg[i] < thresh_val){
                end_offset = i+start_offset+100;
                break;
            }
        }
        
        //printf("good=%lu bad=%lu\n", d_cnt_good, d_cnt_bad);

        // discard burst if we dont have a good ending ...
        if(end_offset <= 100+start_offset*2){
            d_cnt_bad++;
            //printf("WARNING: length_detect_c: discarding burst");
            return;
        }
        d_cnt_good++;
 
        // add stuff to metadata            
        meta = pmt::dict_add(meta, pmt::intern("est_len"), pmt::mp(end_offset));  
        meta = pmt::dict_add(meta, pmt::intern("orig_len"), pmt::mp((int)mag.size()));  

        // set up the output vector and send
        pmt::pmt_t data_out(pmt::init_c32vector(end_offset, input));
        message_port_pub(pmt::intern("cpdus"), pmt::cons(meta, data_out));
    }

    int
    length_detect_c_impl::work(int noutput_items,
			  gr_vector_const_void_star &input_items,
			  gr_vector_void_star &output_items)
    {
        return noutput_items;
    }

  } /* namespace burst */
} /* namespace gr */

