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
#include "interp_cpdu_impl.h"
#include <gnuradio/filter/interp_fir_filter_ccf.h>
#include <gnuradio/filter/firdes.h>
#include <volk/volk.h>

namespace gr {
  namespace burst {

    interp_cpdu::sptr
    interp_cpdu::make(double fs)
    {
      return gnuradio::get_initial_sptr
        (new interp_cpdu_impl(fs));
    }

    interp_cpdu_impl::interp_cpdu_impl(double fs)
      : gr::sync_block("interp_cpdu",
              gr::io_signature::make(0,0,0),
              gr::io_signature::make(0,0,0)),
        d_fs(fs)
    {
        message_port_register_in(pmt::mp("cpdus"));
        message_port_register_out(pmt::mp("cpdus"));
        set_msg_handler(pmt::mp("cpdus"), boost::bind(&interp_cpdu_impl::interp, this, _1));
    }

    interp_cpdu_impl::~interp_cpdu_impl()
    {
    }

    void interp_cpdu_impl::interp(pmt::pmt_t msg)
    {
        // get parameters from the pdu
        std::vector<gr_complex> input(pmt::c32vector_elements(pmt::cdr(msg)));
        pmt::pmt_t meta = pmt::car(msg);
        int interp = pmt::to_long(pmt::dict_ref(meta, pmt::mp("interp"), pmt::PMT_NIL));
        double offset = pmt::to_double(pmt::dict_ref(meta, pmt::mp("freq_offset"), pmt::PMT_NIL));
        if(interp < 2){ throw std::runtime_error("bad interp value"); }
        std::vector<gr_complex> output(interp*input.size());
//        printf("interp = %d, offset = %f\n", interp, offset);
 
        // make filter
        std::vector<float> taps(gr::filter::firdes::low_pass(1.0, d_fs, d_fs/(2*interp), d_fs/(2*interp*4)));
//        printf("fs = %f, cutoff = %f, trans = %f (ntaps = %d)\n", d_fs,  d_fs/(2*interp), d_fs/(2*interp*4), taps.size());
//        printf("in len = %d\n", input.size());
        gr::filter::interp_fir_filter_ccf::sptr filt(gr::filter::interp_fir_filter_ccf::make(interp,taps));

        // zero pad the input
        std::vector<gr_complex> zpad_input(taps.size()-1, gr_complex(0,0));
        zpad_input.insert(zpad_input.end(), input.begin(), input.end());

        // perform filtering
        //gr_vector_const_void_star win(1, (const void*) &input[0] );
        gr_vector_const_void_star win(1, (const void*) &zpad_input[0] );
        gr_vector_void_star wout(1, (void*) &output[0]);
        //int nout = filt->work(input.size() * interp - taps.size() + 1, win, wout);
        int nout_exp = output.size();
//        printf("nout_exp = %d\n", nout_exp);
        int nout = filt->work( nout_exp, win, wout);

        // apply frequency offset
        gr_complex phase(1,0), inc(std::exp(gr_complex(0, 2*M_PI*offset/d_fs)));
        for(int i=0; i<output.size(); i++){
            output[i] *= phase;
            phase *= inc;
        }
        //volk_32fc_s32fc_x2_rotator_32fc( &output[0], &output[0], inc, &phase, output.size());

        // put into new pdu and send
        pmt::pmt_t newvec = pmt::init_c32vector(nout, &output[0]);
        msg = pmt::cons( meta, newvec );
        message_port_pub(pmt::mp("cpdus"), msg);
    }

    int
    interp_cpdu_impl::work(int noutput_items,
			  gr_vector_const_void_star &input_items,
			  gr_vector_void_star &output_items)
    {
        return noutput_items;
    }

  } /* namespace burst */
} /* namespace gr */

