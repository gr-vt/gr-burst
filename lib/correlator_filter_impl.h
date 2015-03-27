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

#ifndef INCLUDED_BURST_THRESHOLD_TRIGGER_IMPL_H
#define INCLUDED_BURST_THRESHOLD_TRIGGER_IMPL_H

#include <burst/correlator_filter.h>

namespace gr {
  namespace burst {

    class correlator_filter_impl : public correlator_filter
    {
     private:
      float d_ma_mu;
      float d_ma_sigma;
      int d_l1;
      float d_i1;
      float d_thresh;
      uint64_t d_last_event;

     public:
      correlator_filter_impl(int l1, float i1, float thresh);
      ~correlator_filter_impl();

      // Where all the action really happens
      int work(int noutput_items,
	       gr_vector_const_void_star &input_items,
	       gr_vector_void_star &output_items);
    };

  } // namespace burst
} // namespace gr

#endif /* INCLUDED_BURST_THRESHOLD_TRIGGER_IMPL_H */

