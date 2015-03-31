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


#ifndef INCLUDED_BURST_SYNCHRONIZER_V4_H
#define INCLUDED_BURST_SYNCHRONIZER_V4_H

#include <burst/api.h>
#include <gnuradio/sync_block.h>

namespace gr {
  namespace burst {

    /*!
     * \brief <+description of block+>
     * \ingroup burst
     *
     */
    class BURST_API synchronizer_v4 : virtual public gr::sync_block
    {
     public:
      typedef boost::shared_ptr<synchronizer_v4> sptr;

      /*!
       * \brief Return a shared_ptr to a new instance of burst::synchronizer_v4.
       *
       * To avoid accidental use of raw pointers, burst::synchronizer_v4's
       * constructor is in a private implementation
       * class. burst::synchronizer_v4::make is the public interface for
       * creating new instances.
       */
      static sptr make(double Fs, int sps, std::vector<unsigned char> preamble_bits, std::vector<int> sym_mapping);
      virtual void enableDebugMode() = 0;
      virtual void handler(pmt::pmt_t msg) = 0;
    };

  } // namespace burst
} // namespace gr

#endif /* INCLUDED_BURST_SYNCHRONIZER_V4_H */

