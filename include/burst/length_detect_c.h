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


#ifndef INCLUDED_BURST_LENGTH_DETECT_C_H
#define INCLUDED_BURST_LENGTH_DETECT_C_H

#include <burst/api.h>
#include <gnuradio/sync_block.h>

namespace gr {
  namespace burst {

    /*!
     * \brief <+description of block+>
     * \ingroup burst
     *
     */
    class BURST_API length_detect_c : virtual public gr::sync_block
    {
     public:
      typedef boost::shared_ptr<length_detect_c> sptr;

      /*!
       * \brief Return a shared_ptr to a new instance of burst::length_detect_c.
       *
       * To avoid accidental use of raw pointers, burst::length_detect_c's
       * constructor is in a private implementation
       * class. burst::length_detect_c::make is the public interface for
       * creating new instances.
       */
      static sptr make();
    };

  } // namespace burst
} // namespace gr

#endif /* INCLUDED_BURST_LENGTH_DETECT_C_H */

