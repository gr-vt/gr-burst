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


#ifndef _QA_SYNCHRONIZER_V4_H_
#define _QA_SYNCHRONIZER_V4_H_

#include <cppunit/extensions/HelperMacros.h>
#include <cppunit/TestCase.h>

#include <gnuradio/gr_complex.h>		// to define gr_complex
#include <vector>

namespace gr {
  namespace burst {

  	class qa_synchronizer_v4 : public CppUnit::TestCase
    {
    public:
      CPPUNIT_TEST_SUITE(qa_synchronizer_v4);
      CPPUNIT_TEST(t1);
      CPPUNIT_TEST_SUITE_END();

    private:
      void t1();
    };

  } /* namespace wifi_detector */
} /* namespace gr */

#endif /* _QA_ES_80211B_HELPERS_H_ */

