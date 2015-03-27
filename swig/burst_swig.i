/* -*- c++ -*- */

#define BURST_API

%include "gnuradio.i"			// the common stuff
%include "mapper/mapper/swig/mapper_swig.i"

//load generated python docstrings
%include "burst_swig_doc.i"

%{
#include "burst/interp_cpdu.h"
#include "burst/correlator_filter.h"
#include "burst/length_detect_c.h"
#include "burst/synchronizer_v4.h"
%}




%include "burst/interp_cpdu.h"
GR_SWIG_BLOCK_MAGIC2(burst, interp_cpdu);
%include "burst/correlator_filter.h"
GR_SWIG_BLOCK_MAGIC2(burst, correlator_filter);
%include "burst/length_detect_c.h"
GR_SWIG_BLOCK_MAGIC2(burst, length_detect_c);
%include "burst/synchronizer_v4.h"
GR_SWIG_BLOCK_MAGIC2(burst, synchronizer_v4);
