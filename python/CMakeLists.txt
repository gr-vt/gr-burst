# Copyright 2011 Free Software Foundation, Inc.
#
# This file is part of GNU Radio
#
# GNU Radio is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 3, or (at your option)
# any later version.
#
# GNU Radio is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with GNU Radio; see the file COPYING.  If not, write to
# the Free Software Foundation, Inc., 51 Franklin Street,
# Boston, MA 02110-1301, USA.

########################################################################
# Include python install macros
########################################################################
include(GrPython)
if(NOT PYTHONINTERP_FOUND)
    return()
endif()

########################################################################
# Install python sources
########################################################################
GR_PYTHON_INSTALL(
    FILES
    __init__.py
    arq.py
    padder.py
    random_drop.py
    burst_scheduler.py
    burst_scheduler2.py
    framer.py
    slicer.py
    tofpdu.py
    deframer.py
    framealign.py
    fsk_time_sync.py
    randomizer.py
    preamble_insert.py
    preamble_correlator.py
    pdulib.py
    unpacker.py
    packer.py
    length_detect.py
    cpdu_matlab_writer.py
    #synchronizer.py
    #synchronizer_v2.py
    synchronizer_v3.py
    DESTINATION ${GR_PYTHON_DIR}/burst
)

########################################################################
# Handle the unit tests
########################################################################
include(GrTest)

set(GR_TEST_TARGET_DEPS gnuradio-burst)
set(GR_TEST_PYTHON_DIRS ${CMAKE_BINARY_DIR}/swig)
#GR_ADD_TEST(qa_gr_messages ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/qa_synchronizer_v2.py)
GR_ADD_TEST(qa_gr_messages ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/qa_synchronizer_v3.py)
