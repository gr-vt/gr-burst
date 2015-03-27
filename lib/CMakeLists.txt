# Copyright 2011,2012 Free Software Foundation, Inc.
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
# Setup library
########################################################################
include(GrPlatform) #define LIB_SUFFIX

include_directories(${Boost_INCLUDE_DIR})
link_directories(${Boost_LIBRARY_DIRS})
list(APPEND burst_sources
    interp_cpdu_impl.cc
    correlator_filter_impl.cc
    length_detect_c_impl.cc
    synchronizer_v4_impl.cc
    qa_helpers.cc )

add_library(gnuradio-burst SHARED ${burst_sources})
target_link_libraries(gnuradio-burst ${Boost_LIBRARIES} ${GNURADIO_ALL_LIBRARIES} ${GNURADIO_EVENTSTREAM_LIBRARIES} ${GSL_LIBRARIES} ${CPPUNIT_LIBRARIES} ${MAPPER_LIBRARIES})
set_target_properties(gnuradio-burst PROPERTIES DEFINE_SYMBOL "gnuradio_burst_EXPORTS")

########################################################################
# Install built library files
########################################################################
install(TARGETS gnuradio-burst
    LIBRARY DESTINATION lib${LIB_SUFFIX} # .so/.dylib file
    ARCHIVE DESTINATION lib${LIB_SUFFIX} # .lib file
    RUNTIME DESTINATION bin              # .dll file
)

########################################################################
# Build and register unit test
########################################################################
include(GrTest)

include_directories(${CPPUNIT_INCLUDE_DIRS})

list(APPEND test_burst_sources
    ${CMAKE_CURRENT_SOURCE_DIR}/test_burst.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/qa_synchronizer_v4.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/qa_burst.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/qa_helpers.cc
)

add_executable(test-burst ${test_burst_sources})

target_link_libraries(
  test-burst
  ${GNURADIO_RUNTIME_LIBRARIES}
  ${Boost_LIBRARIES}
  ${CPPUNIT_LIBRARIES}
  ${MAPPER_LIBRARIES}
  gnuradio-burst
)

GR_ADD_TEST(test_burst test-burst)
