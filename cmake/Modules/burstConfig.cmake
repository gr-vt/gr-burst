INCLUDE(FindPkgConfig)
PKG_CHECK_MODULES(PC_BURST burst)

FIND_PATH(
    BURST_INCLUDE_DIRS
    NAMES burst/api.h
    HINTS $ENV{BURST_DIR}/include
        ${PC_BURST_INCLUDEDIR}
    PATHS ${CMAKE_INSTALL_PREEFIX}/include
          /usr/local/include
          /usr/include
)

FIND_LIBRARY(
    BURST_LIBRARIES
    NAMES gnuradio-burst
    HINTS $ENV{BURST_DIR}/lib
        ${PC_BURST_LIBDIR}
    PATHS ${CMAKE_INSTALL_PREFIX}/lib
          ${CMAKE_INSTALL_PREFIX}/lib64
          /usr/local/lib
          /usr/local/lib64
          /usr/lib
          /usr/lib64
)

INCLUDE(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(BURST DEFAULT_MSG BURST_LIBRARIES BURST_INCLUDE_DIRS)
MARK_AS_ADVANCED(BURST_LIBRARIES BURST_INCLUDE_DIRS)

