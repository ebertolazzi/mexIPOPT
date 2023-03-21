#-----------------------------------------------------------------------#
#  file: CMakeLists.txt                                                 #
#                                                                       #
#  version: 1.0  date March, 21, 2023                                   #
#                                                                       #
#  Copyright (C) 2023                                                   #
#                                                                       #
#      Enrico Bertolazzi                                                #
#      Dipartimento di Ingegneria Industriale                           #
#      Universita` degli Studi di Trento                                #
#      Via Sommarive 9, I-38123, Trento, Italy                          #
#      email: enrico.bertolazzi@unitn.it                                #
#-----------------------------------------------------------------------#

cmake_minimum_required( VERSION 3.14 )

project( MexIPOPT )

include(./cmake_utils/CMakeLists-common.txt)

project( ${PROJECT_NAME} VERSION ${UTILS_PROJECT_VERSION} )

if ( NOT DEFINED CMAKE_BUILD_TYPE )
  set( CMAKE_BUILD_TYPE Release )
endif()

find_package( Matlab REQUIRED )

message( STATUS "Matlab_ROOT_DIR = ${Matlab_ROOT_DIR}" )
message( STATUS "PROJECT_NAME    = ${PROJECT_NAME}" )

if ( (UNIX OR LINUX) AND NOT APPLE )
  set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -static-libstdc++ -static-libgcc -Wl,--no-undefined")
  set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -static-libstdc++ -static-libgcc -Wl,--no-undefined")
endif()

include(./cmake_utils/CMakeLists-cflags.txt)
include(./cmake_utils/CMakeLists-utilities.txt)

# Evaluating the suffix to append to target built by this CMake file
utils_artifacts_suffix(ARTIFACTS_STATIC_SUFFIX TRUE)
if( UTILS_BUILD_SHARED )
  utils_artifacts_suffix(ARTIFACTS_DYNAMIC_SUFFIX FALSE)
endif()

message( STATUS "Compiler used: ${CMAKE_CXX_COMPILER_ID}" )
message( STATUS "BASE SUFFIXES = ${ARTIFACTS_STATIC_SUFFIX} ${ARTIFACTS_DYNAMIC_SUFFIX}" )

# extra include directory for compilation
include_directories( src )

if ( (UNIX OR LINUX) AND NOT APPLE )
  include_directories( src_ipopt_linux/coin-or )
  set( MEXNAME ipopt_linux )
elseif ( APPLE )
  include_directories( src_ipopt_osx/coin-or )
  set(CMAKE_MACOSX_RPATH 1)
  set(CMAKE_INSTALL_RPATH ${CMAKE_CURRENT_SOURCE_DIR}/bin/osx )
  set( MEXNAME ipopt_osx )
else()
  include_directories( src_ipopt_win/coin-or )
  set( MEXNAME ipopt_win )
endif()

utils_final_messages()

matlab_add_mex(
  NAME ${MEXNAME}
  SRC src/ipopt.cc src/IpoptInterfaceCommon.cc
  LINK_TO -L${CMAKE_CURRENT_SOURCE_DIR}/bin/osx
          -lblas
          -ldmumps
          -lgcc_s.1
          -lgcc_s.1.1
          -lgfortran.5
          -lipopt.3
          -lipoptamplinterface.3
          -llapack
          -lmpiseq
          -lmumps_common
          -lopenblas.0
          -lpord
          -lquadmath.0
          -lsipopt.3
          -lstdc++.6
  OUTPUT_NAME ${MEXNAME}
)
