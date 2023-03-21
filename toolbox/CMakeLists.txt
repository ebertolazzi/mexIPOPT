#-----------------------------------------------------------------------#
#  file: CMakeLists.txt                                                 #
#                                                                       #
#  version: 1.0   date 9/3/2022                                         #
#                                                                       #
#  Copyright (C) 2022                                                   #
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
include_directories( src src_ipopt_osx/coin-or )

message( STATUS "SOURCES = ${SOURCES}" )

#   _____                  _
#  |_   _|_ _ _ _ __ _ ___| |_ ___
#    | |/ _` | '_/ _` / -_)  _(_-<
#    |_|\__,_|_| \__, \___|\__/__/
#                |___/
#

utils_final_messages()

message( STATUS "CMAKE_BUILD_TYPE = ${CMAKE_BUILD_TYPE}")
if ( CMAKE_BUILD_TYPE MATCHES Debug )
  message( STATUS "DEBUG MODE\n\n")
  add_compile_definitions(DEBUG)
endif()

set(CMAKE_MACOSX_RPATH 1)
set(CMAKE_INSTALL_RPATH ${CMAKE_CURRENT_SOURCE_DIR}/bin/osx )

matlab_add_mex(
  NAME ipopt_osx
  SRC src/ipopt.cc src/IpoptInterfaceCommon.cc
  LINK_TO -L${CMAKE_CURRENT_SOURCE_DIR}/bin/osx -lblas -ldmumps
  ${CMAKE_CURRENT_SOURCE_DIR}/bin/osx/libgcc_s.1.dylib
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
OUTPUT_NAME ipopt_osx
)
