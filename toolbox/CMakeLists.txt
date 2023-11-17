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
#find_package( Octave REQUIRED )

message( STATUS "Matlab_ROOT_DIR = ${Matlab_ROOT_DIR}" )
message( STATUS "PROJECT_NAME    = ${PROJECT_NAME}" )

if ( (UNIX OR LINUX) AND NOT APPLE )
  set( CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -static-libstdc++ -static-libgcc -Wl,--no-undefined" )
  set( CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -static-libstdc++ -static-libgcc -Wl,--no-undefined" )
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

if ( APPLE )

  execute_process(
    COMMAND brew --prefix
    RESULT_VARIABLE BREW
    OUTPUT_VARIABLE BREW_PREFIX
    OUTPUT_STRIP_TRAILING_WHITESPACE
  )
  if ( BREW EQUAL 0 AND EXISTS "${BREW_PREFIX}" )
    message( STATUS "Found Brew installed at ${BREW_PREFIX}" )
  else()
    message( FATAL_ERROR "homebrew not found, please install it and retry")
  endif()

  include_directories( ${BREW_PREFIX}/include/coin-or )

  find_library( IPOPT ipopt HINTS ${BREW_PREFIX}/lib )

elseif ( UNIX OR LINUX )
  include_directories( src_ipopt_linux/coin-or )
else()
  include_directories( src_ipopt_win/coin-or )
endif()

utils_final_messages()

matlab_add_mex(
  NAME ipopt
  SRC src/ipopt.cc src/IpoptInterfaceCommon.cc
  LINK_TO ${IPOPT}
  OUTPUT_NAME ipopt
)

install( TARGETS ipopt DESTINATION ${CMAKE_SOURCE_DIR}/lib )
