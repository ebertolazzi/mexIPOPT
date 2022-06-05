############################################################################
#                                                                          #
#  file: CMakeLists.txt                                                    #
#  Copyright (C) 2022                                                      #
#      Enrico Bertolazzi                                                   #
#      Dipartimento di Ingegneria Industriale                              #
#      Universita` degli Studi di Trento                                   #
#      email: enrico.bertolazzi@ing.unitn.it                               #
#                                                                          #
############################################################################

cmake_minimum_required(VERSION 3.14)

include(./cmake/CMakeLists-common.txt)

set( HSL_IPOPT_PROJECT_VERSION 1.0 )
project( HSL_FOR_IPOPT VERSION ${HSL_IPOPT_PROJECT_VERSION} )

enable_language(Fortran)

if ( HSL_IPOPT_INSTALL_LOCAL )
  set( CMAKE_INSTALL_PREFIX ${CMAKE_HOME_DIRECTORY} )
endif()

message("HSL_FOR_IPOPT current version = v${HSL_IPOPT_PROJECT_VERSION}")

include(./cmake/CMakeLists-cflags.txt)
include(./cmake/CMakeLists-utilities.txt)

# Evaluating the suffix to append to target built by this CMake file
hsl_for_ipopt_artifacts_suffix(ARTIFACTS_STATIC_SUFFIX TRUE)
if( HSL_IPOPT_BUILD_SHARED )
  hsl_for_ipopt_artifacts_suffix(ARTIFACTS_DYNAMIC_SUFFIX FALSE)
endif()

message( STATUS "Compiler used: ${CMAKE_CXX_COMPILER_ID}" )
message( STATUS "BASE SUFFIXES = ${ARTIFACTS_STATIC_SUFFIX} ${ARTIFACTS_DYNAMIC_SUFFIX}" )

#   ___
#  / __| ___ _  _ _ _ __ ___ ___
#  \__ \/ _ \ || | '_/ _/ -_|_-<
#  |___/\___/\_,_|_| \__\___/__/
#
set( SOURCES )
file(
  GLOB S
  ${CMAKE_CURRENT_SOURCE_DIR}/*/*.f
  ${CMAKE_CURRENT_SOURCE_DIR}/*/*/*.f
  ${CMAKE_CURRENT_SOURCE_DIR}/*/*.f90
  ${CMAKE_CURRENT_SOURCE_DIR}/*/*/*.f90
)
foreach(F ${S})
  # build relative path respect to source directory
  file( RELATIVE_PATH RF ${CMAKE_CURRENT_SOURCE_DIR} "${F}" )
  list( APPEND SOURCES ${RF} )
endforeach(F ${S})
#list(FILTER SOURCES EXCLUDE REGEX ".*hsl_fake.cc$")

#   _  _             _
#  | || |___ __ _ __| |___ _ _ ___
#  | __ / -_) _` / _` / -_) '_(_-<
#  |_||_\___\__,_\__,_\___|_| /__/
#
set( PUBLIC_HEADERS )
file(
  GLOB S
  ${CMAKE_CURRENT_SOURCE_DIR}/*/C/*.h
)
foreach(F ${S})
  # build relative path respect to source directory
  file( RELATIVE_PATH RF ${CMAKE_CURRENT_SOURCE_DIR} "${F}" )
  list( APPEND PUBLIC_HEADERS ${RF} )
endforeach(F ${S})
set( PRIVATE_HEADERS )
# no private headers

#list(FILTER PUBLIC_HEADERS EXCLUDE REGEX ".*.tmpl$")
#list(FILTER PRIVATE_HEADERS EXCLUDE REGEX ".*.tmpl$")

# extra include directory for compilation
#include_directories( ./lib3rd/include )

#       _                           _                 _
#    __| | ___ _ __   ___ _ __   __| | ___ _ __   ___(_) ___  ___
#   / _` |/ _ \ '_ \ / _ \ '_ \ / _` |/ _ \ '_ \ / __| |/ _ \/ __|
#  | (_| |  __/ |_) |  __/ | | | (_| |  __/ | | | (__| |  __/\__ \
#   \__,_|\___| .__/ \___|_| |_|\__,_|\___|_| |_|\___|_|\___||___/
#             |_|
#
set(DEPEND_TARGETS)

#   _____                  _
#  |_   _|_ _ _ _ __ _ ___| |_ ___
#    | |/ _` | '_/ _` / -_)  _(_-<
#    |_|\__,_|_| \__, \___|\__/__/
#                |___/
#
# define target
hsl_for_ipopt_setup_target(
  ${PROJECT_NAME}
  TRUE # build static
  "${CMAKE_CURRENT_SOURCE_DIR}"
  "${SOURCES}"
  "${PUBLIC_HEADERS}"
  "${PRIVATE_HEADERS}"
  TRUE # flatten headers
)
if ( HSL_IPOPT_BUILD_SHARED )
  hsl_for_ipopt_setup_target(
    ${PROJECT_NAME}
    FALSE # build dynamic
    "${CMAKE_CURRENT_SOURCE_DIR}"
    "${SOURCES}"
    "${PUBLIC_HEADERS}"
    "${PRIVATE_HEADERS}"
    TRUE # flatten headers
  )
  #target_link_libraries( "${HSL_IPOPT_NAMESPACE}_${PROJECT_NAME}" ${lapackblas_libraries} )
  target_link_libraries( "${HSL_IPOPT_NAMESPACE}_${PROJECT_NAME}" "-framework Accelerate -L/usr/local/Cellar/metis/5.1.0/lib -lmetis" )
endif()

#   ___         _        _ _
#  |_ _|_ _  __| |_ __ _| | |
#   | || ' \(_-<  _/ _` | | |
#  |___|_||_/__/\__\__,_|_|_|
#
install(
  TARGETS ${TARGETS_TO_BE_INSTALLED}
  # EXPORT (da definire)
  LIBRARY DESTINATION ${INSTALL_PREFIX}/dll
  ARCHIVE DESTINATION ${INSTALL_PREFIX}/lib
  RUNTIME DESTINATION ${INSTALL_PREFIX}/bin
  PUBLIC_HEADER DESTINATION ${PROJECT_PUBLIC_HEADER_DESTINATION}
)

#   ___         _
#  | _ \__ _ __| |____ _ __ _ ___
#  |  _/ _` / _| / / _` / _` / -_)
#  |_| \__,_\__|_\_\__,_\__, \___|
#                       |___/
#

# do not work
#include( ./cmake/CMakeLists-cpack.txt )

hsl_for_ipopt_final_messages()
