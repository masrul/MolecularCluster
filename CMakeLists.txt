######################################################################
# @author      : Masrul Huda (mail2masrul@gmail.com)
# @file        : CMakeLists
# @created     : Thursday Jun 30, 2022 19:52:24 CDT
######################################################################

cmake_minimum_required(VERSION 3.0)
set(CMAKE_C_COMPILER /opt/local/bin/gcc-mp-11)
set(CMAKE_CXX_COMPILER /opt/local/bin/g++-mp-11)

project("cluster")
set(CMAKE_VERBOSE_MAKEFILE ON)

# Optimization flag 
if ("${CMAKE_CXX_COMPILER_ID}" MATCHES "Clang")
    set(CMAKE_CXX_FLAGS "-O2 -Wall  -std=c++11")
    set(CMAKE_C_FLAGS "-O2 -Wno-incompatible-function-pointer-types")
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
    set(CMAKE_CXX_FLAGS " -fopenmp -O2 -Wno-incompatible-function-pointer-types")
    set(CMAKE_C_FLAGS " -fopenmp -O2")
endif()


set(ENV{PKG_CONFIG_PATH} "$ENV{PKG_CONFIG_PATH}:/Users/MasrulHuda/Software/libgmxtraj/lib/pkgconfig")
include(FindPkgConfig)
pkg_search_module(GMXTRAJ REQUIRED libgmxtraj)

include_directories(${GMXTRAJ_INCLUDE_DIRS})
link_directories (${GMXTRAJ_LIBRARY_DIRS})
add_executable(${PROJECT_NAME} src/cluster.cpp)
target_link_libraries(${PROJECT_NAME} ${GMXTRAJ_LIBRARIES})
