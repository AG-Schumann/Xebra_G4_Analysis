# CMakeLists.txt for event package. It creates a library with dictionary and a main program
cmake_minimum_required(VERSION 3.0 FATAL_ERROR)
project(XebraMCProcessor)
set(CMAKE_MACOSX_RPATH 1)

# You need to tell CMake where to find the ROOT installation. This can be done in a number of ways:
#   - ROOT built with classic configure/make use the provided $ROOTSYS/etc/cmake/FindROOT.cmake
#   - ROOT built with CMake. Add in CMAKE_PREFIX_PATH the installation prefix for ROOT
list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS})

#---Locate the ROOT package and defines a number of variables (e.g. ROOT_INCLUDE_DIRS)
set(CMAKE_MODULE_PATH $ENV{ROOTSYS}/etc/cmake)
find_package(ROOT)
#find_package(NEST REQUIRED)
#include_directories(${NEST_INCLUDE_DIRS})

#---Define useful ROOT functions and macros (e.g. ROOT_GENERATE_DICTIONARY)
include_directories(${PROJECT_SOURCE_DIR}/
                    ${ROOT_INCLUDE_DIR}/)

#----------------------------------------------------------------------------
# Locate sources and headers for this project
# NB: headers are included so they will show up in IDEs
#
include(${ROOT_USE_FILE})
include_directories(include)
file(GLOB sources "src/*.cc")
file(GLOB headers "include/*.h")
#file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc)
#file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.hh)



#---Create  a main program using the library
add_executable(XebraMCProcessor XebraMCProcessor.cpp ${sources} ${headers})
target_link_libraries(XebraMCProcessor ${ROOT_LIBRARIES})

install(TARGETS XebraMCProcessor DESTINATION bin)
