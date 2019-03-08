# CMakeLists.txt for event package. It creates a library with dictionary and a main program
cmake_minimum_required(VERSION 3.0 FATAL_ERROR)
project(clustering)
set(CMAKE_MACOSX_RPATH 1)

# You need to tell CMake where to find the ROOT installation. This can be done in a number of ways:
#   - ROOT built with classic configure/make use the provided $ROOTSYS/etc/cmake/FindROOT.cmake
#   - ROOT built with CMake. Add in CMAKE_PREFIX_PATH the installation prefix for ROOT
list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS})

#---Locate the ROOT package and defines a number of variables (e.g. ROOT_INCLUDE_DIRS)
set(CMAKE_MODULE_PATH $ENV{ROOTSYS}/etc/cmake)
find_package(ROOT)

#---Define useful ROOT functions and macros (e.g. ROOT_GENERATE_DICTIONARY)
 include_directories(${PROJECT_SOURCE_DIR}/
                    ${ROOT_INCLUDE_DIR}/)

#----------------------------------------------------------------------------
# Locate sources and headers for this project
# NB: headers are included so they will show up in IDEs
#
include(${ROOT_USE_FILE})

file(GLOB sources ${PROJECT_SOURCE_DIR}/*.C)
file(GLOB headers ${PROJECT_SOURCE_DIR}/*.h)

#---Create  a main program using the library
add_executable(clustering Clustering.cpp ${sources} ${headers})
target_link_libraries(clustering ${ROOT_LIBRARIES})

install(TARGETS clustering DESTINATION bin)
