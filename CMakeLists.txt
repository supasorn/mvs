cmake_minimum_required(VERSION 2.8 FATAL_ERROR)

project(research)

FIND_PACKAGE(OpenCV REQUIRED)
SET(EXTRA_LIBS GL glut GLU /usr/local/lib/libgflags.a)
SET(OPTICALFLOW /home/supasorn/research/OpticalFlow/mex/CVOpticalFlow.a)
SET(CMAKE_CXX_LINK_FLAGS "${CMAKE_CXX_LINK_FLAGS} -lpthread -pthread")
#link_directories(/usr/local/lib)
add_definitions(-fopenmp)
add_definitions(-D__UBUNTU__)
add_definitions(-std=c++0x)
#FIND_PACKAGE(CUDA)

SET(CUDA_NVCC_FLAGS "-arch=sm_20")
SET(OSBASE "/home/supasorn")
message(${CMAKE_SYSTEM_VERSION})
message(${CMAKE_CURRENT_SOURCE_DIR})

FIND_PACKAGE(Ceres 1.8 REQUIRED)
FIND_PACKAGE(OpenGL REQUIRED )


message("OPENCV: ${OpenCV_VERSION}")
#link_directories(/usr/lib64 ${PCL_LIBRARY_DIRS})
include_directories(${INCLUDE_DIRECTORIES} ../research/OpticalFlow/mex ${CERES_INCLUDE_DIRS} ${OPENCV_INCLUDE_DIRS} ../research/utils/ ../research/) 

add_definitions(-D_LINUX_MAC -D_OPENCV -w)

#set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} ")
if (NOT CMAKE_BUILD_TYPE)
  message(STATUS "No build type selected, default to Release")
  set(CMAKE_BUILD_TYPE "Release")
endif()

message("Build type is ${CMAKE_BUILD_TYPE}")
message("CMAKE_CXX_FLAGS_DEBUG is ${CMAKE_CXX_FLAGS_DEBUG}")
message("CMAKE_CXX_FLAGS_RELEASE is ${CMAKE_CXX_FLAGS_RELEASE}")

add_library (pyramid ../research/utils/pyramid.cpp)

add_library (utilities ../research/utils/utilities.cpp)
target_link_libraries(utilities pyramid ${OpenCV_LIBS} ${OPTICALFLOW} ${EXTRA_LIBS})

add_library (libviz ../research/viz.cpp)

add_library (raster ../research/utils/raster.cpp)

add_executable (featuredetect featuredetect.cpp)
target_link_libraries (featuredetect ceres utilities libviz ply.so) 
