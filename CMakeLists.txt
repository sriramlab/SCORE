cmake_minimum_required(VERSION 2.8.12)
project(SCORE)

include_directories(include/)
include_directories(/usr/lib64/boost/)

set(CMAKE_CXX_FLAGS "-std=c++0x -O3")
IF (NOT DEFINED SSE_SUPPORT)
    SET(SSE_SUPPORT 0)
ENDIF()
ADD_DEFINITIONS(-DSSE_SUPPORT=${SSE_SUPPORT})

IF (NOT DEFINED DEBUG)
    SET(DEBUG 0)
ENDIF()
ADD_DEFINITIONS(-DDEBUG=${DEBUG})

add_executable(SCORE  src/score.cpp src/genotype.cpp src/storage.cpp)
