cmake_minimum_required(VERSION 2.8.12)
project(SCORE)

include_directories(include/)
FIND_PACKAGE( Boost COMPONENTS program_options REQUIRED )
INCLUDE_DIRECTORIES( ${Boost_INCLUDE_DIR} )

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
TARGET_LINK_LIBRARIES( SCORE LINK_PUBLIC ${Boost_LIBRARIES} )