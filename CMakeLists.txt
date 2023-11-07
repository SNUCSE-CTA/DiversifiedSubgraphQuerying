cmake_minimum_required(VERSION 3.24)
project(DivDAF)

# Replace '/usr/local/' with the directory where you installed SDSL, if different
set(SDSL_INCLUDE_DIRS "/usr/local/include")
set(SDSL_LIBRARIES "/usr/local/lib")

# Add the SDSL headers
include_directories(${SDSL_INCLUDE_DIRS})

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++20 -O3 -flto -pthread -w -DNDEBUG -DDEGORDER -DFIRSTASC -DSECONDASC")
add_executable(DivDAF DAF_div.cpp)
