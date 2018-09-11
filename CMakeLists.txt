cmake_minimum_required(VERSION 3.9)
project(kaneko)

set(kaneko_VERSION_MAJOR 1)
set(kaneko_VERSION_MINOR 0)

set(CMAKE_CXX_STANDARD 11)

set(CMAKE_CXX_COMPILER /usr/local/Cellar/gcc/7.1.0/bin/gcc-7)
set(CMAKE_CXX_FLAGS "-lstdc++ -o3")

configure_file (
        "src/KanekoConfig.h.in"
        "${PROJECT_BINARY_DIR}/KanekoConfig.h"
)

include_directories("${PROJECT_BINARY_DIR}/headers")

add_executable(kaneko src/main.cpp src/bchCoder.cpp src/Decoder.cpp src/KanekoKernelProcessor.cpp src/dataForPlot.cpp
                        headers/bchCoder.h headers/dataForPlot.h headers/Decoder.h headers/KanekoKernelProcessor.h)