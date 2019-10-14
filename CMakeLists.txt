cmake_minimum_required(VERSION 3.14)
project(NumericalAssignmentPS2)

set(CMAKE_CXX_STANDARD 14)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -I/usr/include/python2.7")

set(SOURCE
        erdos_renyi.h
        erdos_renyi.cpp
        grg.h
        grg.cpp
        prob_stuff.h
        prob_stuff.cpp
        matplotlibcpp.h)

add_executable(NumericalAssignmentPS2 ${SOURCE} main.cpp)

find_library(PythonLibs 2.7)
target_include_directories(NumericalAssignmentPS2 PRIVATE /home/sv2019/.local/lib/python2.7/site-packages/numpy/core/include)
target_link_libraries(NumericalAssignmentPS2 /usr/lib/x86_64-linux-gnu/libpython2.7.so)