cmake_minimum_required(VERSION 3.14)
project(random-graphs-models)

set(CMAKE_CXX_STANDARD 14)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -I/usr/include/python2.7")

set(SOURCE
        conf_model.h
        conf_model.cpp
        erdos_renyi.h
        erdos_renyi.cpp
        grg.h
        grg.cpp
        prob_stuff.h
        prob_stuff.cpp
        matplotlibcpp.h timing.h timing.cpp)

add_executable(random-graphs-models ${SOURCE} main.cpp)

find_library(PythonLibs 2.7)
target_include_directories(random-graphs-models PRIVATE /home/sv2019/.local/lib/python2.7/site-packages/numpy/core/include)
target_link_libraries(random-graphs-models /usr/lib/x86_64-linux-gnu/libpython2.7.so)

find_package(Boost COMPONENTS system filesystem REQUIRED)
target_link_libraries(random-graphs-models ${Boost_FILESYSTEM_LIBRARY} ${Boost_SYSTEM_LIBRARY})