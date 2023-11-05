cmake_minimum_required(VERSION 3.26)
project(tsp)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O0 -std=gnu++17 -static -lrt -Wl,--whole-archive -lpthread -Wl,--no-whole-archive")
#set(CMAKE_VERBOSE_MAKEFILE ON)

add_executable(tsp main.cpp
        util.cpp
        util.h
        lin_kernighan.cpp
        lin_kernighan.h)
