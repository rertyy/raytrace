cmake_minimum_required(VERSION 3.28)
project(raytrace)

set(CMAKE_CXX_STANDARD 20)

add_executable(raytrace main.cpp
        vec3.h
        color.h
        ray.h)
