cmake_minimum_required(VERSION 3.28)
project(raytrace)

set(CMAKE_CXX_STANDARD 20)

add_executable(raytrace src/main.cpp
        include/vec3.h
        include/color.h
        include/ray.h
        include/hittable.h
        include/sphere.h
        include/hittable_list.h
        include/rtweekend.h
        include/interval.h
        include/camera.h
        include/material.h
)
