cmake_minimum_required(VERSION 3.9)
project(GridQuery)

set(CMAKE_CXX_STANDARD 11)

add_executable(GridQuery
        main.cpp
        src/ClassType.cpp
        src/CodeBeauty.cpp
        src/ComAlgorithm.cpp
        ine_src/INEClassType.cpp
        ine_src/INE.cpp
        ine_src/INEBeautyCode.cpp)

include_directories(./include)
include_directories(./ine_include)