cmake_minimum_required(VERSION 3.10)

project(coloringgirgs)

add_subdirectory(extern/pybind11)
pybind11_add_module(coloringgirgs coloringgirgs.cpp)

set(CMAKE_CXX_FLAGS "-O3") 

include_directories("${CMAKE_SOURCE_DIR}/girgs/build/source/girgs/include")
include_directories("${CMAKE_SOURCE_DIR}/girgs/source/girgs/include/girgs")

target_link_libraries(${PROJECT_NAME} PRIVATE "${CMAKE_SOURCE_DIR}/girgs/build/libgirgs.so.1")
