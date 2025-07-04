#!/bin/bash
mkdir extern/pybind11/build
cd extern/pybind11/build
cmake ..
make
cd ../../../
mkdir girgs/build
cd girgs/build
cmake ..
make
cd ../..
mkdir build
cd build
cmake ..
make