#!/bin/bash

rm -r CMakeFiles
rm CMakeCache.txt
rm cmake_install.cmake
rm Makefile
cmake -DCMAKE_BUILD_TYPE=Release CMakeLists.txt && cmake --build . --parallel 8
# mv bloodflow ../cmake
