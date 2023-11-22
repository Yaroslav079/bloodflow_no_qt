cmake_minimum_required(VERSION 3.0)

set(PROJECT bloodflow)
project(${PROJECT})

include_directories(lib)
include_directories(src)
include_directories(ext)

find_package (Eigen3 3.4.0 REQUIRED)
find_package(OpenMP REQUIRED)
find_library(MATH_LIBRARY m)

set (HEADERS
    lib/bayes_runner.h
    lib/calculator.h
    lib/csv_reader.h
    lib/detailed_run.h
    lib/edge.h
    lib/graph.h
    lib/heart_advanced_valves.h
    lib/heart_reg.h
    lib/heart_part.h
    lib/matrix_utils.h
    lib/rcrwindkessel.h
    lib/rescaler.h
    lib/task.h
    lib/ukf.h
    lib/vertex.h
    ext/json.hpp
)

set (SOURCES
    src/bayes_runner.cpp
    src/calculator.cpp
    src/edge.cpp
    src/graph.cpp
    src/heart_part.cpp
    src/main.cpp
    src/matrix_utils.cpp
    src/task.cpp
    src/ukf.cpp
)

add_executable (${PROJECT} ${HEADERS} ${SOURCES})

include(CheckCXXCompilerFlag)
unset(COMPILER_SUPPORTS_MARCH_NATIVE CACHE)
CHECK_CXX_COMPILER_FLAG(-march=native COMPILER_SUPPORTS_MARCH_NATIVE)
if(COMPILER_SUPPORTS_MARCH_NATIVE)
    target_compile_options(${PROJECT_NAME} PUBLIC -march=native)
endif()

target_include_directories(${PROJECT} PUBLIC "${OpenMP_CXX_INCLUDE_DIRS}")
target_link_libraries(${PROJECT} PUBLIC Eigen3::Eigen)
target_link_libraries(${PROJECT} PUBLIC "${OpenMP_CXX_FLAGS}")
target_compile_options(${PROJECT} PUBLIC "${OpenMP_CXX_FLAGS}")