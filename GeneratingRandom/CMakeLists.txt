cmake_minimum_required(VERSION 3.15)
project(GeneratingRandom)

find_package(OpenMP REQUIRED)

set(CMAKE_CXX_STANDARD 14)

add_executable(GeneratingRandom main.c)
target_compile_options(GeneratingRandom PRIVATE -Wall ${OpenMP_CXX_FLAGS})
target_link_libraries(GeneratingRandom PRIVATE ${OpenMP_CXX_FLAGS})
