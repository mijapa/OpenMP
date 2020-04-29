cmake_minimum_required(VERSION 3.15)
project(BucketSort)

find_package(OpenMP REQUIRED)

set(CMAKE_CXX_STANDARD 14)

add_executable(BucketSort main.c)
target_compile_options(BucketSort PRIVATE -Wall ${OpenMP_CXX_FLAGS})
target_link_libraries(BucketSort PRIVATE ${OpenMP_CXX_FLAGS})