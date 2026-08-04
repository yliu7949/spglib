cmake_minimum_required(VERSION 3.25...3.29)
project(Spglib_MATLAB_package_test LANGUAGES CXX)

find_package(Spglib CONFIG REQUIRED)

set(SPGLIB_MATLAB_WITH_TESTS ON CACHE BOOL "" FORCE)
enable_testing()
add_subdirectory("${CMAKE_CURRENT_LIST_DIR}/../../../matlab" matlab)
