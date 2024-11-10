cmake_minimum_required(VERSION 3.25...3.29)

# Set the project name, version, and specify the C++ language
project(Spglib_MATLAB
        VERSION 1.0
        LANGUAGES CXX
)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Set the root path for MATLAB
set(MATLAB_ROOT "/Applications/MATLAB_R2024b.app")
find_package(Matlab REQUIRED)

# Set the root path for the Spglib library
set(SPGLIB_ROOT "../install")
find_package(Spglib REQUIRED)

# Add header file path
include_directories("${SPGLIB_ROOT}/include")

# Set source files
set(SOURCE_FILES _spglib.cpp)

# Use the MEX compiler to generate a MEX file
matlab_add_mex(NAME spglib_mex SRC ${SOURCE_FILES} LINK_TO Spglib::symspg)
target_compile_options(spglib_mex PRIVATE -Wall)

# Set the output properties for the generated MEX file
set_target_properties(spglib_mex PROPERTIES
        OUTPUT_NAME "symspg"
        PREFIX ""
        LIBRARY_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/+spglib"
)

# Ensure the 'install' and 'install/+spglib' directories exist
file(MAKE_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/install")
file(MAKE_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/install/+spglib")

# Define MATLAB .m files to be copied
file(GLOB MATLAB_FILES "${CMAKE_CURRENT_SOURCE_DIR}/*.m")

# Add a custom command to copy the generated MEX file to the install directory after build
add_custom_command(TARGET spglib_mex POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy "$<TARGET_FILE:spglib_mex>" "${CMAKE_CURRENT_SOURCE_DIR}/install/+spglib/"
        COMMENT "Copying MEX file to install directory"
)

# Add a custom command to copy MATLAB .m files to the install directory after build
add_custom_command(TARGET spglib_mex POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy ${MATLAB_FILES} "${CMAKE_CURRENT_SOURCE_DIR}/install/"
        COMMENT "Copying MATLAB .m files to install directory"
)
