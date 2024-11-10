cmake_minimum_required(VERSION 3.25...3.29)

# Set the project name and version, and specify the C++ language
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

# Check if the +spglib folder exists, and create it if it does not
file(MAKE_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/+spglib")

# Set the suffix for the generated file to MATLAB MEX format
set_target_properties(spglib_mex PROPERTIES
        OUTPUT_NAME "symspg"
        PREFIX ""
        LIBRARY_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/+spglib"
)

# Check if the install folder exists, and create it if it does not
file(MAKE_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/install")

# Copy all *.m files from the current directory to the install folder
file(GLOB MATLAB_FILES "${CMAKE_CURRENT_SOURCE_DIR}/*.m")
foreach(file ${MATLAB_FILES})
    configure_file(${file} "${CMAKE_CURRENT_SOURCE_DIR}/install" COPYONLY)
endforeach()

# Copy the +spglib folder and its contents to the install folder
file(COPY "${CMAKE_CURRENT_SOURCE_DIR}/+spglib"
        DESTINATION "${CMAKE_CURRENT_SOURCE_DIR}/install")

# Allow the option to decide whether to compile _spglib_test.cpp
option(BUILD_TEST_EXECUTABLE "Build test executable" ON)

# If the switch is turned on, add a target to compile _spglib_test.cpp
if(BUILD_TEST_EXECUTABLE)
    set(TEST_SOURCE_FILES _spglib_test.cpp)
    add_executable(test_executable ${TEST_SOURCE_FILES})
    set_target_properties(test_executable PROPERTIES
            OUTPUT_NAME "_test"
            PREFIX ""
            RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}"
    )
    target_link_libraries(test_executable PRIVATE Spglib::symspg)
    target_include_directories(test_executable PRIVATE "${SPGLIB_ROOT}/include")
endif()
