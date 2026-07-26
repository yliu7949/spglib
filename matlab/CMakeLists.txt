cmake_minimum_required(VERSION 3.25...3.29)

# Set the project name, version, and specify the C++ language
project(Spglib_MATLAB
        VERSION 1.0
        LANGUAGES CXX
)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Set the full MATLAB installation used to build the MEX file. Override this
# cache entry with -DMatlab_ROOT_DIR=... when building another MATLAB release.
set(Matlab_ROOT_DIR "/Applications/MATLAB_R2026b.app"
        CACHE PATH "MATLAB installation root")
find_package(Matlab REQUIRED)

# _spglib.cpp is C++ code that implements the classic C mexFunction API.
# FindMatlab otherwise treats every C++ source as a C++ MEX API target on macOS
# and injects unresolved mexFunctionAdapter symbols.
set(Matlab_HAS_CPP_API FALSE)

# Set the root path for the Spglib library
set(SPGLIB_ROOT "../install")
set(Spglib_DIR "../install/lib/cmake/Spglib")
find_package(Spglib REQUIRED)

# Add header file path
include_directories("${SPGLIB_ROOT}/include")

# Set source files
set(SOURCE_FILES
    _spglib.cpp
    dataset_converter.cpp
)

# Use only the C MEX API libraries required by this source. FindMatlab may
# otherwise add the C++ Engine and Data Array libraries, which are not needed
# by symspg and are not portable to a MATLAB Runtime deployment.
matlab_add_mex(
        NAME spglib_mex
        SRC ${SOURCE_FILES}
        NO_IMPLICIT_LINK_TO_MATLAB_LIBRARIES
        LINK_TO Matlab::mex Matlab::mx Spglib::symspg
)
target_compile_options(spglib_mex PRIVATE -Wall)

# Set the output properties for the generated MEX file
set_target_properties(spglib_mex PROPERTIES
        OUTPUT_NAME "symspg"
        PREFIX ""
        LIBRARY_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/+spglib"
        SKIP_BUILD_RPATH TRUE
        BUILD_RPATH ""
        INSTALL_RPATH ""
)

# Ensure the 'install' and 'install/+spglib' directories exist
file(MAKE_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/install")
file(MAKE_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/install/+spglib")

# Define MATLAB .m files to be copied
file(GLOB MATLAB_FILES "${CMAKE_CURRENT_SOURCE_DIR}/*.m")

# Add a custom command to copy the generated MEX file and .m files to the install directory after build
add_custom_command(TARGET spglib_mex POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy "$<TARGET_FILE:spglib_mex>" "${CMAKE_CURRENT_SOURCE_DIR}/install/+spglib/"
        COMMAND ${CMAKE_COMMAND} -E copy ${MATLAB_FILES} "${CMAKE_CURRENT_SOURCE_DIR}/install/+spglib/"
        COMMENT "Copying MEX file and MATLAB .m files to install directory"
)

# Allow the option to decide whether to compile _spglib_test.cpp
option(BUILD_TEST_EXECUTABLE "Build test executable" OFF)

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
