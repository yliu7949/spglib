cmake_minimum_required(VERSION 3.25...3.29)

if (CMAKE_SOURCE_DIR STREQUAL CMAKE_CURRENT_SOURCE_DIR)
    project(Spglib_MATLAB VERSION 1.0 LANGUAGES CXX)
else ()
    enable_language(CXX)
endif ()

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

if (DEFINED SPGLIB_WITH_TESTS)
    set(_Spglib_MATLAB_default_tests ${SPGLIB_WITH_TESTS})
else ()
    set(_Spglib_MATLAB_default_tests ON)
endif ()
option(SPGLIB_MATLAB_WITH_TESTS "Spglib: Run MATLAB unit tests"
       ${_Spglib_MATLAB_default_tests})

# Matlab_ROOT_DIR may be supplied on the command line. Otherwise FindMatlab
# discovers an installed MATLAB release using CMake's standard search logic.
if (SPGLIB_MATLAB_WITH_TESTS)
    find_package(Matlab REQUIRED COMPONENTS MAIN_PROGRAM)
else ()
    find_package(Matlab REQUIRED)
endif ()

# _spglib.cpp is C++ code that implements the classic C mexFunction API.
# FindMatlab otherwise treats every C++ source as a C++ MEX API target on macOS
# and injects unresolved mexFunctionAdapter symbols.
set(Matlab_HAS_CPP_API FALSE)

# Reuse the in-tree target when this directory is built from the main project;
# require an installed Spglib package when it is configured standalone.
if (NOT TARGET Spglib::symspg)
    find_package(Spglib CONFIG REQUIRED)
elseif (TARGET Spglib_symspg)
    # Static libraries linked into a MEX module need position-independent code
    # on platforms such as Linux.
    set_property(TARGET Spglib_symspg PROPERTY POSITION_INDEPENDENT_CODE ON)
endif ()

# Set source files
set(SOURCE_FILES
    _spglib.cpp
    dataset_converter.cpp
)

# MATLAB package layout used by both the source tree and the install bundle.
set(MATLAB_PACKAGE_DIR
        "${CMAKE_CURRENT_SOURCE_DIR}/+kssolv/+analysis/+spglib")
set(SPGLIB_MATLAB_OUTPUT_ROOT "${CMAKE_CURRENT_BINARY_DIR}/install"
        CACHE PATH "Root directory of the generated MATLAB package")
set(MATLAB_INSTALL_PACKAGE_DIR
        "${SPGLIB_MATLAB_OUTPUT_ROOT}/+kssolv/+analysis/+spglib")

# Define the MATLAB source and license files explicitly so that unrelated
# top-level .m files are never included in the install bundle.
set(MATLAB_PACKAGE_FILES
        "${MATLAB_PACKAGE_DIR}/Spglib.m"
        "${MATLAB_PACKAGE_DIR}/SpglibError.m"
        "${MATLAB_PACKAGE_DIR}/Spglib.rights"
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
if (MSVC)
    target_compile_options(spglib_mex PRIVATE /W3)
else ()
    target_compile_options(spglib_mex PRIVATE -Wall)
endif ()

# Set the output properties for the generated MEX file
set_target_properties(spglib_mex PROPERTIES
        OUTPUT_NAME "symspg"
        PREFIX ""
        LIBRARY_OUTPUT_DIRECTORY "${MATLAB_INSTALL_PACKAGE_DIR}"
        RUNTIME_OUTPUT_DIRECTORY "${MATLAB_INSTALL_PACKAGE_DIR}"
        SKIP_BUILD_RPATH TRUE
        BUILD_RPATH ""
        INSTALL_RPATH ""
)

# Ensure the generated package directory exists.
file(MAKE_DIRECTORY "${MATLAB_INSTALL_PACKAGE_DIR}")

# Copy the generated MEX file, MATLAB sources, and Spglib license to the
# generated package after each build.
add_custom_command(TARGET spglib_mex POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E make_directory "${MATLAB_INSTALL_PACKAGE_DIR}"
        COMMAND ${CMAKE_COMMAND} -E copy_if_different
                ${MATLAB_PACKAGE_FILES} "${MATLAB_INSTALL_PACKAGE_DIR}/"
        COMMENT "Installing the kssolv.analysis.spglib MATLAB package"
)

if (SPGLIB_MATLAB_WITH_TESTS)
    enable_testing()
    matlab_add_unit_test(
            NAME matlab_spglib
            UNITTEST_FILE "${CMAKE_CURRENT_SOURCE_DIR}/test/SpglibTest.m"
            ADDITIONAL_PATH "${SPGLIB_MATLAB_OUTPUT_ROOT}"
            TIMEOUT 600
    )
endif ()

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
endif()
