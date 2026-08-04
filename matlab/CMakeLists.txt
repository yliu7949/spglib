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

# MATLAB_ROOT_DIR is the public configuration variable used by this project.
# CMake's FindMatlab module expects the mixed-case Matlab_ROOT_DIR spelling, so
# keep that implementation detail local to this directory.
set(MATLAB_ROOT_DIR "" CACHE PATH "Root directory of the MATLAB installation")
if (MATLAB_ROOT_DIR)
    set(Matlab_ROOT_DIR "${MATLAB_ROOT_DIR}")
endif ()

# Without an explicit root, FindMatlab discovers an installed MATLAB release
# using CMake's standard search logic.
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
endif ()

get_target_property(_Spglib_MATLAB_library_type Spglib::symspg TYPE)
if (NOT _Spglib_MATLAB_library_type STREQUAL "STATIC_LIBRARY")
    message(FATAL_ERROR
            "The MATLAB package requires a static Spglib::symspg target; "
            "the configured target is ${_Spglib_MATLAB_library_type}."
    )
endif ()

if (TARGET Spglib_symspg)
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

# Define the MATLAB source, example, and license files explicitly so that
# unrelated files are never included in the install bundle.
set(MATLAB_PACKAGE_FILES
        "${MATLAB_PACKAGE_DIR}/Spglib.m"
        "${MATLAB_PACKAGE_DIR}/SpglibError.m"
        "${MATLAB_PACKAGE_DIR}/Spglib.rights"
        "${MATLAB_PACKAGE_DIR}/SpglibTest.m"
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

# Keep MATLAB sources, the usage example, and the Spglib license synchronized
# even when an incremental build does not need to relink the MEX file.
add_custom_target(spglib_matlab_package_files
        COMMAND ${CMAKE_COMMAND} -E make_directory "${MATLAB_INSTALL_PACKAGE_DIR}"
        COMMAND ${CMAKE_COMMAND} -E copy_if_different
                ${MATLAB_PACKAGE_FILES} "${MATLAB_INSTALL_PACKAGE_DIR}/"
        DEPENDS ${MATLAB_PACKAGE_FILES}
        COMMENT "Synchronizing the kssolv.analysis.spglib MATLAB package files"
)
add_dependencies(spglib_mex spglib_matlab_package_files)

if (SPGLIB_MATLAB_WITH_TESTS)
    enable_testing()
    matlab_add_unit_test(
            NAME matlab_spglib
            UNITTEST_FILE
                    "${CMAKE_CURRENT_SOURCE_DIR}/cmake/run_spglib_package_tests.m"
            CUSTOM_TEST_COMMAND
                    "run_spglib_package_tests('${MATLAB_INSTALL_PACKAGE_DIR}/SpglibTest.m')"
            ADDITIONAL_PATH "${SPGLIB_MATLAB_OUTPUT_ROOT}"
            TIMEOUT 600
    )
    add_test(
            NAME matlab_package_layout
            COMMAND ${CMAKE_COMMAND}
                    -DMATLAB_PACKAGE_DIR=${MATLAB_INSTALL_PACKAGE_DIR}
                    -DMATLAB_MEX_EXTENSION=${Matlab_MEX_EXTENSION}
                    -P ${CMAKE_CURRENT_SOURCE_DIR}/cmake/check_package.cmake
    )
    set_tests_properties(matlab_spglib matlab_package_layout
            PROPERTIES LABELS "matlab;package"
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
