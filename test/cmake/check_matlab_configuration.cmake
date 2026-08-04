if (NOT DEFINED MODE OR NOT DEFINED SOURCE_DIR OR NOT DEFINED BINARY_DIR)
    message(FATAL_ERROR "MODE, SOURCE_DIR, and BINARY_DIR are required")
endif ()

set(configure_command
        "${CMAKE_COMMAND}"
        -S "${SOURCE_DIR}"
        -B "${BINARY_DIR}"
)
if (DEFINED GENERATOR AND NOT GENERATOR STREQUAL "")
    list(APPEND configure_command -G "${GENERATOR}")
endif ()
if (DEFINED GENERATOR_PLATFORM AND NOT GENERATOR_PLATFORM STREQUAL "")
    list(APPEND configure_command -A "${GENERATOR_PLATFORM}")
endif ()
if (DEFINED GENERATOR_TOOLSET AND NOT GENERATOR_TOOLSET STREQUAL "")
    list(APPEND configure_command -T "${GENERATOR_TOOLSET}")
endif ()
if (DEFINED C_COMPILER AND NOT C_COMPILER STREQUAL "")
    list(APPEND configure_command "-DCMAKE_C_COMPILER=${C_COMPILER}")
endif ()

list(APPEND configure_command
        -DSPGLIB_WITH_TESTS=OFF
        -DSPGLIB_INSTALL=OFF
)

if (MODE STREQUAL "disabled")
    list(APPEND configure_command
            -DSPGLIB_WITH_MATLAB=OFF
            -DCMAKE_DISABLE_FIND_PACKAGE_Matlab=TRUE
    )
elseif (MODE STREQUAL "shared")
    list(APPEND configure_command
            -DSPGLIB_WITH_MATLAB=ON
            -DSPGLIB_SHARED_LIBS=ON
    )
else ()
    message(FATAL_ERROR "Unknown MODE: ${MODE}")
endif ()

execute_process(
        COMMAND ${configure_command}
        RESULT_VARIABLE configure_result
        OUTPUT_VARIABLE configure_stdout
        ERROR_VARIABLE configure_stderr
)
set(configure_output "${configure_stdout}\n${configure_stderr}")

if (MODE STREQUAL "disabled")
    if (NOT configure_result EQUAL 0)
        message(FATAL_ERROR
                "Configuring with MATLAB disabled failed:\n${configure_output}"
        )
    endif ()
elseif (configure_result EQUAL 0)
    message(FATAL_ERROR "A shared MATLAB configuration was unexpectedly accepted")
elseif (NOT configure_output MATCHES
        "SPGLIB_WITH_MATLAB requires a static spglib library")
    message(FATAL_ERROR
            "The shared configuration failed for an unexpected reason:\n${configure_output}"
    )
endif ()
