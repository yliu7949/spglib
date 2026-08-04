if (NOT DEFINED MATLAB_PACKAGE_DIR OR NOT DEFINED MATLAB_MEX_EXTENSION)
    message(FATAL_ERROR
            "MATLAB_PACKAGE_DIR and MATLAB_MEX_EXTENSION are required"
    )
endif ()

set(required_package_files
        "Spglib.m"
        "SpglibError.m"
        "Spglib.rights"
        "SpglibTest.m"
        "symspg.${MATLAB_MEX_EXTENSION}"
)
foreach (package_file IN LISTS required_package_files)
    if (NOT EXISTS "${MATLAB_PACKAGE_DIR}/${package_file}")
        message(FATAL_ERROR
                "Generated MATLAB package is missing ${package_file}: "
                "${MATLAB_PACKAGE_DIR}"
        )
    endif ()
endforeach ()
