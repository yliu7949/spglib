set_property(DIRECTORY APPEND
        PROPERTY LABELS c_api
)

Spglib_add_gtest(TARGET Spglib_tests SOURCES
        test_error.cpp
        test_symmetry_search.cpp
        test_spacegroup_type_search.cpp
        test_find_primitive_cell.cpp
        test_refine_cell.cpp
        test_standardization.cpp
        test_dataset_access.cpp
        test_magnetic_symmetry.cpp
        test_magnetic_dataset.cpp
        test_kpoints.cpp
)

# These tests intentionally call deprecated C APIs to preserve compatibility
# coverage. Suppress the expected diagnostics only for the affected sources.
set_source_files_properties(
        "${CMAKE_CURRENT_SOURCE_DIR}/test_symmetry_search.cpp"
        "${CMAKE_CURRENT_SOURCE_DIR}/test_spacegroup_type_search.cpp"
        "${CMAKE_CURRENT_SOURCE_DIR}/test_dataset_access.cpp"
        "${CMAKE_CURRENT_SOURCE_DIR}/test_magnetic_symmetry.cpp"
        DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/../.."
        PROPERTIES COMPILE_OPTIONS
        "$<$<CXX_COMPILER_ID:MSVC>:/wd4996>;$<$<CXX_COMPILER_ID:GNU,Clang,AppleClang,Intel,IntelLLVM>:-Wno-deprecated-declarations>"
)
