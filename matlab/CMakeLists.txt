cmake_minimum_required(VERSION 3.25...3.29)

# 设置项目名称和版本，并指定 C++ 语言
project(Spglib_MATLAB
        VERSION 1.0
        LANGUAGES CXX
)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# 设置 MATLAB 的根路径
set(MATLAB_ROOT "/Applications/MATLAB_R2024b.app")
find_package(Matlab REQUIRED)

# 设置 Spglib 库的根路径
set(SPGLIB_ROOT "/Users/liu/Documents/GitHub/spglib/install")
find_package(Spglib REQUIRED)

# 添加头文件路径
include_directories("${SPGLIB_ROOT}/include")

# 设置源文件
set(SOURCE_FILES _spglib.cpp)

# 使用 MEX 编译器生成 MEX 文件
matlab_add_mex(NAME spglib_mex SRC ${SOURCE_FILES} LINK_TO Spglib::symspg)

# 判断 spglib 文件夹是否存在，若不存在则创建
file(MAKE_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/spglib")

# 设置生成文件的后缀为 MATLAB MEX 格式
set_target_properties(spglib_mex PROPERTIES
        OUTPUT_NAME "symspg"
        PREFIX ""
        LIBRARY_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/spglib"
)
