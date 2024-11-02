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
set(Spglib_DIR "/Users/liu/Documents/GitHub/spglib/install/lib/cmake/Spglib")
find_package(Matlab)
find_package(Spglib)

# 添加头文件路径
include_directories("/Users/liu/Documents/GitHub/spglib/install/include")

# 设置源文件
set(SOURCE_FILES _spglib.cpp)

# 使用 MEX 编译器生成 MEX 文件
matlab_add_mex(NAME spglib_mex SRC ${SOURCE_FILES} LINK_TO Spglib::symspg)

# 设置生成文件的后缀为 MATLAB MEX 格式
set_target_properties(spglib_mex PROPERTIES OUTPUT_NAME "symspg" PREFIX "" SUFFIX ".mexmaca64")
