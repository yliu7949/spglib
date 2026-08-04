# Spglib matlab 模块

## 一、编译 matlab 模块

推荐在 spglib 主工程中直接启用 MATLAB 接口：

```shell
cmake -S . -B build \
    -DSPGLIB_SHARED_LIBS=OFF \
    -DSPGLIB_WITH_MATLAB=ON \
    -DMatlab_ROOT_DIR=/path/to/MATLAB
cmake --build build --config Release
```

如果 CMake 能自动找到 MATLAB，可以省略 `Matlab_ROOT_DIR`。生成的包位于
`build/matlab/install/+kssolv/+analysis/+spglib`。

也可以单独配置 MATLAB 接口，并链接已经安装的 spglib：

```shell
cmake -S matlab -B matlab-build \
    -DSpglib_DIR=/path/to/lib/cmake/Spglib \
    -DMatlab_ROOT_DIR=/path/to/MATLAB
cmake --build matlab-build --config Release
```

将生成的 `install` 目录（即包含 `+kssolv` 的目录）添加到 MATLAB 路径后即可调用。

同时构建并运行 C 和 MATLAB 测试：

```shell
cmake -S . -B build \
    -DSPGLIB_SHARED_LIBS=OFF \
    -DSPGLIB_WITH_MATLAB=ON \
    -DSPGLIB_WITH_TESTS=ON \
    -DMatlab_ROOT_DIR=/path/to/MATLAB
cmake --build build --config Release
ctest --test-dir build --output-on-failure -C Release
```

## 二、使用示例

`test/SpglibTest.m` 文件中包含了许多具体的可供参考的使用示例。

例如，获取版本号：

```matlab
disp(kssolv.analysis.spglib.Spglib.getVersion())
```
