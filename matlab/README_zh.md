# Spglib matlab 模块

## 一、编译 matlab 模块

推荐在 spglib 主工程中直接启用 MATLAB 接口：

```shell
cmake -S . -B build \
    -DSPGLIB_SHARED_LIBS=OFF \
    -DSPGLIB_WITH_MATLAB=ON \
    -DMATLAB_ROOT_DIR=/path/to/MATLAB
cmake --build build --config Release
```

`MATLAB_ROOT_DIR` 必须指向 MATLAB 安装根目录，而不是其中的 `bin` 目录。各平台
的典型值如下：

```text
macOS:   /Applications/MATLAB_R2026b.app
Linux:   /usr/local/MATLAB/R2026b
Windows: C:/Program Files/MATLAB/R2026b
```

例如，macOS 可传入 `-DMATLAB_ROOT_DIR=/Applications/MATLAB_R2026b.app`，
Windows 可传入 `"-DMATLAB_ROOT_DIR=C:/Program Files/MATLAB/R2026b"`。如果
CMake 能自动找到 MATLAB，可以省略 `MATLAB_ROOT_DIR`。生成的包位于
`build/matlab/install/+kssolv/+analysis/+spglib`。

MATLAB 包必须静态链接 spglib。使用 `SPGLIB_SHARED_LIBS=ON` 的配置，以及
独立构建时解析到共享版 `Spglib::symspg` 的配置都会被拒绝。静态链接可使生成
的包保持自包含，并适用于 MATLAB Runtime 部署。

也可以单独配置 MATLAB 接口，并链接已经安装的 spglib：

```shell
cmake -S matlab -B matlab-build \
    -DSpglib_DIR=/path/to/lib/cmake/Spglib \
    -DMATLAB_ROOT_DIR=/path/to/MATLAB
cmake --build matlab-build --config Release
```

将生成的 `install` 目录（即包含 `+kssolv` 的目录）添加到 MATLAB 路径后即可调用。

同时构建并运行 C 和 MATLAB 测试：

```shell
cmake -S . -B build \
    -DSPGLIB_SHARED_LIBS=OFF \
    -DSPGLIB_WITH_MATLAB=ON \
    -DSPGLIB_WITH_TESTS=ON \
    -DMATLAB_ROOT_DIR=/path/to/MATLAB
cmake --build build --config Release
ctest --test-dir build --output-on-failure -C Release
```

只运行 MATLAB 相关的单元测试、配置测试、独立构建测试和包布局测试时，可使用：
`ctest --test-dir build -L matlab --output-on-failure -C Release`。

## 二、使用示例

生成的安装包包含 `SpglibTest.m`，其中提供了许多具体的使用示例，也可以作为
MATLAB 单元测试运行。

例如，获取版本号：

```matlab
disp(kssolv.analysis.spglib.Spglib.getVersion())
```
