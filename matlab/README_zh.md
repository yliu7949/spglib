# Spglib matlab 模块

## 一、编译 matlab 模块

首先需要使用下面的命令静态编译 `spglib` 库：

```shell
cd spglib
mkdir build
cmake -S . -B build -DSPGLIB_SHARED_LIBS=OFF -DCMAKE_BUILD_TYPE=Release
cmake --build build --config Release
cmake --install build --config Release --prefix="./install"
```

然后进入 `matlab` 文件夹内，根据实际情况修改 `CMakeLists.txt` 中的 `MATLAB_ROOT` 的值，执行下面的命令编译：

```shell
mkdir build
cmake -S . -B build
cmake --build build
```

编译结束后会在当前文件夹下出现 `install` 文件夹。将 `install` 文件夹复制到需要使用的地方，将其添加到 MATLAB 路径后即可在 MATLAB 中调用相关函数。

## 二、使用示例

`install` 文件夹下的 SpglibTest.m 文件中包含了许多具体的可供参考的使用示例。

例如，获取版本号：

```matlab
disp(Spglib.getVersion())
```
