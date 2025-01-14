# Spglib matlab

## 1. Compile the MATLAB Module

First, compile the `spglib` library statically using the following commands:

```shell
cd spglib

# Windows/macOS
mkdir build
cmake -S . -B build -DSPGLIB_SHARED_LIBS=OFF -DSPGLIB_WITH_TESTS=OFF -DCMAKE_BUILD_TYPE=Release
cmake --build build --config Release
cmake --install build --config Release --prefix="./install"

# Linux
mkdir build
cmake -S . -B build -DSPGLIB_SHARED_LIBS=OFF -DSPGLIB_WITH_TESTS=OFF -DCMAKE_BUILD_TYPE=Release -DCMAKE_C_FLAGS="-fPIC"
cmake --build build --config Release
cmake --install build --config Release --prefix="./install"
```

Then, go into the `matlab` folder, modify the `MATLAB_ROOT` value in `CMakeLists.txt` according to your actual setup, and execute the following commands to compile:

```shell
# macOS/Linux
mkdir build
cmake -S . -B build
cmake --build build

# Windows
mkdir build
cmake -G "Ninja" -S . -B build -DCMAKE_CXX_COMPILER=clang-cl
cmake --build build
```

After the compilation is complete, an `install` folder will appear in the current directory. Copy the `install` folder to the desired location, and after adding it to the MATLAB path, you can call the related functions in MATLAB.

## 2. Usage Example

The `SpglibTest.m` file in the `install` folder contains many concrete usage examples for reference.

For example, to get the version number:

```matlab
disp(Spglib.getVersion())
```
