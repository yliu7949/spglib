# Spglib matlab

## 1. Compile the MATLAB Module

The recommended build configures the MATLAB interface together with spglib:

```shell
cmake -S . -B build \
    -DSPGLIB_SHARED_LIBS=OFF \
    -DSPGLIB_WITH_MATLAB=ON \
    -DMatlab_ROOT_DIR=/path/to/MATLAB
cmake --build build --config Release
```

If CMake can discover MATLAB automatically, omit `Matlab_ROOT_DIR`. The
generated package is under
`build/matlab/install/+kssolv/+analysis/+spglib`.

The MATLAB interface can also be built against an installed spglib package:

```shell
cmake -S matlab -B matlab-build \
    -DSpglib_DIR=/path/to/lib/cmake/Spglib \
    -DMatlab_ROOT_DIR=/path/to/MATLAB
cmake --build matlab-build --config Release
```

Add the generated `install` directory (the directory containing `+kssolv`) to
the MATLAB path before calling the package.

To build and run the C and MATLAB tests together:

```shell
cmake -S . -B build \
    -DSPGLIB_SHARED_LIBS=OFF \
    -DSPGLIB_WITH_MATLAB=ON \
    -DSPGLIB_WITH_TESTS=ON \
    -DMatlab_ROOT_DIR=/path/to/MATLAB
cmake --build build --config Release
ctest --test-dir build --output-on-failure -C Release
```

## 2. Usage Example

The `test/SpglibTest.m` file contains concrete usage examples for reference.

For example, to get the version number:

```matlab
disp(kssolv.analysis.spglib.Spglib.getVersion())
```
