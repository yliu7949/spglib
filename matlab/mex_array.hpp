#pragma once

#include <cstdint>

#include "mex.h"

namespace mexutil {

inline void requireInputCount(int actual, int expected,
                              char const* function_name) {
    if (actual != expected) {
        mexErrMsgIdAndTxt("Spglib:invalidNumInputs",
                          "Incorrect number of inputs for %s.", function_name);
    }
}

template <typename Value>
mxArray* makeScalar(Value value) {
    return mxCreateDoubleScalar(static_cast<double>(value));
}

inline mxArray* makeString(char const* value) { return mxCreateString(value); }

template <typename Matrix>
mxArray* makeDoubleMatrix(Matrix const& data, mwSize rows, mwSize columns) {
    mxArray* matrix = mxCreateDoubleMatrix(rows, columns, mxREAL);
    double* output = mxGetPr(matrix);
    for (mwSize row = 0; row < rows; ++row) {
        for (mwSize column = 0; column < columns; ++column) {
            output[row + column * rows] = data[row][column];
        }
    }
    return matrix;
}

template <typename Values>
mxArray* makeIntVector(Values const& data, mwSize size) {
    mxArray* array = mxCreateNumericMatrix(size, 1, mxINT32_CLASS, mxREAL);
    auto* output = static_cast<std::int32_t*>(mxGetData(array));
    for (mwSize index = 0; index < size; ++index) {
        output[index] = static_cast<std::int32_t>(data[index]);
    }
    return array;
}

template <typename Matrix>
mxArray* makeDoubleNx3(Matrix const& data, mwSize rows) {
    return makeDoubleMatrix(data, rows, 3);
}

template <typename Tensor>
mxArray* makeRotations(Tensor const& data, mwSize count) {
    mwSize dimensions[3] = {count, 3, 3};
    mxArray* array = mxCreateNumericArray(3, dimensions, mxINT32_CLASS, mxREAL);
    auto* output = static_cast<std::int32_t*>(mxGetData(array));
    for (mwSize operation = 0; operation < count; ++operation) {
        for (mwSize row = 0; row < 3; ++row) {
            for (mwSize column = 0; column < 3; ++column) {
                output[operation + row * count + column * count * 3] =
                    static_cast<std::int32_t>(data[operation][row][column]);
            }
        }
    }
    return array;
}

template <typename Value>
void setScalarField(mxArray* structure, mwIndex index, char const* field_name,
                    Value value) {
    mxSetField(structure, index, field_name, makeScalar(value));
}

inline void setStringField(mxArray* structure, mwIndex index,
                           char const* field_name, char const* value) {
    mxSetField(structure, index, field_name, makeString(value));
}

template <typename Matrix>
void setDoubleMatrixField(mxArray* structure, mwIndex index,
                          char const* field_name, Matrix const& data,
                          mwSize rows, mwSize columns) {
    mxSetField(structure, index, field_name,
               makeDoubleMatrix(data, rows, columns));
}

template <typename Values>
void setIntArrayField(mxArray* structure, mwIndex index, char const* field_name,
                      Values const& data, mwSize size) {
    mxSetField(structure, index, field_name, makeIntVector(data, size));
}

template <typename Matrix>
void setDouble2DArrayField(mxArray* structure, mwIndex index,
                           char const* field_name, Matrix const& data,
                           mwSize rows) {
    mxSetField(structure, index, field_name, makeDoubleNx3(data, rows));
}

template <typename Tensor>
void set3DIntArrayField(mxArray* structure, mwIndex index,
                        char const* field_name, Tensor const& data,
                        mwSize count) {
    mxSetField(structure, index, field_name, makeRotations(data, count));
}

}  // namespace mexutil
