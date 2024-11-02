#include <string>
#include <unordered_map>
#include "mex.h"
#include "spglib.h"

// 定义一个具有静态方法的类
class SpglibFunctions {
   public:
    // 声明静态方法接口
    static void spg_get_dataset_mex(int nlhs, mxArray* plhs[], int nrhs,
                                    mxArray const* prhs[]);
};

// 定义类型别名，用于指向静态方法的函数指针
typedef void (*SpglibFunction)(int, mxArray*[], int, mxArray const*[]);

// mexFunction 主函数
void mexFunction(int nlhs, mxArray* plhs[], int nrhs, mxArray const* prhs[]) {
    if (nrhs < 1) {
        mexErrMsgIdAndTxt("Spglib:invalidNumInputs",
                          "At least one input required.");
    }
    if (!mxIsChar(prhs[0])) {
        mexErrMsgIdAndTxt("Spglib:invalidInput",
                          "First input must be a function name string.");
    }

    char function_name[256];
    mxGetString(prhs[0], function_name, sizeof(function_name));

    // 将函数名称映射到类的静态方法
    std::unordered_map<std::string, SpglibFunction> function_map = {
        {"spg_get_dataset", SpglibFunctions::spg_get_dataset_mex},
    };

    // 查找并调用对应的函数
    auto it = function_map.find(function_name);
    if (it != function_map.end()) {
        it->second(nlhs, plhs, nrhs - 1, prhs + 1);
    } else {
        mexErrMsgIdAndTxt("Spglib:invalidFunction", "Unknown function name.");
    }
}

// 宏定义，用于简化字段赋值

// 设置标量字段
#define SET_SCALAR_FIELD(matlab_struct, index, field_name, value) \
    mxSetField(matlab_struct, index, field_name, mxCreateDoubleScalar(value))

// 设置字符串字段
#define SET_STRING_FIELD(matlab_struct, index, field_name, value) \
    mxSetField(matlab_struct, index, field_name, mxCreateString(value))

// 设置 3x3 double 矩阵字段
#define SET_DOUBLE_MATRIX_FIELD(matlab_struct, index, field_name, data, rows, \
                                cols)                                         \
    do {                                                                      \
        mxArray* matrix = mxCreateDoubleMatrix(rows, cols, mxREAL);           \
        double* ptr = mxGetPr(matrix);                                        \
        for (int i = 0; i < rows; ++i)                                        \
            for (int j = 0; j < cols; ++j) ptr[i + j * rows] = data[i][j];    \
        mxSetField(matlab_struct, index, field_name, matrix);                 \
    } while (0)

// 设置一维整数数组字段
#define SET_INT_ARRAY_FIELD(matlab_struct, index, field_name, data,        \
                            num_elements)                                  \
    do {                                                                   \
        mxArray* array =                                                   \
            mxCreateNumericMatrix(num_elements, 1, mxINT32_CLASS, mxREAL); \
        int32_t* ptr = static_cast<int32_t*>(mxGetData(array));            \
        for (int i = 0; i < num_elements; ++i) ptr[i] = data[i];           \
        mxSetField(matlab_struct, index, field_name, array);               \
    } while (0)

// 设置二维 double 数组字段（Nx3 矩阵）
#define SET_DOUBLE_2D_ARRAY_FIELD(matlab_struct, index, field_name, data, \
                                  rows)                                   \
    do {                                                                  \
        mxArray* array = mxCreateDoubleMatrix(rows, 3, mxREAL);           \
        double* ptr = mxGetPr(array);                                     \
        for (int i = 0; i < rows; ++i)                                    \
            for (int j = 0; j < 3; ++j) ptr[i + j * rows] = data[i][j];   \
        mxSetField(matlab_struct, index, field_name, array);              \
    } while (0)

// 设置旋转矩阵数组（3x3xN int 数组）
#define SET_3D_INT_ARRAY_FIELD(matlab_struct, index, field_name, data,         \
                               num_elements)                                   \
    do {                                                                       \
        mxArray* array =                                                       \
            mxCreateNumericMatrix(3, 3 * num_elements, mxINT32_CLASS, mxREAL); \
        int32_t* ptr = static_cast<int32_t*>(mxGetData(array));                \
        for (int k = 0; k < num_elements; ++k)                                 \
            for (int i = 0; i < 3; ++i)                                        \
                for (int j = 0; j < 3; ++j)                                    \
                    ptr[i + j * 3 + k * 9] = data[k][i][j];                    \
        mxSetField(matlab_struct, index, field_name, array);                   \
    } while (0)

// dataset = symspg('spg_get_dataset', lattice, position, types, num_atom,
// symprec)
void SpglibFunctions::spg_get_dataset_mex(int nlhs, mxArray* plhs[], int nrhs,
                                          mxArray const* prhs[]) {
    /*
     SpglibDataset * spg_get_dataset(const double lattice[3][3],
                           const double position[][3],
                           const int types[],
                           const int num_atom,
                           const double symprec);
     */

    // 验证输入参数数量
    int const expected_number_of_inputs = 5;
    if (nrhs != expected_number_of_inputs) {
        mexErrMsgIdAndTxt("Spglib:invalidNumInputs",
                          "Incorrect number of inputs for spg_get_dataset.");
    }

    // 提取和验证 lattice, position, types 和 symprec 的代码
    double lattice[3][3];
    if (mxGetM(prhs[0]) != 3 || mxGetN(prhs[0]) != 3) {
        mexErrMsgIdAndTxt("Spglib:invalidLattice",
                          "Lattice matrix must be 3x3.");
    }

    double* lattice_ptr = mxGetPr(prhs[0]);
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            lattice[i][j] = lattice_ptr[i + 3 * j];
        }
    }

    int num_atom = mxGetM(prhs[1]);
    double position[num_atom][3];
    if (mxGetN(prhs[1]) != 3) {
        mexErrMsgIdAndTxt("Spglib:invalidPosition",
                          "Position array must be Nx3.");
    }
    double* position_ptr = mxGetPr(prhs[1]);
    for (int i = 0; i < num_atom; i++) {
        for (int j = 0; j < 3; j++) {
            position[i][j] = position_ptr[i + num_atom * j];
        }
    }

    if (mxGetNumberOfElements(prhs[2]) != num_atom) {
        mexErrMsgIdAndTxt("Spglib:invalidTypes",
                          "Types array size must match the number of atoms.");
    }
    int types[num_atom];
    double* types_ptr = mxGetPr(prhs[2]);
    for (int i = 0; i < num_atom; i++) {
        types[i] = static_cast<int>(types_ptr[i]);
    }

    if (!mxIsDouble(prhs[4]) || mxGetNumberOfElements(prhs[4]) != 1) {
        mexErrMsgIdAndTxt("Spglib:invalidSymprec", "Symprec must be a scalar.");
    }
    double symprec = mxGetScalar(prhs[4]);

    // 调用 spg_get_dataset
    SpglibDataset* dataset =
        spg_get_dataset(lattice, position, types, num_atom, symprec);
    if (dataset == nullptr) {
        mexErrMsgIdAndTxt("Spglib:datasetError", "Failed to get dataset.");
    }

    // 定义 MATLAB 结构体的字段名称
    int const field_count = 24;
    char const* field_names[field_count] = {"spacegroup_number",
                                            "hall_number",
                                            "international_symbol",
                                            "hall_symbol",
                                            "choice",
                                            "transformation_matrix",
                                            "origin_shift",
                                            "n_operations",
                                            "rotations",
                                            "translations",
                                            "n_atoms",
                                            "wyckoffs",
                                            "site_symmetry_symbols",
                                            "equivalent_atoms",
                                            "crystallographic_orbits",
                                            "primitive_lattice",
                                            "mapping_to_primitive",
                                            "n_std_atoms",
                                            "std_lattice",
                                            "std_types",
                                            "std_positions",
                                            "std_rotation_matrix",
                                            "std_mapping_to_primitive",
                                            "pointgroup_symbol"};

    // 创建 MATLAB 结构体
    plhs[0] = mxCreateStructMatrix(1, 1, field_count, field_names);

    // 使用宏填充结构体字段
    SET_SCALAR_FIELD(plhs[0], 0, "spacegroup_number",
                     dataset->spacegroup_number);
    SET_SCALAR_FIELD(plhs[0], 0, "hall_number", dataset->hall_number);
    SET_STRING_FIELD(plhs[0], 0, "international_symbol",
                     dataset->international_symbol);
    SET_STRING_FIELD(plhs[0], 0, "hall_symbol", dataset->hall_symbol);
    SET_STRING_FIELD(plhs[0], 0, "choice", dataset->choice);
    SET_DOUBLE_MATRIX_FIELD(plhs[0], 0, "transformation_matrix",
                            dataset->transformation_matrix, 3, 3);

    // origin_shift (3x1 double array)
    mxArray* origin_shift_mx = mxCreateDoubleMatrix(3, 1, mxREAL);
    double* origin_shift_ptr = mxGetPr(origin_shift_mx);
    for (int i = 0; i < 3; ++i) origin_shift_ptr[i] = dataset->origin_shift[i];
    mxSetField(plhs[0], 0, "origin_shift", origin_shift_mx);

    SET_SCALAR_FIELD(plhs[0], 0, "n_operations", dataset->n_operations);
    SET_3D_INT_ARRAY_FIELD(plhs[0], 0, "rotations", dataset->rotations,
                           dataset->n_operations);
    SET_DOUBLE_2D_ARRAY_FIELD(plhs[0], 0, "translations", dataset->translations,
                              dataset->n_operations);
    SET_SCALAR_FIELD(plhs[0], 0, "n_atoms", dataset->n_atoms);
    SET_INT_ARRAY_FIELD(plhs[0], 0, "wyckoffs", dataset->wyckoffs,
                        dataset->n_atoms);

    // 设置 site_symmetry_symbols 为字符数组
    char const* site_symmetry_symbols[dataset->n_atoms];
    for (int i = 0; i < dataset->n_atoms; ++i) {
        site_symmetry_symbols[i] = dataset->site_symmetry_symbols[i];
    }
    mxArray* site_symmetry_symbols_mx =
        mxCreateCharMatrixFromStrings(dataset->n_atoms, site_symmetry_symbols);
    mxSetField(plhs[0], 0, "site_symmetry_symbols", site_symmetry_symbols_mx);
    SET_INT_ARRAY_FIELD(plhs[0], 0, "equivalent_atoms",
                        dataset->equivalent_atoms, dataset->n_atoms);
    SET_INT_ARRAY_FIELD(plhs[0], 0, "crystallographic_orbits",
                        dataset->crystallographic_orbits, dataset->n_atoms);
    SET_DOUBLE_MATRIX_FIELD(plhs[0], 0, "primitive_lattice",
                            dataset->primitive_lattice, 3, 3);
    SET_INT_ARRAY_FIELD(plhs[0], 0, "mapping_to_primitive",
                        dataset->mapping_to_primitive, dataset->n_atoms);
    SET_SCALAR_FIELD(plhs[0], 0, "n_std_atoms", dataset->n_std_atoms);
    SET_DOUBLE_MATRIX_FIELD(plhs[0], 0, "std_lattice", dataset->std_lattice, 3,
                            3);
    SET_INT_ARRAY_FIELD(plhs[0], 0, "std_types", dataset->std_types,
                        dataset->n_std_atoms);
    SET_DOUBLE_2D_ARRAY_FIELD(plhs[0], 0, "std_positions",
                              dataset->std_positions, dataset->n_std_atoms);
    SET_DOUBLE_MATRIX_FIELD(plhs[0], 0, "std_rotation_matrix",
                            dataset->std_rotation_matrix, 3, 3);
    SET_INT_ARRAY_FIELD(plhs[0], 0, "std_mapping_to_primitive",
                        dataset->std_mapping_to_primitive,
                        dataset->n_std_atoms);
    SET_STRING_FIELD(plhs[0], 0, "pointgroup_symbol",
                     dataset->pointgroup_symbol);

    // 清理 spglib 数据集
    spg_free_dataset(dataset);
}
