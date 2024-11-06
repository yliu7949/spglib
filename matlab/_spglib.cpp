#include <string>
#include <unordered_map>
#include "mex.h"
#include "spglib.h"

// 定义一个具有静态方法的类
class SpglibFunctions {
   public:
    // 声明静态方法接口
    static void spg_get_version_mex(int nlhs, mxArray* plhs[], int nrhs,
                                    mxArray const* prhs[]);
    static void spg_get_version_full_mex(int nlhs, mxArray* plhs[], int nrhs,
                                         mxArray const* prhs[]);
    static void spg_get_commit_mex(int nlhs, mxArray* plhs[], int nrhs,
                                   mxArray const* prhs[]);
    static void spg_get_major_version_mex(int nlhs, mxArray* plhs[], int nrhs,
                                          mxArray const* prhs[]);
    static void spg_get_minor_version_mex(int nlhs, mxArray* plhs[], int nrhs,
                                          mxArray const* prhs[]);
    static void spg_get_micro_version_mex(int nlhs, mxArray* plhs[], int nrhs,
                                          mxArray const* prhs[]);
    static void spg_get_dataset_mex(int nlhs, mxArray* plhs[], int nrhs,
                                    mxArray const* prhs[]);
    static void spg_get_magnetic_dataset_mex(int nlhs, mxArray* plhs[],
                                             int nrhs, mxArray const* prhs[]);
    static void spgms_get_magnetic_dataset_mex(int nlhs, mxArray* plhs[],
                                               int nrhs, mxArray const* prhs[]);
    static void spgat_get_dataset_mex(int nlhs, mxArray* plhs[], int nrhs,
                                      mxArray const* prhs[]);
    static void spg_get_dataset_with_hall_number_mex(int nlhs, mxArray* plhs[],
                                                     int nrhs,
                                                     mxArray const* prhs[]);
    static void spgat_get_dataset_with_hall_number_mex(int nlhs,
                                                       mxArray* plhs[],
                                                       int nrhs,
                                                       mxArray const* prhs[]);
    static void spg_get_symmetry_with_collinear_spin_mex(int nlhs,
                                                         mxArray* plhs[],
                                                         int nrhs,
                                                         mxArray const* prhs[]);
    static void spgat_get_symmetry_with_collinear_spin_mex(
        int nlhs, mxArray* plhs[], int nrhs, mxArray const* prhs[]);
    static void spgms_get_symmetry_with_collinear_spin_mex(
        int nlhs, mxArray* plhs[], int nrhs, mxArray const* prhs[]);
    static void spg_get_symmetry_with_site_tensors_mex(int nlhs,
                                                       mxArray* plhs[],
                                                       int nrhs,
                                                       mxArray const* prhs[]);
    static void spgat_get_symmetry_with_site_tensors_mex(int nlhs,
                                                         mxArray* plhs[],
                                                         int nrhs,
                                                         mxArray const* prhs[]);
    static void spgms_get_symmetry_with_site_tensors_mex(int nlhs,
                                                         mxArray* plhs[],
                                                         int nrhs,
                                                         mxArray const* prhs[]);
    static void spg_get_spacegroup_type_from_symmetry_mex(
        int nlhs, mxArray* plhs[], int nrhs, mxArray const* prhs[]);
    static void spg_get_magnetic_spacegroup_type_from_symmetry_mex(
        int nlhs, mxArray* plhs[], int nrhs, mxArray const* prhs[]);
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
        {"spg_get_version", SpglibFunctions::spg_get_version_mex},
        {"spg_get_version_full", SpglibFunctions::spg_get_version_full_mex},
        {"spg_get_commit", SpglibFunctions::spg_get_commit_mex},
        {"spg_get_major_version", SpglibFunctions::spg_get_major_version_mex},
        {"spg_get_minor_version", SpglibFunctions::spg_get_minor_version_mex},
        {"spg_get_micro_version", SpglibFunctions::spg_get_micro_version_mex},
        {"spg_get_dataset", SpglibFunctions::spg_get_dataset_mex},
        {"spg_get_magnetic_dataset",
         SpglibFunctions::spg_get_magnetic_dataset_mex},
        {"spgms_get_magnetic_dataset",
         SpglibFunctions::spgms_get_magnetic_dataset_mex},
        {"spgat_get_dataset", SpglibFunctions::spgat_get_dataset_mex},
        {"spg_get_dataset_with_hall_number",
         SpglibFunctions::spg_get_dataset_with_hall_number_mex},
        {"spgat_get_dataset_with_hall_number",
         SpglibFunctions::spgat_get_dataset_with_hall_number_mex},
        {"spg_get_symmetry_with_collinear_spin",
         SpglibFunctions::spg_get_symmetry_with_collinear_spin_mex},
        {"spgat_get_symmetry_with_collinear_spin",
         SpglibFunctions::spgat_get_symmetry_with_collinear_spin_mex},
        {"spgms_get_symmetry_with_collinear_spin",
         SpglibFunctions::spgms_get_symmetry_with_collinear_spin_mex},
        {"spg_get_symmetry_with_site_tensors",
         SpglibFunctions::spg_get_symmetry_with_site_tensors_mex},
        {"spgat_get_symmetry_with_site_tensors",
         SpglibFunctions::spgat_get_symmetry_with_site_tensors_mex},
        {"spgms_get_symmetry_with_site_tensors",
         SpglibFunctions::spgms_get_symmetry_with_site_tensors_mex},
        {"spg_get_spacegroup_type_from_symmetry",
         SpglibFunctions::spg_get_spacegroup_type_from_symmetry_mex},
        {"spg_get_magnetic_spacegroup_type_from_symmetry",
         SpglibFunctions::spg_get_magnetic_spacegroup_type_from_symmetry_mex},
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

// 设置旋转矩阵数组（Nx3x3 int 数组）
#define SET_3D_INT_ARRAY_FIELD(matlab_struct, index, field_name, data,         \
                               num_elements)                                   \
    do {                                                                       \
        mwSize dims[3] = {static_cast<mwSize>(num_elements), 3, 3};            \
        mxArray* array = mxCreateNumericArray(3, dims, mxINT32_CLASS, mxREAL); \
        int32_t* ptr = static_cast<int32_t*>(mxGetData(array));                \
        for (int k = 0; k < num_elements; ++k)                                 \
            for (int i = 0; i < 3; ++i)                                        \
                for (int j = 0; j < 3; ++j)                                    \
                    ptr[k + i * num_elements + j * num_elements * 3] =         \
                        data[k][i][j];                                         \
        mxSetField(matlab_struct, index, field_name, array);                   \
    } while (0)

// version = symspg('spg_get_version')
void SpglibFunctions::spg_get_version_mex(int nlhs, mxArray* plhs[], int nrhs,
                                          mxArray const* prhs[]) {
    if (nrhs != 0) {
        mexErrMsgIdAndTxt("Spglib:invalidNumInputs",
                          "spg_get_version does not take any input arguments.");
    }

    char const* version = spg_get_version();
    plhs[0] = mxCreateString(version);
}

// version = symspg('spg_get_version_full')
void SpglibFunctions::spg_get_version_full_mex(int nlhs, mxArray* plhs[],
                                               int nrhs,
                                               mxArray const* prhs[]) {
    if (nrhs != 0) {
        mexErrMsgIdAndTxt(
            "Spglib:invalidNumInputs",
            "spg_get_version_full does not take any input arguments.");
    }

    char const* version_full = spg_get_version_full();
    plhs[0] = mxCreateString(version_full);
}

// commit = symspg('spg_get_commit')
void SpglibFunctions::spg_get_commit_mex(int nlhs, mxArray* plhs[], int nrhs,
                                         mxArray const* prhs[]) {
    if (nrhs != 0) {
        mexErrMsgIdAndTxt("Spglib:invalidNumInputs",
                          "spg_get_commit does not take any input arguments.");
    }

    char const* commit = spg_get_commit();
    plhs[0] = mxCreateString(commit);
}

// version = symspg('spg_get_major_version')
void SpglibFunctions::spg_get_major_version_mex(int nlhs, mxArray* plhs[],
                                                int nrhs,
                                                mxArray const* prhs[]) {
    if (nrhs != 0) {
        mexErrMsgIdAndTxt(
            "Spglib:invalidNumInputs",
            "spg_get_major_version does not take any input arguments.");
    }
    int major_version = spg_get_major_version();
    plhs[0] = mxCreateDoubleScalar(static_cast<double>(major_version));
}

// version = symspg('spg_get_minor_version')
void SpglibFunctions::spg_get_minor_version_mex(int nlhs, mxArray* plhs[],
                                                int nrhs,
                                                mxArray const* prhs[]) {
    if (nrhs != 0) {
        mexErrMsgIdAndTxt(
            "Spglib:invalidNumInputs",
            "spg_get_minor_version does not take any input arguments.");
    }
    int minor_version = spg_get_minor_version();
    plhs[0] = mxCreateDoubleScalar(static_cast<double>(minor_version));
}

// version = symspg('spg_get_micro_version')
void SpglibFunctions::spg_get_micro_version_mex(int nlhs, mxArray* plhs[],
                                                int nrhs,
                                                mxArray const* prhs[]) {
    if (nrhs != 0) {
        mexErrMsgIdAndTxt(
            "Spglib:invalidNumInputs",
            "spg_get_micro_version does not take any input arguments.");
    }
    int micro_version = spg_get_micro_version();
    plhs[0] = mxCreateDoubleScalar(static_cast<double>(micro_version));
}

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

// dataset = symspg('spg_get_magnetic_dataset', lattice, position, types,
// tensors, tensor_rank, num_atom, is_axial, symprec)
void SpglibFunctions::spg_get_magnetic_dataset_mex(int nlhs, mxArray* plhs[],
                                                   int nrhs,
                                                   mxArray const* prhs[]) {
    /*
     SpglibMagneticDataset *spg_get_magnetic_dataset(
         double const lattice[3][3],
         double const position[][3],
         int const types[],
         double const *tensors,
         int const tensor_rank,
         int const num_atom,
         int const is_axial,
         double const symprec);
    */

    // 验证输入参数数量
    int const expected_number_of_inputs = 8;
    if (nrhs != expected_number_of_inputs) {
        mexErrMsgIdAndTxt(
            "Spglib:invalidNumInputs",
            "Incorrect number of inputs for spg_get_magnetic_dataset.");
    }

    // 提取和验证 lattice 参数
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

    // 提取和验证 position 参数
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

    // 提取和验证 types 参数
    if (mxGetNumberOfElements(prhs[2]) != num_atom) {
        mexErrMsgIdAndTxt("Spglib:invalidTypes",
                          "Types array size must match the number of atoms.");
    }
    int types[num_atom];
    double* types_ptr = mxGetPr(prhs[2]);
    for (int i = 0; i < num_atom; i++) {
        types[i] = static_cast<int>(types_ptr[i]);
    }

    // 提取和验证 tensors 参数
    double* tensors = mxGetPr(prhs[3]);

    // 提取和验证 tensor_rank 参数
    int tensor_rank = static_cast<int>(mxGetScalar(prhs[4]));

    // 提取和验证 is_axial 参数
    int is_axial = static_cast<int>(mxGetScalar(prhs[6]));

    // 提取和验证 symprec 参数
    if (!mxIsDouble(prhs[7]) || mxGetNumberOfElements(prhs[7]) != 1) {
        mexErrMsgIdAndTxt("Spglib:invalidSymprec", "Symprec must be a scalar.");
    }
    double symprec = mxGetScalar(prhs[7]);

    // 调用 spg_get_magnetic_dataset
    SpglibMagneticDataset* dataset =
        spg_get_magnetic_dataset(lattice, position, types, tensors, tensor_rank,
                                 num_atom, is_axial, symprec);

    if (dataset == nullptr) {
        mexErrMsgIdAndTxt("Spglib:datasetError",
                          "Failed to get magnetic dataset.");
    }

    // 定义 MATLAB 结构体的字段名称
    int const field_count = 18;
    char const* field_names[field_count] = {"uni_number",
                                            "msg_type",
                                            "hall_number",
                                            "tensor_rank",
                                            "n_operations",
                                            "rotations",
                                            "translations",
                                            "time_reversals",
                                            "n_atoms",
                                            "equivalent_atoms",
                                            "transformation_matrix",
                                            "origin_shift",
                                            "n_std_atoms",
                                            "std_lattice",
                                            "std_types",
                                            "std_positions",
                                            "std_tensors",
                                            "std_rotation_matrix"};

    // 创建 MATLAB 结构体
    plhs[0] = mxCreateStructMatrix(1, 1, field_count, field_names);

    // 填充结构体字段
    SET_SCALAR_FIELD(plhs[0], 0, "uni_number", dataset->uni_number);
    SET_SCALAR_FIELD(plhs[0], 0, "msg_type", dataset->msg_type);
    SET_SCALAR_FIELD(plhs[0], 0, "hall_number", dataset->hall_number);
    SET_SCALAR_FIELD(plhs[0], 0, "tensor_rank", dataset->tensor_rank);

    SET_SCALAR_FIELD(plhs[0], 0, "n_operations", dataset->n_operations);
    SET_3D_INT_ARRAY_FIELD(plhs[0], 0, "rotations", dataset->rotations,
                           dataset->n_operations);
    SET_DOUBLE_2D_ARRAY_FIELD(plhs[0], 0, "translations", dataset->translations,
                              dataset->n_operations);

    SET_INT_ARRAY_FIELD(plhs[0], 0, "time_reversals", dataset->time_reversals,
                        dataset->n_operations);

    SET_SCALAR_FIELD(plhs[0], 0, "n_atoms", dataset->n_atoms);
    SET_INT_ARRAY_FIELD(plhs[0], 0, "equivalent_atoms",
                        dataset->equivalent_atoms, dataset->n_atoms);

    SET_DOUBLE_MATRIX_FIELD(plhs[0], 0, "transformation_matrix",
                            dataset->transformation_matrix, 3, 3);

    mxArray* origin_shift_mx = mxCreateDoubleMatrix(3, 1, mxREAL);
    double* origin_shift_ptr = mxGetPr(origin_shift_mx);
    for (int i = 0; i < 3; ++i) origin_shift_ptr[i] = dataset->origin_shift[i];
    mxSetField(plhs[0], 0, "origin_shift", origin_shift_mx);

    SET_SCALAR_FIELD(plhs[0], 0, "n_std_atoms", dataset->n_std_atoms);
    SET_DOUBLE_MATRIX_FIELD(plhs[0], 0, "std_lattice", dataset->std_lattice, 3,
                            3);
    SET_INT_ARRAY_FIELD(plhs[0], 0, "std_types", dataset->std_types,
                        dataset->n_std_atoms);
    SET_DOUBLE_2D_ARRAY_FIELD(plhs[0], 0, "std_positions",
                              dataset->std_positions, dataset->n_std_atoms);

    int tensor_elements =
        dataset->n_std_atoms * std::pow(dataset->tensor_rank, 3);
    mxArray* std_tensors_mx = mxCreateDoubleMatrix(tensor_elements, 1, mxREAL);
    double* std_tensors_ptr = mxGetPr(std_tensors_mx);
    for (int i = 0; i < tensor_elements; ++i)
        std_tensors_ptr[i] = dataset->std_tensors[i];
    mxSetField(plhs[0], 0, "std_tensors", std_tensors_mx);

    SET_DOUBLE_MATRIX_FIELD(plhs[0], 0, "std_rotation_matrix",
                            dataset->std_rotation_matrix, 3, 3);

    // 清理 spglib 磁性数据集
    spg_free_magnetic_dataset(dataset);
}

// dataset = symspg('spgms_get_magnetic_dataset', lattice, position, types,
// tensors, tensor_rank, num_atom, is_axial, symprec, angle_tolerance,
// mag_symprec)
void SpglibFunctions::spgms_get_magnetic_dataset_mex(int nlhs, mxArray* plhs[],
                                                     int nrhs,
                                                     mxArray const* prhs[]) {
    /*
     SpglibMagneticDataset *spgms_get_magnetic_dataset(
         double const lattice[3][3],
         double const position[][3],
         int const types[],
         double const *tensors,
         int const tensor_rank,
         int const num_atom,
         int const is_axial,
         double const symprec,
         double const angle_tolerance,
         double const mag_symprec);
    */

    // 验证输入参数数量
    int const expected_number_of_inputs = 10;
    if (nrhs != expected_number_of_inputs) {
        mexErrMsgIdAndTxt(
            "Spglib:invalidNumInputs",
            "Incorrect number of inputs for spgms_get_magnetic_dataset.");
    }

    // 提取和验证 lattice 参数
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

    // 提取和验证 position 参数
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

    // 提取和验证 types 参数
    if (mxGetNumberOfElements(prhs[2]) != num_atom) {
        mexErrMsgIdAndTxt("Spglib:invalidTypes",
                          "Types array size must match the number of atoms.");
    }
    int types[num_atom];
    double* types_ptr = mxGetPr(prhs[2]);
    for (int i = 0; i < num_atom; i++) {
        types[i] = static_cast<int>(types_ptr[i]);
    }

    // 提取和验证 tensors 参数
    double* tensors = mxGetPr(prhs[3]);

    // 提取和验证 tensor_rank 参数
    int tensor_rank = static_cast<int>(mxGetScalar(prhs[4]));

    // 提取和验证 is_axial 参数
    int is_axial = static_cast<int>(mxGetScalar(prhs[6]));

    // 提取和验证 symprec 参数
    if (!mxIsDouble(prhs[7]) || mxGetNumberOfElements(prhs[7]) != 1) {
        mexErrMsgIdAndTxt("Spglib:invalidSymprec", "Symprec must be a scalar.");
    }
    double symprec = mxGetScalar(prhs[7]);

    // 提取和验证 angle_tolerance 参数
    if (!mxIsDouble(prhs[8]) || mxGetNumberOfElements(prhs[8]) != 1) {
        mexErrMsgIdAndTxt("Spglib:invalidAngleTolerance",
                          "Angle tolerance must be a scalar.");
    }
    double angle_tolerance = mxGetScalar(prhs[8]);

    // 提取和验证 mag_symprec 参数
    if (!mxIsDouble(prhs[9]) || mxGetNumberOfElements(prhs[9]) != 1) {
        mexErrMsgIdAndTxt("Spglib:invalidMagSymprec",
                          "Magnetic symmetry precision must be a scalar.");
    }
    double mag_symprec = mxGetScalar(prhs[9]);

    // 调用 spgms_get_magnetic_dataset
    SpglibMagneticDataset* dataset = spgms_get_magnetic_dataset(
        lattice, position, types, tensors, tensor_rank, num_atom, is_axial,
        symprec, angle_tolerance, mag_symprec);

    if (dataset == nullptr) {
        mexErrMsgIdAndTxt("Spglib:datasetError",
                          "Failed to get magnetic dataset.");
    }

    // 定义 MATLAB 结构体的字段名称
    int const field_count = 18;
    char const* field_names[field_count] = {"uni_number",
                                            "msg_type",
                                            "hall_number",
                                            "tensor_rank",
                                            "n_operations",
                                            "rotations",
                                            "translations",
                                            "time_reversals",
                                            "n_atoms",
                                            "equivalent_atoms",
                                            "transformation_matrix",
                                            "origin_shift",
                                            "n_std_atoms",
                                            "std_lattice",
                                            "std_types",
                                            "std_positions",
                                            "std_tensors",
                                            "std_rotation_matrix"};

    // 创建 MATLAB 结构体
    plhs[0] = mxCreateStructMatrix(1, 1, field_count, field_names);

    // 填充结构体字段
    SET_SCALAR_FIELD(plhs[0], 0, "uni_number", dataset->uni_number);
    SET_SCALAR_FIELD(plhs[0], 0, "msg_type", dataset->msg_type);
    SET_SCALAR_FIELD(plhs[0], 0, "hall_number", dataset->hall_number);
    SET_SCALAR_FIELD(plhs[0], 0, "tensor_rank", dataset->tensor_rank);

    SET_SCALAR_FIELD(plhs[0], 0, "n_operations", dataset->n_operations);
    SET_3D_INT_ARRAY_FIELD(plhs[0], 0, "rotations", dataset->rotations,
                           dataset->n_operations);
    SET_DOUBLE_2D_ARRAY_FIELD(plhs[0], 0, "translations", dataset->translations,
                              dataset->n_operations);

    SET_INT_ARRAY_FIELD(plhs[0], 0, "time_reversals", dataset->time_reversals,
                        dataset->n_operations);

    SET_SCALAR_FIELD(plhs[0], 0, "n_atoms", dataset->n_atoms);
    SET_INT_ARRAY_FIELD(plhs[0], 0, "equivalent_atoms",
                        dataset->equivalent_atoms, dataset->n_atoms);

    SET_DOUBLE_MATRIX_FIELD(plhs[0], 0, "transformation_matrix",
                            dataset->transformation_matrix, 3, 3);

    mxArray* origin_shift_mx = mxCreateDoubleMatrix(3, 1, mxREAL);
    double* origin_shift_ptr = mxGetPr(origin_shift_mx);
    for (int i = 0; i < 3; ++i) origin_shift_ptr[i] = dataset->origin_shift[i];
    mxSetField(plhs[0], 0, "origin_shift", origin_shift_mx);

    SET_SCALAR_FIELD(plhs[0], 0, "n_std_atoms", dataset->n_std_atoms);
    SET_DOUBLE_MATRIX_FIELD(plhs[0], 0, "std_lattice", dataset->std_lattice, 3,
                            3);
    SET_INT_ARRAY_FIELD(plhs[0], 0, "std_types", dataset->std_types,
                        dataset->n_std_atoms);
    SET_DOUBLE_2D_ARRAY_FIELD(plhs[0], 0, "std_positions",
                              dataset->std_positions, dataset->n_std_atoms);

    int tensor_elements =
        dataset->n_std_atoms * std::pow(dataset->tensor_rank, 3);
    mxArray* std_tensors_mx = mxCreateDoubleMatrix(tensor_elements, 1, mxREAL);
    double* std_tensors_ptr = mxGetPr(std_tensors_mx);
    for (int i = 0; i < tensor_elements; ++i)
        std_tensors_ptr[i] = dataset->std_tensors[i];
    mxSetField(plhs[0], 0, "std_tensors", std_tensors_mx);

    SET_DOUBLE_MATRIX_FIELD(plhs[0], 0, "std_rotation_matrix",
                            dataset->std_rotation_matrix, 3, 3);

    // 清理 spglib 磁性数据集
    spg_free_magnetic_dataset(dataset);
}

// dataset = symspg('spgat_get_dataset', lattice, position, types, num_atom,
// symprec, angle_tolerance)
void SpglibFunctions::spgat_get_dataset_mex(int nlhs, mxArray* plhs[], int nrhs,
                                            mxArray const* prhs[]) {
    /*
     SpglibDataset *spgat_get_dataset(
         double const lattice[3][3],
         double const position[][3],
         int const types[],
         int const num_atom,
         double const symprec,
         double const angle_tolerance);
    */

    // 验证输入参数数量
    int const expected_number_of_inputs = 6;
    if (nrhs != expected_number_of_inputs) {
        mexErrMsgIdAndTxt("Spglib:invalidNumInputs",
                          "Incorrect number of inputs for spgat_get_dataset.");
    }

    // 提取和验证 lattice 参数
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

    // 提取和验证 position 参数
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

    // 提取和验证 types 参数
    if (mxGetNumberOfElements(prhs[2]) != num_atom) {
        mexErrMsgIdAndTxt("Spglib:invalidTypes",
                          "Types array size must match the number of atoms.");
    }
    int types[num_atom];
    double* types_ptr = mxGetPr(prhs[2]);
    for (int i = 0; i < num_atom; i++) {
        types[i] = static_cast<int>(types_ptr[i]);
    }

    // 提取和验证 symprec 参数
    if (!mxIsDouble(prhs[4]) || mxGetNumberOfElements(prhs[4]) != 1) {
        mexErrMsgIdAndTxt("Spglib:invalidSymprec", "Symprec must be a scalar.");
    }
    double symprec = mxGetScalar(prhs[4]);

    // 提取和验证 angle_tolerance 参数
    if (!mxIsDouble(prhs[5]) || mxGetNumberOfElements(prhs[5]) != 1) {
        mexErrMsgIdAndTxt("Spglib:invalidAngleTolerance",
                          "Angle tolerance must be a scalar.");
    }
    double angle_tolerance = mxGetScalar(prhs[5]);

    // 调用 spgat_get_dataset
    SpglibDataset* dataset = spgat_get_dataset(
        lattice, position, types, num_atom, symprec, angle_tolerance);

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

// dataset = symspg('spg_get_dataset_with_hall_number', lattice, position,
// types, num_atom, hall_number symprec)
void SpglibFunctions::spg_get_dataset_with_hall_number_mex(
    int nlhs, mxArray* plhs[], int nrhs, mxArray const* prhs[]) {
    /*
     SpglibDataset * spg_get_dataset_with_hall_number(const double
     lattice[3][3], const double position[][3], const int types[], const int
     num_atom, const int hall_number, const double symprec)
     */
    // 验证输入参数数量
    int const expected_number_of_inputs = 6;
    if (nrhs != expected_number_of_inputs) {
        mexErrMsgIdAndTxt(
            "Spglib:invalidNumInputs",
            "Incorrect number of inputs for spg_get_dataset_with_hall_number.");
    }

    // 提取和验证 lattice, position, types, hall_number 和 symprec 的代码
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

    int hall_number = static_cast<int>(mxGetScalar(prhs[4]));

    if (!mxIsDouble(prhs[5]) || mxGetNumberOfElements(prhs[5]) != 1) {
        mexErrMsgIdAndTxt("Spglib:invalidSymprec", "Symprec must be a scalar.");
    }
    double symprec = mxGetScalar(prhs[5]);

    // 调用 spg_get_dataset
    SpglibDataset* dataset = spg_get_dataset_with_hall_number(
        lattice, position, types, num_atom, hall_number, symprec);
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

// dataset = symspg('spgat_get_dataset_with_hall_number', lattice, position,
// types, num_atom, hall_number, symprec, angle_tolerance)
void SpglibFunctions::spgat_get_dataset_with_hall_number_mex(
    int nlhs, mxArray* plhs[], int nrhs, mxArray const* prhs[]) {
    /*
     SpglibDataset *spgat_get_dataset_with_hall_number(
         double const lattice[3][3],
         double const position[][3],
         int const types[],
         int const num_atom,
         int const hall_number,
         double const symprec,
         double const angle_tolerance);
    */

    // 验证输入参数数量
    int const expected_number_of_inputs = 7;
    if (nrhs != expected_number_of_inputs) {
        mexErrMsgIdAndTxt("Spglib:invalidNumInputs",
                          "Incorrect number of inputs for "
                          "spgat_get_dataset_with_hall_number.");
    }

    // 提取和验证 lattice 参数
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

    // 提取和验证 position 参数
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

    // 提取和验证 types 参数
    if (mxGetNumberOfElements(prhs[2]) != num_atom) {
        mexErrMsgIdAndTxt("Spglib:invalidTypes",
                          "Types array size must match the number of atoms.");
    }
    int types[num_atom];
    double* types_ptr = mxGetPr(prhs[2]);
    for (int i = 0; i < num_atom; i++) {
        types[i] = static_cast<int>(types_ptr[i]);
    }

    // 提取和验证 hall_number 参数
    int hall_number = static_cast<int>(mxGetScalar(prhs[4]));

    // 提取和验证 symprec 参数
    if (!mxIsDouble(prhs[5]) || mxGetNumberOfElements(prhs[5]) != 1) {
        mexErrMsgIdAndTxt("Spglib:invalidSymprec", "Symprec must be a scalar.");
    }
    double symprec = mxGetScalar(prhs[5]);

    // 提取和验证 angle_tolerance 参数
    if (!mxIsDouble(prhs[6]) || mxGetNumberOfElements(prhs[6]) != 1) {
        mexErrMsgIdAndTxt("Spglib:invalidAngleTolerance",
                          "Angle tolerance must be a scalar.");
    }
    double angle_tolerance = mxGetScalar(prhs[6]);

    // 调用 spgat_get_dataset_with_hall_number
    SpglibDataset* dataset = spgat_get_dataset_with_hall_number(
        lattice, position, types, num_atom, hall_number, symprec,
        angle_tolerance);

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

// [rotations, translations, equivalent_atoms, n_operations] =
// symspg('spg_get_symmetry_with_collinear_spin', max_size, lattice, position,
// types, spins, num_atom, symprec)
void SpglibFunctions::spg_get_symmetry_with_collinear_spin_mex(
    int nlhs, mxArray* plhs[], int nrhs, mxArray const* prhs[]) {
    /*
     int spg_get_symmetry_with_collinear_spin(
         int rotation[][3][3],
         double translation[][3],
         int equivalent_atoms[],
         int const max_size,
         double const lattice[3][3],
         double const position[][3],
         int const types[],
         double const spins[],
         int const num_atom,
         double const symprec);
    */

    // 验证输入参数数量
    int const expected_number_of_inputs = 7;
    if (nrhs != expected_number_of_inputs) {
        mexErrMsgIdAndTxt("Spglib:invalidNumInputs",
                          "Incorrect number of inputs for "
                          "spg_get_symmetry_with_collinear_spin.");
    }

    // 提取和验证 max_size 参数
    int max_size = static_cast<int>(mxGetScalar(prhs[0]));

    // 提取和验证 lattice 参数
    double lattice[3][3];
    if (mxGetM(prhs[1]) != 3 || mxGetN(prhs[1]) != 3) {
        mexErrMsgIdAndTxt("Spglib:invalidLattice",
                          "Lattice matrix must be 3x3.");
    }
    double* lattice_ptr = mxGetPr(prhs[1]);
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            lattice[i][j] = lattice_ptr[i + 3 * j];
        }
    }

    // 提取和验证 position 参数
    int num_atom = mxGetM(prhs[2]);
    double position[num_atom][3];
    if (mxGetN(prhs[2]) != 3) {
        mexErrMsgIdAndTxt("Spglib:invalidPosition",
                          "Position array must be Nx3.");
    }
    double* position_ptr = mxGetPr(prhs[2]);
    for (int i = 0; i < num_atom; i++) {
        for (int j = 0; j < 3; j++) {
            position[i][j] = position_ptr[i + num_atom * j];
        }
    }

    // 提取和验证 types 参数
    if (mxGetNumberOfElements(prhs[3]) != num_atom) {
        mexErrMsgIdAndTxt("Spglib:invalidTypes",
                          "Types array size must match the number of atoms.");
    }
    int types[num_atom];
    double* types_ptr = mxGetPr(prhs[3]);
    for (int i = 0; i < num_atom; i++) {
        types[i] = static_cast<int>(types_ptr[i]);
    }

    // 提取和验证 spins 参数
    if (mxGetNumberOfElements(prhs[4]) != num_atom) {
        mexErrMsgIdAndTxt("Spglib:invalidSpins",
                          "Spins array size must match the number of atoms.");
    }
    double spins[num_atom];
    double* spins_ptr = mxGetPr(prhs[4]);
    for (int i = 0; i < num_atom; i++) {
        spins[i] = spins_ptr[i];
    }

    // 提取和验证 symprec 参数
    if (!mxIsDouble(prhs[6]) || mxGetNumberOfElements(prhs[6]) != 1) {
        mexErrMsgIdAndTxt("Spglib:invalidSymprec", "Symprec must be a scalar.");
    }
    double symprec = mxGetScalar(prhs[6]);

    // 初始化 rotation, translation, 和 equivalent_atoms 数组
    int rotation[max_size][3][3];
    double translation[max_size][3];
    int equivalent_atoms[num_atom];

    // 调用 spg_get_symmetry_with_collinear_spin
    int n_operations = spg_get_symmetry_with_collinear_spin(
        rotation, translation, equivalent_atoms, max_size, lattice, position,
        types, spins, num_atom, symprec);

    if (n_operations == 0) {
        mexErrMsgIdAndTxt(
            "Spglib:symmetryError",
            "Failed to get symmetry operations with collinear spin.");
    }

    // 创建输出数组并设置字段
    // 输出旋转矩阵 (Nx3x3 int 数组)
    mwSize dims[3] = {static_cast<mwSize>(n_operations), 3, 3};
    plhs[0] = mxCreateNumericArray(3, dims, mxINT32_CLASS, mxREAL);
    auto* rotations_out = static_cast<int32_t*>(mxGetData(plhs[0]));

    for (int k = 0; k < n_operations; ++k) {
        for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < 3; ++j) {
                rotations_out[k + i * n_operations + j * n_operations * 3] =
                    rotation[k][i][j];
            }
        }
    }

    // 输出平移数组 (Nx3 double 数组)
    plhs[1] = mxCreateDoubleMatrix(n_operations, 3, mxREAL);
    double* translations_out = mxGetPr(plhs[1]);
    for (int i = 0; i < n_operations; ++i) {
        for (int j = 0; j < 3; ++j) {
            translations_out[i + j * n_operations] = translation[i][j];
        }
    }

    // 输出等效原子数组 (num_atom int 数组)
    plhs[2] = mxCreateNumericMatrix(num_atom, 1, mxINT32_CLASS, mxREAL);
    int32_t* equivalent_atoms_out = static_cast<int32_t*>(mxGetData(plhs[2]));
    for (int i = 0; i < num_atom; ++i) {
        equivalent_atoms_out[i] = equivalent_atoms[i];
    }

    // 输出操作数量
    plhs[3] = mxCreateDoubleScalar(static_cast<double>(n_operations));
}

// [rotations, translations, equivalent_atoms, n_operations] =
// symspg('spgat_get_symmetry_with_collinear_spin', max_size, lattice, position,
// types, spins, num_atom, symprec, angle_tolerance)
void SpglibFunctions::spgat_get_symmetry_with_collinear_spin_mex(
    int nlhs, mxArray* plhs[], int nrhs, mxArray const* prhs[]) {
    /*
     int spgat_get_symmetry_with_collinear_spin(
         int rotation[][3][3],
         double translation[][3],
         int equivalent_atoms[],
         int const max_size,
         double const lattice[3][3],
         double const position[][3],
         int const types[],
         double const spins[],
         int const num_atom,
         double const symprec,
         double const angle_tolerance);
    */

    // 验证输入参数数量
    int const expected_number_of_inputs = 8;
    if (nrhs != expected_number_of_inputs) {
        mexErrMsgIdAndTxt("Spglib:invalidNumInputs",
                          "Incorrect number of inputs for "
                          "spgat_get_symmetry_with_collinear_spin.");
    }

    // 提取和验证 max_size 参数
    int max_size = static_cast<int>(mxGetScalar(prhs[0]));

    // 提取和验证 lattice 参数
    double lattice[3][3];
    if (mxGetM(prhs[1]) != 3 || mxGetN(prhs[1]) != 3) {
        mexErrMsgIdAndTxt("Spglib:invalidLattice",
                          "Lattice matrix must be 3x3.");
    }
    double* lattice_ptr = mxGetPr(prhs[1]);
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            lattice[i][j] = lattice_ptr[i + 3 * j];
        }
    }

    // 提取和验证 position 参数
    int num_atom = mxGetM(prhs[2]);
    double position[num_atom][3];
    if (mxGetN(prhs[2]) != 3) {
        mexErrMsgIdAndTxt("Spglib:invalidPosition",
                          "Position array must be Nx3.");
    }
    double* position_ptr = mxGetPr(prhs[2]);
    for (int i = 0; i < num_atom; i++) {
        for (int j = 0; j < 3; j++) {
            position[i][j] = position_ptr[i + num_atom * j];
        }
    }

    // 提取和验证 types 参数
    if (mxGetNumberOfElements(prhs[3]) != num_atom) {
        mexErrMsgIdAndTxt("Spglib:invalidTypes",
                          "Types array size must match the number of atoms.");
    }
    int types[num_atom];
    double* types_ptr = mxGetPr(prhs[3]);
    for (int i = 0; i < num_atom; i++) {
        types[i] = static_cast<int>(types_ptr[i]);
    }

    // 提取和验证 spins 参数
    if (mxGetNumberOfElements(prhs[4]) != num_atom) {
        mexErrMsgIdAndTxt("Spglib:invalidSpins",
                          "Spins array size must match the number of atoms.");
    }
    double spins[num_atom];
    double* spins_ptr = mxGetPr(prhs[4]);
    for (int i = 0; i < num_atom; i++) {
        spins[i] = spins_ptr[i];
    }

    // 提取和验证 symprec 参数
    if (!mxIsDouble(prhs[6]) || mxGetNumberOfElements(prhs[6]) != 1) {
        mexErrMsgIdAndTxt("Spglib:invalidSymprec", "Symprec must be a scalar.");
    }
    double symprec = mxGetScalar(prhs[6]);

    // 提取和验证 angle_tolerance 参数
    if (!mxIsDouble(prhs[7]) || mxGetNumberOfElements(prhs[7]) != 1) {
        mexErrMsgIdAndTxt("Spglib:invalidAngleTolerance",
                          "Angle tolerance must be a scalar.");
    }
    double angle_tolerance = mxGetScalar(prhs[7]);

    // 初始化 rotation, translation, 和 equivalent_atoms 数组
    int rotation[max_size][3][3];
    double translation[max_size][3];
    int equivalent_atoms[num_atom];

    // 调用 spgat_get_symmetry_with_collinear_spin
    int n_operations = spgat_get_symmetry_with_collinear_spin(
        rotation, translation, equivalent_atoms, max_size, lattice, position,
        types, spins, num_atom, symprec, angle_tolerance);

    if (n_operations == 0) {
        mexErrMsgIdAndTxt(
            "Spglib:symmetryError",
            "Failed to get symmetry operations with collinear spin.");
    }

    // 创建输出数组并设置字段
    // 输出旋转矩阵 (Nx3x3 int 数组)
    mwSize dims[3] = {static_cast<mwSize>(n_operations), 3, 3};
    plhs[0] = mxCreateNumericArray(3, dims, mxINT32_CLASS, mxREAL);
    auto* rotations_out = static_cast<int32_t*>(mxGetData(plhs[0]));

    for (int k = 0; k < n_operations; ++k) {
        for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < 3; ++j) {
                rotations_out[k + i * n_operations + j * n_operations * 3] =
                    rotation[k][i][j];
            }
        }
    }

    // 输出平移数组 (Nx3 double 数组)
    plhs[1] = mxCreateDoubleMatrix(n_operations, 3, mxREAL);
    double* translations_out = mxGetPr(plhs[1]);
    for (int i = 0; i < n_operations; ++i) {
        for (int j = 0; j < 3; ++j) {
            translations_out[i + j * n_operations] = translation[i][j];
        }
    }

    // 输出等效原子数组 (num_atom int 数组)
    plhs[2] = mxCreateNumericMatrix(num_atom, 1, mxINT32_CLASS, mxREAL);
    int32_t* equivalent_atoms_out = static_cast<int32_t*>(mxGetData(plhs[2]));
    for (int i = 0; i < num_atom; ++i) {
        equivalent_atoms_out[i] = equivalent_atoms[i];
    }

    // 输出操作数量
    plhs[3] = mxCreateDoubleScalar(static_cast<double>(n_operations));
}
// [rotations, translations, equivalent_atoms, n_operations] =
// symspg('spgms_get_symmetry_with_collinear_spin', max_size, lattice, position,
// types, spins, num_atom, symprec, angle_tolerance, mag_symprec)
void SpglibFunctions::spgms_get_symmetry_with_collinear_spin_mex(
    int nlhs, mxArray* plhs[], int nrhs, mxArray const* prhs[]) {
    /*
     int spgms_get_symmetry_with_collinear_spin(
         int rotation[][3][3],
         double translation[][3],
         int equivalent_atoms[],
         int const max_size,
         double const lattice[3][3],
         double const position[][3],
         int const types[],
         double const spins[],
         int const num_atom,
         double const symprec,
         double const angle_tolerance,
         double const mag_symprec);
    */

    // 验证输入参数数量
    int const expected_number_of_inputs = 9;
    if (nrhs != expected_number_of_inputs) {
        mexErrMsgIdAndTxt("Spglib:invalidNumInputs",
                          "Incorrect number of inputs for "
                          "spgms_get_symmetry_with_collinear_spin.");
    }

    // 提取和验证 max_size 参数
    int max_size = static_cast<int>(mxGetScalar(prhs[0]));

    // 提取和验证 lattice 参数
    double lattice[3][3];
    if (mxGetM(prhs[1]) != 3 || mxGetN(prhs[1]) != 3) {
        mexErrMsgIdAndTxt("Spglib:invalidLattice",
                          "Lattice matrix must be 3x3.");
    }
    double* lattice_ptr = mxGetPr(prhs[1]);
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            lattice[i][j] = lattice_ptr[i + 3 * j];
        }
    }

    // 提取和验证 position 参数
    int num_atom = mxGetM(prhs[2]);
    double position[num_atom][3];
    if (mxGetN(prhs[2]) != 3) {
        mexErrMsgIdAndTxt("Spglib:invalidPosition",
                          "Position array must be Nx3.");
    }
    double* position_ptr = mxGetPr(prhs[2]);
    for (int i = 0; i < num_atom; i++) {
        for (int j = 0; j < 3; j++) {
            position[i][j] = position_ptr[i + num_atom * j];
        }
    }

    // 提取和验证 types 参数
    if (mxGetNumberOfElements(prhs[3]) != num_atom) {
        mexErrMsgIdAndTxt("Spglib:invalidTypes",
                          "Types array size must match the number of atoms.");
    }
    int types[num_atom];
    double* types_ptr = mxGetPr(prhs[3]);
    for (int i = 0; i < num_atom; i++) {
        types[i] = static_cast<int>(types_ptr[i]);
    }

    // 提取和验证 spins 参数
    if (mxGetNumberOfElements(prhs[4]) != num_atom) {
        mexErrMsgIdAndTxt("Spglib:invalidSpins",
                          "Spins array size must match the number of atoms.");
    }
    double spins[num_atom];
    double* spins_ptr = mxGetPr(prhs[4]);
    for (int i = 0; i < num_atom; i++) {
        spins[i] = spins_ptr[i];
    }

    // 提取和验证 symprec 参数
    if (!mxIsDouble(prhs[6]) || mxGetNumberOfElements(prhs[6]) != 1) {
        mexErrMsgIdAndTxt("Spglib:invalidSymprec", "Symprec must be a scalar.");
    }
    double symprec = mxGetScalar(prhs[6]);

    // 提取和验证 angle_tolerance 参数
    if (!mxIsDouble(prhs[7]) || mxGetNumberOfElements(prhs[7]) != 1) {
        mexErrMsgIdAndTxt("Spglib:invalidAngleTolerance",
                          "Angle tolerance must be a scalar.");
    }
    double angle_tolerance = mxGetScalar(prhs[7]);

    // 提取和验证 mag_symprec 参数
    if (!mxIsDouble(prhs[8]) || mxGetNumberOfElements(prhs[8]) != 1) {
        mexErrMsgIdAndTxt("Spglib:invalidMagSymprec",
                          "Magnetic symmetry precision must be a scalar.");
    }
    double mag_symprec = mxGetScalar(prhs[8]);

    // 初始化 rotation, translation, 和 equivalent_atoms 数组
    int rotation[max_size][3][3];
    double translation[max_size][3];
    int equivalent_atoms[num_atom];

    // 调用 spgms_get_symmetry_with_collinear_spin
    int n_operations = spgms_get_symmetry_with_collinear_spin(
        rotation, translation, equivalent_atoms, max_size, lattice, position,
        types, spins, num_atom, symprec, angle_tolerance, mag_symprec);

    if (n_operations == 0) {
        mexErrMsgIdAndTxt(
            "Spglib:symmetryError",
            "Failed to get symmetry operations with collinear spin.");
    }

    // 创建输出数组并设置字段
    // 输出旋转矩阵 (Nx3x3 int 数组)
    mwSize dims[3] = {static_cast<mwSize>(n_operations), 3, 3};
    plhs[0] = mxCreateNumericArray(3, dims, mxINT32_CLASS, mxREAL);
    auto* rotations_out = static_cast<int32_t*>(mxGetData(plhs[0]));

    for (int k = 0; k < n_operations; ++k) {
        for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < 3; ++j) {
                rotations_out[k + i * n_operations + j * n_operations * 3] =
                    rotation[k][i][j];
            }
        }
    }

    // 输出平移数组 (Nx3 double 数组)
    plhs[1] = mxCreateDoubleMatrix(n_operations, 3, mxREAL);
    double* translations_out = mxGetPr(plhs[1]);
    for (int i = 0; i < n_operations; ++i) {
        for (int j = 0; j < 3; ++j) {
            translations_out[i + j * n_operations] = translation[i][j];
        }
    }

    // 输出等效原子数组 (num_atom int 数组)
    plhs[2] = mxCreateNumericMatrix(num_atom, 1, mxINT32_CLASS, mxREAL);
    int32_t* equivalent_atoms_out = static_cast<int32_t*>(mxGetData(plhs[2]));
    for (int i = 0; i < num_atom; ++i) {
        equivalent_atoms_out[i] = equivalent_atoms[i];
    }

    // 输出操作数量
    plhs[3] = mxCreateDoubleScalar(static_cast<double>(n_operations));
}

// [rotations, translations, equivalent_atoms, primitive_lattice, spin_flips,
// n_operations] = symspg('spg_get_symmetry_with_site_tensors', max_size,
// lattice, position, types, tensors, tensor_rank, num_atom, with_time_reversal,
// is_axial, symprec)
void SpglibFunctions::spg_get_symmetry_with_site_tensors_mex(
    int nlhs, mxArray* plhs[], int nrhs, mxArray const* prhs[]) {
    /*
     int spg_get_symmetry_with_site_tensors(
         int rotation[][3][3],
         double translation[][3],
         int equivalent_atoms[],
         double primitive_lattice[3][3],
         int *spin_flips,
         int const max_size,
         double const lattice[3][3],
         double const position[][3],
         int const types[],
         double const *tensors,
         int const tensor_rank,
         int const num_atom,
         int const with_time_reversal,
         int const is_axial,
         double const symprec);
    */

    // 验证输入参数数量
    int const expected_number_of_inputs =
        10;  // Updated expected number of inputs
    if (nrhs != expected_number_of_inputs) {
        mexErrMsgIdAndTxt("Spglib:invalidNumInputs",
                          "Incorrect number of inputs for "
                          "spg_get_symmetry_with_site_tensors.");
    }

    // 提取和验证 max_size 参数
    int max_size = static_cast<int>(mxGetScalar(prhs[0]));

    // 提取和验证 lattice 参数
    double lattice[3][3];
    if (mxGetM(prhs[1]) != 3 || mxGetN(prhs[1]) != 3) {
        mexErrMsgIdAndTxt("Spglib:invalidLattice",
                          "Lattice matrix must be 3x3.");
    }
    double* lattice_ptr = mxGetPr(prhs[1]);
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            lattice[i][j] = lattice_ptr[i + 3 * j];
        }
    }

    // 提取和验证 position 参数
    int num_atom = mxGetM(prhs[2]);
    double position[num_atom][3];
    if (mxGetN(prhs[2]) != 3) {
        mexErrMsgIdAndTxt("Spglib:invalidPosition",
                          "Position array must be Nx3.");
    }
    double* position_ptr = mxGetPr(prhs[2]);
    for (int i = 0; i < num_atom; i++) {
        for (int j = 0; j < 3; j++) {
            position[i][j] = position_ptr[i + num_atom * j];
        }
    }

    // 提取和验证 types 参数
    if (mxGetNumberOfElements(prhs[3]) != num_atom) {
        mexErrMsgIdAndTxt("Spglib:invalidTypes",
                          "Types array size must match the number of atoms.");
    }
    int types[num_atom];
    double* types_ptr = mxGetPr(prhs[3]);
    for (int i = 0; i < num_atom; i++) {
        types[i] = static_cast<int>(types_ptr[i]);
    }

    // 提取和验证 tensors 参数
    double* tensors = mxGetPr(prhs[4]);

    // 提取和验证 tensor_rank 参数
    int tensor_rank = static_cast<int>(mxGetScalar(prhs[5]));

    // 提取和验证 with_time_reversal 参数
    int with_time_reversal = static_cast<int>(mxGetScalar(prhs[7]));

    // 提取和验证 is_axial 参数
    int is_axial = static_cast<int>(mxGetScalar(prhs[8]));

    // 提取和验证 symprec 参数
    if (!mxIsDouble(prhs[9]) || mxGetNumberOfElements(prhs[9]) != 1) {
        mexErrMsgIdAndTxt("Spglib:invalidSymprec", "Symprec must be a scalar.");
    }
    double symprec = mxGetScalar(prhs[9]);

    // 初始化 rotation, translation, equivalent_atoms, 和 primitive_lattice 数组
    int rotation[max_size][3][3];
    double translation[max_size][3];
    int equivalent_atoms[num_atom];
    double primitive_lattice[3][3];
    int spin_flips[max_size];

    // 调用 spg_get_symmetry_with_site_tensors
    int n_operations = spg_get_symmetry_with_site_tensors(
        rotation, translation, equivalent_atoms, primitive_lattice, spin_flips,
        max_size, lattice, position, types, tensors, tensor_rank, num_atom,
        with_time_reversal, is_axial, symprec);

    if (n_operations == 0) {
        mexErrMsgIdAndTxt(
            "Spglib:symmetryError",
            "Failed to get symmetry operations with site tensors.");
    }

    // 创建输出数组并设置字段
    // 输出旋转矩阵 (Nx3x3 int 数组)
    mwSize dims[3] = {static_cast<mwSize>(n_operations), 3, 3};
    plhs[0] = mxCreateNumericArray(3, dims, mxINT32_CLASS, mxREAL);
    auto* rotations_out = static_cast<int32_t*>(mxGetData(plhs[0]));

    for (int k = 0; k < n_operations; ++k) {
        for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < 3; ++j) {
                rotations_out[k + i * n_operations + j * n_operations * 3] =
                    rotation[k][i][j];
            }
        }
    }

    // 输出平移数组 (Nx3 double 数组)
    plhs[1] = mxCreateDoubleMatrix(n_operations, 3, mxREAL);
    double* translations_out = mxGetPr(plhs[1]);
    for (int i = 0; i < n_operations; ++i) {
        for (int j = 0; j < 3; ++j) {
            translations_out[i + j * n_operations] = translation[i][j];
        }
    }

    // 输出等效原子数组 (num_atom int 数组)
    plhs[2] = mxCreateNumericMatrix(num_atom, 1, mxINT32_CLASS, mxREAL);
    auto* equivalent_atoms_out = static_cast<int32_t*>(mxGetData(plhs[2]));
    for (int i = 0; i < num_atom; ++i) {
        equivalent_atoms_out[i] = equivalent_atoms[i];
    }

    // 输出 primitive_lattice (3x3 double 数组)
    plhs[3] = mxCreateDoubleMatrix(3, 3, mxREAL);
    double* primitive_lattice_out = mxGetPr(plhs[3]);
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            primitive_lattice_out[i + j * 3] = primitive_lattice[i][j];
        }
    }

    // 输出 spin_flips (num_atom int 数组)
    plhs[4] = mxCreateNumericMatrix(n_operations, 1, mxINT32_CLASS, mxREAL);
    auto* spin_flips_out = static_cast<int32_t*>(mxGetData(plhs[4]));
    for (int i = 0; i < n_operations; ++i) {
        spin_flips_out[i] = static_cast<int32_t>(spin_flips[i]);
    }

    // 输出操作数量
    plhs[5] = mxCreateDoubleScalar(static_cast<double>(n_operations));
}

// [rotations, translations, equivalent_atoms, primitive_lattice, spin_flips,
// n_operations] = symspg('spgat_get_symmetry_with_site_tensors', max_size,
// lattice, position, types, tensors, tensor_rank, num_atom, with_time_reversal,
// is_axial, symprec, angle_tolerance)
void SpglibFunctions::spgat_get_symmetry_with_site_tensors_mex(
    int nlhs, mxArray* plhs[], int nrhs, mxArray const* prhs[]) {
    /*
     int spgat_get_symmetry_with_site_tensors(
         int rotation[][3][3],
         double translation[][3],
         int equivalent_atoms[],
         double primitive_lattice[3][3],
         int *spin_flips,
         int const max_size,
         double const lattice[3][3],
         double const position[][3],
         int const types[],
         double const *tensors,
         int const tensor_rank,
         int const num_atom,
         int const with_time_reversal,
         int const is_axial,
         double const symprec,
         double const angle_tolerance);
    */

    // 验证输入参数数量
    int const expected_number_of_inputs =
        11;  // Adjusted to match inputs without spin_flips
    if (nrhs != expected_number_of_inputs) {
        mexErrMsgIdAndTxt("Spglib:invalidNumInputs",
                          "Incorrect number of inputs for "
                          "spgat_get_symmetry_with_site_tensors.");
    }

    // 提取和验证 max_size 参数
    int max_size = static_cast<int>(mxGetScalar(prhs[0]));

    // 提取和验证 lattice 参数
    double lattice[3][3];
    if (mxGetM(prhs[1]) != 3 || mxGetN(prhs[1]) != 3) {
        mexErrMsgIdAndTxt("Spglib:invalidLattice",
                          "Lattice matrix must be 3x3.");
    }
    double* lattice_ptr = mxGetPr(prhs[1]);
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            lattice[i][j] = lattice_ptr[i + 3 * j];
        }
    }

    // 提取和验证 position 参数
    int num_atom = mxGetM(prhs[2]);
    double position[num_atom][3];
    if (mxGetN(prhs[2]) != 3) {
        mexErrMsgIdAndTxt("Spglib:invalidPosition",
                          "Position array must be Nx3.");
    }
    double* position_ptr = mxGetPr(prhs[2]);
    for (int i = 0; i < num_atom; i++) {
        for (int j = 0; j < 3; j++) {
            position[i][j] = position_ptr[i + num_atom * j];
        }
    }

    // 提取和验证 types 参数
    if (mxGetNumberOfElements(prhs[3]) != num_atom) {
        mexErrMsgIdAndTxt("Spglib:invalidTypes",
                          "Types array size must match the number of atoms.");
    }
    int types[num_atom];
    double* types_ptr = mxGetPr(prhs[3]);
    for (int i = 0; i < num_atom; i++) {
        types[i] = static_cast<int>(types_ptr[i]);
    }

    // 提取和验证 tensors 参数
    double* tensors = mxGetPr(prhs[4]);

    // 提取和验证 tensor_rank 参数
    int tensor_rank = static_cast<int>(mxGetScalar(prhs[5]));

    // 提取和验证 with_time_reversal 参数
    int with_time_reversal = static_cast<int>(mxGetScalar(prhs[7]));

    // 提取和验证 is_axial 参数
    int is_axial = static_cast<int>(mxGetScalar(prhs[8]));

    // 提取和验证 symprec 参数
    if (!mxIsDouble(prhs[9]) || mxGetNumberOfElements(prhs[9]) != 1) {
        mexErrMsgIdAndTxt("Spglib:invalidSymprec", "Symprec must be a scalar.");
    }
    double symprec = mxGetScalar(prhs[9]);

    // 提取和验证 angle_tolerance 参数
    if (!mxIsDouble(prhs[10]) || mxGetNumberOfElements(prhs[10]) != 1) {
        mexErrMsgIdAndTxt("Spglib:invalidAngleTolerance",
                          "Angle tolerance must be a scalar.");
    }
    double angle_tolerance = mxGetScalar(prhs[9]);

    // 初始化 rotation, translation, equivalent_atoms, 和 primitive_lattice 数组
    int rotation[max_size][3][3];
    double translation[max_size][3];
    int equivalent_atoms[num_atom];
    double primitive_lattice[3][3];
    int spin_flips[num_atom];  // Declare spin_flips as an output variable

    // 调用 spgat_get_symmetry_with_site_tensors
    int n_operations = spgat_get_symmetry_with_site_tensors(
        rotation, translation, equivalent_atoms, primitive_lattice, spin_flips,
        max_size, lattice, position, types, tensors, tensor_rank, num_atom,
        with_time_reversal, is_axial, symprec, angle_tolerance);

    if (n_operations == 0) {
        mexErrMsgIdAndTxt(
            "Spglib:symmetryError",
            "Failed to get symmetry operations with site tensors.");
    }

    // 创建输出数组并设置字段
    // 输出旋转矩阵 (Nx3x3 int 数组)
    mwSize dims[3] = {static_cast<mwSize>(n_operations), 3, 3};
    plhs[0] = mxCreateNumericArray(3, dims, mxINT32_CLASS, mxREAL);
    auto* rotations_out = static_cast<int32_t*>(mxGetData(plhs[0]));

    for (int k = 0; k < n_operations; ++k) {
        for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < 3; ++j) {
                rotations_out[k + i * n_operations + j * n_operations * 3] =
                    rotation[k][i][j];
            }
        }
    }

    // 输出平移数组 (Nx3 double 数组)
    plhs[1] = mxCreateDoubleMatrix(n_operations, 3, mxREAL);
    double* translations_out = mxGetPr(plhs[1]);
    for (int i = 0; i < n_operations; ++i) {
        for (int j = 0; j < 3; ++j) {
            translations_out[i + j * n_operations] = translation[i][j];
        }
    }

    // 输出等效原子数组 (num_atom int 数组)
    plhs[2] = mxCreateNumericMatrix(num_atom, 1, mxINT32_CLASS, mxREAL);
    auto* equivalent_atoms_out = static_cast<int32_t*>(mxGetData(plhs[2]));
    for (int i = 0; i < num_atom; ++i) {
        equivalent_atoms_out[i] = equivalent_atoms[i];
    }

    // 输出 primitive_lattice (3x3 double 数组)
    plhs[3] = mxCreateDoubleMatrix(3, 3, mxREAL);
    double* primitive_lattice_out = mxGetPr(plhs[3]);
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            primitive_lattice_out[i + j * 3] = primitive_lattice[i][j];
        }
    }

    // 输出 spin_flips (num_atom int 数组)
    plhs[4] = mxCreateNumericMatrix(num_atom, 1, mxINT32_CLASS, mxREAL);
    auto* spin_flips_out = static_cast<int32_t*>(mxGetData(plhs[4]));
    for (int i = 0; i < num_atom; ++i) {
        spin_flips_out[i] = spin_flips[i];
    }

    // 输出操作数量
    plhs[5] = mxCreateDoubleScalar(static_cast<double>(n_operations));
}

// [rotations, translations, equivalent_atoms, primitive_lattice, spin_flips,
// n_operations] = symspg('spgms_get_symmetry_with_site_tensors', max_size,
// lattice, position, types, tensors, tensor_rank, num_atom, with_time_reversal,
// is_axial, symprec, angle_tolerance, mag_symprec)
void SpglibFunctions::spgms_get_symmetry_with_site_tensors_mex(
    int nlhs, mxArray* plhs[], int nrhs, mxArray const* prhs[]) {
    /*
     int spgms_get_symmetry_with_site_tensors(
         int rotation[][3][3],
         double translation[][3],
         int equivalent_atoms[],
         double primitive_lattice[3][3],
         int *spin_flips,
         int const max_size,
         double const lattice[3][3],
         double const position[][3],
         int const types[],
         double const *tensors,
         int const tensor_rank,
         int const num_atom,
         int const with_time_reversal,
         int const is_axial,
         double const symprec,
         double const angle_tolerance,
         double const mag_symprec);
    */

    // 验证输入参数数量
    int const expected_number_of_inputs =
        12;  // Adjusted to match inputs without spin_flips
    if (nrhs != expected_number_of_inputs) {
        mexErrMsgIdAndTxt("Spglib:invalidNumInputs",
                          "Incorrect number of inputs for "
                          "spgms_get_symmetry_with_site_tensors.");
    }

    // 提取和验证 max_size 参数
    int max_size = static_cast<int>(mxGetScalar(prhs[0]));

    // 提取和验证 lattice 参数
    double lattice[3][3];
    if (mxGetM(prhs[1]) != 3 || mxGetN(prhs[1]) != 3) {
        mexErrMsgIdAndTxt("Spglib:invalidLattice",
                          "Lattice matrix must be 3x3.");
    }
    double* lattice_ptr = mxGetPr(prhs[1]);
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            lattice[i][j] = lattice_ptr[i + 3 * j];
        }
    }

    // 提取和验证 position 参数
    int num_atom = mxGetM(prhs[2]);
    double position[num_atom][3];
    if (mxGetN(prhs[2]) != 3) {
        mexErrMsgIdAndTxt("Spglib:invalidPosition",
                          "Position array must be Nx3.");
    }
    double* position_ptr = mxGetPr(prhs[2]);
    for (int i = 0; i < num_atom; i++) {
        for (int j = 0; j < 3; j++) {
            position[i][j] = position_ptr[i + num_atom * j];
        }
    }

    // 提取和验证 types 参数
    if (mxGetNumberOfElements(prhs[3]) != num_atom) {
        mexErrMsgIdAndTxt("Spglib:invalidTypes",
                          "Types array size must match the number of atoms.");
    }
    int types[num_atom];
    double* types_ptr = mxGetPr(prhs[3]);
    for (int i = 0; i < num_atom; i++) {
        types[i] = static_cast<int>(types_ptr[i]);
    }

    // 提取和验证 tensors 参数
    double* tensors = mxGetPr(prhs[4]);

    // 提取和验证 tensor_rank 参数
    int tensor_rank = static_cast<int>(mxGetScalar(prhs[5]));

    // 提取和验证 with_time_reversal 参数
    int with_time_reversal = static_cast<int>(mxGetScalar(prhs[7]));

    // 提取和验证 is_axial 参数
    int is_axial = static_cast<int>(mxGetScalar(prhs[8]));

    // 提取和验证 symprec 参数
    if (!mxIsDouble(prhs[9]) || mxGetNumberOfElements(prhs[9]) != 1) {
        mexErrMsgIdAndTxt("Spglib:invalidSymprec", "Symprec must be a scalar.");
    }
    double symprec = mxGetScalar(prhs[9]);

    // 提取和验证 angle_tolerance 参数
    if (!mxIsDouble(prhs[10]) || mxGetNumberOfElements(prhs[10]) != 1) {
        mexErrMsgIdAndTxt("Spglib:invalidAngleTolerance",
                          "Angle tolerance must be a scalar.");
    }
    double angle_tolerance = mxGetScalar(prhs[10]);

    // 提取和验证 mag_symprec 参数
    if (!mxIsDouble(prhs[11]) || mxGetNumberOfElements(prhs[11]) != 1) {
        mexErrMsgIdAndTxt("Spglib:invalidMagSymprec",
                          "Magnetic symmetry precision must be a scalar.");
    }
    double mag_symprec = mxGetScalar(prhs[11]);

    // 初始化 rotation, translation, equivalent_atoms, 和 primitive_lattice 数组
    int rotation[max_size][3][3];
    double translation[max_size][3];
    int equivalent_atoms[num_atom];
    double primitive_lattice[3][3];
    int spin_flips[num_atom];  // Declare spin_flips as an output variable

    // 调用 spgms_get_symmetry_with_site_tensors
    int n_operations = spgms_get_symmetry_with_site_tensors(
        rotation, translation, equivalent_atoms, primitive_lattice, spin_flips,
        max_size, lattice, position, types, tensors, tensor_rank, num_atom,
        with_time_reversal, is_axial, symprec, angle_tolerance, mag_symprec);

    if (n_operations == 0) {
        mexErrMsgIdAndTxt(
            "Spglib:symmetryError",
            "Failed to get symmetry operations with site tensors.");
    }

    // 创建输出数组并设置字段
    // 输出旋转矩阵 (Nx3x3 int 数组)
    mwSize dims[3] = {static_cast<mwSize>(n_operations), 3, 3};
    plhs[0] = mxCreateNumericArray(3, dims, mxINT32_CLASS, mxREAL);
    auto* rotations_out = static_cast<int32_t*>(mxGetData(plhs[0]));

    for (int k = 0; k < n_operations; ++k) {
        for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < 3; ++j) {
                rotations_out[k + i * n_operations + j * n_operations * 3] =
                    rotation[k][i][j];
            }
        }
    }

    // 输出平移数组 (Nx3 double 数组)
    plhs[1] = mxCreateDoubleMatrix(n_operations, 3, mxREAL);
    double* translations_out = mxGetPr(plhs[1]);
    for (int i = 0; i < n_operations; ++i) {
        for (int j = 0; j < 3; ++j) {
            translations_out[i + j * n_operations] = translation[i][j];
        }
    }

    // 输出等效原子数组 (num_atom int 数组)
    plhs[2] = mxCreateNumericMatrix(num_atom, 1, mxINT32_CLASS, mxREAL);
    int32_t* equivalent_atoms_out = static_cast<int32_t*>(mxGetData(plhs[2]));
    for (int i = 0; i < num_atom; ++i) {
        equivalent_atoms_out[i] = equivalent_atoms[i];
    }

    // 输出 primitive_lattice (3x3 double 数组)
    plhs[3] = mxCreateDoubleMatrix(3, 3, mxREAL);
    double* primitive_lattice_out = mxGetPr(plhs[3]);
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            primitive_lattice_out[i + j * 3] = primitive_lattice[i][j];
        }
    }

    // 输出 spin_flips (num_atom int 数组)
    plhs[4] = mxCreateNumericMatrix(num_atom, 1, mxINT32_CLASS, mxREAL);
    auto* spin_flips_out = static_cast<int32_t*>(mxGetData(plhs[4]));
    for (int i = 0; i < num_atom; ++i) {
        spin_flips_out[i] = spin_flips[i];
    }

    // 输出操作数量
    plhs[5] = mxCreateDoubleScalar(static_cast<double>(n_operations));
}

// spacegroup_type = symspg('spg_get_spacegroup_type_from_symmetry', rotation,
// translation, num_operations, lattice, symprec)
void SpglibFunctions::spg_get_spacegroup_type_from_symmetry_mex(
    int nlhs, mxArray* plhs[], int nrhs, mxArray const* prhs[]) {
    /*
     SpglibSpacegroupType spg_get_spacegroup_type_from_symmetry(
         int const rotation[][3][3],
         double const translation[][3],
         int const num_operations,
         double const lattice[3][3],
         double const symprec);
    */

    // 验证输入参数数量
    int const expected_number_of_inputs = 5;
    if (nrhs != expected_number_of_inputs) {
        mexErrMsgIdAndTxt("Spglib:invalidNumInputs",
                          "Incorrect number of inputs for "
                          "spg_get_spacegroup_type_from_symmetry.");
    }

    // 提取和验证 rotation 参数
    mwSize const* dims = mxGetDimensions(prhs[0]);
    int num_operations = dims[0];  // N 的值
    if (dims[1] != 3 || dims[2] != 3) {
        mexErrMsgIdAndTxt("Spglib:invalidRotation",
                          "Rotation matrix must have dimensions Nx3x3.");
    }

    int rotation[num_operations][3][3];
    auto* rotation_ptr = static_cast<int32_t*>(mxGetData(prhs[0]));

    for (int k = 0; k < num_operations; k++) {
        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++) {
                rotation[k][i][j] = rotation_ptr[k + i * num_operations +
                                                 j * num_operations * 3];
            }
        }
    }

    // 提取和验证 translation 参数
    double translation[num_operations][3];
    if (mxGetM(prhs[1]) != num_operations || mxGetN(prhs[1]) != 3) {
        mexErrMsgIdAndTxt("Spglib:invalidTranslation",
                          "Translation array must be Nx3.");
    }
    double* translation_ptr = mxGetPr(prhs[1]);
    for (int i = 0; i < num_operations; i++) {
        for (int j = 0; j < 3; j++) {
            translation[i][j] = translation_ptr[i + j * num_operations];
        }
    }

    // 提取和验证 lattice 参数
    double lattice[3][3];
    if (mxGetM(prhs[3]) != 3 || mxGetN(prhs[3]) != 3) {
        mexErrMsgIdAndTxt("Spglib:invalidLattice",
                          "Lattice matrix must be 3x3.");
    }
    double* lattice_ptr = mxGetPr(prhs[3]);
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            lattice[i][j] = lattice_ptr[i + 3 * j];
        }
    }

    // 提取和验证 symprec 参数
    if (!mxIsDouble(prhs[4]) || mxGetNumberOfElements(prhs[4]) != 1) {
        mexErrMsgIdAndTxt("Spglib:invalidSymprec", "Symprec must be a scalar.");
    }
    double symprec = mxGetScalar(prhs[4]);

    // 调用 spg_get_spacegroup_type_from_symmetry
    SpglibSpacegroupType spacegroup_type =
        spg_get_spacegroup_type_from_symmetry(rotation, translation,
                                              num_operations, lattice, symprec);

    // 创建 MATLAB 结构体并设置字段
    char const* field_names[] = {"number",
                                 "international_short",
                                 "international_full",
                                 "international",
                                 "schoenflies",
                                 "hall_number",
                                 "hall_symbol",
                                 "choice",
                                 "pointgroup_international",
                                 "pointgroup_schoenflies",
                                 "arithmetic_crystal_class_number",
                                 "arithmetic_crystal_class_symbol"};
    plhs[0] = mxCreateStructMatrix(1, 1, 12, field_names);

    SET_SCALAR_FIELD(plhs[0], 0, "number", spacegroup_type.number);
    SET_STRING_FIELD(plhs[0], 0, "international_short",
                     spacegroup_type.international_short);
    SET_STRING_FIELD(plhs[0], 0, "international_full",
                     spacegroup_type.international_full);
    SET_STRING_FIELD(plhs[0], 0, "international",
                     spacegroup_type.international);
    SET_STRING_FIELD(plhs[0], 0, "schoenflies", spacegroup_type.schoenflies);
    SET_SCALAR_FIELD(plhs[0], 0, "hall_number", spacegroup_type.hall_number);
    SET_STRING_FIELD(plhs[0], 0, "hall_symbol", spacegroup_type.hall_symbol);
    SET_STRING_FIELD(plhs[0], 0, "choice", spacegroup_type.choice);
    SET_STRING_FIELD(plhs[0], 0, "pointgroup_international",
                     spacegroup_type.pointgroup_international);
    SET_STRING_FIELD(plhs[0], 0, "pointgroup_schoenflies",
                     spacegroup_type.pointgroup_schoenflies);
    SET_SCALAR_FIELD(plhs[0], 0, "arithmetic_crystal_class_number",
                     spacegroup_type.arithmetic_crystal_class_number);
    SET_STRING_FIELD(plhs[0], 0, "arithmetic_crystal_class_symbol",
                     spacegroup_type.arithmetic_crystal_class_symbol);
}

// magnetic_spacegroup_type =
// symspg('spg_get_magnetic_spacegroup_type_from_symmetry', rotation,
// translation, time_reversals, num_operations, lattice, symprec)
void SpglibFunctions::spg_get_magnetic_spacegroup_type_from_symmetry_mex(
    int nlhs, mxArray* plhs[], int nrhs, mxArray const* prhs[]) {
    /*
     SpglibMagneticSpacegroupType
     spg_get_magnetic_spacegroup_type_from_symmetry( int const
     rotations[][3][3], double const translations[][3], int const
     *time_reversals, int const num_operations, double const lattice[3][3],
         double const symprec);
    */

    // 验证输入参数数量
    int const expected_number_of_inputs = 6;
    if (nrhs != expected_number_of_inputs) {
        mexErrMsgIdAndTxt("Spglib:invalidNumInputs",
                          "Incorrect number of inputs for "
                          "spg_get_magnetic_spacegroup_type_from_symmetry.");
    }

    // 提取和验证 rotation 参数
    mwSize const* dims = mxGetDimensions(prhs[0]);
    int num_operations = dims[0];  // N 的值
    if (dims[1] != 3 || dims[2] != 3) {
        mexErrMsgIdAndTxt("Spglib:invalidRotation",
                          "Rotation matrix must have dimensions Nx3x3.");
    }

    int rotation[num_operations][3][3];
    auto* rotation_ptr = static_cast<int32_t*>(mxGetData(prhs[0]));

    for (int k = 0; k < num_operations; k++) {
        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++) {
                rotation[k][i][j] = rotation_ptr[k + i * num_operations +
                                                 j * num_operations * 3];
            }
        }
    }

    // 提取和验证 translation 参数
    double translation[num_operations][3];
    if (mxGetM(prhs[1]) != num_operations || mxGetN(prhs[1]) != 3) {
        mexErrMsgIdAndTxt("Spglib:invalidTranslation",
                          "Translation array must be Nx3.");
    }
    double* translation_ptr = mxGetPr(prhs[1]);
    for (int i = 0; i < num_operations; i++) {
        for (int j = 0; j < 3; j++) {
            translation[i][j] = translation_ptr[i + j * num_operations];
        }
    }

    // 提取和验证 time_reversals 参数
    if (mxGetNumberOfElements(prhs[2]) != num_operations) {
        mexErrMsgIdAndTxt(
            "Spglib:invalidTimeReversals",
            "Time reversals array size must match the number of operations.");
    }
    int const* time_reversals = static_cast<int32_t const*>(mxGetData(prhs[2]));

    // 提取和验证 lattice 参数
    double lattice[3][3];
    if (mxGetM(prhs[4]) != 3 || mxGetN(prhs[4]) != 3) {
        mexErrMsgIdAndTxt("Spglib:invalidLattice",
                          "Lattice matrix must be 3x3.");
    }
    double* lattice_ptr = mxGetPr(prhs[4]);
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            lattice[i][j] = lattice_ptr[i + 3 * j];
        }
    }

    // 提取和验证 symprec 参数
    if (!mxIsDouble(prhs[5]) || mxGetNumberOfElements(prhs[5]) != 1) {
        mexErrMsgIdAndTxt("Spglib:invalidSymprec", "Symprec must be a scalar.");
    }
    double symprec = mxGetScalar(prhs[5]);

    // 调用 spg_get_magnetic_spacegroup_type_from_symmetry
    SpglibMagneticSpacegroupType magnetic_spacegroup_type =
        spg_get_magnetic_spacegroup_type_from_symmetry(
            rotation, translation, time_reversals, num_operations, lattice,
            symprec);

    // 创建 MATLAB 结构体并设置字段
    char const* field_names[] = {"uni_number", "litvin_number", "bns_number",
                                 "og_number",  "number",        "type"};
    plhs[0] = mxCreateStructMatrix(1, 1, 6, field_names);

    SET_SCALAR_FIELD(plhs[0], 0, "uni_number",
                     magnetic_spacegroup_type.uni_number);
    SET_SCALAR_FIELD(plhs[0], 0, "litvin_number",
                     magnetic_spacegroup_type.litvin_number);
    SET_STRING_FIELD(plhs[0], 0, "bns_number",
                     magnetic_spacegroup_type.bns_number);
    SET_STRING_FIELD(plhs[0], 0, "og_number",
                     magnetic_spacegroup_type.og_number);
    SET_SCALAR_FIELD(plhs[0], 0, "number", magnetic_spacegroup_type.number);
    SET_SCALAR_FIELD(plhs[0], 0, "type", magnetic_spacegroup_type.type);
}
