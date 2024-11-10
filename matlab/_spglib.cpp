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
    static void spg_get_error_code_mex(int nlhs, mxArray* plhs[], int nrhs,
                                       mxArray const* prhs[]);
    static void spg_get_error_message_mex(int nlhs, mxArray* plhs[], int nrhs,
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
    static void spg_get_pointgroup_mex(int nlhs, mxArray* plhs[], int nrhs,
                                       mxArray const* prhs[]);
    static void spg_get_symmetry_from_database_mex(int nlhs, mxArray* plhs[],
                                                   int nrhs,
                                                   mxArray const* prhs[]);
    static void spg_get_magnetic_symmetry_from_database_mex(
        int nlhs, mxArray* plhs[], int nrhs, mxArray const* prhs[]);
    static void spg_get_spacegroup_type_mex(int nlhs, mxArray* plhs[], int nrhs,
                                            mxArray const* prhs[]);
    static void spg_get_magnetic_spacegroup_type_mex(int nlhs, mxArray* plhs[],
                                                     int nrhs,
                                                     mxArray const* prhs[]);
    static void spg_standardize_cell_mex(int nlhs, mxArray* plhs[], int nrhs,
                                         mxArray const* prhs[]);
    static void spgat_standardize_cell_mex(int nlhs, mxArray* plhs[], int nrhs,
                                           mxArray const* prhs[]);
    static void spg_find_primitive_mex(int nlhs, mxArray* plhs[], int nrhs,
                                       mxArray const* prhs[]);
    static void spgat_find_primitive_mex(int nlhs, mxArray* plhs[], int nrhs,
                                         mxArray const* prhs[]);
    static void spg_refine_cell_mex(int nlhs, mxArray* plhs[], int nrhs,
                                    mxArray const* prhs[]);
    static void spgat_refine_cell_mex(int nlhs, mxArray* plhs[], int nrhs,
                                      mxArray const* prhs[]);
    static void spg_delaunay_reduce_mex(int nlhs, mxArray* plhs[], int nrhs,
                                        mxArray const* prhs[]);
    static void spg_get_grid_point_from_address_mex(int nlhs, mxArray* plhs[],
                                                    int nrhs,
                                                    mxArray const* prhs[]);
    static void spg_get_dense_grid_point_from_address_mex(
        int nlhs, mxArray* plhs[], int nrhs, mxArray const* prhs[]);
    static void spg_get_ir_reciprocal_mesh_mex(int nlhs, mxArray* plhs[],
                                               int nrhs, mxArray const* prhs[]);
    static void spg_get_dense_ir_reciprocal_mesh_mex(int nlhs, mxArray* plhs[],
                                                     int nrhs,
                                                     mxArray const* prhs[]);
    static void spg_get_stabilized_reciprocal_mesh_mex(int nlhs,
                                                       mxArray* plhs[],
                                                       int nrhs,
                                                       mxArray const* prhs[]);
    static void spg_get_dense_stabilized_reciprocal_mesh_mex(
        int nlhs, mxArray* plhs[], int nrhs, mxArray const* prhs[]);
    static void spg_get_dense_grid_points_by_rotations_mex(
        int nlhs, mxArray* plhs[], int nrhs, mxArray const* prhs[]);
    static void spg_get_dense_BZ_grid_points_by_rotations_mex(
        int nlhs, mxArray* plhs[], int nrhs, mxArray const* prhs[]);
    static void spg_relocate_BZ_grid_address_mex(int nlhs, mxArray* plhs[],
                                                 int nrhs,
                                                 mxArray const* prhs[]);
    static void spg_relocate_dense_BZ_grid_address_mex(int nlhs,
                                                       mxArray* plhs[],
                                                       int nrhs,
                                                       mxArray const* prhs[]);
    static void spg_niggli_reduce_mex(int nlhs, mxArray* plhs[], int nrhs,
                                      mxArray const* prhs[]);
};

// 定义类型别名，用于指向静态方法的函数指针
typedef void (*SpglibFunction)(int, mxArray*[], int, mxArray const*[]);

void show_matrix_3d(double const lattice[3][3]) {
    for (int i = 0; i < 3; i++) {
        mexPrintf("%f %f %f\n", lattice[0][i], lattice[1][i], lattice[2][i]);
    }
}

void show_cell(double const lattice[3][3], double const positions[][3],
               int const types[], int const num_atoms) {
    mexPrintf("num_atoms: %d\n", num_atoms);
    mexPrintf("Lattice parameter:\n");
    show_matrix_3d(lattice);
    mexPrintf("Atomic positions:\n");
    for (int i = 0; i < num_atoms; i++) {
        mexPrintf("%d: %f %f %f\n", types[i], positions[i][0], positions[i][1],
                  positions[i][2]);
    }
}

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
        {"spg_get_error_code", SpglibFunctions::spg_get_error_code_mex},
        {"spg_get_error_message", SpglibFunctions::spg_get_error_message_mex},
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
        {"spg_get_pointgroup", SpglibFunctions::spg_get_pointgroup_mex},
        {"spg_get_symmetry_from_database",
         SpglibFunctions::spg_get_symmetry_from_database_mex},
        {"spg_get_magnetic_symmetry_from_database",
         SpglibFunctions::spg_get_magnetic_symmetry_from_database_mex},
        {"spg_get_spacegroup_type",
         SpglibFunctions::spg_get_spacegroup_type_mex},
        {"spg_get_magnetic_spacegroup_type",
         SpglibFunctions::spg_get_magnetic_spacegroup_type_mex},
        {"spg_standardize_cell", SpglibFunctions::spg_standardize_cell_mex},
        {"spgat_standardize_cell", SpglibFunctions::spgat_standardize_cell_mex},
        {"spg_find_primitive", SpglibFunctions::spg_find_primitive_mex},
        {"spgat_find_primitive", SpglibFunctions::spgat_find_primitive_mex},
        {"spg_refine_cell", SpglibFunctions::spg_refine_cell_mex},
        {"spgat_refine_cell", SpglibFunctions::spgat_refine_cell_mex},
        {"spg_delaunay_reduce", SpglibFunctions::spg_delaunay_reduce_mex},
        {"spg_get_grid_point_from_address",
         SpglibFunctions::spg_get_grid_point_from_address_mex},
        {"spg_get_dense_grid_point_from_address",
         SpglibFunctions::spg_get_dense_grid_point_from_address_mex},
        {"spg_get_ir_reciprocal_mesh",
         SpglibFunctions::spg_get_ir_reciprocal_mesh_mex},
        {"spg_get_dense_ir_reciprocal_mesh",
         SpglibFunctions::spg_get_dense_ir_reciprocal_mesh_mex},
        {"spg_get_stabilized_reciprocal_mesh",
         SpglibFunctions::spg_get_stabilized_reciprocal_mesh_mex},
        {"spg_get_dense_stabilized_reciprocal_mesh",
         SpglibFunctions::spg_get_dense_stabilized_reciprocal_mesh_mex},
        {"spg_get_dense_grid_points_by_rotations",
         SpglibFunctions::spg_get_dense_grid_points_by_rotations_mex},
        {"spg_get_dense_BZ_grid_points_by_rotations",
         SpglibFunctions::spg_get_dense_BZ_grid_points_by_rotations_mex},
        {"spg_relocate_BZ_grid_address",
         SpglibFunctions::spg_relocate_BZ_grid_address_mex},
        {"spg_relocate_dense_BZ_grid_address",
         SpglibFunctions::spg_relocate_dense_BZ_grid_address_mex},
        {"spg_niggli_reduce", SpglibFunctions::spg_niggli_reduce_mex}};

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

// error_code = symspg('spg_get_error_code')
void SpglibFunctions::spg_get_error_code_mex(int nlhs, mxArray* plhs[],
                                             int nrhs, mxArray const* prhs[]) {
    /*
     SpglibError spg_get_error_code(void);
    */

    // 验证输入参数数量
    if (nrhs != 0) {
        mexErrMsgIdAndTxt("Spglib:invalidNumInputs",
                          "No inputs expected for spg_get_error_code.");
    }

    // 调用 spg_get_error_code
    SpglibError error_code = spg_get_error_code();

    // 创建输出标量 error_code
    plhs[0] = mxCreateDoubleScalar(static_cast<double>(error_code));
}

// error_message = symspg('spg_get_error_message', error_code)
void SpglibFunctions::spg_get_error_message_mex(int nlhs, mxArray* plhs[],
                                                int nrhs,
                                                mxArray const* prhs[]) {
    /*
     char *spg_get_error_message(SpglibError spglib_error);
    */

    // 验证输入参数数量
    int const expected_number_of_inputs = 1;
    if (nrhs != expected_number_of_inputs) {
        mexErrMsgIdAndTxt(
            "Spglib:invalidNumInputs",
            "Incorrect number of inputs for spg_get_error_message.");
    }

    // 提取和验证 spglib_error 参数
    SpglibError spglib_error = static_cast<SpglibError>(mxGetScalar(prhs[0]));

    // 调用 spg_get_error_message
    char* error_message = spg_get_error_message(spglib_error);

    // 创建输出字符串
    plhs[0] = mxCreateString(error_message);
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

// [symbol, trans_mat, result] =
// symspg('spg_get_pointgroup', rotations, num_rotations)
void SpglibFunctions::spg_get_pointgroup_mex(int nlhs, mxArray* plhs[],
                                             int nrhs, mxArray const* prhs[]) {
    /*
     int spg_get_pointgroup(char symbol[6], int trans_mat[3][3],
                            int const rotations[][3][3], int const
     num_rotations);
    */

    // 验证输入参数数量
    int const expected_number_of_inputs = 2;
    if (nrhs != expected_number_of_inputs) {
        mexErrMsgIdAndTxt("Spglib:invalidNumInputs",
                          "Incorrect number of inputs for spg_get_pointgroup.");
    }

    // 提取和验证 rotation 参数
    mwSize const* dims = mxGetDimensions(prhs[0]);
    int num_operations = dims[0];
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

    char symbol[6];
    int trans_mat[3][3];

    // 调用 spg_get_pointgroup
    int result =
        spg_get_pointgroup(symbol, trans_mat, rotation, num_operations);

    // 创建输出数组并设置字段
    // 设置 symbol 字段
    plhs[0] = mxCreateString(symbol);

    // 输出变换矩阵 (3x3 int 数组)
    plhs[1] = mxCreateNumericMatrix(3, 3, mxINT32_CLASS, mxREAL);
    int* trans_mat_out = static_cast<int*>(mxGetData(plhs[1]));
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            trans_mat_out[i + j * 3] = trans_mat[i][j];
        }
    }

    // 设置 result 字段
    plhs[2] = mxCreateDoubleScalar(static_cast<double>(result));
}

// [rotations, translations] = symspg('spg_get_symmetry_from_database',
// hall_number)
void SpglibFunctions::spg_get_symmetry_from_database_mex(
    int nlhs, mxArray* plhs[], int nrhs, mxArray const* prhs[]) {
    /*
     int spg_get_symmetry_from_database(int rotations[192][3][3],
                                        double translations[192][3],
                                        int const hall_number);
    */

    // 验证输入参数数量
    int const expected_number_of_inputs = 1;
    if (nrhs != expected_number_of_inputs) {
        mexErrMsgIdAndTxt(
            "Spglib:invalidNumInputs",
            "Incorrect number of inputs for spg_get_symmetry_from_database.");
    }

    // 提取和验证 hall_number 参数
    if (mxGetNumberOfElements(prhs[0]) != 1) {
        mexErrMsgIdAndTxt("Spglib:invalidHallNumber",
                          "Hall number must be a scalar.");
    }
    int hall_number = static_cast<int>(mxGetScalar(prhs[0]));

    // 准备输出数组
    int rotations[192][3][3];
    double translations[192][3];

    // 调用 spg_get_symmetry_from_database
    int num_operations =
        spg_get_symmetry_from_database(rotations, translations, hall_number);

    // 创建输出 rotation 数组 (Nx3x3 int 数组)
    mwSize rotation_dims[3] = {static_cast<mwSize>(num_operations), 3, 3};
    plhs[0] = mxCreateNumericArray(3, rotation_dims, mxINT32_CLASS, mxREAL);
    int* rotations_out = static_cast<int*>(mxGetData(plhs[0]));
    for (int k = 0; k < num_operations; ++k) {
        for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < 3; ++j) {
                rotations_out[k + i * num_operations + j * num_operations * 3] =
                    rotations[k][i][j];
            }
        }
    }

    // 创建输出 translation 数组 (Nx3 double 数组)
    plhs[1] = mxCreateDoubleMatrix(num_operations, 3, mxREAL);
    double* translations_out = mxGetPr(plhs[1]);
    for (int i = 0; i < num_operations; ++i) {
        for (int j = 0; j < 3; ++j) {
            translations_out[i + j * num_operations] = translations[i][j];
        }
    }
}

// [rotations, translations, time_reversals] =
// symspg('spg_get_magnetic_symmetry_from_database', uni_number, hall_number)
void SpglibFunctions::spg_get_magnetic_symmetry_from_database_mex(
    int nlhs, mxArray* plhs[], int nrhs, mxArray const* prhs[]) {
    /*
     int spg_get_magnetic_symmetry_from_database(int rotations[384][3][3],
                                                 double translations[384][3],
                                                 int time_reversals[384],
                                                 int const uni_number,
                                                 int const hall_number);
    */

    // 验证输入参数数量
    int const expected_number_of_inputs = 2;
    if (nrhs != expected_number_of_inputs) {
        mexErrMsgIdAndTxt("Spglib:invalidNumInputs",
                          "Incorrect number of inputs for "
                          "spg_get_magnetic_symmetry_from_database.");
    }

    // 提取和验证 uni_number 参数
    if (mxGetNumberOfElements(prhs[0]) != 1) {
        mexErrMsgIdAndTxt("Spglib:invalidUniNumber",
                          "Uni number must be a scalar.");
    }
    int uni_number = static_cast<int>(mxGetScalar(prhs[0]));

    // 提取和验证 hall_number 参数
    if (mxGetNumberOfElements(prhs[1]) != 1) {
        mexErrMsgIdAndTxt("Spglib:invalidHallNumber",
                          "Hall number must be a scalar.");
    }
    int hall_number = static_cast<int>(mxGetScalar(prhs[1]));

    // 准备输出数组
    int rotations[384][3][3];
    double translations[384][3];
    int time_reversals[384];

    // 调用 spg_get_magnetic_symmetry_from_database
    int num_operations = spg_get_magnetic_symmetry_from_database(
        rotations, translations, time_reversals, uni_number, hall_number);

    // 创建输出 rotation 数组 (Nx3x3 int 数组)
    mwSize rotation_dims[3] = {static_cast<mwSize>(num_operations), 3, 3};
    plhs[0] = mxCreateNumericArray(3, rotation_dims, mxINT32_CLASS, mxREAL);
    int* rotations_out = static_cast<int*>(mxGetData(plhs[0]));
    for (int k = 0; k < num_operations; ++k) {
        for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < 3; ++j) {
                rotations_out[k + i * num_operations + j * num_operations * 3] =
                    rotations[k][i][j];
            }
        }
    }

    // 创建输出 translation 数组 (Nx3 double 数组)
    plhs[1] = mxCreateDoubleMatrix(num_operations, 3, mxREAL);
    double* translations_out = mxGetPr(plhs[1]);
    for (int i = 0; i < num_operations; ++i) {
        for (int j = 0; j < 3; ++j) {
            translations_out[i + j * num_operations] = translations[i][j];
        }
    }

    // 创建输出 time_reversals 数组 (N int 数组)
    plhs[2] = mxCreateNumericMatrix(num_operations, 1, mxINT32_CLASS, mxREAL);
    int* time_reversals_out = static_cast<int*>(mxGetData(plhs[2]));
    for (int i = 0; i < num_operations; ++i) {
        time_reversals_out[i] = time_reversals[i];
    }
}

// spacegroup = symspg('spg_get_spacegroup_type', hall_number)
void SpglibFunctions::spg_get_spacegroup_type_mex(int nlhs, mxArray* plhs[],
                                                  int nrhs,
                                                  mxArray const* prhs[]) {
    /*
     SpglibSpacegroupType spg_get_spacegroup_type(int const hall_number);
    */

    // 验证输入参数数量
    int const expected_number_of_inputs = 1;
    if (nrhs != expected_number_of_inputs) {
        mexErrMsgIdAndTxt(
            "Spglib:invalidNumInputs",
            "Incorrect number of inputs for spg_get_spacegroup_type.");
    }

    // 提取和验证 hall_number 参数
    if (mxGetNumberOfElements(prhs[0]) != 1) {
        mexErrMsgIdAndTxt("Spglib:invalidHallNumber",
                          "Hall number must be a scalar.");
    }
    int hall_number = static_cast<int>(mxGetScalar(prhs[0]));

    // 调用 spg_get_spacegroup_type
    SpglibSpacegroupType spacegroup = spg_get_spacegroup_type(hall_number);

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

    // 设置结构体字段
    mxSetField(plhs[0], 0, "number",
               mxCreateDoubleScalar(static_cast<double>(spacegroup.number)));
    mxSetField(plhs[0], 0, "international_short",
               mxCreateString(spacegroup.international_short));
    mxSetField(plhs[0], 0, "international_full",
               mxCreateString(spacegroup.international_full));
    mxSetField(plhs[0], 0, "international",
               mxCreateString(spacegroup.international));
    mxSetField(plhs[0], 0, "schoenflies",
               mxCreateString(spacegroup.schoenflies));
    mxSetField(
        plhs[0], 0, "hall_number",
        mxCreateDoubleScalar(static_cast<double>(spacegroup.hall_number)));
    mxSetField(plhs[0], 0, "hall_symbol",
               mxCreateString(spacegroup.hall_symbol));
    mxSetField(plhs[0], 0, "choice", mxCreateString(spacegroup.choice));
    mxSetField(plhs[0], 0, "pointgroup_international",
               mxCreateString(spacegroup.pointgroup_international));
    mxSetField(plhs[0], 0, "pointgroup_schoenflies",
               mxCreateString(spacegroup.pointgroup_schoenflies));
    mxSetField(plhs[0], 0, "arithmetic_crystal_class_number",
               mxCreateDoubleScalar(static_cast<double>(
                   spacegroup.arithmetic_crystal_class_number)));
    mxSetField(plhs[0], 0, "arithmetic_crystal_class_symbol",
               mxCreateString(spacegroup.arithmetic_crystal_class_symbol));
}

// spacegroup = symspg('spg_get_magnetic_spacegroup_type', uni_number)
void SpglibFunctions::spg_get_magnetic_spacegroup_type_mex(
    int nlhs, mxArray* plhs[], int nrhs, mxArray const* prhs[]) {
    /*
     SpglibMagneticSpacegroupType spg_get_magnetic_spacegroup_type(int const
     uni_number);
    */

    // 验证输入参数数量
    int const expected_number_of_inputs = 1;
    if (nrhs != expected_number_of_inputs) {
        mexErrMsgIdAndTxt(
            "Spglib:invalidNumInputs",
            "Incorrect number of inputs for spg_get_magnetic_spacegroup_type.");
    }

    // 提取和验证 uni_number 参数
    if (mxGetNumberOfElements(prhs[0]) != 1) {
        mexErrMsgIdAndTxt("Spglib:invalidUniNumber",
                          "Uni number must be a scalar.");
    }
    int uni_number = static_cast<int>(mxGetScalar(prhs[0]));

    // 调用 spg_get_magnetic_spacegroup_type
    SpglibMagneticSpacegroupType spacegroup =
        spg_get_magnetic_spacegroup_type(uni_number);

    // 创建 MATLAB 结构体并设置字段
    char const* field_names[] = {"uni_number", "litvin_number", "bns_number",
                                 "og_number",  "number",        "type"};
    plhs[0] = mxCreateStructMatrix(1, 1, 6, field_names);

    // 设置结构体字段
    mxSetField(
        plhs[0], 0, "uni_number",
        mxCreateDoubleScalar(static_cast<double>(spacegroup.uni_number)));
    mxSetField(
        plhs[0], 0, "litvin_number",
        mxCreateDoubleScalar(static_cast<double>(spacegroup.litvin_number)));
    mxSetField(plhs[0], 0, "bns_number", mxCreateString(spacegroup.bns_number));
    mxSetField(plhs[0], 0, "og_number", mxCreateString(spacegroup.og_number));
    mxSetField(plhs[0], 0, "number",
               mxCreateDoubleScalar(static_cast<double>(spacegroup.number)));
    mxSetField(plhs[0], 0, "type",
               mxCreateDoubleScalar(static_cast<double>(spacegroup.type)));
}

// [lattice, position, types, num_primitive_atom] =
// symspg('spg_standardize_cell', lattice, position, types, num_atom,
// to_primitive, no_idealize, symprec)
void SpglibFunctions::spg_standardize_cell_mex(int nlhs, mxArray* plhs[],
                                               int nrhs,
                                               mxArray const* prhs[]) {
    /*
     int spg_standardize_cell(double lattice[3][3], double position[][3],
                              int types[], int const num_atom,
                              int const to_primitive, int const no_idealize,
                              double const symprec);
    */

    // 验证输入参数数量
    int const expected_number_of_inputs = 7;
    if (nrhs != expected_number_of_inputs) {
        mexErrMsgIdAndTxt(
            "Spglib:invalidNumInputs",
            "Incorrect number of inputs for spg_standardize_cell.");
    }

    // 提取和验证 lattice 参数
    if (mxGetM(prhs[0]) != 3 || mxGetN(prhs[0]) != 3) {
        mexErrMsgIdAndTxt("Spglib:invalidLattice",
                          "Lattice must be a 3x3 matrix.");
    }
    double* lattice_ptr = mxGetPr(prhs[0]);
    double lattice[3][3];
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            lattice[i][j] = lattice_ptr[i + 3 * j];
        }
    }

    // 提取和验证 position 参数
    mwSize num_atom = mxGetM(prhs[1]);
    if (mxGetN(prhs[1]) != 3) {
        mexErrMsgIdAndTxt("Spglib:invalidPosition",
                          "Position must be an Nx3 matrix.");
    }
    double* position_ptr = mxGetPr(prhs[1]);
    double position[num_atom][3];
    for (mwSize i = 0; i < num_atom; ++i) {
        for (int j = 0; j < 3; ++j) {
            position[i][j] = position_ptr[i + j * num_atom];
        }
    }

    // 提取和验证 types 参数
    if (mxGetNumberOfElements(prhs[2]) != num_atom) {
        mexErrMsgIdAndTxt("Spglib:invalidTypes",
                          "Types array size must match the number of atoms.");
    }
    int* types_ptr = static_cast<int*>(mxGetData(prhs[2]));
    int types[num_atom];
    for (mwSize i = 0; i < num_atom; ++i) {
        types[i] = types_ptr[i];
    }

    // 提取和验证其他标量参数
    int to_primitive = static_cast<int>(mxGetScalar(prhs[3]));
    int no_idealize = static_cast<int>(mxGetScalar(prhs[4]));
    double symprec = mxGetScalar(prhs[5]);

    // 调用 spg_standardize_cell
    int num_primitive_atom = spg_standardize_cell(
        lattice, position, types, num_atom, to_primitive, no_idealize, symprec);

    // 创建输出 lattice 数组 (3x3 double 数组)
    plhs[0] = mxCreateDoubleMatrix(3, 3, mxREAL);
    double* lattice_out = mxGetPr(plhs[0]);
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            lattice_out[i + 3 * j] = lattice[i][j];
        }
    }

    // 创建输出 position 数组 (num_primitive_atom x 3 double 数组)
    plhs[1] = mxCreateDoubleMatrix(num_primitive_atom, 3, mxREAL);
    double* position_out = mxGetPr(plhs[1]);
    for (int i = 0; i < num_primitive_atom; ++i) {
        for (int j = 0; j < 3; ++j) {
            position_out[i + j * num_primitive_atom] = position[i][j];
        }
    }

    // 创建输出 types 数组 (num_primitive_atom int 数组)
    plhs[2] =
        mxCreateNumericMatrix(num_primitive_atom, 1, mxINT32_CLASS, mxREAL);
    int* types_out = static_cast<int*>(mxGetData(plhs[2]));
    for (int i = 0; i < num_primitive_atom; ++i) {
        types_out[i] = types[i];
    }

    // 创建输出 num_primitive_atom 标量
    plhs[3] = mxCreateDoubleScalar(static_cast<double>(num_primitive_atom));
}

// [lattice, position, types, num_primitive_atom] =
// symspg('spgat_standardize_cell', lattice, position, types, num_atom,
// to_primitive, no_idealize, symprec, angle_tolerance)
void SpglibFunctions::spgat_standardize_cell_mex(int nlhs, mxArray* plhs[],
                                                 int nrhs,
                                                 mxArray const* prhs[]) {
    /*
     int spgat_standardize_cell(double lattice[3][3], double position[][3],
                                int types[], int const num_atom,
                                int const to_primitive, int const no_idealize,
                                double const symprec, double const
     angle_tolerance);
    */

    // 验证输入参数数量
    int const expected_number_of_inputs = 8;
    if (nrhs != expected_number_of_inputs) {
        mexErrMsgIdAndTxt(
            "Spglib:invalidNumInputs",
            "Incorrect number of inputs for spgat_standardize_cell.");
    }

    // 提取和验证 lattice 参数
    if (mxGetM(prhs[0]) != 3 || mxGetN(prhs[0]) != 3) {
        mexErrMsgIdAndTxt("Spglib:invalidLattice",
                          "Lattice must be a 3x3 matrix.");
    }
    double* lattice_ptr = mxGetPr(prhs[0]);
    double lattice[3][3];
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            lattice[i][j] = lattice_ptr[i + 3 * j];
        }
    }

    // 提取和验证 position 参数
    mwSize num_atom = mxGetM(prhs[1]);
    if (mxGetN(prhs[1]) != 3) {
        mexErrMsgIdAndTxt("Spglib:invalidPosition",
                          "Position must be an Nx3 matrix.");
    }
    double* position_ptr = mxGetPr(prhs[1]);
    double position[num_atom][3];
    for (mwSize i = 0; i < num_atom; ++i) {
        for (int j = 0; j < 3; ++j) {
            position[i][j] = position_ptr[i + j * num_atom];
        }
    }

    // 提取和验证 types 参数
    if (mxGetNumberOfElements(prhs[2]) != num_atom) {
        mexErrMsgIdAndTxt("Spglib:invalidTypes",
                          "Types array size must match the number of atoms.");
    }
    int* types_ptr = static_cast<int*>(mxGetData(prhs[2]));
    int types[num_atom];
    for (mwSize i = 0; i < num_atom; ++i) {
        types[i] = types_ptr[i];
    }

    // 提取和验证其他标量参数
    int to_primitive = static_cast<int>(mxGetScalar(prhs[3]));
    int no_idealize = static_cast<int>(mxGetScalar(prhs[4]));
    double symprec = mxGetScalar(prhs[5]);
    double angle_tolerance = mxGetScalar(prhs[6]);

    // 调用 spgat_standardize_cell
    int num_primitive_atom =
        spgat_standardize_cell(lattice, position, types, num_atom, to_primitive,
                               no_idealize, symprec, angle_tolerance);

    // 创建输出 lattice 数组 (3x3 double 数组)
    plhs[0] = mxCreateDoubleMatrix(3, 3, mxREAL);
    double* lattice_out = mxGetPr(plhs[0]);
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            lattice_out[i + 3 * j] = lattice[i][j];
        }
    }

    // 创建输出 position 数组 (num_primitive_atom x 3 double 数组)
    plhs[1] = mxCreateDoubleMatrix(num_primitive_atom, 3, mxREAL);
    double* position_out = mxGetPr(plhs[1]);
    for (int i = 0; i < num_primitive_atom; ++i) {
        for (int j = 0; j < 3; ++j) {
            position_out[i + j * num_primitive_atom] = position[i][j];
        }
    }

    // 创建输出 types 数组 (num_primitive_atom int 数组)
    plhs[2] =
        mxCreateNumericMatrix(num_primitive_atom, 1, mxINT32_CLASS, mxREAL);
    int* types_out = static_cast<int*>(mxGetData(plhs[2]));
    for (int i = 0; i < num_primitive_atom; ++i) {
        types_out[i] = types[i];
    }

    // 创建输出 num_primitive_atom 标量
    plhs[3] = mxCreateDoubleScalar(static_cast<double>(num_primitive_atom));
}

// [lattice, position, types, num_primitive_atom] = symspg('spg_find_primitive',
// lattice, position, types, num_atom, symprec)
void SpglibFunctions::spg_find_primitive_mex(int nlhs, mxArray* plhs[],
                                             int nrhs, mxArray const* prhs[]) {
    /*
     int spg_find_primitive(double lattice[3][3], double position[][3],
                            int types[], int const num_atom,
                            double const symprec);
    */

    // 验证输入参数数量
    int const expected_number_of_inputs = 5;
    if (nrhs != expected_number_of_inputs) {
        mexErrMsgIdAndTxt("Spglib:invalidNumInputs",
                          "Incorrect number of inputs for spg_find_primitive.");
    }

    // 提取和验证 lattice 参数
    if (mxGetM(prhs[0]) != 3 || mxGetN(prhs[0]) != 3) {
        mexErrMsgIdAndTxt("Spglib:invalidLattice",
                          "Lattice must be a 3x3 matrix.");
    }
    double* lattice_ptr = mxGetPr(prhs[0]);
    double lattice[3][3];
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            lattice[i][j] = lattice_ptr[i + 3 * j];
        }
    }

    // 提取和验证 position 参数
    mwSize num_atom = mxGetM(prhs[1]);
    if (mxGetN(prhs[1]) != 3) {
        mexErrMsgIdAndTxt("Spglib:invalidPosition",
                          "Position must be an Nx3 matrix.");
    }
    double* position_ptr = mxGetPr(prhs[1]);
    double position[num_atom][3];
    for (mwSize i = 0; i < num_atom; ++i) {
        for (int j = 0; j < 3; ++j) {
            position[i][j] = position_ptr[i + j * num_atom];
        }
    }

    // 提取和验证 types 参数
    if (mxGetNumberOfElements(prhs[2]) != num_atom) {
        mexErrMsgIdAndTxt("Spglib:invalidTypes",
                          "Types array size must match the number of atoms.");
    }
    int* types_ptr = static_cast<int*>(mxGetData(prhs[2]));
    int types[num_atom];
    for (mwSize i = 0; i < num_atom; ++i) {
        types[i] = types_ptr[i];
    }

    // 提取和验证 symprec 参数
    double symprec = mxGetScalar(prhs[4]);

    // 调用 spg_find_primitive
    int num_primitive_atom =
        spg_find_primitive(lattice, position, types, num_atom, symprec);

    // 创建输出 lattice 数组 (3x3 double 数组)
    plhs[0] = mxCreateDoubleMatrix(3, 3, mxREAL);
    double* lattice_out = mxGetPr(plhs[0]);
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            lattice_out[i + 3 * j] = lattice[i][j];
        }
    }

    // 创建输出 position 数组 (num_primitive_atom x 3 double 数组)
    plhs[1] = mxCreateDoubleMatrix(num_primitive_atom, 3, mxREAL);
    double* position_out = mxGetPr(plhs[1]);
    for (int i = 0; i < num_primitive_atom; ++i) {
        for (int j = 0; j < 3; ++j) {
            position_out[i + j * num_primitive_atom] = position[i][j];
        }
    }

    // 创建输出 types 数组 (num_primitive_atom int 数组)
    plhs[2] =
        mxCreateNumericMatrix(num_primitive_atom, 1, mxINT32_CLASS, mxREAL);
    int* types_out = static_cast<int*>(mxGetData(plhs[2]));
    for (int i = 0; i < num_primitive_atom; ++i) {
        types_out[i] = types[i];
    }

    // 创建输出 num_primitive_atom 标量
    plhs[3] = mxCreateDoubleScalar(static_cast<double>(num_primitive_atom));
}

// [lattice, position, types, num_primitive_atom] =
// symspg('spgat_find_primitive', lattice, position, types, num_atom, symprec,
// angle_tolerance)
void SpglibFunctions::spgat_find_primitive_mex(int nlhs, mxArray* plhs[],
                                               int nrhs,
                                               mxArray const* prhs[]) {
    /*
     int spgat_find_primitive(double lattice[3][3], double position[][3],
                              int types[], int const num_atom,
                              double const symprec, double const
     angle_tolerance);
    */

    // 验证输入参数数量
    int const expected_number_of_inputs = 6;
    if (nrhs != expected_number_of_inputs) {
        mexErrMsgIdAndTxt(
            "Spglib:invalidNumInputs",
            "Incorrect number of inputs for spgat_find_primitive.");
    }

    // 提取和验证 lattice 参数
    if (mxGetM(prhs[0]) != 3 || mxGetN(prhs[0]) != 3) {
        mexErrMsgIdAndTxt("Spglib:invalidLattice",
                          "Lattice must be a 3x3 matrix.");
    }
    double* lattice_ptr = mxGetPr(prhs[0]);
    double lattice[3][3];
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            lattice[i][j] = lattice_ptr[i + 3 * j];
        }
    }

    // 提取和验证 position 参数
    mwSize num_atom = mxGetM(prhs[1]);
    if (mxGetN(prhs[1]) != 3) {
        mexErrMsgIdAndTxt("Spglib:invalidPosition",
                          "Position must be an Nx3 matrix.");
    }
    double* position_ptr = mxGetPr(prhs[1]);
    double position[num_atom][3];
    for (mwSize i = 0; i < num_atom; ++i) {
        for (int j = 0; j < 3; ++j) {
            position[i][j] = position_ptr[i + j * num_atom];
        }
    }

    // 提取和验证 types 参数
    if (mxGetNumberOfElements(prhs[2]) != num_atom) {
        mexErrMsgIdAndTxt("Spglib:invalidTypes",
                          "Types array size must match the number of atoms.");
    }
    int* types_ptr = static_cast<int*>(mxGetData(prhs[2]));
    int types[num_atom];
    for (mwSize i = 0; i < num_atom; ++i) {
        types[i] = types_ptr[i];
    }

    // 提取和验证 symprec 参数
    double symprec = mxGetScalar(prhs[4]);

    // 提取和验证 angle_tolerance 参数
    double angle_tolerance = mxGetScalar(prhs[5]);

    // 调用 spgat_find_primitive
    int num_primitive_atom = spgat_find_primitive(
        lattice, position, types, num_atom, symprec, angle_tolerance);

    // 创建输出 lattice 数组 (3x3 double 数组)
    plhs[0] = mxCreateDoubleMatrix(3, 3, mxREAL);
    double* lattice_out = mxGetPr(plhs[0]);
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            lattice_out[i + 3 * j] = lattice[i][j];
        }
    }

    // 创建输出 position 数组 (num_primitive_atom x 3 double 数组)
    plhs[1] = mxCreateDoubleMatrix(num_primitive_atom, 3, mxREAL);
    double* position_out = mxGetPr(plhs[1]);
    for (int i = 0; i < num_primitive_atom; ++i) {
        for (int j = 0; j < 3; ++j) {
            position_out[i + j * num_primitive_atom] = position[i][j];
        }
    }

    // 创建输出 types 数组 (num_primitive_atom int 数组)
    plhs[2] =
        mxCreateNumericMatrix(num_primitive_atom, 1, mxINT32_CLASS, mxREAL);
    int* types_out = static_cast<int*>(mxGetData(plhs[2]));
    for (int i = 0; i < num_primitive_atom; ++i) {
        types_out[i] = types[i];
    }

    // 创建输出 num_primitive_atom 标量
    plhs[3] = mxCreateDoubleScalar(static_cast<double>(num_primitive_atom));
}

// [lattice, position, types, num_atom_bravais] = symspg('spg_refine_cell',
// lattice, position, types, num_atom, symprec)
void SpglibFunctions::spg_refine_cell_mex(int nlhs, mxArray* plhs[], int nrhs,
                                          mxArray const* prhs[]) {
    /*
     int spg_refine_cell(double lattice[3][3], double position[][3],
                         int types[], int const num_atom,
                         double const symprec);
    */

    // 验证输入参数数量
    int const expected_number_of_inputs = 5;
    if (nrhs != expected_number_of_inputs) {
        mexErrMsgIdAndTxt("Spglib:invalidNumInputs",
                          "Incorrect number of inputs for spg_refine_cell.");
    }

    // 提取和验证 lattice 参数
    if (mxGetM(prhs[0]) != 3 || mxGetN(prhs[0]) != 3) {
        mexErrMsgIdAndTxt("Spglib:invalidLattice",
                          "Lattice must be a 3x3 matrix.");
    }
    double* lattice_ptr = mxGetPr(prhs[0]);
    double lattice[3][3];
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            lattice[i][j] = lattice_ptr[i + 3 * j];
        }
    }

    // 提取和验证 position 参数
    mwSize num_atom = mxGetM(prhs[1]);
    if (mxGetN(prhs[1]) != 3) {
        mexErrMsgIdAndTxt("Spglib:invalidPosition",
                          "Position must be an Nx3 matrix.");
    }
    double* position_ptr = mxGetPr(prhs[1]);
    double position[4 * num_atom][3];  // 这里必须设定 position
                                       // 数组的第一个维度为 4 * num_atom
    for (mwSize i = 0; i < num_atom; ++i) {
        for (int j = 0; j < 3; ++j) {
            position[i][j] = position_ptr[i + j * num_atom];
        }
    }

    // 提取和验证 types 参数
    if (mxGetNumberOfElements(prhs[2]) != num_atom) {
        mexErrMsgIdAndTxt("Spglib:invalidTypes",
                          "Types array size must match the number of atoms.");
    }
    int* types_ptr = static_cast<int*>(mxGetData(prhs[2]));
    int types[4 *
              num_atom];  // 这里必须设定 types 数组的第一个维度为 4 * num_atom
    for (mwSize i = 0; i < num_atom; ++i) {
        types[i] = types_ptr[i];
    }

    // 提取和验证 symprec 参数
    double symprec = mxGetScalar(prhs[4]);

    // 调用 spg_refine_cell
    int num_atom_bravais =
        spg_refine_cell(lattice, position, types, num_atom, symprec);

    // 创建输出 lattice 数组 (3x3 double 数组)
    plhs[0] = mxCreateDoubleMatrix(3, 3, mxREAL);
    double* lattice_out = mxGetPr(plhs[0]);
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            lattice_out[i + 3 * j] = lattice[i][j];
        }
    }

    // 创建输出 position 数组 (num_atom_bravais x 3 double 数组)
    plhs[1] = mxCreateDoubleMatrix(num_atom_bravais, 3, mxREAL);
    double* position_out = mxGetPr(plhs[1]);
    for (int i = 0; i < num_atom_bravais; ++i) {
        for (int j = 0; j < 3; ++j) {
            position_out[i + j * num_atom_bravais] = position[i][j];
        }
    }

    // 创建输出 types 数组 (num_atom_bravais int 数组)
    plhs[2] = mxCreateNumericMatrix(num_atom_bravais, 1, mxINT32_CLASS, mxREAL);
    int* types_out = static_cast<int*>(mxGetData(plhs[2]));
    for (int i = 0; i < num_atom_bravais; ++i) {
        types_out[i] = types[i];
    }

    // 创建输出 num_atom_bravais 标量
    plhs[3] = mxCreateDoubleScalar(static_cast<double>(num_atom_bravais));
}

// [lattice, position, types, num_atom_bravais] = symspg('spgat_refine_cell',
// lattice, position, types, num_atom, symprec, angle_tolerance)
void SpglibFunctions::spgat_refine_cell_mex(int nlhs, mxArray* plhs[], int nrhs,
                                            mxArray const* prhs[]) {
    /*
     int spgat_refine_cell(double lattice[3][3], double position[][3],
                           int types[], int const num_atom,
                           double const symprec, double const angle_tolerance);
    */

    // 验证输入参数数量
    int const expected_number_of_inputs = 6;
    if (nrhs != expected_number_of_inputs) {
        mexErrMsgIdAndTxt("Spglib:invalidNumInputs",
                          "Incorrect number of inputs for spgat_refine_cell.");
    }

    // 提取和验证 lattice 参数
    if (mxGetM(prhs[0]) != 3 || mxGetN(prhs[0]) != 3) {
        mexErrMsgIdAndTxt("Spglib:invalidLattice",
                          "Lattice must be a 3x3 matrix.");
    }
    double* lattice_ptr = mxGetPr(prhs[0]);
    double lattice[3][3];
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            lattice[i][j] = lattice_ptr[i + 3 * j];
        }
    }

    // 提取和验证 position 参数
    mwSize num_atom = mxGetM(prhs[1]);
    if (mxGetN(prhs[1]) != 3) {
        mexErrMsgIdAndTxt("Spglib:invalidPosition",
                          "Position must be an Nx3 matrix.");
    }
    double* position_ptr = mxGetPr(prhs[1]);
    double position[4 * num_atom][3];  // 这里必须设定 position
                                       // 数组的第一个维度为 4 * num_atom
    for (mwSize i = 0; i < num_atom; ++i) {
        for (int j = 0; j < 3; ++j) {
            position[i][j] = position_ptr[i + j * num_atom];
        }
    }

    // 提取和验证 types 参数
    if (mxGetNumberOfElements(prhs[2]) != num_atom) {
        mexErrMsgIdAndTxt("Spglib:invalidTypes",
                          "Types array size must match the number of atoms.");
    }
    int* types_ptr = static_cast<int*>(mxGetData(prhs[2]));
    int types[4 *
              num_atom];  // 这里必须设定 types 数组的第一个维度为 4 * num_atom
    for (mwSize i = 0; i < num_atom; ++i) {
        types[i] = types_ptr[i];
    }

    // 提取和验证 symprec 参数
    double symprec = mxGetScalar(prhs[4]);

    // 提取和验证 angle_tolerance 参数
    double angle_tolerance = mxGetScalar(prhs[5]);

    // 调用 spgat_refine_cell
    int num_atom_bravais = spgat_refine_cell(lattice, position, types, num_atom,
                                             symprec, angle_tolerance);

    // 创建输出 lattice 数组 (3x3 double 数组)
    plhs[0] = mxCreateDoubleMatrix(3, 3, mxREAL);
    double* lattice_out = mxGetPr(plhs[0]);
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            lattice_out[i + 3 * j] = lattice[i][j];
        }
    }

    // 创建输出 position 数组 (num_atom_bravais x 3 double 数组)
    plhs[1] = mxCreateDoubleMatrix(num_atom_bravais, 3, mxREAL);
    double* position_out = mxGetPr(plhs[1]);
    for (int i = 0; i < num_atom_bravais; ++i) {
        for (int j = 0; j < 3; ++j) {
            position_out[i + j * num_atom_bravais] = position[i][j];
        }
    }

    // 创建输出 types 数组 (num_atom_bravais int 数组)
    plhs[2] = mxCreateNumericMatrix(num_atom_bravais, 1, mxINT32_CLASS, mxREAL);
    int* types_out = static_cast<int*>(mxGetData(plhs[2]));
    for (int i = 0; i < num_atom_bravais; ++i) {
        types_out[i] = types[i];
    }

    // 创建输出 num_atom_bravais 标量
    plhs[3] = mxCreateDoubleScalar(static_cast<double>(num_atom_bravais));
}

// [lattice, result] = symspg('spg_delaunay_reduce', lattice, symprec)
void SpglibFunctions::spg_delaunay_reduce_mex(int nlhs, mxArray* plhs[],
                                              int nrhs, mxArray const* prhs[]) {
    /*
     int spg_delaunay_reduce(double lattice[3][3], double const symprec);
    */

    // 验证输入参数数量
    int const expected_number_of_inputs = 2;
    if (nrhs != expected_number_of_inputs) {
        mexErrMsgIdAndTxt(
            "Spglib:invalidNumInputs",
            "Incorrect number of inputs for spg_delaunay_reduce.");
    }

    // 提取和验证 lattice 参数
    if (mxGetM(prhs[0]) != 3 || mxGetN(prhs[0]) != 3) {
        mexErrMsgIdAndTxt("Spglib:invalidLattice",
                          "Lattice must be a 3x3 matrix.");
    }
    double* lattice_ptr = mxGetPr(prhs[0]);
    double lattice[3][3];
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            lattice[i][j] = lattice_ptr[i + 3 * j];
        }
    }

    // 提取和验证 symprec 参数
    double symprec = mxGetScalar(prhs[1]);

    // 调用 spg_delaunay_reduce
    int result = spg_delaunay_reduce(lattice, symprec);

    // 创建输出 lattice 数组 (3x3 double 数组)
    plhs[0] = mxCreateDoubleMatrix(3, 3, mxREAL);
    double* lattice_out = mxGetPr(plhs[0]);
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            lattice_out[i + 3 * j] = lattice[i][j];
        }
    }

    // 创建输出 result 标量
    plhs[1] = mxCreateDoubleScalar(static_cast<double>(result));
}

// grid_point_index = symspg('spg_get_grid_point_from_address', grid_address,
// mesh)
void SpglibFunctions::spg_get_grid_point_from_address_mex(
    int nlhs, mxArray* plhs[], int nrhs, mxArray const* prhs[]) {
    /*
     int spg_get_grid_point_from_address(int const grid_address[3], int const
     mesh[3]);
    */

    // 验证输入参数数量
    int const expected_number_of_inputs = 2;
    if (nrhs != expected_number_of_inputs) {
        mexErrMsgIdAndTxt(
            "Spglib:invalidNumInputs",
            "Incorrect number of inputs for spg_get_grid_point_from_address.");
    }

    // 提取和验证 grid_address 参数
    if (mxGetNumberOfElements(prhs[0]) != 3) {
        mexErrMsgIdAndTxt("Spglib:invalidGridAddress",
                          "Grid address must be an array of 3 elements.");
    }
    int* grid_address_ptr = static_cast<int*>(mxGetData(prhs[0]));
    int grid_address[3] = {grid_address_ptr[0], grid_address_ptr[1],
                           grid_address_ptr[2]};

    // 提取和验证 mesh 参数
    if (mxGetNumberOfElements(prhs[1]) != 3) {
        mexErrMsgIdAndTxt("Spglib:invalidMesh",
                          "Mesh must be an array of 3 elements.");
    }
    int* mesh_ptr = static_cast<int*>(mxGetData(prhs[1]));
    int mesh[3] = {mesh_ptr[0], mesh_ptr[1], mesh_ptr[2]};

    // 调用 spg_get_grid_point_from_address
    int grid_point_index = spg_get_grid_point_from_address(grid_address, mesh);

    // 创建输出标量 grid_point_index
    plhs[0] = mxCreateDoubleScalar(static_cast<double>(grid_point_index));
}

// dense_grid_point_index = symspg('spg_get_dense_grid_point_from_address',
// grid_address, mesh)
void SpglibFunctions::spg_get_dense_grid_point_from_address_mex(
    int nlhs, mxArray* plhs[], int nrhs, mxArray const* prhs[]) {
    /*
     size_t spg_get_dense_grid_point_from_address(int const grid_address[3], int
     const mesh[3]);
    */

    // 验证输入参数数量
    int const expected_number_of_inputs = 2;
    if (nrhs != expected_number_of_inputs) {
        mexErrMsgIdAndTxt("Spglib:invalidNumInputs",
                          "Incorrect number of inputs for "
                          "spg_get_dense_grid_point_from_address.");
    }

    // 提取和验证 grid_address 参数
    if (mxGetNumberOfElements(prhs[0]) != 3) {
        mexErrMsgIdAndTxt("Spglib:invalidGridAddress",
                          "Grid address must be an array of 3 elements.");
    }
    int* grid_address_ptr = static_cast<int*>(mxGetData(prhs[0]));
    int grid_address[3] = {grid_address_ptr[0], grid_address_ptr[1],
                           grid_address_ptr[2]};

    // 提取和验证 mesh 参数
    if (mxGetNumberOfElements(prhs[1]) != 3) {
        mexErrMsgIdAndTxt("Spglib:invalidMesh",
                          "Mesh must be an array of 3 elements.");
    }
    int* mesh_ptr = static_cast<int*>(mxGetData(prhs[1]));
    int mesh[3] = {mesh_ptr[0], mesh_ptr[1], mesh_ptr[2]};

    // 调用 spg_get_dense_grid_point_from_address
    size_t dense_grid_point_index =
        spg_get_dense_grid_point_from_address(grid_address, mesh);

    // 创建输出标量 dense_grid_point_index
    plhs[0] = mxCreateDoubleScalar(static_cast<double>(dense_grid_point_index));
}

// [grid_address, ir_mapping_table, num_ir_kpoints] =
// symspg('spg_get_ir_reciprocal_mesh', mesh, is_shift, is_time_reversal,
// lattice, position, types, num_atom, symprec)
void SpglibFunctions::spg_get_ir_reciprocal_mesh_mex(int nlhs, mxArray* plhs[],
                                                     int nrhs,
                                                     mxArray const* prhs[]) {
    /*
     int spg_get_ir_reciprocal_mesh(int grid_address[][3], int
     ir_mapping_table[], int const mesh[3], int const is_shift[3], int const
     is_time_reversal, double const lattice[3][3], double const position[][3],
     int const types[], int const num_atom, double const symprec);
    */

    // 验证输入参数数量
    int const expected_number_of_inputs = 8;
    if (nrhs != expected_number_of_inputs) {
        mexErrMsgIdAndTxt(
            "Spglib:invalidNumInputs",
            "Incorrect number of inputs for spg_get_ir_reciprocal_mesh.");
    }

    // 提取和验证 mesh 参数
    if (mxGetNumberOfElements(prhs[0]) != 3) {
        mexErrMsgIdAndTxt("Spglib:invalidMesh",
                          "Mesh must be an array of 3 elements.");
    }
    int* mesh_ptr = static_cast<int*>(mxGetData(prhs[0]));
    int mesh[3] = {mesh_ptr[0], mesh_ptr[1], mesh_ptr[2]};

    // 提取和验证 is_shift 参数
    if (mxGetNumberOfElements(prhs[1]) != 3) {
        mexErrMsgIdAndTxt("Spglib:invalidShift",
                          "is_shift must be an array of 3 elements.");
    }
    int* is_shift_ptr = static_cast<int*>(mxGetData(prhs[1]));
    int is_shift[3] = {is_shift_ptr[0], is_shift_ptr[1], is_shift_ptr[2]};

    // 提取 is_time_reversal 参数
    int is_time_reversal = static_cast<int>(mxGetScalar(prhs[2]));

    // 提取和验证 lattice 参数
    if (mxGetM(prhs[3]) != 3 || mxGetN(prhs[3]) != 3) {
        mexErrMsgIdAndTxt("Spglib:invalidLattice",
                          "Lattice must be a 3x3 matrix.");
    }
    double* lattice_ptr = mxGetPr(prhs[3]);
    double lattice[3][3];
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            lattice[i][j] = lattice_ptr[i + 3 * j];
        }
    }

    // 提取和验证 position 参数
    mwSize num_atom = mxGetM(prhs[4]);
    if (mxGetN(prhs[4]) != 3) {
        mexErrMsgIdAndTxt("Spglib:invalidPosition",
                          "Position must be an Nx3 matrix.");
    }
    double* position_ptr = mxGetPr(prhs[4]);
    double position[num_atom][3];
    for (mwSize i = 0; i < num_atom; ++i) {
        for (int j = 0; j < 3; ++j) {
            position[i][j] = position_ptr[i + j * num_atom];
        }
    }

    // 提取和验证 types 参数
    if (mxGetNumberOfElements(prhs[5]) != num_atom) {
        mexErrMsgIdAndTxt("Spglib:invalidTypes",
                          "Types array size must match the number of atoms.");
    }
    int* types_ptr = static_cast<int*>(mxGetData(prhs[5]));
    int types[num_atom];
    for (mwSize i = 0; i < num_atom; ++i) {
        types[i] = types_ptr[i];
    }

    // 提取 symprec 参数
    double symprec = mxGetScalar(prhs[7]);

    // 创建 grid_address 数组 (最大可能数目的网格点)
    int grid_address[mesh[0] * mesh[1] * mesh[2]][3];
    int ir_mapping_table[mesh[0] * mesh[1] * mesh[2]];

    // 调用 spg_get_ir_reciprocal_mesh
    int num_ir_kpoints = spg_get_ir_reciprocal_mesh(
        grid_address, ir_mapping_table, mesh, is_shift, is_time_reversal,
        lattice, position, types, num_atom, symprec);

    // 创建输出 grid_address 数组 (num_ir_kpoints x 3 double 数组)
    plhs[0] = mxCreateDoubleMatrix(num_ir_kpoints, 3, mxREAL);
    double* grid_address_out = mxGetPr(plhs[0]);
    for (int i = 0; i < num_ir_kpoints; ++i) {
        for (int j = 0; j < 3; ++j) {
            grid_address_out[i + j * num_ir_kpoints] =
                static_cast<double>(grid_address[i][j]);
        }
    }

    // 创建输出 ir_mapping_table 数组 (num_ir_kpoints int 数组)
    plhs[1] = mxCreateNumericMatrix(num_ir_kpoints, 1, mxINT32_CLASS, mxREAL);
    int* ir_mapping_table_out = static_cast<int*>(mxGetData(plhs[1]));
    for (int i = 0; i < num_ir_kpoints; ++i) {
        ir_mapping_table_out[i] = ir_mapping_table[i];
    }

    // 创建输出 num_ir_kpoints 标量
    plhs[2] = mxCreateDoubleScalar(static_cast<double>(num_ir_kpoints));
}

// [grid_address, ir_mapping_table, num_ir_kpoints] =
// symspg('spg_get_dense_ir_reciprocal_mesh', mesh, is_shift, is_time_reversal,
// lattice, position, types, num_atom, symprec)
void SpglibFunctions::spg_get_dense_ir_reciprocal_mesh_mex(
    int nlhs, mxArray* plhs[], int nrhs, mxArray const* prhs[]) {
    /*
     size_t spg_get_dense_ir_reciprocal_mesh(int grid_address[][3], size_t
     ir_mapping_table[], int const mesh[3], int const is_shift[3], int const
     is_time_reversal, double const lattice[3][3], double const position[][3],
     int const types[], int const num_atom, double const symprec);
    */

    // 验证输入参数数量
    int const expected_number_of_inputs = 8;
    if (nrhs != expected_number_of_inputs) {
        mexErrMsgIdAndTxt(
            "Spglib:invalidNumInputs",
            "Incorrect number of inputs for spg_get_dense_ir_reciprocal_mesh.");
    }

    // 提取和验证 mesh 参数
    if (mxGetNumberOfElements(prhs[0]) != 3) {
        mexErrMsgIdAndTxt("Spglib:invalidMesh",
                          "Mesh must be an array of 3 elements.");
    }
    int* mesh_ptr = static_cast<int*>(mxGetData(prhs[0]));
    int mesh[3] = {mesh_ptr[0], mesh_ptr[1], mesh_ptr[2]};

    // 提取和验证 is_shift 参数
    if (mxGetNumberOfElements(prhs[1]) != 3) {
        mexErrMsgIdAndTxt("Spglib:invalidShift",
                          "is_shift must be an array of 3 elements.");
    }
    int* is_shift_ptr = static_cast<int*>(mxGetData(prhs[1]));
    int is_shift[3] = {is_shift_ptr[0], is_shift_ptr[1], is_shift_ptr[2]};

    // 提取 is_time_reversal 参数
    int is_time_reversal = static_cast<int>(mxGetScalar(prhs[2]));

    // 提取和验证 lattice 参数
    if (mxGetM(prhs[3]) != 3 || mxGetN(prhs[3]) != 3) {
        mexErrMsgIdAndTxt("Spglib:invalidLattice",
                          "Lattice must be a 3x3 matrix.");
    }
    double* lattice_ptr = mxGetPr(prhs[3]);
    double lattice[3][3];
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            lattice[i][j] = lattice_ptr[i + 3 * j];
        }
    }

    // 提取和验证 position 参数
    mwSize num_atom = mxGetM(prhs[4]);
    if (mxGetN(prhs[4]) != 3) {
        mexErrMsgIdAndTxt("Spglib:invalidPosition",
                          "Position must be an Nx3 matrix.");
    }
    double* position_ptr = mxGetPr(prhs[4]);
    double position[num_atom][3];
    for (mwSize i = 0; i < num_atom; ++i) {
        for (int j = 0; j < 3; ++j) {
            position[i][j] = position_ptr[i + j * num_atom];
        }
    }

    // 提取和验证 types 参数
    if (mxGetNumberOfElements(prhs[5]) != num_atom) {
        mexErrMsgIdAndTxt("Spglib:invalidTypes",
                          "Types array size must match the number of atoms.");
    }
    int* types_ptr = static_cast<int*>(mxGetData(prhs[5]));
    int types[num_atom];
    for (mwSize i = 0; i < num_atom; ++i) {
        types[i] = types_ptr[i];
    }

    // 提取 symprec 参数
    double symprec = mxGetScalar(prhs[7]);

    // 创建 grid_address 数组 (最大可能数目的网格点)
    size_t num_total_grid_points =
        static_cast<size_t>(mesh[0] * mesh[1] * mesh[2]);
    int grid_address[num_total_grid_points][3];
    size_t ir_mapping_table[num_total_grid_points];

    // 调用 spg_get_dense_ir_reciprocal_mesh
    size_t num_ir_kpoints = spg_get_dense_ir_reciprocal_mesh(
        grid_address, ir_mapping_table, mesh, is_shift, is_time_reversal,
        lattice, position, types, num_atom, symprec);

    // 创建输出 grid_address 数组 (num_total_grid_points x 3 double 数组)
    plhs[0] = mxCreateDoubleMatrix(num_total_grid_points, 3, mxREAL);
    double* grid_address_out = mxGetPr(plhs[0]);
    for (size_t i = 0; i < num_total_grid_points; ++i) {
        for (int j = 0; j < 3; ++j) {
            grid_address_out[i + j * num_total_grid_points] =
                static_cast<double>(grid_address[i][j]);
        }
    }

    // 创建输出 ir_mapping_table 数组 (num_total_grid_points size_t 数组)
    plhs[1] =
        mxCreateNumericMatrix(num_total_grid_points, 1, mxUINT64_CLASS, mxREAL);
    size_t* ir_mapping_table_out = static_cast<size_t*>(mxGetData(plhs[1]));
    for (size_t i = 0; i < num_total_grid_points; ++i) {
        ir_mapping_table_out[i] = ir_mapping_table[i];
    }

    // 创建输出 num_ir_kpoints 标量
    plhs[2] = mxCreateDoubleScalar(static_cast<double>(num_ir_kpoints));
}

// [grid_address, ir_mapping_table, num_ir_kpoints] =
// symspg('spg_get_stabilized_reciprocal_mesh', mesh, is_shift,
// is_time_reversal, num_rot, rotations, num_q, qpoints)
void SpglibFunctions::spg_get_stabilized_reciprocal_mesh_mex(
    int nlhs, mxArray* plhs[], int nrhs, mxArray const* prhs[]) {
    /*
     int spg_get_stabilized_reciprocal_mesh(int grid_address[][3], int
     ir_mapping_table[], int const mesh[3], int const is_shift[3], int const
     is_time_reversal, int const num_rot, int const rotations[][3][3], int const
     num_q, double const qpoints[][3]);
    */

    // 验证输入参数数量
    int const expected_number_of_inputs = 7;
    if (nrhs != expected_number_of_inputs) {
        mexErrMsgIdAndTxt("Spglib:invalidNumInputs",
                          "Incorrect number of inputs for "
                          "spg_get_stabilized_reciprocal_mesh.");
    }

    // 提取和验证 mesh 参数
    if (mxGetNumberOfElements(prhs[0]) != 3) {
        mexErrMsgIdAndTxt("Spglib:invalidMesh",
                          "Mesh must be an array of 3 elements.");
    }
    int* mesh_ptr = static_cast<int*>(mxGetData(prhs[0]));
    int mesh[3] = {mesh_ptr[0], mesh_ptr[1], mesh_ptr[2]};

    // 提取和验证 is_shift 参数
    if (mxGetNumberOfElements(prhs[1]) != 3) {
        mexErrMsgIdAndTxt("Spglib:invalidShift",
                          "is_shift must be an array of 3 elements.");
    }
    int* is_shift_ptr = static_cast<int*>(mxGetData(prhs[1]));
    int is_shift[3] = {is_shift_ptr[0], is_shift_ptr[1], is_shift_ptr[2]};

    // 提取 is_time_reversal 参数
    int is_time_reversal = static_cast<int>(mxGetScalar(prhs[2]));

    // 提取和验证 num_rot 参数
    int num_rot = static_cast<int>(mxGetScalar(prhs[3]));

    // 提取和验证 rotations 参数
    if (mxGetNumberOfElements(prhs[4]) != num_rot * 9) {
        mexErrMsgIdAndTxt("Spglib:invalidRotations",
                          "Rotations must be a num_rot x 3 x 3 array.");
    }
    int* rotations_ptr = static_cast<int*>(mxGetData(prhs[4]));
    int rotations[num_rot][3][3];
    for (int k = 0; k < num_rot; ++k) {
        for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < 3; ++j) {
                rotations[k][i][j] =
                    rotations_ptr[k + i * num_rot + j * num_rot * 3];
            }
        }
    }

    // 提取和验证 num_q 参数
    int num_q = static_cast<int>(mxGetScalar(prhs[5]));

    // 提取和验证 qpoints 参数
    if (mxGetM(prhs[6]) != num_q || mxGetN(prhs[6]) != 3) {
        mexErrMsgIdAndTxt("Spglib:invalidQPoints",
                          "Qpoints must be a num_q x 3 array.");
    }
    double* qpoints_ptr = mxGetPr(prhs[6]);
    double qpoints[num_q][3];
    for (int i = 0; i < num_q; ++i) {
        for (int j = 0; j < 3; ++j) {
            qpoints[i][j] = qpoints_ptr[i + j * num_q];
        }
    }

    // 创建 grid_address 数组 (最大可能数目的网格点)
    int num_total_grid_points = mesh[0] * mesh[1] * mesh[2];
    int grid_address[num_total_grid_points][3];
    int ir_mapping_table[num_total_grid_points];

    // 调用 spg_get_stabilized_reciprocal_mesh
    int num_ir_kpoints = spg_get_stabilized_reciprocal_mesh(
        grid_address, ir_mapping_table, mesh, is_shift, is_time_reversal,
        num_rot, rotations, num_q, qpoints);

    // 创建输出 grid_address 数组 (num_total_grid_points x 3 double 数组)
    plhs[0] = mxCreateDoubleMatrix(num_total_grid_points, 3, mxREAL);
    double* grid_address_out = mxGetPr(plhs[0]);
    for (int i = 0; i < num_total_grid_points; ++i) {
        for (int j = 0; j < 3; ++j) {
            grid_address_out[i + j * num_total_grid_points] =
                static_cast<double>(grid_address[i][j]);
        }
    }

    // 创建输出 ir_mapping_table 数组 (num_total_grid_points int 数组)
    plhs[1] =
        mxCreateNumericMatrix(num_total_grid_points, 1, mxINT32_CLASS, mxREAL);
    int* ir_mapping_table_out = static_cast<int*>(mxGetData(plhs[1]));
    for (int i = 0; i < num_total_grid_points; ++i) {
        ir_mapping_table_out[i] = ir_mapping_table[i];
    }

    // 创建输出 num_ir_kpoints 标量
    plhs[2] = mxCreateDoubleScalar(static_cast<double>(num_ir_kpoints));
}

// [grid_address, ir_mapping_table, num_ir_kpoints] =
// symspg('spg_get_dense_stabilized_reciprocal_mesh', mesh, is_shift,
// is_time_reversal, num_rot, rotations, num_q, qpoints)
void SpglibFunctions::spg_get_dense_stabilized_reciprocal_mesh_mex(
    int nlhs, mxArray* plhs[], int nrhs, mxArray const* prhs[]) {
    /*
     size_t spg_get_dense_stabilized_reciprocal_mesh(int grid_address[][3],
     size_t ir_mapping_table[], int const mesh[3], int const is_shift[3], int
     const is_time_reversal, int const num_rot, int const rotations[][3][3], int
     const num_q, double const qpoints[][3]);
    */

    // 验证输入参数数量
    int const expected_number_of_inputs = 7;
    if (nrhs != expected_number_of_inputs) {
        mexErrMsgIdAndTxt("Spglib:invalidNumInputs",
                          "Incorrect number of inputs for "
                          "spg_get_dense_stabilized_reciprocal_mesh.");
    }

    // 提取和验证 mesh 参数
    if (mxGetNumberOfElements(prhs[0]) != 3) {
        mexErrMsgIdAndTxt("Spglib:invalidMesh",
                          "Mesh must be an array of 3 elements.");
    }
    int* mesh_ptr = static_cast<int*>(mxGetData(prhs[0]));
    int mesh[3] = {mesh_ptr[0], mesh_ptr[1], mesh_ptr[2]};

    // 提取和验证 is_shift 参数
    if (mxGetNumberOfElements(prhs[1]) != 3) {
        mexErrMsgIdAndTxt("Spglib:invalidShift",
                          "is_shift must be an array of 3 elements.");
    }
    int* is_shift_ptr = static_cast<int*>(mxGetData(prhs[1]));
    int is_shift[3] = {is_shift_ptr[0], is_shift_ptr[1], is_shift_ptr[2]};

    // 提取 is_time_reversal 参数
    int is_time_reversal = static_cast<int>(mxGetScalar(prhs[2]));

    // 提取和验证 num_rot 参数
    int num_rot = static_cast<int>(mxGetScalar(prhs[3]));

    // 提取和验证 rotations 参数
    if (mxGetNumberOfElements(prhs[4]) != num_rot * 9) {
        mexErrMsgIdAndTxt("Spglib:invalidRotations",
                          "Rotations must be a num_rot x 3 x 3 array.");
    }
    int* rotations_ptr = static_cast<int*>(mxGetData(prhs[4]));
    int rotations[num_rot][3][3];
    for (int k = 0; k < num_rot; ++k) {
        for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < 3; ++j) {
                rotations[k][i][j] =
                    rotations_ptr[k + i * num_rot + j * num_rot * 3];
            }
        }
    }

    // 提取和验证 num_q 参数
    int num_q = static_cast<int>(mxGetScalar(prhs[5]));

    // 提取和验证 qpoints 参数
    if (mxGetM(prhs[6]) != num_q || mxGetN(prhs[6]) != 3) {
        mexErrMsgIdAndTxt("Spglib:invalidQPoints",
                          "Qpoints must be a num_q x 3 array.");
    }
    double* qpoints_ptr = mxGetPr(prhs[6]);
    double qpoints[num_q][3];
    for (int i = 0; i < num_q; ++i) {
        for (int j = 0; j < 3; ++j) {
            qpoints[i][j] = qpoints_ptr[i + j * num_q];
        }
    }

    // 创建 grid_address 数组 (最大可能数目的网格点)
    size_t num_total_grid_points =
        static_cast<size_t>(mesh[0] * mesh[1] * mesh[2]);
    int grid_address[num_total_grid_points][3];
    size_t ir_mapping_table[num_total_grid_points];

    // 调用 spg_get_dense_stabilized_reciprocal_mesh
    size_t num_ir_kpoints = spg_get_dense_stabilized_reciprocal_mesh(
        grid_address, ir_mapping_table, mesh, is_shift, is_time_reversal,
        num_rot, rotations, num_q, qpoints);

    // 创建输出 grid_address 数组 (num_total_grid_points x 3 double 数组)
    plhs[0] = mxCreateDoubleMatrix(num_total_grid_points, 3, mxREAL);
    double* grid_address_out = mxGetPr(plhs[0]);
    for (size_t i = 0; i < num_total_grid_points; ++i) {
        for (int j = 0; j < 3; ++j) {
            grid_address_out[i + j * num_total_grid_points] =
                static_cast<double>(grid_address[i][j]);
        }
    }

    // 创建输出 ir_mapping_table 数组 (num_total_grid_points size_t 数组)
    plhs[1] =
        mxCreateNumericMatrix(num_total_grid_points, 1, mxUINT64_CLASS, mxREAL);
    size_t* ir_mapping_table_out = static_cast<size_t*>(mxGetData(plhs[1]));
    for (size_t i = 0; i < num_total_grid_points; ++i) {
        ir_mapping_table_out[i] = ir_mapping_table[i];
    }

    // 创建输出 num_ir_kpoints 标量
    plhs[2] = mxCreateDoubleScalar(static_cast<double>(num_ir_kpoints));
}

// rot_grid_points = symspg('spg_get_dense_grid_points_by_rotations',
// address_orig, num_rot, rot_reciprocal, mesh, is_shift)
void SpglibFunctions::spg_get_dense_grid_points_by_rotations_mex(
    int nlhs, mxArray* plhs[], int nrhs, mxArray const* prhs[]) {
    /*
     void spg_get_dense_grid_points_by_rotations(size_t rot_grid_points[], int
     const address_orig[3], int const num_rot, int const rot_reciprocal[][3][3],
                                                 int const mesh[3], int const
     is_shift[3]);
    */

    // 验证输入参数数量
    int const expected_number_of_inputs = 5;
    if (nrhs != expected_number_of_inputs) {
        mexErrMsgIdAndTxt("Spglib:invalidNumInputs",
                          "Incorrect number of inputs for "
                          "spg_get_dense_grid_points_by_rotations.");
    }

    // 提取和验证 address_orig 参数
    if (mxGetNumberOfElements(prhs[0]) != 3) {
        mexErrMsgIdAndTxt("Spglib:invalidAddressOrig",
                          "Address_orig must be an array of 3 elements.");
    }
    int* address_orig_ptr = static_cast<int*>(mxGetData(prhs[0]));
    int address_orig[3] = {address_orig_ptr[0], address_orig_ptr[1],
                           address_orig_ptr[2]};

    // 提取和验证 num_rot 参数
    int num_rot = static_cast<int>(mxGetScalar(prhs[1]));

    // 提取和验证 rot_reciprocal 参数
    if (mxGetNumberOfElements(prhs[2]) != num_rot * 9) {
        mexErrMsgIdAndTxt("Spglib:invalidRotReciprocal",
                          "Rot_reciprocal must be a num_rot x 3 x 3 array.");
    }
    int* rot_reciprocal_ptr = static_cast<int*>(mxGetData(prhs[2]));
    int rot_reciprocal[num_rot][3][3];
    for (int k = 0; k < num_rot; ++k) {
        for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < 3; ++j) {
                rot_reciprocal[k][i][j] =
                    rot_reciprocal_ptr[k + i * num_rot + j * num_rot * 3];
            }
        }
    }

    // 提取和验证 mesh 参数
    if (mxGetNumberOfElements(prhs[3]) != 3) {
        mexErrMsgIdAndTxt("Spglib:invalidMesh",
                          "Mesh must be an array of 3 elements.");
    }
    int* mesh_ptr = static_cast<int*>(mxGetData(prhs[3]));
    int mesh[3] = {mesh_ptr[0], mesh_ptr[1], mesh_ptr[2]};

    // 提取和验证 is_shift 参数
    if (mxGetNumberOfElements(prhs[4]) != 3) {
        mexErrMsgIdAndTxt("Spglib:invalidShift",
                          "is_shift must be an array of 3 elements.");
    }
    int* is_shift_ptr = static_cast<int*>(mxGetData(prhs[4]));
    int is_shift[3] = {is_shift_ptr[0], is_shift_ptr[1], is_shift_ptr[2]};

    // 创建输出 rot_grid_points 数组 (num_rot size_t 数组)
    plhs[0] = mxCreateNumericMatrix(num_rot, 1, mxUINT64_CLASS, mxREAL);
    size_t* rot_grid_points = static_cast<size_t*>(mxGetData(plhs[0]));

    // 调用 spg_get_dense_grid_points_by_rotations
    spg_get_dense_grid_points_by_rotations(
        rot_grid_points, address_orig, num_rot, rot_reciprocal, mesh, is_shift);
}

// rot_grid_points = symspg('spg_get_dense_BZ_grid_points_by_rotations',
// address_orig, num_rot, rot_reciprocal, mesh, is_shift, bz_map)
void SpglibFunctions::spg_get_dense_BZ_grid_points_by_rotations_mex(
    int nlhs, mxArray* plhs[], int nrhs, mxArray const* prhs[]) {
    /*
     void spg_get_dense_BZ_grid_points_by_rotations(size_t rot_grid_points[],
     int const address_orig[3], int const num_rot, int const
     rot_reciprocal[][3][3], int const mesh[3], int const is_shift[3], size_t
     const bz_map[]);
    */

    // 验证输入参数数量
    int const expected_number_of_inputs = 6;
    if (nrhs != expected_number_of_inputs) {
        mexErrMsgIdAndTxt("Spglib:invalidNumInputs",
                          "Incorrect number of inputs for "
                          "spg_get_dense_BZ_grid_points_by_rotations.");
    }

    // 提取和验证 address_orig 参数
    if (mxGetNumberOfElements(prhs[0]) != 3) {
        mexErrMsgIdAndTxt("Spglib:invalidAddressOrig",
                          "Address_orig must be an array of 3 elements.");
    }
    int* address_orig_ptr = static_cast<int*>(mxGetData(prhs[0]));
    int address_orig[3] = {address_orig_ptr[0], address_orig_ptr[1],
                           address_orig_ptr[2]};

    // 提取和验证 num_rot 参数
    int num_rot = static_cast<int>(mxGetScalar(prhs[1]));

    // 提取和验证 rot_reciprocal 参数
    if (mxGetNumberOfElements(prhs[2]) != num_rot * 9) {
        mexErrMsgIdAndTxt("Spglib:invalidRotReciprocal",
                          "Rot_reciprocal must be a num_rot x 3 x 3 array.");
    }
    int* rot_reciprocal_ptr = static_cast<int*>(mxGetData(prhs[2]));
    int rot_reciprocal[num_rot][3][3];
    for (int k = 0; k < num_rot; ++k) {
        for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < 3; ++j) {
                rot_reciprocal[k][i][j] =
                    rot_reciprocal_ptr[k + i * num_rot + j * num_rot * 3];
            }
        }
    }

    // 提取和验证 mesh 参数
    if (mxGetNumberOfElements(prhs[3]) != 3) {
        mexErrMsgIdAndTxt("Spglib:invalidMesh",
                          "Mesh must be an array of 3 elements.");
    }
    int* mesh_ptr = static_cast<int*>(mxGetData(prhs[3]));
    int mesh[3] = {mesh_ptr[0], mesh_ptr[1], mesh_ptr[2]};

    // 提取和验证 is_shift 参数
    if (mxGetNumberOfElements(prhs[4]) != 3) {
        mexErrMsgIdAndTxt("Spglib:invalidShift",
                          "is_shift must be an array of 3 elements.");
    }
    int* is_shift_ptr = static_cast<int*>(mxGetData(prhs[4]));
    int is_shift[3] = {is_shift_ptr[0], is_shift_ptr[1], is_shift_ptr[2]};

    // 提取和验证 bz_map 参数
    size_t* bz_map = static_cast<size_t*>(mxGetData(prhs[5]));

    // 创建输出 rot_grid_points 数组 (num_rot size_t 数组)
    plhs[0] = mxCreateNumericMatrix(num_rot, 1, mxUINT64_CLASS, mxREAL);
    size_t* rot_grid_points = static_cast<size_t*>(mxGetData(plhs[0]));

    // 调用 spg_get_dense_BZ_grid_points_by_rotations
    spg_get_dense_BZ_grid_points_by_rotations(rot_grid_points, address_orig,
                                              num_rot, rot_reciprocal, mesh,
                                              is_shift, bz_map);
}

// [bz_grid_address, bz_map, num_ir_grid_points] =
// symspg('spg_relocate_BZ_grid_address', grid_address, mesh, rec_lattice,
// is_shift)
void SpglibFunctions::spg_relocate_BZ_grid_address_mex(int nlhs,
                                                       mxArray* plhs[],
                                                       int nrhs,
                                                       mxArray const* prhs[]) {
    /*
     int spg_relocate_BZ_grid_address(int bz_grid_address[][3], int bz_map[],
                                      int const grid_address[][3], int const
     mesh[3], double const rec_lattice[3][3], int const is_shift[3]);
    */

    // 验证输入参数数量
    int const expected_number_of_inputs = 4;
    if (nrhs != expected_number_of_inputs) {
        mexErrMsgIdAndTxt(
            "Spglib:invalidNumInputs",
            "Incorrect number of inputs for spg_relocate_BZ_grid_address.");
    }

    // 提取和验证 grid_address 参数
    mwSize num_grid_points = mxGetM(prhs[0]);
    if (mxGetN(prhs[0]) != 3) {
        mexErrMsgIdAndTxt("Spglib:invalidGridAddress",
                          "Grid_address must be an Nx3 array.");
    }
    int* grid_address_ptr = static_cast<int*>(mxGetData(prhs[0]));
    int grid_address[num_grid_points][3];
    for (mwSize i = 0; i < num_grid_points; ++i) {
        for (int j = 0; j < 3; ++j) {
            grid_address[i][j] = grid_address_ptr[i + j * num_grid_points];
        }
    }

    // 提取和验证 mesh 参数
    if (mxGetNumberOfElements(prhs[1]) != 3) {
        mexErrMsgIdAndTxt("Spglib:invalidMesh",
                          "Mesh must be an array of 3 elements.");
    }
    int* mesh_ptr = static_cast<int*>(mxGetData(prhs[1]));
    int mesh[3] = {mesh_ptr[0], mesh_ptr[1], mesh_ptr[2]};

    // 提取和验证 rec_lattice 参数
    if (mxGetM(prhs[2]) != 3 || mxGetN(prhs[2]) != 3) {
        mexErrMsgIdAndTxt("Spglib:invalidRecLattice",
                          "Rec_lattice must be a 3x3 matrix.");
    }
    double* rec_lattice_ptr = mxGetPr(prhs[2]);
    double rec_lattice[3][3];
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            rec_lattice[i][j] = rec_lattice_ptr[i + 3 * j];
        }
    }

    // 提取和验证 is_shift 参数
    if (mxGetNumberOfElements(prhs[3]) != 3) {
        mexErrMsgIdAndTxt("Spglib:invalidShift",
                          "is_shift must be an array of 3 elements.");
    }
    int* is_shift_ptr = static_cast<int*>(mxGetData(prhs[3]));
    int is_shift[3] = {is_shift_ptr[0], is_shift_ptr[1], is_shift_ptr[2]};

    // 创建输出 bz_grid_address 数组 (prod(mesh + 1) x 3 int 数组)
    size_t bz_grid_address_size = (mesh[0] + 1) * (mesh[1] + 1) * (mesh[2] + 1);
    plhs[0] =
        mxCreateNumericMatrix(bz_grid_address_size, 3, mxINT32_CLASS, mxREAL);
    int(*bz_grid_address)[3] = static_cast<int(*)[3]>(mxGetData(plhs[0]));

    // 创建输出 bz_map 数组 (prod(mesh * 2) int 数组)
    size_t bz_map_size = mesh[0] * 2 * mesh[1] * 2 * mesh[2] * 2;
    plhs[1] = mxCreateNumericMatrix(bz_map_size, 1, mxINT32_CLASS, mxREAL);
    int* bz_map = static_cast<int*>(mxGetData(plhs[1]));

    // 调用 spg_relocate_BZ_grid_address
    int num_ir_grid_points = spg_relocate_BZ_grid_address(
        bz_grid_address, bz_map, grid_address, mesh, rec_lattice, is_shift);

    // 创建输出 num_ir_grid_points 标量
    plhs[2] = mxCreateDoubleScalar(static_cast<double>(num_ir_grid_points));
}

// [bz_grid_address, bz_map, num_ir_grid_points] =
// symspg('spg_relocate_dense_BZ_grid_address', grid_address, mesh, rec_lattice,
// is_shift)
void SpglibFunctions::spg_relocate_dense_BZ_grid_address_mex(
    int nlhs, mxArray* plhs[], int nrhs, mxArray const* prhs[]) {
    /*
     size_t spg_relocate_dense_BZ_grid_address(int bz_grid_address[][3], size_t
     bz_map[], int const grid_address[][3], int const mesh[3], double const
     rec_lattice[3][3], int const is_shift[3]);
    */

    // 验证输入参数数量
    int const expected_number_of_inputs = 4;
    if (nrhs != expected_number_of_inputs) {
        mexErrMsgIdAndTxt("Spglib:invalidNumInputs",
                          "Incorrect number of inputs for "
                          "spg_relocate_dense_BZ_grid_address.");
    }

    // 提取和验证 grid_address 参数
    mwSize num_grid_points = mxGetM(prhs[0]);
    if (mxGetN(prhs[0]) != 3) {
        mexErrMsgIdAndTxt("Spglib:invalidGridAddress",
                          "Grid_address must be an Nx3 array.");
    }
    int* grid_address_ptr = static_cast<int*>(mxGetData(prhs[0]));
    int grid_address[num_grid_points][3];
    for (mwSize i = 0; i < num_grid_points; ++i) {
        for (int j = 0; j < 3; ++j) {
            grid_address[i][j] = grid_address_ptr[i + j * num_grid_points];
        }
    }

    // 提取和验证 mesh 参数
    if (mxGetNumberOfElements(prhs[1]) != 3) {
        mexErrMsgIdAndTxt("Spglib:invalidMesh",
                          "Mesh must be an array of 3 elements.");
    }
    int* mesh_ptr = static_cast<int*>(mxGetData(prhs[1]));
    int mesh[3] = {mesh_ptr[0], mesh_ptr[1], mesh_ptr[2]};

    // 提取和验证 rec_lattice 参数
    if (mxGetM(prhs[2]) != 3 || mxGetN(prhs[2]) != 3) {
        mexErrMsgIdAndTxt("Spglib:invalidRecLattice",
                          "Rec_lattice must be a 3x3 matrix.");
    }
    double* rec_lattice_ptr = mxGetPr(prhs[2]);
    double rec_lattice[3][3];
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            rec_lattice[i][j] = rec_lattice_ptr[i + 3 * j];
        }
    }

    // 提取和验证 is_shift 参数
    if (mxGetNumberOfElements(prhs[3]) != 3) {
        mexErrMsgIdAndTxt("Spglib:invalidShift",
                          "is_shift must be an array of 3 elements.");
    }
    int* is_shift_ptr = static_cast<int*>(mxGetData(prhs[3]));
    int is_shift[3] = {is_shift_ptr[0], is_shift_ptr[1], is_shift_ptr[2]};

    // 创建输出 bz_grid_address 数组 (prod(mesh + 1) x 3 int 数组)
    size_t bz_grid_address_size = (mesh[0] + 1) * (mesh[1] + 1) * (mesh[2] + 1);
    plhs[0] =
        mxCreateNumericMatrix(bz_grid_address_size, 3, mxINT32_CLASS, mxREAL);
    int(*bz_grid_address)[3] = static_cast<int(*)[3]>(mxGetData(plhs[0]));

    // 创建输出 bz_map 数组 (prod(mesh * 2) size_t 数组)
    size_t bz_map_size = mesh[0] * 2 * mesh[1] * 2 * mesh[2] * 2;
    plhs[1] = mxCreateNumericMatrix(bz_map_size, 1, mxUINT64_CLASS, mxREAL);
    size_t* bz_map = static_cast<size_t*>(mxGetData(plhs[1]));

    // 调用 spg_relocate_dense_BZ_grid_address
    size_t num_ir_grid_points = spg_relocate_dense_BZ_grid_address(
        bz_grid_address, bz_map, grid_address, mesh, rec_lattice, is_shift);

    // 创建输出 num_ir_grid_points 标量
    plhs[2] = mxCreateDoubleScalar(static_cast<double>(num_ir_grid_points));
}

// [lattice, success] = symspg('spg_niggli_reduce', lattice, symprec)
void SpglibFunctions::spg_niggli_reduce_mex(int nlhs, mxArray* plhs[], int nrhs,
                                            mxArray const* prhs[]) {
    /*
     int spg_niggli_reduce(double lattice[3][3], double const symprec);
    */

    // 验证输入参数数量
    int const expected_number_of_inputs = 2;
    if (nrhs != expected_number_of_inputs) {
        mexErrMsgIdAndTxt("Spglib:invalidNumInputs",
                          "Incorrect number of inputs for spg_niggli_reduce.");
    }

    // 提取和验证 lattice 参数
    if (mxGetM(prhs[0]) != 3 || mxGetN(prhs[0]) != 3) {
        mexErrMsgIdAndTxt("Spglib:invalidLattice",
                          "Lattice must be a 3x3 matrix.");
    }
    double* lattice_ptr = mxGetPr(prhs[0]);
    double lattice[3][3];
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            lattice[i][j] = lattice_ptr[i + 3 * j];
        }
    }

    // 提取 symprec 参数
    double symprec = mxGetScalar(prhs[1]);

    // 调用 spg_niggli_reduce
    int success = spg_niggli_reduce(lattice, symprec);

    // 创建输出 lattice 数组 (3x3 double 数组)
    plhs[0] = mxCreateDoubleMatrix(3, 3, mxREAL);
    double* lattice_out = mxGetPr(plhs[0]);
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            lattice_out[i + 3 * j] = lattice[i][j];
        }
    }

    // 创建输出 success 标量
    plhs[1] = mxCreateDoubleScalar(static_cast<double>(success));
}
