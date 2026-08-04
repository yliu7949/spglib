#include <string>
#include <unordered_map>
#include "dataset_converter.hpp"
#include "mex.h"
#include "mex_array.hpp"
#include "mex_buffer.hpp"
#include "spglib.h"

// Define a class containing static methods
class SpglibFunctions {
   public:
    // Declare the static method interface
    static void spg_get_version_mex(int nlhs, mxArray *plhs[], int nrhs,
                                    mxArray const *prhs[]);
    static void spg_get_version_full_mex(int nlhs, mxArray *plhs[], int nrhs,
                                         mxArray const *prhs[]);
    static void spg_get_commit_mex(int nlhs, mxArray *plhs[], int nrhs,
                                   mxArray const *prhs[]);
    static void spg_get_major_version_mex(int nlhs, mxArray *plhs[], int nrhs,
                                          mxArray const *prhs[]);
    static void spg_get_minor_version_mex(int nlhs, mxArray *plhs[], int nrhs,
                                          mxArray const *prhs[]);
    static void spg_get_micro_version_mex(int nlhs, mxArray *plhs[], int nrhs,
                                          mxArray const *prhs[]);
    static void spg_get_error_code_mex(int nlhs, mxArray *plhs[], int nrhs,
                                       mxArray const *prhs[]);
    static void spg_get_error_message_mex(int nlhs, mxArray *plhs[], int nrhs,
                                          mxArray const *prhs[]);
    static void spg_get_dataset_mex(int nlhs, mxArray *plhs[], int nrhs,
                                    mxArray const *prhs[]);
    static void spg_get_magnetic_dataset_mex(int nlhs, mxArray *plhs[],
                                             int nrhs, mxArray const *prhs[]);
    static void spgms_get_magnetic_dataset_mex(int nlhs, mxArray *plhs[],
                                               int nrhs, mxArray const *prhs[]);
    static void spgat_get_dataset_mex(int nlhs, mxArray *plhs[], int nrhs,
                                      mxArray const *prhs[]);
    static void spg_get_dataset_with_hall_number_mex(int nlhs, mxArray *plhs[],
                                                     int nrhs,
                                                     mxArray const *prhs[]);
    static void spgat_get_dataset_with_hall_number_mex(int nlhs,
                                                       mxArray *plhs[],
                                                       int nrhs,
                                                       mxArray const *prhs[]);
    static void spg_get_symmetry_with_collinear_spin_mex(int nlhs,
                                                         mxArray *plhs[],
                                                         int nrhs,
                                                         mxArray const *prhs[]);
    static void spgat_get_symmetry_with_collinear_spin_mex(
        int nlhs, mxArray *plhs[], int nrhs, mxArray const *prhs[]);
    static void spgms_get_symmetry_with_collinear_spin_mex(
        int nlhs, mxArray *plhs[], int nrhs, mxArray const *prhs[]);
    static void spg_get_symmetry_with_site_tensors_mex(int nlhs,
                                                       mxArray *plhs[],
                                                       int nrhs,
                                                       mxArray const *prhs[]);
    static void spgat_get_symmetry_with_site_tensors_mex(int nlhs,
                                                         mxArray *plhs[],
                                                         int nrhs,
                                                         mxArray const *prhs[]);
    static void spgms_get_symmetry_with_site_tensors_mex(int nlhs,
                                                         mxArray *plhs[],
                                                         int nrhs,
                                                         mxArray const *prhs[]);
    static void spg_get_spacegroup_type_from_symmetry_mex(
        int nlhs, mxArray *plhs[], int nrhs, mxArray const *prhs[]);
    static void spg_get_magnetic_spacegroup_type_from_symmetry_mex(
        int nlhs, mxArray *plhs[], int nrhs, mxArray const *prhs[]);
    static void spg_get_pointgroup_mex(int nlhs, mxArray *plhs[], int nrhs,
                                       mxArray const *prhs[]);
    static void spg_get_symmetry_from_database_mex(int nlhs, mxArray *plhs[],
                                                   int nrhs,
                                                   mxArray const *prhs[]);
    static void spg_get_magnetic_symmetry_from_database_mex(
        int nlhs, mxArray *plhs[], int nrhs, mxArray const *prhs[]);
    static void spg_get_spacegroup_type_mex(int nlhs, mxArray *plhs[], int nrhs,
                                            mxArray const *prhs[]);
    static void spg_get_magnetic_spacegroup_type_mex(int nlhs, mxArray *plhs[],
                                                     int nrhs,
                                                     mxArray const *prhs[]);
    static void spg_standardize_cell_mex(int nlhs, mxArray *plhs[], int nrhs,
                                         mxArray const *prhs[]);
    static void spgat_standardize_cell_mex(int nlhs, mxArray *plhs[], int nrhs,
                                           mxArray const *prhs[]);
    static void spg_find_primitive_mex(int nlhs, mxArray *plhs[], int nrhs,
                                       mxArray const *prhs[]);
    static void spgat_find_primitive_mex(int nlhs, mxArray *plhs[], int nrhs,
                                         mxArray const *prhs[]);
    static void spg_refine_cell_mex(int nlhs, mxArray *plhs[], int nrhs,
                                    mxArray const *prhs[]);
    static void spgat_refine_cell_mex(int nlhs, mxArray *plhs[], int nrhs,
                                      mxArray const *prhs[]);
    static void spg_delaunay_reduce_mex(int nlhs, mxArray *plhs[], int nrhs,
                                        mxArray const *prhs[]);
    static void spg_get_grid_point_from_address_mex(int nlhs, mxArray *plhs[],
                                                    int nrhs,
                                                    mxArray const *prhs[]);
    static void spg_get_dense_grid_point_from_address_mex(
        int nlhs, mxArray *plhs[], int nrhs, mxArray const *prhs[]);
    static void spg_get_ir_reciprocal_mesh_mex(int nlhs, mxArray *plhs[],
                                               int nrhs, mxArray const *prhs[]);
    static void spg_get_dense_ir_reciprocal_mesh_mex(int nlhs, mxArray *plhs[],
                                                     int nrhs,
                                                     mxArray const *prhs[]);
    static void spg_get_stabilized_reciprocal_mesh_mex(int nlhs,
                                                       mxArray *plhs[],
                                                       int nrhs,
                                                       mxArray const *prhs[]);
    static void spg_get_dense_stabilized_reciprocal_mesh_mex(
        int nlhs, mxArray *plhs[], int nrhs, mxArray const *prhs[]);
    static void spg_get_dense_grid_points_by_rotations_mex(
        int nlhs, mxArray *plhs[], int nrhs, mxArray const *prhs[]);
    static void spg_get_dense_BZ_grid_points_by_rotations_mex(
        int nlhs, mxArray *plhs[], int nrhs, mxArray const *prhs[]);
    static void spg_relocate_BZ_grid_address_mex(int nlhs, mxArray *plhs[],
                                                 int nrhs,
                                                 mxArray const *prhs[]);
    static void spg_relocate_dense_BZ_grid_address_mex(int nlhs,
                                                       mxArray *plhs[],
                                                       int nrhs,
                                                       mxArray const *prhs[]);
    static void spg_niggli_reduce_mex(int nlhs, mxArray *plhs[], int nrhs,
                                      mxArray const *prhs[]);
};

// Define the function-pointer type for static methods
typedef void (*SpglibFunction)(int, mxArray *[], int, mxArray const *[]);

namespace {

void throwLastSpglibError(char const *operation) {
    SpglibError const error = spg_get_error_code();
    mexErrMsgIdAndTxt("Spglib:spglibError", "%s: %s", operation,
                      spg_get_error_message(error));
}

void validateNumAtoms(mxArray const *value, mwSize const expected) {
    if (!mxIsNumeric(value) || mxIsComplex(value) ||
        mxGetNumberOfElements(value) != 1 ||
        mxGetScalar(value) != static_cast<double>(expected)) {
        mexErrMsgIdAndTxt(
            "Spglib:invalidNumAtoms",
            "num_atom must be a real numeric scalar matching the number of "
            "position rows.");
    }
}

}  // namespace

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

// Main MEX entry point
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, mxArray const *prhs[]) {
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

    // Map function names to static methods
    static std::unordered_map<std::string, SpglibFunction> const function_map =
        {{"spg_get_version", SpglibFunctions::spg_get_version_mex},
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
         {"spgat_standardize_cell",
          SpglibFunctions::spgat_standardize_cell_mex},
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

    // Find and invoke the requested function
    auto it = function_map.find(function_name);
    if (it != function_map.end()) {
        it->second(nlhs, plhs, nrhs - 1, prhs + 1);
    } else {
        mexErrMsgIdAndTxt("Spglib:invalidFunction", "Unknown function name.");
    }
}

// version = symspg('spg_get_version')
void SpglibFunctions::spg_get_version_mex(int nlhs, mxArray *plhs[], int nrhs,
                                          mxArray const *prhs[]) {
    if (nrhs != 0) {
        mexErrMsgIdAndTxt("Spglib:invalidNumInputs",
                          "spg_get_version does not take any input arguments.");
    }

    char const *version = spg_get_version();
    plhs[0] = mxCreateString(version);
}

// version = symspg('spg_get_version_full')
void SpglibFunctions::spg_get_version_full_mex(int nlhs, mxArray *plhs[],
                                               int nrhs,
                                               mxArray const *prhs[]) {
    if (nrhs != 0) {
        mexErrMsgIdAndTxt(
            "Spglib:invalidNumInputs",
            "spg_get_version_full does not take any input arguments.");
    }

    char const *version_full = spg_get_version_full();
    plhs[0] = mxCreateString(version_full);
}

// commit = symspg('spg_get_commit')
void SpglibFunctions::spg_get_commit_mex(int nlhs, mxArray *plhs[], int nrhs,
                                         mxArray const *prhs[]) {
    if (nrhs != 0) {
        mexErrMsgIdAndTxt("Spglib:invalidNumInputs",
                          "spg_get_commit does not take any input arguments.");
    }

    char const *commit = spg_get_commit();
    plhs[0] = mxCreateString(commit);
}

// version = symspg('spg_get_major_version')
void SpglibFunctions::spg_get_major_version_mex(int nlhs, mxArray *plhs[],
                                                int nrhs,
                                                mxArray const *prhs[]) {
    if (nrhs != 0) {
        mexErrMsgIdAndTxt(
            "Spglib:invalidNumInputs",
            "spg_get_major_version does not take any input arguments.");
    }
    int major_version = spg_get_major_version();
    plhs[0] = mxCreateDoubleScalar(static_cast<double>(major_version));
}

// version = symspg('spg_get_minor_version')
void SpglibFunctions::spg_get_minor_version_mex(int nlhs, mxArray *plhs[],
                                                int nrhs,
                                                mxArray const *prhs[]) {
    if (nrhs != 0) {
        mexErrMsgIdAndTxt(
            "Spglib:invalidNumInputs",
            "spg_get_minor_version does not take any input arguments.");
    }
    int minor_version = spg_get_minor_version();
    plhs[0] = mxCreateDoubleScalar(static_cast<double>(minor_version));
}

// version = symspg('spg_get_micro_version')
void SpglibFunctions::spg_get_micro_version_mex(int nlhs, mxArray *plhs[],
                                                int nrhs,
                                                mxArray const *prhs[]) {
    if (nrhs != 0) {
        mexErrMsgIdAndTxt(
            "Spglib:invalidNumInputs",
            "spg_get_micro_version does not take any input arguments.");
    }
    int micro_version = spg_get_micro_version();
    plhs[0] = mxCreateDoubleScalar(static_cast<double>(micro_version));
}

// error_code = symspg('spg_get_error_code')
void SpglibFunctions::spg_get_error_code_mex(int nlhs, mxArray *plhs[],
                                             int nrhs, mxArray const *prhs[]) {
    /*
     SpglibError spg_get_error_code(void);
    */

    // Validate the number of input arguments
    if (nrhs != 0) {
        mexErrMsgIdAndTxt("Spglib:invalidNumInputs",
                          "No inputs expected for spg_get_error_code.");
    }

    // Call spg_get_error_code
    SpglibError error_code = spg_get_error_code();

    // Create the output error_code scalar
    plhs[0] = mxCreateDoubleScalar(static_cast<double>(error_code));
}

// error_message = symspg('spg_get_error_message', error_code)
void SpglibFunctions::spg_get_error_message_mex(int nlhs, mxArray *plhs[],
                                                int nrhs,
                                                mxArray const *prhs[]) {
    /*
     const char *spg_get_error_message(SpglibError spglib_error);
    */

    // Validate the number of input arguments
    int const expected_number_of_inputs = 1;
    if (nrhs != expected_number_of_inputs) {
        mexErrMsgIdAndTxt(
            "Spglib:invalidNumInputs",
            "Incorrect number of inputs for spg_get_error_message.");
    }

    // Extract and validate the spglib_error argument
    SpglibError spglib_error = static_cast<SpglibError>(mxGetScalar(prhs[0]));

    // Call spg_get_error_message
    char const *error_message = spg_get_error_message(spglib_error);

    // Create the output string
    plhs[0] = mxCreateString(error_message);
}

// dataset = symspg('spg_get_dataset', lattice, position, types, num_atom,
// symprec)
void SpglibFunctions::spg_get_dataset_mex(int nlhs, mxArray *plhs[], int nrhs,
                                          mxArray const *prhs[]) {
    /*
     SpglibDataset * spg_get_dataset(const double lattice[3][3],
                           const double position[][3],
                           const int types[],
                           const int num_atom,
                           const double symprec);
     */

    // Validate the number of input arguments
    int const expected_number_of_inputs = 5;
    if (nrhs != expected_number_of_inputs) {
        mexErrMsgIdAndTxt("Spglib:invalidNumInputs",
                          "Incorrect number of inputs for spg_get_dataset.");
    }

    // Extract and validate the lattice, position, types, and symprec arguments
    double lattice[3][3];
    if (mxGetM(prhs[0]) != 3 || mxGetN(prhs[0]) != 3) {
        mexErrMsgIdAndTxt("Spglib:invalidLattice",
                          "Lattice matrix must be 3x3.");
    }

    double *lattice_ptr = mxGetPr(prhs[0]);
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            lattice[i][j] = lattice_ptr[i + 3 * j];
        }
    }

    int num_atom = mxGetM(prhs[1]);
    mexutil::Buffer2D<double, 3> position(num_atom);
    if (mxGetN(prhs[1]) != 3) {
        mexErrMsgIdAndTxt("Spglib:invalidPosition",
                          "Position array must be Nx3.");
    }
    double *position_ptr = mxGetPr(prhs[1]);
    for (int i = 0; i < num_atom; i++) {
        for (int j = 0; j < 3; j++) {
            position[i][j] = position_ptr[i + num_atom * j];
        }
    }

    if (mxGetNumberOfElements(prhs[2]) != num_atom) {
        mexErrMsgIdAndTxt("Spglib:invalidTypes",
                          "Types array size must match the number of atoms.");
    }
    mexutil::Buffer1D<int> types(num_atom);
    double *types_ptr = mxGetPr(prhs[2]);
    for (int i = 0; i < num_atom; i++) {
        types[i] = static_cast<int>(types_ptr[i]);
    }

    if (!mxIsDouble(prhs[4]) || mxGetNumberOfElements(prhs[4]) != 1) {
        mexErrMsgIdAndTxt("Spglib:invalidSymprec", "Symprec must be a scalar.");
    }
    double symprec = mxGetScalar(prhs[4]);

    // Call spg_get_dataset
    mexutil::DatasetPtr dataset(
        spg_get_dataset(lattice, position, types, num_atom, symprec));
    if (!dataset) {
        throwLastSpglibError("spg_get_dataset failed");
    }
    plhs[0] = mexutil::makeDatasetStruct(*dataset);
}

// dataset = symspg('spg_get_magnetic_dataset', lattice, position, types,
// tensors, tensor_rank, num_atom, is_axial, symprec)
void SpglibFunctions::spg_get_magnetic_dataset_mex(int nlhs, mxArray *plhs[],
                                                   int nrhs,
                                                   mxArray const *prhs[]) {
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

    // Validate the number of input arguments
    int const expected_number_of_inputs = 8;
    if (nrhs != expected_number_of_inputs) {
        mexErrMsgIdAndTxt(
            "Spglib:invalidNumInputs",
            "Incorrect number of inputs for spg_get_magnetic_dataset.");
    }

    // Extract and validate the lattice argument
    double lattice[3][3];
    if (mxGetM(prhs[0]) != 3 || mxGetN(prhs[0]) != 3) {
        mexErrMsgIdAndTxt("Spglib:invalidLattice",
                          "Lattice matrix must be 3x3.");
    }

    double *lattice_ptr = mxGetPr(prhs[0]);
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            lattice[i][j] = lattice_ptr[i + 3 * j];
        }
    }

    // Extract and validate the position argument
    int num_atom = mxGetM(prhs[1]);
    mexutil::Buffer2D<double, 3> position(num_atom);
    if (mxGetN(prhs[1]) != 3) {
        mexErrMsgIdAndTxt("Spglib:invalidPosition",
                          "Position array must be Nx3.");
    }
    double *position_ptr = mxGetPr(prhs[1]);
    for (int i = 0; i < num_atom; i++) {
        for (int j = 0; j < 3; j++) {
            position[i][j] = position_ptr[i + num_atom * j];
        }
    }

    // Extract and validate the types argument
    if (mxGetNumberOfElements(prhs[2]) != num_atom) {
        mexErrMsgIdAndTxt("Spglib:invalidTypes",
                          "Types array size must match the number of atoms.");
    }
    mexutil::Buffer1D<int> types(num_atom);
    double *types_ptr = mxGetPr(prhs[2]);
    for (int i = 0; i < num_atom; i++) {
        types[i] = static_cast<int>(types_ptr[i]);
    }

    // Extract and validate the tensors argument
    double *tensors = mxGetPr(prhs[3]);

    // Extract and validate the tensor_rank argument
    int tensor_rank = static_cast<int>(mxGetScalar(prhs[4]));

    // Extract and validate the is_axial argument
    int is_axial = static_cast<int>(mxGetScalar(prhs[6]));

    // Extract and validate the symprec argument
    if (!mxIsDouble(prhs[7]) || mxGetNumberOfElements(prhs[7]) != 1) {
        mexErrMsgIdAndTxt("Spglib:invalidSymprec", "Symprec must be a scalar.");
    }
    double symprec = mxGetScalar(prhs[7]);

    // Call spg_get_magnetic_dataset
    mexutil::MagneticDatasetPtr dataset(
        spg_get_magnetic_dataset(lattice, position, types, tensors, tensor_rank,
                                 num_atom, is_axial, symprec));
    if (!dataset) {
        throwLastSpglibError("spg_get_magnetic_dataset failed");
    }
    plhs[0] = mexutil::makeMagneticDatasetStruct(*dataset);
}

// dataset = symspg('spgms_get_magnetic_dataset', lattice, position, types,
// tensors, tensor_rank, num_atom, is_axial, symprec, angle_tolerance,
// mag_symprec)
void SpglibFunctions::spgms_get_magnetic_dataset_mex(int nlhs, mxArray *plhs[],
                                                     int nrhs,
                                                     mxArray const *prhs[]) {
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

    // Validate the number of input arguments
    int const expected_number_of_inputs = 10;
    if (nrhs != expected_number_of_inputs) {
        mexErrMsgIdAndTxt(
            "Spglib:invalidNumInputs",
            "Incorrect number of inputs for spgms_get_magnetic_dataset.");
    }

    // Extract and validate the lattice argument
    double lattice[3][3];
    if (mxGetM(prhs[0]) != 3 || mxGetN(prhs[0]) != 3) {
        mexErrMsgIdAndTxt("Spglib:invalidLattice",
                          "Lattice matrix must be 3x3.");
    }

    double *lattice_ptr = mxGetPr(prhs[0]);
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            lattice[i][j] = lattice_ptr[i + 3 * j];
        }
    }

    // Extract and validate the position argument
    int num_atom = mxGetM(prhs[1]);
    mexutil::Buffer2D<double, 3> position(num_atom);
    if (mxGetN(prhs[1]) != 3) {
        mexErrMsgIdAndTxt("Spglib:invalidPosition",
                          "Position array must be Nx3.");
    }
    double *position_ptr = mxGetPr(prhs[1]);
    for (int i = 0; i < num_atom; i++) {
        for (int j = 0; j < 3; j++) {
            position[i][j] = position_ptr[i + num_atom * j];
        }
    }

    // Extract and validate the types argument
    if (mxGetNumberOfElements(prhs[2]) != num_atom) {
        mexErrMsgIdAndTxt("Spglib:invalidTypes",
                          "Types array size must match the number of atoms.");
    }
    mexutil::Buffer1D<int> types(num_atom);
    double *types_ptr = mxGetPr(prhs[2]);
    for (int i = 0; i < num_atom; i++) {
        types[i] = static_cast<int>(types_ptr[i]);
    }

    // Extract and validate the tensors argument
    double *tensors = mxGetPr(prhs[3]);

    // Extract and validate the tensor_rank argument
    int tensor_rank = static_cast<int>(mxGetScalar(prhs[4]));

    // Extract and validate the is_axial argument
    int is_axial = static_cast<int>(mxGetScalar(prhs[6]));

    // Extract and validate the symprec argument
    if (!mxIsDouble(prhs[7]) || mxGetNumberOfElements(prhs[7]) != 1) {
        mexErrMsgIdAndTxt("Spglib:invalidSymprec", "Symprec must be a scalar.");
    }
    double symprec = mxGetScalar(prhs[7]);

    // Extract and validate the angle_tolerance argument
    if (!mxIsDouble(prhs[8]) || mxGetNumberOfElements(prhs[8]) != 1) {
        mexErrMsgIdAndTxt("Spglib:invalidAngleTolerance",
                          "Angle tolerance must be a scalar.");
    }
    double angle_tolerance = mxGetScalar(prhs[8]);

    // Extract and validate the mag_symprec argument
    if (!mxIsDouble(prhs[9]) || mxGetNumberOfElements(prhs[9]) != 1) {
        mexErrMsgIdAndTxt("Spglib:invalidMagSymprec",
                          "Magnetic symmetry precision must be a scalar.");
    }
    double mag_symprec = mxGetScalar(prhs[9]);

    // Call spgms_get_magnetic_dataset
    mexutil::MagneticDatasetPtr dataset(spgms_get_magnetic_dataset(
        lattice, position, types, tensors, tensor_rank, num_atom, is_axial,
        symprec, angle_tolerance, mag_symprec));
    if (!dataset) {
        throwLastSpglibError("spgms_get_magnetic_dataset failed");
    }
    plhs[0] = mexutil::makeMagneticDatasetStruct(*dataset);
}

// dataset = symspg('spgat_get_dataset', lattice, position, types, num_atom,
// symprec, angle_tolerance)
void SpglibFunctions::spgat_get_dataset_mex(int nlhs, mxArray *plhs[], int nrhs,
                                            mxArray const *prhs[]) {
    /*
     SpglibDataset *spgat_get_dataset(
         double const lattice[3][3],
         double const position[][3],
         int const types[],
         int const num_atom,
         double const symprec,
         double const angle_tolerance);
    */

    // Validate the number of input arguments
    int const expected_number_of_inputs = 6;
    if (nrhs != expected_number_of_inputs) {
        mexErrMsgIdAndTxt("Spglib:invalidNumInputs",
                          "Incorrect number of inputs for spgat_get_dataset.");
    }

    // Extract and validate the lattice argument
    double lattice[3][3];
    if (mxGetM(prhs[0]) != 3 || mxGetN(prhs[0]) != 3) {
        mexErrMsgIdAndTxt("Spglib:invalidLattice",
                          "Lattice matrix must be 3x3.");
    }

    double *lattice_ptr = mxGetPr(prhs[0]);
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            lattice[i][j] = lattice_ptr[i + 3 * j];
        }
    }

    // Extract and validate the position argument
    int num_atom = mxGetM(prhs[1]);
    mexutil::Buffer2D<double, 3> position(num_atom);
    if (mxGetN(prhs[1]) != 3) {
        mexErrMsgIdAndTxt("Spglib:invalidPosition",
                          "Position array must be Nx3.");
    }
    double *position_ptr = mxGetPr(prhs[1]);
    for (int i = 0; i < num_atom; i++) {
        for (int j = 0; j < 3; j++) {
            position[i][j] = position_ptr[i + num_atom * j];
        }
    }

    // Extract and validate the types argument
    if (mxGetNumberOfElements(prhs[2]) != num_atom) {
        mexErrMsgIdAndTxt("Spglib:invalidTypes",
                          "Types array size must match the number of atoms.");
    }
    mexutil::Buffer1D<int> types(num_atom);
    double *types_ptr = mxGetPr(prhs[2]);
    for (int i = 0; i < num_atom; i++) {
        types[i] = static_cast<int>(types_ptr[i]);
    }

    // Extract and validate the symprec argument
    if (!mxIsDouble(prhs[4]) || mxGetNumberOfElements(prhs[4]) != 1) {
        mexErrMsgIdAndTxt("Spglib:invalidSymprec", "Symprec must be a scalar.");
    }
    double symprec = mxGetScalar(prhs[4]);

    // Extract and validate the angle_tolerance argument
    if (!mxIsDouble(prhs[5]) || mxGetNumberOfElements(prhs[5]) != 1) {
        mexErrMsgIdAndTxt("Spglib:invalidAngleTolerance",
                          "Angle tolerance must be a scalar.");
    }
    double angle_tolerance = mxGetScalar(prhs[5]);

    // Call spgat_get_dataset
    mexutil::DatasetPtr dataset(spgat_get_dataset(
        lattice, position, types, num_atom, symprec, angle_tolerance));
    if (!dataset) {
        throwLastSpglibError("spgat_get_dataset failed");
    }
    plhs[0] = mexutil::makeDatasetStruct(*dataset);
}

// dataset = symspg('spg_get_dataset_with_hall_number', lattice, position,
// types, num_atom, hall_number symprec)
void SpglibFunctions::spg_get_dataset_with_hall_number_mex(
    int nlhs, mxArray *plhs[], int nrhs, mxArray const *prhs[]) {
    /*
     SpglibDataset * spg_get_dataset_with_hall_number(const double
     lattice[3][3], const double position[][3], const int types[], const int
     num_atom, const int hall_number, const double symprec)
     */
    // Validate the number of input arguments
    int const expected_number_of_inputs = 6;
    if (nrhs != expected_number_of_inputs) {
        mexErrMsgIdAndTxt(
            "Spglib:invalidNumInputs",
            "Incorrect number of inputs for spg_get_dataset_with_hall_number.");
    }

    // Extract and validate the lattice, position, types, hall_number, and
    // symprec arguments
    double lattice[3][3];
    if (mxGetM(prhs[0]) != 3 || mxGetN(prhs[0]) != 3) {
        mexErrMsgIdAndTxt("Spglib:invalidLattice",
                          "Lattice matrix must be 3x3.");
    }

    double *lattice_ptr = mxGetPr(prhs[0]);
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            lattice[i][j] = lattice_ptr[i + 3 * j];
        }
    }

    int num_atom = mxGetM(prhs[1]);
    mexutil::Buffer2D<double, 3> position(num_atom);
    if (mxGetN(prhs[1]) != 3) {
        mexErrMsgIdAndTxt("Spglib:invalidPosition",
                          "Position array must be Nx3.");
    }
    double *position_ptr = mxGetPr(prhs[1]);
    for (int i = 0; i < num_atom; i++) {
        for (int j = 0; j < 3; j++) {
            position[i][j] = position_ptr[i + num_atom * j];
        }
    }

    if (mxGetNumberOfElements(prhs[2]) != num_atom) {
        mexErrMsgIdAndTxt("Spglib:invalidTypes",
                          "Types array size must match the number of atoms.");
    }
    mexutil::Buffer1D<int> types(num_atom);
    double *types_ptr = mxGetPr(prhs[2]);
    for (int i = 0; i < num_atom; i++) {
        types[i] = static_cast<int>(types_ptr[i]);
    }

    int hall_number = static_cast<int>(mxGetScalar(prhs[4]));

    if (!mxIsDouble(prhs[5]) || mxGetNumberOfElements(prhs[5]) != 1) {
        mexErrMsgIdAndTxt("Spglib:invalidSymprec", "Symprec must be a scalar.");
    }
    double symprec = mxGetScalar(prhs[5]);

    // Call spg_get_dataset
    mexutil::DatasetPtr dataset(spg_get_dataset_with_hall_number(
        lattice, position, types, num_atom, hall_number, symprec));
    if (!dataset) {
        throwLastSpglibError("spg_get_dataset_with_hall_number failed");
    }
    plhs[0] = mexutil::makeDatasetStruct(*dataset);
}

// dataset = symspg('spgat_get_dataset_with_hall_number', lattice, position,
// types, num_atom, hall_number, symprec, angle_tolerance)
void SpglibFunctions::spgat_get_dataset_with_hall_number_mex(
    int nlhs, mxArray *plhs[], int nrhs, mxArray const *prhs[]) {
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

    // Validate the number of input arguments
    int const expected_number_of_inputs = 7;
    if (nrhs != expected_number_of_inputs) {
        mexErrMsgIdAndTxt("Spglib:invalidNumInputs",
                          "Incorrect number of inputs for "
                          "spgat_get_dataset_with_hall_number.");
    }

    // Extract and validate the lattice argument
    double lattice[3][3];
    if (mxGetM(prhs[0]) != 3 || mxGetN(prhs[0]) != 3) {
        mexErrMsgIdAndTxt("Spglib:invalidLattice",
                          "Lattice matrix must be 3x3.");
    }

    double *lattice_ptr = mxGetPr(prhs[0]);
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            lattice[i][j] = lattice_ptr[i + 3 * j];
        }
    }

    // Extract and validate the position argument
    int num_atom = mxGetM(prhs[1]);
    mexutil::Buffer2D<double, 3> position(num_atom);
    if (mxGetN(prhs[1]) != 3) {
        mexErrMsgIdAndTxt("Spglib:invalidPosition",
                          "Position array must be Nx3.");
    }
    double *position_ptr = mxGetPr(prhs[1]);
    for (int i = 0; i < num_atom; i++) {
        for (int j = 0; j < 3; j++) {
            position[i][j] = position_ptr[i + num_atom * j];
        }
    }

    // Extract and validate the types argument
    if (mxGetNumberOfElements(prhs[2]) != num_atom) {
        mexErrMsgIdAndTxt("Spglib:invalidTypes",
                          "Types array size must match the number of atoms.");
    }
    mexutil::Buffer1D<int> types(num_atom);
    double *types_ptr = mxGetPr(prhs[2]);
    for (int i = 0; i < num_atom; i++) {
        types[i] = static_cast<int>(types_ptr[i]);
    }

    // Extract and validate the hall_number argument
    int hall_number = static_cast<int>(mxGetScalar(prhs[4]));

    // Extract and validate the symprec argument
    if (!mxIsDouble(prhs[5]) || mxGetNumberOfElements(prhs[5]) != 1) {
        mexErrMsgIdAndTxt("Spglib:invalidSymprec", "Symprec must be a scalar.");
    }
    double symprec = mxGetScalar(prhs[5]);

    // Extract and validate the angle_tolerance argument
    if (!mxIsDouble(prhs[6]) || mxGetNumberOfElements(prhs[6]) != 1) {
        mexErrMsgIdAndTxt("Spglib:invalidAngleTolerance",
                          "Angle tolerance must be a scalar.");
    }
    double angle_tolerance = mxGetScalar(prhs[6]);

    // Call spgat_get_dataset_with_hall_number
    mexutil::DatasetPtr dataset(spgat_get_dataset_with_hall_number(
        lattice, position, types, num_atom, hall_number, symprec,
        angle_tolerance));
    if (!dataset) {
        throwLastSpglibError("spgat_get_dataset_with_hall_number failed");
    }
    plhs[0] = mexutil::makeDatasetStruct(*dataset);
}

// [rotations, translations, equivalent_atoms, n_operations] =
// symspg('spg_get_symmetry_with_collinear_spin', max_size, lattice, position,
// types, spins, num_atom, symprec)
void SpglibFunctions::spg_get_symmetry_with_collinear_spin_mex(
    int nlhs, mxArray *plhs[], int nrhs, mxArray const *prhs[]) {
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

    // Validate the number of input arguments
    int const expected_number_of_inputs = 7;
    if (nrhs != expected_number_of_inputs) {
        mexErrMsgIdAndTxt("Spglib:invalidNumInputs",
                          "Incorrect number of inputs for "
                          "spg_get_symmetry_with_collinear_spin.");
    }

    // Extract and validate the max_size argument
    int max_size = static_cast<int>(mxGetScalar(prhs[0]));

    // Extract and validate the lattice argument
    double lattice[3][3];
    if (mxGetM(prhs[1]) != 3 || mxGetN(prhs[1]) != 3) {
        mexErrMsgIdAndTxt("Spglib:invalidLattice",
                          "Lattice matrix must be 3x3.");
    }
    double *lattice_ptr = mxGetPr(prhs[1]);
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            lattice[i][j] = lattice_ptr[i + 3 * j];
        }
    }

    // Extract and validate the position argument
    int num_atom = mxGetM(prhs[2]);
    mexutil::Buffer2D<double, 3> position(num_atom);
    if (mxGetN(prhs[2]) != 3) {
        mexErrMsgIdAndTxt("Spglib:invalidPosition",
                          "Position array must be Nx3.");
    }
    double *position_ptr = mxGetPr(prhs[2]);
    for (int i = 0; i < num_atom; i++) {
        for (int j = 0; j < 3; j++) {
            position[i][j] = position_ptr[i + num_atom * j];
        }
    }

    // Extract and validate the types argument
    if (mxGetNumberOfElements(prhs[3]) != num_atom) {
        mexErrMsgIdAndTxt("Spglib:invalidTypes",
                          "Types array size must match the number of atoms.");
    }
    mexutil::Buffer1D<int> types(num_atom);
    double *types_ptr = mxGetPr(prhs[3]);
    for (int i = 0; i < num_atom; i++) {
        types[i] = static_cast<int>(types_ptr[i]);
    }

    // Extract and validate the spins argument
    if (mxGetNumberOfElements(prhs[4]) != num_atom) {
        mexErrMsgIdAndTxt("Spglib:invalidSpins",
                          "Spins array size must match the number of atoms.");
    }
    mexutil::Buffer1D<double> spins(num_atom);
    double *spins_ptr = mxGetPr(prhs[4]);
    for (int i = 0; i < num_atom; i++) {
        spins[i] = spins_ptr[i];
    }

    // Extract and validate the symprec argument
    if (!mxIsDouble(prhs[6]) || mxGetNumberOfElements(prhs[6]) != 1) {
        mexErrMsgIdAndTxt("Spglib:invalidSymprec", "Symprec must be a scalar.");
    }
    double symprec = mxGetScalar(prhs[6]);

    // Initialize the rotation, translation, and equivalent_atoms arrays
    mexutil::Buffer3D<int, 3, 3> rotation(max_size);
    mexutil::Buffer2D<double, 3> translation(max_size);
    mexutil::Buffer1D<int> equivalent_atoms(num_atom);

    // Call spg_get_symmetry_with_collinear_spin
    int n_operations = spg_get_symmetry_with_collinear_spin(
        rotation, translation, equivalent_atoms, max_size, lattice, position,
        types, spins, num_atom, symprec);

    if (n_operations == 0) {
        throwLastSpglibError("spg_get_symmetry_with_collinear_spin failed");
    }

    // Create and populate the output arrays
    // Output the rotation matrices (Nx3x3 int array)
    mwSize dims[3] = {static_cast<mwSize>(n_operations), 3, 3};
    plhs[0] = mxCreateNumericArray(3, dims, mxINT32_CLASS, mxREAL);
    auto *rotations_out = static_cast<int32_t *>(mxGetData(plhs[0]));

    for (int k = 0; k < n_operations; ++k) {
        for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < 3; ++j) {
                rotations_out[k + i * n_operations + j * n_operations * 3] =
                    rotation[k][i][j];
            }
        }
    }

    // Output the translations (Nx3 double array)
    plhs[1] = mxCreateDoubleMatrix(n_operations, 3, mxREAL);
    double *translations_out = mxGetPr(plhs[1]);
    for (int i = 0; i < n_operations; ++i) {
        for (int j = 0; j < 3; ++j) {
            translations_out[i + j * n_operations] = translation[i][j];
        }
    }

    // Output the equivalent atoms (num_atom int array)
    plhs[2] = mxCreateNumericMatrix(num_atom, 1, mxINT32_CLASS, mxREAL);
    int32_t *equivalent_atoms_out = static_cast<int32_t *>(mxGetData(plhs[2]));
    for (int i = 0; i < num_atom; ++i) {
        equivalent_atoms_out[i] = equivalent_atoms[i];
    }

    // Output the number of operations
    plhs[3] = mxCreateDoubleScalar(static_cast<double>(n_operations));
}

// [rotations, translations, equivalent_atoms, n_operations] =
// symspg('spgat_get_symmetry_with_collinear_spin', max_size, lattice, position,
// types, spins, num_atom, symprec, angle_tolerance)
void SpglibFunctions::spgat_get_symmetry_with_collinear_spin_mex(
    int nlhs, mxArray *plhs[], int nrhs, mxArray const *prhs[]) {
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

    // Validate the number of input arguments
    int const expected_number_of_inputs = 8;
    if (nrhs != expected_number_of_inputs) {
        mexErrMsgIdAndTxt("Spglib:invalidNumInputs",
                          "Incorrect number of inputs for "
                          "spgat_get_symmetry_with_collinear_spin.");
    }

    // Extract and validate the max_size argument
    int max_size = static_cast<int>(mxGetScalar(prhs[0]));

    // Extract and validate the lattice argument
    double lattice[3][3];
    if (mxGetM(prhs[1]) != 3 || mxGetN(prhs[1]) != 3) {
        mexErrMsgIdAndTxt("Spglib:invalidLattice",
                          "Lattice matrix must be 3x3.");
    }
    double *lattice_ptr = mxGetPr(prhs[1]);
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            lattice[i][j] = lattice_ptr[i + 3 * j];
        }
    }

    // Extract and validate the position argument
    int num_atom = mxGetM(prhs[2]);
    mexutil::Buffer2D<double, 3> position(num_atom);
    if (mxGetN(prhs[2]) != 3) {
        mexErrMsgIdAndTxt("Spglib:invalidPosition",
                          "Position array must be Nx3.");
    }
    double *position_ptr = mxGetPr(prhs[2]);
    for (int i = 0; i < num_atom; i++) {
        for (int j = 0; j < 3; j++) {
            position[i][j] = position_ptr[i + num_atom * j];
        }
    }

    // Extract and validate the types argument
    if (mxGetNumberOfElements(prhs[3]) != num_atom) {
        mexErrMsgIdAndTxt("Spglib:invalidTypes",
                          "Types array size must match the number of atoms.");
    }
    mexutil::Buffer1D<int> types(num_atom);
    double *types_ptr = mxGetPr(prhs[3]);
    for (int i = 0; i < num_atom; i++) {
        types[i] = static_cast<int>(types_ptr[i]);
    }

    // Extract and validate the spins argument
    if (mxGetNumberOfElements(prhs[4]) != num_atom) {
        mexErrMsgIdAndTxt("Spglib:invalidSpins",
                          "Spins array size must match the number of atoms.");
    }
    mexutil::Buffer1D<double> spins(num_atom);
    double *spins_ptr = mxGetPr(prhs[4]);
    for (int i = 0; i < num_atom; i++) {
        spins[i] = spins_ptr[i];
    }

    // Extract and validate the symprec argument
    if (!mxIsDouble(prhs[6]) || mxGetNumberOfElements(prhs[6]) != 1) {
        mexErrMsgIdAndTxt("Spglib:invalidSymprec", "Symprec must be a scalar.");
    }
    double symprec = mxGetScalar(prhs[6]);

    // Extract and validate the angle_tolerance argument
    if (!mxIsDouble(prhs[7]) || mxGetNumberOfElements(prhs[7]) != 1) {
        mexErrMsgIdAndTxt("Spglib:invalidAngleTolerance",
                          "Angle tolerance must be a scalar.");
    }
    double angle_tolerance = mxGetScalar(prhs[7]);

    // Initialize the rotation, translation, and equivalent_atoms arrays
    mexutil::Buffer3D<int, 3, 3> rotation(max_size);
    mexutil::Buffer2D<double, 3> translation(max_size);
    mexutil::Buffer1D<int> equivalent_atoms(num_atom);

    // Call spgat_get_symmetry_with_collinear_spin
    int n_operations = spgat_get_symmetry_with_collinear_spin(
        rotation, translation, equivalent_atoms, max_size, lattice, position,
        types, spins, num_atom, symprec, angle_tolerance);

    if (n_operations == 0) {
        throwLastSpglibError("spgat_get_symmetry_with_collinear_spin failed");
    }

    // Create and populate the output arrays
    // Output the rotation matrices (Nx3x3 int array)
    mwSize dims[3] = {static_cast<mwSize>(n_operations), 3, 3};
    plhs[0] = mxCreateNumericArray(3, dims, mxINT32_CLASS, mxREAL);
    auto *rotations_out = static_cast<int32_t *>(mxGetData(plhs[0]));

    for (int k = 0; k < n_operations; ++k) {
        for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < 3; ++j) {
                rotations_out[k + i * n_operations + j * n_operations * 3] =
                    rotation[k][i][j];
            }
        }
    }

    // Output the translations (Nx3 double array)
    plhs[1] = mxCreateDoubleMatrix(n_operations, 3, mxREAL);
    double *translations_out = mxGetPr(plhs[1]);
    for (int i = 0; i < n_operations; ++i) {
        for (int j = 0; j < 3; ++j) {
            translations_out[i + j * n_operations] = translation[i][j];
        }
    }

    // Output the equivalent atoms (num_atom int array)
    plhs[2] = mxCreateNumericMatrix(num_atom, 1, mxINT32_CLASS, mxREAL);
    int32_t *equivalent_atoms_out = static_cast<int32_t *>(mxGetData(plhs[2]));
    for (int i = 0; i < num_atom; ++i) {
        equivalent_atoms_out[i] = equivalent_atoms[i];
    }

    // Output the number of operations
    plhs[3] = mxCreateDoubleScalar(static_cast<double>(n_operations));
}

// [rotations, translations, equivalent_atoms, n_operations] =
// symspg('spgms_get_symmetry_with_collinear_spin', max_size, lattice, position,
// types, spins, num_atom, symprec, angle_tolerance, mag_symprec)
void SpglibFunctions::spgms_get_symmetry_with_collinear_spin_mex(
    int nlhs, mxArray *plhs[], int nrhs, mxArray const *prhs[]) {
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

    // Validate the number of input arguments
    int const expected_number_of_inputs = 9;
    if (nrhs != expected_number_of_inputs) {
        mexErrMsgIdAndTxt("Spglib:invalidNumInputs",
                          "Incorrect number of inputs for "
                          "spgms_get_symmetry_with_collinear_spin.");
    }

    // Extract and validate the max_size argument
    int max_size = static_cast<int>(mxGetScalar(prhs[0]));

    // Extract and validate the lattice argument
    double lattice[3][3];
    if (mxGetM(prhs[1]) != 3 || mxGetN(prhs[1]) != 3) {
        mexErrMsgIdAndTxt("Spglib:invalidLattice",
                          "Lattice matrix must be 3x3.");
    }
    double *lattice_ptr = mxGetPr(prhs[1]);
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            lattice[i][j] = lattice_ptr[i + 3 * j];
        }
    }

    // Extract and validate the position argument
    int num_atom = mxGetM(prhs[2]);
    mexutil::Buffer2D<double, 3> position(num_atom);
    if (mxGetN(prhs[2]) != 3) {
        mexErrMsgIdAndTxt("Spglib:invalidPosition",
                          "Position array must be Nx3.");
    }
    double *position_ptr = mxGetPr(prhs[2]);
    for (int i = 0; i < num_atom; i++) {
        for (int j = 0; j < 3; j++) {
            position[i][j] = position_ptr[i + num_atom * j];
        }
    }

    // Extract and validate the types argument
    if (mxGetNumberOfElements(prhs[3]) != num_atom) {
        mexErrMsgIdAndTxt("Spglib:invalidTypes",
                          "Types array size must match the number of atoms.");
    }
    mexutil::Buffer1D<int> types(num_atom);
    double *types_ptr = mxGetPr(prhs[3]);
    for (int i = 0; i < num_atom; i++) {
        types[i] = static_cast<int>(types_ptr[i]);
    }

    // Extract and validate the spins argument
    if (mxGetNumberOfElements(prhs[4]) != num_atom) {
        mexErrMsgIdAndTxt("Spglib:invalidSpins",
                          "Spins array size must match the number of atoms.");
    }
    mexutil::Buffer1D<double> spins(num_atom);
    double *spins_ptr = mxGetPr(prhs[4]);
    for (int i = 0; i < num_atom; i++) {
        spins[i] = spins_ptr[i];
    }

    // Extract and validate the symprec argument
    if (!mxIsDouble(prhs[6]) || mxGetNumberOfElements(prhs[6]) != 1) {
        mexErrMsgIdAndTxt("Spglib:invalidSymprec", "Symprec must be a scalar.");
    }
    double symprec = mxGetScalar(prhs[6]);

    // Extract and validate the angle_tolerance argument
    if (!mxIsDouble(prhs[7]) || mxGetNumberOfElements(prhs[7]) != 1) {
        mexErrMsgIdAndTxt("Spglib:invalidAngleTolerance",
                          "Angle tolerance must be a scalar.");
    }
    double angle_tolerance = mxGetScalar(prhs[7]);

    // Extract and validate the mag_symprec argument
    if (!mxIsDouble(prhs[8]) || mxGetNumberOfElements(prhs[8]) != 1) {
        mexErrMsgIdAndTxt("Spglib:invalidMagSymprec",
                          "Magnetic symmetry precision must be a scalar.");
    }
    double mag_symprec = mxGetScalar(prhs[8]);

    // Initialize the rotation, translation, and equivalent_atoms arrays
    mexutil::Buffer3D<int, 3, 3> rotation(max_size);
    mexutil::Buffer2D<double, 3> translation(max_size);
    mexutil::Buffer1D<int> equivalent_atoms(num_atom);

    // Call spgms_get_symmetry_with_collinear_spin
    int n_operations = spgms_get_symmetry_with_collinear_spin(
        rotation, translation, equivalent_atoms, max_size, lattice, position,
        types, spins, num_atom, symprec, angle_tolerance, mag_symprec);

    if (n_operations == 0) {
        throwLastSpglibError("spgms_get_symmetry_with_collinear_spin failed");
    }

    // Create and populate the output arrays
    // Output the rotation matrices (Nx3x3 int array)
    mwSize dims[3] = {static_cast<mwSize>(n_operations), 3, 3};
    plhs[0] = mxCreateNumericArray(3, dims, mxINT32_CLASS, mxREAL);
    auto *rotations_out = static_cast<int32_t *>(mxGetData(plhs[0]));

    for (int k = 0; k < n_operations; ++k) {
        for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < 3; ++j) {
                rotations_out[k + i * n_operations + j * n_operations * 3] =
                    rotation[k][i][j];
            }
        }
    }

    // Output the translations (Nx3 double array)
    plhs[1] = mxCreateDoubleMatrix(n_operations, 3, mxREAL);
    double *translations_out = mxGetPr(plhs[1]);
    for (int i = 0; i < n_operations; ++i) {
        for (int j = 0; j < 3; ++j) {
            translations_out[i + j * n_operations] = translation[i][j];
        }
    }

    // Output the equivalent atoms (num_atom int array)
    plhs[2] = mxCreateNumericMatrix(num_atom, 1, mxINT32_CLASS, mxREAL);
    int32_t *equivalent_atoms_out = static_cast<int32_t *>(mxGetData(plhs[2]));
    for (int i = 0; i < num_atom; ++i) {
        equivalent_atoms_out[i] = equivalent_atoms[i];
    }

    // Output the number of operations
    plhs[3] = mxCreateDoubleScalar(static_cast<double>(n_operations));
}

// [rotations, translations, equivalent_atoms, primitive_lattice, spin_flips,
// n_operations] = symspg('spg_get_symmetry_with_site_tensors', max_size,
// lattice, position, types, tensors, tensor_rank, num_atom, with_time_reversal,
// is_axial, symprec)
void SpglibFunctions::spg_get_symmetry_with_site_tensors_mex(
    int nlhs, mxArray *plhs[], int nrhs, mxArray const *prhs[]) {
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

    // Validate the number of input arguments
    int const expected_number_of_inputs =
        10;  // Updated expected number of inputs
    if (nrhs != expected_number_of_inputs) {
        mexErrMsgIdAndTxt("Spglib:invalidNumInputs",
                          "Incorrect number of inputs for "
                          "spg_get_symmetry_with_site_tensors.");
    }

    // Extract and validate the max_size argument
    int max_size = static_cast<int>(mxGetScalar(prhs[0]));

    // Extract and validate the lattice argument
    double lattice[3][3];
    if (mxGetM(prhs[1]) != 3 || mxGetN(prhs[1]) != 3) {
        mexErrMsgIdAndTxt("Spglib:invalidLattice",
                          "Lattice matrix must be 3x3.");
    }
    double *lattice_ptr = mxGetPr(prhs[1]);
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            lattice[i][j] = lattice_ptr[i + 3 * j];
        }
    }

    // Extract and validate the position argument
    int num_atom = mxGetM(prhs[2]);
    mexutil::Buffer2D<double, 3> position(num_atom);
    if (mxGetN(prhs[2]) != 3) {
        mexErrMsgIdAndTxt("Spglib:invalidPosition",
                          "Position array must be Nx3.");
    }
    double *position_ptr = mxGetPr(prhs[2]);
    for (int i = 0; i < num_atom; i++) {
        for (int j = 0; j < 3; j++) {
            position[i][j] = position_ptr[i + num_atom * j];
        }
    }

    // Extract and validate the types argument
    if (mxGetNumberOfElements(prhs[3]) != num_atom) {
        mexErrMsgIdAndTxt("Spglib:invalidTypes",
                          "Types array size must match the number of atoms.");
    }
    mexutil::Buffer1D<int> types(num_atom);
    double *types_ptr = mxGetPr(prhs[3]);
    for (int i = 0; i < num_atom; i++) {
        types[i] = static_cast<int>(types_ptr[i]);
    }

    // Extract and validate the tensors argument
    double *tensors = mxGetPr(prhs[4]);

    // Extract and validate the tensor_rank argument
    int tensor_rank = static_cast<int>(mxGetScalar(prhs[5]));

    // Extract and validate the with_time_reversal argument
    int with_time_reversal = static_cast<int>(mxGetScalar(prhs[7]));

    // Extract and validate the is_axial argument
    int is_axial = static_cast<int>(mxGetScalar(prhs[8]));

    // Extract and validate the symprec argument
    if (!mxIsDouble(prhs[9]) || mxGetNumberOfElements(prhs[9]) != 1) {
        mexErrMsgIdAndTxt("Spglib:invalidSymprec", "Symprec must be a scalar.");
    }
    double symprec = mxGetScalar(prhs[9]);

    // Initialize the rotation, translation, equivalent_atoms, and
    // primitive_lattice arrays
    mexutil::Buffer3D<int, 3, 3> rotation(max_size);
    mexutil::Buffer2D<double, 3> translation(max_size);
    mexutil::Buffer1D<int> equivalent_atoms(num_atom);
    double primitive_lattice[3][3];
    mexutil::Buffer1D<int> spin_flips(max_size);

    // Call spg_get_symmetry_with_site_tensors
    int n_operations = spg_get_symmetry_with_site_tensors(
        rotation, translation, equivalent_atoms, primitive_lattice, spin_flips,
        max_size, lattice, position, types, tensors, tensor_rank, num_atom,
        with_time_reversal, is_axial, symprec);

    if (n_operations == 0) {
        throwLastSpglibError("spg_get_symmetry_with_site_tensors failed");
    }

    // Create and populate the output arrays
    // Output the rotation matrices (Nx3x3 int array)
    mwSize dims[3] = {static_cast<mwSize>(n_operations), 3, 3};
    plhs[0] = mxCreateNumericArray(3, dims, mxINT32_CLASS, mxREAL);
    auto *rotations_out = static_cast<int32_t *>(mxGetData(plhs[0]));

    for (int k = 0; k < n_operations; ++k) {
        for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < 3; ++j) {
                rotations_out[k + i * n_operations + j * n_operations * 3] =
                    rotation[k][i][j];
            }
        }
    }

    // Output the translations (Nx3 double array)
    plhs[1] = mxCreateDoubleMatrix(n_operations, 3, mxREAL);
    double *translations_out = mxGetPr(plhs[1]);
    for (int i = 0; i < n_operations; ++i) {
        for (int j = 0; j < 3; ++j) {
            translations_out[i + j * n_operations] = translation[i][j];
        }
    }

    // Output the equivalent atoms (num_atom int array)
    plhs[2] = mxCreateNumericMatrix(num_atom, 1, mxINT32_CLASS, mxREAL);
    auto *equivalent_atoms_out = static_cast<int32_t *>(mxGetData(plhs[2]));
    for (int i = 0; i < num_atom; ++i) {
        equivalent_atoms_out[i] = equivalent_atoms[i];
    }

    // Output primitive_lattice (3x3 double array)
    plhs[3] = mxCreateDoubleMatrix(3, 3, mxREAL);
    double *primitive_lattice_out = mxGetPr(plhs[3]);
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            primitive_lattice_out[i + j * 3] = primitive_lattice[i][j];
        }
    }

    // Output spin_flips (num_atom int array)
    plhs[4] = mxCreateNumericMatrix(n_operations, 1, mxINT32_CLASS, mxREAL);
    auto *spin_flips_out = static_cast<int32_t *>(mxGetData(plhs[4]));
    for (int i = 0; i < n_operations; ++i) {
        spin_flips_out[i] = static_cast<int32_t>(spin_flips[i]);
    }

    // Output the number of operations
    plhs[5] = mxCreateDoubleScalar(static_cast<double>(n_operations));
}

// [rotations, translations, equivalent_atoms, primitive_lattice, spin_flips,
// n_operations] = symspg('spgat_get_symmetry_with_site_tensors', max_size,
// lattice, position, types, tensors, tensor_rank, num_atom, with_time_reversal,
// is_axial, symprec, angle_tolerance)
void SpglibFunctions::spgat_get_symmetry_with_site_tensors_mex(
    int nlhs, mxArray *plhs[], int nrhs, mxArray const *prhs[]) {
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

    // Validate the number of input arguments
    int const expected_number_of_inputs =
        11;  // Adjusted to match inputs without spin_flips
    if (nrhs != expected_number_of_inputs) {
        mexErrMsgIdAndTxt("Spglib:invalidNumInputs",
                          "Incorrect number of inputs for "
                          "spgat_get_symmetry_with_site_tensors.");
    }

    // Extract and validate the max_size argument
    int max_size = static_cast<int>(mxGetScalar(prhs[0]));

    // Extract and validate the lattice argument
    double lattice[3][3];
    if (mxGetM(prhs[1]) != 3 || mxGetN(prhs[1]) != 3) {
        mexErrMsgIdAndTxt("Spglib:invalidLattice",
                          "Lattice matrix must be 3x3.");
    }
    double *lattice_ptr = mxGetPr(prhs[1]);
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            lattice[i][j] = lattice_ptr[i + 3 * j];
        }
    }

    // Extract and validate the position argument
    int num_atom = mxGetM(prhs[2]);
    mexutil::Buffer2D<double, 3> position(num_atom);
    if (mxGetN(prhs[2]) != 3) {
        mexErrMsgIdAndTxt("Spglib:invalidPosition",
                          "Position array must be Nx3.");
    }
    double *position_ptr = mxGetPr(prhs[2]);
    for (int i = 0; i < num_atom; i++) {
        for (int j = 0; j < 3; j++) {
            position[i][j] = position_ptr[i + num_atom * j];
        }
    }

    // Extract and validate the types argument
    if (mxGetNumberOfElements(prhs[3]) != num_atom) {
        mexErrMsgIdAndTxt("Spglib:invalidTypes",
                          "Types array size must match the number of atoms.");
    }
    mexutil::Buffer1D<int> types(num_atom);
    double *types_ptr = mxGetPr(prhs[3]);
    for (int i = 0; i < num_atom; i++) {
        types[i] = static_cast<int>(types_ptr[i]);
    }

    // Extract and validate the tensors argument
    double *tensors = mxGetPr(prhs[4]);

    // Extract and validate the tensor_rank argument
    int tensor_rank = static_cast<int>(mxGetScalar(prhs[5]));

    // Extract and validate the with_time_reversal argument
    int with_time_reversal = static_cast<int>(mxGetScalar(prhs[7]));

    // Extract and validate the is_axial argument
    int is_axial = static_cast<int>(mxGetScalar(prhs[8]));

    // Extract and validate the symprec argument
    if (!mxIsDouble(prhs[9]) || mxGetNumberOfElements(prhs[9]) != 1) {
        mexErrMsgIdAndTxt("Spglib:invalidSymprec", "Symprec must be a scalar.");
    }
    double symprec = mxGetScalar(prhs[9]);

    // Extract and validate the angle_tolerance argument
    if (!mxIsDouble(prhs[10]) || mxGetNumberOfElements(prhs[10]) != 1) {
        mexErrMsgIdAndTxt("Spglib:invalidAngleTolerance",
                          "Angle tolerance must be a scalar.");
    }
    double angle_tolerance = mxGetScalar(prhs[9]);

    // Initialize the rotation, translation, equivalent_atoms, and
    // primitive_lattice arrays
    mexutil::Buffer3D<int, 3, 3> rotation(max_size);
    mexutil::Buffer2D<double, 3> translation(max_size);
    mexutil::Buffer1D<int> equivalent_atoms(num_atom);
    double primitive_lattice[3][3];
    mexutil::Buffer1D<int> spin_flips(max_size);

    // Call spgat_get_symmetry_with_site_tensors
    int n_operations = spgat_get_symmetry_with_site_tensors(
        rotation, translation, equivalent_atoms, primitive_lattice, spin_flips,
        max_size, lattice, position, types, tensors, tensor_rank, num_atom,
        with_time_reversal, is_axial, symprec, angle_tolerance);

    if (n_operations == 0) {
        throwLastSpglibError("spgat_get_symmetry_with_site_tensors failed");
    }

    // Create and populate the output arrays
    // Output the rotation matrices (Nx3x3 int array)
    mwSize dims[3] = {static_cast<mwSize>(n_operations), 3, 3};
    plhs[0] = mxCreateNumericArray(3, dims, mxINT32_CLASS, mxREAL);
    auto *rotations_out = static_cast<int32_t *>(mxGetData(plhs[0]));

    for (int k = 0; k < n_operations; ++k) {
        for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < 3; ++j) {
                rotations_out[k + i * n_operations + j * n_operations * 3] =
                    rotation[k][i][j];
            }
        }
    }

    // Output the translations (Nx3 double array)
    plhs[1] = mxCreateDoubleMatrix(n_operations, 3, mxREAL);
    double *translations_out = mxGetPr(plhs[1]);
    for (int i = 0; i < n_operations; ++i) {
        for (int j = 0; j < 3; ++j) {
            translations_out[i + j * n_operations] = translation[i][j];
        }
    }

    // Output the equivalent atoms (num_atom int array)
    plhs[2] = mxCreateNumericMatrix(num_atom, 1, mxINT32_CLASS, mxREAL);
    auto *equivalent_atoms_out = static_cast<int32_t *>(mxGetData(plhs[2]));
    for (int i = 0; i < num_atom; ++i) {
        equivalent_atoms_out[i] = equivalent_atoms[i];
    }

    // Output primitive_lattice (3x3 double array)
    plhs[3] = mxCreateDoubleMatrix(3, 3, mxREAL);
    double *primitive_lattice_out = mxGetPr(plhs[3]);
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            primitive_lattice_out[i + j * 3] = primitive_lattice[i][j];
        }
    }

    // Output spin_flips (n_operations int array)
    plhs[4] = mxCreateNumericMatrix(n_operations, 1, mxINT32_CLASS, mxREAL);
    auto *spin_flips_out = static_cast<int32_t *>(mxGetData(plhs[4]));
    for (int i = 0; i < n_operations; ++i) {
        spin_flips_out[i] = spin_flips[i];
    }

    // Output the number of operations
    plhs[5] = mxCreateDoubleScalar(static_cast<double>(n_operations));
}

// [rotations, translations, equivalent_atoms, primitive_lattice, spin_flips,
// n_operations] = symspg('spgms_get_symmetry_with_site_tensors', max_size,
// lattice, position, types, tensors, tensor_rank, num_atom, with_time_reversal,
// is_axial, symprec, angle_tolerance, mag_symprec)
void SpglibFunctions::spgms_get_symmetry_with_site_tensors_mex(
    int nlhs, mxArray *plhs[], int nrhs, mxArray const *prhs[]) {
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

    // Validate the number of input arguments
    int const expected_number_of_inputs =
        12;  // Adjusted to match inputs without spin_flips
    if (nrhs != expected_number_of_inputs) {
        mexErrMsgIdAndTxt("Spglib:invalidNumInputs",
                          "Incorrect number of inputs for "
                          "spgms_get_symmetry_with_site_tensors.");
    }

    // Extract and validate the max_size argument
    int max_size = static_cast<int>(mxGetScalar(prhs[0]));

    // Extract and validate the lattice argument
    double lattice[3][3];
    if (mxGetM(prhs[1]) != 3 || mxGetN(prhs[1]) != 3) {
        mexErrMsgIdAndTxt("Spglib:invalidLattice",
                          "Lattice matrix must be 3x3.");
    }
    double *lattice_ptr = mxGetPr(prhs[1]);
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            lattice[i][j] = lattice_ptr[i + 3 * j];
        }
    }

    // Extract and validate the position argument
    int num_atom = mxGetM(prhs[2]);
    mexutil::Buffer2D<double, 3> position(num_atom);
    if (mxGetN(prhs[2]) != 3) {
        mexErrMsgIdAndTxt("Spglib:invalidPosition",
                          "Position array must be Nx3.");
    }
    double *position_ptr = mxGetPr(prhs[2]);
    for (int i = 0; i < num_atom; i++) {
        for (int j = 0; j < 3; j++) {
            position[i][j] = position_ptr[i + num_atom * j];
        }
    }

    // Extract and validate the types argument
    if (mxGetNumberOfElements(prhs[3]) != num_atom) {
        mexErrMsgIdAndTxt("Spglib:invalidTypes",
                          "Types array size must match the number of atoms.");
    }
    mexutil::Buffer1D<int> types(num_atom);
    double *types_ptr = mxGetPr(prhs[3]);
    for (int i = 0; i < num_atom; i++) {
        types[i] = static_cast<int>(types_ptr[i]);
    }

    // Extract and validate the tensors argument
    double *tensors = mxGetPr(prhs[4]);

    // Extract and validate the tensor_rank argument
    int tensor_rank = static_cast<int>(mxGetScalar(prhs[5]));

    // Extract and validate the with_time_reversal argument
    int with_time_reversal = static_cast<int>(mxGetScalar(prhs[7]));

    // Extract and validate the is_axial argument
    int is_axial = static_cast<int>(mxGetScalar(prhs[8]));

    // Extract and validate the symprec argument
    if (!mxIsDouble(prhs[9]) || mxGetNumberOfElements(prhs[9]) != 1) {
        mexErrMsgIdAndTxt("Spglib:invalidSymprec", "Symprec must be a scalar.");
    }
    double symprec = mxGetScalar(prhs[9]);

    // Extract and validate the angle_tolerance argument
    if (!mxIsDouble(prhs[10]) || mxGetNumberOfElements(prhs[10]) != 1) {
        mexErrMsgIdAndTxt("Spglib:invalidAngleTolerance",
                          "Angle tolerance must be a scalar.");
    }
    double angle_tolerance = mxGetScalar(prhs[10]);

    // Extract and validate the mag_symprec argument
    if (!mxIsDouble(prhs[11]) || mxGetNumberOfElements(prhs[11]) != 1) {
        mexErrMsgIdAndTxt("Spglib:invalidMagSymprec",
                          "Magnetic symmetry precision must be a scalar.");
    }
    double mag_symprec = mxGetScalar(prhs[11]);

    // Initialize the rotation, translation, equivalent_atoms, and
    // primitive_lattice arrays
    mexutil::Buffer3D<int, 3, 3> rotation(max_size);
    mexutil::Buffer2D<double, 3> translation(max_size);
    mexutil::Buffer1D<int> equivalent_atoms(num_atom);
    double primitive_lattice[3][3];
    mexutil::Buffer1D<int> spin_flips(max_size);

    // Call spgms_get_symmetry_with_site_tensors
    int n_operations = spgms_get_symmetry_with_site_tensors(
        rotation, translation, equivalent_atoms, primitive_lattice, spin_flips,
        max_size, lattice, position, types, tensors, tensor_rank, num_atom,
        with_time_reversal, is_axial, symprec, angle_tolerance, mag_symprec);

    if (n_operations == 0) {
        throwLastSpglibError("spgms_get_symmetry_with_site_tensors failed");
    }

    // Create and populate the output arrays
    // Output the rotation matrices (Nx3x3 int array)
    mwSize dims[3] = {static_cast<mwSize>(n_operations), 3, 3};
    plhs[0] = mxCreateNumericArray(3, dims, mxINT32_CLASS, mxREAL);
    auto *rotations_out = static_cast<int32_t *>(mxGetData(plhs[0]));

    for (int k = 0; k < n_operations; ++k) {
        for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < 3; ++j) {
                rotations_out[k + i * n_operations + j * n_operations * 3] =
                    rotation[k][i][j];
            }
        }
    }

    // Output the translations (Nx3 double array)
    plhs[1] = mxCreateDoubleMatrix(n_operations, 3, mxREAL);
    double *translations_out = mxGetPr(plhs[1]);
    for (int i = 0; i < n_operations; ++i) {
        for (int j = 0; j < 3; ++j) {
            translations_out[i + j * n_operations] = translation[i][j];
        }
    }

    // Output the equivalent atoms (num_atom int array)
    plhs[2] = mxCreateNumericMatrix(num_atom, 1, mxINT32_CLASS, mxREAL);
    int32_t *equivalent_atoms_out = static_cast<int32_t *>(mxGetData(plhs[2]));
    for (int i = 0; i < num_atom; ++i) {
        equivalent_atoms_out[i] = equivalent_atoms[i];
    }

    // Output primitive_lattice (3x3 double array)
    plhs[3] = mxCreateDoubleMatrix(3, 3, mxREAL);
    double *primitive_lattice_out = mxGetPr(plhs[3]);
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            primitive_lattice_out[i + j * 3] = primitive_lattice[i][j];
        }
    }

    // Output spin_flips (n_operations int array)
    plhs[4] = mxCreateNumericMatrix(n_operations, 1, mxINT32_CLASS, mxREAL);
    auto *spin_flips_out = static_cast<int32_t *>(mxGetData(plhs[4]));
    for (int i = 0; i < n_operations; ++i) {
        spin_flips_out[i] = spin_flips[i];
    }

    // Output the number of operations
    plhs[5] = mxCreateDoubleScalar(static_cast<double>(n_operations));
}

// spacegroup_type = symspg('spg_get_spacegroup_type_from_symmetry', rotation,
// translation, num_operations, lattice, symprec)
void SpglibFunctions::spg_get_spacegroup_type_from_symmetry_mex(
    int nlhs, mxArray *plhs[], int nrhs, mxArray const *prhs[]) {
    /*
     SpglibSpacegroupType spg_get_spacegroup_type_from_symmetry(
         int const rotation[][3][3],
         double const translation[][3],
         int const num_operations,
         double const lattice[3][3],
         double const symprec);
    */

    // Validate the number of input arguments
    int const expected_number_of_inputs = 5;
    if (nrhs != expected_number_of_inputs) {
        mexErrMsgIdAndTxt("Spglib:invalidNumInputs",
                          "Incorrect number of inputs for "
                          "spg_get_spacegroup_type_from_symmetry.");
    }

    // Extract and validate the rotation argument
    mwSize const *dims = mxGetDimensions(prhs[0]);
    int num_operations = dims[0];  // Value of N
    if (dims[1] != 3 || dims[2] != 3) {
        mexErrMsgIdAndTxt("Spglib:invalidRotation",
                          "Rotation matrix must have dimensions Nx3x3.");
    }

    mexutil::Buffer3D<int, 3, 3> rotation(num_operations);
    auto *rotation_ptr = static_cast<int32_t *>(mxGetData(prhs[0]));

    for (int k = 0; k < num_operations; k++) {
        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++) {
                rotation[k][i][j] = rotation_ptr[k + i * num_operations +
                                                 j * num_operations * 3];
            }
        }
    }

    // Extract and validate the translation argument
    mexutil::Buffer2D<double, 3> translation(num_operations);
    if (mxGetM(prhs[1]) != num_operations || mxGetN(prhs[1]) != 3) {
        mexErrMsgIdAndTxt("Spglib:invalidTranslation",
                          "Translation array must be Nx3.");
    }
    double *translation_ptr = mxGetPr(prhs[1]);
    for (int i = 0; i < num_operations; i++) {
        for (int j = 0; j < 3; j++) {
            translation[i][j] = translation_ptr[i + j * num_operations];
        }
    }

    // Extract and validate the lattice argument
    double lattice[3][3];
    if (mxGetM(prhs[3]) != 3 || mxGetN(prhs[3]) != 3) {
        mexErrMsgIdAndTxt("Spglib:invalidLattice",
                          "Lattice matrix must be 3x3.");
    }
    double *lattice_ptr = mxGetPr(prhs[3]);
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            lattice[i][j] = lattice_ptr[i + 3 * j];
        }
    }

    // Extract and validate the symprec argument
    if (!mxIsDouble(prhs[4]) || mxGetNumberOfElements(prhs[4]) != 1) {
        mexErrMsgIdAndTxt("Spglib:invalidSymprec", "Symprec must be a scalar.");
    }
    double symprec = mxGetScalar(prhs[4]);

    // Call spg_get_spacegroup_type_from_symmetry
    SpglibSpacegroupType spacegroup_type =
        spg_get_spacegroup_type_from_symmetry(rotation, translation,
                                              num_operations, lattice, symprec);

    // Create and populate the MATLAB struct
    char const *field_names[] = {"number",
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

    mexutil::setScalarField(plhs[0], 0, "number", spacegroup_type.number);
    mexutil::setStringField(plhs[0], 0, "international_short",
                            spacegroup_type.international_short);
    mexutil::setStringField(plhs[0], 0, "international_full",
                            spacegroup_type.international_full);
    mexutil::setStringField(plhs[0], 0, "international",
                            spacegroup_type.international);
    mexutil::setStringField(plhs[0], 0, "schoenflies",
                            spacegroup_type.schoenflies);
    mexutil::setScalarField(plhs[0], 0, "hall_number",
                            spacegroup_type.hall_number);
    mexutil::setStringField(plhs[0], 0, "hall_symbol",
                            spacegroup_type.hall_symbol);
    mexutil::setStringField(plhs[0], 0, "choice", spacegroup_type.choice);
    mexutil::setStringField(plhs[0], 0, "pointgroup_international",
                            spacegroup_type.pointgroup_international);
    mexutil::setStringField(plhs[0], 0, "pointgroup_schoenflies",
                            spacegroup_type.pointgroup_schoenflies);
    mexutil::setScalarField(plhs[0], 0, "arithmetic_crystal_class_number",
                            spacegroup_type.arithmetic_crystal_class_number);
    mexutil::setStringField(plhs[0], 0, "arithmetic_crystal_class_symbol",
                            spacegroup_type.arithmetic_crystal_class_symbol);
}

// magnetic_spacegroup_type =
// symspg('spg_get_magnetic_spacegroup_type_from_symmetry', rotation,
// translation, time_reversals, num_operations, lattice, symprec)
void SpglibFunctions::spg_get_magnetic_spacegroup_type_from_symmetry_mex(
    int nlhs, mxArray *plhs[], int nrhs, mxArray const *prhs[]) {
    /*
     SpglibMagneticSpacegroupType
     spg_get_magnetic_spacegroup_type_from_symmetry( int const
     rotations[][3][3], double const translations[][3], int const
     *time_reversals, int const num_operations, double const lattice[3][3],
         double const symprec);
    */

    // Validate the number of input arguments
    int const expected_number_of_inputs = 6;
    if (nrhs != expected_number_of_inputs) {
        mexErrMsgIdAndTxt("Spglib:invalidNumInputs",
                          "Incorrect number of inputs for "
                          "spg_get_magnetic_spacegroup_type_from_symmetry.");
    }

    // Extract and validate the rotation argument
    mwSize const *dims = mxGetDimensions(prhs[0]);
    int num_operations = dims[0];  // Value of N
    if (dims[1] != 3 || dims[2] != 3) {
        mexErrMsgIdAndTxt("Spglib:invalidRotation",
                          "Rotation matrix must have dimensions Nx3x3.");
    }

    mexutil::Buffer3D<int, 3, 3> rotation(num_operations);
    auto *rotation_ptr = static_cast<int32_t *>(mxGetData(prhs[0]));

    for (int k = 0; k < num_operations; k++) {
        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++) {
                rotation[k][i][j] = rotation_ptr[k + i * num_operations +
                                                 j * num_operations * 3];
            }
        }
    }

    // Extract and validate the translation argument
    mexutil::Buffer2D<double, 3> translation(num_operations);
    if (mxGetM(prhs[1]) != num_operations || mxGetN(prhs[1]) != 3) {
        mexErrMsgIdAndTxt("Spglib:invalidTranslation",
                          "Translation array must be Nx3.");
    }
    double *translation_ptr = mxGetPr(prhs[1]);
    for (int i = 0; i < num_operations; i++) {
        for (int j = 0; j < 3; j++) {
            translation[i][j] = translation_ptr[i + j * num_operations];
        }
    }

    // Extract and validate the time_reversals argument
    if (mxGetNumberOfElements(prhs[2]) != num_operations) {
        mexErrMsgIdAndTxt(
            "Spglib:invalidTimeReversals",
            "Time reversals array size must match the number of operations.");
    }
    int const *time_reversals =
        static_cast<int32_t const *>(mxGetData(prhs[2]));

    // Extract and validate the lattice argument
    double lattice[3][3];
    if (mxGetM(prhs[4]) != 3 || mxGetN(prhs[4]) != 3) {
        mexErrMsgIdAndTxt("Spglib:invalidLattice",
                          "Lattice matrix must be 3x3.");
    }
    double *lattice_ptr = mxGetPr(prhs[4]);
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            lattice[i][j] = lattice_ptr[i + 3 * j];
        }
    }

    // Extract and validate the symprec argument
    if (!mxIsDouble(prhs[5]) || mxGetNumberOfElements(prhs[5]) != 1) {
        mexErrMsgIdAndTxt("Spglib:invalidSymprec", "Symprec must be a scalar.");
    }
    double symprec = mxGetScalar(prhs[5]);

    // Call spg_get_magnetic_spacegroup_type_from_symmetry
    SpglibMagneticSpacegroupType magnetic_spacegroup_type =
        spg_get_magnetic_spacegroup_type_from_symmetry(
            rotation, translation, time_reversals, num_operations, lattice,
            symprec);

    // Create and populate the MATLAB struct
    char const *field_names[] = {"uni_number", "litvin_number", "bns_number",
                                 "og_number",  "number",        "type"};
    plhs[0] = mxCreateStructMatrix(1, 1, 6, field_names);

    mexutil::setScalarField(plhs[0], 0, "uni_number",
                            magnetic_spacegroup_type.uni_number);
    mexutil::setScalarField(plhs[0], 0, "litvin_number",
                            magnetic_spacegroup_type.litvin_number);
    mexutil::setStringField(plhs[0], 0, "bns_number",
                            magnetic_spacegroup_type.bns_number);
    mexutil::setStringField(plhs[0], 0, "og_number",
                            magnetic_spacegroup_type.og_number);
    mexutil::setScalarField(plhs[0], 0, "number",
                            magnetic_spacegroup_type.number);
    mexutil::setScalarField(plhs[0], 0, "type", magnetic_spacegroup_type.type);
}

// [symbol, trans_mat, result] =
// symspg('spg_get_pointgroup', rotations, num_rotations)
void SpglibFunctions::spg_get_pointgroup_mex(int nlhs, mxArray *plhs[],
                                             int nrhs, mxArray const *prhs[]) {
    /*
     int spg_get_pointgroup(char symbol[6], int trans_mat[3][3],
                            int const rotations[][3][3], int const
     num_rotations);
    */

    // Validate the number of input arguments
    int const expected_number_of_inputs = 2;
    if (nrhs != expected_number_of_inputs) {
        mexErrMsgIdAndTxt("Spglib:invalidNumInputs",
                          "Incorrect number of inputs for spg_get_pointgroup.");
    }

    // Extract and validate the rotation argument
    mwSize const *dims = mxGetDimensions(prhs[0]);
    int num_operations = dims[0];
    if (dims[1] != 3 || dims[2] != 3) {
        mexErrMsgIdAndTxt("Spglib:invalidRotation",
                          "Rotation matrix must have dimensions Nx3x3.");
    }

    mexutil::Buffer3D<int, 3, 3> rotation(num_operations);
    auto *rotation_ptr = static_cast<int32_t *>(mxGetData(prhs[0]));

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

    // Call spg_get_pointgroup
    int result =
        spg_get_pointgroup(symbol, trans_mat, rotation, num_operations);

    // Create and populate the output arrays
    // Set the symbol field
    plhs[0] = mxCreateString(symbol);

    // Output the transformation matrix (3x3 int array)
    plhs[1] = mxCreateNumericMatrix(3, 3, mxINT32_CLASS, mxREAL);
    int *trans_mat_out = static_cast<int *>(mxGetData(plhs[1]));
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            trans_mat_out[i + j * 3] = trans_mat[i][j];
        }
    }

    // Set the result field
    plhs[2] = mxCreateDoubleScalar(static_cast<double>(result));
}

// [rotations, translations] = symspg('spg_get_symmetry_from_database',
// hall_number)
void SpglibFunctions::spg_get_symmetry_from_database_mex(
    int nlhs, mxArray *plhs[], int nrhs, mxArray const *prhs[]) {
    /*
     int spg_get_symmetry_from_database(int rotations[192][3][3],
                                        double translations[192][3],
                                        int const hall_number);
    */

    // Validate the number of input arguments
    int const expected_number_of_inputs = 1;
    if (nrhs != expected_number_of_inputs) {
        mexErrMsgIdAndTxt(
            "Spglib:invalidNumInputs",
            "Incorrect number of inputs for spg_get_symmetry_from_database.");
    }

    // Extract and validate the hall_number argument
    if (mxGetNumberOfElements(prhs[0]) != 1) {
        mexErrMsgIdAndTxt("Spglib:invalidHallNumber",
                          "Hall number must be a scalar.");
    }
    int hall_number = static_cast<int>(mxGetScalar(prhs[0]));

    // Prepare the output arrays
    int rotations[192][3][3];
    double translations[192][3];

    // Call spg_get_symmetry_from_database
    int num_operations =
        spg_get_symmetry_from_database(rotations, translations, hall_number);

    // Create the output rotation array (Nx3x3 int array)
    mwSize rotation_dims[3] = {static_cast<mwSize>(num_operations), 3, 3};
    plhs[0] = mxCreateNumericArray(3, rotation_dims, mxINT32_CLASS, mxREAL);
    int *rotations_out = static_cast<int *>(mxGetData(plhs[0]));
    for (int k = 0; k < num_operations; ++k) {
        for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < 3; ++j) {
                rotations_out[k + i * num_operations + j * num_operations * 3] =
                    rotations[k][i][j];
            }
        }
    }

    // Create the output translation array (Nx3 double array)
    plhs[1] = mxCreateDoubleMatrix(num_operations, 3, mxREAL);
    double *translations_out = mxGetPr(plhs[1]);
    for (int i = 0; i < num_operations; ++i) {
        for (int j = 0; j < 3; ++j) {
            translations_out[i + j * num_operations] = translations[i][j];
        }
    }
}

// [rotations, translations, time_reversals] =
// symspg('spg_get_magnetic_symmetry_from_database', uni_number, hall_number)
void SpglibFunctions::spg_get_magnetic_symmetry_from_database_mex(
    int nlhs, mxArray *plhs[], int nrhs, mxArray const *prhs[]) {
    /*
     int spg_get_magnetic_symmetry_from_database(int rotations[384][3][3],
                                                 double translations[384][3],
                                                 int time_reversals[384],
                                                 int const uni_number,
                                                 int const hall_number);
    */

    // Validate the number of input arguments
    int const expected_number_of_inputs = 2;
    if (nrhs != expected_number_of_inputs) {
        mexErrMsgIdAndTxt("Spglib:invalidNumInputs",
                          "Incorrect number of inputs for "
                          "spg_get_magnetic_symmetry_from_database.");
    }

    // Extract and validate the uni_number argument
    if (mxGetNumberOfElements(prhs[0]) != 1) {
        mexErrMsgIdAndTxt("Spglib:invalidUniNumber",
                          "Uni number must be a scalar.");
    }
    int uni_number = static_cast<int>(mxGetScalar(prhs[0]));

    // Extract and validate the hall_number argument
    if (mxGetNumberOfElements(prhs[1]) != 1) {
        mexErrMsgIdAndTxt("Spglib:invalidHallNumber",
                          "Hall number must be a scalar.");
    }
    int hall_number = static_cast<int>(mxGetScalar(prhs[1]));

    // Prepare the output arrays
    int rotations[384][3][3];
    double translations[384][3];
    int time_reversals[384];

    // Call spg_get_magnetic_symmetry_from_database
    int num_operations = spg_get_magnetic_symmetry_from_database(
        rotations, translations, time_reversals, uni_number, hall_number);

    // Create the output rotation array (Nx3x3 int array)
    mwSize rotation_dims[3] = {static_cast<mwSize>(num_operations), 3, 3};
    plhs[0] = mxCreateNumericArray(3, rotation_dims, mxINT32_CLASS, mxREAL);
    int *rotations_out = static_cast<int *>(mxGetData(plhs[0]));
    for (int k = 0; k < num_operations; ++k) {
        for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < 3; ++j) {
                rotations_out[k + i * num_operations + j * num_operations * 3] =
                    rotations[k][i][j];
            }
        }
    }

    // Create the output translation array (Nx3 double array)
    plhs[1] = mxCreateDoubleMatrix(num_operations, 3, mxREAL);
    double *translations_out = mxGetPr(plhs[1]);
    for (int i = 0; i < num_operations; ++i) {
        for (int j = 0; j < 3; ++j) {
            translations_out[i + j * num_operations] = translations[i][j];
        }
    }

    // Create the output time_reversals array (N int array)
    plhs[2] = mxCreateNumericMatrix(num_operations, 1, mxINT32_CLASS, mxREAL);
    int *time_reversals_out = static_cast<int *>(mxGetData(plhs[2]));
    for (int i = 0; i < num_operations; ++i) {
        time_reversals_out[i] = time_reversals[i];
    }
}

// spacegroup = symspg('spg_get_spacegroup_type', hall_number)
void SpglibFunctions::spg_get_spacegroup_type_mex(int nlhs, mxArray *plhs[],
                                                  int nrhs,
                                                  mxArray const *prhs[]) {
    /*
     SpglibSpacegroupType spg_get_spacegroup_type(int const hall_number);
    */

    // Validate the number of input arguments
    int const expected_number_of_inputs = 1;
    if (nrhs != expected_number_of_inputs) {
        mexErrMsgIdAndTxt(
            "Spglib:invalidNumInputs",
            "Incorrect number of inputs for spg_get_spacegroup_type.");
    }

    // Extract and validate the hall_number argument
    if (mxGetNumberOfElements(prhs[0]) != 1) {
        mexErrMsgIdAndTxt("Spglib:invalidHallNumber",
                          "Hall number must be a scalar.");
    }
    int hall_number = static_cast<int>(mxGetScalar(prhs[0]));

    // Call spg_get_spacegroup_type
    SpglibSpacegroupType spacegroup = spg_get_spacegroup_type(hall_number);

    // Create and populate the MATLAB struct
    char const *field_names[] = {"number",
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

    // Populate the struct fields
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
    int nlhs, mxArray *plhs[], int nrhs, mxArray const *prhs[]) {
    /*
     SpglibMagneticSpacegroupType spg_get_magnetic_spacegroup_type(int const
     uni_number);
    */

    // Validate the number of input arguments
    int const expected_number_of_inputs = 1;
    if (nrhs != expected_number_of_inputs) {
        mexErrMsgIdAndTxt(
            "Spglib:invalidNumInputs",
            "Incorrect number of inputs for spg_get_magnetic_spacegroup_type.");
    }

    // Extract and validate the uni_number argument
    if (mxGetNumberOfElements(prhs[0]) != 1) {
        mexErrMsgIdAndTxt("Spglib:invalidUniNumber",
                          "Uni number must be a scalar.");
    }
    int uni_number = static_cast<int>(mxGetScalar(prhs[0]));

    // Call spg_get_magnetic_spacegroup_type
    SpglibMagneticSpacegroupType spacegroup =
        spg_get_magnetic_spacegroup_type(uni_number);

    // Create and populate the MATLAB struct
    char const *field_names[] = {"uni_number", "litvin_number", "bns_number",
                                 "og_number",  "number",        "type"};
    plhs[0] = mxCreateStructMatrix(1, 1, 6, field_names);

    // Populate the struct fields
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
void SpglibFunctions::spg_standardize_cell_mex(int nlhs, mxArray *plhs[],
                                               int nrhs,
                                               mxArray const *prhs[]) {
    /*
     int spg_standardize_cell(double lattice[3][3], double position[][3],
                              int types[], int const num_atom,
                              int const to_primitive, int const no_idealize,
                              double const symprec);
    */

    // Validate the number of input arguments
    int const expected_number_of_inputs = 7;
    if (nrhs != expected_number_of_inputs) {
        mexErrMsgIdAndTxt(
            "Spglib:invalidNumInputs",
            "Incorrect number of inputs for spg_standardize_cell.");
    }

    // Extract and validate the lattice argument
    if (mxGetM(prhs[0]) != 3 || mxGetN(prhs[0]) != 3) {
        mexErrMsgIdAndTxt("Spglib:invalidLattice",
                          "Lattice must be a 3x3 matrix.");
    }
    double *lattice_ptr = mxGetPr(prhs[0]);
    double lattice[3][3];
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            lattice[i][j] = lattice_ptr[i + 3 * j];
        }
    }

    // Extract and validate the position argument
    mwSize num_atom = mxGetM(prhs[1]);
    if (mxGetN(prhs[1]) != 3) {
        mexErrMsgIdAndTxt("Spglib:invalidPosition",
                          "Position must be an Nx3 matrix.");
    }
    validateNumAtoms(prhs[3], num_atom);
    double *position_ptr = mxGetPr(prhs[1]);
    mexutil::Buffer2D<double, 3> position(4 * num_atom);
    for (mwSize i = 0; i < num_atom; ++i) {
        for (int j = 0; j < 3; ++j) {
            position[i][j] = position_ptr[i + j * num_atom];
        }
    }

    // Extract and validate the types argument
    if (mxGetNumberOfElements(prhs[2]) != num_atom) {
        mexErrMsgIdAndTxt("Spglib:invalidTypes",
                          "Types array size must match the number of atoms.");
    }
    int *types_ptr = static_cast<int *>(mxGetData(prhs[2]));
    mexutil::Buffer1D<int> types(4 * num_atom);
    for (mwSize i = 0; i < num_atom; ++i) {
        types[i] = types_ptr[i];
    }

    // Extract and validate the remaining scalar arguments
    int to_primitive = static_cast<int>(mxGetScalar(prhs[4]));
    int no_idealize = static_cast<int>(mxGetScalar(prhs[5]));
    double symprec = mxGetScalar(prhs[6]);

    // Call spg_standardize_cell
    int num_primitive_atom = spg_standardize_cell(
        lattice, position, types, num_atom, to_primitive, no_idealize, symprec);

    // Create the output lattice array (3x3 double array)
    plhs[0] = mxCreateDoubleMatrix(3, 3, mxREAL);
    double *lattice_out = mxGetPr(plhs[0]);
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            lattice_out[i + 3 * j] = lattice[i][j];
        }
    }

    // Create the output position array (num_primitive_atom x 3 double array)
    plhs[1] = mxCreateDoubleMatrix(num_primitive_atom, 3, mxREAL);
    double *position_out = mxGetPr(plhs[1]);
    for (int i = 0; i < num_primitive_atom; ++i) {
        for (int j = 0; j < 3; ++j) {
            position_out[i + j * num_primitive_atom] = position[i][j];
        }
    }

    // Create the output types array (num_primitive_atom int array)
    plhs[2] =
        mxCreateNumericMatrix(num_primitive_atom, 1, mxINT32_CLASS, mxREAL);
    int *types_out = static_cast<int *>(mxGetData(plhs[2]));
    for (int i = 0; i < num_primitive_atom; ++i) {
        types_out[i] = types[i];
    }

    // Create the output num_primitive_atom scalar
    plhs[3] = mxCreateDoubleScalar(static_cast<double>(num_primitive_atom));
}

// [lattice, position, types, num_primitive_atom] =
// symspg('spgat_standardize_cell', lattice, position, types, num_atom,
// to_primitive, no_idealize, symprec, angle_tolerance)
void SpglibFunctions::spgat_standardize_cell_mex(int nlhs, mxArray *plhs[],
                                                 int nrhs,
                                                 mxArray const *prhs[]) {
    /*
     int spgat_standardize_cell(double lattice[3][3], double position[][3],
                                int types[], int const num_atom,
                                int const to_primitive, int const no_idealize,
                                double const symprec, double const
     angle_tolerance);
    */

    // Validate the number of input arguments
    int const expected_number_of_inputs = 8;
    if (nrhs != expected_number_of_inputs) {
        mexErrMsgIdAndTxt(
            "Spglib:invalidNumInputs",
            "Incorrect number of inputs for spgat_standardize_cell.");
    }

    // Extract and validate the lattice argument
    if (mxGetM(prhs[0]) != 3 || mxGetN(prhs[0]) != 3) {
        mexErrMsgIdAndTxt("Spglib:invalidLattice",
                          "Lattice must be a 3x3 matrix.");
    }
    double *lattice_ptr = mxGetPr(prhs[0]);
    double lattice[3][3];
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            lattice[i][j] = lattice_ptr[i + 3 * j];
        }
    }

    // Extract and validate the position argument
    mwSize num_atom = mxGetM(prhs[1]);
    if (mxGetN(prhs[1]) != 3) {
        mexErrMsgIdAndTxt("Spglib:invalidPosition",
                          "Position must be an Nx3 matrix.");
    }
    validateNumAtoms(prhs[3], num_atom);
    double *position_ptr = mxGetPr(prhs[1]);
    mexutil::Buffer2D<double, 3> position(4 * num_atom);
    for (mwSize i = 0; i < num_atom; ++i) {
        for (int j = 0; j < 3; ++j) {
            position[i][j] = position_ptr[i + j * num_atom];
        }
    }

    // Extract and validate the types argument
    if (mxGetNumberOfElements(prhs[2]) != num_atom) {
        mexErrMsgIdAndTxt("Spglib:invalidTypes",
                          "Types array size must match the number of atoms.");
    }
    int *types_ptr = static_cast<int *>(mxGetData(prhs[2]));
    mexutil::Buffer1D<int> types(4 * num_atom);
    for (mwSize i = 0; i < num_atom; ++i) {
        types[i] = types_ptr[i];
    }

    // Extract and validate the remaining scalar arguments
    int to_primitive = static_cast<int>(mxGetScalar(prhs[4]));
    int no_idealize = static_cast<int>(mxGetScalar(prhs[5]));
    double symprec = mxGetScalar(prhs[6]);
    double angle_tolerance = mxGetScalar(prhs[7]);

    // Call spgat_standardize_cell
    int num_primitive_atom =
        spgat_standardize_cell(lattice, position, types, num_atom, to_primitive,
                               no_idealize, symprec, angle_tolerance);

    // Create the output lattice array (3x3 double array)
    plhs[0] = mxCreateDoubleMatrix(3, 3, mxREAL);
    double *lattice_out = mxGetPr(plhs[0]);
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            lattice_out[i + 3 * j] = lattice[i][j];
        }
    }

    // Create the output position array (num_primitive_atom x 3 double array)
    plhs[1] = mxCreateDoubleMatrix(num_primitive_atom, 3, mxREAL);
    double *position_out = mxGetPr(plhs[1]);
    for (int i = 0; i < num_primitive_atom; ++i) {
        for (int j = 0; j < 3; ++j) {
            position_out[i + j * num_primitive_atom] = position[i][j];
        }
    }

    // Create the output types array (num_primitive_atom int array)
    plhs[2] =
        mxCreateNumericMatrix(num_primitive_atom, 1, mxINT32_CLASS, mxREAL);
    int *types_out = static_cast<int *>(mxGetData(plhs[2]));
    for (int i = 0; i < num_primitive_atom; ++i) {
        types_out[i] = types[i];
    }

    // Create the output num_primitive_atom scalar
    plhs[3] = mxCreateDoubleScalar(static_cast<double>(num_primitive_atom));
}

// [lattice, position, types, num_primitive_atom] = symspg('spg_find_primitive',
// lattice, position, types, num_atom, symprec)
void SpglibFunctions::spg_find_primitive_mex(int nlhs, mxArray *plhs[],
                                             int nrhs, mxArray const *prhs[]) {
    /*
     int spg_find_primitive(double lattice[3][3], double position[][3],
                            int types[], int const num_atom,
                            double const symprec);
    */

    // Validate the number of input arguments
    int const expected_number_of_inputs = 5;
    if (nrhs != expected_number_of_inputs) {
        mexErrMsgIdAndTxt("Spglib:invalidNumInputs",
                          "Incorrect number of inputs for spg_find_primitive.");
    }

    // Extract and validate the lattice argument
    if (mxGetM(prhs[0]) != 3 || mxGetN(prhs[0]) != 3) {
        mexErrMsgIdAndTxt("Spglib:invalidLattice",
                          "Lattice must be a 3x3 matrix.");
    }
    double *lattice_ptr = mxGetPr(prhs[0]);
    double lattice[3][3];
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            lattice[i][j] = lattice_ptr[i + 3 * j];
        }
    }

    // Extract and validate the position argument
    mwSize num_atom = mxGetM(prhs[1]);
    if (mxGetN(prhs[1]) != 3) {
        mexErrMsgIdAndTxt("Spglib:invalidPosition",
                          "Position must be an Nx3 matrix.");
    }
    double *position_ptr = mxGetPr(prhs[1]);
    mexutil::Buffer2D<double, 3> position(num_atom);
    for (mwSize i = 0; i < num_atom; ++i) {
        for (int j = 0; j < 3; ++j) {
            position[i][j] = position_ptr[i + j * num_atom];
        }
    }

    // Extract and validate the types argument
    if (mxGetNumberOfElements(prhs[2]) != num_atom) {
        mexErrMsgIdAndTxt("Spglib:invalidTypes",
                          "Types array size must match the number of atoms.");
    }
    int *types_ptr = static_cast<int *>(mxGetData(prhs[2]));
    mexutil::Buffer1D<int> types(num_atom);
    for (mwSize i = 0; i < num_atom; ++i) {
        types[i] = types_ptr[i];
    }

    // Extract and validate the symprec argument
    double symprec = mxGetScalar(prhs[4]);

    // Call spg_find_primitive
    int num_primitive_atom =
        spg_find_primitive(lattice, position, types, num_atom, symprec);

    // Create the output lattice array (3x3 double array)
    plhs[0] = mxCreateDoubleMatrix(3, 3, mxREAL);
    double *lattice_out = mxGetPr(plhs[0]);
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            lattice_out[i + 3 * j] = lattice[i][j];
        }
    }

    // Create the output position array (num_primitive_atom x 3 double array)
    plhs[1] = mxCreateDoubleMatrix(num_primitive_atom, 3, mxREAL);
    double *position_out = mxGetPr(plhs[1]);
    for (int i = 0; i < num_primitive_atom; ++i) {
        for (int j = 0; j < 3; ++j) {
            position_out[i + j * num_primitive_atom] = position[i][j];
        }
    }

    // Create the output types array (num_primitive_atom int array)
    plhs[2] =
        mxCreateNumericMatrix(num_primitive_atom, 1, mxINT32_CLASS, mxREAL);
    int *types_out = static_cast<int *>(mxGetData(plhs[2]));
    for (int i = 0; i < num_primitive_atom; ++i) {
        types_out[i] = types[i];
    }

    // Create the output num_primitive_atom scalar
    plhs[3] = mxCreateDoubleScalar(static_cast<double>(num_primitive_atom));
}

// [lattice, position, types, num_primitive_atom] =
// symspg('spgat_find_primitive', lattice, position, types, num_atom, symprec,
// angle_tolerance)
void SpglibFunctions::spgat_find_primitive_mex(int nlhs, mxArray *plhs[],
                                               int nrhs,
                                               mxArray const *prhs[]) {
    /*
     int spgat_find_primitive(double lattice[3][3], double position[][3],
                              int types[], int const num_atom,
                              double const symprec, double const
     angle_tolerance);
    */

    // Validate the number of input arguments
    int const expected_number_of_inputs = 6;
    if (nrhs != expected_number_of_inputs) {
        mexErrMsgIdAndTxt(
            "Spglib:invalidNumInputs",
            "Incorrect number of inputs for spgat_find_primitive.");
    }

    // Extract and validate the lattice argument
    if (mxGetM(prhs[0]) != 3 || mxGetN(prhs[0]) != 3) {
        mexErrMsgIdAndTxt("Spglib:invalidLattice",
                          "Lattice must be a 3x3 matrix.");
    }
    double *lattice_ptr = mxGetPr(prhs[0]);
    double lattice[3][3];
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            lattice[i][j] = lattice_ptr[i + 3 * j];
        }
    }

    // Extract and validate the position argument
    mwSize num_atom = mxGetM(prhs[1]);
    if (mxGetN(prhs[1]) != 3) {
        mexErrMsgIdAndTxt("Spglib:invalidPosition",
                          "Position must be an Nx3 matrix.");
    }
    double *position_ptr = mxGetPr(prhs[1]);
    mexutil::Buffer2D<double, 3> position(num_atom);
    for (mwSize i = 0; i < num_atom; ++i) {
        for (int j = 0; j < 3; ++j) {
            position[i][j] = position_ptr[i + j * num_atom];
        }
    }

    // Extract and validate the types argument
    if (mxGetNumberOfElements(prhs[2]) != num_atom) {
        mexErrMsgIdAndTxt("Spglib:invalidTypes",
                          "Types array size must match the number of atoms.");
    }
    int *types_ptr = static_cast<int *>(mxGetData(prhs[2]));
    mexutil::Buffer1D<int> types(num_atom);
    for (mwSize i = 0; i < num_atom; ++i) {
        types[i] = types_ptr[i];
    }

    // Extract and validate the symprec argument
    double symprec = mxGetScalar(prhs[4]);

    // Extract and validate the angle_tolerance argument
    double angle_tolerance = mxGetScalar(prhs[5]);

    // Call spgat_find_primitive
    int num_primitive_atom = spgat_find_primitive(
        lattice, position, types, num_atom, symprec, angle_tolerance);

    // Create the output lattice array (3x3 double array)
    plhs[0] = mxCreateDoubleMatrix(3, 3, mxREAL);
    double *lattice_out = mxGetPr(plhs[0]);
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            lattice_out[i + 3 * j] = lattice[i][j];
        }
    }

    // Create the output position array (num_primitive_atom x 3 double array)
    plhs[1] = mxCreateDoubleMatrix(num_primitive_atom, 3, mxREAL);
    double *position_out = mxGetPr(plhs[1]);
    for (int i = 0; i < num_primitive_atom; ++i) {
        for (int j = 0; j < 3; ++j) {
            position_out[i + j * num_primitive_atom] = position[i][j];
        }
    }

    // Create the output types array (num_primitive_atom int array)
    plhs[2] =
        mxCreateNumericMatrix(num_primitive_atom, 1, mxINT32_CLASS, mxREAL);
    int *types_out = static_cast<int *>(mxGetData(plhs[2]));
    for (int i = 0; i < num_primitive_atom; ++i) {
        types_out[i] = types[i];
    }

    // Create the output num_primitive_atom scalar
    plhs[3] = mxCreateDoubleScalar(static_cast<double>(num_primitive_atom));
}

// [lattice, position, types, num_atom_bravais] = symspg('spg_refine_cell',
// lattice, position, types, num_atom, symprec)
void SpglibFunctions::spg_refine_cell_mex(int nlhs, mxArray *plhs[], int nrhs,
                                          mxArray const *prhs[]) {
    /*
     int spg_refine_cell(double lattice[3][3], double position[][3],
                         int types[], int const num_atom,
                         double const symprec);
    */

    // Validate the number of input arguments
    int const expected_number_of_inputs = 5;
    if (nrhs != expected_number_of_inputs) {
        mexErrMsgIdAndTxt("Spglib:invalidNumInputs",
                          "Incorrect number of inputs for spg_refine_cell.");
    }

    // Extract and validate the lattice argument
    if (mxGetM(prhs[0]) != 3 || mxGetN(prhs[0]) != 3) {
        mexErrMsgIdAndTxt("Spglib:invalidLattice",
                          "Lattice must be a 3x3 matrix.");
    }
    double *lattice_ptr = mxGetPr(prhs[0]);
    double lattice[3][3];
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            lattice[i][j] = lattice_ptr[i + 3 * j];
        }
    }

    // Extract and validate the position argument
    mwSize num_atom = mxGetM(prhs[1]);
    if (mxGetN(prhs[1]) != 3) {
        mexErrMsgIdAndTxt("Spglib:invalidPosition",
                          "Position must be an Nx3 matrix.");
    }
    double *position_ptr = mxGetPr(prhs[1]);
    mexutil::Buffer2D<double, 3> position(
        4 * num_atom);  // position must be allocated for
                        // 4 * num_atom rows
    for (mwSize i = 0; i < num_atom; ++i) {
        for (int j = 0; j < 3; ++j) {
            position[i][j] = position_ptr[i + j * num_atom];
        }
    }

    // Extract and validate the types argument
    if (mxGetNumberOfElements(prhs[2]) != num_atom) {
        mexErrMsgIdAndTxt("Spglib:invalidTypes",
                          "Types array size must match the number of atoms.");
    }
    int *types_ptr = static_cast<int *>(mxGetData(prhs[2]));
    mexutil::Buffer1D<int> types(
        4 * num_atom);  // types must be allocated with 4 * num_atom elements
    for (mwSize i = 0; i < num_atom; ++i) {
        types[i] = types_ptr[i];
    }

    // Extract and validate the symprec argument
    double symprec = mxGetScalar(prhs[4]);

    // Call spg_refine_cell
    int num_atom_bravais =
        spg_refine_cell(lattice, position, types, num_atom, symprec);

    // Create the output lattice array (3x3 double array)
    plhs[0] = mxCreateDoubleMatrix(3, 3, mxREAL);
    double *lattice_out = mxGetPr(plhs[0]);
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            lattice_out[i + 3 * j] = lattice[i][j];
        }
    }

    // Create the output position array (num_atom_bravais x 3 double array)
    plhs[1] = mxCreateDoubleMatrix(num_atom_bravais, 3, mxREAL);
    double *position_out = mxGetPr(plhs[1]);
    for (int i = 0; i < num_atom_bravais; ++i) {
        for (int j = 0; j < 3; ++j) {
            position_out[i + j * num_atom_bravais] = position[i][j];
        }
    }

    // Create the output types array (num_atom_bravais int array)
    plhs[2] = mxCreateNumericMatrix(num_atom_bravais, 1, mxINT32_CLASS, mxREAL);
    int *types_out = static_cast<int *>(mxGetData(plhs[2]));
    for (int i = 0; i < num_atom_bravais; ++i) {
        types_out[i] = types[i];
    }

    // Create the output num_atom_bravais scalar
    plhs[3] = mxCreateDoubleScalar(static_cast<double>(num_atom_bravais));
}

// [lattice, position, types, num_atom_bravais] = symspg('spgat_refine_cell',
// lattice, position, types, num_atom, symprec, angle_tolerance)
void SpglibFunctions::spgat_refine_cell_mex(int nlhs, mxArray *plhs[], int nrhs,
                                            mxArray const *prhs[]) {
    /*
     int spgat_refine_cell(double lattice[3][3], double position[][3],
                           int types[], int const num_atom,
                           double const symprec, double const angle_tolerance);
    */

    // Validate the number of input arguments
    int const expected_number_of_inputs = 6;
    if (nrhs != expected_number_of_inputs) {
        mexErrMsgIdAndTxt("Spglib:invalidNumInputs",
                          "Incorrect number of inputs for spgat_refine_cell.");
    }

    // Extract and validate the lattice argument
    if (mxGetM(prhs[0]) != 3 || mxGetN(prhs[0]) != 3) {
        mexErrMsgIdAndTxt("Spglib:invalidLattice",
                          "Lattice must be a 3x3 matrix.");
    }
    double *lattice_ptr = mxGetPr(prhs[0]);
    double lattice[3][3];
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            lattice[i][j] = lattice_ptr[i + 3 * j];
        }
    }

    // Extract and validate the position argument
    mwSize num_atom = mxGetM(prhs[1]);
    if (mxGetN(prhs[1]) != 3) {
        mexErrMsgIdAndTxt("Spglib:invalidPosition",
                          "Position must be an Nx3 matrix.");
    }
    double *position_ptr = mxGetPr(prhs[1]);
    mexutil::Buffer2D<double, 3> position(
        4 * num_atom);  // position must be allocated for
                        // 4 * num_atom rows
    for (mwSize i = 0; i < num_atom; ++i) {
        for (int j = 0; j < 3; ++j) {
            position[i][j] = position_ptr[i + j * num_atom];
        }
    }

    // Extract and validate the types argument
    if (mxGetNumberOfElements(prhs[2]) != num_atom) {
        mexErrMsgIdAndTxt("Spglib:invalidTypes",
                          "Types array size must match the number of atoms.");
    }
    int *types_ptr = static_cast<int *>(mxGetData(prhs[2]));
    mexutil::Buffer1D<int> types(
        4 * num_atom);  // types must be allocated with 4 * num_atom elements
    for (mwSize i = 0; i < num_atom; ++i) {
        types[i] = types_ptr[i];
    }

    // Extract and validate the symprec argument
    double symprec = mxGetScalar(prhs[4]);

    // Extract and validate the angle_tolerance argument
    double angle_tolerance = mxGetScalar(prhs[5]);

    // Call spgat_refine_cell
    int num_atom_bravais = spgat_refine_cell(lattice, position, types, num_atom,
                                             symprec, angle_tolerance);

    // Create the output lattice array (3x3 double array)
    plhs[0] = mxCreateDoubleMatrix(3, 3, mxREAL);
    double *lattice_out = mxGetPr(plhs[0]);
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            lattice_out[i + 3 * j] = lattice[i][j];
        }
    }

    // Create the output position array (num_atom_bravais x 3 double array)
    plhs[1] = mxCreateDoubleMatrix(num_atom_bravais, 3, mxREAL);
    double *position_out = mxGetPr(plhs[1]);
    for (int i = 0; i < num_atom_bravais; ++i) {
        for (int j = 0; j < 3; ++j) {
            position_out[i + j * num_atom_bravais] = position[i][j];
        }
    }

    // Create the output types array (num_atom_bravais int array)
    plhs[2] = mxCreateNumericMatrix(num_atom_bravais, 1, mxINT32_CLASS, mxREAL);
    int *types_out = static_cast<int *>(mxGetData(plhs[2]));
    for (int i = 0; i < num_atom_bravais; ++i) {
        types_out[i] = types[i];
    }

    // Create the output num_atom_bravais scalar
    plhs[3] = mxCreateDoubleScalar(static_cast<double>(num_atom_bravais));
}

// [lattice, result] = symspg('spg_delaunay_reduce', lattice, symprec)
void SpglibFunctions::spg_delaunay_reduce_mex(int nlhs, mxArray *plhs[],
                                              int nrhs, mxArray const *prhs[]) {
    /*
     int spg_delaunay_reduce(double lattice[3][3], double const symprec);
    */

    // Validate the number of input arguments
    int const expected_number_of_inputs = 2;
    if (nrhs != expected_number_of_inputs) {
        mexErrMsgIdAndTxt(
            "Spglib:invalidNumInputs",
            "Incorrect number of inputs for spg_delaunay_reduce.");
    }

    // Extract and validate the lattice argument
    if (mxGetM(prhs[0]) != 3 || mxGetN(prhs[0]) != 3) {
        mexErrMsgIdAndTxt("Spglib:invalidLattice",
                          "Lattice must be a 3x3 matrix.");
    }
    double *lattice_ptr = mxGetPr(prhs[0]);
    double lattice[3][3];
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            lattice[i][j] = lattice_ptr[i + 3 * j];
        }
    }

    // Extract and validate the symprec argument
    double symprec = mxGetScalar(prhs[1]);

    // Call spg_delaunay_reduce
    int result = spg_delaunay_reduce(lattice, symprec);

    // Create the output lattice array (3x3 double array)
    plhs[0] = mxCreateDoubleMatrix(3, 3, mxREAL);
    double *lattice_out = mxGetPr(plhs[0]);
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            lattice_out[i + 3 * j] = lattice[i][j];
        }
    }

    // Create the output result scalar
    plhs[1] = mxCreateDoubleScalar(static_cast<double>(result));
}

// grid_point_index = symspg('spg_get_grid_point_from_address', grid_address,
// mesh)
void SpglibFunctions::spg_get_grid_point_from_address_mex(
    int nlhs, mxArray *plhs[], int nrhs, mxArray const *prhs[]) {
    /*
     int spg_get_grid_point_from_address(int const grid_address[3], int const
     mesh[3]);
    */

    // Validate the number of input arguments
    int const expected_number_of_inputs = 2;
    if (nrhs != expected_number_of_inputs) {
        mexErrMsgIdAndTxt(
            "Spglib:invalidNumInputs",
            "Incorrect number of inputs for spg_get_grid_point_from_address.");
    }

    // Extract and validate the grid_address argument
    if (mxGetNumberOfElements(prhs[0]) != 3) {
        mexErrMsgIdAndTxt("Spglib:invalidGridAddress",
                          "Grid address must be an array of 3 elements.");
    }
    int *grid_address_ptr = static_cast<int *>(mxGetData(prhs[0]));
    int grid_address[3] = {grid_address_ptr[0], grid_address_ptr[1],
                           grid_address_ptr[2]};

    // Extract and validate the mesh argument
    if (mxGetNumberOfElements(prhs[1]) != 3) {
        mexErrMsgIdAndTxt("Spglib:invalidMesh",
                          "Mesh must be an array of 3 elements.");
    }
    int *mesh_ptr = static_cast<int *>(mxGetData(prhs[1]));
    int mesh[3] = {mesh_ptr[0], mesh_ptr[1], mesh_ptr[2]};

    // Call spg_get_grid_point_from_address
    int grid_point_index = spg_get_grid_point_from_address(grid_address, mesh);

    // Create the output grid_point_index scalar
    plhs[0] = mxCreateDoubleScalar(static_cast<double>(grid_point_index));
}

// dense_grid_point_index = symspg('spg_get_dense_grid_point_from_address',
// grid_address, mesh)
void SpglibFunctions::spg_get_dense_grid_point_from_address_mex(
    int nlhs, mxArray *plhs[], int nrhs, mxArray const *prhs[]) {
    /*
     size_t spg_get_dense_grid_point_from_address(int const grid_address[3], int
     const mesh[3]);
    */

    // Validate the number of input arguments
    int const expected_number_of_inputs = 2;
    if (nrhs != expected_number_of_inputs) {
        mexErrMsgIdAndTxt("Spglib:invalidNumInputs",
                          "Incorrect number of inputs for "
                          "spg_get_dense_grid_point_from_address.");
    }

    // Extract and validate the grid_address argument
    if (mxGetNumberOfElements(prhs[0]) != 3) {
        mexErrMsgIdAndTxt("Spglib:invalidGridAddress",
                          "Grid address must be an array of 3 elements.");
    }
    int *grid_address_ptr = static_cast<int *>(mxGetData(prhs[0]));
    int grid_address[3] = {grid_address_ptr[0], grid_address_ptr[1],
                           grid_address_ptr[2]};

    // Extract and validate the mesh argument
    if (mxGetNumberOfElements(prhs[1]) != 3) {
        mexErrMsgIdAndTxt("Spglib:invalidMesh",
                          "Mesh must be an array of 3 elements.");
    }
    int *mesh_ptr = static_cast<int *>(mxGetData(prhs[1]));
    int mesh[3] = {mesh_ptr[0], mesh_ptr[1], mesh_ptr[2]};

    // Call spg_get_dense_grid_point_from_address
    size_t dense_grid_point_index =
        spg_get_dense_grid_point_from_address(grid_address, mesh);

    // Create the output dense_grid_point_index scalar
    plhs[0] = mxCreateDoubleScalar(static_cast<double>(dense_grid_point_index));
}

// [grid_address, ir_mapping_table, num_ir_kpoints] =
// symspg('spg_get_ir_reciprocal_mesh', mesh, is_shift, is_time_reversal,
// lattice, position, types, num_atom, symprec)
void SpglibFunctions::spg_get_ir_reciprocal_mesh_mex(int nlhs, mxArray *plhs[],
                                                     int nrhs,
                                                     mxArray const *prhs[]) {
    /*
     int spg_get_ir_reciprocal_mesh(int grid_address[][3], int
     ir_mapping_table[], int const mesh[3], int const is_shift[3], int const
     is_time_reversal, double const lattice[3][3], double const position[][3],
     int const types[], int const num_atom, double const symprec);
    */

    // Validate the number of input arguments
    int const expected_number_of_inputs = 8;
    if (nrhs != expected_number_of_inputs) {
        mexErrMsgIdAndTxt(
            "Spglib:invalidNumInputs",
            "Incorrect number of inputs for spg_get_ir_reciprocal_mesh.");
    }

    // Extract and validate the mesh argument
    if (mxGetNumberOfElements(prhs[0]) != 3) {
        mexErrMsgIdAndTxt("Spglib:invalidMesh",
                          "Mesh must be an array of 3 elements.");
    }
    int *mesh_ptr = static_cast<int *>(mxGetData(prhs[0]));
    int mesh[3] = {mesh_ptr[0], mesh_ptr[1], mesh_ptr[2]};

    // Extract and validate the is_shift argument
    if (mxGetNumberOfElements(prhs[1]) != 3) {
        mexErrMsgIdAndTxt("Spglib:invalidShift",
                          "is_shift must be an array of 3 elements.");
    }
    int *is_shift_ptr = static_cast<int *>(mxGetData(prhs[1]));
    int is_shift[3] = {is_shift_ptr[0], is_shift_ptr[1], is_shift_ptr[2]};

    // Extract the is_time_reversal argument
    int is_time_reversal = static_cast<int>(mxGetScalar(prhs[2]));

    // Extract and validate the lattice argument
    if (mxGetM(prhs[3]) != 3 || mxGetN(prhs[3]) != 3) {
        mexErrMsgIdAndTxt("Spglib:invalidLattice",
                          "Lattice must be a 3x3 matrix.");
    }
    double *lattice_ptr = mxGetPr(prhs[3]);
    double lattice[3][3];
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            lattice[i][j] = lattice_ptr[i + 3 * j];
        }
    }

    // Extract and validate the position argument
    mwSize num_atom = mxGetM(prhs[4]);
    if (mxGetN(prhs[4]) != 3) {
        mexErrMsgIdAndTxt("Spglib:invalidPosition",
                          "Position must be an Nx3 matrix.");
    }
    double *position_ptr = mxGetPr(prhs[4]);
    mexutil::Buffer2D<double, 3> position(num_atom);
    for (mwSize i = 0; i < num_atom; ++i) {
        for (int j = 0; j < 3; ++j) {
            position[i][j] = position_ptr[i + j * num_atom];
        }
    }

    // Extract and validate the types argument
    if (mxGetNumberOfElements(prhs[5]) != num_atom) {
        mexErrMsgIdAndTxt("Spglib:invalidTypes",
                          "Types array size must match the number of atoms.");
    }
    int *types_ptr = static_cast<int *>(mxGetData(prhs[5]));
    mexutil::Buffer1D<int> types(num_atom);
    for (mwSize i = 0; i < num_atom; ++i) {
        types[i] = types_ptr[i];
    }

    // Extract the symprec argument
    double symprec = mxGetScalar(prhs[7]);

    // Allocate grid_address for the maximum possible number of grid points
    mexutil::Buffer2D<int, 3> grid_address(mesh[0] * mesh[1] * mesh[2]);
    mexutil::Buffer1D<int> ir_mapping_table(mesh[0] * mesh[1] * mesh[2]);

    // Call spg_get_ir_reciprocal_mesh
    int num_ir_kpoints = spg_get_ir_reciprocal_mesh(
        grid_address, ir_mapping_table, mesh, is_shift, is_time_reversal,
        lattice, position, types, num_atom, symprec);

    int const num_total_grid_points = mesh[0] * mesh[1] * mesh[2];

    // Create the output grid_address array (num_total_grid_points x 3 double
    // array)
    plhs[0] = mxCreateDoubleMatrix(num_total_grid_points, 3, mxREAL);
    double *grid_address_out = mxGetPr(plhs[0]);
    for (int i = 0; i < num_total_grid_points; ++i) {
        for (int j = 0; j < 3; ++j) {
            grid_address_out[i + j * num_total_grid_points] =
                static_cast<double>(grid_address[i][j]);
        }
    }

    // Create the output ir_mapping_table array (num_total_grid_points int
    // array)
    plhs[1] =
        mxCreateNumericMatrix(num_total_grid_points, 1, mxINT32_CLASS, mxREAL);
    int *ir_mapping_table_out = static_cast<int *>(mxGetData(plhs[1]));
    for (int i = 0; i < num_total_grid_points; ++i) {
        ir_mapping_table_out[i] = ir_mapping_table[i];
    }

    // Create the output num_ir_kpoints scalar
    plhs[2] = mxCreateDoubleScalar(static_cast<double>(num_ir_kpoints));
}

// [grid_address, ir_mapping_table, num_ir_kpoints] =
// symspg('spg_get_dense_ir_reciprocal_mesh', mesh, is_shift, is_time_reversal,
// lattice, position, types, num_atom, symprec)
void SpglibFunctions::spg_get_dense_ir_reciprocal_mesh_mex(
    int nlhs, mxArray *plhs[], int nrhs, mxArray const *prhs[]) {
    /*
     size_t spg_get_dense_ir_reciprocal_mesh(int grid_address[][3], size_t
     ir_mapping_table[], int const mesh[3], int const is_shift[3], int const
     is_time_reversal, double const lattice[3][3], double const position[][3],
     int const types[], int const num_atom, double const symprec);
    */

    // Validate the number of input arguments
    int const expected_number_of_inputs = 8;
    if (nrhs != expected_number_of_inputs) {
        mexErrMsgIdAndTxt(
            "Spglib:invalidNumInputs",
            "Incorrect number of inputs for spg_get_dense_ir_reciprocal_mesh.");
    }

    // Extract and validate the mesh argument
    if (mxGetNumberOfElements(prhs[0]) != 3) {
        mexErrMsgIdAndTxt("Spglib:invalidMesh",
                          "Mesh must be an array of 3 elements.");
    }
    int *mesh_ptr = static_cast<int *>(mxGetData(prhs[0]));
    int mesh[3] = {mesh_ptr[0], mesh_ptr[1], mesh_ptr[2]};

    // Extract and validate the is_shift argument
    if (mxGetNumberOfElements(prhs[1]) != 3) {
        mexErrMsgIdAndTxt("Spglib:invalidShift",
                          "is_shift must be an array of 3 elements.");
    }
    int *is_shift_ptr = static_cast<int *>(mxGetData(prhs[1]));
    int is_shift[3] = {is_shift_ptr[0], is_shift_ptr[1], is_shift_ptr[2]};

    // Extract the is_time_reversal argument
    int is_time_reversal = static_cast<int>(mxGetScalar(prhs[2]));

    // Extract and validate the lattice argument
    if (mxGetM(prhs[3]) != 3 || mxGetN(prhs[3]) != 3) {
        mexErrMsgIdAndTxt("Spglib:invalidLattice",
                          "Lattice must be a 3x3 matrix.");
    }
    double *lattice_ptr = mxGetPr(prhs[3]);
    double lattice[3][3];
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            lattice[i][j] = lattice_ptr[i + 3 * j];
        }
    }

    // Extract and validate the position argument
    mwSize num_atom = mxGetM(prhs[4]);
    if (mxGetN(prhs[4]) != 3) {
        mexErrMsgIdAndTxt("Spglib:invalidPosition",
                          "Position must be an Nx3 matrix.");
    }
    double *position_ptr = mxGetPr(prhs[4]);
    mexutil::Buffer2D<double, 3> position(num_atom);
    for (mwSize i = 0; i < num_atom; ++i) {
        for (int j = 0; j < 3; ++j) {
            position[i][j] = position_ptr[i + j * num_atom];
        }
    }

    // Extract and validate the types argument
    if (mxGetNumberOfElements(prhs[5]) != num_atom) {
        mexErrMsgIdAndTxt("Spglib:invalidTypes",
                          "Types array size must match the number of atoms.");
    }
    int *types_ptr = static_cast<int *>(mxGetData(prhs[5]));
    mexutil::Buffer1D<int> types(num_atom);
    for (mwSize i = 0; i < num_atom; ++i) {
        types[i] = types_ptr[i];
    }

    // Extract the symprec argument
    double symprec = mxGetScalar(prhs[7]);

    // Allocate grid_address for the maximum possible number of grid points
    size_t num_total_grid_points =
        static_cast<size_t>(mesh[0] * mesh[1] * mesh[2]);
    mexutil::Buffer2D<int, 3> grid_address(num_total_grid_points);
    mexutil::Buffer1D<size_t> ir_mapping_table(num_total_grid_points);

    // Call spg_get_dense_ir_reciprocal_mesh
    size_t num_ir_kpoints = spg_get_dense_ir_reciprocal_mesh(
        grid_address, ir_mapping_table, mesh, is_shift, is_time_reversal,
        lattice, position, types, num_atom, symprec);

    // Create the output grid_address array (num_total_grid_points x 3 double
    // array)
    plhs[0] = mxCreateDoubleMatrix(num_total_grid_points, 3, mxREAL);
    double *grid_address_out = mxGetPr(plhs[0]);
    for (size_t i = 0; i < num_total_grid_points; ++i) {
        for (int j = 0; j < 3; ++j) {
            grid_address_out[i + j * num_total_grid_points] =
                static_cast<double>(grid_address[i][j]);
        }
    }

    // Create the output ir_mapping_table array (num_total_grid_points size_t
    // array)
    plhs[1] =
        mxCreateNumericMatrix(num_total_grid_points, 1, mxUINT64_CLASS, mxREAL);
    size_t *ir_mapping_table_out = static_cast<size_t *>(mxGetData(plhs[1]));
    for (size_t i = 0; i < num_total_grid_points; ++i) {
        ir_mapping_table_out[i] = ir_mapping_table[i];
    }

    // Create the output num_ir_kpoints scalar
    plhs[2] = mxCreateDoubleScalar(static_cast<double>(num_ir_kpoints));
}

// [grid_address, ir_mapping_table, num_ir_kpoints] =
// symspg('spg_get_stabilized_reciprocal_mesh', mesh, is_shift,
// is_time_reversal, num_rot, rotations, num_q, qpoints)
void SpglibFunctions::spg_get_stabilized_reciprocal_mesh_mex(
    int nlhs, mxArray *plhs[], int nrhs, mxArray const *prhs[]) {
    /*
     int spg_get_stabilized_reciprocal_mesh(int grid_address[][3], int
     ir_mapping_table[], int const mesh[3], int const is_shift[3], int const
     is_time_reversal, int const num_rot, int const rotations[][3][3], int const
     num_q, double const qpoints[][3]);
    */

    // Validate the number of input arguments
    int const expected_number_of_inputs = 7;
    if (nrhs != expected_number_of_inputs) {
        mexErrMsgIdAndTxt("Spglib:invalidNumInputs",
                          "Incorrect number of inputs for "
                          "spg_get_stabilized_reciprocal_mesh.");
    }

    // Extract and validate the mesh argument
    if (mxGetNumberOfElements(prhs[0]) != 3) {
        mexErrMsgIdAndTxt("Spglib:invalidMesh",
                          "Mesh must be an array of 3 elements.");
    }
    int *mesh_ptr = static_cast<int *>(mxGetData(prhs[0]));
    int mesh[3] = {mesh_ptr[0], mesh_ptr[1], mesh_ptr[2]};

    // Extract and validate the is_shift argument
    if (mxGetNumberOfElements(prhs[1]) != 3) {
        mexErrMsgIdAndTxt("Spglib:invalidShift",
                          "is_shift must be an array of 3 elements.");
    }
    int *is_shift_ptr = static_cast<int *>(mxGetData(prhs[1]));
    int is_shift[3] = {is_shift_ptr[0], is_shift_ptr[1], is_shift_ptr[2]};

    // Extract the is_time_reversal argument
    int is_time_reversal = static_cast<int>(mxGetScalar(prhs[2]));

    // Extract and validate the num_rot argument
    int num_rot = static_cast<int>(mxGetScalar(prhs[3]));

    // Extract and validate the rotations argument
    if (mxGetNumberOfElements(prhs[4]) != num_rot * 9) {
        mexErrMsgIdAndTxt("Spglib:invalidRotations",
                          "Rotations must be a num_rot x 3 x 3 array.");
    }
    int *rotations_ptr = static_cast<int *>(mxGetData(prhs[4]));
    mexutil::Buffer3D<int, 3, 3> rotations(num_rot);
    for (int k = 0; k < num_rot; ++k) {
        for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < 3; ++j) {
                rotations[k][i][j] =
                    rotations_ptr[k + i * num_rot + j * num_rot * 3];
            }
        }
    }

    // Extract and validate the num_q argument
    int num_q = static_cast<int>(mxGetScalar(prhs[5]));

    // Extract and validate the qpoints argument
    if (mxGetM(prhs[6]) != num_q || mxGetN(prhs[6]) != 3) {
        mexErrMsgIdAndTxt("Spglib:invalidQPoints",
                          "Qpoints must be a num_q x 3 array.");
    }
    double *qpoints_ptr = mxGetPr(prhs[6]);
    mexutil::Buffer2D<double, 3> qpoints(num_q);
    for (int i = 0; i < num_q; ++i) {
        for (int j = 0; j < 3; ++j) {
            qpoints[i][j] = qpoints_ptr[i + j * num_q];
        }
    }

    // Allocate grid_address for the maximum possible number of grid points
    int num_total_grid_points = mesh[0] * mesh[1] * mesh[2];
    mexutil::Buffer2D<int, 3> grid_address(num_total_grid_points);
    mexutil::Buffer1D<int> ir_mapping_table(num_total_grid_points);

    // Call spg_get_stabilized_reciprocal_mesh
    int num_ir_kpoints = spg_get_stabilized_reciprocal_mesh(
        grid_address, ir_mapping_table, mesh, is_shift, is_time_reversal,
        num_rot, rotations, num_q, qpoints);

    // Create the output grid_address array (num_total_grid_points x 3 double
    // array)
    plhs[0] = mxCreateDoubleMatrix(num_total_grid_points, 3, mxREAL);
    double *grid_address_out = mxGetPr(plhs[0]);
    for (int i = 0; i < num_total_grid_points; ++i) {
        for (int j = 0; j < 3; ++j) {
            grid_address_out[i + j * num_total_grid_points] =
                static_cast<double>(grid_address[i][j]);
        }
    }

    // Create the output ir_mapping_table array (num_total_grid_points int
    // array)
    plhs[1] =
        mxCreateNumericMatrix(num_total_grid_points, 1, mxINT32_CLASS, mxREAL);
    int *ir_mapping_table_out = static_cast<int *>(mxGetData(plhs[1]));
    for (int i = 0; i < num_total_grid_points; ++i) {
        ir_mapping_table_out[i] = ir_mapping_table[i];
    }

    // Create the output num_ir_kpoints scalar
    plhs[2] = mxCreateDoubleScalar(static_cast<double>(num_ir_kpoints));
}

// [grid_address, ir_mapping_table, num_ir_kpoints] =
// symspg('spg_get_dense_stabilized_reciprocal_mesh', mesh, is_shift,
// is_time_reversal, num_rot, rotations, num_q, qpoints)
void SpglibFunctions::spg_get_dense_stabilized_reciprocal_mesh_mex(
    int nlhs, mxArray *plhs[], int nrhs, mxArray const *prhs[]) {
    /*
     size_t spg_get_dense_stabilized_reciprocal_mesh(int grid_address[][3],
     size_t ir_mapping_table[], int const mesh[3], int const is_shift[3], int
     const is_time_reversal, int const num_rot, int const rotations[][3][3], int
     const num_q, double const qpoints[][3]);
    */

    // Validate the number of input arguments
    int const expected_number_of_inputs = 7;
    if (nrhs != expected_number_of_inputs) {
        mexErrMsgIdAndTxt("Spglib:invalidNumInputs",
                          "Incorrect number of inputs for "
                          "spg_get_dense_stabilized_reciprocal_mesh.");
    }

    // Extract and validate the mesh argument
    if (mxGetNumberOfElements(prhs[0]) != 3) {
        mexErrMsgIdAndTxt("Spglib:invalidMesh",
                          "Mesh must be an array of 3 elements.");
    }
    int *mesh_ptr = static_cast<int *>(mxGetData(prhs[0]));
    int mesh[3] = {mesh_ptr[0], mesh_ptr[1], mesh_ptr[2]};

    // Extract and validate the is_shift argument
    if (mxGetNumberOfElements(prhs[1]) != 3) {
        mexErrMsgIdAndTxt("Spglib:invalidShift",
                          "is_shift must be an array of 3 elements.");
    }
    int *is_shift_ptr = static_cast<int *>(mxGetData(prhs[1]));
    int is_shift[3] = {is_shift_ptr[0], is_shift_ptr[1], is_shift_ptr[2]};

    // Extract the is_time_reversal argument
    int is_time_reversal = static_cast<int>(mxGetScalar(prhs[2]));

    // Extract and validate the num_rot argument
    int num_rot = static_cast<int>(mxGetScalar(prhs[3]));

    // Extract and validate the rotations argument
    if (mxGetNumberOfElements(prhs[4]) != num_rot * 9) {
        mexErrMsgIdAndTxt("Spglib:invalidRotations",
                          "Rotations must be a num_rot x 3 x 3 array.");
    }
    int *rotations_ptr = static_cast<int *>(mxGetData(prhs[4]));
    mexutil::Buffer3D<int, 3, 3> rotations(num_rot);
    for (int k = 0; k < num_rot; ++k) {
        for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < 3; ++j) {
                rotations[k][i][j] =
                    rotations_ptr[k + i * num_rot + j * num_rot * 3];
            }
        }
    }

    // Extract and validate the num_q argument
    int num_q = static_cast<int>(mxGetScalar(prhs[5]));

    // Extract and validate the qpoints argument
    if (mxGetM(prhs[6]) != num_q || mxGetN(prhs[6]) != 3) {
        mexErrMsgIdAndTxt("Spglib:invalidQPoints",
                          "Qpoints must be a num_q x 3 array.");
    }
    double *qpoints_ptr = mxGetPr(prhs[6]);
    mexutil::Buffer2D<double, 3> qpoints(num_q);
    for (int i = 0; i < num_q; ++i) {
        for (int j = 0; j < 3; ++j) {
            qpoints[i][j] = qpoints_ptr[i + j * num_q];
        }
    }

    // Allocate grid_address for the maximum possible number of grid points
    size_t num_total_grid_points =
        static_cast<size_t>(mesh[0] * mesh[1] * mesh[2]);
    mexutil::Buffer2D<int, 3> grid_address(num_total_grid_points);
    mexutil::Buffer1D<size_t> ir_mapping_table(num_total_grid_points);

    // Call spg_get_dense_stabilized_reciprocal_mesh
    size_t num_ir_kpoints = spg_get_dense_stabilized_reciprocal_mesh(
        grid_address, ir_mapping_table, mesh, is_shift, is_time_reversal,
        num_rot, rotations, num_q, qpoints);

    // Create the output grid_address array (num_total_grid_points x 3 double
    // array)
    plhs[0] = mxCreateDoubleMatrix(num_total_grid_points, 3, mxREAL);
    double *grid_address_out = mxGetPr(plhs[0]);
    for (size_t i = 0; i < num_total_grid_points; ++i) {
        for (int j = 0; j < 3; ++j) {
            grid_address_out[i + j * num_total_grid_points] =
                static_cast<double>(grid_address[i][j]);
        }
    }

    // Create the output ir_mapping_table array (num_total_grid_points size_t
    // array)
    plhs[1] =
        mxCreateNumericMatrix(num_total_grid_points, 1, mxUINT64_CLASS, mxREAL);
    size_t *ir_mapping_table_out = static_cast<size_t *>(mxGetData(plhs[1]));
    for (size_t i = 0; i < num_total_grid_points; ++i) {
        ir_mapping_table_out[i] = ir_mapping_table[i];
    }

    // Create the output num_ir_kpoints scalar
    plhs[2] = mxCreateDoubleScalar(static_cast<double>(num_ir_kpoints));
}

// rot_grid_points = symspg('spg_get_dense_grid_points_by_rotations',
// address_orig, num_rot, rot_reciprocal, mesh, is_shift)
void SpglibFunctions::spg_get_dense_grid_points_by_rotations_mex(
    int nlhs, mxArray *plhs[], int nrhs, mxArray const *prhs[]) {
    /*
     void spg_get_dense_grid_points_by_rotations(size_t rot_grid_points[], int
     const address_orig[3], int const num_rot, int const rot_reciprocal[][3][3],
                                                 int const mesh[3], int const
     is_shift[3]);
    */

    // Validate the number of input arguments
    int const expected_number_of_inputs = 5;
    if (nrhs != expected_number_of_inputs) {
        mexErrMsgIdAndTxt("Spglib:invalidNumInputs",
                          "Incorrect number of inputs for "
                          "spg_get_dense_grid_points_by_rotations.");
    }

    // Extract and validate the address_orig argument
    if (mxGetNumberOfElements(prhs[0]) != 3) {
        mexErrMsgIdAndTxt("Spglib:invalidAddressOrig",
                          "Address_orig must be an array of 3 elements.");
    }
    int *address_orig_ptr = static_cast<int *>(mxGetData(prhs[0]));
    int address_orig[3] = {address_orig_ptr[0], address_orig_ptr[1],
                           address_orig_ptr[2]};

    // Extract and validate the num_rot argument
    int num_rot = static_cast<int>(mxGetScalar(prhs[1]));

    // Extract and validate the rot_reciprocal argument
    if (mxGetNumberOfElements(prhs[2]) != num_rot * 9) {
        mexErrMsgIdAndTxt("Spglib:invalidRotReciprocal",
                          "Rot_reciprocal must be a num_rot x 3 x 3 array.");
    }
    int *rot_reciprocal_ptr = static_cast<int *>(mxGetData(prhs[2]));
    mexutil::Buffer3D<int, 3, 3> rot_reciprocal(num_rot);
    for (int k = 0; k < num_rot; ++k) {
        for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < 3; ++j) {
                rot_reciprocal[k][i][j] =
                    rot_reciprocal_ptr[k + i * num_rot + j * num_rot * 3];
            }
        }
    }

    // Extract and validate the mesh argument
    if (mxGetNumberOfElements(prhs[3]) != 3) {
        mexErrMsgIdAndTxt("Spglib:invalidMesh",
                          "Mesh must be an array of 3 elements.");
    }
    int *mesh_ptr = static_cast<int *>(mxGetData(prhs[3]));
    int mesh[3] = {mesh_ptr[0], mesh_ptr[1], mesh_ptr[2]};

    // Extract and validate the is_shift argument
    if (mxGetNumberOfElements(prhs[4]) != 3) {
        mexErrMsgIdAndTxt("Spglib:invalidShift",
                          "is_shift must be an array of 3 elements.");
    }
    int *is_shift_ptr = static_cast<int *>(mxGetData(prhs[4]));
    int is_shift[3] = {is_shift_ptr[0], is_shift_ptr[1], is_shift_ptr[2]};

    // Create the output rot_grid_points array (num_rot size_t array)
    plhs[0] = mxCreateNumericMatrix(num_rot, 1, mxUINT64_CLASS, mxREAL);
    size_t *rot_grid_points = static_cast<size_t *>(mxGetData(plhs[0]));

    // Call spg_get_dense_grid_points_by_rotations
    spg_get_dense_grid_points_by_rotations(
        rot_grid_points, address_orig, num_rot, rot_reciprocal, mesh, is_shift);
}

// rot_grid_points = symspg('spg_get_dense_BZ_grid_points_by_rotations',
// address_orig, num_rot, rot_reciprocal, mesh, is_shift, bz_map)
void SpglibFunctions::spg_get_dense_BZ_grid_points_by_rotations_mex(
    int nlhs, mxArray *plhs[], int nrhs, mxArray const *prhs[]) {
    /*
     void spg_get_dense_BZ_grid_points_by_rotations(size_t rot_grid_points[],
     int const address_orig[3], int const num_rot, int const
     rot_reciprocal[][3][3], int const mesh[3], int const is_shift[3], size_t
     const bz_map[]);
    */

    // Validate the number of input arguments
    int const expected_number_of_inputs = 6;
    if (nrhs != expected_number_of_inputs) {
        mexErrMsgIdAndTxt("Spglib:invalidNumInputs",
                          "Incorrect number of inputs for "
                          "spg_get_dense_BZ_grid_points_by_rotations.");
    }

    // Extract and validate the address_orig argument
    if (mxGetNumberOfElements(prhs[0]) != 3) {
        mexErrMsgIdAndTxt("Spglib:invalidAddressOrig",
                          "Address_orig must be an array of 3 elements.");
    }
    int *address_orig_ptr = static_cast<int *>(mxGetData(prhs[0]));
    int address_orig[3] = {address_orig_ptr[0], address_orig_ptr[1],
                           address_orig_ptr[2]};

    // Extract and validate the num_rot argument
    int num_rot = static_cast<int>(mxGetScalar(prhs[1]));

    // Extract and validate the rot_reciprocal argument
    if (mxGetNumberOfElements(prhs[2]) != num_rot * 9) {
        mexErrMsgIdAndTxt("Spglib:invalidRotReciprocal",
                          "Rot_reciprocal must be a num_rot x 3 x 3 array.");
    }
    int *rot_reciprocal_ptr = static_cast<int *>(mxGetData(prhs[2]));
    mexutil::Buffer3D<int, 3, 3> rot_reciprocal(num_rot);
    for (int k = 0; k < num_rot; ++k) {
        for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < 3; ++j) {
                rot_reciprocal[k][i][j] =
                    rot_reciprocal_ptr[k + i * num_rot + j * num_rot * 3];
            }
        }
    }

    // Extract and validate the mesh argument
    if (mxGetNumberOfElements(prhs[3]) != 3) {
        mexErrMsgIdAndTxt("Spglib:invalidMesh",
                          "Mesh must be an array of 3 elements.");
    }
    int *mesh_ptr = static_cast<int *>(mxGetData(prhs[3]));
    int mesh[3] = {mesh_ptr[0], mesh_ptr[1], mesh_ptr[2]};

    // Extract and validate the is_shift argument
    if (mxGetNumberOfElements(prhs[4]) != 3) {
        mexErrMsgIdAndTxt("Spglib:invalidShift",
                          "is_shift must be an array of 3 elements.");
    }
    int *is_shift_ptr = static_cast<int *>(mxGetData(prhs[4]));
    int is_shift[3] = {is_shift_ptr[0], is_shift_ptr[1], is_shift_ptr[2]};

    // Extract and validate the bz_map argument
    mwSize const bz_map_size = static_cast<mwSize>(mesh[0] * 2) *
                               static_cast<mwSize>(mesh[1] * 2) *
                               static_cast<mwSize>(mesh[2] * 2);
    if (!mxIsUint64(prhs[5]) || mxGetNumberOfElements(prhs[5]) != bz_map_size) {
        mexErrMsgIdAndTxt(
            "Spglib:invalidBZMap",
            "Bz_map must be a uint64 array with prod(mesh * 2) elements.");
    }
    size_t const *bz_map = static_cast<size_t const *>(mxGetData(prhs[5]));

    // Create the output rot_grid_points array (num_rot size_t array)
    plhs[0] = mxCreateNumericMatrix(num_rot, 1, mxUINT64_CLASS, mxREAL);
    size_t *rot_grid_points = static_cast<size_t *>(mxGetData(plhs[0]));

    // Call spg_get_dense_BZ_grid_points_by_rotations
    spg_get_dense_BZ_grid_points_by_rotations(rot_grid_points, address_orig,
                                              num_rot, rot_reciprocal, mesh,
                                              is_shift, bz_map);
}

// [bz_grid_address, bz_map, num_ir_grid_points] =
// symspg('spg_relocate_BZ_grid_address', grid_address, mesh, rec_lattice,
// is_shift)
void SpglibFunctions::spg_relocate_BZ_grid_address_mex(int nlhs,
                                                       mxArray *plhs[],
                                                       int nrhs,
                                                       mxArray const *prhs[]) {
    /*
     int spg_relocate_BZ_grid_address(int bz_grid_address[][3], int bz_map[],
                                      int const grid_address[][3], int const
     mesh[3], double const rec_lattice[3][3], int const is_shift[3]);
    */

    // Validate the number of input arguments
    int const expected_number_of_inputs = 4;
    if (nrhs != expected_number_of_inputs) {
        mexErrMsgIdAndTxt(
            "Spglib:invalidNumInputs",
            "Incorrect number of inputs for spg_relocate_BZ_grid_address.");
    }

    // Extract and validate the grid_address argument
    mwSize num_grid_points = mxGetM(prhs[0]);
    if (mxGetN(prhs[0]) != 3) {
        mexErrMsgIdAndTxt("Spglib:invalidGridAddress",
                          "Grid_address must be an Nx3 array.");
    }
    int *grid_address_ptr = static_cast<int *>(mxGetData(prhs[0]));
    mexutil::Buffer2D<int, 3> grid_address(num_grid_points);
    for (mwSize i = 0; i < num_grid_points; ++i) {
        for (int j = 0; j < 3; ++j) {
            grid_address[i][j] = grid_address_ptr[i + j * num_grid_points];
        }
    }

    // Extract and validate the mesh argument
    if (mxGetNumberOfElements(prhs[1]) != 3) {
        mexErrMsgIdAndTxt("Spglib:invalidMesh",
                          "Mesh must be an array of 3 elements.");
    }
    int *mesh_ptr = static_cast<int *>(mxGetData(prhs[1]));
    int mesh[3] = {mesh_ptr[0], mesh_ptr[1], mesh_ptr[2]};

    // Extract and validate the rec_lattice argument
    if (mxGetM(prhs[2]) != 3 || mxGetN(prhs[2]) != 3) {
        mexErrMsgIdAndTxt("Spglib:invalidRecLattice",
                          "Rec_lattice must be a 3x3 matrix.");
    }
    double *rec_lattice_ptr = mxGetPr(prhs[2]);
    double rec_lattice[3][3];
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            rec_lattice[i][j] = rec_lattice_ptr[i + 3 * j];
        }
    }

    // Extract and validate the is_shift argument
    if (mxGetNumberOfElements(prhs[3]) != 3) {
        mexErrMsgIdAndTxt("Spglib:invalidShift",
                          "is_shift must be an array of 3 elements.");
    }
    int *is_shift_ptr = static_cast<int *>(mxGetData(prhs[3]));
    int is_shift[3] = {is_shift_ptr[0], is_shift_ptr[1], is_shift_ptr[2]};

    // Create the output bz_grid_address array (prod(mesh + 1) x 3 int array)
    size_t bz_grid_address_size = (mesh[0] + 1) * (mesh[1] + 1) * (mesh[2] + 1);
    plhs[0] =
        mxCreateNumericMatrix(bz_grid_address_size, 3, mxINT32_CLASS, mxREAL);
    int (*bz_grid_address)[3] = static_cast<int (*)[3]>(mxGetData(plhs[0]));

    // Create the output bz_map array (prod(mesh * 2) int array)
    size_t bz_map_size = mesh[0] * 2 * mesh[1] * 2 * mesh[2] * 2;
    plhs[1] = mxCreateNumericMatrix(bz_map_size, 1, mxINT32_CLASS, mxREAL);
    int *bz_map = static_cast<int *>(mxGetData(plhs[1]));

    // Call spg_relocate_BZ_grid_address
    int num_ir_grid_points = spg_relocate_BZ_grid_address(
        bz_grid_address, bz_map, grid_address, mesh, rec_lattice, is_shift);

    // Create the output num_ir_grid_points scalar
    plhs[2] = mxCreateDoubleScalar(static_cast<double>(num_ir_grid_points));
}

// [bz_grid_address, bz_map, num_ir_grid_points] =
// symspg('spg_relocate_dense_BZ_grid_address', grid_address, mesh, rec_lattice,
// is_shift)
void SpglibFunctions::spg_relocate_dense_BZ_grid_address_mex(
    int nlhs, mxArray *plhs[], int nrhs, mxArray const *prhs[]) {
    /*
     size_t spg_relocate_dense_BZ_grid_address(int bz_grid_address[][3], size_t
     bz_map[], int const grid_address[][3], int const mesh[3], double const
     rec_lattice[3][3], int const is_shift[3]);
    */

    // Validate the number of input arguments
    int const expected_number_of_inputs = 4;
    if (nrhs != expected_number_of_inputs) {
        mexErrMsgIdAndTxt("Spglib:invalidNumInputs",
                          "Incorrect number of inputs for "
                          "spg_relocate_dense_BZ_grid_address.");
    }

    // Extract and validate the grid_address argument
    mwSize num_grid_points = mxGetM(prhs[0]);
    if (mxGetN(prhs[0]) != 3) {
        mexErrMsgIdAndTxt("Spglib:invalidGridAddress",
                          "Grid_address must be an Nx3 array.");
    }
    int *grid_address_ptr = static_cast<int *>(mxGetData(prhs[0]));
    mexutil::Buffer2D<int, 3> grid_address(num_grid_points);
    for (mwSize i = 0; i < num_grid_points; ++i) {
        for (int j = 0; j < 3; ++j) {
            grid_address[i][j] = grid_address_ptr[i + j * num_grid_points];
        }
    }

    // Extract and validate the mesh argument
    if (mxGetNumberOfElements(prhs[1]) != 3) {
        mexErrMsgIdAndTxt("Spglib:invalidMesh",
                          "Mesh must be an array of 3 elements.");
    }
    int *mesh_ptr = static_cast<int *>(mxGetData(prhs[1]));
    int mesh[3] = {mesh_ptr[0], mesh_ptr[1], mesh_ptr[2]};

    // Extract and validate the rec_lattice argument
    if (mxGetM(prhs[2]) != 3 || mxGetN(prhs[2]) != 3) {
        mexErrMsgIdAndTxt("Spglib:invalidRecLattice",
                          "Rec_lattice must be a 3x3 matrix.");
    }
    double *rec_lattice_ptr = mxGetPr(prhs[2]);
    double rec_lattice[3][3];
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            rec_lattice[i][j] = rec_lattice_ptr[i + 3 * j];
        }
    }

    // Extract and validate the is_shift argument
    if (mxGetNumberOfElements(prhs[3]) != 3) {
        mexErrMsgIdAndTxt("Spglib:invalidShift",
                          "is_shift must be an array of 3 elements.");
    }
    int *is_shift_ptr = static_cast<int *>(mxGetData(prhs[3]));
    int is_shift[3] = {is_shift_ptr[0], is_shift_ptr[1], is_shift_ptr[2]};

    // Create the output bz_grid_address array (prod(mesh + 1) x 3 int array)
    size_t bz_grid_address_size = (mesh[0] + 1) * (mesh[1] + 1) * (mesh[2] + 1);
    plhs[0] =
        mxCreateNumericMatrix(bz_grid_address_size, 3, mxINT32_CLASS, mxREAL);
    int (*bz_grid_address)[3] = static_cast<int (*)[3]>(mxGetData(plhs[0]));

    // Create the output bz_map array (prod(mesh * 2) size_t array)
    size_t bz_map_size = mesh[0] * 2 * mesh[1] * 2 * mesh[2] * 2;
    plhs[1] = mxCreateNumericMatrix(bz_map_size, 1, mxUINT64_CLASS, mxREAL);
    size_t *bz_map = static_cast<size_t *>(mxGetData(plhs[1]));

    // Call spg_relocate_dense_BZ_grid_address
    size_t num_ir_grid_points = spg_relocate_dense_BZ_grid_address(
        bz_grid_address, bz_map, grid_address, mesh, rec_lattice, is_shift);

    // Create the output num_ir_grid_points scalar
    plhs[2] = mxCreateDoubleScalar(static_cast<double>(num_ir_grid_points));
}

// [lattice, success] = symspg('spg_niggli_reduce', lattice, symprec)
void SpglibFunctions::spg_niggli_reduce_mex(int nlhs, mxArray *plhs[], int nrhs,
                                            mxArray const *prhs[]) {
    /*
     int spg_niggli_reduce(double lattice[3][3], double const symprec);
    */

    // Validate the number of input arguments
    int const expected_number_of_inputs = 2;
    if (nrhs != expected_number_of_inputs) {
        mexErrMsgIdAndTxt("Spglib:invalidNumInputs",
                          "Incorrect number of inputs for spg_niggli_reduce.");
    }

    // Extract and validate the lattice argument
    if (mxGetM(prhs[0]) != 3 || mxGetN(prhs[0]) != 3) {
        mexErrMsgIdAndTxt("Spglib:invalidLattice",
                          "Lattice must be a 3x3 matrix.");
    }
    double *lattice_ptr = mxGetPr(prhs[0]);
    double lattice[3][3];
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            lattice[i][j] = lattice_ptr[i + 3 * j];
        }
    }

    // Extract the symprec argument
    double symprec = mxGetScalar(prhs[1]);

    // Call spg_niggli_reduce
    int success = spg_niggli_reduce(lattice, symprec);

    // Create the output lattice array (3x3 double array)
    plhs[0] = mxCreateDoubleMatrix(3, 3, mxREAL);
    double *lattice_out = mxGetPr(plhs[0]);
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            lattice_out[i + 3 * j] = lattice[i][j];
        }
    }

    // Create the output success scalar
    plhs[1] = mxCreateDoubleScalar(static_cast<double>(success));
}
