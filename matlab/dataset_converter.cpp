#include "dataset_converter.hpp"

#include "mex_array.hpp"
#include "mex_buffer.hpp"

namespace mexutil {

void DatasetDeleter::operator()(SpglibDataset* dataset) const noexcept {
    spg_free_dataset(dataset);
}

void MagneticDatasetDeleter::operator()(
    SpglibMagneticDataset* dataset) const noexcept {
    spg_free_magnetic_dataset(dataset);
}

namespace {

mxArray* makeDoubleVector(double const* values, mwSize size) {
    mxArray* array = mxCreateDoubleMatrix(size, 1, mxREAL);
    double* output = mxGetPr(array);
    for (mwSize index = 0; index < size; ++index) {
        output[index] = values[index];
    }
    return array;
}

}  // namespace

mxArray* makeDatasetStruct(SpglibDataset const& dataset) {
    char const* field_names[] = {
        "spacegroup_number",
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
        "pointgroup_symbol",
    };
    auto const field_count =
        static_cast<int>(sizeof(field_names) / sizeof(field_names[0]));
    mxArray* result = mxCreateStructMatrix(1, 1, field_count, field_names);

    setScalarField(result, 0, "spacegroup_number", dataset.spacegroup_number);
    setScalarField(result, 0, "hall_number", dataset.hall_number);
    setStringField(result, 0, "international_symbol",
                   dataset.international_symbol);
    setStringField(result, 0, "hall_symbol", dataset.hall_symbol);
    setStringField(result, 0, "choice", dataset.choice);
    setDoubleMatrixField(result, 0, "transformation_matrix",
                         dataset.transformation_matrix, 3, 3);
    mxSetField(result, 0, "origin_shift",
               makeDoubleVector(dataset.origin_shift, 3));
    setScalarField(result, 0, "n_operations", dataset.n_operations);
    set3DIntArrayField(result, 0, "rotations", dataset.rotations,
                       dataset.n_operations);
    setDouble2DArrayField(result, 0, "translations", dataset.translations,
                          dataset.n_operations);
    setScalarField(result, 0, "n_atoms", dataset.n_atoms);
    setIntArrayField(result, 0, "wyckoffs", dataset.wyckoffs, dataset.n_atoms);

    Buffer1D<char const*> site_symmetry_symbols(dataset.n_atoms);
    for (int index = 0; index < dataset.n_atoms; ++index) {
        site_symmetry_symbols[index] = dataset.site_symmetry_symbols[index];
    }
    mxSetField(
        result, 0, "site_symmetry_symbols",
        mxCreateCharMatrixFromStrings(dataset.n_atoms, site_symmetry_symbols));

    setIntArrayField(result, 0, "equivalent_atoms", dataset.equivalent_atoms,
                     dataset.n_atoms);
    setIntArrayField(result, 0, "crystallographic_orbits",
                     dataset.crystallographic_orbits, dataset.n_atoms);
    setDoubleMatrixField(result, 0, "primitive_lattice",
                         dataset.primitive_lattice, 3, 3);
    setIntArrayField(result, 0, "mapping_to_primitive",
                     dataset.mapping_to_primitive, dataset.n_atoms);
    setScalarField(result, 0, "n_std_atoms", dataset.n_std_atoms);
    setDoubleMatrixField(result, 0, "std_lattice", dataset.std_lattice, 3, 3);
    setIntArrayField(result, 0, "std_types", dataset.std_types,
                     dataset.n_std_atoms);
    setDouble2DArrayField(result, 0, "std_positions", dataset.std_positions,
                          dataset.n_std_atoms);
    setDoubleMatrixField(result, 0, "std_rotation_matrix",
                         dataset.std_rotation_matrix, 3, 3);
    setIntArrayField(result, 0, "std_mapping_to_primitive",
                     dataset.std_mapping_to_primitive, dataset.n_std_atoms);
    setStringField(result, 0, "pointgroup_symbol", dataset.pointgroup_symbol);
    return result;
}

mxArray* makeMagneticDatasetStruct(SpglibMagneticDataset const& dataset) {
    char const* field_names[] = {
        "uni_number",
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
        "std_rotation_matrix",
    };
    auto const field_count =
        static_cast<int>(sizeof(field_names) / sizeof(field_names[0]));
    mxArray* result = mxCreateStructMatrix(1, 1, field_count, field_names);

    setScalarField(result, 0, "uni_number", dataset.uni_number);
    setScalarField(result, 0, "msg_type", dataset.msg_type);
    setScalarField(result, 0, "hall_number", dataset.hall_number);
    setScalarField(result, 0, "tensor_rank", dataset.tensor_rank);
    setScalarField(result, 0, "n_operations", dataset.n_operations);
    set3DIntArrayField(result, 0, "rotations", dataset.rotations,
                       dataset.n_operations);
    setDouble2DArrayField(result, 0, "translations", dataset.translations,
                          dataset.n_operations);
    setIntArrayField(result, 0, "time_reversals", dataset.time_reversals,
                     dataset.n_operations);
    setScalarField(result, 0, "n_atoms", dataset.n_atoms);
    setIntArrayField(result, 0, "equivalent_atoms", dataset.equivalent_atoms,
                     dataset.n_atoms);
    setDoubleMatrixField(result, 0, "transformation_matrix",
                         dataset.transformation_matrix, 3, 3);
    mxSetField(result, 0, "origin_shift",
               makeDoubleVector(dataset.origin_shift, 3));
    setScalarField(result, 0, "n_std_atoms", dataset.n_std_atoms);
    setDoubleMatrixField(result, 0, "std_lattice", dataset.std_lattice, 3, 3);
    setIntArrayField(result, 0, "std_types", dataset.std_types,
                     dataset.n_std_atoms);
    setDouble2DArrayField(result, 0, "std_positions", dataset.std_positions,
                          dataset.n_std_atoms);

    // spglib supports scalar (rank 0) and vector (rank 1) site tensors.
    int const components_per_tensor = dataset.tensor_rank == 0 ? 1 : 3;
    int const tensor_elements = dataset.n_std_atoms * components_per_tensor;
    mxSetField(result, 0, "std_tensors",
               makeDoubleVector(dataset.std_tensors, tensor_elements));
    setDoubleMatrixField(result, 0, "std_rotation_matrix",
                         dataset.std_rotation_matrix, 3, 3);
    return result;
}

}  // namespace mexutil
