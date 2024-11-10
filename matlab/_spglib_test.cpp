#include <stdio.h>
#include <cstdlib>
#include "spglib.h"

static void show_spg_dataset(double lattice[3][3], double const origin_shift[3],
                             double position[][3], int const num_atom,
                             int const types[]);
void show_cell(double const lattice[3][3], double const positions[][3],
               int const types[], int const num_atoms);

int main() {
    /*
    mexPrintf("mesh: %d, %d, %d\n", mesh[0], mesh[1], mesh[2]);
    mexPrintf("is_shift: %d, %d, %d\n", is_shift[0], is_shift[1], is_shift[2]);
    mexPrintf("is_time_reversal: %d\n", is_time_reversal);
    mexPrintf("symprec: %.6f\n", symprec);
    show_cell(lattice, position, types, num_atom);
     */
    double lattice[3][3] = {{0, 2, 2}, {2, 0, 2}, {2, 2, 0}};

    /* 4 times larger memory space must be prepared. */
    double position[4][3];
    int types[4];

    int num_atom_bravais;
    int num_atom = 1;
    double symprec = 1e-5;

    position[0][0] = 0;
    position[0][1] = 0;
    position[0][2] = 0;
    types[0] = 1;

    show_cell(lattice, position, types, num_atom);

    /* lattice, position, and types are overwritten. */
    num_atom_bravais =
        spg_refine_cell(lattice, position, types, num_atom, symprec);

    printf("\nnum_atom_bravais: %d\n", num_atom_bravais);
    show_cell(lattice, position, types, num_atom_bravais);
}

int main4() {
    /* MAGNDATA #0.1: LaMnO3 */
    /* BNS: Pn'ma' (62.448), MHall: -P 2ac' 2n' (546) */
    int max_size, size, i;
    double lattice[][3] = {{5.7461, 0, 0}, {0, 7.6637, 0}, {0, 0, 5.5333}};
    /* clang-format off */
    double position[][3] = {
        {0.051300, 0.250000, 0.990500}, /* La */
        {0.948700, 0.750000, 0.009500},
        {0.551300, 0.250000, 0.509500},
        {0.448700, 0.750000, 0.490500},
        {0.000000, 0.000000, 0.500000}, /* Mn */
        {0.000000, 0.500000, 0.500000},
        {0.500000, 0.500000, 0.000000},
        {0.500000, 0.000000, 0.000000},
        {0.484900, 0.250000, 0.077700}, /* O1 */
        {0.515100, 0.750000, 0.922300},
        {0.984900, 0.250000, 0.422300},
        {0.015100, 0.750000, 0.577700},
        {0.308500, 0.040800, 0.722700}, /* O2 */
        {0.691500, 0.540800, 0.277300},
        {0.691500, 0.959200, 0.277300},
        {0.308500, 0.459200, 0.722700},
        {0.808500, 0.459200, 0.777300},
        {0.191500, 0.959200, 0.222700},
        {0.191500, 0.540800, 0.222700},
        {0.808500, 0.040800, 0.777300},
    };
    int types[] = {1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3};
    double tensors[] = {
        0, 0, 0, /* La */
        0, 0, 0,
        0, 0, 0,
        0, 0, 0,
        3.87, 0, 0, /* Mn */
        -3.87, 0, 0,
        -3.87, 0, 0,
        3.87, 0, 0,
        0, 0, 0, /* O1 */
        0, 0, 0,
        0, 0, 0,
        0, 0, 0,
        0, 0, 0, /* O2 */
        0, 0, 0,
        0, 0, 0,
        0, 0, 0,
        0, 0, 0,
        0, 0, 0,
        0, 0, 0,
        0, 0, 0,
    };
    /* clang-format on */
    int num_atom = 20;

    int equivalent_atoms[20];
    double primitive_lattice[3][3];
    int(*rotation)[3][3];
    double(*translation)[3];
    int *spin_flips;

    // 96  = 48 * 2
    // = (max number of order of point group) * (spin degrees of freedom)
    max_size = num_atom * 96;
    rotation = (int(*)[3][3])malloc(sizeof(int[3][3]) * max_size);
    translation = (double(*)[3])malloc(sizeof(double[3]) * max_size);
    spin_flips = (int *)malloc(sizeof(int *) * max_size);

    /* Find equivalent_atoms, primitive_lattice, spin_flips */
    size = spg_get_symmetry_with_site_tensors(
        rotation, translation, equivalent_atoms, primitive_lattice, spin_flips,
        max_size, lattice, position, types, tensors, 1, num_atom,
        1 /* with_time_reversal */, 1 /* is_axial */, 1e-5);

    printf("Size: %d\n", size);
    printf("Spin flips:\n");
    for (i = 0; i < size; i++) {
        printf("%d ", spin_flips[i]);
    }
    printf("\n");

    printf("time_reversals: \n");
    int *time_reversals;
    time_reversals = (int *)malloc(sizeof(int *) * size);
    for (i = 0; i < size; i++) {
        time_reversals[i] = (1 - spin_flips[i]) / 2;
        printf("%d ", time_reversals[i]);
    }
    printf("\n");

    SpglibMagneticSpacegroupType msgtype =
        spg_get_magnetic_spacegroup_type_from_symmetry(
            rotation, translation, time_reversals, size, lattice, 1e-5);
    printf("uni_number: %d\n", msgtype.uni_number);

    // Free allocated memory
    free(rotation);
    free(translation);
    free(spin_flips);
    free(time_reversals);

    return 0;
}

int main3() {
    double lattice[3][3] = {{4, 0, 0}, {0, 4, 0}, {0, 0, 4}};
    double position[][3] = {{0, 0, 0}, {0.5, 0.5, 0.5}};
    int types[] = {1, 1};
    int num_atom = 2;
    double symprec = 1e-5;

    SpglibSpacegroupType spgtype;
    SpglibDataset *dataset;

    dataset = NULL;
    dataset = spg_get_dataset(lattice, position, types, num_atom, symprec);

    printf("spacegroup_number = %d is found by spg_get_dataset.\n",
           dataset->spacegroup_number);
    spgtype = spg_get_spacegroup_type_from_symmetry(
        dataset->rotations, dataset->translations, dataset->n_operations,
        lattice, symprec);
    printf("number = %d is found by spg_get_spacegroup_type_from_symmetry.\n",
           spgtype.number);

    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            printf("%d ", dataset->rotations[0][i][j]);
        }
        printf("\n");
    }

    spg_free_dataset(dataset);

    return 0;
}

int main2() {
    printf("*** Example of spg_get_dataset (Rutile two unit cells) ***:\n");
    double lattice[3][3] = {{4, 0, 0}, {0, 4, 0}, {0, 0, 3}};
    double origin_shift[3] = {0.1, 0.1, 0};
    double position[][3] = {
        {0, 0, 0},        {0.5, 0.5, 0.25}, {0.3, 0.3, 0},    {0.7, 0.7, 0},
        {0.2, 0.8, 0.25}, {0.8, 0.2, 0.25}, {0, 0, 0.5},      {0.5, 0.5, 0.75},
        {0.3, 0.3, 0.5},  {0.7, 0.7, 0.5},  {0.2, 0.8, 0.75}, {0.8, 0.2, 0.75}};
    int types[] = {1, 1, 2, 2, 2, 2, 1, 1, 2, 2, 2, 2};
    int num_atom = 12;

    show_spg_dataset(lattice, origin_shift, position, num_atom, types);
    return 0;
}

static void show_spg_dataset(double lattice[3][3], double const origin_shift[3],
                             double position[][3], int const num_atom,
                             int const types[]) {
    SpglibDataset *dataset;
    char ptsymbol[6];
    int pt_trans_mat[3][3];

    char const *wl = "abcdefghijklmnopqrstuvwxyz";

    for (int i = 0; i < num_atom; i++) {
        for (int j = 0; j < 3; j++) {
            position[i][j] += origin_shift[j];
        }
    }

    dataset = spg_get_dataset(lattice, position, types, num_atom, 1e-5);

    printf("International: %s (%d)\n", dataset->international_symbol,
           dataset->spacegroup_number);
    printf("Hall symbol:   %s\n", dataset->hall_symbol);
    spg_get_pointgroup(ptsymbol, pt_trans_mat, dataset->rotations,
                       dataset->n_operations);
    printf("Point group:   %s\n", ptsymbol);
    printf("Transformation matrix:\n");
    for (int i = 0; i < 3; i++) {
        printf("%f %f %f\n", dataset->transformation_matrix[i][0],
               dataset->transformation_matrix[i][1],
               dataset->transformation_matrix[i][2]);
    }
    printf("Wyckoff letters:\n");
    for (int i = 0; i < dataset->n_atoms; i++) {
        printf("%c ", wl[dataset->wyckoffs[i]]);
    }
    printf("\n");
    printf("Equivalent atoms:\n");
    for (int i = 0; i < dataset->n_atoms; i++) {
        printf("%d ", dataset->equivalent_atoms[i]);
    }
    printf("\n");

    for (int i = 0; i < dataset->n_operations; i++) {
        printf("--- %d ---\n", i + 1);
        for (int j = 0; j < 3; j++) {
            printf("%2d %2d %2d\n", dataset->rotations[i][j][0],
                   dataset->rotations[i][j][1], dataset->rotations[i][j][2]);
        }
        printf("%f %f %f\n", dataset->translations[i][0],
               dataset->translations[i][1], dataset->translations[i][2]);
    }

    spg_free_dataset(dataset);
}

void show_matrix_3d(double const lattice[3][3]) {
    for (int i = 0; i < 3; i++) {
        printf("%f %f %f\n", lattice[0][i], lattice[1][i], lattice[2][i]);
    }
}

void show_cell(double const lattice[3][3], double const positions[][3],
               int const types[], int const num_atoms) {
    printf("Lattice parameter:\n");
    show_matrix_3d(lattice);
    printf("Atomic positions:\n");
    for (int i = 0; i < num_atoms; i++) {
        printf("%d: %f %f %f\n", types[i], positions[i][0], positions[i][1],
               positions[i][2]);
    }
}
