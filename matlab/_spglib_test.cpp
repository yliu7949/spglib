#include <stdio.h>
#include <cstdlib>
#include "spglib.h"

static void show_spg_dataset(double lattice[3][3], double const origin_shift[3],
                             double position[][3], int const num_atom,
                             int const types[]);
void show_cell(double const lattice[3][3], double const positions[][3],
               int const types[], int const num_atoms);

int main() {
    double lattice[3][3] = {
        {1.7807733239034182, -3.0843898737640294, 0.0000000000000000},
        {1.7807733239034182, 3.0843898737640294, 0.0000000000000000},
        {0.0000000000000000, 0.0000000000000000, 14.8472256999999992}};
    double position[][3] = {
        {0.3333333333333333, 0.6666666666666666, 0.8708151000000000},
        {0.6666666666666667, 0.3333333333333334, 0.1291849000000000},
        {0.6666666666666667, 0.3333333333333333, 0.3708151000000000},
        {0.3333333333333333, 0.6666666666666667, 0.6291849000000000},
        {0.3333333333333333, 0.6666666666666666, 0.2500000000000000},
        {0.6666666666666667, 0.3333333333333334, 0.7500000000000000}};
    int types[] = {52, 52, 52, 52, 74, 74};
    int num_atom = 6;
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

    show_cell(lattice, position, types, num_atom);
    double origin_shift[] = {0.0, 0.0, 0.0};
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
