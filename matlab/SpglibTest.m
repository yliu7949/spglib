classdef SpglibTest < matlab.unittest.TestCase

    methods (Test)
        function getVersionTest(~)
            version = spglib.Spglib.getVersion();
            disp(version);
        end

        function getVersionFullTest(~)
            versionFull = spglib.Spglib.getVersion("full");
            disp(versionFull);
        end

        function getMajorVersionTest(~)
            version = spglib.Spglib.getVersion("major");
            disp(version);
        end

        function getMinorVersionTest(~)
            version = spglib.Spglib.getVersion("minor");
            disp(version);
        end

        function getMicroVersionTest(~)
            version = spglib.Spglib.getVersion("micro");
            disp(version);
        end

        function getCommitTest(~)
            commit = spglib.Spglib.getCommit();
            disp(commit);
        end

        function getErrorCodeTest(~)
            disp(spglib.Spglib.getErrorCode());
        end

        function getErrorMessageTest(testCase)
            error_message = spglib.Spglib.getErrorMessage(4);
            testCase.assertEqual(error_message, 'too close distance between atoms');
        end

        function getDatasetTest1(testCase)
            lattice = [4, 0, 0; 0, 4, 0; 0, 0, 3];
            origin_shift = [0.1, 0.1, 0];
            position = [
                0,   0,    0;
                0.5, 0.5, 0.25;
                0.3, 0.3,    0;
                0.7, 0.7,    0;
                0.2, 0.8, 0.25;
                0.8, 0.2, 0.25;
                0,   0, 0.50;
                0.5, 0.5, 0.75;
                0.3, 0.3, 0.50;
                0.7, 0.7, 0.50;
                0.2, 0.8, 0.75;
                0.8, 0.2, 0.75
                ];
            types = [1, 1, 2, 2, 2, 2, 1, 1, 2, 2, 2, 2];
            num_atom = 12;
            symprec = 1e-5;

            for i = 1:num_atom
                position(i, :) = position(i, :) + origin_shift;
            end

            dataset = spglib.Spglib.getDataset(lattice, position, types, num_atom, symprec);
            testCase.assertEqual(dataset.spacegroup_number, 136);
        end

        function getDatasetTest2(testCase)
            lattice = [3.111, 0, 0; -1.5555, 2.6942050311733885, 0; 0, 0, 4.988];
            position = [
                1.0 / 3, 2.0 / 3, 0.0;
                2.0 / 3, 1.0 / 3, 0.5;
                1.0 / 3, 2.0 / 3, 0.6181;
                2.0 / 3, 1.0 / 3, 0.1181
                ];
            types = [1, 1, 2, 2];
            num_atom = 4;
            symprec = 1e-5;
            angle_tolerance = -1;

            dataset = spglib.Spglib.getDataset(lattice, position, types, num_atom, symprec, angle_tolerance);
            testCase.assertEqual(dataset.spacegroup_number, 186);
        end

        function getMagneticDatasetTest1(testCase)
            lattice = [5, 0, 0; 0, 5, 0; 0, 0, 3];
            position = [
                0,   0, 0;
                0.5, 0.5, 0.5;
                0.3, 0.3, 0;
                0.7, 0.7, 0;
                0.2, 0.8, 0.5;
                0.8, 0.2, 0.5
                ];
            types = [1, 1, 2, 2, 2, 2];
            spins = [0.3, 0.3, 0, 0, 0, 0];
            num_atom = 6;
            symprec = 1e-5;

            dataset = spglib.Spglib.getMagneticDataset(lattice, position, types, spins, 0, num_atom, 0, symprec);
            testCase.assertEqual(dataset.msg_type, 1);
            testCase.assertEqual(dataset.uni_number, 1155);
        end

        function getMagneticDatasetTest2(testCase)
            lattice = [
                0.00000000, 0.00000000, -4.05000000;
                -5.00000000,  0.00000000,  0.00000000;
                -2.50000000,  4.33012702,  0.00000000
                ];
            positions = [
                0.50000000, 0.33333333, 0.33333333;
                0.50000000, 0.66666667, 0.66666667;
                0.00000000, 0.00000000, 0.00000000;
                0.00000000, 0.00000000, 0.50000000;
                0.00000000, 0.50000000, 0.50000000;
                0.00000000, 0.50000000, 0.00000000
                ];
            types = [1, 1, 1, 2, 2, 2];
            tensors = [
                -0.00200000,  0.00200000,  1.90000000,  0.00200000,  0.00000000;
                1.91100000, -0.00100000,  0.00200000, -2.23300000, -0.00000000;
                -0.00000000, -0.06000000,  0.00000000,  0.00000000, -0.03200000;
                0.00000000, -0.00000000, -0.02900000,  0.00000000,  0.00000000;
                ];

            tensor_rank = 1;
            num_atoms = 6;
            is_axial = true;
            symprec = 1e-5;
            angle_tolerance = -1;
            mag_symprec = 1e-3;

            dataset = spglib.Spglib.getMagneticDataset(lattice, positions, types, tensors, tensor_rank, num_atoms, is_axial, symprec, angle_tolerance, mag_symprec);
            testCase.assertEqual(dataset.msg_type, 1);
        end

        function getDatasetWithHallNumberTest1(testCase)
            lattice = [3.111, 0, 0; -1.5555, 2.6942050311733885, 0; 0, 0, 4.988];
            position = [
                1.0 / 3, 2.0 / 3, 0.0;
                2.0 / 3, 1.0 / 3, 0.5;
                1.0 / 3, 2.0 / 3, 0.6181;
                2.0 / 3, 1.0 / 3, 0.1181
                ];
            types = [1, 1, 2, 2];
            num_atom = 4;
            hall_number = 480;
            symprec = 1e-5;

            dataset = spglib.Spglib.getDatasetWithHallNumber(lattice, position, types, num_atom, hall_number, symprec);
            testCase.assertEqual(dataset.spacegroup_number, 186);
        end

        function getDatasetWithHallNumberTest2(testCase)
            lattice = [3.111, 0, 0; -1.5555, 2.6942050311733885, 0; 0, 0, 4.988];
            position = [
                1.0 / 3, 2.0 / 3, 0.0;
                2.0 / 3, 1.0 / 3, 0.5;
                1.0 / 3, 2.0 / 3, 0.6181;
                2.0 / 3, 1.0 / 3, 0.1181
                ];
            types = [1, 1, 2, 2];
            num_atom = 4;
            hall_number = 480;
            symprec = 1e-5;
            angle_tolerance = -1;

            dataset = spglib.Spglib.getDatasetWithHallNumber(lattice, position, types, num_atom, hall_number, symprec, angle_tolerance);
            testCase.assertEqual(dataset.spacegroup_number, 186);
        end

        function getSymmetryWithCollinearSpinTest1(testCase)
            lattice = [
                4, 0, 0;
                0, 4, 0;
                0, 0, 4
                ];

            position = [
                0, 0, 0;
                0.5, 0.5, 0.5
                ];

            types = [1, 1];
            spins = [1, 1];
            num_atom = 2;
            max_size = 300;

            [rotations, translations, equivalent_atoms, num_operations] = spglib.Spglib.getSymmetryWithCollinearSpin(max_size, lattice, position, types, spins, num_atom, 1e-5);

            testCase.assertEqual([size(rotations)], [96, 3, 3]);
            testCase.assertEqual([size(translations)], [96, 3]);
            testCase.assertEqual(length(equivalent_atoms), 2);
            testCase.assertEqual(num_operations, 96);
        end

        function getSymmetryWithCollinearSpinTest2(testCase)
            lattice = [
                4, 0, 0;
                0, 4, 0;
                0, 0, 4
                ];

            position = [
                0, 0, 0;
                0.5, 0.5, 0.5
                ];

            types = [1, 1];
            spins = [1, 1];
            num_atom = 2;
            max_size = 300;
            angle_tolerance = -1;

            [rotations, translations, equivalent_atoms, num_operations] = spglib.Spglib.getSymmetryWithCollinearSpin(max_size, lattice, position, types, spins, num_atom, 1e-5, angle_tolerance);

            testCase.assertEqual([size(rotations)], [96, 3, 3]);
            testCase.assertEqual([size(translations)], [96, 3]);
            testCase.assertEqual(length(equivalent_atoms), 2);
            testCase.assertEqual(num_operations, 96);
        end

        function getSymmetryWithCollinearSpinTest3(testCase)
            lattice = [
                4, 0, 0;
                0, 4, 0;
                0, 0, 4
                ];

            position = [
                0, 0, 0;
                0.5, 0.5, 0.5
                ];

            types = [1, 1];
            spins = [1, 1];
            num_atom = 2;
            max_size = 300;
            angle_tolerance = -1;
            mag_symprec = 1e-5;

            [rotations, translations, equivalent_atoms, num_operations] = spglib.Spglib.getSymmetryWithCollinearSpin(max_size, lattice, position, types, spins, num_atom, 1e-5, angle_tolerance, mag_symprec);

            testCase.assertEqual([size(rotations)], [96, 3, 3]);
            testCase.assertEqual([size(translations)], [96, 3]);
            testCase.assertEqual(length(equivalent_atoms), 2);
            testCase.assertEqual(num_operations, 96);
        end

        function getSymmetryWithSiteTensorsTest(testCase)
            lattice = [
                5.7461, 0, 0;
                0, 7.6637, 0;
                0, 0, 5.5333
                ];
            position = [
                0.051300, 0.250000, 0.990500;
                0.948700, 0.750000, 0.009500;
                0.551300, 0.250000, 0.509500;
                0.448700, 0.750000, 0.490500;
                0.000000, 0.000000, 0.500000;
                0.000000, 0.500000, 0.500000;
                0.500000, 0.500000, 0.000000;
                0.500000, 0.000000, 0.000000;
                0.484900, 0.250000, 0.077700;
                0.515100, 0.750000, 0.922300;
                0.984900, 0.250000, 0.422300;
                0.015100, 0.750000, 0.577700;
                0.308500, 0.040800, 0.722700;
                0.691500, 0.540800, 0.277300;
                0.691500, 0.959200, 0.277300;
                0.308500, 0.459200, 0.722700;
                0.808500, 0.459200, 0.777300;
                0.191500, 0.959200, 0.222700;
                0.191500, 0.540800, 0.222700;
                0.808500, 0.040800, 0.777300
                ];
            types = [1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3];
            tensors = [
                0, 0, 0;
                0, 0, 0;
                0, 0, 0;
                0, 0, 0;
                3.87, 0, 0;
                -3.87, 0, 0;
                -3.87, 0, 0;
                3.87, 0, 0;
                0, 0, 0;
                0, 0, 0;
                0, 0, 0;
                0, 0, 0;
                0, 0, 0;
                0, 0, 0;
                0, 0, 0;
                0, 0, 0;
                0, 0, 0;
                0, 0, 0;
                0, 0, 0;
                0, 0, 0
                ];
            num_atom = 20;
            max_size = num_atom * 96;
            symprec = 1e-5;

            [~, ~, ~, ~, ~, num_operations]  = spglib.Spglib.getSymmetryWithSiteTensors(max_size, lattice, position, types, tensors, 1, num_atom, true, true, symprec);
            testCase.assertEqual(num_operations, 8);
        end

        function getSpacegroupTypeFromSymmetryTest(testCase)
            lattice = [4, 0, 0; 0, 4, 0; 0, 0, 4];
            position = [
                0, 0, 0;
                0.5, 0.5, 0.5
                ];
            types = [1, 1];
            num_atom = 2;
            symprec = 1e-5;

            dataset = spglib.Spglib.getDataset(lattice, position, types, num_atom, symprec);
            spacegroupType = spglib.Spglib.getSpacegroupTypeFromSymmetry(dataset.rotations, dataset.translations, dataset.n_operations, lattice, symprec);
            testCase.assertEqual(spacegroupType.number, dataset.spacegroup_number);
            testCase.assertEqual(spacegroupType.number, 229);
        end

        function getMagneticSpacegroupTypeFromSymmetryTest(testCase)
            lattice = [
                5.7461, 0, 0;
                0, 7.6637, 0;
                0, 0, 5.5333
                ];
            position = [
                0.051300, 0.250000, 0.990500;
                0.948700, 0.750000, 0.009500;
                0.551300, 0.250000, 0.509500;
                0.448700, 0.750000, 0.490500;
                0.000000, 0.000000, 0.500000;
                0.000000, 0.500000, 0.500000;
                0.500000, 0.500000, 0.000000;
                0.500000, 0.000000, 0.000000;
                0.484900, 0.250000, 0.077700;
                0.515100, 0.750000, 0.922300;
                0.984900, 0.250000, 0.422300;
                0.015100, 0.750000, 0.577700;
                0.308500, 0.040800, 0.722700;
                0.691500, 0.540800, 0.277300;
                0.691500, 0.959200, 0.277300;
                0.308500, 0.459200, 0.722700;
                0.808500, 0.459200, 0.777300;
                0.191500, 0.959200, 0.222700;
                0.191500, 0.540800, 0.222700;
                0.808500, 0.040800, 0.777300
                ];
            types = [1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3];
            tensors = [
                0, 0, 0;
                0, 0, 0;
                0, 0, 0;
                0, 0, 0;
                3.87, 0, 0;
                -3.87, 0, 0;
                -3.87, 0, 0;
                3.87, 0, 0;
                0, 0, 0;
                0, 0, 0;
                0, 0, 0;
                0, 0, 0;
                0, 0, 0;
                0, 0, 0;
                0, 0, 0;
                0, 0, 0;
                0, 0, 0;
                0, 0, 0;
                0, 0, 0;
                0, 0, 0
                ];
            num_atom = 20;
            max_size = num_atom * 96;
            symprec = 1e-5;

            [rotations, translations, ~, ~, spin_flips, num_operations]  = spglib.Spglib.getSymmetryWithSiteTensors(max_size, lattice, position, types, tensors, 1, num_atom, true, true, symprec);
            testCase.assertEqual(num_operations, 8);

            time_reversals = idivide(1 - spin_flips, 2);
            spacegroupType = spglib.Spglib.getMagneticSpacegroupTypeFromSymmetry(rotations, translations, time_reversals, num_operations, lattice, symprec);
            testCase.assertEqual(spacegroupType.uni_number, 546);
        end

        function getPointgroupTest(testCase)
            lattice = [4, 0, 0; 0, 4, 0; 0, 0, 4];
            position = [
                0, 0, 0;
                0.5, 0.5, 0.5
                ];
            types = [1, 1];
            num_atom = 2;
            symprec = 1e-5;

            dataset = spglib.Spglib.getDataset(lattice, position, types, num_atom, symprec);
            [symbol, trans_mat, ~] = spglib.Spglib.getPointgroup(dataset.rotations, dataset.n_operations);
            testCase.assertEqual(symbol, 'm-3m');
            testCase.assertEqual(trans_mat, int32(eye(3)));
        end

        function getSymmetryFromDatabaseTest(~)
            spglib.Spglib.getSymmetryFromDatabase(460);
        end

        function getMagneticSymmetryFromDatabaseTest(testCase)
            rotations = spglib.Spglib.getMagneticSymmetryFromDatabase(1242, 434);
            testCase.assertTrue(size(rotations, 1) > 0);
        end

        function getSpacegroupTypeTest(testCase)
            spacegroup = spglib.Spglib.getSpacegroupType(446);
            testCase.assertEqual(spacegroup.number, 156);
        end

        function getMagneticSpacegroupType(testCase)
            spacegroup = spglib.Spglib.getMagneticSpacegroupType(1279);
            testCase.assertEqual(spacegroup.number, 156);
        end

        function standardizeCellTest(testCase)
            lattice = [
                4.8076344022756095, 0, 0;
                -2.4038172011378047, 4.1635335244786962, 0;
                0, 0, 13.1172699198127543
                ];

            position = [
                0.0000000000000000, 0.0000000000000000, 0.3521850942289043;
                0.6666666666666643, 0.3333333333333357, 0.6855184275622400;
                0.3333333333333357, 0.6666666666666643, 0.0188517608955686;
                0.0000000000000000, 0.0000000000000000, 0.6478149057711028;
                0.6666666666666643, 0.3333333333333357, 0.9811482391044314;
                0.3333333333333357, 0.6666666666666643, 0.3144815724377600;
                0.0000000000000000, 0.0000000000000000, 0.1478149057710957;
                0.6666666666666643, 0.3333333333333357, 0.4811482391044314;
                0.3333333333333357, 0.6666666666666643, 0.8144815724377600;
                0.0000000000000000, 0.0000000000000000, 0.8521850942288972;
                0.6666666666666643, 0.3333333333333357, 0.1855184275622400;
                0.3333333333333357, 0.6666666666666643, 0.5188517608955686;
                0.3061673906454899, 0.0000000000000000, 0.2500000000000000;
                0.9728340573121541, 0.3333333333333357, 0.5833333333333357;
                0.6395007239788255, 0.6666666666666643, 0.9166666666666643;
                0.6938326093545102, 0.0000000000000000, 0.7500000000000000;
                0.3604992760211744, 0.3333333333333357, 0.0833333333333357;
                0.0271659426878458, 0.6666666666666643, 0.4166666666666643;
                0.0000000000000000, 0.3061673906454899, 0.2500000000000000;
                0.6666666666666643, 0.6395007239788255, 0.5833333333333357;
                0.3333333333333357, 0.9728340573121541, 0.9166666666666643;
                0.0000000000000000, 0.6938326093545102, 0.7500000000000000;
                0.6666666666666643, 0.0271659426878458, 0.0833333333333357;
                0.3333333333333357, 0.3604992760211744, 0.4166666666666643;
                0.6938326093545102, 0.6938326093545102, 0.2500000000000000;
                0.3604992760211744, 0.0271659426878458, 0.5833333333333357;
                0.0271659426878458, 0.3604992760211744, 0.9166666666666643;
                0.3061673906454899, 0.3061673906454899, 0.7500000000000000;
                0.9728340573121541, 0.6395007239788255, 0.0833333333333357;
                0.6395007239788255, 0.9728340573121541, 0.4166666666666643
                ];

            types = [ones(1, 12), 2 * ones(1, 18)];
            num_atom = 30;
            to_primitive = 1;
            no_idealize = 1;
            symprec = 1e-5;

            % 调用 Spglib 函数
            [~, ~, ~, num_primitive_atom] = spglib.Spglib.standardizeCell( ...
                lattice, position, types, num_atom, to_primitive, no_idealize, symprec);
            testCase.assertEqual(num_primitive_atom, 10);
        end

        function findPrimitiveCellTest1(~)
            lattice = [
                4, 0, 0;
                0, 4, 0;
                0, 0, 4
                ];

            position = [
                0, 0, 0;
                0.5, 0.5, 0.5;
                0.5, 0.5, 0.5
                ];

            types = [1, 1, 1];
            num_atom = 3;
            symprec = 1e-5;

            [~, ~, ~, num_primitive_atom] = spglib.Spglib.findPrimitive(lattice, position, types, num_atom, symprec);
            if num_primitive_atom == 0
                disp(spglib.Spglib.getErrorMessage(spglib.Spglib.getErrorCode()));
            end
        end

        function findPrimitiveCellTest2(testCase)
            lattice = 4 * eye(3);
            position = [0, 0, 0; 0.5, 0.5, 0.5];
            types = [1, 1];
            num_atom = 2;
            symprec = 1e-5;

            [~, ~, ~, num_primitive_atom] = spglib.Spglib.findPrimitive(lattice, position, types, num_atom, symprec);
            testCase.assertEqual(num_primitive_atom, 1);
        end

        function findPrimitiveCellTest3(testCase)
            lattice = 4 * eye(3);
            position = [0, 0, 0; 0.5, 0.5, 0.5];
            types = [1, 1];
            num_atom = 2;
            symprec = 1e-5;
            angle_tolerance = 1e-5;

            [~, ~, ~, num_primitive_atom] = spglib.Spglib.findPrimitive(lattice, position, types, num_atom, symprec, angle_tolerance);
            testCase.assertEqual(num_primitive_atom, 1);
        end
    
        function refineCellTest(testCase)
            lattice = [0, 2, 2; 2, 0, 2; 2, 2, 0];
            position = [0, 0, 0];
            types = 1;
            num_atom = 1;
            symprec = 1e-5;

            [lattice, ~, ~, num_atom_bravais] = ...
                spglib.Spglib.refineCell(lattice, position, types, num_atom, symprec);
            testCase.assertEqual(lattice, 4 * eye(3));
            testCase.assertEqual(num_atom_bravais, 4);
        end
    
        function delaunayReduceTest(~)
            lattice = [
                3.0, 0.1, 0.2;
                0.1, 3.5, 0.3;
                0.2, 0.3, 4.0
                ];
            symprec = 1e-5;

            spglib.Spglib.delaunayReduce(lattice, symprec);
        end

        function getGridPointFromAddressTest(~)
            grid_address = [1, 2, 3];
            mesh = [10, 10, 10];

            grid_point_index = spglib.Spglib.getGridPointFromAddress(grid_address, mesh);
            disp(['Grid point index: ', num2str(grid_point_index)]);
        end

        function getDenseGridPointFromAddressTest(~)
            grid_address = [2, 3, 5];
            mesh = [10, 10, 10];

            dense_grid_point_index = spglib.Spglib.getDenseGridPointFromAddress(grid_address, mesh);
            disp(['Grid point index: ', num2str(dense_grid_point_index)]);
        end
    
        function getIrReciprocalMeshTest(testCase)
            lattice = [
                4, 0, 0;
                0, 4, 0;
                0, 0, 3
                ];

            position = [
                0, 0, 0;
                0.5, 0.5, 0.5;
                0.3, 0.3, 0;
                0.7, 0.7, 0;
                0.2, 0.8, 0.5;
                0.8, 0.2, 0.5
                ];

            types = [1, 1, 2, 2, 2, 2];
            num_atom = 6;
            mesh = [40, 40, 40];
            is_shift = [1, 1, 1];
            symprec = 1e-5;

            [~, ~, num_ir_kpoints] = spglib.Spglib.getIrReciprocalMesh(mesh, is_shift, 1, lattice, position, types, num_atom, symprec);
            testCase.assertEqual(num_ir_kpoints, 4200);
        end

        function getStabilizedReciprocalMeshTest(testCase)
            lattice = [
                4, 0, 0;
                0, 4, 0;
                0, 0, 3
                ];
            position = [
                0, 0, 0;
                0.5, 0.5, 0.5;
                0.3, 0.3, 0;
                0.7, 0.7, 0;
                0.2, 0.8, 0.5;
                0.8, 0.2, 0.5
                ];
            types = [1, 1, 2, 2, 2, 2];
            num_atom = 6;
            mesh = [40, 40, 40];
            is_shift = [1, 1, 1];
            is_time_reversal = 1;
            symprec = 1e-5;
            num_qpoints = 1;
            qpoints = [0, 0.5, 0.5];

            dataset = spglib.Spglib.getDataset(lattice, position, types, num_atom, symprec);
            testCase.assertNotEmpty(dataset);

            [~, ~, num_ir_kpoints] = spglib.Spglib.getStabilizedReciprocalMesh(mesh, is_shift, is_time_reversal, ...
                dataset.n_operations, dataset.rotations, num_qpoints, qpoints);
            testCase.assertEqual(num_ir_kpoints, 8000);
        end
    
        function getDenseStabilizedReciprocalMeshTest(testCase)
            mesh = [40, 40, 40];
            is_shift = [0, 0, 0];
            is_time_reversal = 1;
            rotations = zeros(1, 3, 3);
            rotations(1,:,:) = eye(3);
            num_rotations = size(rotations, 1);
            num_qpoints = 1;
            qpoints = [0, 0, 0];

            [~, ~, num_ir_kpoints] = spglib.Spglib.getDenseStabilizedReciprocalMesh(mesh, is_shift, is_time_reversal, ...
                num_rotations, rotations, num_qpoints, qpoints);

            testCase.assertTrue(num_ir_kpoints > 0);
        end
    
        function relocateBZGridAddressTest(testCase)
            rec_lattice = [
                -0.17573761, 0.17573761, 0.17573761;
                0.17573761, -0.17573761, 0.17573761;
                0.17573761, 0.17573761, -0.17573761
                ];

            rotations = zeros(1, 3, 3);
            rotations(1,:,:) = eye(3);
            mesh = [40, 40, 40];
            is_shift = [0, 0, 0];
            qpoints = [0, 0, 0];

            [grid_address, ~, num_ir] = ...
                spglib.Spglib.getStabilizedReciprocalMesh(mesh, is_shift, 1, 1, rotations, 1, qpoints);
            testCase.assertTrue(num_ir > 0);

            [~, ~, num_ir_grid_points] = ...
                spglib.Spglib.relocateBZGridAddress(grid_address, mesh, rec_lattice, is_shift);
            testCase.assertEqual(num_ir_grid_points, 65861);
        end

        function relocateDenseBZGridAddressTest(testCase)
            rec_lattice = [
                -0.17573761, 0.17573761, 0.17573761;
                0.17573761, -0.17573761, 0.17573761;
                0.17573761, 0.17573761, -0.17573761
                ];

            rotations = zeros(1, 3, 3);
            rotations(1,:,:) = eye(3);
            mesh = [40, 40, 40];
            is_shift = [0, 0, 0];
            qpoints = [0, 0, 0];

            [grid_address, ~, num_ir] = ...
                spglib.Spglib.getDenseStabilizedReciprocalMesh(mesh, is_shift, 1, 1, rotations, 1, qpoints);
            testCase.assertTrue(num_ir > 0);

            [~, ~, num_ir_grid_points] = ...
                spglib.Spglib.relocateDenseBZGridAddress(grid_address, mesh, rec_lattice, is_shift);
            testCase.assertEqual(num_ir_grid_points, 65861);
        end

        function niggliReduceTest(testCase)
            lattice = [2 0 1; 1 2 0; 0 1 2];
            symprec = 1e-5;

            [~, result] = spglib.Spglib.niggliReduce(lattice, symprec);
            testCase.assertEqual(result, 1);
        end
    end
end
