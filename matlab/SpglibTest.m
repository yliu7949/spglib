classdef SpglibTest < matlab.unittest.TestCase

    methods(TestClassSetup)
        function setupPath(~)
            spglibFolder = fullfile(fileparts(mfilename("fullpath")), "spglib");
            addpath(spglibFolder);
        end
    end

    methods (Test)
        function getVersionTest(~)
            version = Spglib.getVersion();
            disp(version);
        end

        function getVersionFullTest(~)
            versionFull = Spglib.getVersion("full");
            disp(versionFull);
        end

        function getMajorVersionTest(~)
            version = Spglib.getVersion("major");
            disp(version);
        end

        function getMinorVersionTest(~)
            version = Spglib.getVersion("minor");
            disp(version);
        end

        function getMicroVersionTest(~)
            version = Spglib.getVersion("micro");
            disp(version);
        end

        function getCommitTest(~)
            commit = Spglib.getCommit();
            disp(commit);
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

            dataset = Spglib.getDataset(lattice, position, types, num_atom, symprec);
            testCase.assertEqual(dataset.spacegroup_number, 136);
        end

        function getDatasetTest2(testCase)
            lattice = [3.111, -1.5555, 0; 0, 2.6942050311733885, 0; 0, 0, 4.988];
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

            dataset = Spglib.getDataset(lattice, position, types, num_atom, symprec, angle_tolerance);
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

            dataset = Spglib.getMagneticDataset(lattice, position, types, spins, 0, num_atom, 0, symprec);
            testCase.assertEqual(dataset.msg_type, 1);
            testCase.assertEqual(dataset.uni_number, 1155);
        end

        function getMagneticDatasetTest2(testCase)
            lattice = [
                0.00000000, -5.00000000, -2.50000000;
                0.00000000,  0.00000000,  4.33012702;
                -4.05000000,  0.00000000,  0.00000000
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

            dataset = Spglib.getMagneticDataset(lattice, positions, types, tensors, tensor_rank, num_atoms, is_axial, symprec, angle_tolerance, mag_symprec);
            testCase.assertEqual(dataset.msg_type, 1);
        end

        function getDatasetWithHallNumberTest1(testCase)
            lattice = [3.111, -1.5555, 0; 0, 2.6942050311733885, 0; 0, 0, 4.988];
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

            dataset = Spglib.getDatasetWithHallNumber(lattice, position, types, num_atom, hall_number, symprec);
            testCase.assertEqual(dataset.spacegroup_number, 186);
        end

        function getDatasetWithHallNumberTest2(testCase)
            lattice = [3.111, -1.5555, 0; 0, 2.6942050311733885, 0; 0, 0, 4.988];
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

            dataset = Spglib.getDatasetWithHallNumber(lattice, position, types, num_atom, hall_number, symprec, angle_tolerance);
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

            [rotations, translations, equivalent_atoms, num_operations] = Spglib.getSymmetryWithCollinearSpin(max_size, lattice, position, types, spins, num_atom, 1e-5);

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

            [rotations, translations, equivalent_atoms, num_operations] = Spglib.getSymmetryWithCollinearSpin(max_size, lattice, position, types, spins, num_atom, 1e-5, angle_tolerance);

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

            [rotations, translations, equivalent_atoms, num_operations] = Spglib.getSymmetryWithCollinearSpin(max_size, lattice, position, types, spins, num_atom, 1e-5, angle_tolerance, mag_symprec);

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

            [~, ~, ~, ~, ~, num_operations]  = Spglib.getSymmetryWithSiteTensors(max_size, lattice, position, types, tensors, 1, num_atom, true, true, symprec);
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

            dataset = Spglib.getDataset(lattice, position, types, num_atom, symprec);
            spacegroupType = Spglib.getSpacegroupTypeFromSymmetry(dataset.rotations, dataset.translations, dataset.n_operations, lattice, symprec);
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

            [rotations, translations, ~, ~, spin_flips, num_operations]  = Spglib.getSymmetryWithSiteTensors(max_size, lattice, position, types, tensors, 1, num_atom, true, true, symprec);
            testCase.assertEqual(num_operations, 8);

            time_reversals = idivide(1 - spin_flips, 2);
            spacegroupType = Spglib.getMagneticSpacegroupTypeFromSymmetry(rotations, translations, time_reversals, num_operations, lattice, symprec);
            testCase.assertEqual(spacegroupType.uni_number, 546);
        end
    end
end
