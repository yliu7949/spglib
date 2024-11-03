classdef SpglibTest < matlab.unittest.TestCase

    methods(TestClassSetup)
        function setupPath(~)
            spglibFolder =  fullfile(fileparts(mfilename("fullpath")), "spglib");
            addpath(spglibFolder);
        end
    end

    methods (Test)
        function getMajorVersionTest(testCase)
            version = Spglib.getMajorVersion();
            testCase.assertEqual(version, 2);
        end

        function getMinorVersionTest(testCase)
            version = Spglib.getMinorVersion();
            testCase.assertEqual(version, 5);
        end

        function getMicroVersionTest(testCase)
            version = Spglib.getMicroVersion();
            testCase.assertEqual(version, 1);
        end

        function getDatasetTest(testCase)
            clc;
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
            dataset = Spglib.getDataset(lattice, position, types, num_atom, symprec);
            testCase.assertEqual(dataset.spacegroup_number, 186);
        end

        function getDatasetWithHallNumberTest(testCase)
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
    end
end
