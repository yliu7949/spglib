classdef SpglibTest < matlab.unittest.TestCase

    properties

    end

    methods(TestClassSetup)
        function setupPath(~)
            addpath("cmake-build-debug");
        end
    end

    methods(TestMethodSetup)
        
    end

    methods (Test)
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
    end
end
