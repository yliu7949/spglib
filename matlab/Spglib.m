classdef Spglib
    %SPGLIB An interface of spglib library.

    methods(Static)
        function version = getMajorVersion()
            version = symspg('spg_get_major_version');
        end

        function version = getMinorVersion()
            version = symspg('spg_get_minor_version');
        end

        function version = getMicroVersion()
            version = symspg('spg_get_micro_version');
        end

        function dataset = getDataset(lattice, position, types, num_atom, symprec)
            arguments
                lattice (3, 3) double
                position (:, 3) double
                types (:, 1) double
                num_atom (1, 1) uint16
                symprec (1, 1) double
            end

            dataset = symspg('spg_get_dataset', lattice, position, types, num_atom, symprec);
        end

        function dataset = getDatasetWithHallNumber(lattice, position, types, num_atom, hall_number, symprec)
            arguments
                lattice (3, 3) double
                position (:, 3) double
                types (:, 1) double
                num_atom (1, 1) uint16
                hall_number (1, 1) int16
                symprec (1, 1) double
            end

            dataset = symspg('spg_get_dataset_with_hall_number', lattice, position, types, num_atom, hall_number, symprec);
        end
    end
end
