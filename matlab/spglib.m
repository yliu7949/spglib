classdef Spglib 
    %SPGLIB An interface of spglib library.

    methods(Static)
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
    end
end
