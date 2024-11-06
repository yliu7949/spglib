classdef Spglib
    %SPGLIB An interface of spglib library.

    methods(Static)
        function version = getVersion(option)
            arguments
                option (1,1) string {mustBeMember(option, ["full", "major", "minor", "micro", "default"])} = "default"
            end

            switch option
                case "full"
                    version = symspg('spg_get_version_full');
                case "major"
                    version = symspg('spg_get_major_version');
                case "minor"
                    version = symspg('spg_get_minor_version');
                case "micro"
                    version = symspg('spg_get_micro_version');
                otherwise
                    version = symspg('spg_get_version');
            end
        end

        function commit = getCommit()
            commit = symspg('spg_get_commit');
        end

        function dataset = getDataset(lattice, position, types, num_atom, symprec, angle_tolerance)
            arguments
                lattice (3, 3) double
                position (:, 3) double
                types (:, 1) double
                num_atom (1, 1) uint16
                symprec (1, 1) double
                angle_tolerance double = []
            end

            if isempty(angle_tolerance)
                dataset = symspg('spg_get_dataset', lattice, position, types, num_atom, symprec);
            else
                dataset = symspg('spgat_get_dataset', lattice, position, types, num_atom, symprec, angle_tolerance);
            end
        end

        function dataset = getMagneticDataset(lattice, position, types, tensors, tensor_rank, num_atom, is_axial, symprec, angle_tolerance, mag_symprec)
            arguments
                lattice (3, 3) double
                position (:, 3) double
                types (:, 1) double
                tensors double
                tensor_rank (1, 1) uint16
                num_atom (1, 1) uint16
                is_axial (1, 1) logical
                symprec (1, 1) double
                angle_tolerance double = []
                mag_symprec double = []
            end

            tensors = reshape(tensors.', [], 1);
            if isempty(angle_tolerance)
                dataset = symspg('spg_get_magnetic_dataset', lattice, position, types, tensors, tensor_rank, num_atom, is_axial, symprec);
            else
                if isempty(mag_symprec)
                    error('Spglib:getMagneticDataset', "'mag_symprec' must be used together with 'angle_tolerance'");
                else
                    dataset = symspg('spgms_get_magnetic_dataset', lattice, position, types, tensors, tensor_rank, num_atom, is_axial, symprec, angle_tolerance, mag_symprec);
                end
            end
        end

        function dataset = getDatasetWithHallNumber(lattice, position, types, num_atom, hall_number, symprec, angle_tolerance)
            arguments
                lattice (3, 3) double
                position (:, 3) double
                types (:, 1) double
                num_atom (1, 1) uint16
                hall_number (1, 1) int16
                symprec (1, 1) double
                angle_tolerance double = []
            end

            if isempty(angle_tolerance)
                dataset = symspg('spg_get_dataset_with_hall_number', lattice, position, types, num_atom, hall_number, symprec);
            else
                dataset = symspg('spgat_get_dataset_with_hall_number', lattice, position, types, num_atom, hall_number, symprec, angle_tolerance);
            end
        end

        function [rotations, translations, equivalent_atoms, num_operations] = getSymmetryWithCollinearSpin(max_size, lattice, position, types, spins, num_atom, symprec, angle_tolerance, mag_symprec)
            arguments
                max_size (1, 1) uint16
                lattice (3, 3) double
                position (:, 3) double
                types (:, 1) double
                spins double
                num_atom (1, 1) uint16
                symprec (1, 1) double
                angle_tolerance double = []
                mag_symprec double = []
            end

            if isempty(angle_tolerance)
                [rotations, translations, equivalent_atoms, num_operations] = symspg('spg_get_symmetry_with_collinear_spin', max_size, lattice, position, types, spins, num_atom, symprec);
            else
                if isempty(mag_symprec)
                    [rotations, translations, equivalent_atoms, num_operations] = symspg('spgat_get_symmetry_with_collinear_spin', max_size, lattice, position, types, spins, num_atom, symprec, angle_tolerance);
                else
                    [rotations, translations, equivalent_atoms, num_operations] = symspg('spgms_get_symmetry_with_collinear_spin', max_size, lattice, position, types, spins, num_atom, symprec, angle_tolerance, mag_symprec);
                end
            end
        end

        function [rotations, translations, equivalent_atoms, primitive_lattice, spin_flips, num_operations] = getSymmetryWithSiteTensors(max_size, lattice, position, types, tensors, tensor_rank, num_atom, with_time_reversal, is_axial, symprec, angle_tolerance, mag_symprec)
            arguments
                max_size (1, 1) uint16
                lattice (3, 3) double
                position (:, 3) double
                types (:, 1) double
                tensors (:, :) double
                tensor_rank (1, 1) uint16
                num_atom (1, 1) uint16
                with_time_reversal (1, 1) logical
                is_axial (1, 1) logical
                symprec (1, 1) double
                angle_tolerance double = []
                mag_symprec double = []
            end

            tensors = reshape(tensors.', [], 1);
            if isempty(angle_tolerance)
                [rotations, translations, equivalent_atoms, primitive_lattice, spin_flips, num_operations] = symspg('spg_get_symmetry_with_site_tensors', max_size, lattice, position, types, tensors, tensor_rank, num_atom, with_time_reversal, is_axial, symprec);
            else
                if isempty(mag_symprec)
                    [rotations, translations, equivalent_atoms, primitive_lattice, spin_flips, num_operations] = symspg('spgat_get_symmetry_with_site_tensors', max_size, lattice, position, types, tensors, tensor_rank, num_atom, with_time_reversal, is_axial, symprec, angle_tolerance);
                else
                    [rotations, translations, equivalent_atoms, primitive_lattice, spin_flips, num_operations] = symspg('spgms_get_symmetry_with_site_tensors', max_size, lattice, position, types, tensors, tensor_rank, num_atom, with_time_reversal, is_axial, symprec, angle_tolerance, mag_symprec);
                end
            end
        end

        function spacegroup_type = getSpacegroupTypeFromSymmetry(rotation, translation, num_operations, lattice, symprec)
            arguments
                rotation (:, 3, 3) int32
                translation (:, 3) double
                num_operations (1, 1) int32
                lattice (3, 3) double
                symprec (1, 1) double
            end

            spacegroup_type = symspg('spg_get_spacegroup_type_from_symmetry', rotation, translation, num_operations, lattice, symprec);
        end

        function magnetic_spacegroup_type = getMagneticSpacegroupTypeFromSymmetry(rotation, translation, time_reversals, num_operations, lattice, symprec)
            arguments
                rotation (:, 3, 3) int32
                translation (:, 3) double
                time_reversals (:, 1) int32
                num_operations (1, 1) int32
                lattice (3, 3) double
                symprec (1, 1) double
            end

            magnetic_spacegroup_type = symspg('spg_get_magnetic_spacegroup_type_from_symmetry', rotation, translation, time_reversals, num_operations, lattice, symprec);
        end
    end
end
