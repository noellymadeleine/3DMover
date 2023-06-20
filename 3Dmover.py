#!/usr/bin/env python3

"""
3Dmover: A Tool for Manipulating 3D Protein Objects

3Dmover is a versatile tool designed to manipulate 3D protein objects in space. It provides functionality to
center the protein at the origin (0, 0, 0) or perform translations to move the protein in 3D space.

Usage: python 3Dmover.py <pdb_file>

The tool reads a PDB file containing the atom coordinates of a protein object, applies the desired manipulation
(transformation), and saves the updated coordinates in a new PDB file.

Example:
    python 3Dmover.py protein.pdb

    This will center the protein object in the specified PDB file at the origin (0, 0, 0) and save the updated
    coordinates in a new PDB file with the suffix '_centered' appended to the original filename.

Note: The PDB file must follow the standard format for atomic coordinate data.

For more information and options, refer to the documentation.
"""


import sys


def center_pdb(pdb_file):
    """
    Centers a protein object represented in a PDB file.

    Given a PDB file containing atom coordinates of a protein object, this function calculates the center of mass
    and shifts the coordinates to center the protein at (0, 0, 0) in 3D space. It updates the PDB file with the
    centered coordinates and returns the path of the updated PDB file.

    Args:
        pdb_file (str): The path to the PDB file.

    Returns:
        str: The path of the updated PDB file with centered coordinates.

    Raises:
        FileNotFoundError: If the provided PDB file does not exist.

    """
    # Read the PDB file and extract the atom coordinates
    with open(pdb_file, 'r') as file:
        pdb_lines = file.readlines()

    atoms = []
    for line in pdb_lines:
        if line.startswith('ATOM'):
            x = float(line[30:38])
            y = float(line[38:46])
            z = float(line[46:54])
            atoms.append((x, y, z))

    # Calculate the center of mass
    total_atoms = len(atoms)
    sum_x = sum(atom[0] for atom in atoms)
    sum_y = sum(atom[1] for atom in atoms)
    sum_z = sum(atom[2] for atom in atoms)

    center_x = sum_x / total_atoms
    center_y = sum_y / total_atoms
    center_z = sum_z / total_atoms

    # Shift the coordinates to center the protein
    centered_atoms = []
    for atom in atoms:
        x = atom[0] - center_x
        y = atom[1] - center_y
        z = atom[2] - center_z
        centered_atoms.append((x, y, z))

    # Update the PDB file with centered coordinates
    centered_pdb_lines = []
    atom_counter = 0  # Counter for lines starting with 'ATOM'
    for line in pdb_lines:
        if line.startswith('ATOM'):
            line = line[:30] + f"{centered_atoms[atom_counter][0]:8.3f}" + \
                   f"{centered_atoms[atom_counter][1]:8.3f}" + f"{centered_atoms[atom_counter][2]:8.3f}" + line[54:]
            atom_counter += 1
        centered_pdb_lines.append(line)

    # Write the updated PDB file
    centered_pdb_file = pdb_file.replace('.pdb', '_centered.pdb')
    with open(centered_pdb_file, 'w') as file:
        file.writelines(centered_pdb_lines)

    return centered_pdb_file


if __name__ == '__main__':
    center_pdb(sys.argv[1])


