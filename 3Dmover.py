#!/usr/bin/env python3

"""
3Dmover: A Tool for Manipulating 3D Protein Objects

3Dmover is a versatile tool designed to manipulate 3D protein objects in space. It provides functionality to
center the protein at the origin (0, 0, 0) or perform translations to move the protein in 3D space.

Usage: python 3Dmover.py <pdb_file> [-c] [-t X Y Z]

The tool reads a PDB file containing the atom coordinates of a protein object, applies the desired manipulation
(transformation), and saves the updated coordinates in a new PDB file.

Arguments:
    pdb_file         Path to the PDB file

Options:
    -c, --center     Center the protein at the origin (0, 0, 0)
    -t, --translate  Translate the protein by the specified amounts in X, Y, and Z directions

Example:
    python 3Dmover.py protein.pdb -c

    This will center the protein object in the specified PDB file at the origin (0, 0, 0) and save the updated
    coordinates in a new PDB file with the suffix '_centered' appended to the original filename.

    python 3Dmover.py protein.pdb -t 1.0 2.0 3.0

    This will translate the protein object in the specified PDB file by the amounts (1.0, 2.0, 3.0) in the X, Y,
    and Z directions, and save the updated coordinates in a new PDB file with the suffix '_translated' appended
    to the original filename.

Note: The PDB file must follow the standard format for atomic coordinate data.

For more information and options, refer to the documentation.
"""


import sys
import argparse


def parse_arguments():
    """
    Parse the command-line arguments for 3Dmover tool.

    Returns:
        argparse.Namespace: An object containing the parsed command-line arguments.

    Raises:
        argparse.ArgumentError: If the command-line arguments are invalid.
    """
    parser = argparse.ArgumentParser(description='3Dmover: A Tool for Manipulating 3D Protein Objects')

    parser.add_argument('pdb_file', type=str, help='Path to the PDB file')

    parser.add_argument('-c', '--center', action='store_true', help='Center the protein at the origin')
    parser.add_argument('-t', '--translate', nargs=3, type=float, metavar=('X', 'Y', 'Z'),
                        help='Translate the protein by the specified amounts in X, Y, and Z directions')

    return parser.parse_args()


def move_pdb(pdb_file, center=True, translation=None):
    """
    Read a PDB file and modify the atom coordinates by centering and/or translating them.

    Args:
        pdb_file (str): The path to the PDB file.
        center (bool, optional): Flag indicating whether to center the atom coordinates. Defaults to True.
        translation (tuple or None, optional): The translation vector as a tuple of three floats (x, y, z).
            If provided, the atom coordinates will be translated by adding the translation vector to each coordinate.
            Defaults to None.

    Returns:
        str: The path of the updated PDB file with modified coordinates.

    Raises:
        FileNotFoundError: If the specified PDB file does not exist.
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

    if center:
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

        centered_atoms = atoms

    if translation is not None:
        translation_vector = translation
        translated_atoms = []
        for atom in atoms:
            x = atom[0] + translation_vector[0]
            y = atom[1] + translation_vector[1]
            z = atom[2] + translation_vector[2]
            translated_atoms.append((x, y, z))

        centered_atoms = translated_atoms

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
    centered_pdb_file = pdb_file.replace('.pdb', '_modif.pdb')
    with open(centered_pdb_file, 'w') as file:
        file.writelines(centered_pdb_lines)

    return centered_pdb_file


if __name__ == '__main__':
    args = parse_arguments()
    move_pdb(args.pdb_file, center=args.center, translation=args.translate)
