# 3Dmover

3Dmover is a versatile tool designed to manipulate 3D protein objects in space. It provides functionality to center the protein at the origin (0, 0, 0) or perform translations to move the protein in 3D space. The tool reads a PDB file containing the atom coordinates of a protein object, applies the desired manipulation (transformation), and saves the updated coordinates in a new PDB file.
Features

- Center the protein at the origin (0, 0, 0).
- Perform translations to move the protein in 3D space.
- Update the PDB file with the manipulated coordinates.


## Installation

1. Clone the repository:

`git clone https://github.com/your_username/3Dmover.git`

2. Change into the project directory:

`cd 3Dmover`

3. Install the required dependencies:

`pip install -r requirements.txt`

## Usage

`python 3Dmover.py <pdb_file> [-c] [-t X Y Z]`

- <pdb_file>: Path to the PDB file containing the protein object coordinates.
- -c, --center: Center the protein at the origin (0, 0, 0).
- -t X Y Z, --translate X Y Z: Translate the protein by the specified amounts in the X, Y, and Z directions.