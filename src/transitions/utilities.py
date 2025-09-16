import numpy as np
import math as math

def parse_atom_coordinates(file_path):
    """
    Reads a text file containing atom coordinates and returns them as a list of numpy arrays.
    
    Args:
        file_path (str): Path to the text file containing atom coordinates.
    
    Returns:
        list: list of numpy arrays. n rows x 3 columns, each row containing the [x, y, z] coordinates of an atom.
    """
    atoms = []
    with open(file_path, 'r') as file:
        for line in file:
            # Skip empty lines
            if not line.strip():
                continue
            # Split the line into components
            parts = line.split()
            # The first part is the atom type, the rest are coordinates
            if len(parts) >= 4:  # Ensure we have atom type + x,y,z coordinates
                try:
                    # Convert string coordinates to floats
                    coords = np.array([float(parts[1]), float(parts[2]), float(parts[3])])
                    atoms.append(coords)
                except (ValueError, IndexError):
                    # Skip lines that don't have valid coordinates
                    continue
    return atoms

def align_molecule_to_z_axis(atoms_list, origin_atom):
    """
    Aligns a molecule's atoms to the z-axis by rotating around the x-axis.
    
    Args:
        atoms_list (list): List of numpy arrays containing atom coordinates.
        origin_atom (numpy array): The coordinates of the atom that will be moved to the origin [0, 0, 0].
    
    Returns:
        list: List of atom coordinates regarding the aligned molecule.
    """    
    # Se lleva el origin_atom a las coordenadas [0, 0, 0]
    for i in range(len(atoms_list)):
        atoms_list[i] = atoms_list[i]-origin_atom

    atoms_list = np.array(atoms_list)

    # Calcular el ángulo de rotación necesario con respecto al eje x
    tetax = -(math.atan(atoms_list[0, 2]/atoms_list[0, 1])) # atan retorna valores de -PI/2 a PI/2
    # Calcular la matriz de rotación respecto al eje x
    matriz_rot_x = np.array([[1, 0, 0],
                            [0, math.cos(tetax), -math.sin(tetax)],
                            [0, math.sin(tetax), math.cos(tetax)]])

    # Aplicar la matriz de rotación a cada átomo de la molécula sp3
    for i in range(len(atoms_list)):
        atoms_list[i] = np.dot(matriz_rot_x, atoms_list[i])
    
    return atoms_list
