import numpy as np
import math as math

def parse_atom_coordinates(file_path):
    """
    Reads a text file containing atom coordinates and returns them as a list of numpy arrays.
    
    Args:
        file_path (str): Path to the text file containing atom coordinates
    
    Returns:
        list: List of numpy arrays, each containing [x, y, z] coordinates for an atom
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

def generate_gjf(atoms_array, n_transitions, num_atoms, template_file_path):
    """
    Generates a .gjf document for each list of atoms given.

    Args: 
        final_atom_list (numpy array): An array of arrays of atoms coordinates, each subarray represents a transition.
        n_transitions (int): Number of transitions.
        num_atoms (int): Number of atoms.
        template_file_path (str): Path of the .gjf template used to create the new files.
    """
    for n in range(n_transitions):
        # Insertar las coordenadas en la plantilla .gjf
        with open(template_file_path, 'r') as file:
            template = file.read()
            # Número de átomos
            for x in range(num_atoms):
                # Número de dimensiones para el vector del átomo
                for y in range(3):
                    template = template.replace(f" {x+1}[n][{y}]", f" {atoms_array[n,x,y]}")

        new_file_path = f"molecules\LI-7\gjf_files\{n+1}transition.gjf"

        with open(new_file_path, 'w') as file:
            file.write(template)
    
    print(f"File '{new_file_path}' created successfully.")
