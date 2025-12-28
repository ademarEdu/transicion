import numpy as np
import math as math
from transitions.molecule import Atom

def get_atoms_dict(file_path):
    """
    Reads a text file containing atom coordinates and returns them as a dictionary of Atom objects.
    
    Args:
        file_path (str): Path to the text file containing atom coordinates.

    Returns:
        dict: A dictionary of Atom objects, keyed by atom index(starts from 1).
    """
    atoms = {}
    with open(file_path, 'r') as file:
        i = 1
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
                    atoms[i] = Atom(parts[0], coords)
                except (ValueError, IndexError):
                    # Skip lines that don't have valid coordinates
                    continue
            i += 1
            
    return atoms
