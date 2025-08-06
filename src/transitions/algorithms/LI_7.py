import math
import numpy as np

# Definir función de magnitud para tomar en cuenta el cambio de distancia entre átomos 
def magnitude(vector): 
    return math.sqrt(sum(pow(element, 2) for element in vector))

# -----------------------------------------------

def run_algorithm(n, sp3_atoms_list, sp2_atoms_list, difvecatomslist):
    """
    Runs the algorithm to calculate the coordinates of each transition. It takes into account the bond distances and angles that change during the transitions.
    
    Args:
        n (int): Number of transitions to calculate.
        sp3_atoms_list (list): List of numpy arrays containing atom coordinates of the sp3 molecule.
        sp2_atoms_list (list): List of numpy arrays containing atom coordinates of the sp2 molecule.
        difvecatomslist (list): List of numpy arrays containing the difference vectors between the coordinates of the sp3 and sp2 molecules divided by the number of transitions.
    Returns:
        numpy array: An array of arrays of atoms coordinates, each subarray represents a transition.
    """
    partialgeoms = [] # Lista para almacenar las geometrías parciales de cada transición

    # Pasar de sp2 a sp3 solo sumando difvecatomslist a las coordenadas de los átomos de la molécula sp3
    for i in range(n):
        partialgeom = []
        for j in range(len(sp3_atoms_list)):
            # Se suman las diferencias a las coordenadas de los átomos de la molécula sp3
            vec = (difvecatomslist[j]*(i+1)) + sp3_atoms_list[j]
            partialgeom.append(np.array(vec))
        # Se agrega la geometría parcial
        partialgeom = np.array(partialgeom)
        partialgeoms.append(partialgeom)

    print(np.array(partialgeoms))
    return np.array(partialgeoms)