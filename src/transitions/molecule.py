import numpy as np
import math

sqrd = lambda x : x*x
magnitude = lambda x : math.sqrt(sqrd(x[0]) + sqrd(x[1]) + sqrd(x[2]))

class Molecule:
    """
    A molecule is a group of atoms bonded together.

    Args:
        name (str): The name of the molecule.
        n_atoms (int): The number of atoms in the molecule.
    """
    def __init__(self, name, atoms_list):
        self.name = name
        self.number_atoms = len(atoms_list)
        self.atoms = {i:atom for i, atom in zip(range(1, self.number_atoms+1), atoms_list)}

class TriCationicMolecule:
    def __init__(self, origin_atom_index,central_atoms_indexes, c_cation_indexes, l_cation_indexes, r_cation_indexes, l_bridge_indexes, r_bridge_indexes, sp3_coords, sp2_coords):
        """
        Initialize an aligned TriCationicMolecule in both states sp3 and sp2.
        
        Args:
            origin_atom_index (int): Index of the atom that will be moved to the origin [0, 0, 0].

            central_atoms_indexes (list): Contains the coordinates of top atom of the central cation and the right adjacent atom.

            c_cation_indexes (list): List of indexes of the central cation atoms.

            l_cation_indexes (list): List of indexes of the left cation atoms.

            r_cation_indexes (list): List of indexes of the right cation atoms.
            
            l_bridge_indexes (list): List of indexes of the left bridge atoms.

            r_bridge_indexes (list): List of indexes of the right bridge atoms.

            sp3_coords (list): List of arrays containing the coordinates of the atoms in the sp3 state.

            sp2_coords (list): List of arrays containing the coordinates of the atoms in the sp2 state.
        """
        # Inicializar atributos posibles (independientes)
        self.number_atoms = len(sp3_coords)
        # Átomo que se llevará a la coordenada [0, 0, 0]
        origin_atom = sp3_coords[origin_atom_index]

        # Llevar el origin_atom a las coordenadas [0, 0, 0]
        for i in range(self.number_atoms):
            sp3_coords[i] = sp3_coords[i]-origin_atom
            sp2_coords[i] = sp2_coords[i]-origin_atom

        # Encontrar las matrices de alineaminto del estado sp3 y sp2
        a_sp3 = self.get_alignment_matrix(sp3_coords[central_atoms_indexes[0]], sp3_coords[central_atoms_indexes[1]])
        a_sp2 = self.get_alignment_matrix(sp2_coords[central_atoms_indexes[0]], sp2_coords[central_atoms_indexes[1]])

        # Alinear coordenadas de los estados sp3 y sp2 a plano xy
        for i in range(self.number_atoms):
            sp3_coords[i] = np.dot(a_sp3, sp3_coords[i])
            sp2_coords[i] = np.dot(a_sp2, sp2_coords[i])

        self.sp3 = sp3_coords

        # Inicializar los atributos que contienen las coordenadas alineadas       
        self.central_cation = {i: sp3_coords[i] for i in central_atoms_indexes}

        self.left_cation = {i: sp3_coords[i] for i in l_cation_indexes}

        self.right_cation = {i: sp3_coords[i] for i in r_cation_indexes}

        self.left_bridge = {i: sp3_coords[i] for i in l_bridge_indexes}

        self.right_bridge = {i: sp3_coords[i] for i in r_bridge_indexes}

        # Inicializar los atributos que contienen las diferencias entre los vectores de sp3 y sp2
        self.l_cation_differences = {i: (sp2_coords[i] - sp3_coords[i]) for i in central_atoms_indexes}

        self.r_cation_differences = {i: (sp2_coords[i] - sp3_coords[i]) for i in r_cation_indexes}

        self.l_bridge_differences = {i: (sp2_coords[i] - sp3_coords[i]) for i in l_bridge_indexes}

        self.r_bridge_differences = {i: (sp2_coords[i] - sp3_coords[i]) for i in r_bridge_indexes}

    def get_alignment_matrix(self, w_1, w_2):
        """
        Finds the matrix A that aligns the central vectors of the molecule to the xy-plane.

        Args:
            w_1 (numpy array): Contains the coordinates of top atom of the central cation.
            w_2 (numpy array): Contains the coordinates of the right adjacent atom.
        """
        # Vectores unitarios desalineados|
        j = w_1/magnitude(w_1)
        i = w_2 + np.dot(j, magnitude(w_2 * np.sin(np.deg2rad(35.5))))
        i = i/magnitude(i)
        k = np.cross(i, j)
        a = np.array([
            [i[0], j[0], k[0]],
            [i[1], j[1], k[1]],
            [i[2], j[2], k[2]]
        ])

        return np.linalg.inv(np.array(a))

    def generate_transitions(self, n_transitions):
        """
        Calculate the coordinates of each transition.

        Args:
            n_transitions (int): Number of transitions that will be calculated.
        
        Returns:
            numpy array: A n_transitions x self.number_atoms x 3 array, each self.number_atoms x 3 array represents a transition.
        """
        self.transitions = np.zeros((n_transitions, self.number_atoms, 3))

        for n in range(n_transitions):
            self.transitions[n] = np.array(self.sp3) # Iniciar la transición con las coordenadas del estado sp3

    def generate_gjf(self, n_transitions, template_file_path):
        """
        Generates a .gjf document for each list of atoms given.

        Args:
            n_transitions (int): Number of transitions.
            template_file_path (str): Path of the .gjf template used to create the new files.
        """
        for n in range(n_transitions):
            # Insertar las coordenadas en la plantilla template.txt y crear los archivos .gjf
            with open(template_file_path, 'r') as file:
                template = file.read()
                # Número de átomos
                for x in range(self.number_atoms):
                    # Número de dimensiones para el vector del átomo
                    for y in range(3):
                        template = template.replace(f" {x+1}[n][{y}]", f" {self.transitions[n,x,y]}")

            new_file_path = f"molecules\LI-7\gjf_files\{n+1}transition.gjf"

            with open(new_file_path, 'w') as file:
                file.write(template)
        
        print(f"File '{new_file_path}' created successfully.")