import numpy as np
import math

sqrd = lambda x : x*x
magnitude = lambda x : math.sqrt(sqrd(x[0]) + sqrd(x[1]) + sqrd(x[2]))

class Molecule:
    """
    A molecule is a group of atoms bonded together.

    Args:
        name (str): The name of the molecule.
        atoms_dict (dict): A dictionary of numpy arrays, each containing the index of an atom starting from 1 as keys, and its coordinates as values.
        origin_atom_index (int): Index of the atom that will be moved to the origin [0, 0, 0].
    """
    def __init__(self, name, atoms_dict, origin_atom_index):
        self.name = name
        self.number_atoms = len(atoms_dict)
        self.atoms = atoms_dict
        # Átomo que se llevará a la coordenada [0, 0, 0]
        self.origin_atom_index = origin_atom_index
        self.origin_atom = self.atoms[origin_atom_index]
        # Llevar el origin_atom a las coordenadas [0, 0, 0]
        for i in self.atoms:
            self.atoms[i] = self.atoms[i] - self.origin_atom

    def __sub__(self, other):
        """
        Subtracts the coordinates of two molecules.

        Args:
            other (Molecule): The molecule to be subtracted from self.

        Returns:
            numpy array: A numpy array containing the difference between the coordinates of the atoms of the two molecules.
        """
        if self.number_atoms != other.number_atoms:
            raise ValueError("Both molecules must have the same number of atoms to perform subtraction.")
        
        diff = {n: self.atoms[n] - other.atoms[n] for n in self.atoms}

        return diff


class TriCationicMolecule(Molecule):
    def __init__(
            self,
            name,
            atoms_dict,
            origin_atom_index,
            central_atoms_indexes,
            left_cation_indexes,
            right_cation_indexes
        ):
        """
        Initialize an aligned TriCationicMolecule in both states sp3 and sp2.
        
        Args:
            origin_atom_index (int): Index of the atom that will be moved to the origin [0, 0, 0].
            central_atoms_indexes (list): Contains the coordinates of top atom of the central cation and the right adjacent atom.
        """
        super().__init__(name, atoms_dict, origin_atom_index)
        # Inicializar atributos posibles (independientes)
        # Encontrar la matriz de alineamiento
        a = self.get_alignment_matrix(self.atoms[central_atoms_indexes[0]], self.atoms[central_atoms_indexes[1]])

        # Alinear coordenadas al plano xy
        for i in range(1, self.number_atoms+1):
            self.atoms[i] = np.dot(a, self.atoms[i])

        # Inicializar atributos que contendrán información calculada después de haber creado el objeto
        self.transitions = None # Se inicializa el atributo que contendrá las coordenadas de cada transición
        self.left_cation = Molecule(
            name = self.name+"_left_cation",
            atoms_dict = {i:self.atoms[i] for i in left_cation_indexes},
            origin_atom_index = left_cation_indexes[0]
        )
        self.right_cation = Molecule(
            name = self.name+"_right_cation",
            atoms_dict = {i:self.atoms[i] for i in right_cation_indexes},
            origin_atom_index = right_cation_indexes[0]
        )
        self.central_cation = Molecule(
            name = self.name+"_central_cation",
            atoms_dict = {i:self.atoms[i] for i in self.atoms if i not in left_cation_indexes and i not in right_cation_indexes},
            origin_atom_index = central_atoms_indexes[0]
        )

        # self.left_bridge = None
        # self.right_bridge = None

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

    def generate_gjf(self, n_transitions, template_file_path, output_directory):
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

            new_file_path = f"{output_directory}/{n+1}transition.gjf"

            with open(new_file_path, 'w') as file:
                file.write(template)
        
        print(f"File '{new_file_path}' created successfully.")