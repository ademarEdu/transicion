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

class TriCationicMolecule(Molecule):
    def __init__(
            self,
            name,
            atoms_list,
            origin_atom_index,
            central_atoms_indexes
        ):
        """
        Initialize an aligned TriCationicMolecule in both states sp3 and sp2.
        
        Args:
            origin_atom_index (int): Index of the atom that will be moved to the origin [0, 0, 0].
            central_atoms_indexes (list): Contains the coordinates of top atom of the central cation and the right adjacent atom.
        """
        super().__init__(name, atoms_list)
        # Inicializar atributos posibles (independientes)
        # Átomo que se llevará a la coordenada [0, 0, 0]
        origin_atom = self.atoms[origin_atom_index]

        # Llevar el origin_atom a las coordenadas [0, 0, 0]
        for i in range(1, self.number_atoms+1):
            self.atoms[i] = self.atoms[i] - origin_atom

        # Encontrar la matriz de alineamiento
        a = self.get_alignment_matrix(self.atoms[central_atoms_indexes[0]], self.atoms[central_atoms_indexes[1]])

        # Alinear coordenadas al plano xy
        for i in range(1, self.number_atoms+1):
            self.atoms[i] = np.dot(a, self.atoms[i])

        # Inicializar atributos que contendrán información calculada después de haber creado el objeto
        self.transitions = None # Se inicializa el atributo que contendrá las coordenadas de cada transición
        # self.left_cation = None
        # self.right_cation = None
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