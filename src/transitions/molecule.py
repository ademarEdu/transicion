import numpy as np
import math
import tempfile
import subprocess
import os

sqrd = lambda x : x*x
magnitude = lambda x : math.sqrt(sqrd(x[0]) + sqrd(x[1]) + sqrd(x[2]))

class Atom:
    """
    An object that represents an atom.
    
    Args:
        name (str): The name of the atom.
        coordinates (numpy array): A numpy array containing the x, y, z coordinates of the atom.
    """

    def __init__(self, name, coordinates):
        self.name = name
        self.coordinates = coordinates

    def __sub__(self, other):
        return Atom(self.name, self.coordinates - other.coordinates)

    def __add__(self, other):
        return Atom(self.name, self.coordinates + other.coordinates)
    
    def __mul__(self, scalar):
        return Atom(self.name, self.coordinates * scalar)

    def __truediv__(self, scalar):
        return Atom(self.name, self.coordinates / scalar)

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
        for i in atoms_dict:
            self.atoms[i] -= self.origin_atom

    def __copy__(self):
        return Molecule(
            name = self.name+"_copy",
            atoms_dict = {i:self.atoms[i] for i in self.atoms},
            origin_atom_index = self.origin_atom_index
        )

    def __sub__(self, other):
        """
        Subtracts the coordinates of two molecules.

        Args:
            other (Molecule): The molecule to be subtracted from self.
        """
        if self.number_atoms != other.number_atoms:
            raise ValueError("Both molecules must have the same number of atoms to perform subtraction.")
        diff = self.__copy__()
        for i in self.atoms:
            diff.atoms[i] -= other.atoms[i]
        return diff
    
    def __add__(self, other):
        """
        Adds the coordinates of two molecules.

        Args:
            other (Molecule): The molecule to be added to self.
        """
        if self.number_atoms != other.number_atoms:
            raise ValueError("Both molecules must have the same number of atoms to perform addition.")
        add = self.__copy__()
        for i in add.atoms:
            add.atoms[i] += other.atoms[i]
        return add

    def __mul__(self, scalar):
        """
        Multiplies the coordinates of the molecule by a scalar.

        Args:
            scalar (float): The scalar to multiply the coordinates by.
        """
        mult = self.__copy__()
        for i in self.atoms:
            mult.atoms[i] *= scalar
        return mult

    def __truediv__(self, scalar):
        """
        Divides the coordinates of the molecule by a scalar.

        Args:
            scalar (float): The scalar to divide the coordinates by.
        """
        div = self.__copy__()
        for i in self.atoms:
            div.atoms[i] /= scalar
        return div

    def gjf_file_content(self, header_file):
        """
        Generates a .gjf document for the molecule.

        Args:
            output_directory (str): Path of the directory where the .gjf file will be created.
        
        Returns:
            list: A list of strings, each string is a line in the .gjf file. The 8th line contains the atoms of the molecule separated by a newline(\n).
        """
        # Obtener los datos de la cabecera
        with open(header_file, 'r') as template_file:
            lines = template_file.readlines()
        # Generar la líneas que contienen los átomos
        atoms_str = ""
        for i in self.atoms:
            atom = self.atoms[i]
            atoms_str += f"{atom.name}    {atom.coordinates[0]:.6f}    {atom.coordinates[1]:.6f}    {atom.coordinates[2]:.6f}\n"
        lines[7] = atoms_str
        
        return lines

    def show_gjf(self, header_file):
        """
        Shows a temporary .gjf file of the molecule with GaussView.

        Args:
            header_file (str): Path of the header file to be used in the .gjf file. Tipically "molecules/header.txt".
        """
        gaussview_path = "C:/Users/dell/Desktop/GaussView_6.0.16+Win64/GaussView_6.0.16_win64/gview.exe"
        # Generar el archivo temporal .gjf
        with tempfile.NamedTemporaryFile(mode="w+", suffix=".gjf", delete=False) as temp_file:
            temp_file.writelines(self.gjf_file_content(header_file))
            temp_path = temp_file.name  # Path to the temp file
        # Abrir el archivo temporal con GaussView
        subprocess.run([gaussview_path, temp_path])
        # Eliminar el archivo temporal
        os.remove(temp_path)

    def save_gjf(self, header_file, new_file_path):
        """
        Saves a .gjf file of the molecule in the specified directory.
        """
        # Generar el archivo .gjf
        with open(f"{new_file_path}", "w") as gjf_file:
            gjf_file.writelines(self.gjf_file_content(header_file))

class Cation(Molecule):
    def __init__(self, name, atoms_dict, origin_atom_index):
        """
        Initialize a Cation molecule that belongs to a TriCationicMolecule.
        
        Args:
            name (str): The name of the cation.
            atoms_dict (dict): A dictionary of numpy arrays, each containing the index of an atom starting from 1 as keys, and its coordinates as values.
            origin_atom_index (int): Index of the atom that will be moved to the origin [0, 0, 0].
        """
        super().__init__(name, atoms_dict, origin_atom_index)

    def move_to_origin(self):
        """
        Moves self.origin_atom of the cation to the origin [0, 0, 0].
        """
        for i in self.atoms:
            self.atoms[i] -= self.origin_atom

    def back_to_original(self):
        """
        Moves the cation back to its original position in the TriCationicMolecule.

        Args:
            origin_coordinates (numpy array): The coordinates of the atom that was moved to the origin.
        """
        for i in self.atoms:
            self.atoms[i] += self.origin_atom

class TriCationicMolecule(Molecule):
    def __init__(
            self,
            name,
            atoms_dict,
            origin_atom_index,
            central_atoms_indexes,
            left_cation_indexes,
            right_cation_indexes,
            left_bridge_indexes,
            right_bridge_indexes
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
        a = self.get_alignment_matrix(self.atoms[central_atoms_indexes[0]].coordinates, self.atoms[central_atoms_indexes[1]].coordinates)

        # Alinear coordenadas al plano xy
        for i in range(1, self.number_atoms+1):
            self.atoms[i].coordinates = np.dot(a, self.atoms[i].coordinates)

        # Inicializar atributos que contendrán información calculada después de haber creado el objeto
        self.left_cation = Cation(
            name = self.name+"_left_cation",
            atoms_dict = {i:self.atoms[i] for i in left_cation_indexes},
            origin_atom_index = left_cation_indexes[0]
        )
        self.right_cation = Cation(
            name = self.name+"_right_cation",
            atoms_dict = {i:self.atoms[i] for i in right_cation_indexes},
            origin_atom_index = right_cation_indexes[0]
        )
        self.left_bridge = Cation(
            name = self.name+"_left_bridge",
            atoms_dict = {i:self.atoms[i] for i in left_bridge_indexes},
            origin_atom_index = left_bridge_indexes[0]
        )
        self.right_bridge = Cation(
            name = self.name+"_right_bridge",
            atoms_dict = {i:self.atoms[i] for i in right_bridge_indexes},
            origin_atom_index = right_bridge_indexes[0]
        )
        self.central_cation = Cation(
            name = self.name+"_central_cation",
            atoms_dict = {i:self.atoms[i] for i in self.atoms if i not in left_cation_indexes and i not in right_cation_indexes and i not in left_bridge_indexes and i not in right_bridge_indexes},
            origin_atom_index = central_atoms_indexes[0]
        )

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
    
    def build_molecule(self, left_cation, right_cation, central_cation, left_bridge, right_bridge):
        """
        Returns a TriCationicMolecule containing the specified parts.

        Args:
            left_cation (bool): Indicates if the left cation of the molecule is included.
            right_cation (bool): Indicates if the right cation of the molecule is included.
            central_cation (bool): Indicates if the central cation of the molecule is included.
            left_bridge (bool): Indicates if the left bridge of the molecule is included.
            right_bridge (bool): Indicates if the right bridge of the molecule is included.
        """
        # # Mover los cationes a su posición original
        # left_cation.back_to_original()
        # right_cation.back_to_original()
        # central_cation.back_to_original()
        # left_bridge.back_to_original()
        # right_bridge.back_to_original()

        # # Unir los átomos de las partes en un solo diccionario
        # self.atoms = {**left_cation.atoms, **right_cation.atoms, **central_cation.atoms, **left_bridge.atoms, **right_bridge.atoms}