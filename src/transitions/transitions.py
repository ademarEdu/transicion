import numpy as np
import transitions.utilities as utils
from transitions.molecule import Molecule, TriCationicMolecule
import sqlite3
import json

def get_data(molecule_name):
    """
    Get the necessary data from the database to calculate the transitions of a specific molecule.
    
    Args:
        molecule_name (str): Name of the molecule to be studied. It has to be the same as the folder name in the molecules directory.

    Returns:
        dict: A dictionary containing the necessary data to calculate the transitions of the molecule.
    """
    conn = sqlite3.connect("database/molecules.db")
    c = conn.cursor()
    c.execute(f"SELECT * FROM molecules WHERE name='{molecule_name}'")
    data = c.fetchone()
    data = {
        "name": data[0],
        "origin_atom_index_sp3": data[1],
        "central_atoms_indexes_sp3": json.loads(data[2]),
        "origin_atom_index_sp2": data[3],
        "central_atoms_indexes_sp2": json.loads(data[4]),
        "left_cation_indexes": json.loads(data[5]),
        "right_cation_indexes": json.loads(data[6]),
        "left_bridge_indexes": json.loads(data[7]),
        "right_bridge_indexes": json.loads(data[8]),
        "left_bridge_correspondance": json.loads(data[9]),
        "right_bridge_correspondance": json.loads(data[10])
        }

    conn.close()
    return data

class Transitions:
    def __init__(self, molecule_name):
        """
        This class represents the transitions generator. It contains the necessary methods to generate the transitions of a molecule from sp3 to sp2 hybridization.

        Args:
            molecule_name (str): Name of the molecule to be studied. It has to be the same as the folder name in the molecules directory.
        """
        # Rutas de los archivos que contienen las coordenadas de los átomos de las moléculas
        sp3_path = f"molecules/{molecule_name}/sp3-coordinates.txt"
        sp2_path = f"molecules/{molecule_name}/sp2-coordinates.txt"
        self.header_file = f"molecules/header.txt"

        # Obtener los datos necesarios para calcular las transiciones
        self.mol_data = get_data(molecule_name)

        # Inicializar las moléculas
        self.sp3_molecule = TriCationicMolecule(
            name = molecule_name+"_sp3",
            atoms_dict = utils.get_atoms_dict(sp3_path),
            origin_atom_index = self.mol_data["origin_atom_index_sp3"],
            central_atoms_indexes= self.mol_data["central_atoms_indexes_sp3"],
            left_cation_indexes= self.mol_data["left_cation_indexes"],
            right_cation_indexes= self.mol_data["right_cation_indexes"],
            left_bridge_indexes= self.mol_data["left_bridge_indexes"],
            right_bridge_indexes= self.mol_data["right_bridge_indexes"]
        )

        self.sp2_molecule = TriCationicMolecule(
            name = molecule_name+"_sp2",
            atoms_dict = utils.get_atoms_dict(sp2_path),
            origin_atom_index = self.mol_data["origin_atom_index_sp2"],
            central_atoms_indexes= self.mol_data["central_atoms_indexes_sp2"],
            left_cation_indexes= self.mol_data["left_cation_indexes"],
            right_cation_indexes= self.mol_data["right_cation_indexes"],
            left_bridge_indexes= self.mol_data["left_bridge_indexes"],
            right_bridge_indexes= self.mol_data["right_bridge_indexes"]
        )

        # Ejemplos del uso de algunos métodos de la clase Molecule y TriCationicMolecule

        # Ejemplo de cómo mostrar o guardar un archivo .gjf
        # sp2_molecule.show_gjf(header_file)
        # sp2_molecule.save_gjf(header_file, f"molecules/{molecule_name}/gjf_files/trial.gjf")

        # # Ejemplo de cómo ver los componentes de una molécula tricatiónica
        # sp3_molecule.left_cation.back_to_original()
        # sp3_molecule.left_cation.show_gjf(header_file)
        # sp3_molecule.right_cation.show_gjf(header_file)
        # sp3_molecule.central_cation.show_gjf(header_file)
        # sp3_molecule.right_bridge.show_gjf(header_file)
        # sp3_molecule.left_bridge.show_gjf(header_file)

        # # Ejemplo de cómo generar una nueva molécula a partir de las partes de una molécula tricatiónica
        # new_sp3_molecule = sp3_molecule.build_molecule(
        #     left_cation=False,
        #     right_cation=False,
        #     left_bridge=True,
        #     right_bridge=True
        # )

        # --------------------------------------

    def gen_bridges_data(self, indexes_dictionaries):
        """
        Generate the values of r, r_prime and delta_r to adjust the bond legths of the bridges during the transitions. You can see more information about what these variables represent and how they are calculated in the notebooks/AdjustingBondLengths.ipynb notebook.

        Args:
            indexes_dictionaries (list): A list of two dictionaries, each one contains the indexes of the atoms involved in the bridges of the molecule. indexes_dictionaries[0] corresponds to the left bridge and indexes_dictionaries[1] corresponds to the right bridge. Each dictionary has the following structure:

            {
            'c': (index_c_0, index_c_1),
            'a': (index_a_0, index_a_1),
            'b': (index_b_0, index_b_1)
            }
        
        Returns:
            dict: A dictionary with the keys 'left' and 'right'. Each key points to another dictionary with the keys 'c_1', 'a_0' and 'b_0' that are related to the respective r, r_prime and delta_r values of each atom.
        """
        # Calcular datos necesarios para ajustar longitudes de enlace
        bridges_data = {'left': {}, 'right': {}}

        # c[0]: el átomo de carbono estacionario
        # c[1]: el átomo de carbono que se moverá durante la Transición
        # a: átomos de hidrógeno que se quedarán en la molécula durante la Transición
        # b: átomos de hidrógeno que saldrán de la molécula para dejar espacio al doble enlace

        # Recorreremos todas las llaves de indexes dictionaries (c, a, b)
        for i in ['c', 'a', 'b']:
            for indexes_dict, side in zip(indexes_dictionaries, ['left', 'right']):
                # Si i es 'c', el átomo con el que tratamos es c[1]; de lo contrario, el átomo con el que tratamos es a[0] o b[0]
                x = indexes_dict[i][1] if i == 'c' else indexes_dict[i][0]
                c_0 = indexes_dict['c'][0] 

                # v_1x = x - c_0
                # Donde x es la posición del átomo en el estado sp3 
                v1 = self.sp3_molecule.atoms[x] - self.sp3_molecule.atoms[c_0]

                # v'_x = x - c_0
                # Donde x es la posición del átomo en el estado sp2
                v_prime = self.sp2_molecule.atoms[x] - self.sp2_molecule.atoms[c_0]

                bridges_data[side][i+('_1' if i == 'c' else '_0')] = {
                    'r': v1.magnitude(),
                    'r_prime': v_prime.magnitude(),
                    'delta_r': (v_prime.magnitude() - v1.magnitude())/self.n
                }

        return bridges_data

    def adjust_bridges(self, n, indexes_dictionaries, bridges_data):
        """
        Adjust the bond lengths of the bridges at the nth transition.

        Args:
            n (int): The number of the current transition.
            indexes_dictionaries (list): A list of dictionaries containing the indexes of the atoms involved in the bridges of the molecule.
            bridges_data (dict): A dictionary containing the necessary data to adjust the bond lengths of the bridges during the transitions.

        Returns:
            Molecule: A new molecule with the bond lengths of the bridges adjusted at the nth transition.
        """
        new_sp3_molecule = self.sp3_molecule.__copy__()

        for i in ['c', 'a', 'b']:
            for side, indexes_dict in zip(['left', 'right'], indexes_dictionaries):
                # Si i es 'c', el átomo con el que tratamos es c[1]; de lo contrario, el átomo con el que tratamos es a[0] o b[0]
                x = indexes_dict[i][1] if i == 'c' else indexes_dict[i][0]
                c_0 = indexes_dict['c'][0]

                # v_2x = x - c_0
                # Donde x es la posición del átomo en el iésimo paso 
                v2 = new_sp3_molecule.atoms[x] - new_sp3_molecule.atoms[c_0]
                
                # nombre del átomo con el que tratamos (i.e. 'c_1', 'a_0' o 'b_0')
                name_atom = i+('_1' if i == 'c' else '_0')
                # ||v_1x||
                r = bridges_data[side][name_atom]['r']
                # delta_r
                delta_r = bridges_data[side][name_atom]['delta_r']
                
                # v_3x = v_2x/ ||v_2x|| (||v_1x|| - ||v_2x|| + delta_r*i)
                # Donde i es el número de paso actual
                v3 = (v2/v2.magnitude()) * (r - v2.magnitude() + delta_r * n )
                
                # Si i es 'c', el átomo con el que tratamos es c[1] y solo se ajusta la posición de ese átomo
                if i == 'c':
                    new_sp3_molecule.atoms[x] += v3
                # De lo contrario, el átomo con el que tratamos es a[0] o b[0] y se ajusta la posición de ese átomo y del átomo simétrico con respecto al átomo c[0], es decir, el átomo a[1] o b[1]
                else:
                    # Ajustar la posición del átomo con el que tratamos
                    new_sp3_molecule.atoms[x] += v3

                    # Ajustar la posición del átomo simétrico con respecto a c
                    # x_1 = c_1 - (x_0 - c_0)
                    # Para mantener los átomos adecuados,reformularemos la ecuación
                    # x_1 = - x_0 + c_0 + c_1
                    new_sp3_molecule.atoms[indexes_dict[i][1]] = -1*(new_sp3_molecule.atoms[indexes_dict[i][0]]) + new_sp3_molecule.atoms[c_0] + new_sp3_molecule.atoms[indexes_dict['c'][1]]

        return new_sp3_molecule

    def gen_transitions(self, n_transitions):
        """
        Generate the transitions of the molecule from sp3 to sp2 hybridization. The transitions are stored as .gjf files in the corresponding folder of the molecule (i.e. molecules/molecule_name/gjf_files).

        Args:
            n_transitions (int): Number of transitions to be calculated.
        """
        # Mostrar la molécula en su estado sp3
        self.sp3_molecule.show_gjf(self.header_file)

        # Generar los datos necesarios para la transición
        self.n = n_transitions
        diff_molecule = self.sp2_molecule - self.sp3_molecule
        diff_molecule /= n_transitions

        # Generar los datos necesarios para ajustar las longitudes de enlace de los puentes durante la Transición
        indexes_dictionaries = [self.mol_data["left_bridge_correspondance"], self.mol_data["right_bridge_correspondance"]]
        bridges_data = self.gen_bridges_data(indexes_dictionaries)

        # Generar cada transición ajustando las longitudes de enlace
        for n in range(1, n_transitions+1):
            self.sp3_molecule += diff_molecule
            adjusted_molecule = self.adjust_bridges(n, indexes_dictionaries, bridges_data)
            adjusted_molecule.show_gjf(self.header_file)

        # # Paso 4. Generar archivos .gjf (opcional)
        # self.sp3_molecule.generate_gjf(n_transitions, template_file_path, f"molecules/{molecule_name}/gjf_files")