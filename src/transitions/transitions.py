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

def get_transitions(molecule_name, n_transitions):
    """
    Calculate the coordinates of each transition.

    Args:
        molecule_name (str): Name of the molecule to be studied. It has to be the same as the folder name in the molecules directory.
        n_transitions (int): Number of transitions to be calculated.
        gjf (bool): If True, .gjf files will be generated in the specified folder.
    """
    # Obtener los datos necesarios para calcular las transiciones
    data = get_data(molecule_name)
    
    # Rutas de los archivos que contienen las coordenadas de los átomos de las moléculas
    sp3_path = f"molecules/{molecule_name}/sp3-coordinates.txt"
    sp2_path = f"molecules/{molecule_name}/sp2-coordinates.txt"
    header_file = f"molecules/header.txt"

    left_cation_origin_index = data["left_cation_indexes"][0]
    right_cation_origin_index = data["right_cation_indexes"][0]

    # Inicializar las moléculas
    sp3_molecule = TriCationicMolecule(
        name = molecule_name+"_sp3",
        atoms_dict = utils.get_atoms_dict(sp3_path),
        origin_atom_index = data["origin_atom_index_sp3"],
        central_atoms_indexes= data["central_atoms_indexes_sp3"],
        left_cation_indexes= data["left_cation_indexes"],
        right_cation_indexes= data["right_cation_indexes"],
        left_bridge_indexes= data["left_bridge_indexes"],
        right_bridge_indexes= data["right_bridge_indexes"]
    )

    sp2_molecule = TriCationicMolecule(
        name = molecule_name+"_sp2",
        atoms_dict = utils.get_atoms_dict(sp2_path),
        origin_atom_index = data["origin_atom_index_sp2"],
        central_atoms_indexes= data["central_atoms_indexes_sp2"],
        left_cation_indexes= data["left_cation_indexes"],
        right_cation_indexes= data["right_cation_indexes"],
        left_bridge_indexes= data["left_bridge_indexes"],
        right_bridge_indexes= data["right_bridge_indexes"]
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

    # Generar transiciones
    diff_molecule = sp2_molecule - sp3_molecule
    diff_molecule /= n_transitions

    # Calcular datos necesarios para ajustar longitudes de enlace
    indexes_dictionaries = [data["left_bridge_correspondance"], data["right_bridge_correspondance"]]
    bridges_data = {'left': {}, 'right': {}}

    # c: el átomo de carbono del puente
    # a: el átomo de hidrógeno que se quedará al transicionar a sp2
    # b: el átomo de hidrógeno que se moverá para dejar espacio al doble enlace

    for i in ['c', 'a', 'b']:
        for indexes_dict, side in zip(indexes_dictionaries, ['left', 'right']):
            if i == 'c':
                v1 = sp3_molecule.atoms[indexes_dict['c'][1]] - sp3_molecule.atoms[indexes_dict['c'][0]]
                r = v1.magnitude()
                r_prime = (sp2_molecule.atoms[indexes_dict['c'][1]] - sp2_molecule.atoms[indexes_dict['c'][0]]).magnitude()
                delta_r = (r_prime - r)/n_transitions
                bridges_data[side][i] = {
                    'v1': v1,
                    'r': r,
                    'r_prime': r_prime,
                    'delta_r': delta_r
                }
            else:
                v1 = sp3_molecule.atoms[indexes_dict[i][0]] - sp3_molecule.atoms[indexes_dict['c'][0]]
                r = v1.magnitude()
                r_prime = (sp2_molecule.atoms[indexes_dict[i][0]] - sp2_molecule.atoms[indexes_dict['c'][0]]).magnitude()
                delta_r = (r_prime - r)/n_transitions
                bridges_data[side][i] = {
                    'v1': v1,
                    'r': r,
                    'r_prime': r_prime,
                    'delta_r': delta_r
                }

    # Generar cada transición ajustando las longitudes de enlace
    for n in range(1, n_transitions+1):
        sp3_molecule += diff_molecule
        new_sp3_molecule = sp3_molecule.__copy__()
        for i in ['c', 'a', 'b']:
            for side, indexes_dict in zip(['left', 'right'], indexes_dictionaries):
                if i == 'c':
                    v2 = new_sp3_molecule.atoms[indexes_dict['c'][1]] - new_sp3_molecule.atoms[indexes_dict['c'][0]]
                    v3 = ((bridges_data[side]['c']['r'] - v2.magnitude() + bridges_data[side]['c']['delta_r'] * n ) / v2.magnitude()) * v2
                    new_sp3_molecule.atoms[indexes_dict['c'][1]] += v3
                else:
                    v2 = new_sp3_molecule.atoms[indexes_dict[i][0]] - new_sp3_molecule.atoms[indexes_dict['c'][0]]
                    v3 = ((bridges_data[side][i]['r'] - v2.magnitude() + bridges_data[side][i]['delta_r'] * n ) / v2.magnitude()) * v2
                    new_sp3_molecule.atoms[indexes_dict[i][0]] += v3
                    new_sp3_molecule.atoms[indexes_dict[i][1]] = (-1)*(new_sp3_molecule.atoms[indexes_dict[i][0]]) + new_sp3_molecule.atoms[indexes_dict['c'][1]] + new_sp3_molecule.atoms[indexes_dict['c'][0]]
        new_sp3_molecule.show_gjf(header_file)

    # # Paso 4. Generar archivos .gjf (opcional)
    # sp3_molecule.generate_gjf(n_transitions, template_file_path, f"molecules/{molecule_name}/gjf_files")