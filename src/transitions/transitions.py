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
        "right_bridge_indexes": json.loads(data[8])
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

    # sp2_molecule.show_gjf(header_file)
    # sp3_molecule.save_gjf(header_file, f"molecules/{molecule_name}/gjf_files/trial.gjf")

    # sp3_molecule.left_bridge.show_gjf(header_file)
    # sp2_molecule.left_bridge.show_gjf(header_file)

    # Generar transiciones
    diff_molecule = sp2_molecule - sp3_molecule
    diff_molecule /= n_transitions
    for i in range(1, n_transitions+1):
        sp3_molecule += diff_molecule
        sp3_molecule.show_gjf(header_file)

    # Generar Transiciones
    # Paso 3. Calcular las coordenadas de cada transición
    # El objetivo de todos los algoritmos debe de ser:
    #   - Hacer los cambios de ángulo en cada transición sumando difvecatomslist a las coordenadas de los átomos de la molécula sp3.
    #   - Tomar en cuenta los cambios de distancia de enlace que ocurren durante las transiciones.
    
    # # Iterar para cada transición
    # # Se procesa cada catión por separado
    # for n in range(1, n_transitions+1):
    #     # Iterar para cada átomo
    #     for i in sp3_molecule.atoms:
    #         if i in sp3_molecule.left_cation.atoms.keys():
    #             sp3_molecule.transitions[n-1][i-1] = sp3_molecule.atoms[i] + (diff_left_cation[i] * n/(n_transitions))
    #         elif i in sp3_molecule.right_cation.atoms.keys():
    #             sp3_molecule.transitions[n-1][i-1] = sp3_molecule.atoms[i] + (diff_right_cation[i] * n/(n_transitions))
    #         else:
    #             sp3_molecule.transitions[n-1][i-1] = sp3_molecule.atoms[i] + (diff_molecule[i] * n/(n_transitions))
    #     # Pegar los cationes a su átomo origen transicionando
    #     for i in [left_cation_origin_index, right_cation_origin_index]:    
    #         if i == left_cation_origin_index:
    #             for j in sp3_molecule.left_cation.atoms.keys():
    #                 sp3_molecule.transitions[n-1][j-1] += diff_molecule[i] * n/(n_transitions)
    #             sp3_molecule.transitions[n-1][i-1] = sp3_molecule.atoms[i] + (diff_molecule[i] * n/(n_transitions))
    #         elif i == right_cation_origin_index:
    #             for j in sp3_molecule.right_cation.atoms.keys():
    #                 sp3_molecule.transitions[n-1][j-1] += diff_molecule[i] * n/(n_transitions)
    #             sp3_molecule.transitions[n-1][i-1] = sp3_molecule.atoms[i] + (diff_molecule[i] * n/(n_transitions))    

    # # Paso 4. Generar archivos .gjf (opcional)
    # sp3_molecule.generate_gjf(n_transitions, template_file_path, f"molecules/{molecule_name}/gjf_files")