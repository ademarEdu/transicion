from turtle import left
import numpy as np
import transitions.utilities as utils
from transitions.molecule import Molecule, TriCationicMolecule

# Datos necesarios para calcular las transiciones
# Para más información sobre qué representan estos datos, consulta el archivo README.md

# Rutas de los archivos que contienen las coordenadas de los átomos de las moléculas
molecule_name = "LI-7" # Editar este valor para cambiar la molécula a estudiar
sp3_path = f"molecules/{molecule_name}/sp3-coordinates.txt"
sp2_path = f"molecules/{molecule_name}/sp2-coordinates.txt"
template_file_path = f"molecules/{molecule_name}/template.txt"

# Este es el átomo que se llevará a la coordenada [0, 0, 0]
origin_atom_index_sp3 = 5 # Editar este valor para cambiar el átomo que se llevará al origen
central_atoms_indexes_sp3 = [6, 15] # Editar este valor para cambiar los índices de los átomos centrales
origin_atom_index_sp2 = 5
central_atoms_indexes_sp2 = [6, 15]

# Índices de las partes de la molécula
left_cation_indexes = [28, 29, 30, 31, 32, 33, 34, 36, 41, 42, 43, 44]
right_cation_indexes = [21, 22, 23, 24, 25, 26, 27, 35, 37, 38, 39, 40]
left_cation_origin_index = left_cation_indexes[0]
right_cation_origin_index = right_cation_indexes[0]
left_bridge_indexes = [11, 12, 13, 14, 19, 20]
right_bridge_indexes = [7, 8, 9, 10, 17, 18]

# Número de transiciones que se calcularán
# Para más información sobre qué representa una transición, consulta el archivo README.md
n_transitions = 4 # Editar este valor para cambiar el número de transiciones que se calcularán

# Generar archivos .gjf
gjf = True # Editar este valor para generar los archivos .gjf o no

# --------------------------------------------------

def get_transitions():
    """
    Calculate the coordinates of each transition.
    """
    # ----------------------------------------------
    # Paso 1. Obtener coordenadas y alinear las moleculas

    # Obtener la lista de coordenadas de los átomos de la molécula
    sp3_atoms_list = utils.parse_atom_coordinates(sp3_path)
    sp2_atoms_list = utils.parse_atom_coordinates(sp2_path)

    sp3_atoms_dict = {i:atom for i, atom in zip(range(1, len(sp3_atoms_list)+1), sp3_atoms_list)}
    sp2_atoms_dict = {i:atom for i, atom in zip(range(1, len(sp2_atoms_list)+1), sp2_atoms_list)}

    # Inicializar las moléculas
    sp3_molecule = TriCationicMolecule(
        name = molecule_name+"_sp3",
        atoms_dict = sp3_atoms_dict,
        origin_atom_index = origin_atom_index_sp3,
        central_atoms_indexes= central_atoms_indexes_sp3,
        left_cation_indexes= left_cation_indexes,
        right_cation_indexes= right_cation_indexes
    )

    sp2_molecule = TriCationicMolecule(
        name = molecule_name+"_sp2",
        atoms_dict = sp2_atoms_dict,
        origin_atom_index = origin_atom_index_sp2,
        central_atoms_indexes= central_atoms_indexes_sp2,
        left_cation_indexes= left_cation_indexes,
        right_cation_indexes= right_cation_indexes
    )

    # sp3_molecule.left_bridge = {i:sp3_molecule.atoms[i] for i in left_bridge_indexes}
    # sp3_molecule.right_bridge = {i:sp3_molecule.atoms[i] for i in right_bridge_indexes}

    # Generar Transiciones
    get_all_atoms = lambda x: np.array(list(x.atoms.values()))
    # Paso 2.- Obtener la diferencia entre las coordenadas de los átomos de la molécula sp3 y sp2
    diff_right_cation = sp2_molecule.right_cation - sp3_molecule.right_cation
    diff_left_cation = sp2_molecule.left_cation - sp3_molecule.left_cation
    diff_molecule = sp2_molecule - sp3_molecule
    
    sp3_molecule.transitions = np.zeros((n_transitions, sp3_molecule.number_atoms, 3))

    # Paso 3. Calcular las coordenadas de cada transición
    # El objetivo de todos los algoritmos debe de ser:
    #   - Hacer los cambios de ángulo en cada transición sumando difvecatomslist a las coordenadas de los átomos de la molécula sp3.
    #   - Tomar en cuenta los cambios de distancia de enlace que ocurren durante las transiciones.
    
    # Iterar para cada transición
    for n in range(1, n_transitions+1):
        # Iterar para cada átomo
        for i in sp3_molecule.atoms:
            if i in sp3_molecule.left_cation.atoms.keys():
                sp3_molecule.transitions[n-1][i-1] = sp3_molecule.atoms[i] + (diff_left_cation[i] * n/(n_transitions))
            elif i in sp3_molecule.right_cation.atoms.keys():
                sp3_molecule.transitions[n-1][i-1] = sp3_molecule.atoms[i] + (diff_right_cation[i] * n/(n_transitions))
            else:
                sp3_molecule.transitions[n-1][i-1] = sp3_molecule.atoms[i] + (diff_molecule[i] * n/(n_transitions))
        # Pegar los cationes a su átomo origen transicionando
        for i in [left_cation_origin_index, right_cation_origin_index]:    
            if i == left_cation_origin_index:
                for j in sp3_molecule.left_cation.atoms.keys():
                    sp3_molecule.transitions[n-1][j-1] += diff_molecule[i] * n/(n_transitions)
                sp3_molecule.transitions[n-1][i-1] = sp3_molecule.atoms[i] + (diff_molecule[i] * n/(n_transitions))
            elif i == right_cation_origin_index:
                for j in sp3_molecule.right_cation.atoms.keys():
                    sp3_molecule.transitions[n-1][j-1] += diff_molecule[i] * n/(n_transitions)
                sp3_molecule.transitions[n-1][i-1] = sp3_molecule.atoms[i] + (diff_molecule[i] * n/(n_transitions))

    # Paso 4. Generar archivos .gjf (opcional)
    if gjf:
        sp3_molecule.generate_gjf(n_transitions, template_file_path, f"molecules/{molecule_name}/gjf_files")