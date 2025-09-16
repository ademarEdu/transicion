import numpy as np
import transitions.utilities as utils
from transitions.molecule import TriCationicMolecule

# Datos necesarios para calcular las transiciones
# Para más información sobre qué representan estos datos, consulta el archivo README.md

# Rutas de los archivos que contienen las coordenadas de los átomos de las moléculas
molecule_name = "LI-7" # Editar este valor para cambiar la molécula a estudiar
sp3_path = "molecules/LI-7/sp3-coordinates.txt"
sp2_path = "molecules/LI-7/sp2-coordinates.txt"
template_file_path = "molecules/LI-7/template.txt"

# Este es el átomo que se llevará a la coordenada [0, 0, 0]
origin_atom_index = 11 # Editar este valor para cambiar el átomo que se llevará al origen
origin_atom_index -= 1

# Número de transiciones que se calcularán
# Para más información sobre qué representa una transición, consulta el archivo README.md
n_transitions = 3 # Editar este valor para cambiar el número de transiciones que se calcularán
n_transitions += 1 # Se suma 1 para incluir la transición final (molécula sp2)

# Generar archivos .gjf
gjf = True # Editar este valor para generar los archivos .gjf o no

# Algoritmo a utilizar
import transitions.algorithms.LI_7 as algorithm # Editar este valor si se quiere calcular las transiciones de otra molecula

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

    # sp3__molecule = TriCationicMolecule(
    #     name = molecule_name+"_sp3",
    #     atoms_list = sp3_atoms_list,
    #     origin_atom_index = 5-1,
    #     central_atoms_indexes=[i-1 for i in [6, 15]]
    # )

    # sp2_molecule = TriCationicMolecule(
    #     name = molecule_name+"_sp2",
    #     atoms_list = sp2_atoms_list,
    #     origin_atom_index = 5-1,
    #     central_atoms_indexes=[i-1 for i in [6, 15]]
    # )
    
    molecule = TriCationicMolecule(
        origin_atom_index= 5-1,
        central_atoms_indexes=[i-1 for i in [6, 15]],
        c_cation_indexes=[i-1 for i in [1, 2, 3, 4, 5, 6, 15, 16]],
        l_cation_indexes=[i-i for i in [28, 29, 30, 31, 32, 33, 34, 36, 41, 42, 43, 44]],
        r_cation_indexes=[i-1 for i in [21, 22, 23, 24, 25, 26, 27, 35, 37, 38, 39, 40]],
        l_bridge_indexes=[i-1 for i in [11, 12, 13, 14, 19, 20]],
        r_bridge_indexes=[i-1 for i in [7, 8, 9, 10, 17, 18]],
        sp3_coords=sp3_atoms_list,
        sp2_coords=sp2_atoms_list
    )

    molecule.generate_transitions(n_transitions)

    molecule.generate_gjf(n_transitions, template_file_path)

    # # ----------------------------------------------
    # # Paso 2. Obtener diferencias entre las coordenadas de los átomos de la molécula sp3 y sp2

    # difvecatomslist = sp2_atoms_list-sp3_atoms_list
    # difvecatomslist = difvecatomslist/n_transitions-1 # Se resta 1 para ajustar la suma de 1 del principio

    # # ----------------------------------------------
    # # Paso 3. Calcular las coordenadas de cada transición (Iniciar el algoritmo)
    # # El objetivo de todos los algoritmos debe de ser:
    # #   - Hacer los cambios de ángulo en cada transición sumando difvecatomslist a las coordenadas de los átomos de la molécula sp3.
    # #   - Tomar en cuenta los cambios de distancia de enlace que ocurren durante las transiciones.
    
    # transitions = algorithm.run_algorithm(n_transitions, sp3_atoms_list, sp2_atoms_list, difvecatomslist)

    # # ----------------------------------------------
    # # Paso 4. Generar archivos .gjf (opcional)
    # if gjf:
    #     utils.generate_gjf(transitions, n_transitions, len(sp3_atoms_list), template_file_path)
