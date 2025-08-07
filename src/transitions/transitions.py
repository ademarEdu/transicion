import numpy as np
import transitions.utilities as utils

# Datos necesarios para calcular las transiciones
# Para más información sobre qué representan estos datos, consulta el archivo README.md

# Rutas de los archivos que contienen las coordenadas de los átomos de las moléculas
sp3_path = "molecules/LI-7/sp2-coordinates.txt"
sp2_path = "molecules/LI-7/sp3-coordinates.txt"
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
    # Alinear las moléculas al eje z
    sp3_atoms_list = utils.align_molecule_to_z_axis(sp3_atoms_list, sp3_atoms_list[origin_atom_index])
    sp2_atoms_list = utils.align_molecule_to_z_axis(sp2_atoms_list, sp2_atoms_list[origin_atom_index])
   
    # ----------------------------------------------
    # Paso 2. Obtener diferencias entre las coordenadas de los átomos de la molécula sp3 y sp2

    difvecatomslist = sp2_atoms_list-sp3_atoms_list
    difvecatomslist = difvecatomslist/n_transitions-1 # Se resta 1 para ajustar la suma de 1 del principio

    # ----------------------------------------------
    # Paso 3. Calcular las coordenadas de cada transición (Iniciar el algoritmo)
    # El objetivo de todos los algoritmos debe de ser:
    #   - Hacer los cambios de ángulo en cada transición sumando difvecatomslist a las coordenadas de los átomos de la molécula sp3.
    #   - Tomar en cuenta los cambios de distancia de enlace que ocurren durante las transiciones.
    
    transitions = algorithm.run_algorithm(n_transitions, sp3_atoms_list, sp2_atoms_list, difvecatomslist)

    # ----------------------------------------------
    # Paso 4. Generar archivos .gjf (opcional)
    if gjf:
        utils.generate_gjf(transitions, n_transitions, len(sp3_atoms_list), template_file_path)
