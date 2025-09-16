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
origin_atom_index = 5 # Editar este valor para cambiar el átomo que se llevará al origen
central_atoms_indexes = [6, 15] # Editar este valor para cambiar los índices de los átomos centrales

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

    sp3_molecule = TriCationicMolecule(
        name = molecule_name+"_sp3",
        atoms_list = sp3_atoms_list,
        origin_atom_index = origin_atom_index,
        central_atoms_indexes= central_atoms_indexes
    )

    sp2_molecule = TriCationicMolecule(
        name = molecule_name+"_sp2",
        atoms_list = sp2_atoms_list,
        origin_atom_index = origin_atom_index,
        central_atoms_indexes= central_atoms_indexes
    )
    
    # Generar Transiciones
    get_all_atoms = lambda x: np.array(list(x.atoms.values()))
    
    # Paso 2.- Obtener la diferencia entre las coordenadas de los átomos de la molécula sp3 y sp2
    diffvecatomslist = get_all_atoms(sp2_molecule) - get_all_atoms(sp3_molecule)
    sp3_molecule.transitions = np.zeros((n_transitions, sp3_molecule.number_atoms, 3))

    # Paso 3. Calcular las coordenadas de cada transición
    # El objetivo de todos los algoritmos debe de ser:
    #   - Hacer los cambios de ángulo en cada transición sumando difvecatomslist a las coordenadas de los átomos de la molécula sp3.
    #   - Tomar en cuenta los cambios de distancia de enlace que ocurren durante las transiciones.
    for n in range(1, n_transitions+1):
        sp3_molecule.transitions[n-1] = get_all_atoms(sp3_molecule) + (diffvecatomslist * n/(n_transitions))

    # Paso 4. Generar archivos .gjf (opcional)
    if gjf:
        sp3_molecule.generate_gjf(n_transitions, template_file_path)