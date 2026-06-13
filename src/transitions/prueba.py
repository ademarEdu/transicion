# Paso 1.- Obtener los datos
import pandas as pd
import numpy as np

with open('../../molecules/LI-7-A/sp3-coordinates.txt', 'r') as f:
    lines = [' '.join(line.split()) for line in f if line.strip()]  # removes extra spaces and empty lines

# Then read the cleaned content
from io import StringIO
df = pd.read_csv(StringIO('\n'.join(lines)), delimiter=' ', names=['Atom', 'X', 'Y', 'Z'])

from molecule import Atom, Molecule
atoms_dict = {i+1:Atom(df['Atom'][i], np.array([df['X'][i], df['Y'][i], df['Z'][i]])) for i in df.index}

# -----------------------------------------------------

# Crear la molécula
molecule = Molecule(name='LI-7-A', atoms_dict=atoms_dict, origin_atom_index=5)
header_file = f"../../molecules/header.txt" # ruta que se usará para crear el archivo .gjf

# Mover el origen al centro de la molécula
c = molecule.atoms[23].magnitude() # distancia del átomo 23 al origen
h = c / (2*np.sin(np.radians(36))) # altura del triángulo equilátero formado por los átomos 5, 23 y 24
print("el valor de h es:", h)

# Mover todos los átomos
for atom in molecule.atoms.values():
    atom.coordinates += np.array([0, h, 0]) # mover la molécula hacia arriba para que el átomo 23 quede en el plano XY

# molecule.show_gjf(header_file)

# ínices de los aniones
upper_anion_i = [25]
left_anion_i = [16, 18, 21, 22, 46, 47, 48, 49]
right_anion_i = [15, 17, 19, 20, 50, 51, 52, 53]

# Vectores que definen la dirección
a = molecule.atoms[25] / molecule.atoms[25].magnitude() # anión superior
b = molecule.atoms[18] / molecule.atoms[18].magnitude() # anión izquierdo
c = molecule.atoms[17] / molecule.atoms[17].magnitude() # anión derecho

# Distancias máxima del vector en amstrongs
paso_amstrongs = 0.25

# Mover los aniones
for j in [-2, -1, 0, 1, 2, 3, 4, 5]:
    new_molecule = molecule.__copy__()
    for i in upper_anion_i:
        new_molecule.atoms[i] += a * j * paso_amstrongs
    for i in left_anion_i:
        new_molecule.atoms[i] += b * j * paso_amstrongs
    for i in right_anion_i:
        new_molecule.atoms[i] += c * j * paso_amstrongs
    print(f"The magnitude of atom[25] at step {j} is: {new_molecule.atoms[25].magnitude()}")
    print(f"The magnitude of atom[18] at step {j} is: {new_molecule.atoms[18].magnitude()}")
    print(f"The magnitude of atom[17] at step {j} is: {new_molecule.atoms[17].magnitude()}")
    new_molecule.save_gjf(header_file=header_file, new_file_path=f"../../prueba/{j+3}transition.gjf")
    print(f"The {j} file was generated.")