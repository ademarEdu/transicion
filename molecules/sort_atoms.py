# Este script se debe correr si los átomos no corresponden en el estado sp2 y sp3

# Especificar la ruta de la molécula
# path = "./LI-7-B-D-F"
# correspondance = {
#     1:1,
#     2:2,
#     3:3,
#     4:4,
#     5:5,
#     6:6,
#     7:7,
#     8:18,
#     9:9,
#     10:8,
#     11:11,
#     12:20,
#     13:13,
#     14:14,
#     15:15,
#     16:16,
#     17:17,
#     18:10,
#     19:19,
#     20:12,
#     21:21,
#     22:22,
#     23:23,
#     24:24,
#     25:25,
#     26:26,
#     27:27,
#     28:28,
#     29:29,
#     30:30,
#     31:31,
#     32:32,
#     33:33,
#     34:34,
#     35:35,
#     36:36,
#     37:37,
#     38:38,
#     39:39,
#     40:40,
#     41:41,
#     42:42,
#     43:47,
#     44:50,
#     45:48,
#     46:49,
#     47:43,
#     48:45,
#     49:46,
#     50:44,
#     51:51,
#     52:52,
#     53:53,
#     54:54,
#     55:55,
#     56:56,
#     57:57,
#     58:58
# }

# Especificar la ruta de la molécula
path = "./LI-7-G"
correspondance = {
    1:1,
    2:2,
    3:3,
    4:4,
    5:5,
    6:6,
    7:7,
    8:16,
    9:9,
    10:8,
    11:10,
    12:12,
    13:19,
    14:11,
    15:13,
    16:14,
    17:15,
    18:17,
    19:18,
    20:20,
    21:21,
    22:22,
    23:23,
    24:24,
    25:25,
    26:26,
    27:27,
    28:28,
    29:29,
    30:30,
    31:31,
    32:32,
    33:33,
    34:34,
    35:35,
    36:36,
    37:37,
    38:38,
    39:39,
    40:40,
    41:41,
    42:42,
    43:43,
    44:44,
    45:45,
    46:46,
    47:47,
    48:48,
    49:49,
    50:50
}

def sort_atoms(path, correspondance_d):
    """
    This function sorts the atoms of the molecule in the sp2 state given its correspondace with the sp3 state. The sorting happens automatically in the file of the molecule in sp2 state.

    Args:
        path (str): The path to the desired molecule directory.
        correspondance_d (dict): A distionary that contains the indexes of the atoms in sp3 as keys and the corresponding indexes in sp2 as values.
    """
    sp2_path = path + "/sp2-coordinates.txt"
    with open(sp2_path, mode="r", encoding="utf-8") as file:
        lines = file.readlines()
    new_lines = [""] * len(lines)
    for atom in correspondance_d:
        new_lines[atom - 1] = lines[correspondance_d[atom] - 1]
    with open(sp2_path, mode="w", encoding="utf-8") as file:
        file.writelines(new_lines)

    print(f"Atoms of {path + '/sp2-coordinates.txt'} sorted successfully.")
        
sort_atoms(path, correspondance)