# Especificar la ruta de la molécula
path = "./LI-7-B-D-F1"

def generate_template(molecule_path):
    """
    Adjusts the template.txt file values.

    The molecule directory has to contain the sp3 and sp2-coordinates.txt files before running this function.

    Args:
        path (str): The path where the template.txt file will be created.
    """
    # Ruta del archivo sp3-coordinates.txt
    sp3_path = molecule_path + "/sp3-coordinates.txt"
    # Almacenar los renglones totalmente nuevos del archivo template.txt
    new_lines = [
        "%nprocshared=4",
        "%mem=512MB",
        "#p b3lyp/6-31g polar=(dcshg,cubic) cphf=rdfreq",
        "",
        "Title Card Required",
        "",
        "3 1",
        "",
        "0.0828424  0.0723227  0.05695416  0.0536039  0.042812964  0.0347811  0.0293956  ",
        "",
        "",
        ""
    ]

    # Leer las líneas del archivo sp3-coordinates.txt
    with open(sp3_path, mode="r", encoding="utf-8") as file:
        lines = file.readlines()

    # Modificar las líneas según sea necesario
    i = 0
    for line in lines:
        lines[i] = line[:1] +  f"                 {i+1}[n][0]   {i+1}[n][1]    {i+1}[n][2]"
        i += 1
        
    # Insertar las nuevas líneas al principio y al final
    lines = new_lines[:7] + lines + new_lines[7:]
    # Crear un nuevo archivo template.txt con las líneas modificadas
    with open(molecule_path + "/template.txt", mode="w", encoding="utf-8") as file:
        file.writelines("\n".join(lines))

# Generar el archivo template.txt en la ruta especificada
generate_template(path)

