# Transición
El programa calcula las coordenadas de los átomos de una molécula al transicionar a otra. El resultado es una serie de archivos .gjf que podrán ser visualizados en el programa GaussView.

## Instalación
En proceso.

## Guía de uso

### Agregar una nueva molécula
Los datos de entrada del programa los tomaremos de dos documentos .gjf que representen una molécula en el software GaussView 6.0. El primer y segundo documento representan la molécula en estado sp3 y sp2 respectivamente.

<div align="center">

| Documento en sp3 | Documento en sp2 |
| ------------ | ------------ |
| <img src="images\readme\img1.png" alt="Mi imagen" width="300"/> | <img src="images\readme\img2.png" alt="Mi imagen" width="300"/> |

</div>

Almacenaremos las coordenadas de la molécula en cada estado en 2 nuevos documentos: **sp2-coordinates.txt** y **sp3-coordinates.txt**. Se deben almacenar en la carpeta molecules\\**nombre**, donde **nombre** corresponde al nombre de la molécula. Para obtener las coordenadas adecuadas —dentro de GaussView— accedemos a **Calculate\GaussianCalculationSetup** (Ctrl+G) y presionamos la pestaña **Preview** donde encontraremos las coordenadas de todos los átomos de nuestra molécula. Dentro de los nuevos documentos deben estar únicamente las coordenadas de cada átomo.

<img src="images\readme\img3.png" alt="Mi imagen" width="400"/>

Se debe notar que al transicionar del estado sp3 a sp2 se **pierden 4 átomos** de hidrogeno debido a la deformación del puente C<sub>2</sub>H<sub>6</sub>. Para poder correr este programa se tienen que **agregar estos 4 átomos a sp2-coordinates.txt**. Las coordenadas de los 4 átomos de hidrógeno quedarán sobrepuestas con las de los 4 átomos de carbono hacia los que se dirigen durante la transición. 

En la molécula LI-7 las parejas de átomos que quedarán sobrepuestos son: (14, 11), (13, 19), (10,7) y (8, 17).

<div align="center">

| Átomos que quedarán sobrepuestos | Átomos agregados a sp2-coordinates.txt |
| ------------ | ------------ |
| <img src="images\readme\img4.png" alt="Mi imagen" width="300"/> | <img src="images\readme\img6.png" alt="Mi imagen" width="300"/> |

</div>

A continuación debemos asegurarnos que los átomos tienen el mismo índice en ambos estados. Por ejemplo, el átomo con índice 1 en sp3 debe corresponder al mismo átomo con índice 1 en sp2. Si los índices no corresponden, se añaden los 4 átomos de hidrógeno sobrepuestos a los 4 de carbono en cualquier lugar y se usa el script **molecules/sort_atoms.py** para reordenar los átomos en sp2.

<div align="center">

| Documento en sp3 | Documento en sp2 |
| ------------ | ------------ |
| <img src="images\readme\img12.png" alt="Mi imagen" width="300"/> | <img src="images\readme\img11.png" alt="Mi imagen" width="300"/> |

</div>

Después se debe agregar el documento **template.txt** a la carpeta molecules\\**nombre** corriendo el script **molecules/generate_template.py**. El documento debe tener la siguiente estructura:

```
%nprocshared=4
%mem=512MB
#p b3lyp/6-31g polar=(dcshg,cubic) cphf=rdfreq

Title Card Required

3 1
Átomo 1                1[n][0] 1[n][1] 1[n][2]        
Átomo 2                2[n][0] 2[n][1] 2[n][2]       
Átomo 3                3[n][0] 3[n][1] 3[n][2]
   .                               .
   .                               .
   .                               .
Átomo i                i[n][0] i[n][1] i[n][2]

0.0828424  0.0723227  0.05695416  0.0536039  0.042812964  0.0347811  0.0293956  



```

Donde:
- **i** es el número de átomos de la molécula.
- **Átomo i** es el elemento correspondiente. Ejemplo: C, H, N.

Además se debe agregar una carpeta llamada **gjf_files** donde se almacenarán todas las coordenadas de cada progresión de la transición. La estructura final de molecules\\**nombre** debe ser:

```
nombre/
├── gjf_files/
├── sp2-coordinates.txt
├── sp3-coordinates.txt
└── template.txt
```

Por último, se debe añadir la nueva molécula a la base de datos database/molecules.db. Para esto se debe agregar el siguiente código al archivo database/database.py reemplazando el código correspondiente antes de la línea ```c.execute("SELECT * FROM molecules")```:

```
# LI-7
c_atoms_i_sp3 = json.dumps([6, 15])
c_atoms_i_sp2 = json.dumps([6, 15])
l_cation_i = json.dumps([28, 29, 30, 31, 32, 33, 34, 36, 41, 42, 43, 44])
r_cation_i = json.dumps([21, 22, 23, 24, 25, 26, 27, 35, 37, 38, 39, 40])
l_bridge_i = json.dumps([11, 12, 13, 14, 19, 20])
r_bridge_i = json.dumps([7, 8, 9, 10, 17, 18])
c.execute("""INSERT INTO molecules VALUES ('LI-7', 5, ?, 5, ?, ?, ?, ?, ?)""", (c_atoms_i_sp3, c_atoms_i_sp2, l_cation_i, r_cation_i, l_bridge_i, r_bridge_i))
```

Este código añade añande una nueva fila a la tabla molecules con los datos de la moléucla LI-7. A continución se explica cada línea del código:

1. ```c_atoms_i_sp3 = json.dumps([6, 15])```: la función **json.dumps()** convierte una lista de python a una cadena de texto en formato JSON. En este caso, los índices 6 y 15 corresponden a los átomos w<sub>1</sub> y w<sub>2</sub> en el estado sp3.

<div align="center">
<img src="images\readme\img8.png" alt="Mi imagen" width="300"/>
</div>

2. ```c_atoms_i_sp2 = json.dumps([6, 15])```: similar a la línea anterior, pero para el estado sp2.

3. ```l_cation_i = json.dumps([28, 29, 30, 31, 32, 33, 34, 36, 41, 42, 43, 44])```: Contiene los índices de los átomos que conforman el catión izquierdo de la molécula LI-7.

<div align="center">
<img src="images\readme\img9.png" alt="Mi imagen" width="300"/>
</div>

4. ```r_cation_i = json.dumps([21, 22, 23, 24, 25, 26, 27, 35, 37, 38, 39, 40])```: Contiene los índices de los átomos que conforman el catión derecho de la molécula LI-7.

5. ```l_bridge_i = json.dumps([11, 12, 13, 14, 19, 20])```: Contiene los índices de los átomos que conforman el puente C<sub>2</sub>H<sub>6</sub> izquierdo de la molécula LI-7.

<div align="center">
<img src="images\readme\img10.png" alt="Mi imagen" width="300"/>
</div>

6. ```r_bridge_i = json.dumps([7, 8, 9, 10, 17, 18])```: Contiene los índices de los átomos que conforman el puente C<sub>2</sub>H<sub>6</sub> derecho de la molécula LI-7.

7. ```c.execute("""INSERT INTO molecules VALUES ('LI-7', 5, ?, 5, ?, ?, ?, ?, ?)""", (c_atoms_i_sp3, c_atoms_i_sp2, l_cation_i, r_cation_i, l_bridge_i, r_bridge_i))```: esta línea inserta una nueva fila en la tabla molecules con los datos de la molécula LI-7. Los 3 valores: LI-7, 5 y 5 corresponden al nombre de la molécula, el átomo central en sp3 y átomo central sp2; los signos de interrogación (?) son marcadores de posición que serán reemplazados por los valores proporcionados en la tupla al final de la línea.

## Base de datos
Esta base de datos almacena información que el programa usará para generar las transiciones. La base de datos se encuentra en database/molecules.db y sus columnas son:

   - **name** (TEXT): Nombre de la molécula. Debe coincidir con el nombre de la carpeta en molecules.
   - **origin_atom_i_sp3** (INTEGER): Índice del átomo que se llevará a las coordenadas [0,0,0] en el estado sp3. 
   - **central_atoms_i_sp3** (TEXT): índices de los 2 átomos w<sub>1</sub> y w<sub>2</sub>, los cuales se encuentran enla estructura pentagonal central de la molécula. Se puede ver más información sobre estos átomos en notebooks/**align_to_xy.ipynb**.
   - **origin_atom_i_sp2** (INTEGER): Índice del átomo que se llevará a las coordenadas [0,0,0] en el estado sp2.
   - **central_atoms_i_sp2** (TEXT): índices de los 2 átomos w<sub>1</sub> y w<sub>2</sub> en el estado sp2.
   - **left_cation_i** (TEXT): índices de los átomos que conforman el catión izquierdo.
   - **right_cation_i** (TEXT): índices de los átomos que conforman el catión derecho.
   - **left_bridge_i** (TEXT): índices de los átomos que conforman el puente C<sub>2</sub>H<sub>6</sub> izquierdo.
   - **right_bridge_i** (TEXT): índices de los átomos que conforman el puente C<sub>2</sub>H<sub>6</sub> derecho.

## Licensia
MIT LICENSE

