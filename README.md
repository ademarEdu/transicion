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

Por úlimo se debe agregar el documento **template.txt** a la carpeta molecules\\**nombre**, el cual debe tener la siguiente estructura:

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

Por último, se deben añadir datos importantes sobre la molécula a la base de datos ubicada en database/molecules.db. Para esto 

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

