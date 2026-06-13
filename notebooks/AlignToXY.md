# ¿Cómo alinear las moleculas al plano xy?

Como se puede apreciar en esta imágen, las moléculas pueden llegar a estar desalineadas con el plano $x y$. Para correr este programa se tuvo que idear el método **get_alignment_matrix()** perteneciente a la clase **molecule**, de esa manera encontramos una transformación lineal que ubica el estado sp3 y sp2 en la misma posición al alinearlos en el plano xy.

| Molécula no alineada al plano xy | Molécula no alineada al plano xy |
|----------------------------------|----------------------------------|
|![image.png](attachment:image.png)|![image-2.png](attachment:image-2.png)|

La función **get_alignment_matrix()** debe encontrar una transformación lineal representada por una matriz $A$ que al ser aplicada al conjunto de vectores $M_{desalineada}$ (contiene todos los átomos de la molécula representados en vectores) alinea la molécula.

$$ M_{desalineada} = {v_{1}, v_{2}, ..., v_{n}} $$

$$ A(M_{desalineada}) = {Av_{1}, Av_{2}, ..., Av_{n}} = M_{alineada}$$

Utilizaremos los vectores unitarios alineados:

$$A\hat{\textbf{\i}} = [1, 0, 0]$$

$$A\hat{\textbf{\j}} = [0, 1, 0]$$

$$A\hat{\textbf{k}} = [0, 0, 1]$$

Esto nos será de utilidad cuando hayamos encontrado los vectores $U_{desalineados} = ( \hat{\textbf{\i}}$, $\hat{\textbf{\j}}, \hat{\textbf{k}} )$ pertenecientes a $M_{desalineada}$.

Ya que:

$$ U_{alineados} =
\begin{bmatrix}
1 & 0 & 0 \\
0 & 1 & 0 \\
0 & 0 & 1
\end{bmatrix}
$$

$$ U_{alineados} = AU_{desalineados} $$

$$ A^{-1}U_{alineados} = U_{desalineados} $$

$$ A^{-1} = U_{desalineados} $$

$$ A^{-1} = 
\begin{bmatrix}
\hat{\textbf{\i}} &
\hat{\textbf{\j}} &
\hat{\textbf{k}}
\end{bmatrix}
$$

podemos usar $A^{-1}$ para encontrar la transformación lineal que buscamos: $A = \begin{bmatrix} \hat{\textbf{\i}} & \hat{\textbf{\j}} & \hat{\textbf{k}} \end{bmatrix}^{-1}$ 

## ¿Cómo encontrar la matriz A de forma computacional?

Al analizar el catión central de las moléculas, observaremos el carbono superior de la estructura pentagonal. De esta manera podemos representarlos en un plano $xy$ dándonos a conocer su posición tras la transformación lineal $A$.

![image-2.png](attachment:image-2.png)

Podemos deducir que:

$$
\begin{align*}
\hat{\textbf{\i}} &= \frac{w_{2} - r}{||w_{2} - r||} \\
\hat{\textbf{\j}} &= \frac{w_1}{||w_1||} \\
\hat{\textbf{k}} &= \hat{\textbf{\i}} \times \hat{\textbf{\j}}
\end{align*}
$$

Donde:

$$ 
\begin{align*}
r &= -\hat{\textbf{\j}} \cdot ||w_{2} \cdot sin\theta|| \\
\theta &= 35.5\degree
\end{align*}
$$

Podemos sustituir para representar A en función de $w_1$, $w_2$ y $\theta$:

$$
\begin{align*}
\hat{\textbf{\i}} &= \frac{w_{2} + \frac{w_1}{||w_1||} \cdot ||w_{2} \cdot sin\theta||}{|| w_{2} + (\frac{w_1}{||w_1||} \cdot ||w_{2} \cdot sin\theta||) ||} \\
\hat{\textbf{\j}} &= \frac{w_1}{||w_1||} \\
\hat{\textbf{k}} &= \hat{\textbf{\i}} \times \hat{\textbf{\j}} \\
A &=
\begin{bmatrix} 
\hat{\textbf{\i}} &
\hat{\textbf{\j}} &
\hat{\textbf{k}}
\end{bmatrix}^{-1}
\end{align*}
$$