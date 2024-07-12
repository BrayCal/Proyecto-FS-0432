# Planteamiento del problema

La dinámica de un estado puro para un sistema cuántico aislado se rige bajo la ecuación Schrödinger ($\hbar = 1$)

\begin{equation}
    \frac{\partial \ket{\psi(t)}}{\partial t} = -i \hat{H} \ket{\psi(t)},
\label{Ecuación 1}
\end{equation}

cuya solución formal está dada por

\begin{equation}
    \ket{\psi(t)} = e^{-i\hat{H}(t-t_o)} \ket{\psi(t=t_o)}.
\label{Ecuación }
\end{equation}

Es decir, la solución involucra resolver de manera numérica la ecuación diferencial o evaluar de alguna forma la exponencial de la matriz. La idea del proyecto es evaluar la dinámica del modelo de Ising empezando de algún estado inicial.

\begin{equation}
    \hat{H} = -J \sum_{i=1}^N \hat{\sigma}_i ^z \hat{\sigma}_{i+1} ^z -g \sum_{i=1}^N \hat{\sigma}_i ^x,
\label{Hamiltoniano}
\end{equation}

donde $J$ es una escala energética que determina la interacción ferromagnética, $g$ el parámetro energético del campo transversal y $\hat{\sigma}_i ^\alpha$ ($\alpha =x, y, z$) son las matrices de Pauli para el espín $i$.  

# Solución desarrollada en Python

La solución implementada involucra el uso de las bibliotecas `numpy`, que permite, entre otras cosas, trabajar con arreglos de matrices, `matplotlib.pyplot` para poder graficar los resultados y `time` para medir tiempos de ejecución.  

Primeramente, la solución requiere de la construcción de la matriz Hamiltoniana a partir de la ecuación, que es precisamente lo que hace el método `hamiltonian` presentado en esta documentación. Además de este, como se puede comprobar, se desarrollan otros métodos que ejecutan pasos necesarios para llegar a los resultados deseados. Un método para evolucionar en el tiempo el estado cuántico inicial, otro para calcular la derivada temporal de la ecuación de Schrödinger, métodos que utilizan Runge-Kutta de cuarto orden, entre otros.  

Con todos los métodos ya definidos, se procede a su aplicación, declarando antes las matrices de Pauli y los parámetros constantes a usar. Se implementan los métodos de RK4 y de diagonalización evaluando los tiempos de ejecución respectivos; a su vez, se grafican los resultados obtenidos.  

