# Modelo de Ising Cuántico Unidimensional en una Grilla de N Espines: Dinámica de Muchos Cuerpos

**Coordinador**: Marlon Brenes

####Estudiantes####: 
 - Daniel Lacayo Zúñiga
 - Jorge Salas Rodríguez
 - Jorge Ramsés Jiménez Ramírez
 - Bryton Ramírez Calderón
 - Josué David Gómez Castillo

La dinámica de un estado puro para un sistema cuántico aislado se rige bajo la ecuación Schrödinger (\( \hbar = 1 \)):

$$
\frac{\partial |\psi(t)\rangle}{\partial t} = -i\hat{H} |\psi(t)\rangle,
$$

cuya solución formal está dada por

$$
|\psi(t)\rangle = e^{-i\hat{H}(t-t_0)} |\psi(t = t_0)\rangle.
$$

Es decir, la solución involucra resolver de manera numérica la ecuación diferencial Eq. (3) o evaluar de alguna forma la exponencial de la matriz Eq. (4). La idea del proyecto es evaluar la dinámica del modelo de Ising Eq. (2) empezando de algún estado inicial.

## Milestones:

- Con base en argumentos numéricos, evaluar cuál de las dos metodologías es mejor implementar para la solución numérica
- Construir la matriz Hamiltoniana mediante productos tensoriales
- Resolver el sistema utilizando algún estado inicial y visualizar su dinámica
- Implementar la solución en **Python**
- Implementar la solución en **C++**
- Encontrar una forma de paralelizar el algoritmo y evaluar la aceleración

