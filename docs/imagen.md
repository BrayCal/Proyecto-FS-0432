# Gráfica obtenida a partir de Python:

La gráfica obtenida representa cómo cambia en el tiempo el valor esperado del operador $\sigma_{z}$ para cada uno de los cuatro espines en una cadena, usando el modelo de Ising con un campo transversal. Se empleó el método de Runge-Kutta de cuarto orden (RK4) para calcular la evolución temporal del sistema.

![](img/C1ramses.png)

# Gráfica obtenida a partir de Python Diagonalization:
En las siguientes gráficas, observaremos las diferencias entre el método RK4 y el método de diagonalización. Utilizaremos una grilla de 1001 puntos y el valor de N será inicialmente 10 y luego la aumentamos a 11.
Veremos las gráficas obtenidas con ambos métodos: RK4 y diagonalización, y analizaremos los tiempos de evolución correspondientes.

Grafica con el metodo RK4 (N=10), tiempo obtenido 62.48s 
![](img/10nrk4.png)

Grafica para el metodos de diagonalización (N=10), tiempo obtenido: 31.794s
![](img/10ndiag.png)

###Ahora veamos el para N = 11.

Grafica con el metodo RK4, tiempo obtenido: 337.792s 
![](img/11nrk4.png)

Grafica con el metodo de diagonalización, tiempo obtenido: 137.4689s

![](img/11ndiag.png)

## Conclusión:

Con base en los resultados, se determino que para N pequeños el método de diagonalización tiende a ser mejor que RK4, sin embargo, conforme crece el valor de N asignado RK4 se acerca en cuanto al tiempo a la diagonalización. A partir de esto se concluye que computacionalmente la diagonalización es  “mejor” método para resolver el problema.
