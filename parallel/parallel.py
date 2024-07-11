# Código de Simulación Cuántica

import numpy as np
import matplotlib.pyplot as plt
import time

# Definición de las matrices de Pauli y otras constantes
sx = np.array([[0, 1], [1, 0]], dtype=complex)
sz = np.array([[1, 0], [0, -1]], dtype=complex)
iden = np.eye(2, dtype=complex)

def operacion(orden):
    """Realiza el producto tensorial de una lista de operadores.

    Args:
        orden (list): Lista de matrices (operadores).

    Returns:
        numpy.ndarray: Producto tensorial de los operadores en la lista.
    """
    producto = orden[0]
    for i in orden[1:]:
        producto = np.kron(producto, i)
    return producto

def hamiltonian(J, g, N):
    """Construye el Hamiltoniano del sistema.

    Args:
        J (float): Constante de acoplamiento.
        g (float): Fuerza del campo externo.
        N (int): Número de espines.

    Returns:
        numpy.ndarray: Hamiltoniano del sistema.
    """
    spins = np.zeros((2**N, 2**N), dtype=complex)
    field = np.zeros((2**N, 2**N), dtype=complex)
    for i in range(N):
        j = (i + 1) % N
        orden = [iden] * N
        orden[i] = sz
        orden[j] = sz
        elemento = operacion(orden)
        spins += elemento
    for i in range(N):
        orden = [iden] * N
        orden[i] = sx
        elemento = operacion(orden)
        field += elemento
    H = -J * spins - g * field
    return H

def evolucion_cuántica(psi_0_transformed, t, eigenvalues, eigenvectors):
    """Evoluciona el estado cuántico en el tiempo usando la diagonalización del Hamiltoniano.

    Args:
        psi_0_transformed (numpy.ndarray): Estado inicial transformado.
        t (float): Tiempo de evolución.
        eigenvalues (numpy.ndarray): Valores propios del Hamiltoniano.
        eigenvectors (numpy.ndarray): Vectores propios del Hamiltoniano.

    Returns:
        numpy.ndarray: Estado cuántico evolucionado.
    """
    return eigenvectors @ (np.exp(-1j * eigenvalues * t) * psi_0_transformed)

def schrodinger_derivative(H, psi):
    """Calcula la derivada del estado cuántico según la ecuación de Schrödinger.

    Args:
        H (numpy.ndarray): Hamiltoniano del sistema.
        psi (numpy.ndarray): Estado cuántico.

    Returns:
        numpy.ndarray: Derivada del estado cuántico.
    """
    return -1j * H @ psi

def rk4_step(H, psi, t, dt):
    """Realiza un paso del método de Runge-Kutta de cuarto orden (RK4).

    Args:
        H (numpy.ndarray): Hamiltoniano del sistema.
        psi (numpy.ndarray): Estado cuántico.
        t (float): Tiempo inicial.
        dt (float): Paso de tiempo.

    Returns:
        numpy.ndarray: Estado cuántico actualizado.
    """
    k1 = dt * schrodinger_derivative(H, psi)
    k2 = dt * schrodinger_derivative(H, psi + 0.5 * k1)
    k3 = dt * schrodinger_derivative(H, psi + 0.5 * k2)
    k4 = dt * schrodinger_derivative(H, psi + k3)
    return psi + (k1 + 2*k2 + 2*k3 + k4) / 6

def expectation_value(state, observable):
    """Calcula el valor de expectación de un observable en un estado cuántico dado.

    Args:
        state (numpy.ndarray): Estado cuántico.
        observable (numpy.ndarray): Observable cuántico.

    Returns:
        float: Valor de expectación.
    """
    return np.vdot(state, observable @ state).real

def operador_z_n_espin(i, N):
    """Crea el operador Z para el i-ésimo espín en un sistema de N espines.

    Args:
        i (int): Índice del espín.
        N (int): Número de espines.

    Returns:
        numpy.ndarray: Operador Z para el i-ésimo espín.
    """
    operadores = [iden] * N
    operadores[i] = sz
    return operacion(operadores)

# Parámetros del sistema
N = 12
J = 1.0
g = 0.5
psi_0 = np.zeros(2**N, dtype=complex)
psi_0[0] = 1
espines = list(range(N))

# Construcción del Hamiltoniano
H = hamiltonian(J, g, N)

# Diagonalización del Hamiltoniano y cálculo de psi_0'
start_time_diag = time.time()
eigenvalues, eigenvectors = np.linalg.eigh(H)
psi_0_transformed = eigenvectors.conj().T @ psi_0
diag_total_time = time.time() - start_time_diag

# Tiempo de evolución y paso de tiempo
times = np.linspace(0, 10, 100)
dt = times[1] - times[0]

# Evolución diagonalizada
start_time_diag = time.time()
for i in espines:
    observable = operador_z_n_espin(i, N)
    expectation_values = np.zeros(len(times))
    for idx, t in enumerate(times):
        psi_t = evolucion_cuántica(psi_0_transformed, t, eigenvalues, eigenvectors)
        expectation_values[idx] = expectation_value(psi_t, observable)
    plt.plot(times, expectation_values, label=f'Diagonalización, espín {i}')
diag_total_time = time.time() - start_time_diag

# Evolución con RK4
start_time_rk4 = time.time()
for i in espines:
    observable = operador_z_n_espin(i, N)
    expectation_values = np.zeros(len(times))
    psi = psi_0.copy()
    for idx, t in enumerate(times):
        expectation_values[idx] = expectation_value(psi, observable)
        psi = rk4_step(H, psi, t, dt)
    plt.plot(times, expectation_values, '--', label=f'RK4, espín {i}')
rk4_total_time = time.time() - start_time_rk4

# Imprimir tiempos
print(f"Tiempo de diagonalización y evolución: {diag_total_time:.4f} s")
print(f"Tiempo de evolución con RK4: {rk4_total_time:.4f} s")

# Graficar resultados de evolución temporal
plt.xlabel('Tiempo')
plt.ylabel('Valor de expectación')
plt.title('Comparación de la evolución temporal del valor de expectación')
plt.legend()
plt.show()
