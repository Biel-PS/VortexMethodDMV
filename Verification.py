import numpy as np
import matplotlib.pyplot as plt
import Parameters as par
import Vortex_Iteration as vi

def main():
    # Definir parámetros
    par.Parameters_definition()

    # Inicializar variables
    start = 1
    finish = 400

    # Almacenar resultados
    N_values = np.arange(start, finish + 1)
    Cl_values = []

    # Calcular Cl para cada valor de N
    for N in N_values:
        # Crear nodos y calcular resultados
        coord = np.zeros((N + 1, 2))
        vi.Calc_coord_Cosinus(coord, par.p, N, par.xh)

        # Inicializar variables para almacenar resultados
        Cl = 0

        # Calcular resultados para un ángulo de ataque de 4 grados
        par.alfa = 4 * (np.pi / 180)
        infoMatrix = vi.Calc_panel(coord, N)
        coefMatrix, RHSmatrix = vi.Iteration_Process(infoMatrix, N)
        Circulation = vi.Circuilation_Calc(coefMatrix, RHSmatrix)
        Cl = vi.Lift_Coeficient(Circulation)*2

        # Almacenar resultados
        Cl_values.append(Cl)


    # Graficar Cl respecto a N
    plt.plot(N_values, Cl_values, color='black', linestyle='dashed', linewidth=1,
             marker='o', markerfacecolor='black', markersize=4)

    plt.title('Coeficiente de Sustentación (Cl) vs Número de Paneles')
    plt.xlabel('Número de Paneles (N)')
    plt.ylabel('Cl')
    plt.grid(True)
    plt.show()
