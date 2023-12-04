import Parameters as par
import Vortex_Iteration as vi
import matplotlib.pyplot as plt
import numpy as np

N = par.M  # Nombre de puntos
coord = np.zeros((N+1, 2))  # filas, columnas; x y



# Lista de valores de x_h para iterar
x_h_values = [0.7, 0.75, 0.8, 0.85]
print('|Angle [deg]|', '|Cl perfil|', '|Delta Cl flap|', '|Cmle|', '|Cmxh|', '|xh|')
for x_h in x_h_values: #recorre los diferentes valores de xh
    alfa = 0.0014682949688161284
    cont = 0
    eta_inicial = 0
    eta_final = 20
    step = 1
    lenght = np.abs(eta_inicial/step) + np.abs(eta_final/step) + 1
    angle = np.zeros((int(lenght), 1))
    hinge_pos = x_h
    Cl = np.zeros((int(lenght), 1))  # Cl con flap del ala
    Cl_flap = np.zeros((int(lenght), 1))  # Cl del flap
    Cmle = np.zeros((int(lenght), 1))  # coef de momentos del perfil respecto borde de ataque
    Cmxh = np.zeros((int(lenght), 1))  # coef de momentos del flap respecto eje de charnela
    for i in range(eta_inicial, eta_final + step, step):
        par.Parameters_definition()  # definimos los parámetros
        par.eta = i * (np.pi / 180)  # convertimos el ángulo a radianes
        vi.Calc_coord_Cosinus(coord, par.p, N, x_h, i)
        # calculamos con el ángulo de 0
        infoMatrix = vi.Calc_panel(coord, N)  # VECTOR NORMAL, VECTOR TANGENTE, X LUMPED VORTEX, X CONTROL POINT
        coefMatrix, RHSmatrix = vi.Iteration_Process(infoMatrix, N)
        Circulation = vi.Circuilation_Calc(coefMatrix, RHSmatrix)
        Cl[cont], Cl_flap[cont] = vi.Lift_Coeficient(Circulation, infoMatrix)  # Cl de el perfil completo i flap
        Cmle[cont], Cmxh[cont] = vi.MomentLE_Coeficient(Circulation, infoMatrix)
        angle[cont] = i
        print(angle[cont], Cl[cont], Cl_flap[cont], Cmle[cont], Cmxh[cont], hinge_pos)
        cont += 1

    # Aquí puedes hacer lo que necesites con los resultados específicos de x_h
    # Por ejemplo, podrías graficar los resultados para cada valor de x_h
    plt.plot(angle, Cmle, label=f'x_h = {x_h}')

# Configuración del gráfico
plt.title('Coeficiente de Sustentación (Cl) vs. Ángulo para diferentes valores de x_h')
plt.xlabel('Ángulo (grados)')
plt.ylabel('Cl')
plt.legend()
plt.grid(True)
plt.show()

