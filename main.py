#Python project regarding the resolution of the discrete vortex metod aplied in a NACA airfol
import Parameters as par
import Vortex_Iteration as vi
import matplotlib.pyplot as plt
import numpy as np

N = par.M #Nombre de punts

coord = np.zeros((N+1,2)) #files columnes; x y




#vi.Calc_coord_Cosinus(coord,par.p,N)
#print(coord)
cont = 0
start =0
finish = 1
step = 1
lenght = np.abs(start/step)+np.abs(finish/step) + 1

par.Parameters_definition()

angle = np.zeros((int(lenght),1))
Cl = np.zeros((int(lenght),1))#Cl con flap del ala
Cl_flap = np.zeros((int(lenght),1))#Cl del flap
Cmle = np.zeros((int(lenght),1)) #coef de momentos del perfil
Cmxh = np.zeros((int(lenght),1))

vi.Calc_coord_UNIFORM(coord, par.p, par.M, par.xh , par.eta)

infoMatrix = vi.Calc_panel(coord, N)  # VECTOR NORMAL, VECTOR TANGENTE, X LUMPED VORTEX, X CONTROL POINT
coefMatrix, RHSmatrix = vi.Iteration_Process(infoMatrix, N)
Circulation = vi.Circuilation_Calc(coefMatrix, RHSmatrix)
Cl[cont], Cl_flap[cont] = vi.Lift_Coeficient(Circulation, infoMatrix)  # Cl de el perfil completo i flap
Cmle[cont], Cmxh[cont] = vi.MomentLE_Coeficient(Circulation, infoMatrix)
print(angle[cont], Cl[cont], Cl_flap[cont], Cmle[cont], Cmxh[cont], par.xh)

print(coord)
print(*zip(*coord))
plt.scatter(*zip(*coord))
plt.show()
"""
plt.plot(angle,Cl, color='black', linestyle='dashed', linewidth = 1,
         marker='o', markerfacecolor='black', markersize=4,label = 'cmle')
plt.plot(angle,Cmxh, color='blue', linestyle='dashed', linewidth = 1,
         marker='o', markerfacecolor='black', markersize=4,label = 'cmxh')
plt.grid(color='black', linestyle='--', linewidth=0.5)

plt.legend()
plt.title('Cl vs atack angle')
plt.xlabel('Atack angle (deg)')
plt.ylabel('Cl')
plt.show()""""""
#CODI PER IMPRIMIR PER PANTALLA EL GRAFIC DE CMXH I CL
fig, ax1 = plt.subplots()

color = 'tab:red'
ax1.set_xlabel('flap deflection angle (s)')
ax1.set_ylabel('Cl', color=color)
ax1.plot(angle, Cl, color=color)
ax1.tick_params(axis='y', labelcolor=color)

ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis

color = 'tab:blue'
ax2.set_ylabel('Cmxh', color=color)  # we already handled the x-label with ax1
ax2.plot(angle, Cmxh, color=color)
ax2.tick_params(axis='y', labelcolor=color)

fig.tight_layout()  # otherwise the right y-label is slightly clipped
plt.show()
"""
"""print(f"infomatrix: ",infoMatrix[0])
print(f"A: ", coefMatrix)
print(f"Circulacion: ",Circulation)
"""
#funcionaa
