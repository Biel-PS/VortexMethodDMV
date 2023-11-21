#Python project regarding the resolution of the discrete vortex metod aplied in a NACA airfol
import Parameters as par
import Vortex_Iteration as vi
import matplotlib.pyplot as plt
import numpy as np

N = par.M #Nombre de punts
coord = np.zeros((N+1,2)) #files columnes; x y
par.Parameters_definition()
vi.Calc_coord_Cosinus(coord,par.p,N,par.xh)
#vi.Calc_coord_Cosinus(coord,par.p,N)
#print(coord)
cont = 0

start =-10
finish = 10
step = 1
lenght = np.abs(start/step)+np.abs(finish/step) + 1
angle = np.zeros((int(lenght),1))
Cl = np.zeros((int(lenght),1))

for i in range(start,finish+step,step):
    par.alfa = i*(np.pi/180)
    infoMatrix = vi.Calc_panel(coord,N)#VECTOR NORMAL, VECTOR TANGENTE,X LUMPED VORTEX, X CONTROL POINT

    coefMatrix,RHSmatrix = vi.Iteration_Process(infoMatrix,N)

    Circulation = vi.Circuilation_Calc(coefMatrix,RHSmatrix)

    """print(infoMatrix)
    print(coefMatrix)
    print(RHSmatrix)
    print(Circulation)
"""
    Cl[cont]=(vi.Lift_Coeficient(Circulation))
    angle[cont] = i
    #print(angle[cont],Cl[cont])
    cont += 1

"""print(coord)
print(*zip(*coord))
plt.scatter(*zip(*coord))
plt.show()"""

plt.plot(angle,Cl, color='black', linestyle='dashed', linewidth = 1,
         marker='o', markerfacecolor='black', markersize=4)

plt.title('Cl vs atack angle')
plt.xlabel('Atack angle (deg)')
plt.ylabel('Cl')
plt.show()

"""print(f"infomatrix: ",infoMatrix[0])
print(f"A: ", coefMatrix)
print(f"Circulacion: ",Circulation)
"""
#proba
