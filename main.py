#Python project regarding the resolution of the discrete vortex metod aplied in a NACA airfol
import Parameters as par
import Vortex_Iteration as vi
import matplotlib.pyplot as plt
import numpy as np

N = par.M + 1 #Nombre de punts
coord = np.zeros((N+1,2)) #files columnes; x y
par.Parameters_definition()
vi.Calc_coord_Cosinus(coord,par.p,N)
#vi.Calc_coord_Cosinus(coord,par.p,N)
#print(coord)
cont = 0
start =-16
finish = 16
step = 2
lenght = np.abs(start)+np.abs(finish)+np.abs(step)
angle = np.zeros((lenght,1))
Cl = np.zeros((lenght,1))

for i in range(start,finish+step,step):
    par.alfa = i*(np.pi/180)
    infoMatrix = vi.Calc_panel(coord,N)#VECTOR NORMAL, VECTOR TANGENTE,X LUMPED VORTEX, X CONTROL POINT

    coefMatrix,RHSmatrix = vi.Iteration_Process(infoMatrix,N)

    Circulation = vi.Circuilation_Calc(coefMatrix,RHSmatrix)

    Cl[cont]=(vi.Lift_Coeficient(Circulation))
    angle[cont] = i
    print(angle[cont],Cl[cont])
    cont += 1



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
#coment
#AAAAA
