#Python project regarding the resolution of the discrete vortex metod aplied in a NACA airfol
import Parameters as par
import Vortex_Iteration as vi
import matplotlib.pyplot as plt
import numpy as np

N = par.M #Nombre de punts

coord = np.zeros((N+1,2)) #files columnes; x y
coord_flap = np.zeros((N+1,2)) #files columnes; x y

par.Parameters_definition()

vi.Calc_coord_Cosinus(coord,par.p,N,par.xh,0) #sense flap
vi.Calc_coord_Cosinus(coord_flap,par.p,N,par.xh,par.eta) #sense flap


#vi.Calc_coord_Cosinus(coord,par.p,N)
#print(coord)
cont = 0

start =-10
finish = 10
step = 1
lenght = np.abs(start/step)+np.abs(finish/step) + 1


angle = np.zeros((int(lenght),1))
Cl = np.zeros((int(lenght),1))#Cl con flap
Cl_flap = np.zeros((int(lenght),1))#Cl sin flap
Cmle = np.zeros((int(lenght),1)) #coef de momentos del perfil
Cmxh = np.zeros((int(lenght),1))
print('|Angle [deg]|','|Cl perfil|', '|Delta Cl flap|','|Cmle|','|Cmxh|')
for i in range(start,finish+step,step):
    par.alfa = i*(np.pi/180)

    #calulo con el angulo de 0
    infoMatrix = vi.Calc_panel(coord,N)#VECTOR NORMAL, VECTOR TANGENTE,X LUMPED VORTEX, X CONTROL POINT
    coefMatrix,RHSmatrix = vi.Iteration_Process(infoMatrix,N)
    Circulation = vi.Circuilation_Calc(coefMatrix,RHSmatrix)
    #parametros calculados para un angulo de falp establecido en par.
    infoMatrix_flap =  vi.Calc_panel(coord_flap,N)#VECTOR NORMAL, VECTOR TANGENTE,X LUMPED VORTEX, X CONTROL POINT
    coefMatrix_flap, RHSmatrix_flap = vi.Iteration_Process(infoMatrix_flap, N)
    Circulation_flap = vi.Circuilation_Calc(coefMatrix_flap, RHSmatrix_flap)

    """print(infoMatrix)
    print(coefMatrix)
    print(RHSmatrix)
    print(Circulation)
"""
    Cl[cont]=(vi.Lift_Coeficient(Circulation_flap)) #Cl de el perfil completo
    Cl_flap[cont] = (vi.Lift_Coeficient(Circulation_flap - Circulation)) #delta de Cl a causa del flap

    Cmle[cont] = (vi.MomentLE_Coeficient(Circulation_flap,infoMatrix))
    Cmxh[cont] = (vi.MomentXH_Coeficient_OnlyFlap(Cl_flap[cont],Cmle[cont],vi.MomentLE_Coeficient(Circulation,infoMatrix)))
    angle[cont] = i
    print(angle[cont],Cl[cont],Cl_flap[cont],Cmle[cont],Cmxh[cont])
    cont += 1

"""print(coord)
print(*zip(*coord))
plt.scatter(*zip(*coord))
plt.show()"""

plt.plot(angle,Cmle, color='black', linestyle='dashed', linewidth = 1,
         marker='o', markerfacecolor='black', markersize=4,label = 'cmle')
plt.plot(angle,Cl, color='blue', linestyle='dashed', linewidth = 1,
         marker='o', markerfacecolor='black', markersize=4,label = 'cl')
plt.grid(color='black', linestyle='--', linewidth=0.5)

plt.legend()
plt.title('Cl vs atack angle')
plt.xlabel('Atack angle (deg)')
plt.ylabel('Cl')
plt.show()

"""print(f"infomatrix: ",infoMatrix[0])
print(f"A: ", coefMatrix)
print(f"Circulacion: ",Circulation)
"""
#funcionaa
