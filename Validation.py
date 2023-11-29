import Parameters as par
import Vortex_Iteration as vi
import matplotlib.pyplot as plt
import numpy as np


N = par.M #Nombre de punts

coord = np.zeros((N+1,2)) #files columnes; x y

par.Parameters_definition()


#vi.Calc_coord_Cosinus(coord,par.p,N)
#print(coord)
cont = 0

start =-10
finish = 10
step = 1
lenght = np.abs(start/step)+np.abs(finish/step) + 1


angle  = np.zeros((int(lenght),1))
print (angle)
Cl = np.zeros((int(lenght),1))#Cl con flap del ala
Cl_flap = np.zeros((int(lenght),1))#Cl del flap
Cmle = np.zeros((int(lenght),1)) #coef de momentos del perfil
Cmxh = np.zeros((int(lenght),1))
print('|Angle [deg]|','|Cl perfil|', '|Delta Cl flap|','|Cmle|','|Cmxh|')
for i in range(start,finish+step,step):
    par.alfa = i*(np.pi/180) #CANVIAR par.eta PER par.alfa SI ES VOL FER ANALÍSI D'ANGLE D'ATAC!!
    vi.Calc_coord_Cosinus(coord, par.p, N, par.xh, par.eta)
    #calulo con el angulo de 0
    infoMatrix = vi.Calc_panel(coord,N)#VECTOR NORMAL, VECTOR TANGENTE, X LUMPED VORTEX, X CONTROL POINT
    coefMatrix,RHSmatrix = vi.Iteration_Process(infoMatrix,N)
    Circulation = vi.Circuilation_Calc(coefMatrix,RHSmatrix)


    Cl[cont],Cl_flap[cont] =(vi.Lift_Coeficient(Circulation,infoMatrix)) #Cl de el perfil completo i flap


    Cmle[cont],Cmxh[cont] = (vi.MomentLE_Coeficient(Circulation,infoMatrix))

    angle[cont] = i
    print(angle[cont],Cl[cont],Cl_flap[cont],Cmle[cont],Cmxh[cont])
    cont += 1

print(np.transpose(angle))
angle_array = np.array(angle).flatten()
cl_array = np.array(Cl).flatten()


Cl_slope = np.polyfit(angle_array, cl_array, 1)

print(Cl_slope)



plt.plot(angle,Cl, color='black', linestyle='dashed', linewidth = 1,
         marker='o', markerfacecolor='black', markersize=4,label = 'cl')
plt.plot(angle,Cmle, color='blue', linestyle='dashed', linewidth = 1,
         marker='o', markerfacecolor='black', markersize=4,label = 'cmle')
plt.grid(color='black', linestyle='--', linewidth=0.5)

plt.legend()
plt.title('Cl vs atack angle')
plt.xlabel('Atack angle (deg)')
plt.ylabel('Cl')
plt.show()
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


