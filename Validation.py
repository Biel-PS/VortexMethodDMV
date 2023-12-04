import Parameters as par
import Vortex_Iteration as vi
import matplotlib.pyplot as plt
import numpy as np

def MomentAC_Coeficient (circulation,infopanel):
    cmac = 0
    cmxh = 0
    for i in range (0,len(circulation)):
        cmac += circulation[i]*(infopanel[i][3][0]-0.25)
        if infopanel[i][2][0]>=par.xh:
            cmxh+= circulation[i]*(infopanel[i][3][0]-par.xh)
    return -2*cmac*np.cos(par.alfa) , -2*cmxh*np.cos(par.alfa)

N = par.M #Nombre de punts
coord = np.zeros((N+1,2)) #files columnes; x y

par.Parameters_definition()


#vi.Calc_coord_Cosinus(coord,par.p,N)
#print(coord)
cont = 0

start =-5
finish = 10
step = 0.5
lenght = np.abs(start/step)+np.abs(finish/step) + 1


angle = np.zeros((int(lenght),1))
Cl = np.zeros((int(lenght),1))#Cl con flap del ala
Cl_flap = np.zeros((int(lenght),1))#Cl del flap
Cmle = np.zeros((int(lenght),1)) #coef de momentos del perfil
Cmxh = np.zeros((int(lenght),1))
Cmac = np.zeros((int(lenght),1)) #coef cm0

print('|Angle [deg]|','|Cl perfil|', '|Delta Cl flap|','|Cmle|','|Cmxh|')
for i in np.arange(start,finish+step,step):
    par.alfa = i*(np.pi/180) #CANVIAR par.eta PER par.alfa SI ES VOL FER ANALÍSI D'ANGLE D'ATAC!!
    vi.Calc_coord_Cosinus(coord, par.p, N, par.xh, par.eta)
    #calulo con el angulo de 0
    infoMatrix = vi.Calc_panel(coord,N)#VECTOR NORMAL, VECTOR TANGENTE, X LUMPED VORTEX, X CONTROL POINT
    coefMatrix,RHSmatrix = vi.Iteration_Process(infoMatrix,N)
    Circulation = vi.Circuilation_Calc(coefMatrix,RHSmatrix)


    Cl[cont],Cl_flap[cont] =(vi.Lift_Coeficient(Circulation,infoMatrix)) #Cl de el perfil completo i flap


    Cmle[cont],Cmxh[cont] = (vi.MomentLE_Coeficient(Circulation,infoMatrix))
    Cmac[cont], Cmxh[cont] = (MomentAC_Coeficient(Circulation, infoMatrix))

    angle[cont] = i
    print(angle[cont],Cl[cont],Cl_flap[cont],Cmle[cont],Cmxh[cont])
    cont += 1


angle_array = np.array(angle).flatten()
cl_array = np.array(Cl).flatten()
cmac_array = np.array(Cmac).flatten()

"Cl slope"
Cl_slope = np.polyfit(angle_array, cl_array, 1)
Cmac_slope = np.polyfit(cl_array,cmac_array,1)

print("Cl slope:",Cl_slope[0])
print(Cmac_slope)
print("Cmo:",Cmac_slope[0]*cl_array[int((start+finish)/2)]+Cmac_slope[1])
"Cmac = Cmac_slope[0]*Cl + Cmac_slope[1]"

"Alpha_0"
alfa_0= -Cl_slope[1]/Cl_slope[0]
print("alfa zero:",alfa_0)




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

plt.plot(Cl,Cmac, color='black', linestyle='dashed', linewidth = 1,
         marker='o', markerfacecolor='black', markersize=4,label = 'cm0')
plt.grid(color='black', linestyle='--', linewidth=0.5)

plt.ylim(-1.25, 1.25)
plt.legend()
plt.title('Cm0 vs Cl')
plt.xlabel('Cl')
plt.ylabel('Cm0')
plt.show()


#CODI PER IMPRIMIR PER PANTALLA EL GRAFIC DE CMXH I CL

"""
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

cont = 0
cont_xh = 0

par.xh=1
start_xh =1
finish_xh = 0.5
step_xh = -0.05
lenght_xh = np.abs(start_xh/step_xh)-np.abs(finish_xh/step_xh)+1

xh = np.zeros((int(lenght_xh),1))
angle = np.zeros((int(lenght),1))
Cl = np.zeros((int(lenght),1))#Cl con flap del ala
Cl_flap = np.zeros((int(lenght),1))#Cl del flap
Cmle = np.zeros((int(lenght),1)) #coef de momentos del perfil
Cmxh = np.zeros((int(lenght),1))
Cmac = np.zeros((int(lenght),1)) #coef cm0
alfa_0 = np.zeros((int(lenght_xh),1))
alfa_flap_0 = np.zeros((int(lenght_xh),1))
increment_alfa_0 = np.zeros((int(lenght_xh),1))


print('|Angle [deg]|','|Cl perfil|', '|Delta Cl flap|','|Cmle|','|Cmxh|')
for k in np.arange(start_xh,finish_xh+step_xh,step_xh):
    print("Xh", cont_xh, ": -------------------------")
    for i in np.arange(start,finish+step,step):
        par.alfa = i*(np.pi/180) #CANVIAR par.eta PER par.alfa SI ES VOL FER ANALÍSI D'ANGLE D'ATAC!!
        vi.Calc_coord_Cosinus(coord, par.p, N, par.xh, par.eta)
        #calulo con el angulo de 0
        infoMatrix = vi.Calc_panel(coord,N)#VECTOR NORMAL, VECTOR TANGENTE, X LUMPED VORTEX, X CONTROL POINT
        coefMatrix,RHSmatrix = vi.Iteration_Process(infoMatrix,N)
        Circulation = vi.Circuilation_Calc(coefMatrix,RHSmatrix)


        Cl[cont],Cl_flap[cont] =(vi.Lift_Coeficient(Circulation,infoMatrix)) #Cl de el perfil completo i flap


        Cmle[cont],Cmxh[cont]= (vi.MomentLE_Coeficient(Circulation,infoMatrix))
        Cmac[cont], Cmxh[cont] = (MomentAC_Coeficient(Circulation, infoMatrix))

        angle[cont] = i
        print(angle[cont],Cl[cont],Cl_flap[cont],Cmle[cont],Cmxh[cont])
        cont += 1

    angle_array = np.array(angle).flatten()
    cl_array = np.array(Cl).flatten()
    cl_flap_array = np.array(Cl_flap).flatten()

    "Cl slope"
    Cl_slope = np.polyfit(angle_array, cl_array, 1)
    Cl_flap_slope = np.polyfit(angle_array, cl_flap_array, 1)

    print("Cl slope:", Cl_slope[0])

    "Alpha_0"
    alfa_0[cont_xh] = -Cl_slope[1] / Cl_slope[0]
    "alfa_flap_0[cont_xh] = -Cl_flap_slope[1] / Cl_flap_slope[0]"
    print("alfa zero:", alfa_0)

    increment_alfa_0[cont_xh]=abs(alfa_0[cont_xh]-alfa_0[0])

    cont = 0
    angle = np.zeros((int(lenght),1))
    par.xh += step_xh
    xh[cont_xh] = 1-k
    cont_xh += 1

print(alfa_0)
print(xh)
print(increment_alfa_0)

plt.plot(xh,increment_alfa_0, color='black', linestyle='dashed', linewidth = 1,
         marker='o', markerfacecolor='black', markersize=4,label = 'increment alpha 0')
plt.grid(color='black', linestyle='--', linewidth=0.5)

plt.legend()
plt.title('Delta Alpha_0 - flap position')
plt.xlabel('Flap position')
plt.ylabel('Delta alpha_0')
plt.show()
