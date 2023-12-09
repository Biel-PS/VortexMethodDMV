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
step = 1
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


"Cl vs alpha"

plt.plot(angle,Cl, color='black', linestyle='dashed', linewidth = 1,
         marker='o', markerfacecolor='black', markersize=4,label = 'Cl')

plt.grid(color='black', linestyle='--', linewidth=0.5)

plt.legend()
plt.title('Linear range of lift coefficient')
plt.xlabel('Angle of attack (deg)')
plt.ylabel('Cl')
plt.show()

"Cmle vs alpha"
plt.plot(angle,Cmle, color='blue', linestyle='dashed', linewidth = 1,
         marker='o', markerfacecolor='black', markersize=4,label = 'Cmle')
plt.grid(color='black', linestyle='--', linewidth=0.5)

plt.legend()
plt.title('Linear range of moment coefficient at leading edge')
plt.xlabel('Angle of attack (deg)')
plt.ylabel('Cmle')
plt.show()

"Cm0 vs Cl"
plt.plot(Cl,Cmac, color='black', linestyle='dashed', linewidth = 1,
         marker='o', markerfacecolor='black', markersize=4,label = 'Cm0')
plt.grid(color='black', linestyle='--', linewidth=0.5)

plt.ylim(-1.25, 1.25)
plt.legend()
plt.title('Moment coefficient at aerodynamic center')
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
cont_eta = 0


start_xh = 0.85
finish_xh = 0.70
step_xh = -0.05
lenght_xh = np.abs(start_xh/step_xh)-np.abs(finish_xh/step_xh)+1
par.xh = start_xh

par.eta=0
start_eta =0
finish_eta = 10
step_eta = 1
lenght_eta = np.abs(start_eta/step_eta)+np.abs(finish_eta/step_eta)+1
eta = np.zeros((int(lenght_eta),1))

xh = np.zeros((int(lenght_xh),1))
angle = np.zeros((int(lenght),1))
Cl = np.zeros((int(lenght),1))#Cl con flap del ala
Cl_flap = np.zeros((int(lenght),1))#Cl del flap
Cmle = np.zeros((int(lenght),1)) #coef de momentos del perfil
Cmxh = np.zeros((int(lenght),1))
Cmac = np.zeros((int(lenght),1)) #coef cm0
alfa_0 = np.zeros((int(lenght_xh),int(lenght_eta)))
alfa_flap_0 = np.zeros((int(lenght_xh),1))
increment_alfa_0 = np.zeros((int(lenght_xh),1))

effectiveness = np.zeros((int(lenght_xh),1))



print('|Angle [deg]|','|Cl perfil|', '|Delta Cl flap|','|Cmle|','|Cmxh|')
for k in np.arange(start_xh,finish_xh,step_xh):
    print("Xh", cont_xh, ": -------------------------")
    for j in np.arange(start_eta,finish_eta+step_eta,step_eta):
        par.eta = j * (np.pi/180)
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
        alfa_0[cont_xh][cont_eta] = -Cl_slope[1] / Cl_slope[0]
        "alfa_flap_0[cont_xh] = -Cl_flap_slope[1] / Cl_flap_slope[0]"
        print("alfa zero:", alfa_0)

        "increment_alfa_0[cont_xh]=abs(alfa_0[cont_xh]-alfa_0[0])"

        cont = 0
        angle = np.zeros((int(lenght),1))

        eta[cont_eta]=par.eta * (180/np.pi)
        cont_eta += 1

    cont_eta=0
    par.eta=0
    alfa_0_array = np.array(alfa_0[cont_xh]).flatten()
    eta_array = np.array(eta).flatten()

    effectiveness_slope = np.polyfit(eta_array, alfa_0_array, 1)
    effectiveness[cont_xh] = abs(effectiveness_slope[0])

    par.xh += step_xh
    xh[cont_xh] = 1-k
    cont_xh += 1

"effectiveness[0]=0"
print(alfa_0)
print(xh)
print(effectiveness)

effectiveness_Experimental = [0.39,0.475,0.54,0.595]

plt.plot(xh,effectiveness, color='black', linestyle='dashed', linewidth = 1,
         marker='o', markerfacecolor='black', markersize=4,label = 'DVM')
plt.plot(xh,effectiveness_Experimental, color='blue', linestyle='dashed', linewidth = 1,
         marker='o', markerfacecolor='black', markersize=4,label = 'Experimental')
plt.plot(xh,effectiveness*0.85, color='red', linestyle='dashed', linewidth = 1,
         marker='o', markerfacecolor='black', markersize=4,label = 'DVM Corrected')
plt.grid(color='black', linestyle='--', linewidth=0.5)

plt.legend()
plt.title('Flap effectiveness - Flap position')
plt.xlabel('Flap position')
plt.ylabel('Flap effectiveness')
plt.show()
