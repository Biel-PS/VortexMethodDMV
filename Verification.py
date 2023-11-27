import numpy as np
import matplotlib.pyplot as plt
import time
from scipy import integrate
import Parameters as par
import Vortex_Iteration as vi

par.Parameters_definition()

th_p = np.arccos(1 - 2 * par.p)

#Definimos las funciones que se integraran mas adelante
def dz(th):
    if th < th_p:
        z = (par.f / par.p**2) * (2 * par.p - 1 + np.cos(th))
    elif th_p <= th:
        z = (par.f / (1 - par.p)**2) * (2 * par.p - 1 + np.cos(th))
    return z

def int1(th):
    return dz(th) * np.cos(th)

def int2(th):
    return dz(th) * np.cos(2*th)

def int3(th):
    return dz(th) * (np.cos(th) - 1)

# Calcular la integral definida de 0 a pi de las diferentes funcionesCULO
result0, error0 = integrate.quad(dz, 0, np.pi)
result1, error1 = integrate.quad(int1, 0, np.pi)
result2, error2 = integrate.quad(int2, 0, np.pi)
result3, error3 = integrate.quad(int3, 0, np.pi)

# Encontrar coef A0,1,2
A0 = par.alfa - (1 / np.pi) * result0
A1 = (2 / np.pi) * result1
A2= (2 / np.pi) * result2

print(f"A0: ",A0)
print(f"A1: ",A1)
print(f"A2: ",A2)

#Encontrar alfa_lo_tat, CL_tat, CM0_tat (sin flap)
alfa_lo_tat=-1/np.pi*result3
CL_tat=(2*A0+A1)*np.pi
CM0_tat=(A2-A1)*np.pi/4

#Con flap
if par.eta>0:
    th_h=np.arccos(1-2*par.xh)
    alfa_lo_tat=alfa_lo_tat-(par.eta/np.pi)*(np.pi-th_h+np.sin(th_h))
    CL_tat=CL_tat+2*(np.pi-th_h+np.sin(th_h))*par.eta
    CM0_tat=CM0_tat-np.sin(th_h)*(1-np.cos(th_h))*par.eta/2

CM_tat=CM0_tat-CL_tat/4

print(CL_tat)
print(CM_tat)


####################

#Ara trobem valors de Cl i CM0 mitjançant DVM per N panells (entre 1 i 400)



par.Parameters_definition()

cont = 0

start =1
finish = 400
step = 1
lenght = np.abs(start/step)+np.abs(finish/step) + 1


N = np.zeros((int(lenght),1))

Cl_flap = np.zeros((int(lenght),1))#Cl del flap

Cmxh = np.zeros((int(lenght),1))
print('|Panels number|','|Cl perfil|','|Cmle|')


data = np.zeros((3, 400))

for i in range(start,finish+step,step):
    start_time = time.time()
    coord = np.zeros((i + 1, 2))  # files columnes; x y
    par.alfa = 4*(np.pi/180) #CANVIAR par.eta PER par.alfa SI ES VOL FER ANALÍSI D'ANGLE D'ATAC!!
    vi.Calc_coord_Cosinus(coord, par.p, i, par.xh, par.eta)
    #calulo con el angulo de 0
    infoMatrix = vi.Calc_panel(coord,i)#VECTOR NORMAL, VECTOR TANGENTE, X LUMPED VORTEX, X CONTROL POINT
    coefMatrix,RHSmatrix = vi.Iteration_Process(infoMatrix,i)
    Circulation = vi.Circuilation_Calc(coefMatrix,RHSmatrix)

    """print(infoMatrix)
    print(coefMatrix)
    print(RHSmatrix)
    print(Circulation)
"""

    data[0, i - 1] = time.time() - start_time


    data[1, i-1], Cl_flap[cont] =(vi.Lift_Coeficient(Circulation,infoMatrix)) #Cl de el perfil completo i flap


    data[2, i - 1], Cmxh[cont] = (vi.MomentLE_Coeficient(Circulation,infoMatrix))

    N[cont] = i
    print(N[cont],data[1, i-1],data[2, i-1])
    cont += 1

error = np.zeros((2, 400))
error[0, :] = np.abs((data[1, i-1] - CL_tat) / CL_tat) * 100
error[1, :] = np.abs((data[2, i-1] - CM_tat) / CM_tat) * 100


# Plots: elapsed time, CL, CM
plt.figure(figsize=(10, 6))

# Elapsed time plot
plt.subplot(3, 1, 1)
plt.plot(data[0, :], label='Elapsed Time')
plt.xlabel('Panel Number')
plt.ylabel('Elapsed Time (s)')
plt.legend()

# CL convergence plot
plt.subplot(3, 1, 2)
plt.plot(data[1, :], label='CL')
plt.xlabel('Panel Number')
plt.ylabel('CL')
plt.legend()

# CM convergence plot
plt.subplot(3, 1, 3)
plt.plot(data[2, :], label='CM')
plt.xlabel('Panel Number')
plt.ylabel('CM')
plt.legend()

plt.tight_layout()
plt.show()




"""
def DVM(i):
    # Replace this with your DVM implementation
    CL = 0.7
    CM = 0.3
    return CL, CM


M = 400  # Maximum number of panels
data = np.zeros((3, M))


for i in range(1, M + 1):
    start_time = time.time()
    CL = vi.Lift_Coeficient(circulation, infopanel)[0]
    CM = vi.MomentLE_Coeficient(circulation, infopanel)[0]
    data[0, i - 1] = CL
    data[1, i - 1] = CM
    data[2, i - 1] = time.time() - start_time

error = np.zeros((2, M))
error[0, :] = np.abs((data[0, :] - CL_tat) / CL_tat) * 100
error[1, :] = np.abs((data[1, :] - CM_tat) / CM_tat) * 100

# Plots: elapsed time, CL, CM
plt.figure(figsize=(10, 6))

# Elapsed time plot
plt.subplot(3, 1, 1)
plt.plot(data[2, :], label='Elapsed Time')
plt.xlabel('Panel Number')
plt.ylabel('Elapsed Time (s)')
plt.legend()

# CL convergence plot
plt.subplot(3, 1, 2)
plt.plot(data[0, :], label='CL')
plt.xlabel('Panel Number')
plt.ylabel('CL')
plt.legend()

# CM convergence plot
plt.subplot(3, 1, 3)
plt.plot(data[1, :], label='CM')
plt.xlabel('Panel Number')
plt.ylabel('CM')
plt.legend()

plt.tight_layout()
plt.show()
"""