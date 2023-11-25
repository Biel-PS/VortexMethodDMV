import numpy as np
from scipy import integrate
import Parameters as par
import Vortex_Iteration as vi

th_p = np.arccos(1 - 2 * par.p)

#Definimos las funciones que se integraran mas adelante
def dz(th):
    if th < th_p:
        z = (par.f / par.p**2) * (2 * par.p - 1 + np.cos(th))
    elif th_p <= th:
        z = (par.f / (1 - par.p**2)) * (2 * par.p - 1 + np.cos(th))
    return z

def int1(th):
    return dz(th) * np.cos(th)

def int2(th):
    return dz(th) * np.cos(2*th)

def int3(th):
    return dz(th) * (np.cos(th) - 1)

# Calcular la integral definida de 0 a pi de las diferentes funciones
result0, error0 = integrate.quad(dz, 0, np.pi)
result1, error1 = integrate.quad(int1, 0, np.pi)
result2, error2 = integrate.quad(int2, 0, np.pi)
result3, error3 = integrate.quad(int3, 0, np.pi)

# Encontrar coef A0,1,2
A0 = par.alfa * np.pi / 180 - (1 / np.pi) * result0
A1 = 2 / np.pi * result1
A2=2 / np.pi * result2

print(A0)
print(A1)
print(A2)

#Encontrar alfa_lo_tat, CL_tat, CM0_tat (sin flap)
alfa_lo_tat=-1/np.pi*result3
CL_tat=(2*A0+A1)*np.pi
CM0_tat=(A2-A1)*np.pi/4

#Con flap
if par.eta>0:
    th_h=np.arccos(1-2*par.xh)
    alfa_lo_tat=alfa_lo_tat-(par.eta/np.pi)*(np.pi-th_h+np.sin(th_h))
    CL_tat=CL_tat+2*(np.pi-th_h+np.sin(th_h))*par.eta
    CM0_tat=CM0_tat-np.sin(th_h)*(1-cos(th_h))*par.eta/2

CM_tat=CM0_tat-CL_tat/4

