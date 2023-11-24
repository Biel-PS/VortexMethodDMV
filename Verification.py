import numpy as np
import matplotlib.pyplot as plt
import Parameters as par
import Vortex_Iteration as vi
from scipy import integrate

th_p=np.arccos(1-2*par.p)

def dz(th):
    if th < th_p:
        z = (par.f/par.p**2)*(2*par.p-1+np.cos(th))
    elif th_p<=th<=np.pi:
        z = (par.f/(1-par.p**2))*(2*par.p-1+np.cos(th))
    return z


# Definir la función que quieres integrar
#def f(x):
    #return dz/dx

# Calcular la integral definida de 0 a 2 de la función f(x)
result0, error = integrate.quad(dz, 0, np.pi)

#vull agafar els diferents valors de x,z del vi.Calc coord cosinus del biel per poder fer la integral.


#trobar coef A0,1,2
A0=par.alfa-(1/np.pi)*result0
print(A0)

################################
"""
def g(x):
    return (dz/dx)*np.cos(z)

# Calcular la integral definida de 0 a 2 de la función f(x)
#result1, error = integrate.quad(g, 0, np.pi)

#A1=2/np.pi*result1

################################

def h(x):
    return (dz/dx)*np.cos(2*z)

# Calcular la integral definida de 0 a 2 de la función f(x)
#result3, error = integrate.quad(h, 0, np.pi)

#A2=2/np.pi*result3

################################
def j(x):
    return (dz/dx)*(np.cos(z)-1)

# Calcular la integral definida de 0 a 2 de la función f(x)
#result4, error = integrate.quad(g, 0, np.pi)

#alfa_lo_tat=-1/np.pi*result4

################################

#CL_tat=(2*A0+A1)*np.pi
#CM0_tat=(A2-A1)*np.pi/4

#Amb flap

#if par.eta>0:
    #th_h=np.arccos(1-2*par.xh)
    #alfa_lo_tat=alfa_lo_tat-(par.eta/np.pi)*(np.pi-th_h+np.sin(th_h))
    #CL_tat=CL_tat+2*(np.pi-th_h+np.sin(th_h))*par.eta
    #CM0_tat=CM0_tat-np.sin(th_h)*(1-cos(th_h))*par.eta/2

#CM_tat=CM0_tat-CL_tat/4


"""