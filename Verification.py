import numpy as np
import matplotlib.pyplot as plt
import time
from scipy import integrate
import Parameters as par
import Vortex_Iteration as vi

par.Parameters_definition()
par.alfa = 4*(np.pi/180)
par.eta = 0
par.xh = 1
#Límite de integración
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

# Calcular la integral definida de 0 a pi de las diferentes funciones
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

#Encontrar alfa_lo_tat, CL_tat, CM0_tat, CM_LE_tat (sin flap)
alfa_lo_tat=-1/np.pi*result3
CL_tat=(2*A0+A1)*np.pi
CM0_tat=(A2-A1)*np.pi/4
CM_LE_tat=CM0_tat-CL_tat/4

print(CL_tat)
print(CM_LE_tat)



#Ara trobem valors de Cl i CM0 mitjançant DVM per N panells (entre 1 i 400)


cont = 0

start =1
finish = 200
step = 1
lenght = 200


N = np.zeros(int(lenght))
Cl_flap = np.zeros(int(lenght))#Cl del flap
Cmxh = np.zeros(int(lenght))

print('|Panels number|','|Cl perfil|','|Cmle|', '|Error_CL|', '|Error_CM|')


data = np.zeros((3, int(lenght)))
error = np.zeros((2, int(lenght)))

for i in range(start,finish+step,step):
    start_time = time.time()
    coord = np.zeros((i + 1, 2))  # files columnes; x y
    vi.Calc_coord_Cosinus(coord, par.p, i, par.xh, par.eta)
    #calulo con el angulo de 0
    infoMatrix = vi.Calc_panel(coord,i)#VECTOR NORMAL, VECTOR TANGENTE, X LUMPED VORTEX, X CONTROL POINT
    coefMatrix,RHSmatrix = vi.Iteration_Process(infoMatrix,i)
    Circulation = vi.Circuilation_Calc(coefMatrix,RHSmatrix)


    data[0, cont] = time.time() - start_time
    data[1, cont], Cl_flap[cont] = vi.Lift_Coeficient(Circulation,infoMatrix) #Cl de el perfil completo i flap
    data[2, cont], Cmxh[cont] = vi.MomentLE_Coeficient(Circulation,infoMatrix)

    error[0, cont] = np.abs((data[1, cont] - CL_tat) / CL_tat) * 100
    error[1, cont] = np.abs((data[2, cont] - CM_LE_tat) / CM_LE_tat) * 100

    N[cont] = i




    print(N[cont],data[1, cont],data[2, cont], error[0,cont], error[1,cont])

    cont += 1


# Elapsed time plot
plt.figure(figsize=(8, 6))
plt.plot(N, data[0, :], label='Elapsed Time')
plt.xlabel('Panel Number')
plt.ylabel('Elapsed Time (s)')
plt.legend()
plt.title('Elapsed Time Plot')
plt.show()

# ...
# ...

# CL convergence plot with error
fig, ax1 = plt.subplots(figsize=(8, 6))

color = 'tab:red'
ax1.set_xlabel('Panel Number')
ax1.set_ylabel('CL', color=color)
ax1.plot(N, data[1, :], label='CL', color=color)
ax1.axhline(y=CL_tat, linestyle='--', color='gray', label='CL_tat')  # Línea constante
ax1.tick_params(axis='y', labelcolor=color)
ax1.legend(loc='upper left')

ax2 = ax1.twinx()
color = 'tab:blue'
ax2.set_ylabel('Error CL (%)', color=color)
ax2.plot(N, error[0, :], label='Error CL', color=color)
ax2.tick_params(axis='y', labelcolor=color)
ax2.legend(loc='upper right')

plt.title('CL Convergence Plot with Error')
plt.show()

# CM convergence plot with error
fig, ax1 = plt.subplots(figsize=(8, 6))

color = 'tab:green'
ax1.set_xlabel('Panel Number')
ax1.set_ylabel('CM', color=color)
ax1.plot(N, data[2, :], label='CM', color=color)
ax1.axhline(y=CM_LE_tat, linestyle='--', color='gray', label='CM_LE_tat')  # Línea constante
ax1.tick_params(axis='y', labelcolor=color)
ax1.legend(loc='upper left')

ax2 = ax1.twinx()
color = 'tab:orange'
ax2.set_ylabel('Error CM (%)', color=color)
ax2.plot(N, error[1, :], label='Error CM', color=color)
ax2.tick_params(axis='y', labelcolor=color)
ax2.legend(loc='upper right')

plt.title('CM Convergence Plot with Error')
plt.show()


######################################################################

M = 10
th_pp = []
"""for i in range(1, M + 1):
    print(i)
    th_pp[i-1] = np.arccos(1 - 2 * p_range[i-1])
    print(th_pp[i-1])
print(th_pp)
"""
valvec = np.zeros(10)
mi_vector = np.zeros(10)
alfa_vector = np.zeros((4, 10))
for j in range(4):
    f_bucle = j * 0.02
    print(f_bucle)
    for i in range(10):
        valor = 0.1 + i * 0.0588888
        th_pp = np.arccos(1 - 2 * valor)
        print(i)
        def dz(th):
            if th < th_pp:
                z = (f_bucle / valor ** 2) * (2 * valor - 1 + np.cos(th))
            elif th_pp <= th:
                z = (f_bucle / (1 - valor) ** 2) * (2 * valor - 1 + np.cos(th))
            return z
        result3, error3 = integrate.quad(int3, 0, np.pi)
        alfa_lo_tat = -1 / np.pi * result3
        mi_vector[i] = alfa_lo_tat * 180 / np.pi
        valvec[i] = valor
    alfa_vector[j] = mi_vector
print("alfa", alfa_vector)

"""
plt.figure(figsize=(10, 6))
plt.plot(mi_vector, label='p')
plt.xlabel('P')
plt.ylabel('alfa')
"""




# Supongamos que tienes un vector de valores


# Crear el gráfico de líneas
aaa = 0
for fila in alfa_vector:
    aaa = aaa + 1
    plt.plot(valvec, fila, linestyle='-', label=f'f = {0.02*(aaa-1)}')
#plt.plot(valvec, alfa_vector, marker='o', linestyle='-', color='b', label='Mi Vector')

# Añadir etiquetas a los ejes y al gráfico
plt.xlabel('p')
plt.ylabel('alfa_lo')
plt.title('Angle zero lift en funció de p i f')
plt.legend()  # Mostrar la leyenda si se especifican etiquetas en la función plot
plt.grid(True)
plt.xlim(0.1, 0.6)
# Mostrar el gráfico
plt.show()
######################################MIQUEL#######################################
######################################################################

valvec = np.zeros(10)
mi_vector = np.zeros(10)
alfa_vector = np.zeros((4, 10))
for j in range(4):
    f_bucle = j * 0.02
    print(f_bucle)
    for i in range(10):
        valor = 0.1 + i * 0.0588888
        th_pp = np.arccos(1 - 2 * valor)
        print(i)
        def dz(th):
            if th < th_pp:
                z = (f_bucle / valor ** 2) * (2 * valor - 1 + np.cos(th))
            elif th_pp <= th:
                z = (f_bucle / (1 - valor) ** 2) * (2 * valor - 1 + np.cos(th))
            return z
        result3, error3 = integrate.quad(int3, 0, np.pi)
        alfa_lo_tat = -1 / np.pi * result3
        mi_vector[i] = alfa_lo_tat * 180 / np.pi
        valvec[i] = valor
    alfa_vector[j] = mi_vector
print("alfa", alfa_vector)



aaa = 0
for fila in alfa_vector:
    aaa = aaa + 1
    plt.plot(valvec, fila, linestyle='-', label=f'f = {0.02*(aaa-1)}')
plt.xlabel('p')
plt.ylabel('alfa_lo')
plt.title('Angle zero lift en funció de p i f')
plt.legend()
plt.grid(True)
plt.xlim(0.1, 0.6)
plt.show()

#cmo
valvec = np.zeros(10)
vector_var = np.zeros(10)
cmo_vector = np.zeros((4, 10))
for j in range(4):
    f_bucle = j * 0.02
    print(f_bucle)
    for i in range(10):
        valor = 0.1 + i * 0.0588888
        th_pp = np.arccos(1 - 2 * valor)
        print(i)
        def dz(th):
            if th < th_pp:
                z = (f_bucle / valor ** 2) * (2 * valor - 1 + np.cos(th))
            elif th_pp <= th:
                z = (f_bucle / (1 - valor) ** 2) * (2 * valor - 1 + np.cos(th))
            return z


        result1, error1 = integrate.quad(int1, 0, np.pi)
        result2, error2 = integrate.quad(int2, 0, np.pi)
        A1 = (2 / np.pi) * result1
        A2 = (2 / np.pi) * result2
        CM0_tat = (A2 - A1) * np.pi / 4
        vector_var[i] = CM0_tat
        valvec[i] = valor
    cmo_vector[j] = vector_var

#plot cmo
aaa = 0
for fila in cmo_vector:
    aaa = aaa + 1
    plt.plot(valvec, fila, linestyle='-', label=f'f = {0.02*(aaa-1)}')

plt.xlabel('p')
plt.ylabel('Cmo')
plt.title('Cmo en funció de p i f')
plt.legend()
plt.grid(True)
plt.xlim(0.1, 0.6)

plt.show()