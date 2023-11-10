#Python project regarding the resolution of the discrete vortex metod aplied in a NACA airfol
import Parameters as par
import Vortex_Iteration as vi
import numpy as np

N = par.M + 1 #Nombre de punts
N = 5
coord = np.zeros((N+1,2)) #files columnes; x y

par.Parameters_definition()
vi.Calc_coord_UNIFORM(coord,par.f,N)
print(coord)
print("next")
#vi.Calc_coord_Cosinus(coord,par.f,N)
#print(coord)
infoMatrix = vi.Calc_panel(coord,N)

print(infoMatrix) #VECTOR NORMAL, VECTOR TANGENTE,X LUMPED VORTEX, X CONTROL POINT
"""print("general")
print(infoMatrix[0])
print("concreto")
print(infoMatrix[0][2])"""
print("neeext")
coefMatrix,RHSmatrix = vi.Iteration_Process(infoMatrix,N)
print(coefMatrix)
print(RHSmatrix)
