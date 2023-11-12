#Python project regarding the resolution of the discrete vortex metod aplied in a NACA airfol
import Parameters as par
import Vortex_Iteration as vi
import numpy as np

N = par.M + 1 #Nombre de punts
coord = np.zeros((N+1,2)) #files columnes; x y
par.Parameters_definition()
vi.Calc_coord_Cosinus(coord,par.f,N)
#vi.Calc_coord_Cosinus(coord,par.f,N)
#print(coord)
infoMatrix = vi.Calc_panel(coord,N)#VECTOR NORMAL, VECTOR TANGENTE,X LUMPED VORTEX, X CONTROL POINT

coefMatrix,RHSmatrix = vi.Iteration_Process(infoMatrix,N)

Circulation = vi.Circuilation_Calc(coefMatrix,RHSmatrix)

Cl = vi.Lift_Coeficient(Circulation)

"""print(f"infomatrix: ",infoMatrix[0])
print(f"A: ", coefMatrix)
print(f"Circulacion: ",Circulation)"""
print(f"Cl: ",Cl[0])
