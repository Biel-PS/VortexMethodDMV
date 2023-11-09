#Python project regarding the resolution of the discrete vortex metod aplied in a NACA airfol
import Parameters as par
import Vortex_Iteration as vi
import numpy as np

N = par.M + 1 #Nombre de punts
N = 10
coord = np.zeros((N+1,2)) #files columnes; x y

par.Parameters_definition()
vi.Calc_coord_Cosinus(coord,par.f,N)
print(coord)
print("next")
#vi.Calc_coord_Cosinus(coord,par.f,N)
#print(coord)
infoMatrix = vi.Calc_panel(coord,N)
print(infoMatrix[0]) #VECTOR NORMAL, VECTOR TANGENTE,X LUMPED VORTEX, X CONTROL POINT
