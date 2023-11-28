import Parameters as par
import Vortex_Iteration as vi
import matplotlib.pyplot as plt
import numpy as np

N = par.M #Nombre de punts
coord = np.zeros((N+1,2)) #files columnes; x y
par.Parameters_definition() #definim els parametres

cont = 0
eta_inicial = -10
eta_final = 20
step = 0.1
lenght = np.abs(eta_inicial/step) + np.abs(eta_final/step)


