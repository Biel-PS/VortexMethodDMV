import Parameters as par
import Vortex_Iteration as vi
import matplotlib.pyplot as plt
import numpy as np

N = par.M; #Nombre de punts

#Loads calculation
CL = 0.0;
CMLE = 0.0;

for i in range (1,par.M,1):
