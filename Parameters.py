#Project parameters

import numpy as np
f,p,alfa,xh,eta = 0.,0.,0.,0.,0.
def Parameters_definition():
    global f,p,alfa,xh,eta
    naca = (0,0,0,0) #input("Insert NACA series:")
    f = int(naca[0])/100 #Maximum camber
    p = int(naca[1])/10 #position of maximum camber
    alfa = 4 *(np.pi/180)

    """ alfa = int(input("Insert angle of attack in degrees: ")) #angle of attack of the analysis, degrees
     xh = int(input("Insert the hinge position in % of chord: "))#position of the hinge
    eta = int(input("Insert the flap's angle of defletion in degrees")) #flap deflection angle degree """""

M = 100 #numer of panels
