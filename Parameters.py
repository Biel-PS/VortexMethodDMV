#Project parameters

import numpy as np
f,p,alfa,xh,eta = 0.02,0.4,4,0,5
def Parameters_definition():
    global f,p,alfa,xh,eta
    naca = (4,4,1,2) #input("Insert NACA series:")
    f = int(naca[0])/100 #Maximum camber
    p = int(naca[1])/10 #position of maximum camber
    angle_ATACK = -4.12217  #atack angle in degrees
    angle_FLAP = 0 #flap angle in degrees
    xh = 0.7  #Cordinate of the hinge in function of the chord


    eta = angle_FLAP * (np.pi/180)
    alfa = angle_ATACK * (np.pi/180) #from degrees to rad
    """ alfa = int(input("Insert angle of attack in degrees: ")) #angle of attack of the analysis, degrees
     xh = int(input("Insert the hinge position in % of chord: "))#position of the hinge
    eta = int(input("Insert the flap's angle of defletion in degrees")) #flap deflection angle degree """""

M = 110 #numer of panels (TO OPTIMIZE)
