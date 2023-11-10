import numpy as np

import Parameters as par

def Calc_coord_UNIFORM (cordMatrix,f,N):
    cont = 0
    step = 1/N
    for x in np.arange(0.0,1.000001,step):
        if x < f:
            cordMatrix[cont,0] = x
            cordMatrix[cont,1] = (par.f/par.p**2)*(2*par.p*x-x**2)
        elif x >= f:
            cordMatrix[cont,0] = x
            cordMatrix[cont,1] = (par.f/(1-par.p**2))*((1-2*par.p)+2*par.p*x-x**2)
        cont +=1

def Calc_coord_Cosinus (cordMatrix,f,N):
    cont = 0
    Control = True
    while Control:
        x = 0.5*(1-np.cos(np.pi*(cont)/(N)))
        if x < f:
            cordMatrix[cont,0] = x
            cordMatrix[cont,1] = (par.f/par.p**2)*(2*par.p*x-x**2)
        elif x >= f and x<0.99999999:
            cordMatrix[cont,0] = x
            cordMatrix[cont,1] = (par.f/(1-par.p**2))*((1-2*par.p)+2*par.p*x-x**2)
        else:
            cordMatrix[cont, 0] = 1
            cordMatrix[cont, 1] = (par.f / (1 - par.p ** 2)) * ((1 - 2 * par.p) -1 + 2 * par.p )
            Control = False
        cont += 1

def Calc_panel (cordMatrix,N):
    cont_i = 0
    cont_k = 0
    cache_x = 0
    cache_z = 0
    property_panel = []
    for i in cordMatrix:
        if cont_k != 0:
            for k in i:
                if cont_i == 0:
                    delta_x = cordMatrix[cont_k,cont_i] - cache_x
                    xo = cache_x
                    cache_x = cordMatrix[cont_k,cont_i]
                    cont_i = 1
                else:
                    delta_z = cordMatrix[cont_k, cont_i] - cache_z
                    zo = cache_z
                    cache_z = cordMatrix[cont_k, cont_i]
                    cont_i = 0

            panel_chord = np.sqrt(delta_z**2 + delta_x**2)
            vec_normal = (-delta_z/panel_chord,delta_x/panel_chord)
            vec_tangent = (delta_x/panel_chord,delta_z/panel_chord)

            pos_lumpedVortex = xo + (0.25*panel_chord)*np.transpose(vec_tangent)
            pos_controlpoint = xo + (0.75 * panel_chord) * np.transpose(vec_tangent)
            conjunto = np.array([np.transpose(vec_normal),np.transpose(vec_tangent),pos_lumpedVortex,pos_controlpoint])
            property_panel.append([cont_k,conjunto])

        cont_k += 1
    return property_panel

