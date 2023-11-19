import numpy as np

import Parameters as par

def Calc_coord_UNIFORM (cordMatrix,p,N): #Create the node points using a linal aprixmation
    cont = 0
    step = 1/N
    for x in np.arange(0.0,1.000001,step):
        if x < p:
            cordMatrix[cont,0] = x
            cordMatrix[cont,1] = (par.f/par.p**2)*(2*par.p*x-x**2)
        elif x >= p:
            cordMatrix[cont,0] = x
            cordMatrix[cont,1] = (par.f/(1-par.p**2))*((1-2*par.p)+2*par.p*x-x**2)
        cont +=1

def Calc_coord_Cosinus (cordMatrix,p,N): #Create the node points using an angular cosinus aproximation
    cont = 0
    Control = True
    while Control:
        x = 0.5*(1-np.cos(np.pi*(cont)/(N)))
        if x < p: #First parabola deffined by the naca 4 digit standard
            cordMatrix[cont,0] = x
            cordMatrix[cont,1] = (par.f/par.p**2)*(2*par.p*x-x**2)
        elif x >= p and x<0.99999999:#second parabola deffined by the naca 4 digit standard
            cordMatrix[cont,0] = x
            cordMatrix[cont,1] = (par.f/(1-par.p**2))*((1-2*par.p)+2*par.p*x-x**2)
        else: #Last node is TE of the airfoil
            cordMatrix[cont, 0] = 1
            cordMatrix[cont, 1] = (par.f / (1 - par.p ** 2)) * ((1 - 2 * par.p) -1 + 2 * par.p )
            Control = False
        cont += 1

def Calc_panel (cordMatrix,N): #Define the panel as a plain segmen between two nodes
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

            vec_normal = (-delta_z/panel_chord,delta_x/panel_chord) #vector normal to the segment

            vec_tangent = (delta_x/panel_chord,delta_z/panel_chord)#vector tangent to the segment

            pos_lumpedVortex = [xo + (0.25*panel_chord)*vec_tangent[0],vec_tangent[1]+zo] #postion of the vortex in one segment (0,25 of the panel lenght)
            pos_controlpoint = [xo + (0.75 * panel_chord)*vec_tangent[0],vec_tangent[1]+zo] #postion of the control point of the segment (0,75 of the panel lenght)

            conjunto = np.array([np.transpose(vec_normal),np.transpose(vec_tangent),pos_lumpedVortex,pos_controlpoint])
            property_panel.append(conjunto) #create a matrix with all the paramenters of every panel

        cont_k += 1
    return property_panel #retrun property matrix

def Iteration_Process(panelMatrix, N): #proces to obtain the matrixes of the circulation equation
    a = np.zeros((len(panelMatrix),len(panelMatrix)))
    RHS = np.zeros((len(panelMatrix),1))
    angle = [np.sin(par.alfa),np.cos(par.alfa)]
    #print(angle)
    #print(f"angle: ", angle)
    for i in range(0,len(panelMatrix)):

        for j in range(0,len(panelMatrix)):
            r2 = (panelMatrix[i][3][0]-panelMatrix[j][2][0])**2+(panelMatrix[i][3][1]-panelMatrix[j][2][1])**2 #r^2 = (xi-xj)^2 +(zi-zj)^2
            u = (panelMatrix[i][3][1]-panelMatrix[j][2][1])/(2*np.pi*r2)
            w = -(panelMatrix[i][3][0]-panelMatrix[j][2][0])/(2*np.pi*r2)
            velocity = np.array([u,w])
            a[i][j] = np.dot(velocity,panelMatrix[i][0])#*np.pi*(1/N) #fila i ,columna j. si aquest fragment es multiplica, es compleix l'exemple de la diapo 12
        RHS[i] = -np.dot(angle,panelMatrix[i][0]) *np.sin(par.alfa)
        #print(f"normal: ", panelMatrix[i][0])

    return a,RHS #return matrix of parameters and result
def Circuilation_Calc (a,RHS): #obtain the circulation matrix solving the matricial equation
   # print(np.matmul(np.linalg.inv(a), RHS))
    return np.matmul(np.linalg.inv(a),RHS)

def Lift_Coeficient (circulation): #Calcule the cl using the circulation matrix
    Cl = 0
    for i in range(0,len(circulation)):
        Cl += 2*circulation[i]
    return Cl
