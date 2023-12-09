import numpy as np

import Parameters as par

def Calc_coord_UNIFORM (cordMatrix,p,N,xh,eta): #Create the node points using a linal aprixmation
    cont = 0
    step = 1/N
    angleMatrix = [[np.cos(par.eta), np.sin(par.eta)], [-np.sin(par.eta), np.cos(par.eta)]]
    for x in np.arange(0.0,1.000001,step):
        if x < p: #first half of the naca parabola
            cordMatrix[cont,0] = x
            cordMatrix[cont,1] = (par.f/par.p**2)*(2*par.p*x-x**2)
        elif x >= p: #second half of the naca parabola
            cordMatrix[cont,0] = x
            cordMatrix[cont,1] = (par.f/(1-par.p)**2)*(1-2*par.p+2*par.p*x-x**2)
        cont +=1
    for i in range(0,len(cordMatrix)): #We rotate the points from the flap the deflection angle
        if cordMatrix[i][0]>= xh and eta !=0:
            cordMatrix[i] = np.matmul(angleMatrix,
                                      [cordMatrix[i][0] - xh, cordMatrix[i][1]]) + [xh, 0]
def Calc_coord_Cosinus (cordMatrix,p,N,xh,eta): #Creem els nodes del perfil mirjançant l'aproximació angular
    cont,Control=0,True
    angleMatrix = [[np.cos(eta), np.sin(eta)], [-np.sin(eta), np.cos(eta)]]
    if par.p != 0:
        Fa1,Fa2 = (par.f/ par.p**2),(par.f/(1-par.p) ** 2 )
    else:
        Fa1,Fa2 = 0,0
    while Control:
        x = 0.5*(1-np.cos(np.pi*((cont)/(N))))
        if x < p: #Definim la primera part de la parabola seguint l'estandard NACA 4-digit
            cordMatrix[cont,0] = x
            cordMatrix[cont,1] =  Fa1 * (2*par.p*x - x**2)
        elif x >= p and x!=1:#Definim la segona paràbola
            cordMatrix[cont,0] = x
            cordMatrix[cont,1] =  Fa2 * ((1 - 2 * par.p) + 2 * par.p * x - x**2)
        else: #Imposem la posició de l'últim node
            cordMatrix[cont, 0] = 1
            cordMatrix[cont, 1] = Fa2 * ((1 - 2 * par.p) - 1 + 2 * par.p )
            Control = False
        cont += 1
    for i in range(0,len(cordMatrix)):#Rotem els punts respecte Xhinge (en cas de eta =! 0)
        if cordMatrix[i][0] > xh and eta !=0:
            cordMatrix[i] = np.matmul(angleMatrix,[cordMatrix[i][0]-xh,cordMatrix[i][1]]) + [xh,0]

def Calc_panel (cordMatrix,N): #Discretizem en panells el perfil
    cont_i, cont_k, cache_x, cache_z = 0,0,0,False
    panelMatrix = []
    for i in cordMatrix: #Definim la geometria de cada panell
        if cont_k != 0: # El primer punt del LE és l'origen
            for k in i:
                if cont_i == False: #Calcul de la component X
                    delta_x = k - cache_x
                    xo, cache_x, cont_i = cache_x, k, True
                else: #Calcul de la component Z
                    delta_z = k - cache_z
                    zo, cache_z, cont_i = cache_z, k, False

            panel_chord = np.hypot(delta_z,delta_x) # Definim la longitud del panell i com chord=(sqrt(deltaz^2 + delta_x^2))
            vec_normal = (-delta_z/panel_chord,delta_x/panel_chord) #vector normal al panell
            vec_tangent = (delta_x/panel_chord,delta_z/panel_chord)#vector tangent al panell
            pos_lumpedVortex = [xo + (0.25*panel_chord)*vec_tangent[0],zo + (0.25*panel_chord)*vec_tangent[1]] #Posició de l'origen de la ciruclació (0,25 de la corda del panell)
            pos_controlpoint = [xo + (0.75 * panel_chord)*vec_tangent[0],zo + (0.75*panel_chord)*vec_tangent[1]] #Posició del punt de control (0,75 de la corda del panell)
            conjunto = np.array([np.transpose(vec_normal),np.transpose(vec_tangent),pos_lumpedVortex,pos_controlpoint]) #Incorporem tota la informació de cada panell a un vector
            panelMatrix.append(conjunto) #Creem una matriu amb tota la informació de tots els panells
            #property_panel[i][j][k] per j = 0,1,2,3 => Vector Normal, Vector Tangent, Posició Vortex, Posició punt de control; respectivament (k entre 0 (x) i 1 (z))
        cont_k += 1
    return panelMatrix #retrun property matrix

def Iteration_Process(panelMatrix, N): #Càlcul de les matrius de coeficients
    a,RHS= np.zeros((len(panelMatrix),len(panelMatrix))),np.zeros((len(panelMatrix),1)) #Definim la longitud de les matrius
    angle = [np.cos(par.alfa),np.sin(par.alfa)] #Definim la matriu d'angle d'atac
    for i in range(0,len(panelMatrix)): #Iterem per tots els punts de control
        for j in range(0,len(panelMatrix)): #Iterem per tots els vortex
            r2 = (panelMatrix[i][3][0]-panelMatrix[j][2][0])**2+\
                 (panelMatrix[i][3][1]-panelMatrix[j][2][1])**2 # r^2 = (xi-xj)^2 +(zi-zj)^2
            u = (panelMatrix[i][3][1]-panelMatrix[j][2][1])/(2*np.pi*r2) #Calculem les component de velocitat
            w = -(panelMatrix[i][3][0]-panelMatrix[j][2][0])/(2*np.pi*r2)
            velocity = np.array([u,w]) #Definim el vector de velocitat
            a[i][j] = np.dot(velocity,panelMatrix[i][0])#Càlcul de la matriu de coeficients [a] adimensional
        RHS[i] = -np.dot(angle,panelMatrix[i][0]) #Càlcul de la matriu de resultats [RHS] adimensional
    return a,RHS #Retornem les matrius
def Circuilation_Calc (a,RHS): #Càlcul de la circulació
    return np.matmul(np.linalg.inv(a),RHS)

def Lift_Coeficient (circulation,infopanel): #Càlcul del Cl del perfil i del flap amb la condició de Kutta-Joukowski
    Cl,Cl_flap,Lift_flap  = 0,0,0
    for i in range(0,len(circulation)):#Càlcul de la component discreta de Cl
        Cl += circulation[i]
        if infopanel[i][3][0] > par.xh: #Càlcul de la component discreta del Clflap
            Cl_flap += circulation[i]
    if par.xh == 1: #Per prevenir indeterminació en Clflap
        Lift_flap = 0
    else:
        Lift_flap = 2 * (Cl_flap / (1 - par.xh))#Clflap en funció de la corda del flap
    return 2*Cl, Lift_flap #Cl del perfil i Cl del flap
def MomentLE_Coeficient (circulation,infopanel):#Càlcul de Cmle i Cmxh del perfil i del flap amb la condició de Kutta-Joukowski
    cmle,cmxh,Mxh = 0,0,0
    for i in range(0, len(circulation)): #Càlcul de les components discretes dels coeficients
        cmle += circulation[i] * (infopanel[i][3][0])
        if infopanel[i][3][0] >= par.xh:
            cmxh += circulation[i] * (infopanel[i][3][0] - par.xh)
    if par.xh == 1:#Condició d'indeterminació
        Mxh = 0
    else:
        Mxh = -2 * cmxh * np.cos(par.alfa) / ((1 - par.xh)**2)
    return -2 * cmle * np.cos(par.alfa), Mxh #Retorna Cmle i Cmxh del perfil i el flap respectivament

#hola