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
def Calc_coord_Cosinus (cordMatrix,p,N,xh,eta): #Create the node points using an angular cosinus aproximation
    angleMatrix = [[np.cos(par.eta),np.sin(par.eta)],[-np.sin(par.eta),np.cos(par.eta)]]
    cont = 0
    Control = True

    Fa1 = (par.f/par.p**2)
    Fa2 = (par.f/  (1-par.p) ** 2 )

    while Control:
        x = 0.5*(1-np.cos(np.pi*((cont)/(N))))
        if x < p: #First parabola deffined by the naca 4 digit standard
            cordMatrix[cont,0] = x
            cordMatrix[cont,1] =  Fa1 * (2*par.p*x - x**2)
        elif x >= p and x<0.99999999:#second parabola deffined by the naca 4 digit standard
            cordMatrix[cont,0] = x
            cordMatrix[cont,1] =  Fa2 * ((1 - 2 * par.p) + 2 * par.p * x - x**2)
        else: #Last node is TE of the airfoil
            cordMatrix[cont, 0] = 1
            cordMatrix[cont, 1] = Fa2 * ((1 - 2 * par.p) - 1 + 2 * par.p )
            Control = False
        cont += 1
    for i in range(0,len(cordMatrix)):#We rotate the points from the flap the deflection angle
        if cordMatrix[i][0] >= xh and eta !=0:
            cordMatrix[i] = np.matmul(angleMatrix,[cordMatrix[i][0]-xh,cordMatrix[i][1]]) + [xh,0]

def Calc_panel (cordMatrix,N): #Define the panel as a plain segmen between two nodes
    cont_i = 0 # we define the cache that will be used in the definition of the panels
    cont_k = 0
    cache_x = 0
    cache_z = 0
    property_panel = []
    for i in cordMatrix: #we define the geometry of each panel
        if cont_k != 0: # we dont need to operate with the first point of the extradros
            for k in i:
                if cont_i == 0: #for the x coordinate of a point
                    delta_x = k - cache_x
                    xo = cache_x
                    cache_x = k
                    cont_i = 1
                else: #for the z coordinate of a point
                    delta_z = k - cache_z
                    zo = cache_z
                    cache_z = k
                    cont_i = 0

            panel_chord = np.hypot(delta_z,delta_x) # We define the length of the panel i functon of the total airfoil chord (sqrt(deltaz^2 + delta_x^2))

            vec_normal = (-delta_z/panel_chord,delta_x/panel_chord) #vector normal to the segment

            vec_tangent = (delta_x/panel_chord,delta_z/panel_chord)#vector tangent to the segment

            pos_lumpedVortex = [xo + (0.25*panel_chord)*vec_tangent[0],zo + (0.25*panel_chord)*vec_tangent[1]] #postion of the vortex in one segment (0,25 of the panel lenght)
            pos_controlpoint = [xo + (0.75 * panel_chord)*vec_tangent[0],zo + (0.75*panel_chord)*vec_tangent[1]] #postion of the control point of the segment (0,75 of the panel lenght)

            conjunto = np.array([np.transpose(vec_normal),np.transpose(vec_tangent),pos_lumpedVortex,pos_controlpoint]) #we set all the data to an array
            property_panel.append(conjunto) #create a matrix with all the paramenters of every panel appending the created array

        cont_k += 1
    return property_panel #retrun property matrix

def Iteration_Process(panelMatrix, N): #proces to obtain the matrixes of the circulation equation
    a = np.zeros((len(panelMatrix),len(panelMatrix))) #define all the matrixes that will be used
    RHS = np.zeros((len(panelMatrix),1))
    angle = [np.cos(par.alfa),np.sin(par.alfa)] #define the atack angle vector

    for i in range(0,len(panelMatrix)): #iterate all the control points

        for j in range(0,len(panelMatrix)): #iterate all the vortices points
            r2 = (panelMatrix[i][3][0]-panelMatrix[j][2][0])**2+(panelMatrix[i][3][1]-panelMatrix[j][2][1])**2 #r^2 = (xi-xj)^2 +(zi-zj)^2
            u = (panelMatrix[i][3][1]-panelMatrix[j][2][1])/(2*np.pi*r2)
            w = -(panelMatrix[i][3][0]-panelMatrix[j][2][0])/(2*np.pi*r2)
            velocity = np.array([u,w]) #define an array with the tangent and normal velocity
            a[i][j] = np.dot(velocity,panelMatrix[i][0])#*np.pi*(1/N) #fila i ,columna j. si aquest fragment es multiplica, es compleix l'exemple de la diapo 12
        RHS[i] = -np.dot(angle,panelMatrix[i][0])
        #print(f"normal: ", panelMatrix[i][0])

    return a,RHS #return matrix of parameters and result
def Circuilation_Calc (a,RHS): #obtain the circulation matrix solving the matricial equation

    return np.matmul(np.linalg.inv(a),RHS)

def Lift_Coeficient (circulation,infopanel): #Calculate the cl using the circulation matrix
    Cl = 0
    Cl_flap = 0
    for i in range(0,len(circulation)): #add all the Cl using the adimensional kutta condition.
        Cl += circulation[i]
        if infopanel[i][2][0] > par.xh:
            Cl_flap += circulation[i]

    return 2*Cl,2*Cl_flap
def MomentLE_Coeficient (circulation,infopanel):
    cmle = 0
    cmxh = 0
    for i in range (0,len(circulation)):
        cmle += circulation[i]*(infopanel[i][3][0])
        if infopanel[i][2][0]>=par.xh:
            cmxh+= circulation[i]*(infopanel[i][3][0]-par.xh)
    return -2*cmle*np.cos(par.alfa) , -2*cmxh*np.cos(par.alfa)

