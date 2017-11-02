# Viscous Burgers Equation
# Jason Turner - Fall 2017
#
# NOTE :
# To do:
#______________________________________________________________________________#

import numpy as np
import matplotlib.pyplot as plt

def highFidelityModel (Re, dx, dt, tstep, L, positionVector,initialConditionFlag):

    #Construct Stencil Matrix
    aCoefficient = 1/(Re*dx**2) - c/(2*dx)
    bCoefficient = -2/(dx**2)
    cCoefficient = 1/(Re*dx**2)+c/(2*dx)

    a = np.ones((np.size(positionVector)-1,1)) * aCoefficient    #off diagonal
    b = np.ones((np.size(positionVector),1)) * bCoefficient
    c = np.ones((np.size(positionVector)-1,1)) * cCoefficient    #off diagonal

    A = np.zeros((np.size(positionVector),np.size(positionVector)))

    for i in range(np.size(positionVector)):
        if i > 0 and i < np.size(positionVector)-1:
            A[i,i] = bCoefficient
            A[i,i+1] = aCoefficient
            A[i,i-1] = cCoefficient
        elif i == 0 :
            A[i,i] = bCoefficient
            A[i,i+1] = aCoefficient
        elif i == np.size(positionVector)-1:
            A[i,i-1] = cCoefficient
            A[i,i] = bCoefficient
        else:
            print(i, 'erRawr')

    A[0,np.size(positionVector)-2] = cCoefficient
    A[np.size(positionVector)-1, 1] = aCoefficient
    I = np.eye(np.size(positionVector))   # NxN identity Matrix
    B = (I - dt * A)                # Evolution Matrix - see page 31

    # Define initial conditions as functions called in the dictionary
    # lookup executed below
    def sinCondition():
        return( np.sin(positionVector) )
    def cosCondition():
        return( np.cos(positionVector) )
    def sawToothCondition():
        u_init = positionVector
        u_init[ 0:int(0.4*L/dx) ] = 0
        u_init[ int(0.4*L/dx):int(0.5*L/dx) ] = ( 5/np.pi *
            ( positionVector[int(0.4*L/dx):int(0.5*L/dx)] - 0.4*L ) )
        u_init[ int(0.5*L/dx):int(0.6*L/dx) ] = ( -5/np.pi *
            ( positionVector[int(0.5*L/dx):int(0.6*L/dx)] - 0.5*L) + 1 )
        u_init[ int(0.6*L/dx):int(L/dx) ] = 0
        return( u_init )
    def boxCondition():
        # This initial condition will not work without upwinding, or a different
        # Discretization Scheme
        u_init = positionVector
        u_init[ 0:int(0.4*L/dx) ] = 0
        u_init[ int(0.4*L/dx):int(0.6*L/dx) ] = 1
        u_init[ int(0.6*L/dx):int(L/dx) ] = 0
        return( u_init )

    # Dictionary style lookup table for initial conditions
    initialConditionType = { 'sin' : sinCondition,
                                'saw' : sawToothCondition,
                                'cos' : cosCondition,
                                'box' : boxCondition,
    }

    u_init = initialConditionType[initialConditionFlag]()

    u = np.zeros(( np.size(u_init) , tstep ))
    u[:,0] = u_init
    #plt.plot(u_init)
    C = np.linalg.inv(B)

    for i in range(1, tstep):
        u[:,i] = C.dot(u[:,i-1])

    return (u, C);

def main():
    #pg 47-## - Research Notebook 1
    # input parameters
    Re = 1.0                        # Reynold's Number - inverse weighting for viscosity
    dx = 0.01*(np.pi)               # Spacial Discretization Value                        # Convection Parameter
    dt = 0.01                       # Time Step size
    tstep = 1000                    # Number of Time Steps
    L = 2*(np.pi)                   # Length of the Domain
#    initialConditionFlag = 'sin'    # Sets Dictionary lookup for initial condtion type

    positionVector = np.linspace(0,L, int(L/dx))
    (X, C) = highFidelityModel(Re, dx, c, dt, tstep, L, positionVector, initialConditionFlag);
    U_ROM = dataDriven_ROM(X, C, tstep, positionVector)
    # (u_FD_ROM, A) = continuous_ROM(Re, c, L)
    # continuous_ROM(Re, c, L)
    plt.plot( positionVector, U_ROM[:,0], positionVector, U_ROM[:,99], positionVector, U_ROM[:,199], positionVector, U_ROM[:,999] )
    plt.show()

main()
