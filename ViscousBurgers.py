# Viscous Burgers Equation
# Jason Turner - Fall 2017
#
# NOTE :
# To do:
#______________________________________________________________________________#

import numpy as np
import matplotlib.pyplot as plt

def highFidelityModel (Re, dx, dt, tstep, L, positionVector,initialConditionFlag, EPSILON):

    # Initialize Matrices
    N = np.size(positionVector)
    A = np.zeros( (N,N) )
    B = A

    # Construct Stencil Matrix A - Convection Stencil
    scalar = 1/(2*dx)
    a = 1
    b = -1

    for i in range(N):
        if i > 0 and i < N-1:
            A[i,i+1] = a
            A[i,i-1] = b
        elif i == 0 :
            A[i,N-2] = b
            A[i,i+1] = a
        elif i == (N - 1):
            A[i,i-1] = b
            A[i,1] = a
        else:
            print(i, 'erRawr')
    A = scalar * A

    # Construct Stencil Matrix B - Diffusion Stencil In Progress!!!
    scalar = 1/(Re * dx**2)
    a = 1
    b = -2

    for i in range(N):
        if i > 0 and i < N-1:
            B[i,i+1] = a
            B[i,i] = b
            B[i,i-1] = a
        elif i == 0 :
            B[i,N-2] = a
            B[i,i] = b
            B[i,i+1] = a
        elif i == (N - 1):
            B[i,i-1] = a
            B[i,i] = b
            B[i,1] = a
        else:
            print(i, 'erRawr')
    B = scalar * B

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


    for i in range(1, tstep):
        # Compute Residual

        # Find root
        while eps > EPSILON

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
