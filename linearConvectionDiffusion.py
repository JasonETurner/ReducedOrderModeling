# Solves Linear Convection Diffusion Equation with 2nd Order finite differencing
# in space and 1st order in time. Then Reduced Order Modeling (ROM) techniques
# are used on the solutions of the FD Scheme. A Snapshot matrix X is generated
# and Proper Orthogonal Decomposition is executed on this to generate a reduced
# basis. Continuous ROM is executed at the PDE level. Currently, only two basis
# functions are used, but these basis functions can be evaluated analytically.
# Jason Turner - Fall 2017
#
# NOTE : The continuous ROM is broken at the 2nd to last stage of computation.
# It probably requires hacking the data structure over to numbpy away from sympy.
# I don't think sympy can handle the vectorized ODE solve command.
# To do: Build Test function... although, isn't the result kind of a test function
# anyways? Maybe try and add in an analytical solution to get "Truth".
#______________________________________________________________________________#

import numpy as np
import matplotlib.pyplot as plt

def highFidelityModel (Re, dx, c, dt, tstep, L, positionVector,initialConditionFlag):

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

def dataDriven_ROM(X, A, tstep, positionVector):
    # USE SVD for data driven ROM approach
    num_modes = 5

    (U, sigma, V) = np.linalg.svd(X)
    Ur = U[:,0:(num_modes-1)]
    (m,n) = np.shape(Ur)

    a = np.zeros((n, tstep))
    a[:,0]= (Ur.transpose()).dot(X[:,0])
    Ar = Ur.transpose().dot(A).dot(Ur)

    for i in range(0,tstep-1):
        a[:,i+1] = Ar.dot(a[:,i])
    u_hat = Ur.dot(a)

    return(u_hat)

def continuous_ROM(Re, c, L):
    import sympy as sp

    x = sp.symbols('x')
    f, g, phi1, phi2 = sp.symbols('f g phi1 phi2', cls=sp.Function)

    a = sp.Matrix([[f], [g]])

    # prescribe Fourier Modes - if you want to do more than two basis functions,
    # should probably figure out a way to use a for loop with the basis
    # functions held in either a list or an n-tuple

    phi1 = sp.sin(1*x)
    phi2 = sp.sin(2*x) # Using just 2 basis functions so far

    # a11 = sp.integrate(phi1*phi1, (x,0,L) )
    # a12 = sp.integrate(phi2*phi1, (x,0,L) )
    # a21 = a12                           # Symmetric
    # a22 = sp.integrate(phi2*phi2, (x,0,L) )

    A = sp.Matrix( [[ sp.Integral( phi1*phi1, (x,0,L) ) ,
                        sp.Integral( phi2*phi1, (x,0,L) ) ] ,
                        [ sp.Integral( phi2*phi1, (x,0,L) ) ,
                        sp.Integral( phi2*phi2, (x,0,L) ) ]] )

    # A = lambdify ( , )

    b11 = (  c * sp.Integral( sp.diff(phi1,x) * phi1, (x,0,L) ) -
            1/Re * sp.Integral( sp.diff(phi1,x,x) * phi1, (x,0,L) )  )

    b12 = (  c * sp.Integral( sp.diff(phi2,x) * phi1, (x,0,L) ) -
            1/Re * sp.Integral( sp.diff(phi2,x,x) * phi1, (x,0,L) )  )

    b21 = (  c * sp.Integral( sp.diff(phi1,x) * phi2, (x,0,L) ) -
            1/Re * sp.Integral( sp.diff(phi1,x,x) * phi2, (x,0,L) )  )

    b22 = (  c * sp.Integral( sp.diff(phi2,x) * phi2, (x,0,L) ) -
            1/Re * sp.Integral( sp.diff(phi2,x,x) * phi2, (x,0,L) )  )
    B = sp.Matrix([[b11, b12], [b21, b22]])
    A = A.doit()                # In theory I should be able to go without eval
    B = B.doit()                # -uating these integrals
    g = -1*(A**-1) * B
#    g = g*a Going to need lambdify() here.
#    diffeq = sp.Eq( a(x).diff(x), g )
#    sp.dsolve(  a(x).diff(x)+ 1*(A**-1) * B * a(x), a )
#    This is the broken line
#    return (a(x))  # This is wrong. Should be Basis matrix multiplication first
#    I only want to return the fully analytic function

def test():
    print('this isn\'t done yet.')
    return

def printing(positionVector, u):
    plt.plot( positionVector, u[:,0] )
    plt.show()
    print('this isn\'t done yet.')
    return

def main():
    #pg 30-31 - Research Notebook 1
    #input parameters
    Re = 1.0                        # Reynold's Number - inverse weighting for viscosity
    dx = 0.01*(np.pi)               # Spacial Discretization Value
    c = 1.0                         # Convection Parameter
    dt = 0.01                       # Time Step size
    tstep = 1000                    # Number of Time Steps
    L = 2*(np.pi)                   # Length of the Domain
    initialConditionFlag = 'sin'    # Sets Dictionary lookup for initial condtion type

    positionVector = np.linspace(0,L, int(L/dx))
    (X, C) = highFidelityModel(Re, dx, c, dt, tstep, L, positionVector, initialConditionFlag);
    U_ROM = dataDriven_ROM(X, C, tstep, positionVector)
    # (u_FD_ROM, A) = continuous_ROM(Re, c, L)
    # continuous_ROM(Re, c, L)
    plt.plot( positionVector, U_ROM[:,0], positionVector, U_ROM[:,99], positionVector, U_ROM[:,199], positionVector, U_ROM[:,999] )
    plt.show()

main()
