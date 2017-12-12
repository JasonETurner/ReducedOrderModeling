# Viscous Burgers Equation
# Jason Turner - Fall 2017
#
# NOTE :
# To do:
#______________________________________________________________________________#

import numpy as np
import matplotlib.pyplot as plt

def highFidelityModel (N, Re, dx, dt, tstep, L, positionVector,initialConditionFlag, EPSILON, MAXCOUNT, dd):
    # Need to add convection speed c

    def NewtonRaphson( u_n, u_nMinus1):

#______________________________________________________________________________#
        def computeResidual( u_n, u_nMinus1):

            #   This function is a subcomponent of the broader highFidelityModel
            #   function. This function's purpose is to compute the residual of a
            #   Finite Difference approximation to the viscous Burger's Equation.

            #   Initialize Stencil Matrices

            R = np.zeros(N)
            Rinertia = np.zeros(N)
            Rdiffusion = np.zeros(N)
            Rconvection = np.zeros(N)

            Rinertia = 1/dt * (u_n-u_nMinus1)

            print('\n Residual inertia term\n = ', Rinertia)

            Rdiffusion = -(1/(Re * dx**2) ) * B.dot(u_n)

            print('\n Residual Diffusion term\n = ', Rdiffusion)

            #Rconvection = (1/(2*dx)) * u_n * A.dot(u_n)

            print('\n Residual Convection term\n = ', Rconvection)

            R = Rinertia + Rdiffusion

            #   print('\n\n The residual is:  R  = ',R)
            #   Newton's method should converge by roughly double the digits per iteration,
            #   Secant method should converge at the golden ration 1.618x digits per iteration
            R_norm = np.linalg.norm(R)
            # print('\nResidual =   ', R)
            # print('\nResidual Norm=   ', R_norm)
            return( R, R_norm )

#______________________________________________________________________________#
        def computeJacobian(u_n):
            # maybe need u_old here. I need the timestep before the timestep I'm computing.
            Convection = np.zeros( (N,N) )
            for ii in range(N):
                if ii > 0 and ii < (N-1) :
                    Convection[ii,ii+1] = -u_n[ii]
                    Convection[ii,ii] = u_n[ii+1] - u_n[ii-1]
                    Convection[ii,ii-1] = -u_n[ii]

                elif ii == 0 :
                    Convection[ii,N-2] = -u_n[0]
                    Convection[ii,ii] = u_n[1]-u_n[N-2]
                    Convection[ii,ii+1] = u_n[0]

                elif ii == (N-1):
                    Convection[ii,ii-1] = -u_n[N-1]
                    Convection[ii,ii] = u_n[1]-u_n[N-2]
                    Convection[ii,1] = u_n[N-1]
                else:
                    print(i, 'erRawr')

            J = (1/dt) * np.eye(N) + ( 1/(2*dx) ) * Convection + ( 1/(Re*(dx**2)) ) * Viscous
            return(J)

#______________________________________________________________________________#

        def computeJacobianFD( u_n, u_nMinus1 ):
            # Compute jacobian by perturbing Residual in varying directions

            DELTA = np.zeros(N) #  Perturbation vector
            Jfd = np.zeros((N,N))

            for ii in range(N):
                delta = np.zeros(N) # This manually zeroing out at every iteration is quite tedious.
                delta[ii] = dd
                # print('\nComputing Residual for Forward Jacobian FD')
                (Rplus, RnormPlus) = computeResidual(u_n + delta, u_nMinus1)
                # print('\nComputing Residual for Backward Jacobian FD')
                (Rminus, RnormMinus) = computeResidual(u_n - delta, u_nMinus1)
                Jfd[:,ii] = (Rplus - Rminus) * (1/(2*dd))

            return(Jfd)

#______________________________________________________________________________#
        #   Start of computation in NRM
        #   Computes next vector in time series using the Newton Raphson method
        #   See pg 47-57

        #   Compute Stencil Matrices for first and second spatial derivatives
        #   This can be precomputed outside of the function call
        #   A is the stencil matrix for the first spatial derivative
        #   B is the stencil matrix for the second spatial derivatives
        #   ii is the count variable.
        A = np.zeros( (N,N) )
        B = np.zeros( (N,N) )
        Viscous = np.zeros( (N,N) )

        for ii in range(N):
            if ii == 0 :
                A[ii,ii+1] = 1
                A[ii,N-2] = -1

                B[ii,ii] = -2
                B[ii,ii+1] = 1
                B[ii,N-2] = 1

                Viscous[ii,N-2] = -1
                Viscous[ii,ii] = 2
                Viscous[ii,ii+1] = -1
            elif (ii > 0) and (ii < (N-1)) :
                A[ii,ii-1] = -1
                A[ii,ii+1] = 1

                B[ii,ii] = -2
                B[ii,ii-1] = 1
                B[ii,ii+1] = 1

                Viscous[ii,ii+1] = -1
                Viscous[ii,ii] = 2
                Viscous[ii,ii-1] = -1
            elif ii == (N-1) :
                A[ii,ii-1] = -1
                A[ii,1] = 1

                B[ii,1] = 1
                B[ii,ii-1] = 1
                B[ii,ii] = -2

                Viscous[ii,ii-1] = -1
                Viscous[ii,ii] = 2
                Viscous[ii,1] = -1
            else:
                print(ii,'\n\n Error')

        #   Future interesting challenge could be computing everything by vector
        #   as opposed to the current Array computation at the end.
#______________________________________________________________________________#

        count = 0
        R_norm = 1

        while (R_norm > EPSILON) and (count < MAXCOUNT):
            print('\n\n\nComputing Jacobian for NRM. Count =',count)
            if count == 0 :
                print('\n Residual Inertia term should be zeros')

            # J = computeJacobian(u_n)
            Jfd = computeJacobianFD(u_n, u_nMinus1)
            (R, R_norm) = computeResidual( u_n, u_nMinus1)
            u_n = np.linalg.solve(Jfd,-R) + u_nMinus1
            # u_nMinus1 = u_n
            count = count + 1

        if count == MAXCOUNT :
            print('\n\n*********** NEWTON RAPHSON METHOD FAILED TO CONVERGE!!!!!! *********')


        print('\nResidual Norm=   ', R_norm, '\n iteration Number = ', count)
        # print('\n\nDifferences in Jacobian analytic and FD',np.linalg.norm(J-Jfd))
        plt.plot(positionVector, u_n, positionVector, R)
        plt.show()

        return(u_n)

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

    X = np.zeros((N, tstep ))
    X[:,0] = u_init
    print('\n\n\n\nX = ', X)
    for jj in range(tstep-1):
        print('\ntime step =   ',jj,'\nComputing vector at timestep =   ',jj+1,)
        X[:,jj+1] = NewtonRaphson(X[:,jj], X[:,jj])
        # X[:,jj+1] = NewtonRaphson(X[:,jj] - 0.1* np.ones(N), X[:,jj])

    return (X);
#______________________________________________________________________________#

def main():
    #pg 47-## - Research Notebook 1
    # input parameters
    Re = 1                      # Reynold's Number - inverse weighting for viscosity
    dx = 0.02*(np.pi)               # Spacial Discretization Value                        # Convection Parameter
    dt = 0.01                        # Time Step size
    tstep = 2                      # Number of Time Steps
    L = 2*(np.pi)                   # Length of the Domain
    initialConditionFlag = 'sin'    # Sets Dictionary lookup for initial condtion type
    EPSILON = 0.0001
    MAXCOUNT = 200
    dd = 0.0001
    positionVector = np.linspace(0,L, int(L/dx))
    N = np.size(positionVector)

    X = highFidelityModel(N, Re, dx, dt, tstep, L, positionVector,initialConditionFlag, EPSILON, MAXCOUNT, dd);
    print('\n\n\nX final =',X)
    #plt.plot( positionVector, X[:,0], positionVector, X[:,tstep-1])
    #plt.show()
#______________________________________________________________________________#

main()

#                    print('ii  =   ',ii,'\nN-1 = ',N-1,'\nN-2  =   ',N-2,'\n\n')
#            print('\ncount  = ',count)
#            print('\n\n\nPASS\n\n\n')
#            print('Residual Norm =    ',Rnorm,'\n')
