# Authors: Kyle Fridberg, Kevin Howarth, Hari Raval, Aditi Memani, Taro Spirig  #
# Course: AM205 / CS205                                                         #
# File: pyStepper.py                                                            #
# Description: Primary algorithms for Electrostatic 2D                          #
# Particle-in-Cell simulations. This file actually runs the simulations by      #
# leveraging the classes created in pyClasses.py                                #
#################################################################################

import numpy as np
from pyClasses import Particle2D, Grid2D, C
import time
import PICAlgorithm 

def ChargeAssignmentStep(g, debug):

    """

    Method to perform dimensionless Charge Assignment (pg. 34 of Hockney & Eastwood)

    Parameters
    ----------
    g: Grid2D object to perform charge assignments in
    debug : boolean variable indicating whether to enter debug mode to step through simulation diagnostics carefully

    Raises
    -------
    None

    Returns
    -------
    None

    """

    # initialize charge density accumulators
    g.Charge = g.C["qBackground"] * np.ones_like(g.Charge)
    # accumulate scaled charge density
    for j in range(len(g.Particles)):
        # locate nearest mesh point
        px = int(np.round(g.Particles[j].x_0[-1]))
        py = int(np.round(g.Particles[j].x_1[-1]))

        if debug:
            try:
                # increment charge density
                g.Charge[py, px] += g.C["delChg"]
            except IndexError:
                print("IndexError for particle with position:", g.Particles[j].x[-1])
        else:
            g.Charge[py, px] += g.C["delChg"]

    # print sum of charges, if debug mode is activated
    if debug:
        print("Check Charge Neutrality. Sum of charges:", np.sum(g.Charge))


def get2DCenteredDifferenceMatrix(N, h, bc):

    """

    Helper function used by PoissonStep(). Function for constructing finite difference operator matrices. 
    For a 2D grid, Poisson's equation can be discretized as:

        u_(i-1,j) + u_(i+1,j) - 4u_(i,j) + u_(i,j-1) + u_(i,j-1)
        -------------------------------------------------------- = -f_(i,j)
                                  h^2

    Since i,j = [1, 2, ... N], there are N^2 grid points. This creates a system of N^2 algebraic equations representing
    the discretized Poisson equation. This can be written as a large, sparse linear system of the form Ax = b. The A
    matrix contains the numerical coefficients of the second-order centered difference stencil. This function constructs
    this matrix for a Poisson grid of size N by N with either zero Dirichlet or periodic Dirichlet boundary conditions.

    Parameters
    ----------
    N:  Scalar int specifying the dimension of the computational grid
    h:  Scalar float specifying the mesh spacing of the computational grid
    bc:  String specifying which boundary condition to use. Available cases include
                - "Poisson2DZeroBoundary" - zero Dirichlet boundary conditions
                - "Poisson2DPeriodic" - periodic Dirichlet boundary conditions

    Raises
    -------
    ValueError:  ValueError is raised if an invalid boundary condition is supplied

    Returns
    -------
    A:  an N^2 by N^2 2D array containing the sparse finite differencing matrix.

    """

    if bc == "Poisson2DZeroBoundary":
        # For the zero boundary conditions we can use a nice trick involving kronecker products to construct the matrix
        A = -2 * np.eye(N) + np.diag(np.ones(N - 1), k=1) + np.diag(np.ones(N - 1), k=-1)
        A = np.kron(np.eye(N), A) + np.kron(A, np.eye(N))
    elif bc == "Poisson2DPeriodic":
        # In the periodic case, the matrix is still very sparse, but each row always has five nonzero entries, as a
        # point on the grid always has four orthogonal neighbors due to periodic wrapping.
        A = np.zeros((N**2, N**2))

        # mapping function from 2D indexing to 1D indexing
        def M(x, y):
            return y * N + x

        # loop through 2D coordinates and populate finite difference A matrix
        for i in range(N):
            for j in range(N):
                Ai = M(i, j)
                A[Ai, M(i, j - 1) % (N ** 2)] = 1.
                A[Ai, M(i - 1, j) % (N ** 2)] = 1.
                A[Ai, M(i, j)] = -4.
                A[Ai, M(i, j + 1) % (N ** 2)] = 1.
                A[Ai, M(i + 1, j) % (N ** 2)] = 1.
    else:
        raise ValueError(("Please provide a valid boundary condition case. Valid cases include 'Poisson2DZeroBoundary'"
                          " and 'Poisson2DPeriodic'."))

    return (1 / h ** 2) * A


def PoissonStep(g):
  
    """

    Method to solve Poisson's equation to get potential at every point on the grid (pg. 35 of Hockney & Eastwood)

    ** NOTE 1: the reference potential is chosen such that the potential is 0 at the 0th grid point **

    ** NOTE 2: we use 0-based indexing so that we are representing mesh points [1 to Ng-1] (unlike Fortran)

    Parameters
    ----------
    g: Grid2D object to perform poisson steps on

    Raises
    -------
    None

    Returns
    -------
    None

    """

    # finite difference matrix for 2D Poisson equation
    if not g.K:
        g.K = get2DCenteredDifferenceMatrix(g.Ng, 1, bc="Poisson2DPeriodic")

    # compute the updated potential by directly solving the linear system
    rho = g.Charge.flatten()
    phi = np.linalg.solve(g.K, rho)
    g.Potential = np.reshape(phi, (g.Ng, g.Ng))


def EFieldStep(g):

    """

    Method to calculate the electric field at every point on the grid using known potentials

    Parameters
    ----------
    g: Grid2D object to perform electric field calculations on

    Raises
    -------
    None

    Returns
    -------
    None

    """

    # calculate electric field at every point on the grid using known potentials (Eq. 2-34, pg 32 of Hockney & Eastwood)
    for i in range(g.Ng):
        for j in range(g.Ng):
            g.EField_0[i, j] = g.Potential[i, (j + 1) % g.Ng] - g.Potential[i, (j - 1) % g.Ng]

    # calculate electric field at every point on the grid using known potentials (Eq. 2-34, pg 32 of Hockney & Eastwood)
    for i in range(g.Ng):
        for j in range(g.Ng):
            g.EField_1[i, j] = g.Potential[(i + 1) % g.Ng, j] - g.Potential[(i - 1) % g.Ng, j]


def ForceInterpStep(g):

    """

    Method to calculate the dimensionless force interpolation for each particle

    Parameters
    ----------
    g: Grid2D object to perform force interpolation calculations on

    Raises
    -------
    None

    Returns
    -------
    None

    """

    # iterate over each particle in the grid to compute the force interpolation
    for prt in g.Particles:
        # get the nearest mesh point
        px = int(np.round(prt.x_0[-1]))
        py = int(np.round(prt.x_1[-1]))

        # extract the electric field at that nearest mesh point
        prt.a_0.append(g.EField_0[py, px])
        prt.a_1.append(g.EField_1[py, px])


def vStep(g):

    """

    Method to calculate the time-step velocity for each particle

    Parameters
    ----------
    g: Grid2D object to perform force interpolation calculations on

    Raises
    -------
    None

    Returns
    -------
    None

    """

    # iterate over each particle in the grid to compute the time-step velocity
    for prt in g.Particles:
        prt.v_0.append(prt.v_0[-1] + prt.a_0[-1])
        prt.v_1.append(prt.v_1[-1] + prt.a_1[-1])


def xStep(g):

    """

    Method to calculate the time-step position for each particle

    Parameters
    ----------
    g: Grid2D object to perform force interpolation calculations on

    Raises
    -------
    None

    Returns
    -------
    None

    """

    # iterate over each particle in the grid to compute the time-step position
    for prt in g.Particles:
        # enforce periodicityA
        prt.x_0.append(prt.x_0[-1] + prt.v_0[-1])
        prt.x_1.append(prt.x_1[-1] + prt.v_1[-1])

        # translate the position outside the ends of the grid
        prt.x_0[-1] = (prt.x_0[-1] + 1 / 2) % g.Ng - 1 / 2
        prt.x_1[-1] = (prt.x_1[-1] + 1 / 2) % g.Ng - 1 / 2


def DiscreteModelStep(g, debug):

    """

    Method to perform a single step through the discrete model on Grid g

    Parameters
    ----------
    g: Grid2D object to perform force interpolation calculations on
    debug: boolean variable indicating whether to enter debug mode to step through simulation diagnostics carefully

    Raises
    -------
    None

    Returns
    -------
    None

    """

    # print charge diagnostics if debug is set to TRUE
    if debug:
        print("Begin Charge Assignment")
    # perform charge assignments
    ChargeAssignmentStep(g, debug)

    # print field equation (poisson) diagnostics if debug is set to TRUE
    if debug:
        print("Begin Poisson Step")
        start = time.perf_counter()
    # perform poisson computations
    PoissonStep(g)
    if debug:
        end = time.perf_counter()
        print("Finished Poisson Step in "+str(end-start)+" s.")

    # print field equation (electric field) diagnostics if debug is set to TRUE
    if debug:
        end = time.perf_counter()
        print("Finished Poisson Step in " + str(end - start) + " s.")
    if debug:
        print("Begin EField Step")
    # perform electric field computations
    EFieldStep(g)

    # print force interpolation diagnostics if debug is set to TRUE
    if debug:
        print("Begin Force Interpolation")
    # perform force interpolation computations
    ForceInterpStep(g)

    # print equations of motion (velocity) diagnostics if debug is set to TRUE
    if debug:
        print("Begin vStep")
    # perform velocity computations
    vStep(g)

    # print equations of motion (position) diagnostics if debug is set to TRUE
    if debug:
        print("Begin xStep")
    # perform position computations
    xStep(g)


def RunDiscreteModel(g, debug=False):

    """

    Method to perform a time-step of the discrete model from 0 to g.T

    Parameters
    ----------
    g: Grid2D object to perform force interpolation calculations on
    debug: boolean variable indicating whether to enter debug mode to step through simulation diagnostics carefully

    Raises
    -------
    None

    Returns
    -------
    None

    """

    # initialize time to time 0
    t = 0

    # perform the time-steps until the end simulation time is met
    while t * g.dt < g.T:
        # print diagnostics, if in debug-mode
        if debug:
            print()
            print("Start step", t)
        DiscreteModelStep(g, debug)
        t += 1

        
def run_simulations(L, Ng, dt, T, particles, simMode):

    """

    Method to run a 2D particle-in-cell simulation according to user-specified inputs

    Parameters
    ----------
    L: length of simulation grid
    Ng: number of cells in grid
    dt: time step size
    T: simulation end time
    particles: list of Particle2D objects to be simulated
    simMode: int specifying whether to run python code or C++ code for the simulation algorithm
                0 - run python simulation code
                1 - run C++ simulation code

    Raises
    -------
    None

    Returns
    -------
    None

    """

    # initialize a grid with inputted conditions
    G = Grid2D(L=L, Ng=Ng, dt=dt, T=T)

    # add all particles to the grid
    for particle in particles:
        G.addParticle(particle)

    # Ensure that enough particles are present to resolve Debye shielding
    if len(G.Particles) < G.L:
        print(
            "WARNING: There may not be enough particles present to resolve Debye shielding. "
            "\n If attempting to represent plasma waves, try using more particles.")

    # run the discrete model
    start, end = None, None
    if simMode == 0:
        # run the python code - no data manipulation required
        print("Data reading and processing complete. Running Python simulation code.")
        start = time.time()
        RunDiscreteModel(G)
        end = time.time()

        print(f"2D SIMULATION COMPLETED SUCCESSFULLY in {np.around(end - start, 2)} seconds!")

        # plot the charge
        G.plotCharge()

        # plot the potential
        G.plotPotential()

        # plot the electric field
        G.plotEField()

        # plot the state of the particles
        G.plotState()

    elif simMode == 1:
        # run the C++ code - must manipulate data into proper input forms
        # eg. initial conditions must be in 1D vectors

        P_xx = np.zeros(len(G.Particles), dtype=float)
        P_xy = np.zeros(len(G.Particles), dtype=float)
        P_vx = np.zeros(len(G.Particles), dtype=float)
        P_vy = np.zeros(len(G.Particles), dtype=float)

        for p in range(len(particles)):
            P_xx[p] = particles[p].x_0[0] 
            P_xy[p] = particles[p].x_1[0] 
            P_vx[p] = particles[p].v_0[0] 
            P_vy[p] = particles[p].v_1[0] 

        print("Data reading and processing complete. Running C++ simulation code.")

        start = time.time()
        result = PICAlgorithm.runSimulation(dt,                  # timestep
                                            T,                   # end time
                                            G.C["qBackground"],  # background charge value
                                            G.C["delChg"],       # charge increment per particle
                                            Ng,                  # number of grid points
                                            1,                   # step size, normalized via units to be = 1
                                            1,                   # interpolation method, value of 1 uses M4 kernel
                                            0,                   # output specifier, value of 0 returns particle positions and velocities at all timesteps
                                            P_xx, P_xy, P_vx, P_vy) # particle initial conditions

        end = time.time()

        print(f"2D SIMULATION COMPLETED SUCCESSFULLY in {np.around(end - start, 2)} seconds!")

        # reshape 1D array to 3D array
        result = result.reshape(len(particles), -1, 4)
        print("Simulation used",result.shape[0],"particles,",Ng,"grid points, and",result.shape[1],"time steps.")
        print("Now processing output data.")

        # now copy result values into grid particles
        for p in range(len(particles)):  # for each particle...
            for ti in range(result.shape[1]):  # for each timestep...
                G.Particles[p].x_0.append(result[p,ti,0])
                G.Particles[p].x_1.append(result[p,ti,1])
                G.Particles[p].v_0.append(result[p,ti,2])
                G.Particles[p].v_1.append(result[p,ti,3])

        # delete memory we are done with now:
        del(P_xx)
        del(P_xy)
        del(P_vx)
        del(P_vy)
        del(result)

        print("Output data processed successfully.")

        # any further analysis needed specifically for the C++ code can be done here
