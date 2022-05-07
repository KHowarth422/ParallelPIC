# Authors: Kyle Fridberg, Kevin Howarth, Hari Raval, Aditi Memani, Taro Spirig  #
# Course: AM205 / CS205                                                         #
# File: pyClasses.py                                                            #
# Description: Class definitions for data structures necesary for               #
# Electrostatic 2D Particle-in-Cell simulations. The file contains the          #
# Particle2D class which represents a particle with some mass, charge,          #
# and kinematic signature and a Grid2D class which represents the 2D grid in    #
# which the particles live                                                      #
#################################################################################

import numpy as np
import matplotlib.pyplot as plt

# dictionary of constants for the 2D two-particle case
C = {
    "Kb": 1.0,  # boltzmann Constant
    "eChg": -1.0,  # electron charge
    "eMass": 1.0,  # electron mass
    "eps0": 1.0,  # vacuum permittivity
    "rho0": 1.0,  # background particle density
    "T0": 1.0,  # particle distribution temperature
    "vd": 1.0  # particle distribution drift velocity
}

# update dictionary of constants to include additional three terms
C.update({"debyeLength": np.sqrt(C["eps0"] * C["Kb"] / (C["rho0"] * C["eChg"] ** 2))})
C.update({"plasmaFreq": np.sqrt(C["rho0"] * C["eChg"] ** 2 / (C["eps0"] * C["eMass"]))})
C.update({"vTh": C["debyeLength"] * C["plasmaFreq"]})


class Particle2D:

    """

    Class containing the definition of a particle in 2 dimensions; represents a
    particle with mass, charge, and kinematic signature

    Instance Variables
    ----------
    ID: string identifier / unique label of an individual particle
    X0: initial x position of particle in dimensionless units of intervals (default value 0)
    X1: initial y position of particle in dimensionless units of intervals (default value 0)
    v0: initial x velocity of particle in dimensionless units of cell widths per time-step (default value 0)
    v1: initial y velocity of particle in dimensionless units of cell widths per time-step (default value 0)
    a0: initial x acceleration of particle in dimensionless units of cell widths per time-step^2
    a1: initial y acceleration of particle in dimensionless units of cell widths per time-step^2


    Returns
    -------
    Particle2D object

    """

    def __init__(self, ID, x0=0, x1=0, v0=0, v1=0):
        # string identifier for a particle
        self.ID = ID
        # x position in dimensionless units of intervals [L * H^-1] = [m * m^-1]
        self.x_0 = [x0]
        # y position in dimensionless units of intervals [L * H^-1] = [m * m^-1]
        self.x_1 = [x1]
        # x velocity in dimensionless units of cell widths per time-step
        self.v_0 = [v0]
        # y velocity in dimensionless units of cell widths per time-step
        self.v_1 = [v1]
        # x acceleration in dimensionless units of cell widths per time-step^2
        self.a_0 = []
        # y acceleration in dimensionless units of cell widths per time-step^2
        self.a_1 = []

    def __copy__(self):

        """

        Method to create a deep copy of a Particle2D object

        Parameters
        ----------
        None

        Raises
        -------
        None

        Returns
        -------
        pNew: Particle2D object containing the same attributes which self contained

        """

        # create a deep copy of a particle
        pNew = Particle2D(ID=self.ID + "_copy")
        pNew.x_0 = self.x_0[:]
        pNew.x_1 = self.x_1[:]
        pNew.v_0 = self.v_0[:]
        pNew.v_1 = self.v_1[:]
        pNew.a_0 = self.a_0[:]
        pNew.a_1 = self.a_1[:]

        return pNew


class Grid2D:

    """

    Class containing the definition of a Grid2D class which represents the 2D grid in
    which a particle lives

    Instance Variables
    ----------
    L: length of grid for simulation
    Ng: grid spacing for simulation
    dt: time-step size
    T: ending simulation time
    Particles: list of all Particles in the grid
    Charge: dimensionless charge at all grid points
    Potential: dimensionless potential at all grid points
    EField_0: dimensionless electric field in horizontal direction at all grid points
    EField_1: dimensionless electric field in vertical direction at all grid points
    PE: potential energy at each time-step
    C: dictionary of physical values

    Returns
    -------
    Grid2D object

    """

    def __init__(self, L, Ng, dt, T):
        # grid length
        self.L = L
        # grid spacing
        self.H = L / Ng
        # number of grid points
        self.Ng = Ng
        # time-step size (plasma frequencies)
        self.dt = dt
        # ending simulation time
        self.T = T
        # list of all Particles in the grid
        self.Particles = np.array([], dtype=Particle2D)
        # dimensionless charge at all grid points
        self.Charge = np.zeros((Ng, Ng))
        # dimensionless potential at all grid points
        self.Potential = np.zeros((Ng, Ng))
        # dimensionless electric field in horizontal direction at all grid points
        self.EField_0 = np.zeros((Ng, Ng))
        # dimensionless electric field in vertical direction at all grid points
        self.EField_1 = np.zeros((Ng, Ng))
        # placeholder variable for storing finite difference matrix so it doesn't need to be constructed each timestep
        self.K = None

        # populate the dictionary with the correct values for plasmaFreqDT and qBackground
        self.C = C.copy()
        self.C.update({"plasmaFreqDT": self.dt * self.C["plasmaFreq"]})
        self.C.update({"qBackground": -self.C["plasmaFreqDT"] ** 2 / 2.})

    def addParticle(self, p):

        """

        Method to add a particle to the grid. If the particle position is outside the grid,
        adjust the position until it is on the periodic image inside the grid

        ** NOTE: The valid range of positions is -0.5 <= x < Ng - 0.5, so that the nearest
        integer to any particle is a valid grid point index **

        Parameters
        ----------
        p: particle to be added to the grid

        Raises
        -------
        None

        Returns
        -------
        None

        """

        # adjust the particle position by moving to the left/right, if necessary
        p.x_0[0] = (p.x_0[0] + 1 / 2) % self.Ng - 1 / 2
        # adjust the particle position by moving to the up/down, if necessary
        p.x_1[0] = (p.x_1[0] + 1 / 2) % self.Ng - 1 / 2

        self.Particles = np.append(self.Particles, p)

        # after adding the particle, update the related parameters in the dictionary
        self.C.update({"avgParticlesPerCell": len(self.Particles) / self.Ng ** 2})
        self.C.update(
            {"delChg": self.C["plasmaFreqDT"] ** 2 * self.Ng ** 2 / (2. * len(self.Particles))})

