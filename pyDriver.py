# Authors: Kyle Fridberg, Kevin Howarth, Hari Raval, Aditi Memani, Taro Spirig  #
# Course: AM205 / CS205                                                         #
# File: pyDriver.py                                                             #
# Description: Run a 2D particle in cell simulation for a desired number of     #
# cells, inputted initial conditions, and a custom grid size. This file is the  #
# only one which the user should interact with to run their simulations. The    #
# required files to run this driver script are pyClasses.py and pyStepper.py.   #
#################################################################################
from multiprocessing.sharedctypes import Value
from pyClasses import Particle2D
from pyStepper import run_simulations
import csv
import time
import os


def initiate_simulations():

    """

    Pass in input necessary to run particle in cell simulation and visualize simulation results

    ** NOTE: User is required to pass in the following inputs as requested by the program prompts **

    L: length of simulation grid (integer)
    Ng: number of cells in grid (integer)
    dt: time step size (integer or float)
    T: end simulation time
    particles: list of Particle1D objects to be simulated (list)

    Parameters
    ----------
    None

    Raises
    -------
    None

    Returns
    -------
    None

    Example Usage
    -------
    See README.md for example usage.

    """

    # read in L, Ng, dt, and T
    L = int(input("Please enter the desired grid length (integer): "))
    Ng = int(input("Please enter the desired number of grid points: "))
    dt = float(input("Please enter the desired time-step: "))
    T = float(input("Please enter the desired simulation time: "))
    simMode = int(input("Please specify the simulation mode. 0 to run Python code, 1 to run C++ code: "))
    if not (simMode == 0 or simMode == 1):
        raise ValueError("ERROR: Must specify simMode as either 0 or 1.")

    if simMode == 1:
        ompNumThreads = int(input("C++ simulation mode selected. Please enter desired number of OpenMP threads (between 1 and 32): "))

        if not (1 <= ompNumThreads and ompNumThreads <= 32):
            raise ValueError("ERROR: Must specify between 1 and 32 threads for OpenMP.")

        os.environ["OMP_NUM_THREADS"] = str(ompNumThreads)

    file_path = str(input("Please enter the filepath to the input file: "))
    print("Input file specified. Now reading in and validating input data.")

    particles = []

    # generate a list of input particles from the user formatted CSV
    with open(file_path, 'r') as csv_file:
        csv_reader = csv.reader(csv_file, delimiter=',')
        for line_count, row in enumerate(csv_reader):
            if line_count == 0:
                pass
            else:
                particles.append(Particle2D(ID=str(row[0]), x0=float(row[1]), x1=float(row[2]),
                                            v0=float(row[3]), v1=float(row[4])))

    # run the simulations
    run_simulations(L, Ng, dt, T, particles, simMode)


initiate_simulations()
