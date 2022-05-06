#------------------------------------------------------------------------------------------------------#
# performanceBenchmark.py                                                                              #
# Authors: Kevin Howarth, Aditi Memani, Hari Raval, Taro Spirig                                        #
# Description: File for PIC algorithm benchmarking based on user-friendly input queries.               #
#              Initializes particles randomly and uses dummy values for charge.                        #
#------------------------------------------------------------------------------------------------------#

import numpy as np
import os
import PICAlgorithm


# query user for perfomance-related inputs
N_p = int(input("Please enter desired number of particles: "))
N_g = int(input("Please enter desired number of grid points: "))
dt = float(input("Please enter the desired time-step: "))
T = float(input("Please enter the desired simulation time: "))
ompNumThreads = int(input("Please enter desired number of OpenMP threads (between 1 and 32): "))

if not (1 <= ompNumThreads and ompNumThreads <= 32):
    raise ValueError("ERROR: Must specify between 1 and 32 threads for OpenMP.")

os.environ["OMP_NUM_THREADS"] = str(ompNumThreads)

# initialize dummy values for physical constants
h, bgkChg, delChg = 1, 1, -0.01

# use M4 kernel and set flag to output time information
method, output = 1, 1

# randomly initialize particle data
P_xx = np.random.uniform(-0.5, N_g - 0.5, N_p)
P_xy = np.random.uniform(-0.5, N_g - 0.5, N_p)
P_vx = np.random.uniform(-2, 2, N_p)
P_vy = np.random.uniform(-2, 2, N_p)

# run PIC simulation and return runtime information
res = PICAlgorithm.runSimulation(dt, T, bgkChg, delChg, N_g, h, method, output, P_xx, P_xy, P_vx, P_vy)
T = float(np.sum(res))

print("Total runtime [s]: ", T)
print("Runtime of individual algorithm steps (all in [s])")
print("Step 1: ", res[0])
print("Step 2: ", res[1])
print("Step 3: ", res[2])
print("Step 4: ", res[3])