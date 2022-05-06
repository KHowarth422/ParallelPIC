#------------------------------------------------------------------------------------------------------#
# performanceBenchmark.py                                                                              #
# Authors: Kevin Howarth, Aditi Memani, Hari Raval, Taro Spirig                                        #
# Description: File for PIC algorithm validation. Runs a short simulation with a small number of       #
#              particles and plots their trajectories. This is to confirm correctness to visual        #
#              inspection.                                                                             #
#------------------------------------------------------------------------------------------------------#

import numpy as np
import matplotlib.pyplot as plt
import os
import PICAlgorithm

# run with 8 threads
os.environ["OMP_NUM_THREADS"] = "8"

# set input parameters
N_p = 4
N_g = 400
dt = 0.1
T = 5

# initialize dummy values for physical constants
h, bgkChg, delChg = 1, 1, -0.1

# use M4 kernel and output particle trajectory information
method, output = 1, 0

# randomly initialize particle data
P_xx = np.random.uniform(-0.5, N_g - 0.5, N_p)
P_xy = np.random.uniform(-0.5, N_g - 0.5, N_p)
P_vx = np.random.uniform(-1, 1, N_p)
P_vy = np.random.uniform(-2, 2, N_p)

# run PIC simulation and return runtime information
result = PICAlgorithm.runSimulation(dt, T, bgkChg, delChg, N_g, h, method, output, P_xx, P_xy, P_vx, P_vy)
result = result.reshape(N_p, -1, 4)

# plot trajectories
for p in range(N_p):
    plt.plot(result[p,0,0], result[p,0,1], 'rx', markersize=8, markeredgewidth=5)
    plt.plot(result[p,1:,0], result[p,1:,1], 'o', markersize=3)

plt.title("Particle Trajectories")
plt.xlabel("X")
plt.ylabel("Y")
plt.grid()
plt.xlim(-0.5, N_g-0.5)
plt.ylim(-0.5, N_g-0.5)
plt.savefig('particleTrajectories.png', bbox_inches='tight')