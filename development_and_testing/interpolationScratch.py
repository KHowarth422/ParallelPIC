#------------------------------------------------------------------------------------------------------#
# interpolationScratch.py                                                                              #
# Authors: Kevin Howarth, Aditi Memani, Hari Raval, Taro Spirig                                        #
# Description: Python file to test the implementation of interpolation schemes in order to             #
# enable ease of plotting the methods with matplotlib. C++ equivalent can be found in                  # 
# interpolationMethods.cpp                                                                             #
#------------------------------------------------------------------------------------------------------#
import numpy as np
import math
import matplotlib.pyplot as plt


#------------------------------------------------------------------------------------------------------#
# helper function to find the nearest gridpoint (x, y) and update its charge dq in the chargeGrid
def nearest_grid_point_P2M(x, y, dq, grid):
    # round to the nearest grid point
    i, j = int(np.round(x)), int(np.round(y))
    # update the charge at the specified location
    grid[j, i] += dq
    # return the updated grid
    return grid
#------------------------------------------------------------------------------------------------------#


#------------------------------------------------------------------------------------------------------#
# interpolate from electric field on grid to electric field at Particle's location based on nearest neighbor 
def nearest_grid_point_M2P(x, y, grid_Ex, grid_Ey):
    # round to the nearest grid point
    i, j = int(np.round(x)), int(np.round(y))
    # update the charge at the specified location
    particle_Ex = grid_Ex[j,i]
    particle_Ey = grid_Ey[j,i]

    # return the updated grid
    return particle_Ex, particle_Ey
#------------------------------------------------------------------------------------------------------#


#------------------------------------------------------------------------------------------------------#
# perform interpolation from particle charges back to grid based on M4 kernel 
def m4_P2M(x, y, dq, grid):
    n = grid.shape[0]
    # account for wrap-around of points
    start_x, start_y = (math.ceil(x) - 2 + n) % n, (math.ceil(y) - 2 + n) % n
    end_x, end_y = (start_x + 4 + n) % n, (start_y + 4 + n) % n

    # perform the interpolation at each grid location
    for i in range(start_x, end_x):
        for j in range(start_y, end_y):
            grid[j, i] += (dq * m4_kernel(abs(x - i)) * m4_kernel(abs(y - j)))

    return grid
#------------------------------------------------------------------------------------------------------#


#------------------------------------------------------------------------------------------------------#
# perform interpolation from electric field on grid to electric field at particle's location based on 
# M4 kernel.
def m4_M2P(x, y, grid_Ex, grid_Ey):
    n = grid.shape[0]
    # account for wrap-around of points
    start_x, start_y = (math.ceil(x) - 2 + n) % n, (math.ceil(y) - 2 + n) % n
    end_x, end_y = (start_x + 4 + n) % n, (start_y + 4 + n) % n
    particle_Ex, particle_Ey = 0, 0

    for i in range(start_x, end_x):
        for j in range(start_y, end_y):
            particle_Ex += grid_Ex[j, i] * m4_kernel(abs(x - i)) * m4_kernel(abs(y - j))
            particle_Ey += grid_Ey[j, i] * m4_kernel(abs(x - i)) * m4_kernel(abs(y - j))

    return particle_Ex, particle_Ey
#------------------------------------------------------------------------------------------------------#


#------------------------------------------------------------------------------------------------------#
# helper function to compute the M4_kernel 
def m4_kernel(u):
    if 0 <= u <= 1:
        return 1 - (5 * u ** 2) / 2 + (3 * u ** 2) / 2
    elif 1 < u <= 2:
        return ((2 - u) ** 2 * (1 - u)) / 2

    return 0
#------------------------------------------------------------------------------------------------------#


#------------------------------------------------------------------------------------------------------#
# function to test the correctness of the Python interpolation methods 
def main():
    # initialize some conditions
    N = 30
    qBackground = 1
    grid_nearest = np.ones((N,N), dtype = float) * qBackground
    grid_M4 = np.ones((N, N), dtype=float) * qBackground

    # initialize particle conditions
    dq = -0.1
    prt_x = np.random.uniform(low=0, high=N-1, size=(20,))
    prt_y = np.random.uniform(low=0, high=N-1, size=(20,))

    # perform the interpolation
    for i in range(20):
        grid_nearest = nearest_grid_point_P2M(prt_x[i], prt_y[i], dq, grid_nearest)
        grid_M4 = m4_interpol(prt_x[i], prt_y[i], dq, grid_M4)

    # plot the results
    X = np.arange(0, N, 1)
    Y = np.arange(0, N, 1)
    X, Y = np.meshgrid(X, Y)

    # plot the results of the interpolation
    f, (ax1, ax2) = plt.subplots(1, 2,  subplot_kw={"projection": "3d"})
    ax1.plot_surface(X, Y, grid_nearest, cmap='viridis', edgecolor='green')
    ax1.set_title('Nearest grid point interp')
    ax2.plot_surface(X, Y, grid_M4, cmap='viridis', edgecolor='green')
    ax2.set_title('M4 interp')
    plt.show()

    prtE = interp(x, y, )
    prt_a.append(..)
#------------------------------------------------------------------------------------------------------#


if __name__ == "__main__":
    main()


