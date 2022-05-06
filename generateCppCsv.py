#------------------------------------------------------------------------------------------------------#
# generateCppCsv.py                                                                                  #
# Authors: Kevin Howarth, Aditi Memani, Hari Raval, Taro Spirig                                        #
# Description: Helper file to generate large csv files as inputs. In principle a user can make a       #
#              custom input file to use their specified particle configurations. This is just an       #
#              example function which generates random initial conditions as proof of concept.         #
#------------------------------------------------------------------------------------------------------#

import csv
import numpy as np

N_p = int(input("Please enter desired number of particles:"))
N_g = int(input("Please enter desired number of grid points:"))
fname = str(input("Please enter desired filename (omitting extension):"))
fname = fname + ".csv"
header = ["ParticleID", "Position X", "Position Y", "Velocity X", "Velocity Y"]

with open(fname, 'w', encoding='UTF8', newline='') as f:
    writer = csv.writer(f)

    # write the header
    writer.writerow(header)

    # write all particles
    # particles are placed at random locations within the grid
    # with velocities randomly chosen between -2 and 2
    for p in range(1, N_p+1):
        writer.writerow([p,      
                         np.random.uniform(-0.5, N_g - 0.5), 
                         np.random.uniform(-0.5, N_g - 0.5), 
                         np.random.uniform(-2, 2), 
                         np.random.uniform(-2, 2)])