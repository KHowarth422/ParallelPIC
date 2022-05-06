#------------------------------------------------------------------------------------------------------#
# runPICAlgorithm.py                                                                                   #
# Authors: Kevin Howarth, Aditi Memani, Hari Raval, Taro Spirig                                        #
# Description: Development file for testing PIC algorithm performance. Not intended for use by users.  #
#------------------------------------------------------------------------------------------------------#
import PICAlgorithm
import numpy as np
import matplotlib.pyplot as plt
import time

#------------------------------------------------------------------------------------------------------#
# Python driver function to launch the PIC algorithm with data and output kinematic data and/or 
# time statistics
def main():

    # create some particle initial conditions and parameters to run the code with 
    dt, T = 0.1, 3
    bgkChg, delChg = 1, -0.5
    N_g, h = int(200), 1

    # method to use for interpolation (0 --> nearest-neigbor, 1 --> M4)
    method = 1 

    # 0 to output particle kinematic data, 1 to output proportionate time information, 
    # 2 to output raw time information
    output = 0 
    N_p = int(6)

    # randomly initialize particle data
    P_xx = np.random.uniform(-0.5, N_g - 0.5, N_p)
    P_xy = np.random.uniform(-0.5, N_g - 0.5, N_p)
    P_vx = np.random.uniform(-2, 2, N_p)
    P_vy = np.random.uniform(-2, 2, N_p)

    # run the PIC algorithm "normally" without processing time information
    if output == 0:

        #N_p = 5
        #P_xx = np.array([0.25, 0.4, 0.5, 0.6, 0.75]) * N_g
        #P_xy = np.array([0.25, 0.25, 0.5, 0.6, 0.6]) * N_g 
        #P_vx = np.array([-1.5, -3., 2.0, 0.1, 1])
        #P_vy = np.array([1., -0.1, 2.0, 0.1, -0.5])

        print("Now calling PICAlgorithm.runSimulation().")
        start = time.time()
        result = PICAlgorithm.runSimulation(dt, T, bgkChg, delChg, N_g, h, method, output, P_xx, P_xy, P_vx, P_vy)
        end = time.time()
        print("Calling PICAlgorithm.runSimulation() was successful.")
        print("Runtime [s]:", end-start)
        print("type(result):",type(result))
        print("result.shape:",result.shape)
        #print("before reshape:")
        #print(result)
        result = result.reshape(N_p, -1, 4)
        #print("particle 0 data:")
        #print(result[0])
        

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

    # run the PIC algorithm to extract proportion of time for each step of the algorithm
    elif output == 1:
        # use 100 timesteps to get a good average
        T = 10 

        N_g1, N_p1 = int(4e2), int(1e6)
        # randomly initialize particle data
        P_xx = np.random.rand(N_p1) * (N_g1 - 0.5)
        P_xy = np.random.rand(N_p1) * (N_g1 - 0.5)
        P_vx = np.random.rand(N_p1) * 3
        P_vy = np.random.rand(N_p1) * 3
        # run PIC algorithm 
        res1 = PICAlgorithm.runSimulation(dt, T, bgkChg, delChg, N_g1, h, method, output, P_xx, P_xy, P_vx, P_vy)
        T1 = float(np.sum(res1))
        # normalize to the proportion of total time
        res1 = res1 / T1  

        N_g2, N_p2 = int(1e2), int(1e6)
        # randomly initialize particle data
        P_xx = np.random.rand(N_p2) * (N_g2 - 0.5)
        P_xy = np.random.rand(N_p2) * (N_g2 - 0.5)
        P_vx = np.random.rand(N_p2) * 3
        P_vy = np.random.rand(N_p2) * 3
        # run PIC algorithm 
        res2 = PICAlgorithm.runSimulation(dt, T, bgkChg, delChg, N_g2, h, method, output, P_xx, P_xy, P_vx, P_vy)
        T2 = float(np.sum(res2))
         # normalize to the proportion of total time
        res2 = res2 / T2 

        N_g3, N_p3 = int(2e2), int(1e5)
        # randomly initialize particle data
        P_xx = np.random.rand(N_p3) * (N_g - 0.5)
        P_xy = np.random.rand(N_p3) * (N_g - 0.5)
        P_vx = np.random.rand(N_p3) * 3
        P_vy = np.random.rand(N_p3) * 3
        # run PIC algorithm
        res3 = PICAlgorithm.runSimulation(dt, T, bgkChg, delChg, N_g3, h, method, output, P_xx, P_xy, P_vx, P_vy)
        T3 = float(np.sum(res3))
        # normalize to the proportion of total time
        res3 = res3 / T3 

        # labeled data and time information to generate proportions of time plot
        labels = ["$N_g$ = "+format(N_g1,'.0E')+ ",\n$N_p$ = "+format(N_p1, '.0E'), 
        "$N_g$ = "+format(N_g2,'.0E')+",\n$N_p$ = "+format(N_p2, '.0E'),
         "$N_g$ = "+format(N_g3, '.0E')+",\n$N_p$ = "+format(N_p3, '.0E')]
        step1times = np.array([res1[0], res2[0], res3[0]], dtype=float)
        step2times = np.array([res1[1], res2[1], res3[1]], dtype=float)
        step3times = np.array([res1[2], res2[2], res3[2]], dtype=float)
        step4times = np.array([res1[3], res2[3], res3[3]], dtype=float)
        width = 0.7

        # generate proportion of time plot 
        fig, ax = plt.subplots()
        ax.bar(labels, step1times, width, label="Step 1")
        ax.bar(labels, step2times, width, bottom=step1times, label="Step 2")
        ax.bar(labels, step3times, width, bottom=step2times+step1times, label="Step 3")
        rects = ax.bar(labels, step4times, width, bottom=step3times+step2times+step1times, label="Step 4")

        # write times above each bar chart in the plot 
        times = [T1, T2, T3]
        for i in range(len(rects)):
            height = rects[i].get_height()
            ax.text(rects[i].get_x() + rects[i].get_width()/2., 1.06,
                    '%.2f s' % times[i],
                    ha='center', va='top')

        ax.set_ylabel("Proportion of Time")
        ax.set_ylim(0, 1.08)
        ax.set_title("Bottleneck Analysis Based on Inputs")
        ax.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
        plt.rcParams.update({'font.size': 20})
        plt.savefig('timeProportion.png', bbox_inches='tight')

    # run the PIC algorithm to extract raw time for each step of the algorithm
    elif output == 2:
        # use 100 timesteps to get a good average
        T = 10  

        N_g1, N_p1 = int(1.5e3),  int(3e7)
        # randomly initialize particle data
        P_xx = np.random.rand(N_p1) * (N_g1 - 0.5)
        P_xy = np.random.rand(N_p1) * (N_g1 - 0.5)
        P_vx = np.random.rand(N_p1) * 3
        P_vy = np.random.rand(N_p1) * 3
        # run PIC algorithm
        res1 = PICAlgorithm.runSimulation(dt, T, bgkChg, delChg, N_g1, h, method, 1, P_xx, P_xy, P_vx, P_vy)
        T1 = float(np.sum(res1))

        print("Total runtime: ", T1)
        print("Step 1: ", res1[0])
        print("Step 2: ", res1[1])
        print("Step 3: ", res1[2])
        print("Step 4: ", res1[3])
#------------------------------------------------------------------------------------------------------#


if __name__ == '__main__':
    main()
