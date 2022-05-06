#------------------------------------------------------------------------------------------------------#
# iterationConvergencePlotter.py                                                                       #
# Authors: Kevin Howarth, Aditi Memani, Hari Raval, Taro Spirig                                        #
# Description: Python file to plot the particle iteratons                                              #
#------------------------------------------------------------------------------------------------------#
import matplotlib.pyplot as plt
import sys

#------------------------------------------------------------------------------------------------------#
# function to generate plot of the particle iterations 
def main():
    # load in the command line arguments list 
    argList = sys.argv

    # validate the command line arguments depending on number inputted 
    if not len(argList) == 2:
        print("Must specify whether to plot the results of a -method test or a -init test.")
        print("Example usage: 'python iterationConvergencePlotter.py -method.")
        sys.exit()
    elif (not argList[1] == "-method") and (not argList[1] == "-init"):
        print("Must specify whether to plot the results of a -method test or a -init test.")
        print("Example usage: 'python iterationConvergencePlotter.py -method.")
        sys.exit()

    if argList[1] == '-init':
        file = open("iterationConvergenceInit.txt", "r")
    elif argList[1] == '-method':
        file = open("iterationConvergenceMethod.txt", "r")

    # read in the iteration values and create lists of iteration and error values
    f = file.readlines()
    n = f[0]
    iters = [float(f[i].split()[0]) for i in range(1, len(f))]
    errGS = [float(f[i].split()[1]) for i in range(1, len(f))]
    errSOR = [float(f[i].split()[2]) for i in range(1, len(f))]

    # generate a plot according to the requested input 
    if argList[1] == '-method':
        plt.semilogy(iters, errGS, '--', label="Gauss-Seidel")
        plt.semilogy(iters, errSOR, '--', label="SOR")
    elif argList[1] == '-init':
        plt.semilogy(iters, errGS, '--', label="uApprox")
        plt.semilogy(iters, errSOR, '--', label="uZeroes")

    # create the plot labels and show the iterations plot 
    plt.xlabel("Iterations")
    plt.ylabel("Sum of Squared Errors")
    plt.title("Iteration Convergence, n = "+n)
    plt.grid()
    plt.legend()
    plt.show()
#------------------------------------------------------------------------------------------------------#


if __name__ == '__main__':
    main()
    