 # **Parallel Particle-in-Cell Algorithm**
#### Team 06: Kevin Howarth, Aditi Memani, Hari Raval, Taro Spirig 
#### Course: CS205
#### Professor: Fabian Wermelinger
#### Team Supervisor: Yuxin (Iris) Ye
 
## **1. Abstract:**
 
 The Particle-In-Cell (PIC) algorithm is widely used to simulate the physics of plasmas because of its O(N) scaling with the number of particles N compared to the O(N^2) runtime of a direct N-body simulation. However, as the number of particles in plasmas can exceed 10^18 particles per cubic meter, a serial implementation of the PIC algorithm becomes exceedingly time consuming and unrealistic. The PIC algorithm functions by interpolating particles to and from a grid on which it solves the Poisson equation. Each of the steps involved in this process are naturally parallelizable. We show that a parallel implementation of the PIC algorithm with 32 threads in a shared memory model can achieve up to a 21x speedup compared to a serial implementation when considering a large number of particles and grid points. We request resources from an HPC architecture to use our parallel implementation to simulate a two-stream instability in a plasma. This computation requires 217.8 core hours.
 
## **2. File Descriptions:**
 
* **PICAlgorithm.cpp:** Implements the four steps of the Particle-in-Cell algorithm along with helper functions and hook to PyBind module.
* **interpolationMethods.cpp:** Implements interpolation methods for use in PIC algorithm interpolation.
* **interpolationMethods.h:** Header file that declares the interpolation functions for interpolationMethods.cpp.
* **iterativePoissonSolvers.cpp:** Implements iterative Poisson solvers for PIC algorithm Step 2.
* **iterativePoissonSolvers.h:** Header file that declares the SOR method for iterativePoissonSolvers.cpp.
* **sample_inputs/cpp/bigCSVInput.csv:** Input file to test our algorithm in C++ using custom initial particle conditions for 10,000 particles.
* **sample_inputs/cpp/miniCSVInput.csv:** Input file to test our algorithm in C++ using custom initial particle conditions for 3 particles.
* **particle.h:** Header files that creates the custom particle data structure.
* **performanceBenchmark.py:** File for PIC algorithm benchmarking based on user-friendly input queries. Initializes particles randomly and uses dummy values for charge.
* **Makefile:** Declares compilers and flags necessary to compile and run both the sequential and the parallel version of our program.
* **generateCppCsv.py:** Helper file to generate large csv files as inputs. In principle a user can make a custom input file to use their specified particle configurations. This is just an example function which generates random initial conditions as proof of concept.
* **LICENSE:** States our MIT License for our program.
* **pyClasses.py:** Class definitions for AM 205 data structures necessary for Python Electrostatic 2D Particle-in-Cell simulations.
* **pyDriver.py:** Driver file to run either Python or C++ 2D Particle-in-Cell simulation for a desired number of particles, inputted initial conditions, and a custom grid size.
* **pyStepper.py:** Primary algorithms for Python Electrostatic 2D Particle-in-Cell simulations.
* **trajectoryPlotter.py:** File for PIC algorithm validation. Runs a short simulation with a small number of particles and plots their trajectories. This is to confirm correctness to visual inspection.
* **development_and_testing/interpolationScratch.py:** Python file to test the implementation of interpolation schemes in order to enable ease of plotting the methods with matplotlib. C++ equivalent can be found in interpolationMethods.cpp.
* **development_and_testing/iterationConvergencePlotter.py:** Python file to plot the particle iterations.
* **development_and_testing/Roofline Analysis.ipynb:** Python file to plot the Roofline Analysis.       
* **development_and_testing/runPICAlgorithm.py:** Development file for testing PIC algorithm performance. Not intended for use by users.
* **sample_inputs/python/particle_2D_input_ex_1.csv:** Input file to test our algorithm in Python using custom initial particle conditions for 3 particles.
* **sample_inputs/python/particle_2D_input_ex_swirl.csv:** Input file to test our algorithm in Python using custom initial particle conditions for 8 particles.
* **milestone/:**  All files within this directory pertain to milestones assigned over the course of this project.
  
 
 
 
 
 
 
## **3. Build Environment:**
 
Start a new session on the cluster using the Intel Broadwell node:
```sh
$ salloc -N1 -c32 broadwell -t 1:00:00
```
Load the modules for Python and Anaconda:
```sh
$ module load python/3.8.5-fasrc01
```
Create the ~/.condarc file in your home directory and add the following content:
```
pkgs_dirs:
 - ~/.conda/pkgs
envs_dirs:
 - ~/.conda/envs
 ```
Create a new environment:
```sh
$ conda create --name new_env
```
Activate your environment:
```sh
$ conda activate new_env
```
Reload Python module and install the pybind11 library:
```sh
module load python/3.8.5-fasrc01
pip install pybind11
```
Load the GCC compiler, OpenMP, and Eigen Library modules:
```sh
module load gcc/10.2.0-fasrc01 openmpi/4.1.1-fasrc01
module load eigen/3.3.7-fasrc01
```
 
 
## **4. Performance Benchmarking Analysis:**
Utilizing our code in this format will allow the developer to compare performance of the sequential version of our algorithm to variations of the parallel version of our algorithm. In this format, particle data is randomly generated by our code. If a user would like to input their own initial conditions, please see 'Custom Particle Input'.
 
 
### 4.1 Sequential Version:
To compile and run the program please run the following commands:
```sh
$ make 
$ python3 performanceBenchmark.py
```
User can then enter the number of particles, grid points, time-step size, and total simulation time.
For a small test case use the following parameters:
  * Particles: 1000000
  * Grid Points: 400
  * Time-step Size: 0.1
  * Total Time: 10
 
The expected runtime for these inputs is between 2 and 3 minutes. Upon simulation completion, a message will print to the terminal stating the total simulation runtime as well as the total amount of time taken by each individual step of the algorithm. Note that these times are computed only for the algorithm steps themselves; operations like I/O and pybind11 data handling are not included. Hence, the time between issuing and completing the simulation command may be slightly larger than the times printed.
 
### 4.2 Parallel Version with the Eigen Library:
To compile and run the program please run the following commands:
```sh
$ make ompNoNUMA
$ python3 performanceBenchmark.py
```
User can then enter the number of particles, grid points, time-step size, total simulation time, and the number of threads.
For a small test case use the following parameters:
  * Particles: 1000000
  * Grid Points: 400
  * Time-step Size: 0.1
  * Total Time: 10
  * OMP_NUM_THREADS = 16
 
The expected runtime for these inputs is between 20 and 30 seconds. As mentioned before, a message will print stating the simulation runtime and total runtime of each individual step.
 
### 4.3 Parallel Version with the Eigen Library and NUMA First Touch:
To compile and run the program please run the following commands:
```sh
$ make omp
$ python3 performanceBenchmark.py
```
User can then enter the number of particles, grid points, time-step size, total simulation time, and the number of threads.
For a large test case, which was also used for our Roofline Analysis, please use the following parameters:
  * Particles: 30000000
  * Grid Points: 1500
  * Time-step Size: 0.1
  * Total Time: 10
  * OMP_NUM_THREADS = 32
 
The expected runtime for these inputs is around 5 minutes.
 
## 5. **Custom Particle Input:**
Utilizing our code in this format will allow the user to input a .csv file with custom initial conditions for each particle. We have provided some sample .csv files to illustrate the input format located in the "sample_inputs" directory. We have also provided a script, "generateCppCsv.py" allowing generation of a .csv file with a user-specified number of grid points and particles. In principle, a user could create a custom .csv file in this style with particle initial conditions specifically constructed to illustrate some physical phenomenon.
 
 
### 5.1 Sequential Version:
To compile and run the program please run the following commands:
```sh
$ make 
$ python3 pyDriver.py
```
User can then enter the length of the grid, number of grid points, time-step size, total simulation time, run the code using Python or C++, file path to .csv file, and number of threads.
For a sample test case use the following parameters:
  * Grid Length: 100
  * Grid Points: 100
  * Time-step Size: 0.1
  * Total Time: 10
  * 0 to run Python code, 1 to run C++ code: 1
  * OMP_NUM_THREADS = 1
  * file path: ./sample_inputs/cpp/miniCSVInput.csv
 
 
The expected runtime of this case is around 10 seconds.
 
### 5.2 Parallel Version with the Eigen Library:
To compile and run the program please run the following commands:
```sh
$ make ompNoNUMA
$ python3 pyDriver.py
```
User can then enter the length of the grid, number of grid points, time-step size, total simulation time, run the code using Python or C++, file path to .csv file, and number of threads.
For a sample test case use the following parameters:
  * Grid Length: 400
  * Grid Points: 400
  * Time-step Size: 0.1
  * Total Time: 10
  * 0 to run Python code, 1 to run C++ code: 1
  * OMP_NUM_THREADS = 32
  * file path: ./sample_inputs/cpp/bigCSVInput.csv
 
The expected runtime of this case is around 2 minutes. Note that the simulation will complete in much less time than this, but the runtime is now dominated by the rudimentary I/O implementation.
 
### 5.3 Parallel Version with the Eigen Library and NUMA First Touch:
To compile and run the program please run the following commands:
```sh
$ make omp
$ python3 pyDriver.py
```
User can then enter the length of the grid, number of grid points, time-step size, total simulation time, run the code using Python or C++, file path to .csv file, and number of threads.
For a sample test case use the following parameters:
  * Grid Length: 400
  * Grid Points: 400
  * Time-step Size: 0.1
  * Total Time: 10
  * 0 to run Python code, 1 to run C++ code: 1
  * file path: ./sample_inputs/cpp/bigCSVInput.csv
  * OMP_NUM_THREADS = 32
 
The expected runtime of this case is also around 2 minutes.
 
## **6. Run Sample Trajectory Simulation:**
The is a simulation that allows the user to visualize particle trajectories using a small test case.
 
To compile and run this simulation:
```sh
$ make omp
$ python3 trajectoryPlotter.py
```
This will produce an image of the trajectories in your current directory, saved as 'particleTrajectories.png' :
![alt text](https://code.harvard.edu/CS205/team06/blob/5ea52294d1177ab2750f7b396869d1d176d73241/sample_trajectory.jpeg)
 

