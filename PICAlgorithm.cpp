/*-----------------------------------------------------------------------------------------------------*/
/* PICAlgorithm.cpp                                                                                    */
/* Authors: Kevin Howarth, Aditi Memani, Hari Raval, Taro Spirig                                        */
/* Description: Implements the four steps of the Particle in Cell algorithm along with helper          */
/* functions and hook to PyBind module                                                                 */
/*-----------------------------------------------------------------------------------------------------*/

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <iostream>
#include <string>
#include <math.h>
#include <chrono>
#include <omp.h>
#include <stdexcept>
#include <Eigen/Dense>
#include "particle.h"
#include "interpolationMethods.h"
#include "iterativePoissonSolvers.h"

namespace py = pybind11;

/*-----------------------------------------------------------------------------------------------------*/
// helper function to populate the Exgrid and Eygrid with the Electric field vector field components 
// given the potentialGrid containing the scalar potential field at the current timestep 
void getElectricFields(rowMajorMatrixf& potentialGrid, rowMajorMatrixf& ExGrid, rowMajorMatrixf& EyGrid) {
    const int N_g = potentialGrid.rows();

// parallelize with collapsing due to more complex inner body loop work 
#pragma omp parallel for collapse(2)
    for (int i = 0; i < N_g; i++) {
        for (int j = 0; j < N_g; j++) {
            // apply the second order-centered difference formula 
            // note that grid spacing "h" is normalized out in units 
            ExGrid(i,j) = potentialGrid(i,(j + N_g + 1) % N_g) - potentialGrid(i,(j + N_g - 1) % N_g);
            EyGrid(i,j) = potentialGrid((i + N_g + 1) % N_g,j) - potentialGrid((i + N_g - 1) % N_g,j);
        }
    }
}
/*-----------------------------------------------------------------------------------------------------*/


/*-----------------------------------------------------------------------------------------------------*/
// helper function to timestep the velocity and position of each particle using the leapfrog method. 
// Takes in the vector of all Particles populated with electric field at current particle location and timestep
void timeStepParticles(Particle* P, const int ti, const int N_g, const int N_p) {
    float xNew, yNew;

// parallelize the loop with no collapsing due to easy compiler vectorization in loop body and with dynamic
// scheduling due to possible imbalanced workload per thread when enforcing periodic boundary condition
#pragma omp parallel for schedule(dynamic, 500) private(xNew, yNew)
    for (unsigned int p = 0; p < N_p; p++) {
        // timestep the velocity
        P[p].vx[ti + 1] = P[p].vx[ti] + P[p].Ex[ti];
        P[p].vy[ti + 1] = P[p].vy[ti] + P[p].Ey[ti];

        // timestep the position and enforce any periodic wraparounds
        xNew = P[p].xx[ti] + P[p].vx[ti + 1];
        yNew = P[p].xy[ti] + P[p].vy[ti + 1];

        // enforce periodic boundary condition in x direction for each particle
        while (xNew > N_g - 0.5) {
            xNew -= N_g;
        }
        // enforce periodic boundary condition in x direction for each particle
        while (xNew <= -0.5) {
            xNew += N_g;
        }

        // enforce periodic boundary condition in y direction for each particle
        while (yNew > N_g - 0.5) {
            yNew -= N_g;
        }
        // enforce periodic boundary condition in y direction for each particle
        while (yNew <= -0.5) {
            yNew += N_g;
        }

        // update the particle after periodic boundary condition is enforced 
        P[p].xx[ti + 1] = xNew;
        P[p].xy[ti + 1] = yNew;
    }
}
/*-----------------------------------------------------------------------------------------------------*/


/*-----------------------------------------------------------------------------------------------------*/
// primary function which carries out the PIC algorithm and is the only function which is exposed to Python

// PARAMETERS: 
// const float dt (timestep size), const float T (final time to simulate to), 
// const float bkgChg (background charge value), const float delChg (charge of single particle), 
// const int N_g (number of grid points), const float h (grid spacing),  
// const int method (flag specifying which interpolation method to use. 0 for nearest neighbor, 1 for M4),
// const int output (flag specifying what to output. 0 for particle kinematic data, 1 for wall clock time information), 
// py::array_t<float>& P_xx (initial x-coordinate of all particles),  
// py::array_t<float>& P_xy (initial y-coordinate of all particles), 
// py::array_t<float>& P_vx (initial vx-coordinate of all particles), 
// py::array_t<float>& P_vy) (initial vy-coordinate of all particles)
py::array_t<float> runSimulation(const float dt, const float T, const float bkgChg, const float delChg, 
                                const int N_g, const float h, const int method, const int output,      
                                py::array_t<float>& P_xx, py::array_t<float>& P_xy, py::array_t<float>& P_vx,   
                                py::array_t<float>& P_vy) {

#ifdef _OPENMP
    // ensure an even number of grid-points if running in parallel (required for Gauss-Seidel)
    if ((N_g % 2) != 0) {
        throw std::invalid_argument("Must specify an even number of grid points when running with parallel support.");
    }
    std::cout << "Running OpenMP with " << omp_get_max_threads() << " max threads." << std::endl;
#endif

    // request buffers and pointers for particle inputs
    py::buffer_info buf_xx = P_xx.request();
    py::buffer_info buf_xy = P_xy.request();
    py::buffer_info buf_vx = P_vx.request();
    py::buffer_info buf_vy = P_vy.request();

    // allocate points for particle data
    float* ptr_xx = (float*)buf_xx.ptr;
    float* ptr_xy = (float*)buf_xy.ptr;
    float* ptr_vx = (float*)buf_vx.ptr;
    float* ptr_vy = (float*)buf_vy.ptr;

    // get total number of timesteps and particles
    int N_t = floor(T/dt);
    int N_p = buf_xx.shape[0];

    // initialize particles and load in particles
    Particle* P = new Particle[N_p]();

    // parallelize particle value updates since each particle data fill-in is independent 
    // this also enables NUMA-aware data initialization, leading to improvements in the algorithm performance
#pragma omp parallel for
    for (int i = 0; i < N_p; i++) {
        P[i] = Particle(N_t);
        P[i].xx[0] = ptr_xx[i];
        P[i].xy[0] = ptr_xy[i];
        P[i].vx[0] = ptr_vx[i];
        P[i].vy[0] = ptr_vy[i];
    }

    // initialize Eigen matrices for the charge, potential, and electric field grids.
    rowMajorMatrixf chargeGrid = rowMajorMatrixf::Constant(N_g, N_g, bkgChg);
    rowMajorMatrixf potentialGrid = rowMajorMatrixf::Zero(N_g, N_g);
    rowMajorMatrixf ExGrid = rowMajorMatrixf::Zero(N_g, N_g);
    rowMajorMatrixf EyGrid = rowMajorMatrixf::Zero(N_g, N_g);

    // initialize timestep index and number of interations to perform for the first timestep
    int ti = 0;
    int iter = 5000; 

    // set up wall clock timer variables to time each step of the algorithm 
    std::chrono::high_resolution_clock::time_point start, end;
    float times[4];
    times[0] = 0;
    times[1] = 0;
    times[2] = 0;
    times[3] = 0;
    
    while (ti < N_t - 1) {
        
        if (ti != 0) {
            chargeGrid = rowMajorMatrixf::Constant(N_g, N_g, bkgChg);
        }
        
        if (output == 1) {
            start = std::chrono::high_resolution_clock::now();
        }

        // PIC ALGORITHM STEP 1: INTERPOLATE CHARGE TO GRID 
        P2M(P, ti, delChg, chargeGrid, method, N_p, bkgChg);

        if (output == 1) {
            end = std::chrono::high_resolution_clock::now();
            times[0] += std::chrono::duration<float>(end - start).count();
            start = std::chrono::high_resolution_clock::now();
        }

        // PIC ALGORITHM STEP 2A: GET POTENTIAL FIELD
        SORIterRowMajorData(potentialGrid, chargeGrid, h, iter);

        // PIC ALGORITHM STEP 2B: GET ELECTRIC VECTOR FIELD
        getElectricFields(potentialGrid, ExGrid, EyGrid);

        if (output == 1) {
            end = std::chrono::high_resolution_clock::now();
            times[1] += std::chrono::duration<float>(end - start).count();
            start = std::chrono::high_resolution_clock::now();
        }

       // PIC ALGORITHM STEP 3: INTERPOLATE E-FIELD TO PARTICLES
        M2P(P, ti, ExGrid, EyGrid, method, N_p);

        if (output == 1) {
            end = std::chrono::high_resolution_clock::now();
            times[2] += std::chrono::duration<float>(end - start).count();
            start = std::chrono::high_resolution_clock::now();
        }

        // PIC ALGORITHM STEP 4: TIMESTEP PARTICLES
        timeStepParticles(P, ti, N_g, N_p);

        if (output == 1) {
            end = std::chrono::high_resolution_clock::now();
            times[3] += std::chrono::duration<float>(end - start).count();
        }

        // NOTE: Use 1000 iter for STEP 2A on subsequent timesteps (after 1st timestep, potentialGrid
        // will be very close to true solution and solution won't change too much betewen timesteps).
        if (ti == 0) {
            iter = 1000;
        }

        // increment timestep counter
        ti += 1;
    }

    if (output == 0) { 

        // allocate memory dynamically for particleData which contains list of particle positions
        // and velocities at each timestep
        /*
        float*** particleData = new float**[N_p];
        for (int p = 0; p < N_p; p++){
            particleData[p] = new float*[4]; 
            for (int pp = 0; pp < 4; pp++){
                // data is formatted as N_p x 4 x N_t py::array_t<float> so
                // Python receives it as an np.array
                particleData[p][pp] = new float[N_t];
            }
        }
        */
        //float* particleData = new float[int(N_p * N_t * 4)];

        py::array_t<float> particleArr = py::array_t<float>(N_p * N_t * 4);
        py::buffer_info particleDataBuffer = particleArr.request();
        float* particleData = (float*)particleDataBuffer.ptr;

        // copy data to particleData array allocated above
#pragma omp parallel for
        for (int p = 0; p < N_p; p++) {
            for (int ti = 0; ti < N_t; ti++) {
                particleData[int(p*N_t*4 + 4*ti)] = P[p].xx[ti];
                particleData[int(p*N_t*4 + 4*ti) + 1] = P[p].xy[ti];
                particleData[int(p*N_t*4 + 4*ti) + 2] = P[p].vx[ti];
                particleData[int(p*N_t*4 + 4*ti) + 3] = P[p].vy[ti];
            }
        }
      
        // delete array of particles
#pragma omp parallel for
        for (int p = 0; p < N_p; p++) {
            // clean individual particles. this must be done explicitly to avoid issues with the destructor
            // being called for a temporary local object and freeing memory that should be copied to a permanent
            // object in the particle vector.
            P[p].clean();
        }
        delete[] P;

        // return as pybindarray. the pybind11 return value policy ensures that this is done by reference
        // and that python will take care of freeing the dynamically allocated memory of particleData.
        // return py::array_t<float>(std::vector<ptrdiff_t>{int(N_p * N_t * 4)}, particleData);
        return particleArr;
    } 
    else if (output == 1) {
        // only output information about runtime
        py::array_t<float> result = py::array_t<float>(4);
        py::buffer_info bufResult = result.request();
        float* ptrResult = (float*)bufResult.ptr;
        ptrResult[0] = times[0];
        ptrResult[1] = times[1];
        ptrResult[2] = times[2];
        ptrResult[3] = times[3];
    
        // delete pointers of particles
        for (int p = 0; p < N_p; p++) {
            P[p].clean();
        }
        delete[] P;

	    return result;
    }

    // free Particle* memory in the case of bad value for "output" is given (not needed under correct usage)
    for (int p = 0; p < N_p; p++) {
        P[p].clean();
    }
    delete[] P;

    // return empty array if bad value for "output" given
    return py::array_t<float>();
}
/*-----------------------------------------------------------------------------------------------------*/


/*-----------------------------------------------------------------------------------------------------*/
// function which uses PyBind11 to expose runSimulation() C++ function to Python
PYBIND11_MODULE(PICAlgorithm, m) {
    // PYBIND11_MODULE() macro creates a function that is called for an import statement in python.

    // module docstring
    m.doc() = "Module for running a PIC simulation using C++ code."; 

    // module_::def() generates binding code exposing runSimulation() to python (note: take_ownership
    // return policy enables output to be passed by reference to Python and Python handles deleting the object)
    m.def("runSimulation", &runSimulation, py::return_value_policy::take_ownership);
}
/*-----------------------------------------------------------------------------------------------------*/
