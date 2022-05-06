/*-----------------------------------------------------------------------------------------------------*/
/* iterativePoissonSolvers.cpp                                                                         */
/* Authors: Kevin Howarth, Aditi Memani, Hari Raval, Taro Spirig                                       */
/* Description: Implements iterative poisson solvers for PIC algorithm Step 2                          */
/*-----------------------------------------------------------------------------------------------------*/

#include <omp.h>
#include <iostream>
#include <fstream>
#include <string.h>
#include <math.h>
#include <random>
#include <Eigen/Dense>
#include <chrono>
#include <string>
#include "particle.h"

#define _USE_MATH_DEFINES

/*-----------------------------------------------------------------------------------------------------*/
// implement the Successive Over-Relaxation (SOR) method to solve a discretization of Poisson's equation 
// with periodic boundary conditions. Note that this implementation uses the red-black ordering scheme for easy 
// parallelization free of race conditions. However, note also that this requires n to be divisible by 2.

// PARAMETERS:
// x:  float** specifying an array of size n by n containing the memory location for the solution to be stored at
// b:  float** specifying an array of size n by n containing the source term for the Poisson equation
// h:  scalar float specifying the mesh spacing of the computational grid
// n:  scalar int specifying the number of grid points in the computational grid
// maxItr:  scalar float specifying the number of iterations to use
void SORIterRowMajorData(rowMajorMatrixf& x, const rowMajorMatrixf &b, const float h, const int maxItr){

    // get the relaxation parameters and grid size
    const int n = x.rows();
    float rho = 0.999; 
    float w = 1.0;
    int k, itr;

    // The SOR method improves upon the Gauss-Seidel method by using the heuristic that continues in the direction 
    // from x_k to x_k+1 if it is a good direction to move in. The SOR parameter omega (notated w) is chosen using 
    // a formula given in (Demmel, 1997). Additionally, Chebyshev Acceleration as described in (Hockney and Eastwood, 
    // 1988) is used to update omega at each iteration such that the initial increase in error often observed with 
    // SOR is mitigated.

    // perform the iterations until max iterations is reached
    for (itr = 0; itr < maxItr; itr++) {
       
// update "BLACK" nodes in "red-black" ordering scheme using old information
#pragma omp parallel for schedule(static) private(k) 
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n / 2; j++) {
                // re-index to j so that horizontal index doesn't depend on outer loop counter variable
                k = 2*j + (i % 2);
                *(x.data() + n*i + k) = (1 - w) * (*(x.data() + n*i + k)) + w * 
                                        (*(x.data() + n*((i - 1 + n) % n) + k) + *(x.data() + n*((i + 1 + n) % n) + k)
                                          + *(x.data() + n*i + ((k - 1 + n) % n)) + *(x.data() + n*i + ((k + 1 + n) % n))
                                          + h * h * (*(b.data() + n*i + k))) * 0.25;
            }
        }

        if (itr == 0) {
            w = 1 / (1 - 0.5 * rho * rho);
        } else {
            w = 1 / (1 - 0.25 * w * rho * rho);
        }

// update "RED" nodes in "red-black" ordering scheme using new information
#pragma omp parallel for schedule(static) private(k)
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n / 2; j += 2) {
                 // re-index to j so that horizontal index doesn't depend on outer loop counter variable
                k = 2 * j + 1 - (i % 2);
                *(x.data() + n*i + k) = (1 - w) * (*(x.data() + n*i + k)) + w * 
                                        (*(x.data() + n*((i - 1 + n) % n) + k) + *(x.data() + n*((i + 1 + n) % n) + k)
                                          + *(x.data() + n*i + ((k - 1 + n) % n)) + *(x.data() + n*i + ((k + 1 + n) % n))
                                          + h * h * (*(b.data() + n*i + k))) * 0.25;
            }
        }

        w = 1 / (1 - 0.25 * w * rho * rho);
    } 
    // pass by reference the solution estimate of Ax = b obtained via Gauss-Seidel iteration
    return;
}
/*-----------------------------------------------------------------------------------------------------*/


/*-----------------------------------------------------------------------------------------------------*/
// create an example source term for Poisson's Equation. The mathematical function chosen is periodic 
// on x, y in [0, 1], which accurately represents a potential source in 2D space.
float sourceFunction(float x, float y) {

    // output the value of the source term at the point (x, y)
    return sin(2 * M_PI * x) + cos(2 * M_PI * y);
}
/*-----------------------------------------------------------------------------------------------------*/


/*-----------------------------------------------------------------------------------------------------*/
// create an exact solution term for Poisson's equation. The mathematical function chosen is periodic 
// on x, y in [0, 1], which accurately represents a potential field in 2D space.
float fieldFunction(float x, float y) {

    // value of the potential field at the point (x, y)
    return (sin(2 * M_PI * x) + cos(2 * M_PI * y)) / (4 * M_PI * M_PI);
}
/*-----------------------------------------------------------------------------------------------------*/


/*-----------------------------------------------------------------------------------------------------*/
// map an index point on the computational grid (xi) to a point in real space, where h specifies the
// mesh spacing of the computational grid. This function assumes that the real-space coordinates being 
// represented lie in [0, 1]. 
float gridIdxToRealSpace(int xi, float h) {

    // real-space location of grid point xi
    return float(xi) * h;
}
/*-----------------------------------------------------------------------------------------------------*/


/*-----------------------------------------------------------------------------------------------------*/
// test function to ensure iterative poisson solver function compiles and executes correctly
// NOTE: must include '-D MAIN' when compiling in order to run this function since main() should only be
// compiled when debugging iterativePoissonSolvers.cpp
#ifdef MAIN
int main(int argc, char* argv[]) {
    // set n, h, and iter for test case
    int n = 500;
    float h = 1 / float(n);
    int iter = 50000;

    // initialize matrices using selected size n 
    rowMajorMatrixf b = rowMajorMatrixf::Zero(n,n);
    rowMajorMatrixf uExact = rowMajorMatrixf::Zero(n,n);
    rowMajorMatrixf xCol = rowMajorMatrixf::Zero(n,n);
    rowMajorMatrixf xRow = rowMajorMatrixf::Zero(n,n);
    rowMajorMatrixf xColData = rowMajorMatrixf::Zero(n,n);
    rowMajorMatrixf xRowData = rowMajorMatrixf::Zero(n,n);

    // compute b(i,j) and uExact(i,j) for each element of the grid
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            b(i,j) = sourceFunction(gridIdxToRealSpace(i, h), gridIdxToRealSpace(j, h));
            uExact(i,j) = fieldFunction(gridIdxToRealSpace(i, h), gridIdxToRealSpace(j, h));
        }
    }

    // variables to hold the start and end time of the SORIterRowMajorData function
    std::chrono::high_resolution_clock::time_point start, end;
    // wall clock time
    float time;

    // time the SORIterRowMajorData function
    start = std::chrono::high_resolution_clock::now();
    SORIterRowMajorData(xRowData, b, h, iter);
    end = std::chrono::high_resolution_clock::now();
    time = std::chrono::duration<float>(end - start).count();

    std::cout << "SOR Function Wallclock time: " << time << " s" << std::endl;

    return 0;
}
#endif
/*-----------------------------------------------------------------------------------------------------*/
