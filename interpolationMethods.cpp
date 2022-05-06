/*-----------------------------------------------------------------------------------------------------*/
/* interpolationMethods.cpp                                                                            */
/* Authors: Kevin Howarth, Aditi Memani, Hari Raval, Taro Spirig                                       */
/* Description: Implements interpolation methods for use in PIC algorithm interpolation                */
/*-----------------------------------------------------------------------------------------------------*/

#include <math.h>
#include <iostream>
#include <string>
#include <Eigen/Dense>
#include "particle.h"

#define _USE_MATH_DEFINES

/*-----------------------------------------------------------------------------------------------------*/
// declare a custom reduction which performs addition of Eigen::Matrix objects
#pragma omp declare reduction(MatrixPlus: rowMajorMatrixf: omp_out=omp_out+omp_in) \
                    initializer(omp_priv=rowMajorMatrixf::Zero(omp_orig.rows(),omp_orig.cols()))
/*-----------------------------------------------------------------------------------------------------*/


/*-----------------------------------------------------------------------------------------------------*/
// helper function to find the nearest gridpoint (x, y) and update its charge dq in the chargeGrid
void nearestGridPoint_P2M(float x, float y, float dq, rowMajorMatrixf& chargeGrid) {

    // round to nearest grid point
    int i = int(round(x));
    int j = int(round(y));

    // update the charge at the specified location
    const int N_g = chargeGrid.rows();
    *(chargeGrid.data() + N_g*j + i) += dq; 
}
/*-----------------------------------------------------------------------------------------------------*/


/*-----------------------------------------------------------------------------------------------------*/
// interpolate from electric field on grid to electric field at Particle's location based on nearest neighbor. 
// Note: 'E' is a vector field and hence has two components (returned by reference in Ex_particle and Ey_particle)
void nearestGridPoint_M2P(float x, float y, rowMajorMatrixf& ExGrid, rowMajorMatrixf& EyGrid, float &Ex_particle, 
                         float &Ey_particle) {

    // round to nearest grid point
    int i = int(round(x));
    int j = int(round(y));

    // update the charge at the specified location
    const int N_g = ExGrid.rows();
    Ex_particle = *(ExGrid.data() + N_g*j + i);
    Ey_particle = *(EyGrid.data() + N_g*j + i);
}
/*-----------------------------------------------------------------------------------------------------*/


/*-----------------------------------------------------------------------------------------------------*/
// helper function to compute the M4_kernel 
float M4_kernel(float u) {
    // M4 kernel operation
    if (0 <= u && u <= 1) {
        return 1 - 0.5 * (5 * u * u) + 0.5 * (3 * u * u);
    } 
    else if (1 < u && u <= 2) {
        return 0.5 * ((2 - u) * (2 - u) * (1 - u));
    }

    return 0;
}
/*-----------------------------------------------------------------------------------------------------*/


/*-----------------------------------------------------------------------------------------------------*/
// perform interpolation from particle charges back to grid based on M4 kernel 
void M4_P2M(float x, float y, float dq, rowMajorMatrixf& chargeGrid) {

    // find grid size
    int N_g = chargeGrid.rows();

    // perform interpolation to each of the 16 surrounding grid points
    int i_wrap, j_wrap;
    int x_ceil = int(ceil(x));
    int y_ceil = int(ceil(y));

    // iterate through each valid location and update the chargeGrid
    for (int i = x_ceil - 2; i < x_ceil + 2; i++) {
        for (int j = y_ceil - 2; j < y_ceil + 2; j++) {
            // ensure periodic wraparounds for grid points beyond the edges
            i_wrap = (i + N_g) % N_g;
            j_wrap = (j + N_g) % N_g;
            *(chargeGrid.data() + N_g*j_wrap + i_wrap) += dq * M4_kernel(abs(x - i)) * M4_kernel(abs(y - j));
        }
    }
}
/*-----------------------------------------------------------------------------------------------------*/


/*-----------------------------------------------------------------------------------------------------*/
// perform interpolation from electric field on grid to electric field at particle's location based on 
// M4 kernel. Note: 'E' is a vector field and hence has two components (returned by reference in Ex_particle 
// and Ey_particle)
void M4_M2P(float x, float y, rowMajorMatrixf& ExGrid, rowMajorMatrixf& EyGrid, float &Ex_particle, float &Ey_particle) {

    // find grid size
    int N_g = ExGrid.rows();

    // perform interpolation to each of the 16 surrounding grid points
    float sum_Ex = 0;
    float sum_Ey = 0;
    int i_wrap, j_wrap;
    int x_ceil = int(ceil(x));
    int y_ceil = int(ceil(y));

    // iterate through each valid location and accumlate the electric field amounts
    for (int i = x_ceil - 2; i < x_ceil + 2; i++) {
        for (int j = y_ceil - 2; j < y_ceil + 2; j++) {
            // ensure periodic wraparounds for grid points beyond the edges when indexing
            i_wrap = (i + N_g) % N_g;
            j_wrap = (j + N_g) % N_g;
            sum_Ex += *(ExGrid.data() + N_g * j_wrap + i_wrap) * M4_kernel(abs(x - i)) * M4_kernel(abs(y - j));
            sum_Ey += *(EyGrid.data() + N_g * j_wrap + i_wrap) * M4_kernel(abs(x - i)) * M4_kernel(abs(y - j));
        }
    }
    // return electric field vector components by reference
    Ex_particle = sum_Ex;
    Ey_particle = sum_Ey;
}
/*-----------------------------------------------------------------------------------------------------*/


/*-----------------------------------------------------------------------------------------------------*/

// perform interpolation of all particles to chargeGrid based on the current timestep (ti) and the
// interpolation kernel to use (method = 0 --> nearest neighbor vs method = 1 --> M4)
void P2M(Particle* P, int ti, const float dq, rowMajorMatrixf& chargeGrid, const int method, 
        const int N_p, const float bkgChg) {

// perform the custom reduction to add Eigen::Matrix objects
#pragma omp parallel for reduction(MatrixPlus: chargeGrid)
    for (unsigned int p = 0; p < N_p; p++) {
        switch(method) {
            case 0:
                // nearest neighbor
                nearestGridPoint_P2M(P[p].xx[ti], P[p].xy[ti], dq, chargeGrid);
                break;
            case 1:
                // M4
                M4_P2M(P[p].xx[ti], P[p].xy[ti], dq, chargeGrid);
                break;
        }
    }
}
/*-----------------------------------------------------------------------------------------------------*/


/*-----------------------------------------------------------------------------------------------------*/
// interpolate electric field to all particles based on the current timestep (ti) and the
// interpolation kernel to use (method = 0 --> nearest neighbor vs method = 1 --> M4)
void M2P(Particle* P, int ti, rowMajorMatrixf& ExGrid, rowMajorMatrixf& EyGrid, const int method, const int N_p) {

// parallelize without collapsing due to easy compiler vectorization in loop body
#pragma omp parallel for schedule(static)
    for (unsigned int p = 0; p < N_p; p++) {
        switch(method) {
            case 0:
                // nearest neighbor
                nearestGridPoint_M2P(P[p].xx[ti], P[p].xy[ti], ExGrid, EyGrid, P[p].Ex[ti], P[p].Ey[ti]);
                break;
            case 1:
                // M4 kernel
                M4_M2P(P[p].xx[ti], P[p].xy[ti], ExGrid, EyGrid, P[p].Ex[ti], P[p].Ey[ti]);
                break;
        }
    }
}
/*-----------------------------------------------------------------------------------------------------*/


/*-----------------------------------------------------------------------------------------------------*/
// test function to ensure interpolation methods compile and execute correctly
// NOTE: must include '-D MAIN' when compiling in order to run this function since main() should only be
// compiled when debugging interpolationMethods.cpp
#ifdef MAIN
int main(void) {

    // pick x, y, and dq values to use for test case
    float x = 2.1; 
    float y = 2.1;
    float dq = 0.1;

    // construct source term and exact solution
    int n = 100;
    float** b = new float*[n];
    float** Ex = new float*[n];
    float** Ey = new float*[n];

    // fill im x and y electric field values and b values
    for (int i = 0; i < n; i++) {
        b[i] = new float[n];
        Ex[i] = new float[n];
        Ey[i] = new float[n];
        for (int j = 0; j < n; j++) {
            b[i][j] = 1.0;
            Ex[i][j] = 1.0;
            Ey[i][j] = 1.0;
        }
    }

    // perform P2M using nearest grid point approach 
    nearestGridPoint_P2M(x, y, dq, b);
    // perform M2P using nearest grid point approach 
    float Ex_particle, Ey_particle;
    nearestGridPoint_M2P(x, y, Ex, Ey, Ex_particle, Ey_particle);

    std::cout << "Nearest Neighbor:\n"
              << "Ex_particle: " << std::to_string(Ex_particle) << std::endl 
              << "Ey_particle: " << std::to_string(Ey_particle) << "\n\n";

    // perform P2M using M4 kernel 
    M4_P2M(x, y, dq, b, n);
    // perform M2P using M4 kernel 
    M4_M2P(x, y, Ex, Ey, n, Ex_particle, Ey_particle);

    std::cout << "M4:\n"
              << "Ex_particle: " << std::to_string(Ex_particle) << std::endl
              << "Ey_particle: " << std::to_string(Ey_particle) << "\n\n";
    
    return 0;
}
#endif
/*-----------------------------------------------------------------------------------------------------*/