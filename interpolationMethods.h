/*-----------------------------------------------------------------------------------------------------*/
/* interpolationMethods.h                                                                              */
/* Authors: Kevin Howarth, Aditi Memani, Hari Raval, Taro Spirig                                       */
/* Description: Header file that declares the interpolation functions for interpolationMethods.cpp     */
/*-----------------------------------------------------------------------------------------------------*/

#ifndef INTERP
#define INTERP

#include <Eigen/Dense>

/*-----------------------------------------------------------------------------------------------------*/
void nearestGridPoint_P2M(float x, float y, float dq, rowMajorMatrixf& chargeGrid);
/*-----------------------------------------------------------------------------------------------------*/


/*-----------------------------------------------------------------------------------------------------*/
void nearestGridPoint_M2P(float x, float y, rowMajorMatrixf& ExGrid, rowMajorMatrixf& EyGrid,
                          float &Ex_particle, float &Ey_particle);
/*-----------------------------------------------------------------------------------------------------*/


/*-----------------------------------------------------------------------------------------------------*/
void M4_P2M(float x, float y, float dq, rowMajorMatrixf& chargeGrid);
/*-----------------------------------------------------------------------------------------------------*/


/*-----------------------------------------------------------------------------------------------------*/
void M4_M2P(float x, float y, rowMajorMatrixf& ExGrid, rowMajorMatrixf& EyGrid, float &Ex_particle, 
            float &Ey_particle);
/*-----------------------------------------------------------------------------------------------------*/


/*-----------------------------------------------------------------------------------------------------*/
void P2M(Particle* P, int ti, const float dq, rowMajorMatrixf& chargeGrid, const int method, 
        const int N_p, const float bkgChg); 
/*-----------------------------------------------------------------------------------------------------*/


/*-----------------------------------------------------------------------------------------------------*/
void M2P(Particle* P, int ti, rowMajorMatrixf& ExGrid, rowMajorMatrixf& EyGrid, const int method, 
        const int N_p);
/*-----------------------------------------------------------------------------------------------------*/

#endif
