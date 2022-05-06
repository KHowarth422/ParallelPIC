/*-----------------------------------------------------------------------------------------------------*/
/* iterativePoissonSolvers.h                                                                           */
/* Authors: Kevin Howarth, Aditi Memani, Hari Raval, Taro Spirig                                       */
/*-----------------------------------------------------------------------------------------------------*/

#ifndef POISSON
#define POISSON

#include <Eigen/Dense>
#include "particle.h"

/*-----------------------------------------------------------------------------------------------------*/
void SORIterRowMajorData(rowMajorMatrixf& x, const rowMajorMatrixf &b, const float h, const int maxItr);
/*-----------------------------------------------------------------------------------------------------*/

#endif
