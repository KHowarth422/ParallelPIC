/*-----------------------------------------------------------------------------------------------------*/
/* particle.h                                                                                          */
/* Authors: Kevin Howarth, Aditi Memani, Hari Raval, Taro Spirig                                       */
/*-----------------------------------------------------------------------------------------------------*/

#ifndef PARTICLE
#define PARTICLE
#include <vector>

/*-----------------------------------------------------------------------------------------------------*/
// a typedef for Eigen matrices that should be seen by every other file
typedef Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> rowMajorMatrixf;
/*-----------------------------------------------------------------------------------------------------*/


/*-----------------------------------------------------------------------------------------------------*/
// class that represents a Particle object 
class Particle {
public:

    int Nt;

// compiler flag to switch between float* array Particles (NUMA) and std::vector Particles (No NUMA)
#ifdef PARTICLENONUMA

    // Particle default constructor for No NUMA Particles
    Particle() {
        Nt = 1;
        xx = std::vector<float>(1);
        xy = std::vector<float>(1);
        vx = std::vector<float>(1);
        vy = std::vector<float>(1);
        Ex = std::vector<float>(1);
        Ey = std::vector<float>(1);
    }

    // Particle input-dependent constructor for No NUMA Particles
    Particle(int Nt_in) {
        Nt = Nt_in;
        xx = std::vector<float>(Nt);
        xy = std::vector<float>(Nt);
        vx = std::vector<float>(Nt);
        vy = std::vector<float>(Nt);
        Ex = std::vector<float>(Nt);
        Ey = std::vector<float>(Nt);
    }

    // empty delete (no dynamic memory for Non-NUMA version of Particle)
    ~Particle() { }

    // empty clean() (no dynamic memory for Non-NUMA version of Particle)
    void clean() { }

    // vectors of position, velocity, and E-field associated with each Particle object
    std::vector<float> xx;  
    std::vector<float> xy;
    std::vector<float> vx;
    std::vector<float> vy;
    std::vector<float> Ex;
    std::vector<float> Ey;

#else
    // Particle default constructor for NUMA Particles 
    Particle() {
        Nt = 1;
        xx = NULL;
        xy = NULL;
        vx = NULL;
        vy = NULL;
        Ex = NULL;
        Ey = NULL;
    }

     // free the memory on the heap (called automatically when Particle object is deleted)
    ~Particle() {
        Nt = -1;
        xx = NULL;
        xy = NULL;
        vx = NULL;
        vy = NULL;
        Ex = NULL;
        Ey = NULL;
    }

    // Particle input-dependent constructor for NUMA Particles 
    Particle(int Nt_in) {
        Nt = Nt_in;
        xx = new float[Nt];
        xy = new float[Nt];
        vx = new float[Nt];
        vy = new float[Nt];
        Ex = new float[Nt];
        Ey = new float[Nt];
    }

    // free the dynamically allocated memory 
    void clean() {
        delete[] xx;
        delete[] xy;
        delete[] vx; 
        delete[] vy; 
        delete[] Ex;
        delete[] Ey;
    }

    // Particle operator for equal
    Particle& operator=(const Particle& rhs) {
        xx = rhs.xx;
        xy = rhs.xy;
        vx = rhs.vx;
        vy = rhs.vy;
        Ex = rhs.Ex;
        Ey = rhs.Ey;
        Nt = rhs.Nt;
        return *this;
    }

    // arrays of position, velocity, and E-field associated with each NUMA Particle object
    float* xx;
    float* xy;
    float* vx;
    float* vy;
    float* Ex;
    float* Ey;
#endif

};
/*-----------------------------------------------------------------------------------------------------*/

#endif