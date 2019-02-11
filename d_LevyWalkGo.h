#ifndef d_LevyWalkGO_GUARD
#define d_LevyWalkGO_GUARD

#include "curand_kernel.h"



__device__ void nextStep(const double oldTime, const double oldPosition, double& newTime,
    double& newPosition, double gamma,  double nu, double eta, double t0,
    double c, curandState_t &state);

    __device__ double measurePosition(double measurementTime, double oldTime,
        double oldPosition, double newTime, double newPosition,
        double nu, double eta,  double c);

__global__ void d_LevyWalkGo(double* d_MSD, int nMSD,
          int nParticles, int seed, double gamma,  double nu, double eta, double t0,
          double c, double* d_times, const int nTimes,
          const double* d_agingTimes, const int nAgingTimes);



#endif
