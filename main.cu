#include <iostream>
#include <cuda.h>
#include <curand_kernel.h>
#include <vector>
#include <time.h>
#include <cmath>
#include <cstdio>

#include "LevyWalkSimulation.h"
#include "d_LevyWalkGo.h"
#include "cHelper.h"
#include "vector_calculus.h"




using namespace std;
typedef vector<double> vec;



int main(void)
{

  /* %%%%%%%%%%%%%%%% Configure Levy Walk Simulation %%%%%%%%%%%%%%%%%%%%%%%% */
  LevyWalkSimulation LW;
  //Aging and observation times
  const int tMax(1000), tMin(1), taMax(0),taMin(0), nTimes(10000), nAgingTimes(1) ;
  LW.calculateMeasurementTimes(tMax, tMin, nTimes, taMax, taMin, nAgingTimes);
  //Model Parameters
  LW.t0 = 1;
  LW.c = 1;
  LW.gamma = 0.98;
  LW.nu = 1.5;
  LW.eta = 1;
  //Simulation parameters
  LW.nParticles    = 100000; //Total size of the Essemble
  LW.maxNParticles = 10000; //memory Breaks down if maxNparticles*nTimes*nAgingTimes = 10^9
  LW.blocksize = 256; //For the kernel call, must bu multiple of 32;
  /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */

  //Start Simulation
  LW.LevyWalkGo();
  LW.fitMSD(30);
  LW.safeResult("Results/","auto","MSD");
  cout << LW.MSDFitParameters[0] << endl;

  return 1;
}
