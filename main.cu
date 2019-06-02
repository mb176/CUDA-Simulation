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
#include "fitting.h"





using namespace std;
typedef vector<double> vec;



int main(void)
{

  /* %%%%%%%%%%%%%%%% Configure Levy Walk Simulation %%%%%%%%%%%%%%%%%%%%%%%% */
  //Aging and observation times
  const int tMax(20000), tMin(0), taMax(0),taMin(0), nTimes(10), nAgingTimes(1) ; //needs at least one aging and one measurement time
  //Histogram parameters
  uint nBins = 10000;
  double histogramRange = pow(10,11);
  //Initialise Simulation:
  LevyWalkSimulation LW(tMax, tMin, nTimes, taMax, taMin, nAgingTimes, nBins, histogramRange);
  //Model Parameters
  LW.t0 = 1; //Characteristic timescale
  LW.c = 1; //Constant Velocity prefactor
  //Simulation parameters
  LW.nParticles    = pow(10,9); //Total size of the Essemble
  LW.maxNParticles = 2000; //memory Breaks down if maxNparticles*nTimes*nAgingTimes = 10^9
  LW.blocksize = 256; //For the kernel call, must bu multiple of 32;
  string targetFolder = "Results/Histogram_06_02/";
  /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */

  //Varying Parameters
  // LW.gamma = 1.5;
  // LW.nu = 1.5;
  // LW.eta = 1.5;
  vec nu= {0.9,1.2,1.5};
  vec gamma = {0.9};

  {
  ofstream output(targetFolder+"MSDTable.txt");
  output << "Simulation mit " << LW.nParticles << " Teilchen; tMax = " << tMax << "; taMax = " << taMax << endl;
  for(int nuIdx = 0; nuIdx != nu.size(); nuIdx++ ){
    for(int gammaIdx = 0; gammaIdx != gamma.size(); gammaIdx++ ){
  //
      //Set parameters Dynamically
      LW.nu = nu[nuIdx];
      LW.eta = nu[nuIdx];
      LW.gamma = gamma[gammaIdx];

      // Set histogram parameters Dynamically
      histogramRange = pow(pow(10,3),analyticPredictionOrdinary(LW.nu, LW.eta, LW.gamma));
      LW.initialiseHistogram(nBins, histogramRange);

      //Start Simulation
      LW.LevyWalkGo();
      // LW.fitMSD(1); //gives number of values to skip
      LW.safeHistogram(0,targetFolder,"auto");
      // LW.fitMSDaging(1); //gives number of values to skip
      // LW.safeResult(targetFolder,"auto","MSD");
      // LW.safeResult(targetFolder,"auto","MSDaging");
      // output << "Gamma: " << LW.gamma << "; nu: " << LW.nu << "; t: "
      // << LW.MSDFitParameters[0]  << "; ta: " << LW.MSDAgingFitParameters[0] << endl;
      cout << "Walk done!" << endl;

    } //end nu loop
  } //end gamma loop
  }



  return 1;
}
