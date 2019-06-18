#include <iostream>
#include <cuda.h>
#include <curand_kernel.h>
#include <vector>
#include <time.h>
#include <cmath>
#include <cstdio>
#include <sstream>

#include "LevyWalkSimulation.h"
#include "d_LevyWalkGo.h"
#include "cHelper.h"
#include "vector_calculus.h"
#include "fitting.h"





using namespace std;
typedef vector<double> vec;



int main(void)
{

  /* %%%%%%%%%%%%%%%% Configure Static parameters %%%%%%%%%%%%%%%%%%%%%%%% */
  //Aging and observation times
  const int tMax(20000), tMin(0), taMax(0),taMin(0), nTimes(10), nAgingTimes(1) ; //needs at least one aging and one measurement time
  //Histogram parameters
  uint nBins = 10000;
  double histogramRange = pow(10,8); //doesn't matter, reset before start of the simulation
  //Initialise Simulation:
  LevyWalkSimulation LW(tMax, tMin, nTimes, taMax, taMin, nAgingTimes, nBins, histogramRange);
  //Model Parameters
  LW.t0 = 1; //Characteristic timescale
  LW.c = 1; //Constant Velocity prefactor
  //Simulation parameters
  LW.nParticles    = pow(10,9); //Total size of the Essemble
  LW.maxNParticles = 2000; //memory Breaks down if maxNparticles*nTimes*nAgingTimes = 10^9
  LW.blocksize = 256; //For the kernel call, must be multiple of 32;
  string targetFolder = "Results/Histogram_06_12/";
  /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */

  //Varying Parameters

  vec gamma = {1.3};
  vec nu= {1.3};
  vec eta = {1.1,1.3,1.5};

  {
  stringstream MetaFileName;
  MetaFileName << "Metadata_"<< "gamma_" <<LW.gamma << "_nu_" <<LW.nu << "_eta_" << LW.eta;
  ofstream output(targetFolder+MetaFileName.str());
  output << "Simulation mit  "<< LW.nParticles << " Teilchen; tMax = " << tMax << "; taMax = " << taMax << endl;
  for(int nuIdx = 0; nuIdx != nu.size(); nuIdx++ ){
    for(int gammaIdx = 0; gammaIdx != gamma.size(); gammaIdx++ ){
      for(int etaIdx = 0; etaIdx != eta.size(); etaIdx++){

        //Set parameters Dynamically
        LW.nu = nu[nuIdx];
        LW.eta = eta[etaIdx];
        LW.gamma = gamma[gammaIdx];



        // Find histogramRange appropriate for theses Parameters
        // histogramRange = LW.c/2 * pow(10000,analyticPredictionOrdinary(LW.nu, LW.eta, LW.gamma));
        double maxDistance = LW.maximalDistance(); //How for do 5000 particles get with this setup?

        histogramRange = 1.5*maxDistance;
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
      } // end eta loop
    } //end nu loop
  } //end gamma loop
  }



  return 1;
}
