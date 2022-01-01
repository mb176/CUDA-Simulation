#include <iostream>
#include <cuda.h>
#include <curand_kernel.h>
#include <vector>
#include <time.h>
#include <cmath>
#include <cstdio>
#include <sstream>
#include <chrono>
#include <iomanip>
#include <algorithm>
#include <thread>

#include "LevyWalkSimulation.h"
#include "d_LevyWalkGo.h"
#include "cHelper.h"
#include "vector_calculus.h"
#include "fitting.h"

using namespace std;
typedef vector<double> vec;

int main(void)
{

  //Loop over parameters
  vec gamma = {1.1};
  vec nu= {1.3};
  vec eta = {0,0.2,-0.2};

  //Output Folder and metaData filename
  string targetFolder = "Results/Histogram_aged/";
  ofstream metaFile;

  for(int nuIdx = 0; nuIdx != nu.size(); nuIdx++ ){
    for(int gammaIdx = 0; gammaIdx != gamma.size(); gammaIdx++ ){
      for(int etaIdx = 0; etaIdx != eta.size(); etaIdx++){
/* %%%%%%%%%%%%%%%% Configure simulation parameters %%%%%%%%%%%%%%%%%%%%%%%% */
        //Aging and observation times
        const int tMax(1000), tMin(0), taMax(20000),taMin(0), nTimes(20), nAgingTimes(3); //needs at least one aging and one measurement time
        //Initialise Simulation:
        LevyWalkSimulation LW(tMax, tMin, nTimes, taMax, taMin, nAgingTimes);
        //Model Parameters
        LW.t0 = 1; //Characteristic timescale
        LW.c = 1; //Constant Velocity prefactor
        LW.nu = nu[nuIdx];
        LW.eta = nu[nuIdx]+eta[etaIdx];
        LW.gamma = gamma[gammaIdx];
        //Simulation parameters
        LW.nParticles = pow(10,9); //Total size of the Essemble
        LW.maxNParticles = 2048; //GPU memory breaks down if maxNparticles*nTimes*nAgingTimes = 3*10^9
        LW.blocksize = 256; //For the kernel call, must be multiple of 32;
        //Histogram parameters
        uint nBins = 1000;
        double histogramRange = 1; //Estimate HistgramRange (needs nu, eta, gamma, etc.)
        LW.initialiseHistogram(nBins, histogramRange); //just to test
        histogramRange = 1.5*LW.maximalDistance(); //How for do 5000 particles travel with this setup?
        LW.initialiseHistogram(nBins, histogramRange);
  /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
        //Create Meta file
        stringstream MetaFileName;
        MetaFileName << "metaData" << "_gamma_"<<LW.gamma<<"_nu_"<<LW.nu<< "_eta_"<<LW.eta<<".txt";
        ofstream metaFile(targetFolder+MetaFileName.str());

        //Write parameters into file
          metaFile << "Simulation mit  "<< LW.nParticles << " Teilchen; tMax = " << tMax << "; taMax = " << taMax << "; nBins = " << nBins
                    << " ; gamma: " << LW.gamma << " ;nu: " << LW.nu << endl;


        //Start Simulation
        auto start = chrono::steady_clock::now();
        LW.LevyWalkGo();
        auto end = chrono::steady_clock::now();
        cout << "Simulation duration: " << chrono::duration_cast<chrono::seconds>(end-start).count() << " s" <<endl;

        // Safe Histogram from last aging time (1 if no aging)
        LW.safeHistogram(LW.agingTimes.size()-1,targetFolder,"auto");


      } // end eta loop
    } //end nu loop
  } //end gamma loop

  return 1;
}
