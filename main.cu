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
  vec gamma = {1.5};
  vec nu= {0.4,0.8,1.5};
  vec eta = {-0.2,0,0.2};

  //Output Folder and metaData filename
  string targetFolder = "Results/Simulation_06_26/";
  string MetaFileName = "MetaData2.txt";
  ofstream metaFile(targetFolder+MetaFileName);

  for(int nuIdx = 0; nuIdx != nu.size(); nuIdx++ ){
    for(int gammaIdx = 0; gammaIdx != gamma.size(); gammaIdx++ ){
      for(int etaIdx = 0; etaIdx != eta.size(); etaIdx++){
/* %%%%%%%%%%%%%%%% Configure simulation parameters %%%%%%%%%%%%%%%%%%%%%%%% */
        //Aging and observation times
        const int tMax(1000), tMin(0), taMax(100000),taMin(0), nTimes(20), nAgingTimes(20) ; //needs at least one aging and one measurement time
        //Initialise Simulation:
        LevyWalkSimulation LW(tMax, tMin, nTimes, taMax, taMin, nAgingTimes);
        //Model Parameters
        LW.t0 = 1; //Characteristic timescale
        LW.c = 1; //Constant Velocity prefactor
        LW.nu = nu[nuIdx];
        LW.eta = nu[nuIdx]+eta[etaIdx];
        LW.gamma = gamma[gammaIdx];
        //Simulation parameters
        LW.nParticles    = pow(10,5); //Total size of the Essemble
        LW.maxNParticles = 2048; //GPU memory breaks down if maxNparticles*nTimes*nAgingTimes = 3*10^9
        LW.blocksize = 256; //For the kernel call, must be multiple of 32;
        //Histogram parameters
        uint nBins = 1;
        double maxDistance =3;//LW.maximalDistance(); //How for do 5000 particles travel with this setup?
        double histogramRange = 1.5*maxDistance; //Estimate HistgramRange (needs nu, eta, gamma, etc.)
        //Initialise Histogram
        LW.initialiseHistogram(nBins, histogramRange);
  /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */

        //Write parameters into file
        if(nuIdx == 0 && gammaIdx == 0 && etaIdx==0){
          metaFile << "Simulation mit  "<< LW.nParticles << " Teilchen; tMax = " << tMax << "; taMax = " << taMax << "; nBins = " << nBins<< endl;
        }
        //Start Simulation
        auto start = chrono::steady_clock::now();
        LW.LevyWalkGo();
        auto end = chrono::steady_clock::now();
        cout << "Simulation duration: " << chrono::duration_cast<chrono::seconds>(end-start).count() << " s" <<endl;

        // Histogram
        // LW.safeHistogram(0,targetFolder,"auto")

        //Find analytical predictions for the t and ta dependence
        double tExponentPrediction = 0;
        double taExponentPrediction = 0;
        if(taMax > tMax){
          std::vector<double> vecTemp = LW.analyticPredictionAged();
          tExponentPrediction = vecTemp[0];
      # #Run executable directly
# ./$ex# #Run executable directly
# ./$executableName
# # Create detached screen session on computer and run it there
# screen -d -m -S simulation ./$executableNameecutableName
# # Create detached screen session on computer and run it there
# screen -d -m -S simulation ./$executableName    taExponentPrediction = vecTemp[1];
        } else {
          tExponentPrediction = LW.analyticPredictionOrdinary();
        }

        //Time dependence of the MSD (needs ntimes>1)
        LW.fitMSD(3); //gives number of values to skip
        LW.safeResult(targetFolder,"auto","MSD");
        metaFile << std::setprecision(3) << "Gamma: " << LW.gamma << "; nu: "
                 << LW.nu << "; eta: " << LW.eta << ": t-exp: "<< LW.MSDFitParameters[0]
                 << " (" << tExponentPrediction << ",Chi^2=" << LW.MSDFitParameters[2] <<") ";

        //Aging time dependence of the MSD (needs nAgingTimestimes>1 and taMax/nAgingTimes>tMax)
        LW.fitMSDaging(3); //gives number of values to skip
        LW.safeResult(targetFolder,"auto","MSDaging");
        metaFile << "; ta-exp: " << LW.MSDAgingFitParameters[0]
                 << "( "<< taExponentPrediction << ",chi=" << LW.MSDAgingFitParameters[2] <<")";


        metaFile << endl;
      } // end eta loop
    } //end nu loop
  } //end gamma loop

  return 1;
}
