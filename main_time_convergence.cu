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
  vec gamma = {0.4,0.6,0.8};
  vec nu= {1.2};
  vec eta = {0};
  vec times = {2,5,7,10,15,20,30,40,50,75,100};


  //Output Folder and metaData filename
  string targetFolder = "Results/MSD_time_convergence/"; //don't forget the final "/"!
  ofstream metaFile;


  for(int nuIdx = 0; nuIdx != nu.size(); nuIdx++ ){
    for(int gammaIdx = 0; gammaIdx != gamma.size(); gammaIdx++ ){
      for(int etaIdx = 0; etaIdx != eta.size(); etaIdx++){
        //Save deviation from prediction to fit
        vec deviation={};

        // Create output file
        stringstream MetaFileName;
        MetaFileName << "convergence_"<< "gamma_" << gamma[gammaIdx] <<"_nu_" << nu[nuIdx] << "_eta_" << nu[nuIdx]+eta[etaIdx]<<".txt" ;
        ofstream metaFile(targetFolder+MetaFileName.str());

        for(int tIdx =0; tIdx != times.size(); tIdx ++){
/* %%%%%%%%%%%%%%%% Configure simulation parameters %%%%%%%%%%%%%%%%%%%%%%%% */
        //Aging and observation times
        const int tMax(times[tIdx]), tMin(0), taMax(0),taMin(0), nTimes(20), nAgingTimes(1) ; //needs at least one aging and one measurement time
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
        if(tIdx ==0){
          metaFile << "Simulation mit  "<< LW.nParticles << " Teilchen; tMax = " << tMax << "; taMax = " << taMax << "; nBins = " << nBins << endl;
        }
        //Start Simulation
        auto start = chrono::steady_clock::now();
        LW.LevyWalkGo();
        auto end = chrono::steady_clock::now();
        cout << "Simulation duration: " << chrono::duration_cast<chrono::seconds>(end-start).count() << " s" <<endl;

        //Find analytical predictions for the t and ta dependence
        double tExponentPrediction = 0;
        double taExponentPrediction = 0;
        if(taMax > tMax){
          std::vector<double> vecTemp = LW.analyticPredictionAged();
          tExponentPrediction = vecTemp[0];
          taExponentPrediction = vecTemp[1];
        } else {
          tExponentPrediction = LW.analyticPredictionOrdinary();
        }

        //Time dependence of the MSD (needs ntimes>1)
        LW.fitMSD(3); //gives number of values to skip
        deviation.push_back(abs(LW.MSDFitParameters[0]-tExponentPrediction));

        } // end times loop

        //Write results into metaFile
        for(int tIdx =0; tIdx != times.size(); tIdx++){ metaFile << times[tIdx] << " "; }
        metaFile << "# duration of observation"<< endl;

        for(int tIdx =0; tIdx != times.size(); tIdx++){ metaFile << deviation[tIdx] << " "; }
        metaFile << "# deviation of the fitted exponent from the true value"<< endl;

        //Fit deviation from prediction:
        vec parameters;
        parameters = exponent_fit(times,deviation);
        metaFile << parameters[0] << " "<< parameters[1] << " # Exponent and prefactor of the fit of the deviation" << endl;
      } // end eta loop
    } //end nu loop
  } //end gamma loop


  return 1;
}
