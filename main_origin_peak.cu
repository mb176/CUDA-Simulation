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
  vec gamma = {0.8};
  vec nu= {0.8,0.5,1.4};
  vec eta = {0}; //this is actually the difference to nu, i.e. eta-nu

  //Output Folder and metaData filename
  string targetFolder = "Results/Histogram_origin_peak/";
  ofstream metaFile;

  for(int nuIdx = 0; nuIdx != nu.size(); nuIdx++ ){
    for(int gammaIdx = 0; gammaIdx != gamma.size(); gammaIdx++ ){
      for(int etaIdx = 0; etaIdx != eta.size(); etaIdx++){
/* %%%%%%%%%%%%%%%% Configure simulation parameters %%%%%%%%%%%%%%%%%%%%%%%% */
        //Aging and observation times
        const int tMax(20000), tMin(0), taMax(0),taMin(0), nTimes(20), nAgingTimes(1) ; //needs at least one aging and one measurement time
        //Initialise Simulation:
        LevyWalkSimulation LW(tMax, tMin, nTimes, taMax, taMin, nAgingTimes);
        //Model Parameters
        LW.t0 = 1; //Characteristic timescale
        LW.c = 1; //Constant Velocity prefactor
        LW.nu = nu[nuIdx];
        LW.eta = nu[nuIdx]+eta[etaIdx];
        LW.gamma = gamma[gammaIdx];
        //Simulation parameters
        LW.nParticles    = pow(10,6); //Total size of the Essemble
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
        MetaFileName << "origin_peak" << "_gamma_"<<LW.gamma<<"_nu_"<<LW.nu<< "_eta_"<<LW.eta<<".txt";
        ofstream metaFile(targetFolder+MetaFileName.str());

        //Write parameters into file
          metaFile << "Simulation mit  "<< LW.nParticles << " Teilchen; tMax = " << tMax << "; taMax = " << taMax << "; nBins = " << nBins
                    << " ; gamma: " << LW.gamma << " ;nu: " << LW.nu << endl;


        //Start Simulation
        auto start = chrono::steady_clock::now();
        LW.LevyWalkGo();
        auto end = chrono::steady_clock::now();
        cout << "Simulation duration: " << chrono::duration_cast<chrono::seconds>(end-start).count() << " s" <<endl;

        // Histogram
        LW.safeHistogram(0,targetFolder,"auto");

        const uint skip = 1;
        // Write time steps
        vec time = {};
        for(int tIdx = skip; tIdx != LW.times.size(); tIdx++){
          time.push_back(LW.times[tIdx]);
          metaFile  << *(time.end()-1) << " ";
        }
        metaFile << " # time" << endl;

        // Write probability density at origin over time
        vec pOrigin= {};
        for(int tIdx = skip; tIdx != LW.times.size(); tIdx++){
          pOrigin.push_back(LW.histograms[tIdx][0]/double(LW.nParticles));
          metaFile << *(pOrigin.end()-1) << " ";
        }
        metaFile << " # probability density at origin" << endl;

        // Write probability density at peak over time
        vec pPeak = {};
        for(int tIdx = skip; tIdx != LW.times.size(); tIdx++){
          //Find Cut off for empty bins
          int cutOff = 0;
          for(int idx = LW.nBins-1; idx!= 0; idx--){
            if (LW.histograms[tIdx][idx]!=0){
              if(cutOff<idx){ cutOff = idx;}
              break;
            }
          }
          //Find Peak
          pPeak.push_back(*max_element(LW.histograms[tIdx].begin()+int(3*cutOff/4),LW.histograms[tIdx].begin()+int(cutOff+1))/double(LW.nParticles));
          metaFile << *(pPeak.end()-1) << " ";
        }
        metaFile << " # probability density at peak" <<endl;

        //Fitting
        vec fitPOrigin, fitPPeak;
        fitPOrigin = exponent_fit(time, pOrigin);
        fitPPeak = exponent_fit(time, pPeak);
        metaFile  << fitPOrigin[0] << " " << fitPOrigin[1] << " " << fitPPeak[0] << " " << fitPPeak[1] << " "
        << fitPOrigin[2] << " " << fitPPeak[2] << " # exponent origin, offset origin, exponent peak, offset peak, chi2 orign, chi2 peak" << endl;

        //Hand over some parameters
        metaFile  << LW.gamma << " " << LW.nu << " " << LW.eta << " # gamma, nu, eta" <<endl;

      } // end eta loop
    } //end nu loop
  } //end gamma loop

  return 1;
}
