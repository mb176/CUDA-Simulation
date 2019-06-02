#ifndef LevyWalkSimulation_GUARD
#define LevyWalkSimulation_GUARD

#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <vector>



class LevyWalkSimulation{
public:
  //Parameters
  double nu, eta, gamma, t0, c; //Model Parameters
  int blocksize, maxNParticles, nParticles;  //Simulation Parameters, maxNParticles defines how big an essemble can get
  std::vector<double> times, agingTimes; //Measurement times

  //Constructers;
  LevyWalkSimulation(){create();};
  LevyWalkSimulation(double tMax, double tMin, int nTimes, double taMax, double taMin, int nAgingTimes, uint nBins, double histogramRange){
    create();
    calculateMeasurementTimes( tMax,  tMin,  nTimes,  taMax,  taMin,  nAgingTimes);
    initialiseHistogram(nBins, histogramRange);
  };
  void calculateMeasurementTimes(double tMax, double tMin, int nTimes, double taMax, double taMin, int nAgingTimes);
  void initialiseHistogram(uint nBins, uint histogramRange);
  void create();



  //Main routine
  void LevyWalkGo();

  // MSD
  std::vector<double> MSD, MSDaging;
  std::vector<double> MSDFitParameters; //Set by fitMSD, [slope, offset, error]
  std::vector<double> MSDAgingFitParameters; //Set by fitMSDaging, [slope, offset, error]
  void clearResults();

  // Histogram
  uint nBins;
  double histogramRange;
  std::vector<std::vector<int>>  histograms;


  //Fitting
  double fitMSD(const int skipValues);
  double fitMSDaging(const int skipValues);

  //Output
  std::ofstream safeResult(std::string path, std::string filename, std::string type);
    // type is MSD or MSD aging, depending on what you want to safe
    // for filename = "auto" the name gets generated automatically
  std::ofstream safeHistogram(int agingTimeIdx, std::string path, std::string filename);

  private:
      int uninitialisedValue = -1; //To check if user has set the parameters
      bool parameterTestSuccessful(void); //To test parameters
      std::vector<double> createParameterVector(std::string type); //To communicate with plotScript.py via safeResult
};

double analyticPredictionOrdinary(double nu, double eta, double gamma);

#endif
