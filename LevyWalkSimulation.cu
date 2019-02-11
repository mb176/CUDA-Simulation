#include <iostream>
#include <cstdio>
#include <time.h>
#include <string>
#include <iomanip>
#include <fstream>
#include <cmath>
#include <sstream>
#include <vector>
#include <algorithm>
#include <numeric>

#include "vector_calculus.h"
#include "fitting.h"
#include "cHelper.h"
#include "d_LevyWalkGo.h"

#define cudaError(ans) { gpuAssert((ans), __FILE__, __LINE__); }
inline void gpuAssert(cudaError_t code, const char *file, int line, bool abort=true)
{
   if (code != cudaSuccess)
   {
      fprintf(stderr,"GPUassert: %s %s %d\n", cudaGetErrorString(code), file, line);
      if (abort) exit(code);
   }
}

using namespace std;
typedef vector<double> vec;



class LevyWalkSimulation{
public:
  //Parameters
  double nu, eta, gamma, t0, c; //Model Parameters
  int blocksize, maxNParticles, nParticles;  //Simulation Parameters, maxNParticles defines how big an essemble can get
  std::vector<double> times, agingTimes; //Measurement times

  //Constructers;
  LevyWalkSimulation(){create();};
  LevyWalkSimulation(double tMax, double tMin, int nTimes, double taMax, double taMin, int nAgingTimes){
    create();
    calculateMeasurementTimes( tMax,  tMin,  nTimes,  taMax,  taMin,  nAgingTimes);
  };
  void calculateMeasurementTimes(double tMax, double tMin, int nTimes, double taMax, double taMin, int nAgingTimes);
  void create();

  //Main routine
  void LevyWalkGo();

  //Results
  std::vector<double> MSD, MSDaging;
  std::vector<double> MSDFitParameters; //Set by fitMSD, [slope, offset, error]
  std::vector<double> MSDAgingFitParameters; //Set by fitMSDaging, [slope, offset, error]
  void clearResults();

  //Fitting
  double fitMSD(const int skipValues);
  double fitMSDaging(const int skipValues);

  //Output
  std::ofstream safeResult(std::string path, std::string filename, std::string type);
    // type is "MSD" or "MSDaging", depending on what you want to safe
    // for filename = "auto" the name gets generated automatically


private:
    int uninitialisedValue = -1; //To check if user has set the parameters
    bool parameterTestSuccessful(void); //To test parameters
    std::vector<double> createParameterVector(std::string type); //To communicate with plotScript.py via safeResult

};

void LevyWalkSimulation::create(){
    nu = uninitialisedValue;
    eta = uninitialisedValue;
    gamma = uninitialisedValue;
    t0 = uninitialisedValue;
    c = uninitialisedValue;
    nParticles = uninitialisedValue;
    blocksize = 256;
    maxNParticles = 1000000;
}

void LevyWalkSimulation::calculateMeasurementTimes(double tMax, double tMin, int nTimes, double taMax, double taMin, int nAgingTimes){
    if(tMax > (taMax-taMin)/nAgingTimes && nAgingTimes>1){
        cerr << "The intervalls between aging times are smaller then the measurement duration!" << endl;
        return;
    }
    times.resize(nTimes);
    range(times, tMin, tMax);
    agingTimes.resize(nAgingTimes);
    range(agingTimes, taMin , taMax);
}

bool LevyWalkSimulation::parameterTestSuccessful(void){
    //Check Parameters
    if(nu <= 0 || eta<= 0 || gamma <= 0||t0 <= 0 || c <= 0 ){
        cerr << "Error: Model parameters not properly initialised" << endl;
        return false;
    } else if(nParticles <= 0) {
        cerr << "Error: Number of particles not set/ <= 0" << endl;
        return false;
    }
    //Check Times and agingTimes vectors
    //Not empty?
    if(times.size() == 0 ){
        cerr << "Error: times can't be empty" << endl;
        return false;
    } else if (agingTimes.size() == 0 ){
        cerr << "Error: agingTimes can't be empty" << endl;
        return false;
    }
    //Sorted?
    vec test = times; sort(test.begin(),test.end());
    if(test!=times){
        cerr << "Error: Times not sorted" << endl;
        return false;
    }
    test = agingTimes; sort(test.begin(),test.end());
    if(test!=agingTimes){
        cerr << "Error: agingTimes not sorted" << endl;
        return false;
    }
    //agingTimes steps bigger than largest element of times?
    if(agingTimes.size() > 1){
    adjacent_difference(agingTimes.begin(), agingTimes.end(), test.begin());
        if(*min_element(test.begin()+1,test.end()) < times.back()){
            cerr << "The intervalls between aging times are smaller then the measurement duration!" << endl;
            return false;
        }
    }

    //All good
    return true;
}

void LevyWalkSimulation::LevyWalkGo(){
    //Barricades
    if(parameterTestSuccessful()!=true){
        return;
    }

    //Constants
    const int nMeasurements = times.size()*agingTimes.size();

    //Vectors to Save MSD
    vec subtotalSD(nMeasurements,0);
    vec totalMSD(nMeasurements,0);

    //Split the Essemble into parts that the memory can handle
    int nEssembles = (nParticles-1)/maxNParticles + 1;

    //Loop over all Essembles
    {
    int essembleSize;
    int nBlocks;
    int nEntries;

    for(int essembleIdx = 0; essembleIdx != nEssembles; essembleIdx++){
        //Set Essemble size
        if(essembleIdx == nEssembles-1){
            essembleSize = nParticles - (nEssembles-1)*maxNParticles;
        } else {
            essembleSize = maxNParticles;
        }

        //Local Constants
        nBlocks = (essembleSize-1)/blocksize+1;
        nEntries = essembleSize * nMeasurements;

        //copy agingTimes and measurementTimes to device:
        double * d_agingTimes = vectorToDevice(&agingTimes[0], agingTimes.size());
        double * d_times = vectorToDevice(&times[0], times.size());

        //Create squared displacement (SD) Vectors
        vec SD(nEntries,0);
        double* d_SD = vectorToDevice(&SD[0],nEntries);

        //Let essembleSize Walkers walk and record their SD
        d_LevyWalkGo<<<(essembleSize-1)/blocksize+1, blocksize>>>
            (d_SD, nEntries, essembleSize, time(NULL), gamma, nu, eta, t0, c, d_times, times.size(),
            d_agingTimes, agingTimes.size()) ;



        //Retrieve SD from device
        cudaError(cudaDeviceSynchronize());
        cudaError(cudaMemcpy(&SD[0], d_SD, nEntries*sizeof(double) , cudaMemcpyDeviceToHost));

        //Sum over the current Enssemble and store it in subtotalSD
        {
        double  sum;
        int measurementIdx;
        for(int taIdx = 0; taIdx != agingTimes.size(); taIdx++){
          for(int timeIdx = 0; timeIdx != times.size(); timeIdx++){
            measurementIdx = taIdx * times.size() + timeIdx;
            sum = 0;
            for(int blockIdx=0; blockIdx!=nBlocks; blockIdx++){
              sum+=SD[measurementIdx*essembleSize+blockIdx*blocksize];
            }
          subtotalSD[measurementIdx]=sum;
          }
        }
        }

        //Free memory
        cudaFree(d_agingTimes);
        cudaFree(d_times);
        cudaFree(d_SD);

        //Sum subtotalSD to totalMSD
        add(totalMSD, subtotalSD, totalMSD);

    }//End loop over essembles
    }

    //Average over nParticles
    mult(totalMSD, 1.0/nParticles);

    //Save the t dependence of MSD for longest ta:
    MSD.clear();
    for(int idx = 0; idx != times.size(); idx++){
        MSD.push_back(totalMSD[ times.size()*(agingTimes.size()-1) + idx ] );
    }

    //Save ta dependence of MSD for longest t:
    MSDaging.clear();
    for(int idx = 0; idx != agingTimes.size(); idx++){
        MSDaging.push_back(totalMSD[ idx * times.size()+times.size()-1 ] );
    }
    }

double LevyWalkSimulation::fitMSD(const int skipValues){
    if(MSD.size()!= times.size() ){
        throw domain_error("Can't fit; MSDaging not yet calculated");
    }
    if(*min_element(times.begin(),times.end())<=0){
        throw domain_error("Can't fit; Nonpositive measurement times");
    } else if (*min_element(MSD.begin(),MSD.end())<=0){
        throw domain_error("Can't fit; MSD has nonpositive values");
    }
    vec x,y;
    for(int idx = skipValues; idx!=times.size(); idx++){
      x.push_back(times[idx]);
      y.push_back(MSD[idx]);
    }
    MSDFitParameters = exponent_fit(x,y);
    return MSDFitParameters[0];
}

double LevyWalkSimulation::fitMSDaging(const int skipValues){
    if(MSDaging.size()!= agingTimes.size() ){
        throw domain_error("MSDaging not yet calculated");
    }
    if(*min_element(agingTimes.begin(),agingTimes.end())<=0){
        throw domain_error("Can't fit; Nonpositive aging times");
    } else if (*min_element(MSDaging.begin(),MSDaging.end())<=0){
        throw domain_error("Can't fit; MSDaging has nonpositive values");
    }
    vec x,y;
    for(int idx = skipValues; idx!=times.size(); idx++){
      x.push_back(agingTimes[idx]);
      y.push_back(MSDaging[idx]);
    }
    MSDAgingFitParameters = exponent_fit(x, y);
    return MSDAgingFitParameters[0];
}

void LevyWalkSimulation::clearResults(){
  MSD.clear();
  MSDaging.clear();
  MSDFitParameters.clear(); //Set by fitMSD, [slope, offset, error]
  MSDAgingFitParameters.clear(); //Set by fitMSDaging, [slope, offset, error]
}

std::ofstream LevyWalkSimulation::safeResult(string path, string filename, string type)
{// type is MSD or MSD aging, depending on what you want to safe
 // for filename = "auto" the name gets generated automatically

    // Write data matrix [time,MSD, parameters]
    std::vector<std::vector<double> > data;
    if(type == "MSD"){
        if(MSD.size()!= times.size() ){
            throw domain_error("Can't Plot: MSD not yet calculated");
        }
        if(MSDFitParameters.size()== 0 ){
            throw domain_error("MSD not yet fitted");
        }
        //Write times and MSD into Data
        data.push_back(times);
        data.push_back(MSD);
        data.push_back(createParameterVector("MSD"));

    } else if (type == "MSDaging"){
        if(MSDaging.size()!= agingTimes.size() ){
            throw domain_error("Can't Plot: MSDaging not yet calculated");
        }
        if(MSDAgingFitParameters.size()== 0 ){
            throw domain_error("MSDaging not yet fitted");
        }
        //Write times and MSD into Data
        data.push_back(agingTimes);
        data.push_back(MSDaging);
        data.push_back(createParameterVector("MSDaging"));
    }

    //Write Path
    string completePath;
    if(filename == "auto"){
        std::ostringstream name;
        name << std::setprecision(2) << type << "_gamma_"<<gamma<<"_nu_"<<nu<<"_ta_"<<agingTimes.back()<<".txt";//<<"_exp_"<<MSDFitParameters[0]
        completePath = path + name.str();
    } else {
        completePath = path + filename;
    }

    return write(completePath, data);
}

std::vector<double> LevyWalkSimulation::createParameterVector(string type){
    //Write Parameters into data
    std::vector<double> parameters;
    if(type == "MSD"){
        parameters.resize(MSD.size());
        if(MSD.size() <6){
            throw domain_error("Can't save: MSD has less then 6 entries");
        }
        parameters[0] = 1; //1 For MSD, 2 for MSD aging
        parameters[1] = MSDFitParameters[0];
        parameters[2] = MSDFitParameters[1];
    } else if (type == "MSDaging"){
        if(MSDaging.size() <6){
            throw domain_error("Can't save: MSDaging has less then 6 entries");
        }
        parameters.resize(MSDaging.size());
        parameters[0] = 2; //1 For MSD, 2 for MSD aging
        parameters[1] = MSDAgingFitParameters[0];
        parameters[2] = MSDAgingFitParameters[1];
    }
    parameters[3] = gamma;
    parameters[4] = nu;
    parameters[5] = eta;
    return parameters;
}
