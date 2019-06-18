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
  int blocksize; // #number of threads per block in call of d_LevyWalkGo
  int maxNParticles; //defines how big an essemble can get before it is split up
  int nParticles; //total number of particles to be simulated
  std::vector<double> times, agingTimes; //For every agingTime the position will be measured for every time;

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
  // General Idea: First we split all the particles into essembles that the
  // device memory can handle; We then hand a vector SD to the device, where every thread has
  // spots to write down its squared displacement for every measurement time (for every
  // agingTime for every time). Every thread simulates one particle and writes down
  // its SD. Each block performs a reduction so that the sum
  // of the MSDs is in its first entry. On the host we add up these entries for
  // every essemble and then over all ensemble, divide by nParticles and obtain
  // the MSD.
  //The threads also indicate their position in the bin vector for Histograms

  // MSD
  std::vector<double> MSD, MSDaging;
  std::vector<double> MSDFitParameters; //Set by fitMSD, [slope, offset, error]
  std::vector<double> MSDAgingFitParameters; //Set by fitMSDaging, [slope, offset, error]
  void clearResults();

  // Histogram
  uint nBins;
  double histogramRange;
  double maximalDistance();
  std::vector<std::vector<int>>  histograms; //histograms[measurementTime][bin]

  //Fitting
  double fitMSD(const int skipValues);
  double fitMSDaging(const int skipValues);

  //Output
  std::ofstream safeResult(std::string path, std::string filename, std::string type);
    // type is "MSD" or "MSDaging", depending on what you want to safe
    // for filename = "auto" the name gets generated automatically
  std::ofstream safeHistogram(int agingTimeIdx, std::string path, std::string filename);


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

double analyticPredictionOrdinary(double nu, double eta, double gamma){
  if(nu<=0 || eta <= 0 || gamma <= 0){
    throw domain_error("Can't give Prediction: Parameters are not positive.");
  }

  if(gamma <= 2*(nu-eta)) {
    cout << "MSD does not converge for these parameters" << endl;
  }
  if(gamma < 1){
    if(2*nu < gamma){
      return gamma;
    } else if ( 2*nu >= gamma){
      return 2*nu;
    }
  } else if (gamma >= 1){
    if(2*nu < gamma){
      return 1;
    } else if ( 2*nu >= gamma){
      return 2*nu-gamma+1;
    }
  }
  return -1; //To check if all cases are covered
}

void LevyWalkSimulation::initialiseHistogram(uint numberOfBins, uint histogramrange){
  if(times.size()==0){
    cerr << "Can't initialse histogram: times vector not set." << endl;
    return;
  }
  if(agingTimes.size()==0){
    cerr << "Can't initialse histogram: agingTimes vector not set." << endl;
    return;
  }
  //Delete previous histograms
  histograms = {};
  nBins = numberOfBins;
  histogramRange = histogramrange;
  //Construct histogram matrix:
  vector<int> histogram(nBins,0);
  histograms.resize(times.size()*agingTimes.size(),histogram);
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
    //Histograms initialised?
    if(histograms.size()==0 || histograms[0].size()==0){
      cerr << "Error: Histograms was not initialised" << endl;
      return false;
    }

    //All good
    return true;
}

void LevyWalkSimulation::LevyWalkGo(){


    //Barricades
    if(parameterTestSuccessful()!=true){
        cerr << "Can't start Levy Walk: Parameter test unsuccessful." << endl;
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
        //Set current essemble size
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

        //Create squared displacement (SD) and bins vectors
        vec SD(nEntries,1);
        double* d_SD = vectorToDevice(&SD[0],nEntries);
        vector<int> bins(nEntries,0);
        int * d_bins;
        cudaError(cudaMalloc((void**)&d_bins, nEntries*sizeof(int)) );
        cudaError(cudaMemcpy(d_bins, &bins[0], nEntries*sizeof(int) ,cudaMemcpyHostToDevice));

        //Let essembleSize Walkers walk and record their SD
        d_LevyWalkGo<<<(essembleSize-1)/blocksize+1, blocksize>>>
            (d_SD, nEntries, essembleSize, time(NULL), gamma, nu, eta, t0, c, d_times, times.size(),
            d_agingTimes, agingTimes.size(), d_bins, nBins, histogramRange) ;

        //Retrieve SD and bins from device
        cudaError(cudaDeviceSynchronize());
        cudaError(cudaMemcpy(&SD[0], d_SD, nEntries*sizeof(double) , cudaMemcpyDeviceToHost));
        cudaError(cudaMemcpy(&bins[0], d_bins, nEntries*sizeof(int) , cudaMemcpyDeviceToHost));



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

        //Sum subtotalSD to totalMSD
        add(totalMSD, subtotalSD, totalMSD);

        //Insert bins into the histograms
        {
        int binIdx, measurementIdx;
        for(int taIdx = 0; taIdx!=agingTimes.size(); taIdx++){
          for(int timeIdx = 0; timeIdx!=times.size(); timeIdx++){
            measurementIdx = taIdx * times.size() + timeIdx;
            for(int particleIdx = 0; particleIdx!=essembleSize; particleIdx++){
              binIdx = bins[measurementIdx*essembleSize+particleIdx];
              histograms[measurementIdx][binIdx]++;
            }
          }
        }
        }

        //Free memory
        cudaFree(d_agingTimes);
        cudaFree(d_times);
        cudaFree(d_SD);
        cudaFree(d_bins);



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
    if(*min_element(times.begin()+skipValues,times.end())<=0){
        throw domain_error("Can't fit; Nonpositive measurement times");
    } else if (*min_element(MSD.begin()+skipValues,MSD.end())<=0){
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
        throw domain_error("Can't fit; MSDaging not yet calculated");
    }
    if(*min_element(agingTimes.begin()+skipValues,agingTimes.end())<=0){
        throw domain_error("Can't fit; Nonpositive aging times");
    } else if (*min_element(MSDaging.begin()+skipValues,MSDaging.end())<=0){
        throw domain_error("Can't fit; MSDaging has nonpositive values");
    }
    vec x,y;
    for(int idx = skipValues; idx!=agingTimes.size(); idx++){
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
{// type is "MSD" or "MSDaging", depending on what you want to safe
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

std::ofstream LevyWalkSimulation::safeHistogram(int agingTimeIdx, std::string path, std::string filename){
  if(agingTimeIdx>=agingTimes.size()){
    throw domain_error("Can't save histogram: agingTimeIdx out of bounds.");
  }
  //Create complete path
  string completePath;
  if(filename == "auto"){
      std::ostringstream name;
      name << std::setprecision(2) << "histogram" << "_gamma_"<<gamma<<"_nu_"<<nu<< "_eta_"<<eta<<"_ta_"<<agingTimes.back()<<".txt";//<<"_exp_"<<MSDFitParameters[0]
      completePath = path + name.str();
  } else {
      completePath = path + filename;
  }

  ofstream outfile(completePath);

  //Cut off empty bins
  int cutOff = 0;
  for(int timeIdx=0; timeIdx!=times.size(); timeIdx++){
    for(int idx = nBins-1; idx!= 0; idx--){
      if (histograms[agingTimeIdx*times.size()+timeIdx][idx]!=0){
        if(cutOff<idx){ cutOff = idx;}
        break;
      }
    }
  }

  //Write x values (use middle of each bin)
  for(int i = 0; i!= cutOff+1; i++){
    outfile << i*histogramRange/nBins+histogramRange/nBins/2 << " ";
  }
  outfile << endl;

  //Write normalised histogramm for each time out of times:
  for(int timeIdx=0; timeIdx!=times.size(); timeIdx++){
    for(int i = 0; i!= cutOff+1; i++){
      outfile << double(histograms[agingTimeIdx*times.size()+timeIdx][i])/nParticles<< " ";
    }

    outfile << endl;
  }
  return outfile;
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

double LevyWalkSimulation::maximalDistance(){
  //runs levyWalkGo() for fewer particles and measures the maximum Distance totalMax;
  int nParticlesBAK = nParticles;
  nParticles = 5000;
  int blocksizeBAK = blocksize;
  blocksize =1; //To avoid reduction on GPU
  double totalMax=0;

  //Barricades
  if(parameterTestSuccessful()!=true){
      cerr << "Can't start Levy Walk: Parameter test unsuccessful." << endl;
      return 0;
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
      //Set current essemble size
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

      //Create squared displacement (SD) and bins vectors
      vec SD(nEntries,1);
      double* d_SD = vectorToDevice(&SD[0],nEntries);
      vector<int> bins(nEntries,0);
      int * d_bins;
      cudaError(cudaMalloc((void**)&d_bins, nEntries*sizeof(int)) );
      cudaError(cudaMemcpy(d_bins, &bins[0], nEntries*sizeof(int) ,cudaMemcpyHostToDevice));

      //Let essembleSize Walkers walk and record their SD
      d_LevyWalkGo<<<(essembleSize-1)/blocksize+1, blocksize>>>
          (d_SD, nEntries, essembleSize, time(NULL), gamma, nu, eta, t0, c, d_times, times.size(),
          d_agingTimes, agingTimes.size(), d_bins, nBins, histogramRange) ;

      //Retrieve SD and bins from device
      cudaError(cudaDeviceSynchronize());
      cudaError(cudaMemcpy(&SD[0], d_SD, nEntries*sizeof(double) , cudaMemcpyDeviceToHost));
      cudaError(cudaMemcpy(&bins[0], d_bins, nEntries*sizeof(int) , cudaMemcpyDeviceToHost));



      //Find biggest squared displacement in the essemble and store it in subtotalSD
      {
      double  essembleMax=0;
      int measurementIdx;
      for(int taIdx = 0; taIdx != agingTimes.size(); taIdx++){
        for(int timeIdx = 0; timeIdx != times.size(); timeIdx++){
          measurementIdx = taIdx * times.size() + timeIdx;
          for(int blockIdx=0; blockIdx!=nBlocks; blockIdx++){
            essembleMax= max(essembleMax,SD[measurementIdx*essembleSize+blockIdx*blocksize]);
          }
          subtotalSD[measurementIdx]=essembleMax;
        }
      }
      }

      //Find biggest SD for all ensembles at final time and agingTime
      uint finalIdx = agingTimes.size()*times.size()-1;
      totalMax = max(subtotalSD[finalIdx],totalMax);

      //Free memory
      cudaFree(d_agingTimes);
      cudaFree(d_times);
      cudaFree(d_SD);
      cudaFree(d_bins);
  }//End loop over essembles
  }

  //Reset all the changes in LevyWalkSimulation
  nParticles = nParticlesBAK;
  blocksize = blocksizeBAK;
  clearResults();

  totalMax = pow(totalMax,0.5);
  return totalMax;
}
