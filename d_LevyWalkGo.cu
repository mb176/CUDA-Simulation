#include<cmath>
#include<iostream>
#include<vector>
#include<time.h>

#include "curand_kernel.h"
#include "cHelper.h"
#include "d_LevyWalkGo.h"

using namespace std;

__device__ void nextStep(const double oldTime, const double oldPosition, double& newTime,
    double& newPosition, double gamma,  double nu, double t0,
    double c, curandState_t &state){


    //Get length, duration and direction of the next step;
    double randomNumber = curand_uniform_double(&state);
    double duration = t0*(pow(randomNumber,-1/gamma)-1);
    double length = c*pow(duration,nu);
    randomNumber = curand_uniform_double(&state);
    int direction = (randomNumber+1.5)/2;
    direction = 2*(direction-0.5); //+1 or -1


    //Update positions
    newPosition = oldPosition + direction*length;
    newTime = oldTime + duration;
}

__device__ double measurePosition(double measurementTime, double oldTime,
    double oldPosition, double newTime, double newPosition,
    double nu, double eta,  double c){
  //Calculate Times
  double elapsedTime = measurementTime - oldTime;
  double totalTime = newTime - oldTime;

  //Calculate position at measurement time
  double length = c*pow(elapsedTime,eta)*pow(totalTime, nu-eta);
  double direction = newPosition-oldPosition;
  direction = direction/abs(direction);

  return oldPosition+direction*length;
}

__device__ int positionToBin(double position, uint nBins, double histogramRange){
  int binIdx;
  if(abs(position)>histogramRange){
    binIdx = nBins-1;
  } else {
    //Take absoulute value of position since it is symmetric
    binIdx = floor(abs(position)/(histogramRange/nBins));
  }
  return binIdx;
}

__global__ void d_LevyWalkGo(double* d_SD, int nSD,
  int nParticles, int seed, double gamma,  double nu, double eta, double t0,
  double c, double* d_times, const int nTimes, const double* d_agingTimes,
  const int nAgingTimes, int *d_bins, uint nBins, double histogramRange){
  /* Input explainations:
  -d_SD: is a vecor containing the squared displacements SDs for all Essemble elements and all measurement
  times, so [<SDs at t0>, <SDs at t1>,...]
  */

  //Find thread Index
  const int threadIndex = blockDim.x*blockIdx.x+threadIdx.x;

  //Is this thread part of the Essemble?
  if(threadIndex < nParticles){

  //Initialize RNG
  curandState_t state;
  curand_init(seed,threadIndex,0,&state);

  //Measure SD for every measurent time and every aging time
  {//Loop variables:
  int measurementIdx, SDIdx, measurementTime;
  double position;
  double oldTime=0;
  double newTime=0;
  double oldPosition = 0;
  double newPosition = 0;
  double startPosition = 0; //Set at beginning of each observation
  int stepCount=0;
  for(int taIdx = 0; taIdx != nAgingTimes; taIdx++){
    for(int timeIdx = 0; timeIdx != nTimes; timeIdx++){
      measurementTime = d_agingTimes[taIdx]+d_times[timeIdx];

      //Catch up to measurement time
      while(newTime <= measurementTime ){
        oldTime = newTime; oldPosition = newPosition;
        //Get new Position and time:
        nextStep(oldTime, oldPosition, newTime, newPosition, gamma, nu, t0, c, state);
        stepCount++;
      }

      //Find position at measurement time
      position = measurePosition(measurementTime, oldTime,
           oldPosition, newTime, newPosition,
           nu, eta, c);

      if(timeIdx==0){//Start of observation?
        startPosition = position;
      }


      //Add new SD and bin value
      measurementIdx = taIdx * nTimes+timeIdx;
      SDIdx = measurementIdx*nParticles+threadIndex;
      d_SD[SDIdx] = (position-startPosition)*(position-startPosition);
      d_bins[SDIdx] = positionToBin(position-startPosition, nBins, histogramRange);
    }// end loop time
  }// end loop tAge
  }

  //Reduction d_SD
  {//Loop variables:
  int measurementIdx;
  for(int taIdx = 0; taIdx != nAgingTimes; taIdx++){
    for(int timeIdx = 0; timeIdx != nTimes; timeIdx++){
      measurementIdx = taIdx * nTimes + timeIdx;
      //Do Reduction
      __syncthreads();
      for(unsigned int s = blockDim.x/2; s>0; s>>=1 ){
        if(threadIdx.x<s && threadIndex+s <nParticles){
            d_SD[measurementIdx*nParticles+threadIndex] +=
              d_SD[measurementIdx*nParticles+threadIndex+s];
        }
        __syncthreads();
      }
    }
  }
  }//End Reduction

} //endif (threadIndex < nParticles)
}
