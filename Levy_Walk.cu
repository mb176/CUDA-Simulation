#include<cmath>
#include<iostream>
#include<vector>
#include<time.h>

#include "curand_kernel.h"

using namespace std;



  //Levy_Walk::go()
    //Barricades
      //intervalls between t_age < biggest t_measure
      //only positive times
    //get vector pointers for v_t_age,  v_t_measure
    // allocate vectors for v_t_age*v_t_measure MSD
    //load vectors to cuda
    //call __global__ Levy_Walk_go<<<(n_essemble+blocksize-1)/blocksize,blocksize>>>
    //get vectors back from device
    //return vectors

__global__ void LevyWalkGo(double* d_test, int n, int seed, double gamma,
  double nu, double eta, double t0, double c, int dim){
  /*( uint n_essemble,int *d_v_t_age, int v_t_age_size,  int *d_v_t_measure, int v_t_measure_size,
double** d_MSD) */
  //Testing
  const int tidx = blockDim.x*blockIdx.x+threadIdx.x;

  //Initialize RNG
  curandState_t state;
  curand_init(seed,tidx,0,&state);
  for(int i =0; i!= n; i++){
    double number = curand_uniform_double(&state);
    double step_duration = t0*(pow(number,-1/gamma)-1);
    d_test[i] = step_duration;
  }


  //Reduction
  for(unsigned int s = blockDim.x/2; s>0; s>>=1 ){
    if(threadIdx.x<s && tidx+s <n){
        d_test[tidx] += d_test[tidx+s];
    }
    __syncthreads();
  }




  //Save times and positions of the particle at recent turning points:
  // double* old_position; //set to 0; for
  // double* new_position;
  // double* step; //vector of current step
  // double old_t = 0; // set to zero
  // double new_t = 0; //set to zero

  //Save time and position at the measurement times:
  // double *start_position; //start position of current measurement
  // double* position;
  // double start_time; //start time of current measurement
  // double t=0;

  // Allocate shared space to write output for the entire block to avoid
  // access conflict:
  //__shared__ double **tmp_output

  //allocate Memory for vectors

  //allocate Memory for tmp_output

  //Invariant: The measurements for idx_t_age t_ages are complete
  //for(idx_t_age=0; idx_t_age!=v_t_age_size; idx_t_age++){

    //Invariant: The measurements for idx_t_measure t_measures are complete
    //for(idx_t_measure=0; idx_t_measure !=v_t_age_size; idx_t_measure++){

      //time of next measurement:
      //t = d_v_t_age[idx_t_age]+d_v_t_measure[idx_t_measure]);

      //while(new_t <= t){
        //Invariant: the particle has not progressed past the current measurement time
        //generate one step
        //new_t = old_t + t_step;
        // new_position = old_position + step; //as vectors
      //}

    //find_position()
    //if(idx_t_measure == 0){ //set start for new measurement
      //start_time = t;
      //start_position = position;
    //}

    //find MSD

    //write to MSD to tmp_output;

    //} //end idx_t_measure
  //}//end idx_t_age

  //Last kernel writes to global output;
  //if(threadIdx.x = blockDim.x){
    //make sure all threads are done:
    //__syncthreads();
    //write to output
  //}

  //Free device memory:
  //cudaFree(old_position);
  //cudaFree(new_position);
  //cudaFree(step);
  //cudaFree(start_position);
  //cudaFree(position);
  //cudaFree(tmp_output);
}

  // }

// __device__ find_position(double old_t, double new_t, new_position, old_position
// )
