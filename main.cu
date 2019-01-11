#include <iostream>
#include <cuda.h>
#include <curand_kernel.h>
#include <vector>
#include <time.h>
#include <cmath>

#include "Levy_Walk.h"

using namespace std;
typedef vector<double> vec;
//#include "test.h"

int main(void)
{
  const int blocksize = 256;
  LevyWalk LW;
  //Set Physical Parameters
  LW.nu = 1;
  LW.gamma = 0.5;
  LW.eta = LW.nu;

  //experiment with distribution

  vec v(pow(10,6),1);
  double *d_v=&v[0];
  cudaMalloc((void**)&d_v, v.size()*sizeof(double));

  cudaMemcpy(d_v, &v[0], v.size()*sizeof(double) , cudaMemcpyHostToDevice);

  LevyWalkGo<<<(v.size()-1)/blocksize+1, blocksize>>>
  (d_v,v.size(), time(NULL), LW.gamma, LW.nu, LW.eta, LW.t0, LW.c, LW.dim);
  
  cudaMemcpy(&v[0], d_v, v.size()*sizeof(double) , cudaMemcpyDeviceToHost);
  //Sum up all blocks
  for(int i =1; i!= (v.size()+blocksize-1)/blocksize; i++){//(n+blockDim.x-1)/blockDim.x
    v[0] += v[i*blocksize];
  }

  // for(int i= v.size()-1;i>=0; i--){
  //   cout << v[i]<< endl;
  // }
  cout << "avr " << v[0]/v.size()<< endl;




  cudaFree(d_v);



  //Testing section
  // 0) test if member functions can be called globally
  // 1) Levy_Walk.go()
  // 2) Fitting

  //Initialize Simulation Block
  //name
  //dim
  //Run over every element of v_t_measure for each v_t_age
  //t_max, n_measure, v_t_measure
  //v_t_age
  //v_nu
  //v_gamma
  //eta
  //t_0
  //c
  //n_essemble
  //n_avr

  //Call Sim1.rum()
}
