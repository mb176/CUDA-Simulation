#ifndef LEVY_WALK_GUARD
#define LEVY_WALK_GUARD

struct LevyWalk{
  //Simulation Parameters
  std::vector<double> tMeasure, tAge;
  unsigned int n_essemble;

  //Physical Parameters
  double gamma, eta, nu;
  int dim = 1; //No proper implementation of probability distribution for d>1
  double t0 = 1;
  double c = 1;
};


__global__ void LevyWalkGo(double * , int, int, double gamma,
  double nu, double eta, double t0, double c, int dim);


#endif
