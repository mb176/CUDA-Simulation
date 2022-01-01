#ifndef cHelper_GUARD
#define cHelper_GUARD

__host__ __device__ double** createMatrix(const int, const int);

__host__ __device__ void freeMatrix(double** );

__host__  double** createMatrixOnDevice(const int , const int );

__host__ void freeMatrixOnDevice(double**);

__host__ double* vectorToDevice(double* vector, int size);

#endif
