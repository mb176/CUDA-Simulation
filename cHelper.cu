#include <stdio.h>

using namespace std;

__host__ __device__ double** createMatrix(const int nVectors, const int vectorLength){
    if(nVectors<=0){
       printf("Error: positive number of vecors needed");
     }
     if(nVectors<=0){
       printf("Error: positive number of vecors element needed");
     }

    const int nEntries = nVectors*vectorLength;
    double *bigVector; //Stores all value linearly
    bigVector = (double*)malloc(sizeof(double)*nEntries);

    double** matrix;
    matrix = (double**)malloc(sizeof(double*)*nVectors);
    for(int idx = 0; idx!=nVectors; idx++){
        matrix[idx]=&bigVector[idx*vectorLength];
    }
    return matrix;
}

__host__ __device__ void freeMatrix(double** matrix){
  //Delete all cells
  free(matrix[0]);
  //Delete the pointer field
  free(matrix);
}



//Matrix doesnt work on Device :(
/*
 __host__  double** createMatrixOnDevice(const int nVectors, const int vectorLength){
    if(nVectors<=0){
       printf("Error: positive number of vecors needed");
     }
     if(nVectors<=0){
       printf("Error: positive number of vecors element needed");
     }

    const int nEntries = nVectors*vectorLength;
    double *d_bigVector; //Stores all value linearly
    cudaMalloc((void**)&d_bigVector, sizeof(double)*nEntries);
    double** d_matrix;
    cudaMalloc((void**)&d_matrix, sizeof(double*)*nVectors);

    double** matrix;
    matrix = (double**)malloc(sizeof(double*)*nVectors);
    for(int idx = 0; idx!=nVectors; idx++){
        matrix[idx]=&d_bigVector[idx*vectorLength];
    }
    cudaMemcpy(d_matrix, matrix, nVectors*sizeof(double*) ,cudaMemcpyHostToDevice);
    return d_matrix;
}

__host__ void freeMatrixOnDevice(double** matrix){
  //Delete all cells
  cudaFree(matrix[0]);
  //Delete the pointer field
  cudaFree(matrix);
}
*/

__host__ double* vectorToDevice(double* vector, int size){
  double* d_vector;
  cudaError(cudaMalloc((void**)&d_vector, size*sizeof(double)) );
  cudaError(cudaMemcpy(d_vector, vector, size*sizeof(double) ,cudaMemcpyHostToDevice));
  return d_vector;
}
