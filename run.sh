#!/bin/bash

files_cu="main.cu LevyWalkSimulation.cu d_LevyWalkGo.cu cHelper.cu vector_calculus.cpp fitting.cpp"
files_h="LevyWalkSimulation.h d_LevyWalkGo.h cHelper.h vector_calculus.h fitting.h"

scp $files_cu $files_h mbothe@pool12.physik.hu-berlin.de:/users/stud/mbothe/MA/CUDA_Simulation_RW
ssh mbothe@pool12.physik.hu-berlin.de << EOF #Liste der Shell-Befehle
cd /users/stud/mbothe/MA/CUDA_Simulation_RW
module load cuda/9.0  #To avoid -fcpreprocessed error and allow c++11
nvcc --std c++11 -o main --relocatable-device-code true ${files_cu}
echo "End of compilation"
./main
EOF
