#!/bin/bash

#look up GPU: nvidia-smi
#gdong2 has Tesla M2050, corresponds to sm_50 -arch sm_50
files_cu="main.cu LevyWalkSimulation.cu d_LevyWalkGo.cu cHelper.cu vector_calculus.cpp fitting.cpp"
files_h="LevyWalkSimulation.h d_LevyWalkGo.h cHelper.h vector_calculus.h fitting.h"
computer="pool4"
executableName="main3";
scp $files_cu $files_h mbothe@$computer.physik.hu-berlin.de:/users/stud/mbothe/MA/CUDA_Simulation_RW
ssh mbothe@$computer.physik.hu-berlin.de << EOF #Liste der Shell-Befehle
cd /users/stud/mbothe/MA/CUDA_Simulation_RW
module load cuda/9.0  #To avoid -fcpreprocessed error and allow c++11
nvcc --std c++11 -o $executableName --relocatable-device-code true ${files_cu}
#if [$? -eq 0 ]; then # $? is the exit status of the last command
echo "Compilation done"
#./$executableName
#else
#  echo "Compilation failed"
#fi
EOF
#Independent Simulation: ssh, screen, execute main, strg+a, strg+d, exit ssh
