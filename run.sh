#!/bin/bash

files_cu="/users/stud/mbothe/MA/CUDA_Simulation_RW/main.cu /users/stud/mbothe/MA/CUDA_Simulation_RW/LevyWalkSimulation.cu"
files_h="/users/stud/mbothe/MA/CUDA_Simulation_RW/main.h /users/stud/mbothe/MA/CUDA_Simulation_RW/LevyWalkSimulation.h"
scp $files_cu $files_h mbothe@gdong4.physik.hu-berlin.de:/users/stud/mbothe/MA/CUDA_Simulation_RW

ssh  mbothe@gdong4.physik.hu-berlin.de "nvcc -o main ${files_cu}; MA/CUDA_Simulation_RW/main"
