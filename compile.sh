#!/bin/bash

#look up GPU: nvidia-smi
#gdong2 has Tesla M2050, corresponds to sm_50 -arch sm_50

computer="pool3"
executableName="PDF_origin_peak2";
files_cu="main_origin_peak.cu LevyWalkSimulation.cu d_LevyWalkGo.cu cHelper.cu vector_calculus.cpp fitting.cpp"
files_h="LevyWalkSimulation.h d_LevyWalkGo.h cHelper.h vector_calculus.h fitting.h"

# send the program files to HU server
rsync -v ~/MA/CUDA_Simulation_RW/* mbothe@$computer.physik.hu-berlin.de:MA/CUDA_Simulation_RW/

# Log into HU server to compile and execute the programs
ssh mbothe@$computer.physik.hu-berlin.de << EOF #Liste der Shell-Befehle
cd /users/stud/mbothe/MA/CUDA_Simulation_RW
module load cuda/9.0  #To avoid -fcpreprocessed error and allow c++11
nvcc --std c++11 -o $executableName --relocatable-device-code true ${files_cu}
echo "Compilation done"
# #Run executable directly
# ./$executableName
# echo "Simulation started"
# Create detached screen session on computer and run it there
screen -d -m -S simulation ./$executableName
echo "Simulation started in screen session"
EOF

#Copy Results back from server
rsync -auv mbothe@$computer.physik.hu-berlin.de:MA/CUDA_Simulation_RW/Results/ ~/MA/CUDA_Simulation_RW/Results/


#Independent Simulation: ssh, screen, execute main, strg+a, strg+d, exit ssh
