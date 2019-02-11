#!/bin/bash
scp mbothe@pool12.physik.hu-berlin.de:/users/stud/mbothe/MA/CUDA_Simulation_RW/Results/*.txt /home/marius/MA/CUDA_Simulation_RW/Results/
files=Results/*.txt
python plotScript.py $files
