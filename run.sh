#!/bin/bash
nvcc -o main main.cu Levy_Walk.cu
echo "Compilation done!"
./main
