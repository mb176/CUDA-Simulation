#!/bin/bash
#scp mbothe@pool12.physik.hu-berlin.de:/users/stud/mbothe/MA/CUDA_Simulation_RW/Results/*.txt /home/marius/MA/CUDA_Simulation_RW/Results/

#Aged histograms
folder=MA/CUDA_Simulation_RW/Results/Histogram_aged
files=histogram_gamma_0.6_nu_1.3_eta_1.3_ta_1e+05.txt
plotter=histogramAgedPlotterPresi.py

# #Ordinary histograms
# folder=MA/CUDA_Simulation_RW/Results/Histograms
# files=histogram_gamma_0.6_nu_1.3_eta_1.5_ta_0.txt
# plotter=histogramPlotterPresi.py

#PDF at origin fit
# folder=MA/CUDA_Simulation_RW/Results/Histogram_origin_peak
# files=origin*.txt
# plotter=originPeakPlotter.py

# #Fit of semi-analytical inverse laplace
# folder=MA/CUDA_Simulation_RW/Results/InverseLaplace
# files=histogram*.txt
# plotter=inverseLaplacePlotter.py

python Plotter/$plotter /home/marius/$folder/$files
# # send new pictures to university
# rsync -v /home/marius/$folder/* mbothe@pool10.physik.hu-berlin.de:$folder
