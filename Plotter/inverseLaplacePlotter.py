import sys
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc;
from matplotlib import rcParams
rcParams.update({'figure.autolayout': True})

# Skript should be executed on the .txt files that are created by safeHistogram().
# It plots histograms at different points in time

nPictures = len(sys.argv)
for picIdx in range(1,nPictures):
    # Read data
    # These files are created by copying the the inverse laplace values from
    # Mathematica to the bottom of the histogram files
    path = sys.argv[picIdx];
    data = np.genfromtxt(path,skip_header=1, skip_footer=2)
    params = np.genfromtxt(path,skip_header=0, skip_footer=len(data)+2)
    laplaceData = np.genfromtxt(path,skip_header=len(data)+1, skip_footer=0);

    #Read out Parameters:
    gamma =  params[0];
    nu = params[1];
    eta = params[2];
    tMax = params[3];
    nTimes = len(data)-1;

    # Parameters to rescale the Laplace Transform (handcopyed into the .txt file)
    # stretch = params[4];

    #Plot values
    skipValues=0;
    position = data[0,:]
    fontsize = 24;

    fig = plt.figure(figsize=(15/np.sqrt(2),10))
    #plt.xscale('log')

    # Plot histogram from simulation
    ax = plt.subplot(1,1,1,yscale ='log')
    ax.grid(True)
    #ax.axes.ticklabel_format(style='sci',scilimits = (0,0))
    pdf = data[-1,:] #Plot last time step
    ax.plot(position,pdf,label='Simulation')
    # ax.fill_between(position,pdf)

    # Plot the numerical inverse Laplace Transform of the analytic result

    # rescale y axis
    laplaceData[1,:] =  laplaceData[1,:] *pdf[0]/ laplaceData[1,0]

    # Cutoff values below the accuracy of the simulation
    cutoff = 0.3*10**(-7);
    indices =  laplaceData[1,:] > cutoff;
    analyticPosition = laplaceData[0,indices];
    analyticPdf =  laplaceData[1,indices];

    # Rescale x axis
    print(np.max(position)/np.max(analyticPosition))
    if(len(params)==5):
        stretch = params[4];
    else:
        stretch = position[-1]/analyticPosition[-1]
    analyticPosition = analyticPosition*stretch;
    ax.plot(analyticPosition,analyticPdf, label='Analytical result')


    #Labels
    plt.xticks(fontsize = fontsize);
    plt.yticks(fontsize = fontsize);
    rc('font', size=fontsize)
    plt.xlabel('Displacement ($c t_0^{\\nu }$)',axes=ax,fontsize=fontsize)
    plt.ylabel('Probability density',axes=ax,fontsize=fontsize)

    # plt.title(r't = %s'%time,axes=ax)
    plt.legend()
    fig.savefig(path[0:-3]+'png')
    plt.close()
