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
    path = sys.argv[picIdx];
    data = np.genfromtxt(path,skip_header=1, skip_footer=0)
    params = np.genfromtxt(path,skip_header=0, skip_footer=len(data))

    #Read out Parameters:
    gamma =  params[0];
    nu = params[1];
    eta = params[2];
    tMax = params[3];
    nTimes = len(data)-1;

    #How many plots (called axes) do you want?
    xAxes = 1;
    yAxes = 1;
    nAxes = xAxes*yAxes; #number of plots in the figure


    #Plot values
    skipValues=0;
    position = data[0,0:-1]
    fontsize = 50;
    linewidth = 3.5; # for vertical lines

    fig = plt.figure(figsize=(25/np.sqrt(2),15))#(15/np.sqrt(2),10))
    plt.xticks(fontsize = fontsize);
    plt.yticks(fontsize = fontsize);
    rc('font', size=fontsize)
    #plt.xscale('log')

    # fig.subplots(yAxes,xAxes,sharx=True);
    for axesIdx in range(1,nAxes+1):
        ax = plt.subplot(yAxes,xAxes,axesIdx,yscale ='log')

        #Axes style
        ax.grid(True)
        ax.axes.ticklabel_format(axis='x',style='sci',scilimits = (0,0))

        # Calculate the times at which we show the histogram
        tIdxSteps = float(nTimes-1)/float(nAxes);
        tIdx = 1 + int((axesIdx) * tIdxSteps); #+1 beacause the first row of data are positions
        pdf = data[tIdx,0:-1]
        ax.plot(position,pdf,label='PDF of the walker')
        ax.fill_between(position,pdf)

        # boundry: plot vertical line indicating how far a particle can go in a single jump
        time = tMax/(nTimes-1)*(tIdx-1);
        boundry = time**nu;
        ax.axvline(boundry,color='r',label = 'Step with duration t', linestyle = '--',lw=linewidth)

        #Labels
        plt.xlabel('Displacement ($c t_0^{\\nu }$)',axes=ax)
        plt.ylabel('Relative number of particles',axes=ax)
        #plt.title(r't = %.2e $t_0$'%time,axes=ax,fontsize = fontsize)
        plt.legend()
    fig.savefig(path[0:-4]+'_presi'+'.png')
    plt.close()
