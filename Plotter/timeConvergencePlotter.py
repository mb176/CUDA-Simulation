import sys
import numpy as np
import matplotlib.pyplot as plt

# Skript should be executed on the .txt files that are created by main_histogram.cu. It reads out the
# parameters and plots the data as well as the fits

nPictures = len(sys.argv)
for picIdx in range(1,nPictures):
    path = sys.argv[picIdx];
    data = np.genfromtxt(path,skip_header=1, skip_footer=2)
    fits = np.genfromtxt(path,skip_header=4, skip_footer=1)
    params = np.genfromtxt(path,skip_header=5, skip_footer=0)

    #Read out Parameters:
    gamma =  params[0];
    nu = params[1];
    eta = params[2];
    originExp = fits[0]; #obtained by fitting it
    originOffset = fits[1];
    peakExp = fits[2];
    peakOffset = fits[3];

    #Plot values
    skipValues=0;
    time =data[0,0:-1]
    pOrigin =data[1,0:-1];
    pPeak = data[2,0:-1];
    plt.figure()
    plt.xscale('log')
    plt.yscale('log')
    plt.grid(True)
    plt.xlabel('Time')
    plt.ylabel('Probability density')
    plt.gcf().subplots_adjust(left=0.15)
    plt.plot(time,pOrigin,'bo',label='Probability at origin')
    # plt.plot(time,pPeak,'ro',label='Probability at delta peak')
    #Plot fit
    plt.plot(time[skipValues:], time[skipValues:]**originExp * np.exp(originOffset),'b',label='Fit with exponent %s '%(originExp));
    # plt.plot(time[skipValues:], time[skipValues:]**peakExp * np.exp(peakOffset),'r',label='Fit with exponent %s '%(peakExp));
    #Labels

    plt.title(r'Levy Walk ($\gamma =$%s, $\nu=$ %s, $\eta=$ %s)'%(gamma,nu,eta))
    plt.legend()
    plt.savefig(path[0:-3]+'png')
    plt.close()
