import sys
import numpy as np
import matplotlib.pyplot as plt

nPictures = len(sys.argv)
for picIdx in range(1,nPictures):
    path = sys.argv[picIdx];
    data = np.loadtxt(path,skiprows=0)


    #Read out Parameters:
    type = data[2,0]; #1 for MSD, 2 for MSDaging
    a = data[2,1];
    b = data[2,2];
    gamma = data[2,3];
    nu = data[2,4]
    eta = data[2,5]

    #Plot values
    skipValues=0;
    x=data[0,skipValues:-1]
    y=data[1,skipValues:-1]
    plt.figure()
    plt.xscale('log')
    plt.yscale('log')
    plt.grid(True)
    plt.plot(x,y,'bo',label='Data points')
    #Plot fit
    plt.plot(x, x**a * np.exp(b),label='Fit with exponent %s'%a);
    #Labels
    if type==1:
        plt.xlabel('Time')
    elif type==2:
        plt.xlabel('Aging Time')
    else:
        print('Error: unknown image type')
    plt.ylabel('Mean Squared Displacement')
    plt.title(r'Levy Walk ($\gamma =$%s, $\nu=$ %s, $\eta=$ %s)'%(gamma,nu,eta))
    plt.legend()
    plt.savefig(path[0:-3]+'png')
    plt.close()
