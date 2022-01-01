import sys
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc;
from matplotlib import rcParams
rcParams.update({'figure.autolayout': True})

# Skript should be executed on the .txt files that are created by safeResults (ie MSD data). It reads out the
# parameters and plots the data as well as the fits


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
    prediction = data[2,6];



    #Plot values
    skipValues=2;
    x=data[0,0:-1]
    y=data[1,0:-1]

    plt.figure()

    #Style
    # size = 3.5
    # fig = plt.figure(figsize=(4*size/np.sqrt(2),2*size))
    plt.xscale('log')
    plt.yscale('log')
    # plt.grid(True)
    fontsize=14
    plt.xticks(fontsize = fontsize);
    plt.yticks(fontsize = fontsize);
    rc('font', size=fontsize)

    plt.plot(x[0:],y[0:],'bo',label='Data points')
    #Plot fit
    plt.plot(x[skipValues:], x[skipValues:]**a * np.exp(b),label='f(x) = a $x^b$, b = %.6s (prediction: %s)'%(a,prediction));
    #Labels
    if type==1:
        plt.xlabel('Time ($t_0$)', fontsize = fontsize)
    elif type==2:
        plt.xlabel('Aging time ($t_0$)', fontsize = fontsize)
    else:
        print('Error: unknown image type')
    plt.ylabel("MSD ($c^2 t_0^{2 \\nu }$)", fontsize = fontsize)
    #plt.title(r'Levy walk ($\gamma =$%s, $\nu=$ %s, $\eta=$ %s)'%(gamma,nu,eta))
    plt.legend()
    plt.savefig(path[0:-3]+'png')
    plt.close()
