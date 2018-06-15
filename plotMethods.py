###
###  This file contains the functions to plot various quantities of interest during the simulation.
###
import numpy as np
import matplotlib.pyplot as plt


def show_CDF_Plot(data, temp):

    ''' Plots the empirical cumulative distribution function (ECDF), 
    empirical remainder, and the "flat histogram" reference (straight line F_0). '''

    plt.plot(data[:,0],data[:,1],label = '$\overline{F}(E)$')
    plt.plot(data[:,0],temp,label = '$F_0(E)$')
    plt.plot(data[:,0],data[:,2],label = '$\overline{R}(E)$', color='green')
    plt.title('Empirical Cumulative Distribution Function',fontsize = 16)
    plt.xlabel('$E$',fontsize = 16)
    plt.ylabel('Cumulative Probability', fontsize = 16)
    plt.legend(prop={'size': 16})
    plt.show()



def show_DOS_Plot(iteration, lngCoeff, lngBasis, ymin, ymax):

    ''' Plot of the natural log of the density of states '''

    plt.title('Density of States, Iteration %d' % (iteration-1) ,fontsize = 16)
    plt.xlabel('$E$',fontsize = 16)
    plt.ylabel('$ln$ $g(E)$',fontsize = 16)
    x = np.linspace(ymin,ymax,1000)
    y = np.zeros(1000)
    i = 0
    while (i < 1000):
        y[i] = np.dot(np.asarray(lngCoeff),np.asarray(lngBasis(x[i])))
        i = i + 1
    plt.plot(x,y)
    plt.show()



def show_FT_Plot(dataSet, lncCoeff, lncBasis, ymin, ymax):

    ''' Plot of the empirical remainder vs. the analytical remainder '''

    plt.title('Empirical and Analytic Remainders' ,fontsize = 16)
    plt.xlabel('$E$',fontsize = 16)
    plt.ylabel('Remainder',fontsize = 16)
    numPoints = np.size(dataSet[:,0])
    lnC = np.zeros(numPoints)
    i = 0
    while (i < numPoints):
        lnC[i] = np.dot(np.asarray(lncCoeff),np.asarray(lncBasis(dataSet[i,0])))
        i = i + 1
    plt.plot(dataSet[:,0],lnC, label = '$R(E)$', color='red')
    plt.plot(dataSet[:,0],dataSet[:,2], label = '$\overline{R}(E)$',color='green')
    plt.legend(prop={'size': 16})
    plt.show()
