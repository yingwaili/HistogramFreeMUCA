import numpy as np
import sympy
from plotMethods import *


def construct_basisSet(maxBasisTerms, ymax, ymin):
    '''
    Constructs symbolic expressions and lambda functions of the basis set and corresponding derivative.

    Parameters
    ----------
    maxBasisTerms :
        - maximum number of terms allowed in basis set expansion

    ymax, ymin :
        - min. and max. of f(x) within the integration domain

    Returns
    -------
    lncBasis_sym :
        - a list that stores the basis set for ln(c) in SymPy symbolic expressions
    lncBasis :     
        - a list that stores the basis set for ln(c) in SymPy lambda functions 
    lngBasis_sym : 
        - a list that stores the basis set for ln(g) in SymPy symbolic expressions
    lngBasis :
        - a list that stores the basis set for ln(g) in SymPy lambda functions 
    '''

    i  = 0
    q  = sympy.Symbol('q')                 # Define 'q' as a symbol, not a parameter
    lncBasis  = []
    lngBasis = []

    while (i < maxBasisTerms):
        lncBasis.append("sin(%f*%f*((q-%f)/(%f-%f)))" %(i, np.pi, ymin, ymax, ymin))
        i = i + 1
   
    lncBasis_sym = sympy.sympify(lncBasis)
    lncBasis = sympy.lambdify(q,lncBasis_sym)    # Applies expression 'lncBasis_sym' onto argument 'q'
    # ln(g) basis set - derivative of ln(c) basis
    i = 0
    while (i < maxBasisTerms):
        lngBasis.append(str(sympy.diff(lncBasis_sym[i], q)))
        i = i + 1
    lngBasis_sym = sympy.sympify(lngBasis)
    lngBasis = sympy.lambdify(q, lngBasis_sym)

    return lncBasis_sym, lncBasis, lngBasis_sym, lngBasis



def fill_dataSet(function, lngCoeff, lngBasis, points_per_dataSet, xmin, xmax, ymin, ymax, plotCDF):
    '''
    Generates data set; calculates the cumulative distribution function (CDF) and reminder.

    Parameters
    ----------
    function :
        - symbolic expression of the integrand, f(x)
    lngCoeff :
        - coefficients of the basis terms for the natural log of the instantaneous density of states, g
    lngBasis :
        - natural log of the instantaneous density of states, g
        - it acts as the sampling weight for determining acceptance
    points_per_dataSet :
        - number of data points in the data set
        - determines the upper bound of the data array
    xmin, xmax :
        - integration domain
    ymin, ymax :
        - min. and max. of f(x) within the integration domain
    plotCDF :
        - a flag to plot a graph of the CDF; 0 = off, 1 = on

    Returns
    -------
    data :
        - an array with the following indexing scheme:
           [data point number, data values]
          Data values as follows:  [x,0]: sorted data set
                                   [x,1]: CDF
                                   [x,2]: remainder (CDF - straight line)
    '''

    data = np.zeros((points_per_dataSet, 3))

    # Create Lambda function from input function
    x = sympy.Symbol("x")
    f_sym = sympy.sympify(function)
    f = sympy.lambdify(x, f_sym)

    # Fill data set: Accept or reject data points based on acceptance criterion
    j = 0

    rand_xinit = np.random.uniform(xmin, xmax)
    yold = f(rand_xinit)
    while (j < points_per_dataSet):            
        rand_x = np.random.uniform(xmin, xmax)
        y = f(rand_x)
        g_new = np.exp(np.dot(np.asarray(lngCoeff),np.asarray(lngBasis(y))))
        g_old = np.exp(np.dot(np.asarray(lngCoeff),np.asarray(lngBasis(yold))))

        if (g_new < g_old):
            data[j,0] = y
            yold = y
        elif ((g_old/g_new) > np.random.uniform(0,1)): 
            data[j,0] = y
            yold = y
        else:
            data[j,0] = yold
        j = j + 1

    # Sort data set
    data[:,0] = np.sort(data[:,0])

    # Compute CDF of data set
    j = 0
    while (j < points_per_dataSet):
            data[j,1] = float ((j+1.0)/points_per_dataSet)
            j = j + 1

    # Compute remainder
    j = 0
    temp = np.zeros((points_per_dataSet))
    while (j < points_per_dataSet):
        data[j,2] = (data[j,1] - ((data[j,0] - data[0,0]) / 
                                      (data[points_per_dataSet - 1,0] - data[0,0])))
        temp[j] = ((data[j,0] - data[0,0]) /
                                      (data[points_per_dataSet - 1,0] - data[0,0]))
        j = j + 1

    # Plot ECDF, F_0, and empirical remainder
    if (plotCDF == 1):
        show_CDF_Plot(data,temp)

    return data
