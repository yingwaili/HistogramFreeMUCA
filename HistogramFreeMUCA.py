###
###  This is the main program of performing binless multicanonical sampling.
###
import numpy as np
from dataSet_methods import *
from correction import *
from integrate import *
from getExtrema import *
import matplotlib.pyplot as plt
import sys

#from pyinstrument import Profiler


# Read variables from input file
if (len(sys.argv) != 2):
    print('Error: Please Provide Input File as Command Line Argument')
    sys.exit()
else: 
    exec(open(sys.argv[1]).read())
    if (seed == -1):
        seed = np.random.randint(0,100000)
    print('Random Number Seed:', seed)
    np.random.seed(seed)

# Initial parameter, this should not be changed
iteration = 1

# Find global minima and maxima over defined interval
ymin, ymax = getExtrema(function, xmin, xmax)

# Temporary arrays for coefficients from an iteration and overall lng coefficients
lncCoeff = np.zeros(maxBasisTerms)
lngCoeff = np.zeros(maxBasisTerms)
lncBasis_sym, lncBasis, lngBasis_sym, lngBasis  = construct_basisSet(maxBasisTerms, ymax, ymin)

# Main loop of sampling
while (iteration <= maxIter):

    # Generate the data set D
    D = fill_dataSet(function, lngCoeff, lngBasis, points_per_dataSet, xmin, xmax, ymin, ymax, plotCDF)

    # Trim length of lncCoeff and lngCoeff arrays after first iteration
    # this fixes the number of terms in the expansion to numTerms
    # also need to reconstruct basis set, etc. with correct number of terms
    if (iteration == 1):
        
#        profiler = Profiler()
#        profiler.start()

        lncCoeff = correction(D, ymin, ymax, KS_cutoff, maxBasisTerms, lncBasis, plotFT)
        lncCoeff = np.trim_zeros(lncCoeff, 'b')
        numTerms = len(lncCoeff)
        lngCoeff = np.zeros(numTerms)
        lncBasis_sym, lncBasis, lngBasis_sym, lngBasis = construct_basisSet(numTerms, ymax, ymin)

#        profiler.stop()
#        print(profiler.output_text(unicode=True, color=True))

#        profiler.open_in_browser()
#        profiler.output_html()

    else:
        if (improved == 0):
            lncCoeff = correction(D, ymin, ymax, KS_cutoff, numTerms, lncBasis, plotFT)
        elif (improved == 1):
            lncCoeff = correction_improved(D, ymin, ymax, KS_cutoff, numTerms, lncBasis, s, plotFT)

    # Update remainder by adding correction coefficients from current iteration
    lngCoeff = lngCoeff + lncCoeff
    integral = integrate(lngCoeff, lngBasis, xmin, xmax, ymin, ymax, 0.01)
    print ('Iteration:', iteration, '\t Integral Estimate: %5.9f' % integral)
    iteration = iteration + 1

    # Plot DOS after each iteration
    if (plotDOS_iter == 1):
        show_DOS_Plot(iteration, lngCoeff, lngBasis, ymin, ymax)


print("Finished!")

# Plot DOS after maxIter Iterations
if (plotDOS_end == 1):
    show_DOS_Plot(iteration, lngCoeff, lngBasis, ymin, ymax)
