import sympy

def getExtrema(function, xmin, xmax):

    ''' return min and max of function f(x) over the closed interval [xmin, xmax] '''

    # convert input string to symbolic and lambda expressions
    x = sympy.Symbol("x")
    f_sym = sympy.sympify(function)
    f = sympy.lambdify(x, f_sym)

    # find maxima and minima, store in array "temp_extrema"
    df = sympy.diff(f_sym)
    temp_extrema = sympy.solve(df)

    # check if these values are within integration interval
    extrema = []
    i = 0
    j = 0
    while (i < len(temp_extrema)):
        while (j < len(sympy.solve(f_sym-temp_extrema[i]))):
            if ((float(sympy.re(sympy.solve(f_sym-temp_extrema[i])[j])) > xmin) and (float(sympy.re(sympy.solve(f_sym-temp_extrema[i])[j])) < xmax)):
                extrema.append(temp_extrema[i])
            j = j + 1
        i = i + 1
        j = 0

    # find value of function at beginning/end of interval
    extrema.append(f(xmin))
    extrema.append(f(xmax))

    # find min/max over interval
    ymin = float(min(extrema))
    ymax = float(max(extrema))

    return ymin, ymax
