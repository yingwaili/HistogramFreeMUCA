import numpy as np

def integrate(lngCoeff, lngBasis, xmin, xmax, ymin, ymax, dy):
    '''
            Numerical integration by a trapezoidal Riemann sum.

    Parameters
    ----------
    lngCoeff :
        - lngCoeff of the basis terms for the natural log of the instantaneous density of states, g
    lngBasis :
        - natural log of the density of states, g
    ymin, ymax :
        - defines the range of values that y can take
        - analogous to "energy range" in multicanonical or Wang-Landau schemes
    dy :
        - bin width of the Riemann sum

    Returns
    -------
    result of the integral

    '''

    integral = 0.0
    norm     = 0.0

    y = ymin
    while (y <= ymax-dy):
        y1 = y
        y2 = y + dy
        norm = norm + dy*(((np.exp(np.dot(np.asarray(lngCoeff),np.asarray(lngBasis(y1))))) + (np.exp(np.dot(np.asarray(lngCoeff),np.asarray(lngBasis(y2))))))/2.0)

        y = y + dy

    norm = (xmax-xmin)/norm
    y = ymin
    while (y <= ymax-dy):
        y1 = y
        y2 = y + dy
        I1 = (y*np.exp(np.dot(np.asarray(lngCoeff),np.asarray(lngBasis(y1)))))
        I2 = ((y+dy)*np.exp(np.dot(np.asarray(lngCoeff),np.asarray(lngBasis(y2)))))
        integral += dy*(I1+I2)/2.0
        y = y + dy
    return (integral*norm)
