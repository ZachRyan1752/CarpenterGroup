import numpy as np

## One dimensional:
### Polynomial
def Linear(x, m, b):
    y = m * x + b
    return y

def Polynomial_1N(x, m, b):
    y = m * x + b
    return y

def Polynomial_2N(x, m1, m2, b):
    y = m1 * x**2 + m2 * x + b
    return y

def Polynomial_3N(x, m1, m2, m3, b):
    y = m1 * x**3 + m2 * x**2 + m3 * x + b
    return y
    
def Polynomial_4N(x, m1, m2, m3, m4, b):
    y = m1 * x**4 + m2 * x**3 + m3 * x**2 + m4 * x + b
    return y

def Polynomial_5N(x, m1, m2, m3, m4, m5, b):
    y = m1 * x**5 + m2 * x**4 + m3 * x**3 + m4 * x**2 + m5 * x + b
    return y
    
### Exponential
def Exponential_1N(x, a, b, c):
    y = a * b**x + c
    return y


## Two Dimensional:
def twoD_Gaussian(xy, amplitude, xo, yo, sigma_x, sigma_y, theta, offset): # Taken from: https://stackoverflow.com/questions/21566379/fitting-a-2d-gaussian-function-using-scipy-optimize-curve-fit-valueerror-and-m
    x, y = xy
    xo = float(xo)
    yo = float(yo)    
    a = (np.cos(theta)**2)/(2*sigma_x**2) + (np.sin(theta)**2)/(2*sigma_y**2)
    b = -(np.sin(2*theta))/(4*sigma_x**2) + (np.sin(2*theta))/(4*sigma_y**2)
    c = (np.sin(theta)**2)/(2*sigma_x**2) + (np.cos(theta)**2)/(2*sigma_y**2)
    g = offset + amplitude*np.exp( - (a*((x-xo)**2) + 2*b*(x-xo)*(y-yo) 
                            + c*((y-yo)**2)))
    return g.ravel()

def Double_twoD_Gaussian(xy, amplitude, xo, yo, sigma_x, sigma_y, theta, offset, amplitude2, xo2, yo2, sigma_x2, sigma_y2, theta2, offset2): # Taken from: https://stackoverflow.com/questions/21566379/fitting-a-2d-gaussian-function-using-scipy-optimize-curve-fit-valueerror-and-m
    x, y = xy
    xo = float(xo)
    yo = float(yo)    
    a = (np.cos(theta)**2)/(2*sigma_x**2) + (np.sin(theta)**2)/(2*sigma_y**2)
    b = -(np.sin(2*theta))/(4*sigma_x**2) + (np.sin(2*theta))/(4*sigma_y**2)
    c = (np.sin(theta)**2)/(2*sigma_x**2) + (np.cos(theta)**2)/(2*sigma_y**2)
    g = offset + amplitude*np.exp( - (a*((x-xo)**2) + 2*b*(x-xo)*(y-yo) 
                            + c*((y-yo)**2)))
    
    xo2 = float(xo2)
    yo2 = float(yo2)    
    a2 = (np.cos(theta2)**2)/(2*sigma_x2**2) + (np.sin(theta2)**2)/(2*sigma_y2**2)
    b2 = -(np.sin(2*theta2))/(4*sigma_x2**2) + (np.sin(2*theta2))/(4*sigma_y2**2)
    c2 = (np.sin(theta2)**2)/(2*sigma_x2**2) + (np.cos(theta2)**2)/(2*sigma_y2**2)
    g2 = offset2 + amplitude2*np.exp( - (a2*((x-xo2)**2) + 2*b2*(x-xo2)*(y-yo2) 
                            + c2*((y-yo2)**2)))
    
    g3 = g + g2
    return g3.ravel()
