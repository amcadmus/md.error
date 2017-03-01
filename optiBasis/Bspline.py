#!/usr/bin/env python3

import sys
import numpy as np
import matplotlib.pyplot as plt

def test () :
    print ("here")

def _make_bspline_coeff (order) :    
    if (order != 4): 
        raise RuntimeError ("order other than 4 not implemented!")
    mat = np.zeros ((order, order), dtype = float)    
    mat[0][0] = 0.00000000000000000e+00
    mat[0][1] = 0.00000000000000000e+00
    mat[0][2] = 0.00000000000000000e+00
    mat[0][3] = 1.66666666666666657e-01
    
    mat[1][0] = 6.66666666666666630e-01
    mat[1][1] = -2.00000000000000000e+00
    mat[1][2] = 2.00000000000000000e+00
    mat[1][3] = -5.00000000000000000e-01
    
    mat[2][0] = -7.33333333333333304e+00
    mat[2][1] = 1.00000000000000000e+01
    mat[2][2] = -4.00000000000000000e+00
    mat[2][3] = 5.00000000000000000e-01
    
    mat[3][0] = 1.06666666666666661e+01
    mat[3][1] = -8.00000000000000000e+00
    mat[3][2] = 2.00000000000000000e+00
    mat[3][3] = -1.66666666666666657e-01
    return mat

def _make_bspline_d_coeff (order) :
    if (order != 4): 
        raise RuntimeError ("order other than 4 not implemented!")
    mat = np.zeros ((order, order-1), dtype = float)
    mat[0][0] = 0.00000000000000000e+00;
    mat[0][1] = 0.00000000000000000e+00;
    mat[0][2] = 5.00000000000000000e-01;

    mat[1][0] = -2.00000000000000000e+00;
    mat[1][1] = 4.00000000000000000e+00;
    mat[1][2] = -1.50000000000000000e+00;
    
    mat[2][0] = 1.00000000000000000e+01;
    mat[2][1] = -8.00000000000000000e+00;
    mat[2][2] = 1.50000000000000000e+00;
    
    mat[3][0] = -8.00000000000000000e+00;
    mat[3][1] = 4.00000000000000000e+00;
    mat[3][2] = -5.00000000000000000e-01;
    return mat


class Bspline (object) :
    def __init__ (self,
                  order = 4) :
        self.order = int(order)
        self.holder = self.order * 0.5
        self.coeff = _make_bspline_coeff (self.order)
        self.coeff_d = _make_bspline_d_coeff (self.order)

    def _element (self,
                  xx) :
        if (xx <= - self.holder or xx >= self.holder) : 
            return 0
        nx = xx + self.holder
        tmpa = self.coeff[int(nx)]
        buff = tmpa[self.order-1]
        for ii in range (self.order - 2, -1, -1) :
            buff = buff * nx + tmpa[ii]
        return buff

    def _element_d (self,
                    xx) :
        if (xx <= - self.holder or xx >= self.holder) : 
            return 0
        nx = xx + self.holder
        tmpa = self.coeff_d[int(nx)]
        buff = tmpa[self.order-2]
        for ii in range (self.order - 3, -1, -1) :
            buff = buff * nx + tmpa[ii]
        return buff
    
    def __call__ (self, 
                  xx) : 
        result = np.zeros (np.prod(xx.shape))
        count = 0
        for ii in np.nditer (xx) :
            result[count] = self._element(ii)
            count = count + 1
        result.reshape (xx.shape)
        return result

    def deriv (self, 
               xx) : 
        result = np.zeros (np.prod(xx.shape))
        count = 0
        for ii in np.nditer (xx) :
            result[count] = self._element_d(ii)
            count = count + 1
        result.reshape (xx.shape)
        return result
    
if __name__ == "__main__" :
    bspline = Bspline(4)
    x = np.arange (-5, 5, 0.01)
    y = bspline(x)
    print (x.shape)
    print (y.shape)
    plt.plot (x, y)
    plt.show()
