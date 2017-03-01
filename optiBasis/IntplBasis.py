#!/usr/bin/env python3

import sys
import numpy as np
from scipy.interpolate import interp1d

class IntplBasis (object) :
    def __init__ (self,
                  xx,
                  vv,
                  CC) :
        self.func = interp1d (xx, vv, axis=0, kind="cubic")
        self.CC = CC
        
    def __call__ (self,
                  xx) :
        result = np.zeros (xx.shape)
        absxx = np.fabs(xx)
        for ii in range(len(xx)) :
            if (absxx[ii] < self.CC) : 
                result[ii] = self.func (absxx[ii])
        return result
