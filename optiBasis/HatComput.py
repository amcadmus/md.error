#!/usr/bin/env python3

import sys
import numpy as np
from Bspline import Bspline
import matplotlib.pyplot as plt

def _naive_hat_bspline4 (mm) :
    if mm == 0: return 1
    tmp = np.pi * mm
    return np.power (np.sin(tmp)/tmp, 4)

def _naive_test_bspline4 (mm, KK) :
    return _naive_hat_bspline4(mm / float(KK)) / float(KK)

class HatComput (object) :
    def __init__ (self,
                  basis,
                  KK = 32,
                  over_smpl = 30,
                  ) :
        self.basis = basis
        self.KK = KK
        self.over_smpl = over_smpl
        self.MM = self.KK * self.over_smpl
        
        xx = np.arange (0, self.KK, float(self.KK) / float(self.MM))
        # print (len(xx))
        for ii in range (len(xx)) :
            if xx[ii] > float(self.KK / 2.) : 
                xx[ii] = xx[ii] - self.KK

        fx = self.basis (xx)
        ffx = np.fft.fft (fx) / len(fx)
        self.phim = np.real(ffx)

        # # naive check with Bspline 4
        # upper=64
        # mm = np.arange (0, upper)
        # naive_hat = np.zeros (upper)
        # for ii in range (upper) : 
        #     naive_hat[ii] = _naive_test_bspline4 (ii, self.KK) / 2        
        # plt.plot (np.real(ffx[0:upper]))
        # plt.plot (np.imag(ffx[0:upper]))
        # plt.plot (naive_hat)
        # plt.plot ((naive_hat - np.real(ffx[0:upper])))
        # plt.show()
        # plt.plot (np.arange (0, self.KK, float(self.KK) / float(self.MM)), fx)
        # plt.show ()

    def __call__ (self,
                  mm) :
        if mm < 0 : mm = -mm
        return self.phim[mm]
    

if __name__ == "__main__" : 
    bs = Bspline(4)
    hc = HatComput (bs, KK = 32, over_smpl = 300)
