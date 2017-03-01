#!/usr/bin/env python3

import sys
import numpy as np
from Bspline import Bspline
from Hermite import symm_hermite
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
    
class HermiteBasisHatComput (object) :
    def __init__ (self, 
                  CC,
                  n_bin,
                  KK = 32,
                  over_smpl = 30,
                  ) :
        self.CC = CC
        self.n_bin = n_bin
        self.hh = self.CC / self.n_bin
        self.KK = KK
        self.over_smpl = int(over_smpl)
        self.MM = self.KK * self.over_smpl
        
        xx = np.arange (0, self.KK, float(self.KK) / float(self.MM))
        # print (len(xx))
        for ii in range (len(xx)) :
            if xx[ii] > float(self.KK / 2.) : 
                xx[ii] = xx[ii] - self.KK
        
        self.hathv = []
        self.hathd = []

        # at knot 0
        [hv, hd] = symm_hermite (xx, 0, self.hh)
        hd = np.zeros (hd.shape)
        fhv = np.fft.fft (hv) / len(hv)
        # fhd = np.fft.fft (hd) / len(hd)
        self.hathv = [np.real(fhv)]
        # np.append (self.hathd, np.real(fhd))
        # plt.plot (np.arange (0, self.KK, float(self.KK) / float(self.MM)), hv)
        # plt.plot (np.arange (0, self.KK, float(self.KK) / float(self.MM)), hd)
        # plt.show()

        # at "normal" knot
        for ii in range (1, self.n_bin) :
            [hv, hd] = symm_hermite (xx, ii, self.hh)
            fhv = np.fft.fft (hv) / len(hv)
            fhd = np.fft.fft (hd) / len(hd)
            if len(self.hathv) == 0 : self.hathv = [np.real(fhv)] 
            else : self.hathv = np.append (self.hathv, [np.real(fhv)], axis = 0)
            if len(self.hathd) == 0 : self.hathd = [np.real(fhd)] 
            else : self.hathd = np.append (self.hathd, [np.real(fhd)], axis = 0)
            # plt.plot (np.arange (0, self.KK, float(self.KK) / float(self.MM)), hv)
            # plt.plot (np.arange (0, self.KK, float(self.KK) / float(self.MM)), hd)
            # plt.show()

        # at boundary knot
        [hv, hd] = symm_hermite (xx, self.n_bin, self.hh)
        for ii in range (len(xx)) : 
            if xx[ii] > CC : 
                hv[ii] = 0
                hd[ii] = 0
            elif xx[ii] < -CC : 
                hv[ii] = 0
                hd[ii] = 0
        fhv = np.fft.fft (hv) / len(hv)
        fhd = np.fft.fft (hd) / len(hd)
        if len(self.hathv) == 0 : self.hathv = [np.real(fhv)] 
        else : self.hathv = np.append (self.hathv, [np.real(fhv)], axis = 0)
        if len(self.hathd) == 0 : self.hathd = [np.real(fhd)] 
        else : self.hathd = np.append (self.hathd, [np.real(fhd)], axis = 0)
        # plt.plot (np.arange (0, self.KK, float(self.KK) / float(self.MM)), hv)
        # plt.plot (np.arange (0, self.KK, float(self.KK) / float(self.MM)), hd)
        # plt.show()
        
    def set_value (self, vi, di) :
        self.vi = vi
        self.di = di
        assert (len(self.hathv) == len(self.vi))
        assert (len(self.hathd) == len(self.di))
        
    def __call__ (self, mm) :
        if mm < 0 : mm = -mm
        result = 0
        result = result + np.dot (self.vi, self.hathv[:,mm])
        result = result + np.dot (self.di, self.hathd[:,mm])
        # for ii in range(len(self.vi)) :
        #     result = result + self.vi[ii] * self.hathv[ii][mm]
        # for ii in range(len(self.di)) :
        #     result = result + self.di[ii] * self.hathd[ii][mm]
        return result

    def basis_value (self, mm) :
        return np.append (self.hathv[:,mm], self.hathd[:,mm])
        

if __name__ == "__main__" : 
    bs = Bspline(4)
    # hc = HatComput (bs, KK = 32, over_smpl = 300)
    
    CC = 2.
    nbin = 100
    KK = 32
    hhc = HermiteBasisHatComput (CC, nbin, KK, over_smpl = 80)

    hh = CC / float(nbin)
    xx = np.arange (0, CC + hh, hh)
    vv = bs (xx)
    dd_ = bs.deriv (xx)
    dd = dd_[1:len(dd_)]
    print_mat = [xx]
    print_mat = np.append (print_mat, [vv], axis = 0)
    print_mat = np.append (print_mat, [dd_], axis = 0)
    np.savetxt ('bspline.out', print_mat.T)

    hhc.set_value (vv, dd)
    
    # naive check with Bspline 4
    upper=64
    mm = np.arange (0, upper)
    naive_hat = np.zeros (upper)
    cmp_hat = np.zeros (upper)
    for ii in range (upper) : 
        naive_hat[ii] = _naive_test_bspline4 (ii, KK) 
        cmp_hat[ii] = hhc(ii)
    # plt.plot (cmp_hat)
    # plt.plot (naive_hat)
    plt.plot (naive_hat - cmp_hat)
    plt.show()
    
    # print (xx[0:10])
    # print (vv[0:10])
    # print (dd_[0:10])
    # plt.plot(xx, vv)
    # plt.plot(xx, dd_)
    # plt.show()
