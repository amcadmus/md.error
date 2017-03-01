#!/usr/bin/env python3

import numpy as np
from Hermite import symm_hermite
import matplotlib.pyplot as plt

def plot_hermite (CC,
                  hh, 
                  vv, 
                  dd ) :
    sample_h = hh / float(100)
    xx = np.arange (-CC, CC + sample_h * 0.5, sample_h)
    value = np.zeros (len(xx))

    assert (len(vv) == len(dd))
    for ii in range(len(vv)) :
        [hv, hd] = symm_hermite (xx, ii, hh)
        # plt.plot (xx, hv)
        # plt.plot (xx, hd)
        # plt.show()
        value = value + vv[ii] * hv + dd[ii] * hd
    
    return [xx, value]


if __name__ == '__main__' :
    # vv = np.loadtxt ('spline.v.out')
    # dd = np.loadtxt ('spline.d.out')
    # vv = np.loadtxt ('tmp.v.out')
    # dd = np.loadtxt ('tmp.d.out')
    vv = np.loadtxt ('basis.1.out')
    dd = np.loadtxt ('deriv.1.out')
    
    CC = 2
    nbin = len(vv) - 1
    hh = CC / float(nbin)
    print ("nbin is %d", nbin)
    print ("hh   is %f", hh)
    [sx, sv] = plot_hermite (CC, hh, vv, dd)
    
    print_mat = [sx]
    print_mat = np.append (print_mat, [sv], axis=0)
    np.savetxt ('new.out', print_mat.T)
