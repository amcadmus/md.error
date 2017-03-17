#!/usr/bin/env python3

import numpy as np
import argparse
# import matplotlib.pyplot as plt
from Hermite import symm_hermite
from Hermite import symm_hermite_deriv

def _parse_argument():
    parser = argparse.ArgumentParser(
        description="*** Print a smooth Hermite interpolated basis. ***")
    parser.add_argument('input_file', type=str, 
                        help='the basis file')
    parser.add_argument('-s','--sampling-points', type=int, default = 100,
                        help='the number of sampling points in each interval')
    parser.add_argument('-o','--output', type=str, default = "sampled.out",
                        help='the ouptput file')
    args = parser.parse_args()
    return args

def plot_hermite (CC,
                  hh, 
                  vv, 
                  dd ,
                  sampling) :
    sample_h = hh / float(sampling)
    xx = np.arange (-CC, CC + sample_h * 0.5, sample_h)
    value = np.zeros (len(xx))
    deriv = np.zeros (len(xx))

    assert (len(vv) == len(dd))
    for ii in range(len(vv)) :
        [hv, hd] = symm_hermite (xx, ii, hh)
        value = value + vv[ii] * hv + dd[ii] * hd
        [hv, hd] = symm_hermite_deriv (xx, ii, hh)
        deriv = deriv + vv[ii] * hv + dd[ii] * hd
    
    return [xx, value, deriv]


if __name__ == '__main__' :
    args = _parse_argument ()

    data = np.loadtxt (args.input_file)
    data = data.T
    xx = data[0]
    vv = data[1]
    dd = data[2]
    
    CC = 2
    nbin = len(vv) - 1
    hh = CC / float(nbin)
    print ("nbin is %d", nbin)
    print ("hh   is %f", hh)
    [sx, sv, sd] = plot_hermite (CC, hh, vv, dd, args.sampling_points)
    
    print_mat = [sx]
    print_mat = np.append (print_mat, [sv], axis=0)
    print_mat = np.append (print_mat, [sd], axis=0)
    np.savetxt (args.output, print_mat.T)
