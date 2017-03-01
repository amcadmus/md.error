#!/usr/bin/env python3

import sys
import numpy as np
import scipy as sp
import argparse
from Region import Region
from Bspline import Bspline
from IkError import HermiteLossFunc

def _parse_argument():
    parser = argparse.ArgumentParser(
        description="*** Compute the optimal basis for SPME. Error estimated for water like homogeneous system. ***")
    parser.add_argument('-b', '--beta', type=float, default = 3.5,
                        help='the splitting parameter')
    parser.add_argument('-k', '--numb-grid', type=int, nargs = '+', default = [32],
                        help='number of grid points')
    parser.add_argument('-l', '--box-size', type=float, nargs = '+', default = [3.72412],
                        help='box size')
    parser.add_argument('-c', '--cut-off', type=int, default = 2,
                        help='cut-off radius of basis')
    parser.add_argument('-n', '--numb-bin', type=int, default = 10,
                        help='number of bins for piecewise cubic Hermite')
    args = parser.parse_args()
    return args

def _main () :
    args = _parse_argument()
        
    if len(args.box_size) == 1 : 
        box_size = [args.box_size[0], args.box_size[0], args.box_size[0]]
    else :
        box_size = args.box_size
    assert (len(box_size) == 3)

    if len(args.numb_grid) == 1 : 
        numb_grid = [args.numb_grid[0], args.numb_grid[0], args.numb_grid[0]]
    else :
        numb_grid = args.numb_grid
    assert (len(numb_grid) == 3)
    
    KK = numb_grid
    print (KK)
    LL = np.array(box_size)
    beta = args.beta
    region = Region (LL)    
    CC = args.cut_off
    nbin = args.numb_bin    

    # make a water like system
    q2 = 33.456 * np.prod(LL) * (-0.8476 * -0.8476 + 0.4238 * 0.4238 * 2)
    natoms = 33.456 * np.prod(LL) * 3

    lossfunc = HermiteLossFunc (CC, nbin, beta, KK, q2, natoms, region)
    
    bstep = CC / float(nbin)
    bs = Bspline(CC*2)    
    bx = np.arange (0, CC + .5 * bstep, bstep)
    bv_ = bs(bx)
    bd_ = bs.deriv (bx)
    scale = bv_[0]
    bv = bv_[1:]
    bd = bd_[1:]
    bv = bv / scale
    bd = bd / scale
    init_vv = np.append (bv, bd)

    aa = sp.optimize.minimize (lossfunc.value, init_vv, jac = lossfunc.deriv, method='BFGS', options={'disp': True}, tol=1e-4)
    

if __name__ == '__main__':
    _main()
