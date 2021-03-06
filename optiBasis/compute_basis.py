#!/usr/bin/env python3

import sys
import numpy as np
import scipy as sp
from scipy import optimize
import argparse
from Region import Region
from Bspline import Bspline
from IkError import HermiteLossFunc
from IkError import HermiteLossFunc_Bound0
import basis_common as bc

def _parse_argument():
    parser = argparse.ArgumentParser(
        description="*** Compute the optimal basis for SPME. Error estimated for water like homogeneous system. ***")
    parser.add_argument('-n', '--numb-bin', type=int, default = 10,
                        help='number of bins for piecewise cubic Hermite')
    parser.add_argument('-c', '--cut-off', type=int, default = 2,
                        help='cut-off radius of basis')
    parser.add_argument('-b', '--beta', type=float, default = 3.5,
                        help='the splitting parameter')
    parser.add_argument('-k', '--numb-grid', type=int, nargs = '+', default = [16],
                        help='number of grid points')
    parser.add_argument('-l', '--box-size', type=float, nargs = '+', default = [1.8621],
                        help='box size')
    parser.add_argument('-t', '--tolerence', type=float, default = 1e-2,
                        help='the tolerence of convergence')
    parser.add_argument('--l-cut', type=int, default = 0,
                        help='the cut-off of l sum. Use 0 for an estimate')
    parser.add_argument('-d', '--discontinuous-bd', action = 'store_true',
                        help='the boundary of the basis can be discontinuous')
    parser.add_argument('-o', '--output', type=str, default = "basis.out",
                        help='the output file of basis')
    args = parser.parse_args()
    return args

def _main () :
    args = _parse_argument()
        
    box_size = bc.make_dim3_vec (args.box_size)
    numb_grid = bc.make_dim3_vec (args.numb_grid)
    
    KK = numb_grid
    LL = np.array(box_size)
    beta = args.beta
    region = Region (LL)    
    CC = args.cut_off
    nbin = args.numb_bin    
    l_cut = args.l_cut
    if args.discontinuous_bd :        
        vanish_boundary = False
    else :
        vanish_boundary = True
    if (l_cut <= 0) :
        l_cut = bc.estimate_l_cut (CC, nbin, vanish_boundary)
        # print ("# estimated l_cut is %d" % l_cut)

    print ("# numb bins         %d" % nbin)
    print ("# cut-off           %d" % CC)
    print ("# beta              %f" % beta)
    print ("# box               %s" % LL)
    print ("# K                 %s" % KK)
    print ("# h                 %s" % (LL/KK))
    print ("# l_cut             %d" % l_cut)
    print ("# vanish bd         %s" % vanish_boundary)
    print ("# ")
    

    # make a water like system
    q2 = 33.456 * np.prod(LL) * (-0.8476 * -0.8476 + 0.4238 * 0.4238 * 2)
    natoms = 33.456 * np.prod(LL) * 3

    if vanish_boundary :
        lossfunc = HermiteLossFunc_Bound0 (CC, nbin, beta, KK, q2, natoms, region, l_cut = l_cut)
    else :
        lossfunc = HermiteLossFunc (CC, nbin, beta, KK, q2, natoms, region, l_cut = l_cut)
    
    
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
    if vanish_boundary :
        bv = np.delete (bv, -1)
        bd = np.delete (bd, -1)
    init_vv = bc.pack_v_d (bv, bd)
    tolerence = args.tolerence

    opt_res = sp.optimize.minimize (lossfunc.value, 
                                    init_vv, 
                                    jac = lossfunc.deriv, 
                                    method='BFGS', 
                                    options={'disp': True, 'gtol': 1e-3}, 
                                    tol=tolerence)

    [rv, rd] = bc.unpack_v_d (opt_res.x)
    if vanish_boundary: 
        pmat = bc.make_print_matrix_bound0 (CC, rv, rd)
    else :
        pmat = bc.make_print_matrix (CC, rv, rd)
    np.savetxt (args.output, pmat, header="%f %f %e"%(beta, LL[0]/KK[0], opt_res.fun))

if __name__ == '__main__':
    _main()
