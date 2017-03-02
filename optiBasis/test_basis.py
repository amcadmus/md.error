#!/usr/bin/env python3

import sys
import os
import numpy as np
import scipy as sp
import argparse
from Region import Region
from IkError import IkError
from IkError import HermiteLossFunc
from basis_common import make_dim3_vec
from basis_common import make_working_vec
from HatComput import HermiteBasisHatComput

def _parse_argument():
    parser = argparse.ArgumentParser(
        description="*** Error estimate of a basis. ***")
    parser.add_argument('input_file', type=str, 
                        help='the basis file')
    parser.add_argument('-b', '--beta', type=float, default = 3.5,
                        help='the splitting parameter')
    parser.add_argument('-k', '--numb-grid', type=int, nargs = '+', default = [32],
                        help='number of grid points')
    parser.add_argument('-l', '--box-size', type=float, nargs = '+', default = [3.72412],
                        help='box size')
    args = parser.parse_args()
    return args

def _main () :
    args = _parse_argument()

    if not os.path.isfile(args.input_file)  :
        raise RuntimeError ("no file " + args.input_file)
        
    box_size = make_dim3_vec (args.box_size)
    numb_grid = make_dim3_vec (args.numb_grid)

    KK = numb_grid
    LL = np.array(box_size)
    beta = args.beta
    region = Region (LL)    
    basis = np.loadtxt (args.input_file)
    [CC,vv,dd] = make_working_vec (basis)
    nbin = len(dd)
    over_sampling = 100 * (nbin / CC)
        
    # make a water like system
    q2 = 33.456 * np.prod(LL) * (-0.8476 * -0.8476 + 0.4238 * 0.4238 * 2)
    natoms = 33.456 * np.prod(LL) * 3
    
    tmp0 = HermiteBasisHatComput (CC, nbin, KK[0], over_smpl = over_sampling)
    tmp1 = HermiteBasisHatComput (CC, nbin, KK[1], over_smpl = over_sampling)
    tmp2 = HermiteBasisHatComput (CC, nbin, KK[2], over_smpl = over_sampling)
    hhc = [tmp0, tmp1, tmp2]
    esti = IkError (beta, KK, hhc, over_cmpt_ratio = 5)    

    hhc[0].set_value (vv,dd)
    hhc[1].set_value (vv,dd)
    hhc[2].set_value (vv,dd)
    error = esti.estimate (q2, natoms, region)
    print ("# estimated error is")
    print (error)


if __name__ == '__main__' :
    _main()
