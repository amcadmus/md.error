#!/usr/bin/env python3

import sys
import numpy as np
import scipy as sp
import argparse
from Region import Region
from Bspline import Bspline
from HatComput import HermiteBasisHatComput_Norm1
from basis_common import unpack_v_d
from basis_common import make_print_matrix
from basis_common import make_dim3_vec
from basis_common import make_working_vec

def _parse_argument():
    parser = argparse.ArgumentParser(
        description="*** Print the Fourier transform of the basis. ***")
    parser.add_argument('input_basis', type=str, default = "basis.out",
                        help='the input file of basis')
    parser.add_argument('-k', '--numb-grid', type=int, default = 32,
                        help='the K value')
    parser.add_argument('-m', '--m-up', type=int, default = 32,
                        help='the max mode to be printed')
    args = parser.parse_args()
    return args

def _main () :
    args = _parse_argument()
        
    p_mat = np.loadtxt (args.input_basis)
    KK = args.numb_grid
    
    [CC, vv, dd] = make_working_vec (p_mat)
    nbins = len(dd)
    over_sampling = 400 * (nbins / CC)

    hat = HermiteBasisHatComput_Norm1 (CC, nbins, KK, over_sampling)
    hat.set_value (vv, dd)
    
    for ii in range (args.m_up) :
        print ("%d  %e" % (ii, hat(ii)))


if __name__ == '__main__':
    _main()
