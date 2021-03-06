#!/usr/bin/env python3

import os
import numpy as np
import subprocess as sp
import argparse
from compress_basis import file_name
from compress_basis import global_basic_sys_size

def _parse_argument():
    parser = argparse.ArgumentParser(
        description="*** Compute the optimal basis for SPME. Error estimated for water like homogeneous system. ***")
    parser.add_argument('-F', '--force-scheme', type=str, default = "ik",
                        help='the force scheme')
    parser.add_argument('-n', '--bin-dens', type=int, default = 10,
                        help='number of bins divided by cut-off for piecewise cubic Hermite')
    parser.add_argument('-c', '--cut-off', type=int, nargs = '+', default = [2],
                        help='cut-off radius of basis')
    parser.add_argument('-b', '--beta', type=float, nargs = '+', default = [1.0],
                        help='the splitting parameter')
    parser.add_argument('-m', '--box-multiple', type=int, default = 4,
                        help='box size is the multiple of basic size %f' % global_basic_sys_size)
    parser.add_argument('-k', '--numb-grid', type=int, nargs = '+', default = [32],
                        help='number of grid points')
    parser.add_argument('-T', '--numb-threads', type=int, default = 2,
                        help='number of threads')
    parser.add_argument('-o', '--output-dir', type=str, default = "dir.basis",
                        help='the output dir of basis')
    args = parser.parse_args()
    return args

def _main () : 
    args = _parse_argument()
    
    force_scheme = args.force_scheme
    bin_dens     = args.bin_dens
    CC           = args.cut_off
    beta         = args.beta
    mul          = args.box_multiple
    KK           = args.numb_grid
    output       = args.output_dir
    numb_threads = args.numb_threads
    if not os.path.exists (output): 
        os.mkdir (output)

    for i_cc in CC:
        nbins = i_cc * bin_dens
        for i_beta in beta:
            for i_kk in KK :
                LL = global_basic_sys_size * mul
                ofile = file_name (force_scheme, nbins, i_cc, i_beta, mul, i_kk)
                if os.path.exists (output+"/"+ofile) : 
                    print ("# existing file %s , do nothing" % (output+"/"+ofile))
                    continue
                command = "time c.src/compute_basis " +  \
                          " -F %s" % force_scheme + \
                          " -n %d" % nbins + \
                          " -c %d" % i_cc + \
                          " -b %.16e" % i_beta + \
                          " -l %.16e" % LL + \
                          " -k %d" % i_kk + \
                          " -T %d" % numb_threads + \
                          " -v" + \
                          " -t 1e-8" + \
                          " -o %s" % (output+"/"+ofile)
                print ("\n# run with command " + command)
                sp.check_call (command, shell = True)
                
if __name__ == '__main__':
    _main()

