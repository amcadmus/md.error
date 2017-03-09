#!/usr/bin/env python3

import os
import sys
import numpy as np
import argparse

global_version_str = "1.1.0"
global_basic_sys_size = 0.93103

def file_name (nbins, CC, beta, multi, KK) : 
    p_nbin = "n%03d" % nbins
    p_cc = "C%02d" % CC
    p_beta = "beta%.2f" % beta
    p_mul = "mul%03d" % multi
    p_kk = "K%03d" % KK
    input_file = "basis." +         \
                 p_nbin + "." +     \
                 p_cc + "." +       \
                 p_beta + "." +     \
                 p_mul + "." +     \
                 p_kk + ".out"    
    return input_file               

def _parse_argument():
    parser = argparse.ArgumentParser(
        description="*** Compress the basis files. ***")
    parser.add_argument('-n', '--bin-dens', type=int, default = 10,
                        help='number of bins divided by cut-off for piecewise cubic Hermite')
    parser.add_argument('-c', '--cut-off', type=int, nargs = '+', default = [2],
                        help='cut-off radius of basis')
    parser.add_argument('-b', '--beta', type=float, nargs = '+', default = [3.5],
                        help='the splitting parameter')
    parser.add_argument('-m', '--box-multiple', type=int, default = 4,
                        help='box size is the multiple of basic size %f' % global_basic_sys_size)
    parser.add_argument('-k', '--numb-grid', type=int, nargs = '+', default = [32],
                        help='number of grid points')
    parser.add_argument('-i', '--input-dir', type=str, default = "dir.basis",
                        help='the input dir of bases')
    parser.add_argument('-o', '--output', type=str, default = "basis.data",
                        help='the output file of basis database')
    args = parser.parse_args()
    return args

def _main () : 
    args = _parse_argument()
    
    bin_dens     = args.bin_dens
    CC           = args.cut_off
    beta         = args.beta
    mul          = args.box_multiple
    KK           = args.numb_grid
    input_dir    = args.input_dir
    output       = args.output

    CC = np.sort (CC)
    beta = np.sort (beta)
    tmpKK = np.sort(KK)
    for ii in range (len(tmpKK)) :
        KK[ii] = tmpKK[-(ii+1)]
#    KK = np.flip(np.sort(KK), 0)        
    print (KK)
    
    # loop order: CC, beta, hh
    ofp = open (output, 'w')

    # version
    ofp.write ("# tabulated PM interpl. basis version: %s\n" % global_version_str)

    # number of bins
    ofp.write ("# bin density\n")
    ofp.write ("%d" % bin_dens)
    ofp.write ("\n")

    # cut-off
    ofp.write ("# CC\n")
    for i_cc in CC:
        ofp.write ("%d" % i_cc)
    ofp.write ("\n")

    # beta
    ofp.write ("# beta\n")
    for i_beta in beta:
        ofp.write ("%.2f " % i_beta)
    ofp.write ("\n")

    # grid size
    ofp.write ("# hh\n")
    for i_kk in KK:
        ofp.write ("%.16e " % (global_basic_sys_size * mul/float(i_kk)))
    ofp.write ("\n")
    
    # opti value
    ofp.write ("# optimized value\n")
    for i_cc in CC:
        nbins = i_cc * bin_dens
        for i_beta in beta:
            for i_kk in KK :
                input_file = file_name (nbins, i_cc, i_beta, mul, i_kk)
                input_file = input_dir + "/" + input_file
#                print (input_file)
                fp = open (input_file, "r")
                value = float(fp.readline().split()[-1])
                ofp.write ("%e " % value)
    ofp.write ("\n")
                
    # opti value
    ofp.write ("# optimized value\n")
    for i_cc in CC:
        nbins = i_cc * bin_dens
        for i_beta in beta:
            for i_kk in KK :
                input_file = file_name (nbins, i_cc, i_beta, mul, i_kk)
                input_file = input_dir + "/" + input_file
                data = np.loadtxt (input_file)
                assert (data.shape == (nbins+1, 3))
                ofp.write ("# from %s\n" % input_file)
                for ii in range (data.shape[0]) :
                    for jj in range (data.shape[1]) :
                        ofp.write ("%19.16e " % data[ii][jj])
                    ofp.write ("\n")
                

if __name__ == '__main__':
    _main()
               

    
