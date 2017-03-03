#!/usr/bin/env python3

import numpy as np

def make_dim3_vec (args_box_size) :
    if len(args_box_size) == 1 : 
        box_size = [args_box_size[0], args_box_size[0], args_box_size[0]]
    else :
        box_size = args_box_size
    assert (len(box_size) == 3)
    return box_size

def unpack_v_d (vv) :
    nbin = int(len(vv) / 2)
    assert (len(vv) == nbin * 2)
    hv = vv[0:nbin]
    # hv = np.insert (hv, 0, 1.)
    hd = vv[nbin:nbin * 2]
    return [hv, hd]

def pack_v_d (hv, hd) :
    return np.append (hv, hd)

def make_print_matrix (CC, vv, dd) :
    nbin = len(dd)
    assert (nbin == len(vv))
    vv1 = np.insert (vv, 0, 1)
    dd1 = np.insert (dd, 0, 0)
    hh = CC / float(nbin)
    print_matrix = [np.arange (0, CC + hh/2., hh)]
    print_matrix = np.append (print_matrix, [vv1], axis = 0)
    print_matrix = np.append (print_matrix, [dd1],axis = 0)
    return print_matrix.T

def make_print_matrix_bound0 (CC, vv, dd) :
    nbin = len(dd) + 1
    assert (nbin == len(vv) + 1)
    vv1 = np.insert (vv, 0, 1)
    dd1 = np.insert (dd, 0, 0)
    vv1 = np.append (vv1, 0)
    dd1 = np.append (dd1, 0)
    hh = CC / float(nbin)
    print_matrix = [np.arange (0, CC + hh/2., hh)]
    print_matrix = np.append (print_matrix, [vv1], axis = 0)
    print_matrix = np.append (print_matrix, [dd1], axis = 0)
    return print_matrix.T

def make_working_vec (print_matrix) :
    print_matrix = print_matrix.T
    xx = print_matrix[0]
    nbin = len(xx) - 1
    hh = xx[1] - xx[0]
    CC = int((xx[-1] + 1e-6))
    vv1 = print_matrix[1]
    dd1 = print_matrix[2]
    vv = np.delete (vv1, 0)
    dd = np.delete (dd1, 0)
    return [CC, vv, dd]

def estimate_l_cut (CC, nbins, vanish_boundary) :
    esti = 0
    if vanish_boundary :
        esti = int (1 * nbins / CC)
    else :
        esti = int (4 * nbins / CC)
    if esti < 2 :
        esti = 2
    return esti
