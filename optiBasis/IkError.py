#!/usr/bin/env python3

import sys
import numpy as np
import scipy as sp
from Region import Region
from Bspline import Bspline
from HatComput import HatComput
from HatComput import HermiteBasisHatComput
from HatComput import HermiteBasisHatComput_Norm1
from HatComput import HermiteBasisHatComput_Norm1_Bound0
from IntplBasis import IntplBasis
from scipy.interpolate import interp1d
from basis_common import unpack_v_d
from basis_common import make_print_matrix
from basis_common import make_print_matrix_bound0

class IkError (object) :
    def __init__ (self, 
                  beta,
                  KK, 
                  hat_comput,
                  l_cut = 1,
    ) :
        self.conversion = 138.93545756169981341199
        self.beta = beta
        self.KK = KK
        self.hat_comput = hat_comput
        self.l_cut = l_cut
        assert len(self.KK) == 3
        
    def prepare_sum (self,
                     region) :
        VV = region.volume()
        sum_o1 = 0
        bd0 = int(self.KK[0]/2)
        bd1 = int(self.KK[1]/2)
        bd2 = int(self.KK[2]/2)

        for m0 in range (-bd0, bd0+1) :
            for m1 in range (-bd1, bd1+1) :
                for m2 in range (-bd2, bd2+1) :
                    if (abs(m0) + abs(m1) + abs(m2) == 0) : 
                        continue
                    gm = self._compute_G (m0, m1, m2, region)
                    gm2 = np.dot(gm, gm)
                    o1e = 0
                    for ll in range (-self.l_cut, self.l_cut+1) : 
                        if (ll == 0) : 
                            continue
                        tmpz = self.hat_comput[0] (m0 + ll * self.KK[0]) / self.hat_comput[0] (m0)
                        o1e = o1e + tmpz * tmpz
                        tmpz = self.hat_comput[1] (m1 + ll * self.KK[1]) / self.hat_comput[1] (m1)
                        o1e = o1e + tmpz * tmpz
                        tmpz = self.hat_comput[2] (m2 + ll * self.KK[2]) / self.hat_comput[2] (m2)
                        o1e = o1e + tmpz * tmpz
                    sum_o1 = sum_o1 + 2. * gm2 * o1e

        sum_o1 = sum_o1 / (2. * np.pi * VV * 2. * np.pi * VV)
        self.sum_o1 = sum_o1

    def prepare_sum_deriv (self, region) :
        VV = region.volume()
        bd0 = int(self.KK[0]/2)
        bd1 = int(self.KK[1]/2)
        bd2 = int(self.KK[2]/2)
        numb_basis = self.hat_comput[0].basis_value(0).shape[0]
        assert (numb_basis == self.hat_comput[1].basis_value(0).shape[0])
        assert (numb_basis == self.hat_comput[2].basis_value(0).shape[0])
        sum_o1 = np.zeros (numb_basis)

        for m0 in range (-bd0, bd0+1) :
            for m1 in range (-bd1, bd1+1) :
                for m2 in range (-bd2, bd2+1) :
                    if (abs(m0) + abs(m1) + abs(m2) == 0) : 
                        continue
                    gm = self._compute_G (m0, m1, m2, region)
                    gm2 = np.dot(gm, gm)
                    o1e = np.zeros (numb_basis)
                    idxmm = [m0, m1, m2]
                    for ll in range (-self.l_cut, self.l_cut+1) : 
                        if (ll == 0) : 
                            continue
                        for dd in range (3) :
                            hat_phil = self.hat_comput[dd] (idxmm[dd] + ll * self.KK[dd])
                            hat_phi0 = self.hat_comput[dd] (idxmm[dd] )
                            hat_basisl = self.hat_comput[dd].basis_value (idxmm[dd] + ll * self.KK[dd])
                            hat_basis0 = self.hat_comput[dd].basis_value (idxmm[dd])
                            o1e = o1e \
                                  + 2 * hat_phil / (hat_phi0 * hat_phi0) * hat_basisl  \
                                  - 2 * hat_phil * hat_phil / (hat_phi0 * hat_phi0 * hat_phi0) * hat_basis0                            
                    sum_o1 = sum_o1 + 2. * gm2 * o1e

        sum_o1 = sum_o1 / (2. * np.pi * VV * 2. * np.pi * VV)
        self.sum_deriv = sum_o1
        

    def _compute_G (self, m0, m1, m2, region) :
        if (m0 == 0 and m1 == 0 and m2 == 0) : 
            return np.zeros (3)
        invbox = region.rec_box()
        mm = m0 * invbox[0] + m1 * invbox[1] + m2 * invbox[2]
        mm2 = np.dot (mm, mm)
        expp = np.exp (- np.pi * np.pi / (self.beta * self.beta) * mm2) / mm2
        return -4. * np.pi * mm * expp

    def estimate (self,
                  q2,
                  natoms,
                  region) :
        self.prepare_sum (region)
        return np.sqrt (self.sum_o1 * q2 * q2 / float(natoms)) * self.conversion        

    def estimate_deriv (self,
                        q2,
                        natoms,
                        region) :
        self.prepare_sum (region)
        self.prepare_sum_deriv (region)
        return self.sum_deriv * 0.5 / np.sqrt (self.sum_o1) * np.sqrt(q2 * q2 / float(natoms)) * self.conversion 

class LossFunc (object) :
    def __init__ (self,
                  CC,
                  xx, 
                  beta,
                  KK,
                  q2, 
                  natoms, 
                  region,
    ) :
        self.CC = CC
        self.xx = xx
        self.beta = beta
        self.KK = KK
        self.q2 = q2
        self.natoms = natoms
        self.region = region
        
    def __call__ (self,
                  vv) :        
        basis = IntplBasis (self.xx, vv, self.CC)
        hat_basis = HatComput (basis, self.KK, over_smpl = 200)
        err_basis = IkError (self.beta, [self.KK, self.KK, self.KK], hat_basis, l_cut = 5)
        error = err_basis.estimate (self.q2, self.natoms, self.region)
        np.savetxt ('basis.out', vv)
        print ("returned %e" % error)
        return error

class HermiteLossFunc (object) :
    def __init__ (self,
                  CC,
                  nbin, 
                  beta,
                  KK,
                  q2, 
                  natoms, 
                  region,
                  l_cut = 10
    ) :
        self.CC = CC
        self.nbin = nbin
        self.beta = beta
        self.KK = KK
        self.q2 = q2
        self.natoms = natoms
        self.region = region
        self.l_cut = l_cut
        over_sampling = 400 * (self.nbin / self.CC)
        tmp0 = HermiteBasisHatComput_Norm1 (self.CC, self.nbin, self.KK[0], over_smpl = over_sampling)
        tmp1 = HermiteBasisHatComput_Norm1 (self.CC, self.nbin, self.KK[1], over_smpl = over_sampling)
        tmp2 = HermiteBasisHatComput_Norm1 (self.CC, self.nbin, self.KK[2], over_smpl = over_sampling)
        self.hhc = [tmp0, tmp1, tmp2]
        self.err_basis = IkError (self.beta, self.KK, self.hhc, l_cut = self.l_cut)

    def value (self, vv) :
        assert (len(vv) == self.nbin * 2)
        [hv, hd] = unpack_v_d (vv)
        self.hhc[0].set_value (hv, hd)
        self.hhc[1].set_value (hv, hd)
        self.hhc[2].set_value (hv, hd)
        error = self.err_basis.estimate (self.q2, self.natoms, self.region)
        print ("# returned %e" % error, end = "   \r")
        print_matrix = make_print_matrix (self.CC, hv, hd)
        np.savetxt ("basis.step.out", print_matrix)
        return error

    def deriv (self, vv) :
        assert (len(vv) == self.nbin * 2)
        [hv, hd] = unpack_v_d (vv)
        self.hhc[0].set_value (hv, hd)
        self.hhc[1].set_value (hv, hd)
        self.hhc[2].set_value (hv, hd)
        deriv = self.err_basis.estimate_deriv (self.q2, self.natoms, self.region)
        # return deriv[1:]
        return deriv        

class HermiteLossFunc_Bound0 (object) :
    def __init__ (self,
                  CC,
                  nbin, 
                  beta,
                  KK,
                  q2, 
                  natoms, 
                  region,
                  l_cut = 10
    ) :
        self.CC = CC
        self.nbin = nbin
        self.beta = beta
        self.KK = KK
        self.q2 = q2
        self.natoms = natoms
        self.region = region
        self.l_cut = l_cut
        over_sampling = 400 * (self.nbin / self.CC)
        tmp0 = HermiteBasisHatComput_Norm1_Bound0 (self.CC, self.nbin, self.KK[0], over_smpl = over_sampling)
        tmp1 = HermiteBasisHatComput_Norm1_Bound0 (self.CC, self.nbin, self.KK[1], over_smpl = over_sampling)
        tmp2 = HermiteBasisHatComput_Norm1_Bound0 (self.CC, self.nbin, self.KK[2], over_smpl = over_sampling)
        self.hhc = [tmp0, tmp1, tmp2]
        self.err_basis = IkError (self.beta, self.KK, self.hhc, l_cut = self.l_cut)

    def value (self, vv) :
        assert (len(vv) == self.nbin * 2 - 2)
        [hv, hd] = unpack_v_d (vv)
        self.hhc[0].set_value (hv, hd)
        self.hhc[1].set_value (hv, hd)
        self.hhc[2].set_value (hv, hd)
        error = self.err_basis.estimate (self.q2, self.natoms, self.region)
        print ("# returned %e" % error, end = "   \r")
        print_matrix = make_print_matrix_bound0 (self.CC, hv, hd)
        np.savetxt ("basis.step.out", print_matrix)
        return error

    def deriv (self, vv) :
        assert (len(vv) == self.nbin * 2 - 2)
        [hv, hd] = unpack_v_d (vv)
        self.hhc[0].set_value (hv, hd)
        self.hhc[1].set_value (hv, hd)
        self.hhc[2].set_value (hv, hd)
        deriv = self.err_basis.estimate_deriv (self.q2, self.natoms, self.region)
        # return deriv[1:]
        return deriv        

if __name__ == "__main__" : 
    q2 = (-0.8476 * -0.8476 + 0.4238 * 0.4238 * 2) * 1728
    natoms = 3 * 1728
    KK = 32
    beta = 3.5
    region = Region ([3.72412,3.72412,3.72412])
    bs = Bspline(4)

    # hat_cmpt = HatComput (bs, KK, over_smpl = 50)
    # err = IkError (beta, [KK, KK, KK], hat_cmpt, l_cut = 2)
    # error = err.estimate (q2, natoms, region)
    # CC = 2
    # bstep = 0.05
    # bx = np.arange (0, CC + bstep, bstep)
    # lossfunc = LossFunc (CC, bx, beta, KK, q2, natoms, region)
    # # bv = bs(bx)
    # bx0 = np.arange (0, CC + 0.1, 0.1)
    # bv0 = np.array([  6.65332244e-01,   6.56613040e-01,   6.30384087e-01,
    #                   5.91167006e-01,   5.40009415e-01,   4.80249323e-01,
    #                   4.15509569e-01,   3.48863540e-01,   2.83952069e-01,
    #                   2.21532255e-01,   1.67361006e-01,   1.21794403e-01,
    #                   8.52340358e-02,   5.62708564e-02,   3.48644458e-02,
    #                   1.97821619e-02,   9.75321961e-03,   3.77135561e-03,
    #                   6.34660481e-04,   2.68299590e-04,  -8.06167151e-04])
    # func0 = interp1d (bx0, bv0, axis = 0)
    # bv = func0 (bx)
    # print (lossfunc(bv))
    # aa = sp.optimize.minimize (lossfunc, bv, method='Nelder-Mead', options={'disp': True})
    # print (aa)
    # np.savetxt ('basis.out', aa.x)

    CC = 2
    nbin = 10
    bstep = CC / float(nbin)
    lossfunc = HermiteLossFunc (CC, nbin, beta, KK, q2, natoms, region)

    bx = np.arange (0, CC + .5 * bstep, bstep)
    bv_ = bs(bx)
    bd_ = bs.deriv (bx)
    scale = bv_[0]
    bv = bv_[1:]
    bd = bd_[1:]
    bv = bv / scale
    bd = bd / scale
    init_vv = np.append (bv, bd)
    print (lossfunc.value(init_vv))
    # aa = sp.optimize.minimize (lossfunc, init_vv, method='Nelder-Mead', options={'disp': True})

    aa = sp.optimize.minimize (lossfunc.value, init_vv, jac = lossfunc.deriv, method='BFGS', options={'disp': True}, tol=1e-4)
    
    
    
