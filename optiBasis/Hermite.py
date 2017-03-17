#!/usr/bin/env python3

import sys
import numpy as np

def _h00 (xx) : 
    return (1 + 2 * xx) * (1 - xx) * (1 - xx)

def _h10 (xx) : 
    return xx * (1 - xx) * (1 - xx)

def _h01 (xx) :
    return xx * xx * (3. - 2. * xx)

def _h11 (xx) :
    return xx * xx * (xx - 1)

def _h00p (xx) : 
    return 6 * xx * (xx - 1)

def _h10p (xx) : 
    return (3 * xx - 1) * (xx - 1)

def _h01p (xx) :
    return 6 * xx * (1 - xx)

def _h11p (xx) :
    return xx * (3 * xx - 2)

class HermiteV (object) :
    def __init__ (self,
                  knot, 
                  hh) :
        self.knot = knot
        self.hh = hh
        
    def _ele (self, 
              xx) :
        if xx < self.knot * self.hh :
            x1 = xx / self.hh - float(self.knot - 1)
            return _h01 (x1)
        else :
            x1 = xx / self.hh - float(self.knot)
            return _h00 (x1)

    def _ele_deriv (self, 
                    xx) :
        if xx < self.knot * self.hh :
            x1 = xx / self.hh - float(self.knot - 1)
            return _h01p (x1) / self.hh
        else :
            x1 = xx / self.hh - float(self.knot)
            return _h00p (x1) / self.hh

    def __call__ (self,
                  xx) :
        result = np.zeros (len(xx))
        for ii in range (len(xx)) :
            if xx[ii] / self.hh > (self.knot - 1) and xx[ii] / self.hh < (self.knot + 1) :
                result[ii] = self._ele (xx[ii])
        return result                    

    def deriv (self,
               xx) :
        result = np.zeros (len(xx))
        for ii in range (len(xx)) :
            if xx[ii] / self.hh > (self.knot - 1) and xx[ii] / self.hh < (self.knot + 1) :
                result[ii] = self._ele_deriv (xx[ii])
        return result                    
    
class HermiteD (object) :
    def __init__ (self,
                  knot, 
                  hh) :
        self.knot = knot
        self.hh = hh
        
    def _ele (self, 
              xx) :
        if xx < self.knot * self.hh :
            x1 = xx / self.hh - float(self.knot - 1)
            return _h11 (x1) * self.hh
        else :
            x1 = xx / self.hh - float(self.knot)
            return _h10 (x1) * self.hh

    def _ele_deriv (self, 
                    xx) :
        if xx < self.knot * self.hh :
            x1 = xx / self.hh - float(self.knot - 1)
            return _h11p (x1)
        else :
            x1 = xx / self.hh - float(self.knot)
            return _h10p (x1)

    def __call__ (self,
                  xx) :
        result = np.zeros (len(xx))
        for ii in range (len(xx)) :
            if xx[ii] / self.hh > (self.knot - 1) and xx[ii] / self.hh < (self.knot + 1) :
                result[ii] = self._ele (xx[ii])
        return result
        
    def deriv (self,
               xx) :
        result = np.zeros (len(xx))
        for ii in range (len(xx)) :
            if xx[ii] / self.hh > (self.knot - 1) and xx[ii] / self.hh < (self.knot + 1) :
                result[ii] = self._ele_deriv (xx[ii])
        return result
        
def symm_hermite (xx,
                  knot, 
                  hh) : 
    ii = knot
    if ii == 0 :
        her_v = HermiteV (ii, hh)
        her_d = HermiteD (ii, hh)
        hv = her_v (xx)
        hd = her_d (xx)
    else :
        her_v = HermiteV (ii, hh)
        her_d = HermiteD (ii, hh)
        hv = her_v (xx)
        hd = her_d (xx)
        her_v = HermiteV (-ii, hh)
        her_d = HermiteD (-ii, hh)
        hv = hv + her_v (xx)
        hd = hd - her_d (xx)            
    return [hv, hd]
    
def symm_hermite_deriv (xx,
                        knot, 
                        hh) : 
    ii = knot
    if ii == 0 :
        her_v = HermiteV (ii, hh)
        her_d = HermiteD (ii, hh)
        hv = her_v.deriv (xx)
        hd = her_d.deriv (xx)
    else :
        her_v = HermiteV (ii, hh)
        her_d = HermiteD (ii, hh)
        hv = her_v.deriv (xx)
        hd = her_d.deriv (xx)
        her_v = HermiteV (-ii, hh)
        her_d = HermiteD (-ii, hh)
        hv = hv + her_v.deriv (xx)
        hd = hd - her_d.deriv (xx)            
    return [hv, hd]
    
