#!/usr/bin/env python3

import sys
import numpy as np

class Region (object) :
    def __init__ (self, 
                  matrix = [1, 1, 1]
                  ) :
        self.mat = np.zeros ((3,3))
        if len(matrix) == 3 :
            for ii in range (3) :
                self.mat[ii][ii] = matrix[ii]
        elif matrix.shape == (3,3) :
            self.mat = matrix
        else :
            raise RuntimeError ("The region matrix should be of shape (3,) or (3,3)")
        self.invmat = np.linalg.inv (self.mat)

    def volume (self) :
        return np.linalg.det (self.mat)    
    
    def rec_box (self) :
        return self.invmat
