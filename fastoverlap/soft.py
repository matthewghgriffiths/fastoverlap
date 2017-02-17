# -*- coding: utf-8 -*-
"""
Created on Fri Feb 17 12:45:15 2017

@author: mg542
"""

import numpy as np

from numpy import sin, cos, sqrt, pi
from scipy.special import eval_jacobi

try:
    from scipy.special.basic import factorial
except ImportError:
    from scipy.misc import factorial

class SOFT(object):
    """
    Class for calculating the 
    Special Orthogonal group (SO(3)) Fourier Transform (SOFT)
    This algorithm is based on the C++ library written by 
    
    http://www.cs.dartmouth.edu/~geelong/soft/
    """
    def __init__(self, bw):
        self.bw = bw
        self.bws = np.arange(0, bw*2)
        self.n = 2*bw
        self.Jmax = bw-1
        self.weights = self.makeweights(bw)
        self.a = pi/bw * self.bws
        self.b = pi/4/bw * (2*self.bws+1)
        self.y = pi/bw * self.bws
        self.Js, self.m1s, self.m2s = map(np.array, zip(*[(l,m1,m2) for l in xrange(self.Jmax+1) 
                                         for m1 in xrange(-l, l+1) 
                                         for m2 in xrange(-l, l+1)]))
        self.Ds = self.calcWignerMatrices()
        self.indFactor = np.array([2*pi/self.n, pi/self.n, 2*pi/self.n])
    ##
    @classmethod
    def makeweights(cls, bw):
        j = np.arange(0,2*bw).astype(float)[:,None]
        k = np.arange(0,bw).astype(float)[None,:]
        fudge = pi/4/bw
        weights = (2/(2*k+1) * 
                   sin((2*j+1)*(2*k+1)*fudge) * 
                   sin((2*j+1)*fudge)/bw).sum(1)
        return weights
    ##
    def calcWignerMatrices(self):
        bw = self.bw
        bws = self.bws
        beta = self.b
        sinb2 = sin(beta/2)
        cosb2 = cos(beta/2)
        cosb = cos(beta)
        Js, m1s, m2s = self.Js, self.m1s, self.m2s
        Dfull = np.zeros((bw, 2*bw-1, 2*bw-1, 2*bw))
        mu = abs(m1s-m2s)
        nu = abs(m1s+m2s)
        s = Js - (mu + nu)/2
        xi = np.ones_like(s)
        xi[m2s<m1s] = (-1)**(m2s-m1s)[m2s<m1s]
        factor = sqrt(factorial(s)*factorial(s+mu+nu)/
                      factorial(s+mu)/factorial(s+nu)) * xi
        jac = eval_jacobi(s[:,None], mu[:,None], nu[:,None],
                     cosb[None,:])
        Dfull[Js[:,None],m1s[:,None],m2s[:,None],
              bws[None,:]] = (factor[:,None] * jac *
                              sinb2[None,:] ** mu[:,None] * 
                              cosb2[None,:] ** nu[:,None])
        return Dfull
    ##
    def SOFT(self, data):
        """
        Compute forward transform
        """
        Jmax=self.Jmax
        bw=self.bw
        data = np.asanyarray(data)
        assert all(n==self.n for n in data.shape)
        S1 =  np.fft.fft(data, axis=0)
        S2 = np.fft.fft(S1, axis=2) * (2.*bw)**-2
        flmm = np.zeros((bw, bw*2 - 1, bw*2 - 1), np.complex128)
        for m1 in xrange(-Jmax,Jmax+1):
            for m2 in xrange(-Jmax,Jmax+1):
                l = max(np.abs([m1,m2]))
                flmm[l:,m1,m2] = self.Ds[l:,m1,m2].dot(self.weights*S2[m1,:,m2])
        return flmm
    ##
    def iSOFT(self, flmm):
        Jmax=self.Jmax
        bw=self.bw
        assert flmm.shape == (bw, bw*2 - 1, bw*2 - 1)
        S2out = np.zeros((2*bw,2*bw,2*bw), np.complex128)
        for m1 in xrange(-Jmax,Jmax+1):
            for m2 in xrange(-Jmax,Jmax+1):
                l = max(np.abs([m1,m2]))
                S2out[m1,:,m2] = self.Ds[l:,m1,m2].T.dot(flmm[l:,m1,m2])
        S1out = np.fft.ifft(S2out, axis=2)
        return np.fft.ifft(S1out, axis=0) * (2*bw)**2
    ##
    def indtoEuler(self, ind):
        R = self.indFactor * np.atleast_2d(ind)
        R[:,1] += 0.5 * pi / self.n
        return R.squeeze()
