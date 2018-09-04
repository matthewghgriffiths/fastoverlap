# -*- coding: utf-8 -*-

import numpy as np
from itertools import izip

from scipy import (sin, cos, sqrt, pi, exp, arange, zeros, float32)
from numpy.linalg import norm
from scipy.special import sph_harm, gamma, jn_zeros, jn#, spherical_jn
try:
    from scipy.special.basic import factorial
except ImportError:
    from scipy.misc import factorial
from scipy.optimize import minimize, brentq

from utils import find_best_permutation, EulerM, coeffs_harmonicBasis,\
findMax, findPeaks, findrotation, eval_grad_jacobi, calcThetaPhiR
    
from soft import SOFT

from sphericalAlignment import BaseSphericalAlignment, SphericalAlign,\
    SphericalHarmonicAlign

def spherical_jn(n, z):
    return sqrt(pi/z/2) * jn(n+0.5,z)

def spherical_jn_zeros(n,nt):
  zerosj = zeros((n+1, nt), dtype=float32)
  zerosj[0] = arange(1,nt+1)*pi
  points = arange(1,nt+n+1)*pi
  racines = zeros(nt+n, dtype=float32)
  Jn = lambda r, n: spherical_jn(n, r)
  for i in xrange(1,n+1):
    for j in xrange(nt+n-i):
      foo = brentq(Jn, points[j], points[j+1], (i,))
      racines[j] = foo
    points = racines
    zerosj[i][:nt] = racines[:nt]
  return (zerosj)

class SphericalFourierAlign(BaseSphericalAlignment):
    def __init__(self, scale, rcut=None, nmax=15, Jmax=15, perm=None):
        self.scale = scale
        self.nmax = nmax
        self.Jmax = Jmax
        self.rcut = rcut
        if self.rcut is not None:
            self.setCoeffs(nmax, Jmax, rcut)
        self.perm = perm
    
    def setCoeffs(self, nmax=None, Jmax=None, rcut=None):
        
        if nmax is None:
            nmax = self.nmax
        else:
            self.nmax = nmax
            
        if Jmax is None:
            Jmax = self.Jmax
        else:
            self.setJ(Jmax)
            
        if rcut is None:
            rcut = self.rcut
        else:
            self.rcut = rcut
        
        # zeros of the spherical bessel function
        self.kln = spherical_jn_zeros(self.Jmax, self.nmax)
        self.jn1 = spherical_jn(np.arange(1,self.nmax+1)[None,:], self.kln)
        self.normf = (2*self.kln*self.rcut**-3*self.jn1**-2)**0.5
        self.expksr = exp(- self.kln**2 * self.scale**2 / self.rcut**2 / 2)
        
    def calcFourierCoeffs(self, pos):
        
        theta, phi, r = calcThetaPhiR(pos)
        Y = self.sphHarm(theta, phi)
        
        Clmn = ((2*self.scale**2)**(3/2) * 
                self.kln[:,None,:] * self.expksr[:,None,:] * 
                (spherical_jn(
                        np.arange(self.Jmax + 1)[:,None,None], 
                        self.kln[:,:,None] * 
                        r[None,None,:] / self.rcut )[:,None,:,:] * 
                 Y[:,:,None,:]).sum(3) )
        return Clmn
        
    def calcSO3Coeffs(self, pos1, pos2):
        c1lmn = self.calcFourierCoeffs(pos1)
        c2lmn = self.calcFourierCoeffs(pos2)
        return np.einsum("lmn,lon->lmo", c1lmn, c2lmn.conj())

    def calcSO3Harm(self, c1nlms, c2nlms, invert=False):
        if invert:
            c2nlms = (c2*(-1)**self.J[None,:,None] for c2 in c2nlms)
        return sum( np.einsum("nlm,nlo->lmo", c1, c2.conj())
                    for c1, c2 in izip(c1nlms, c2nlms))
    
    def alignGroup(self, coords, keepCoords=False):
        n = len(coords)
        if keepCoords:
            aligned = np.empty((2, n, n) + coords[0].shape)
        coeffs = [[self.calcHarmCoeff(p)] for p in coords]
        dists = np.zeros((n, n))
        for i, (pos1, c1) in enumerate(izip(coords, coeffs)):
            for j, (pos2, c2) in enumerate(izip(coords, coeffs)):
                calcCoeffs = lambda b: self.calcSO3Harm(c1, c2, b)
                dist, x1, x2 = self.align(pos1, pos2, calcCoeffs=calcCoeffs)[:3]
                if keepCoords:
                    aligned[0,i,j,...] = x1
                    aligned[1,i,j,...] = x2
                dists[i,j] = dist
        if keepCoords:
            return dists, aligned
        else:
            return dists
        
        

import os
import csv

import matplotlib as mpl
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

import seaborn as sns

plt.ion()

datafolder = "../examples/LJ38"
def readFile(filename):
    with open(filename, 'rb') as f:
        reader = csv.reader(f, delimiter=' ')
        dist = [map(float, row) for row in reader]
    return np.array(dist)

pos1 = readFile(os.path.join(datafolder, 'coords'))
pos2 = readFile(os.path.join(datafolder, 'finish'))

natoms = 38
scale = 0.3
Jmax = 14
Nmax = 20
harmscale = 1.0

soap = SphericalAlign(scale, Jmax)
harm = SphericalHarmonicAlign(scale, harmscale, Nmax, Jmax)

rcut = 4.0

fourier = SphericalFourierAlign(scale, rcut, 30, Jmax)
self = fourier

rs = np.linspace(0,5,256)

#for l in xrange(self.Jmax):
#    for n in xrange(self.nmax):
#        plt.plot(rs, spherical_jn(l, self.kln[l,n] * rs / rcut))

Ilmm = soap.calcSO3Coeffs(pos1, pos2)
hIlmm = harm.calcSO3Coeffs(pos1, pos2)

c1lmn = self.calcFourierCoeffs(pos1)
c2lmn = self.calcFourierCoeffs(pos2)
fIlmm = np.einsum("lmn,lon->lmo", c1lmn, c2lmn.conj())

fIlmm *= norm(Ilmm)/norm(fIlmm)

print norm(Ilmm-hIlmm)

