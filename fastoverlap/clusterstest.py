# -*- coding: utf-8 -*-


from itertools import izip

import numpy as np
from numpy.linalg import norm
from numpy import sin, cos, sqrt, pi, exp
from scipy.special import eval_jacobi, gamma
from scipy.special.orthogonal import genlaguerre, eval_genlaguerre
from scipy.special import sph_harm, iv, eval_jacobi, gamma, hyp1f1, gammaln

from scipy.special import lpmv


import matplotlib.pyplot as plt


from utils import EulerM
from sphericalAlignment import SphericalHarmonicAlign, SphericalAlign

#import integrals_f as intf
#fast = intf.clusterfastoverlap

import fastclusters
commons = fastclusters.commons
clus = fastclusters.clusterfastoverlap
#from clusters import clusterfastoverlap as clus

import os
import csv
datafolder = "../examples/LJ38"
def readFile(filename):
    with open(filename, 'rb') as f:
        reader = csv.reader(f, delimiter=' ')
        dist = [map(float, row) for row in reader]
    return np.array(dist)

pos1 = readFile(os.path.join(datafolder, 'coords'))
pos2 = -readFile(os.path.join(datafolder, 'finish'))
p1, p2 = pos1.flatten(), pos2.flatten()

natoms = 38

max(norm(pos1,axis=1))

def rollaxes(a, axes=0):
    try:
        for axis in axes:
            a = np.roll(a, -(a.shape[axis]-1)/2, axis)
    except TypeError:
        a = np.roll(a, -(a.shape[axes]-1)/2, axes)
    return a


print 'test SO3 Coeffs'

n=30
l=14
r0 = 1.0
sigma = 0.3

c1nml = clus.harmoniccoeffsnml(pos1.flatten(),n,l,r0,sigma)
c2nml = clus.harmoniccoeffsnml(pos2.flatten(),n,l,r0,sigma)
#c1nlm = clus.harmoniccoeffsnlm(pos1.flatten(),n,l,r0,sigma)
#c2nlm = clus.harmoniccoeffsnlm(pos2.flatten(),n,l,r0,sigma)

#plt.semilogy(norm(c1nml,axis=1))
#plt.semilogy(norm(c2nml,axis=1))
#plt.semilogy(norm(c1nml,axis=1).T)
#plt.semilogy(norm(c2nml,axis=1).T)
Imml = clus.dotharmoniccoeffsnml(c1nml,c2nml)
norm(Imml)
#plt.semilogy(norm(c1nlm,axis=2))

sph = SphericalHarmonicAlign(sigma, r0, n, l)

Imml = clus.dotharmoniccoeffsnml(c1nml,c2nml)
Imml = clus.fouriercoeffsmml(p2, p1, l, sigma)
Ilmm0 = np.rollaxis(rollaxes(Imml, (0,1)),-1)
Ilmm = sph.calcSO3Coeffs(pos1, pos2)

print norm(Ilmm - Ilmm0)



self = sph
Ilmm = sph.calcSO3Coeffs(pos1, pos2)
Rs, amp, m, s, f = self.findRotations(Ilmm)

overlap, ilmm0 = clus.calcoverlap(Imml)

from utils import findPeaks
p,a,m,s,f = findPeaks(overlap, npeaks=10)

nrotations = np.array(10)
peaks, amplitude = fastclusters.fastoverlaputils.findpeaks(overlap, 10)

nrotations = np.array(10)
angles, rotms = clus.findrotations(overlap, nrotations)

p1, p2 = pos1.flatten(), pos2.flatten()

nrot = np.array(10)
angles, rotms, d, d2, rmat = clus.rotations(p1,p2,Imml,True,nrot)



commons.setlogic()
commons.initialise(natoms, np.array([np.arange(natoms)+1]), [1])
print 'aligncoeffs'

nrotations = np.array(10)
clus.aligncoeffs(p2, p1, Imml,1,nrotations)

#
debug=True
d, d2 = np.empty(2)
rot = np.empty((3,3))




#fastclusters.minpermdist(p2, p1, debug, 0,0,0,False,False,d,d2,False,rot)

#clus.aligncoeffs(p2, p1, Imml,True,10)
#fastclusters.findrotations()

#plt.semilogy((abs(Ilmm-Ilmm0).flatten()))
#
#plt.semilogy((norm(Ilmm-Ilmm0, axis=2)))
#plt.semilogy((norm(Ilmm-Ilmm0, axis=2)).T)
#
#plt.semilogy(np.sort(abs(Ilmm).flatten()))
#plt.semilogy(np.sort(abs(Imml).flatten()))
#
#plt.semilogy((abs(Ilmm).flatten()))
#plt.semilogy((abs(Ilmm0).flatten()))



#plt.semilogy((abs(Ilmm0).flatten()))

#Ilmm /= Ilmm[0,0,0]
#Ilmm0 /= Ilmm0[0,0,0]



#n=31
#l=31
#r0 = 1.0
#sigma = 0.3
#sph = SphericalHarmonicAlign(sigma, r0, n, l)
#
#pos1 = (np.random.random((natoms,3))-0.5)*3
#pos2 = (np.random.random((natoms,3))-0.5)*3
#
#c1nlm = clus.harmoniccoeffs(pos1.flatten(),n,l,r0,sigma)
#c2nlm = clus.harmoniccoeffs(pos2.flatten(),n,l,r0,sigma)
#Ilmm = clus.dotharmoniccoeffs(c1nlm,c2nlm)
#
#c1nml = clus.harmoniccoeffsnml(pos1.flatten(),n,l,r0,sigma)
#c2nml = clus.harmoniccoeffsnml(pos2.flatten(),n,l,r0,sigma)
#Imml = clus.dotharmoniccoeffsnml(c1nml,c2nml)
#Ilmm0 = np.rollaxis(rollaxes(Imml, (0,1)),-1)
#print norm(sph.calcSO3Coeffs(pos1, pos2) - clus.dotharmoniccoeffs(c1nlm,c2nlm))



#print norm(sph.calcSO3Coeffs(pos1, pos2) - fast.fouriercoeffs(pos2.flatten(), pos1.flatten(), l,sigma))
#Ilmm = clus.fouriercoeffs(pos2.flatten(), pos1.flatten(), l,sigma)

#for i in xrange(30):
#    r1, ylm = clus.rylm(pos1[i], 5)
#    r2, yml = clus.ryml(pos1[i], 5)
#    print norm(rollaxes(yml).swapaxes(0,1) - ylm),
#
#
#c1nlm = clus.harmoniccoeffs(pos1.flatten(),n,l,r0,sigma)
#c1nml = clus.harmoniccoeffsnml(pos1.flatten(),n,l,r0,sigma)
#c1nlm0 = rollaxes(c1nml,(1,)).swapaxes(1,2)
#print norm(c1nlm-c1nlm0)
#
#
#c1nml = clus.harmoniccoeffsnml(pos1.flatten(),n,l,r0,sigma)
#c2nml = clus.harmoniccoeffsnml(pos2.flatten(),n,l,r0,sigma)
#Imml = clus.dotharmoniccoeffsnml(c1nml,c2nml)
#
#Ilmm0 = np.rollaxis(rollaxes(Imml, (0,1)),-1)
#
#print norm(sph.calcSO3Coeffs(pos1, pos2)-Ilmm0)
#
#Imml= clus.fouriercoeffsmml(pos2.flatten(), pos1.flatten(), l, sigma)
#Imml, ilf, yml2, yml1 = clus.fouriercoeffsmem(pos2.flatten(), pos1.flatten(), l, sigma)
#
#Ilmm0 = np.rollaxis(rollaxes(Imml, (0,1)),-1)
#print norm(sph.calcSO3Coeffs(pos1, pos2)-Ilmm0)
#
#overlap, ilmm = clus.calcoverlap(Imml)
#
#
#import clusters
#
#fout = clusters.dsoft.isoft(Ilmm0)
#
#
#
#self=sph
#from utils import calcThetaPhiR
#
#pos1 = np.atleast_2d(pos1)
#pos2 = np.atleast_2d(pos2)
#assert pos1.shape == pos2.shape
#theta1, phi1, r1 = calcThetaPhiR(pos1)
#Y1 = self.sphHarm(theta1, phi1)
###
#theta2, phi2, r2 = calcThetaPhiR(pos2)
#Y2 = self.sphHarm(theta2, phi2).conj()
#r1r2 = r1[:,None]*r2[None,:]/self.scale**2/2
#il = iv(np.arange(self.Jmax+1)[:,None,None] + 0.5, r1r2[None,...])
#il *= np.sqrt(pi/2/r1r2[None,...])
#il *= exp(-(r1[None,:,None]**2+r2[None,None,:]**2)/4/self.scale**2)
#Ilmm = np.einsum("ijk,ilj,imk->ilm",il,Y1,Y2) * 4 * np.pi**2.5 * self.scale**3
#
#
#norm(rollaxes(yml2).swapaxes(0,1)-Y1)
#
#commons.setlogic()
#commons.initialise(natoms, np.array([np.arange(natoms)+1]), [1])
#
#p1 = pos1.flatten()
#p2 = pos2.flatten()
#
#debug=True
#clus.align(p2, p1, debug, 10)


#clus.rylm(pos1[26],l)[1]
#theta2, phi2, r2 = calcThetaPhiR(pos1[26])
#Y2 = self.sphHarm(theta2, phi2).conj()
#timeit clus.fouriercoeffs(pos2.flatten(), pos1.flatten(), l,sigma)


