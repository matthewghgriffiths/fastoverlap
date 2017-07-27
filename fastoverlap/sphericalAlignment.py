# -*- coding: utf-8 -*-
"""
Created on Sun Feb 12 16:15:14 2017

@author: mg542
"""

import numpy as np
from itertools import izip

from numpy import sin, cos, sqrt, pi, exp
from numpy.linalg import norm
from scipy.special import sph_harm, iv, eval_jacobi, gamma, hyp1f1
try:
    from scipy.special.basic import factorial
except ImportError:
    from scipy.misc import factorial
from scipy.optimize import minimize

from utils import find_best_permutation, EulerM, coeffs_harmonicBasis,\
findMax, findPeaks, findrotation, eval_grad_jacobi, calcThetaPhiR
    
from soft import SOFT

import f90
if f90.have_fastclusters:
    fastclusters = f90.fastclusters

class BaseSphericalAlignment(object):
    def calcSO3Coeffs(self, pos1, pos2):
        raise NotImplementedError
    ##
    def setJ(self, Jmax):
        self.Jmax = Jmax
        self.J = np.arange(Jmax+1)
        ## Needed for Spherical Harmonic Calcualtion
        self.Jinds = np.array([j for j in xrange(Jmax+1) for _ in xrange(-j, j+1)])
        self.Minds = np.array([m for j in xrange(Jmax+1) for m in xrange(-j, j+1)])
        self.soft = SOFT(Jmax+1)
        ## Needed for Wigner Matrix Calculation
        self.Js, self.m1s, self.m2s = map(np.array, 
                                          zip(*[(l,m1,m2) for l in xrange(Jmax+1) 
                                                for m1 in xrange(-l, l+1) 
                                                for m2 in xrange(-l, l+1)]))
    ##
    def sphHarm(self, theta, phi):
        theta = np.atleast_1d(theta)
        phi = np.atleast_1d(phi)
        assert theta.shape == phi.shape
        Y = np.zeros((self.Jmax+1, 2*self.Jmax+1)+theta.shape, np.complex128)
        Y[self.Jinds,self.Minds] = sph_harm(self.Minds[:,None], 
                                            self.Jinds[:,None],
                                            phi[None,:], theta[None,:])
        return Y
    ##
    def calcWignerMatrices(self, rot):
        a,b,y = rot
        sinb2 = sin(b/2)
        cosb2 = cos(b/2)
        cosb = cos(b)
        sinb = sin(b)
        Jmax = self.Jmax
        Js, m1s, m2s = self.Js, self.m1s, self.m2s
        Ds = np.zeros((Jmax+1, 2*Jmax+1, 2*Jmax+1), np.complex128)
        grad = np.zeros((3, Jmax+1, 2*Jmax+1, 2*Jmax+1), np.complex128)
        mu = abs(m1s-m2s)
        nu = abs(m1s+m2s)
        s = Js - (mu + nu)/2
        xi = np.ones_like(s)
        xi[m2s<m1s] = (-1)**(m1s-m2s)[m2s<m1s]
        factor = sqrt(factorial(s)*factorial(s+mu+nu)/
                          factorial(s+mu)/factorial(s+nu)) * xi
        jac = eval_jacobi(s, mu, nu, cosb)
        d = (factor * jac * sinb2 ** mu * cosb2 ** nu)
        ##
        gradjac = eval_grad_jacobi(s, mu, nu, cosb) * - sinb
        gradd = (factor * gradjac * sinb2 ** mu * cosb2 ** nu +
                     factor * jac * sinb2 ** (mu-1) * cosb2 ** (nu+1) * mu/2 - 
                     factor * jac * sinb2 ** (mu+1) * cosb2 ** (nu-1) * nu/2)
        ##
        Ds[Js, m1s, m2s] = exp(-1j*m1s*a) * d * exp(-1j*m2s*y)
        grad[0,Js, m1s, m2s] = -1j*m1s * Ds[Js, m1s, m2s]
        grad[1,Js, m1s, m2s] = exp(-1j*m1s*a) * gradd * exp(-1j*m2s*y)
        grad[2,Js, m1s, m2s] = -1j*m2s * Ds[Js, m1s, m2s]
        return Ds, grad
    ##
    def getEnergyGradient(self, rot, Ilmm):
        D, gradD = self.calcWignerMatrices(rot)
        return -(Ilmm*D).real.sum(), -(Ilmm[None,...]*gradD).real.sum((1,2,3))
    ##
    def Hungarian(self, pos1, pos2, perm=None):
        if perm is None:
            if self.perm is None:
                perm = [np.arange(len(pos1))]
            else:
                perm = self.perm
        return find_best_permutation(pos1, pos2, permlist=perm)
    ##
    def maxOverlap(self, R, Ilmm):
        res = minimize(self.getEnergyGradient, R, jac=True, args=(Ilmm,), 
                       method='L-BFGS-B')
        return res.x, res
    ##
    def rotate(self, X, R):
        return X.dot(EulerM(*R))
    ##
    def refine(self, X1, X2, R, permlist=None):
        if permlist is None:
            if self.perm is None:
                permlist = [np.arange(len(X1))]
            else:
                permlist = self.perm         
        X2 = np.dot(X2, EulerM(*R))
        _, perm = self.Hungarian(X1, X2, permlist)
        dist, MR = findrotation(X1, X2[perm])
        return dist, X1, X2[perm].dot(MR.T)
    ##
    def COM_shift(self, pos1, pos2):
        """
        Center the centroid of the particle positions at the origin
        """
        X1 = np.array(pos1, float)
        X2 = np.array(pos2, float)
        X1 -= X1.mean(axis=0)[None,:]
        X2 -= X2.mean(axis=0)[None,:]
        return X1, X2
    ##
    def _align(self, pos1, pos2, perm=None, invert=True):
        X1, X2 = self.COM_shift(pos1, pos2)
        if perm is None:
            if self.perm is None:
                perm = [np.arange(len(pos1))]
            else:
                perm = self.perm            
        Ilmm = sum(self.calcSO3Coeffs(X1[p], X2[p]) for p in perm)
        R, res = self.findRotation(Ilmm)
        if invert:
            Ilmm = sum(self.calcSO3Coeffs(X1[p], -X2[p]) for p in perm)
            invR, invres = self.findRotation(Ilmm)
            if invres < res:
                R = invR
                X2 *= -1
        return self.refine(X1, X2, R, perm)
    ##
    def align(self, pos1, pos2, perm=None, invert=True, calcCoeffs=None):
        X1, X2 = self.COM_shift(pos1, pos2)
        if perm is None:
            if self.perm is None:
                perm = [np.arange(len(pos1))]
            else:
                perm = self.perm
        if calcCoeffs is None:
            Ilmm = sum(self.calcSO3Coeffs(X1[p], X2[p]) for p in perm)
        else:
            Ilmm = calcCoeffs(False)
        R, res = self.findRotation(Ilmm)
        if invert:
            if calcCoeffs is None:
                Ilmm = sum(self.calcSO3Coeffs(X1[p], -X2[p]) for p in perm)
            else:
                Ilmm = calcCoeffs(True)
            invR, invres = self.findRotation(Ilmm)
            if invres < res:
                R = invR
                X2 *= -1
        return self.refine(X1, X2, R, perm)
    ##
    def findRotation(self, Ilmm):
        overlap = self.soft.iSOFT(Ilmm)
        R = self.soft.indtoEuler(findMax(overlap))
        R, res = self.maxOverlap(R, Ilmm.conj())
        return R, res.fun
    ##
    def findRotations(self, Ilmm, nrot=10, width=2):
        overlap = self.soft.iSOFT(Ilmm).real
        peaks = []
        while len(peaks)==0:
            peaks, amplitude, mean, sigma, f = findPeaks(overlap, 
                                                         npeaks=nrot, 
                                                         width=width)
            width += 1
        return np.atleast_2d(self.soft.indtoEuler(peaks)), amplitude, mean, sigma, f
    ##
    def malign(self, pos1, pos2, perm=None, invert=True, calcCoeffs=None, nrot=10):
        X1, X2 = self.COM_shift(pos1, pos2)
        if perm is None:
            if self.perm is None:
                perm = [np.arange(len(pos1))]
            else:
                perm = self.perm
        if calcCoeffs is None:
            Ilmm = sum(self.calcSO3Coeffs(X1[p], X2[p]) for p in perm)
        else:
            Ilmm = calcCoeffs(False)
        Rs = self.findRotations(Ilmm, nrot)[0]
        dist, rX1, rX2 = min((self.refine(X1, X2, R, perm) for R in Rs),
                             key=lambda x: x[0])
        if invert:
            if calcCoeffs is None:
                Ilmm = sum(self.calcSO3Coeffs(X1[p], -X2[p]) for p in perm)
            else:
                Ilmm = calcCoeffs(True)
            Rs = self.findRotations(Ilmm, nrot)[0]
            idist, iX1, iX2 = min((self.refine(X1, X2, R, perm) for R in Rs),
                                  key=lambda x: x[0])
            if idist<dist:
                return idist, iX1, iX2
        return dist, rX1, rX2
    ##
    def __call__(self, pos1, pos2, perm=None, invert=True, calcCoeffs=None, nrot=10):
        dist, X1, X2 = self.align(pos1, pos2, perm, invert, calcCoeffs)
        if (norm(X1-X2,axis=1)>self.scale).sum() > len(X1)/3:
            try:
                mdist, mX1, mX2 = self.malign(pos1, pos2, perm, invert, calcCoeffs, nrot)
                if mdist < dist:
                    return mdist, mX1, mX2
            except:
                pass
        return dist, X1, X2
    

class SphericalAlign(BaseSphericalAlignment):
    def __init__(self, scale, Jmax=15, perm=None):
        self.scale = scale
        self.setJ(Jmax)
        self.perm = perm
    ##
    def calcSO3Coeffs(self, pos1, pos2):
        pos1 = np.atleast_2d(pos1)
        pos2 = np.atleast_2d(pos2)
        assert pos1.shape == pos2.shape
        theta1, phi1, r1 = calcThetaPhiR(pos1)
        Y1 = self.sphHarm(theta1, phi1)
        ##
        theta2, phi2, r2 = calcThetaPhiR(pos2)
        Y2 = self.sphHarm(theta2, phi2).conj()
        r1r2 = r1[:,None]*r2[None,:]/self.scale**2/2
        il = iv(np.arange(self.Jmax+1)[:,None,None] + 0.5, r1r2[None,...])
        il *= np.sqrt(pi/2/r1r2[None,...])
        il *= exp(-(r1[None,:,None]**2+r2[None,None,:]**2)/4/self.scale**2)
        return np.einsum("ijk,ilj,imk->ilm",il,Y1,Y2) * 4 * np.pi**2.5 * self.scale**3
       
     
class SphericalHarmonicAlign(BaseSphericalAlignment):
    def __init__(self, scale, harmscale=1.0, nmax=15, Jmax=15, perm=None):
        self.scale = scale
        self.harmscale = harmscale
        self.setCoeffs(nmax, Jmax, harmscale)
        self.perm = perm
    ##
    @classmethod
    def radialIntegralHarmonic(cls, n,l,rj,sigma,r0):
        """
        Returns the result of the integral,
        \int_0^{\infty} \exp{\left(-\frac{r^2+{r_j}^2}{2\sigma^2}\right)} 
                        \exp{\left(-\frac{r^2}{2r_0^2}\right)}
                        i_l \left( \frac{r r_{j}}{\sigma^2} \right) r^n r^2
                        \; \mathrm{d}r
        """
        return ((2.**(3+n-l)*pi**3)**0.5 * rj**l * 
                (sigma**-2 + r0**-2)**(-0.5*(3+n+l)) * sigma**(-2*l) *
                hyp1f1(0.5*(3+n+l),1.5+l,0.5* rj**2 *r0**2/(r0**2 * sigma**2 +sigma**4)) * 
                gamma(0.5*(3 + l + n))/gamma(1.5+l)) * exp(-0.5*rj**2*sigma**-2)
    ##
    @classmethod
    def HarmInt(cls, n, l, r0):
        """
        Gives the result of the integral
        \int_0^{\infty}
        N_{nl} r^l \exp{\left(-\frac{r^2}{2r^2_0}\right)} 
        L^{l+1/2}_n\left(\frac{r^2}{r^2_0}\right)
        \; \mathrm{d}r
        """
        coeffs = coeffs_harmonicBasis(n, l, r0)
        ns = np.arange(2*n+l+1)
        intn =  gamma(0.5*(ns+3)) * 2**(0.5*(ns+1)) * r0**(3+ns)
        return coeffs.dot(intn)
    ##
    @classmethod
    def HarmCoeffs(cls, nmax, lmax, r0):
        coeffs = np.zeros((nmax+1, lmax+1, 2*nmax+lmax+1))
        for n in xrange(nmax+1):
            for l in xrange(lmax+1):
                coeffs[n,l,:2*n+l+1] = coeffs_harmonicBasis(n,l,r0)
        return coeffs
    ##
    @classmethod
    def HarmInts(cls, nmax, lmax, r0):
        coeffs = cls.HarmCoeffs(nmax, lmax, r0)
        ns = np.arange(2*nmax+lmax+1)
        intn =  gamma(0.5*(ns+3)) * 2**(0.5*(ns+1)) * r0**(3+ns)
        return coeffs.dot(intn)
    ##
    def setCoeffs(self, nmax=None, Jmax=None, harmscale=None):
        if nmax is None:
            nmax = self.nmax
        else:
            self.nmax = nmax
        if Jmax is None:
            Jmax = self.Jmax
        else:
            self.setJ(Jmax)
        if harmscale is None:
            harmscale = self.harmscale
        else:
            self.harmscale = harmscale
        self.coeffs = self.HarmCoeffs(nmax, Jmax, harmscale)
    ##
    def calcHarmCoeff(self, pos):
        """
        Returns the harmonic basis coefficients for a set of gaussians
        at positions specified by pos, with widths set by self.scale and the 
        lengthscale of the harmonic basis set by self.harmonicscale
        """ 
        theta, phi, r = calcThetaPhiR(pos)
        Y = self.sphHarm(theta, phi)
        dslj = self.radialIntegralHarmonic(np.arange(2*self.nmax+
                                                     self.Jmax+1)[:,None,None],
                                           self.J[None,:,None],
                                           r[None,None,:],
                                           self.scale,
                                           self.harmscale)
        dnlj = np.einsum("nls,slj->nlj",self.coeffs, dslj)
        cnlm = np.einsum("nlj,lmj->nlm", dnlj, Y.conj())
        return cnlm
    ##
    def calcSO3Coeffs(self, pos1, pos2):
        c1nlm = self.calcHarmCoeff(pos1)
        c2nlm = self.calcHarmCoeff(pos2)
        return np.einsum("nlm,nlo->lmo", c1nlm.conj(), c2nlm)
    ##
    def calcSO3Harm(self, c1nlms, c2nlms, invert=False):
        if invert:
            c2nlms = (c2*(-1)**self.J[None,:,None] for c2 in c2nlms)
        return sum( np.einsum("nlm,nlo->lmo", c1.conj(), c2)
                    for c1, c2 in izip(c1nlms, c2nlms))
    ##
    def _align(self, pos1, pos2, perm=None, invert=True, cnlms=None):
        X1, X2 = self.COM_shift(pos1, pos2)
        if perm is None:
            if self.perm is None:
                perm = [np.arange(len(pos1))]
            else:
                perm = self.perm
        if cnlms is None:
            c1nlms = [self.calcHarmCoeff(X1[p]) for p in perm]
            c2nlms = [self.calcHarmCoeff(X2[p]) for p in perm]
        else:
            c1nlms, c2nlms = cnlms
        Ilmm = self.calcSO3Harm(c1nlms, c2nlms)
        R, res = self.findRotation(Ilmm)
        if invert:
            c2nlms = [c2*(-1)**self.J[None,:,None] for c2 in c2nlms]
            Ilmm = self.calcSO3Harm(c1nlms, c2nlms)
            invR, invres = self.findRotation(Ilmm)
            if invres < res:
                R = invR
                X2 *= -1
        return self.refine(X1, X2, R, perm)
    ##
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
        
class VarSphericalHarmonicAlign(SphericalHarmonicAlign):
    def __init__(self, scale, harmscale=1.0, nmax=15, Jmax=7, perm=None):
        self.scale = scale
        self.harmscale = harmscale
        self.setCoeffs(nmax, Jmax, harmscale)
        self.perm = perm
    ##
    def calcHarmCoeff(self, pos):
        """
        Returns the harmonic basis coefficients for a set of gaussians
        at positions specified by pos, with widths set by self.scale and the 
        lengthscale of the harmonic basis set by self.harmonicscale
        """ 
        theta, phi, r = calcThetaPhiR(pos)
        Y = self.sphHarm(theta, phi)
        scale = self.scale * r
        dslj = self.radialIntegralHarmonic(np.arange(2*self.nmax+
                                                     self.Jmax+1)[:,None,None],
                                           self.J[None,:,None],
                                           r[None,None,:],
                                           scale[None,None,:],
                                           self.harmscale) ** scale[None,None,:]**-0.5
        cnlm = np.einsum("nls,slj,lmj->nlm",self.coeffs, dslj, Y.conj())
        return cnlm
        
class SphericalAlignFortran(BaseSphericalAlignment):
    """ Class for aligning two isolated structures, wrapper for FORTRAN 
    subroutines to do the numerical heavy-lifting
    
    Parameters:
    ----------
    scale : float, optional
        The width of the Gaussian kernels.
    Jmax : int, optional
        Specifies the maximum angular momentum degree up to which SO(3) coefficients
        will be calculated
    permlist : sequence of arrays, optional
        Each array in perm represents a different permutation group
    Natoms : int, optional
        Number of atoms
        
    Notes:
    ---------
    
    Be very careful about setting perm/Natoms, as if this is not properly set 
    then you will get segfaults. Also be careful about using this object in
    tandem 
    """
    def __init__(self, scale=0.3, Jmax=15, perm=None, Natoms=None):
        self.scale=scale
        self.Jmax=Jmax
        self.fast = fastclusters
        self.fast.clusterfastoverlap.setcluster()
        self.Natoms = Natoms
        if perm is not None:
            self.setPerm(perm)
        elif Natoms is not None:
            self.setPerm(perm)
        else:
            self.perm = perm
        self.malign = self.__call__
    ##
    def setPerm(self, perm):
        self.Natoms = sum(map(len,perm))
        self.perm = perm
        self.nperm = len(perm)
        self.npermsize = map(len, perm)
        self.permgroup = np.concatenate([np.asanyarray(p)+1 for p in perm])
        self.fast.fastoverlaputils.setperm(self.Natoms, self.permgroup, self.npermsize)
    ##
#    def malign(self, pos1, pos2, perm=None, invert=False, nrot=10, debug=False):
#        return self(pos1, pos2, perm, invert, nrot, debug)
    ##
    def align(self, pos1, pos2, perm=None, invert=True, debug=False):
        return self(pos1, pos2, perm, invert, 1, debug)
    ##
    def __call__(self, pos1, pos2, perm=None, invert=True, nrot=10, debug=False):
        """ Aligns two isolated structures, aligns and permutes pos2 to match
        pos1 as closely as possible. Returns the distance between the aligned 
        structures, and the rotated and centred coordinates of the two structures
        
        Parameters
        ----------
        pos1 : (Natoms, 3) array_like
            Coordinate array.
        pos2 : (Natoms, 3) array_like
            Coordinate array to align with pos1.
        nrot : int, optional
            Number of different displacements to test
        perm : sequence of arrays, optional
            Each array in perm represents a different permutation group
        invert : bool
            If true tests inverted configurations as well
        debug : bool, optional
            If true prints debug information
        
        Returns
        -------
        distance : float
            Euclidean distance between X1 and X2
        X1 : (Natoms, 3) array
            pos1 centred on the origin
        X2 : (Natoms, 3) array
            Aligned coordinates of pos2
        rmatbest : (3,3) array
            the rotation matrix that maps centred pos2 onto X2
        """
        if perm is not None:
            self.setPerm(perm)
        elif len(pos1) != self.Natoms:
            self.Natoms = len(pos1)
            self.setPerm([np.arange(self.Natoms)])
            
        coordsb = np.asanyarray(pos1).flatten()
        coordsa = np.asanyarray(pos2).flatten()
        self.fast.clusterfastoverlap.setcluster()
        self.fast.commons.perminvopt = invert
        args = (coordsb,coordsa,debug,self.Jmax,self.scale,nrot)
        dist, _, rmatbest = self.fast.clusterfastoverlap.align(*args)
        return dist, coordsb.reshape(self.Natoms,3), coordsa.reshape(self.Natoms,3), rmatbest
        
class SphericalHarmonicAlignFortran(BaseSphericalAlignment):
    """ Class for aligning two isolated structures, wrapper for FORTRAN 
    subroutines to do the numerical heavy-lifting
    
    Parameters:
    ----------
    scale : float, optional
        The width of the Gaussian kernels.
    Jmax : int, optional
        Specifies the maximum angular momentum degree up to which SO(3) coefficients
        will be calculated
    permlist : sequence of arrays, optional
        Each array in perm represents a different permutation group
    Natoms : int, optional
        Number of atoms
    """
    def __init__(self, scale=0.3, Jmax=15, harmscale=1.0, nmax=20, perm=None, Natoms=None):
        self.scale=scale
        self.Jmax=Jmax
        self.harmscale=harmscale
        self.nmax=nmax
        self.fast = fastclusters
        self.clus = fastclusters.clusterfastoverlap
        self.Natoms = Natoms
        if perm is not None:
            self.setPerm(perm)
        elif Natoms is not None:
            self.setPerm(perm)
        else:
            self.perm = perm
        self.malign = self.__call__
    ##
    def setPerm(self, perm):
        self.Natoms = sum(map(len,perm))
        self.perm = perm
        self.nperm = len(perm)
        self.npermsize = map(len, perm)
        self.permgroup = np.concatenate([np.asanyarray(p)+1 for p in perm])
        self.fast.fastoverlaputils.setperm(self.Natoms, self.permgroup, self.npermsize)
    ##
    def align(self, pos1, pos2, perm=None, invert=True, debug=False):
        return self(pos1, pos2, perm, invert, 1, debug)
    ##
    def __call__(self, pos1, pos2, perm=None, invert=True, nrot=10, debug=False):
        """ Aligns two isolated structures, aligns and permutes pos2 to match
        pos1 as closely as possible. Returns the distance between the aligned 
        structures, and the rotated and centred coordinates of the two structures
        
        Parameters
        ----------
        pos1 : (Natoms, 3) array_like
            Coordinate array.
        pos2 : (Natoms, 3) array_like
            Coordinate array to align with pos1.
        nrot : int, optional
            Number of different displacements to test
        perm : sequence of arrays, optional
            Each array in perm represents a different permutation group
        invert : bool
            If true tests inverted configurations as well
        debug : bool, optional
            If true prints debug information
        
        Returns
        -------
        distance : float
            Euclidean distance between X1 and X2
        X1 : (Natoms, 3) array
            pos1 centred on the origin
        X2 : (Natoms, 3) array
            Aligned coordinates of pos2
        rmatbest : (3,3) array
            the rotation matrix that maps centred pos2 onto X2
        """
        if perm is not None:
            self.setPerm(perm)
        elif len(pos1) != self.Natoms:
            self.Natoms = len(pos1)
            self.setPerm([np.arange(self.Natoms)])
            
        coordsb = np.asanyarray(pos1).flatten()
        coordsa = np.asanyarray(pos2).flatten()
        self.fast.clusterfastoverlap.setcluster()
        self.fast.commons.perminvopt = invert
        args = (coordsb,coordsa,debug,self.nmax,self.Jmax,self.harmscale,self.scale,nrot)
        dist, _, rmatbest = self.clus.alignharm(*args)
        return dist, coordsb.reshape(self.Natoms,3), coordsa.reshape(self.Natoms,3), rmatbest
    ##
    def compareList(self, poslist, perm=None):
        """ Calculates the maximum and average overlap of a list of coordinates
        
        Parameters
        ----------
        poslist : (nlist, Natoms, 3) array_like
            list of coordinates that are being compared
        perm : sequence of arrays, optional
            Each array in perm represents a different permutation group
            
        Returns
        -------
        avgoverlap : (nlist,nlist) array
            Array storing the array of average overlaps of the poslist
        maxoverlap : (nlist,nlist) array
            Array storing the array of maximum overlaps of the poslist
        navgoverlap : (nlist,nlist) array
            geometric mean normalised avgoverlap so,
            navgoverlap[i,j] = avgoverlap[i,j]/sqrt(avgoverlap[i,i]*avgoverlap[j,j])
        nmaxoverlap : (nlist,nlist) array
            geometric mean normalised avgoverlap so,
            nmaxoverlap[i,j] = maxoverlap[i,j]/sqrt(maxoverlap[i,i]*maxoverlap[j,j])
        """
        coords = np.array(poslist)
        nlist, Natoms, dim = coords.shape
        assert dim == 3
        if perm is None:
            if Natoms != self.Natoms:
                self.Natoms = Natoms
                self.setPerm([np.arange(self.Natoms)])
        else:
            self.setPerm(perm)
        # Center coordinates
        coords -= coords.mean(1)[:,None,:]
        coordslist = np.rollaxis(coords.reshape(nlist,-1),-1)
        avgoverlap, maxoverlap = self.clus.calcoverlapmatrices(
            coordslist,self.nmax,self.Jmax,self.harmscale,self.scale)
        diagavg = avgoverlap.diagonal()
        diagmax = maxoverlap.diagonal()
        navgoverlap = avgoverlap / sqrt(diagavg[:,None]*diagavg[None,:])
        nmaxoverlap = maxoverlap / sqrt(diagmax[:,None]*diagmax[None,:])
        return avgoverlap, maxoverlap, navgoverlap, nmaxoverlap
    
    
    
if __name__ == "__main__":
    import os
    import csv
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
    if f90.have_fortran:
        soapf = SphericalAlignFortran(scale, Jmax)
        harmf = SphericalHarmonicAlignFortran(scale, Jmax=Jmax, harmscale=harmscale, nmax=Nmax)
    
    print 'testing alignment on example LJ38 data, distance should = 1.4767'
    print 'SphericalAlign                Alignment:', soap(pos1, pos2)[0]
    print 'SphericalHarmonicAlign        Alignment:', harm(pos1, pos2)[0]
    if f90.have_fortran:
        print 'SphericalAlignFortran         Alignment:', soapf(pos1, pos2)[0]
        print 'SphericalHarmonicAlignFortran Alignment:', harmf(pos1, pos2)[0]
        
    print ''
    print 'Checking inversion isomers'
    
    print 'SphericalAlign                Alignment:', soap(pos1, -pos2)[0]
    print 'SphericalHarmonicAlign        Alignment:', harm(pos1, -pos2)[0]
    if f90.have_fortran:
        print 'SphericalAlignFortran         Alignment:', soapf(pos1, -pos2)[0]
        print 'SphericalHarmonicAlignFortran Alignment:', harmf(pos1, -pos2)[0]

    ###########################################################################
       
    print ''
    print 'testing alignment on randomly generated data, distance should ~ 0'
    # Testing on synthetic data
    rot = np.random.random((3,)) * np.array([2*pi, pi, 2*pi])
    rotM = EulerM(*rot)
    pos3 = np.random.normal(size=(300,3))*2
    pos3 -= pos3.mean(0)[None,:]
    pos4 = pos3.dot(rotM.T)
    
    print 'SphericalAlign                Alignment:', soap(pos3, pos4)[0]
    print 'SphericalHarmonicAlign        Alignment:', harm(pos3, pos4)[0]
    if f90.have_fortran:
        print 'SphericalAlignFortran         Alignment:', soapf(pos3, pos4)[0]
        print 'SphericalHarmonicAlignFortran Alignment:', harmf(pos3, pos4)[0]
    
    print ''
    print 'Checking inversion isomers'
    
    print 'SphericalAlign                Alignment:', soap(pos3, -pos4)[0]
    print 'SphericalHarmonicAlign        Alignment:', harm(pos3, -pos4)[0]
    if f90.have_fortran:
        print 'SphericalAlignFortran         Alignment:', soapf(pos3, -pos4)[0]
        print 'SphericalHarmonicAlignFortran Alignment:', harmf(pos3, -pos4)[0]
        
    print ''
    print 'Comparing rotation angles, Euler angles:'
    print rot
    print 'Rotation matrix:'
    print rotM.T
    soft = SOFT(Jmax+1)
    res = soap.calcSO3Coeffs(pos3, pos4)
    fout = soft.iSOFT(res).real
    R = soft.indtoEuler(findMax(fout))
    print 'Euler angles from obtained from SOFT'
    print R
    if f90.have_fortran:
        print 'Rotation matrix from minpermdist'
        print soapf(pos3, -pos4)[3]
        
#if  __name__ == '__main__':
#    scale = 0.3
#    Jmax = 7
#    Nmax = 15
#    rot = np.random.random((3,)) * np.array([2*pi, pi, 2*pi])
#    pos1 = np.random.normal(size=(60,3))*2
#    pos1 -= pos1.mean(0)[None,:]
#    pos2 = pos1.dot(EulerM(*rot).T)
#    #
#    soap = SphericalAlign(scale, Jmax)
#    harm = SphericalHarmonicAlign(scale, 1.0, Nmax, Jmax)
#    if have_fortran:
#        soapf = SphericalAlignFortran(scale, Jmax)
#        harmf = SphericalHarmonicAlignFortran(scale, Jmax=Jmax, harmscale=1.0, nmax=Nmax)
#    soft = SOFT(Jmax+1)
#    # res = soap.calcCoeff(pos1, pos2)
#    # fout = soft.iSOFT(res.conj()).real
#    res = soap.calcSO3Coeffs(pos1, pos2)
#    fout = soft.iSOFT(res).real
#    R = soft.indtoEuler(findMax(fout))
#    print rot, R
#    #
#    print soap(pos1, pos2)[0], harm(pos1, pos2)[0]
#    print soap(pos1, -pos2)[0], harm(pos1, -pos2)[0]
