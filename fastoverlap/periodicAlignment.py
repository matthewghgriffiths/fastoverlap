'''
Created on 19 Jan 2017

@author: mg542
'''

from itertools import product, combinations_with_replacement, chain, izip

import numpy as np
from numpy.linalg import norm
from scipy import median

from utils import find_best_permutation, findMax, _next_fast_len

try:
    import fastoverlap_f as fast
    have_fortran=True
except ImportError:
    have_fortran=False

class BasePeriodicAlignment(object):
    def findDisps(self, pos1, pos2):
        raise NotImplementedError
    ##
    def align(self, pos1, pos2):
        disps = self.findDisps(pos1, pos2)
        return self.refine(pos1, pos2, disps)
    ##
    def refine(self, x, y, disps):
        """
        Given a displacement vector between two non-permutationally aligned
        configurations, this will refine the displacement vector by 
        performing a permutational alignment followed by a calculation of the median
        displacement, followed by a second permutational alignment, and
        then finally setting the mean displacement to 0.
        """
        disps = np.atleast_2d(disps)
        distperm = [self.Hungarian(x, y - disp[None,:]) + (disp,)
                    for disp in disps]
        dist, perm, disp = min(distperm, key=lambda x: x[0])
        dxs = self.get_disp(x, (y-disp)[perm])
        mediandx = median(dxs, axis=0)
        disp += mediandx
        perm = self.Hungarian(x, y - disp[None,:])[1]
        pos1 = self.periodic(x, True)
        pos2 = (y - disp)[perm]
        ## subtract mean displacement vector
        dx = self.get_disp(pos1, pos2).mean(0)
        disp -= dx
        pos2 += dx[None,:]
        self.periodic(pos2)
        dist = self.get_dist(pos1, pos2)
        return dist, pos1, pos2, perm, disp
    ##
    def periodic(self, x, copy=False):
        if copy:
            x = x.copy()
        x -= np.round(x / self.boxvec) * self.boxvec
        return x
    ##
    def get_disp(self, X1, X2):
        return self.periodic(X1 - X2)
    ##
    def get_dist(self, X1, X2):
        return norm(self.get_disp(X1, X2))
    ##
    def cost_matrix(self, X1, X2):
        """
        Calculating this matrix is most of the computational cost associated
        with this algorithm for large databases, might be worth reimplementing 
        this in c++/fortran?
        """
        disps = X1[None,:,:] - X2[:,None,:]
        disps -= np.round(disps / self.boxvec[None,None,:]) * self.boxvec[None,None,:]
        return norm(disps, axis=2)
    ##
    def Hungarian(self, X1, X2):
        _, permlist = find_best_permutation(X1, X2, self.perm, user_cost_matrix=self.cost_matrix)
        dist = self.get_dist(X1, X2[permlist])
        return dist, permlist
    ##
    def __call__(self, pos1, pos2):
        return self.align(pos1, pos2)

class FourierAlign(BasePeriodicAlignment):
    """
    Alignment procedure based on a maximum likelihood method. 
    It's probably better to use the PeriodicGaussian class to align 
    periodic systems, this is included for completeness.
    
    Parameters
    ----------
    Natoms : int 
    boxvec : array like floats
        defines periodicity of system being aligned
    permlist : optional
        list of allowed permutations. If nothing is given, all atoms will be
        considered as permutable. For no permutations give an empty list []
    
    
    For two structures $\vec{R}^0_j$ and $\vec{R}^1_j$
    
    We want to find the optimal alignment such that we have a permutation matrix
    $\vec{P}$, global displacement $\vec{d}$ such that the local displacement
    $\sum |\vec{e}_j|^2$ is minimized.
    
    \begin{equation}
    \vec{R}^0_j = \vec{P}(\vec{R}^1_j + \vec{e}_j + \vec{d}
    \end{equation}

    We can define fourier coefficients for both structures as follows,     
    
    \begin{align}
    \tilde{C}_\vec{k}^0 = \frac{1}{L^3} \sum_j e^{-i\vec{k}\cdot\vec{R}_j}
    \tilde{C}_\vec{k}^1 = \frac{1}{L^3} \sum_j e^{-i\vec{k}\cdot\vec{R}'_j}
    \end{align}
    
    The ratio of these coefficients will be, 
    
    \begin{equation}
    \frac{\tilde{C}_{\vec{k}}^1}{\tilde{C}_{\vec{k}}^0}=
    e^{i\vec{k}\cdot\vec{d}} \frac
    {\sum_{j=1}^N e^{-i\vec{k}\cdot(\vec{R}_j+\vec{e}_j)}}
    {\sum_{j=1}^N e^{-i\vec{k}\cdot\vec{R}_j}}.
    \end{equation}
    
    Taking a first order approximation of the local displacements we find
    
    \begin{equation}
    \vec{k}\cdot\vec{d} = -i \textrm{log}
    \left(\frac{\tilde{C}_{\vec{k}}^1}{\tilde{C}_{\vec{k}}^0} \right) 
    + i \log{\left( 1 - \frac{i\vec{k}\cdot\tilde{\vec{e}}_{\vec{k}}}{\tilde{C}_{\vec{k}}^0} \right)}
    +  2\pi n_\vec{k}
    \end{equation}
    
    To find $\vec{d}$ we define $\vec{d}=\vec{d}_0+\vec{d}_1$ so
    
    \begin{align}
    \theta_{\vec{k}}^{\vec{d}_0} &=
    \Re{\left[
    - i
    \textrm{log}
    \left(
    \frac{\tilde{C}_{\vec{k}}^1}{\tilde{C}_{\vec{k}}^0} e^{i\vec{k}\cdot\vec{d}_0} 
    \right)\right]} \\
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    \frac{|\vec{k}|}{|\tilde{C}_{\vec{k}}^0|} \hat{e}_{\vec{k}}
    &=
    \Re{\left[
    \frac{\vec{k}\cdot\tilde{\vec{e}}_{\vec{k}}}
    {\tilde{C}_{\vec{k}}^0}
    \right]},
    \end{align}

    This means we can assume that $2 \pi n_\vec{k}$ is 0. 
    If we assume that the standard deviation of $\hat{e}_{\vec{k}}$ 
    can be estimated to be, 
    $\sigma(\hat{e}_{\vec{k}}) \approx F(1+|\vec{p}|^{-3})$, 
    and we define $\sigma_\vec{k} = \sigma(\hat{e}_{\vec{k}}) 
    |\vec{k}|\left/|\tilde{C}_{\vec{k}}^0|\right.$ 
    then the log-likelihood of $\vec{d}$ can be calculated.
	
    \begin{align}
    \log{\Pr(\vec{d}_1|\tilde{C}_{\vec{k}}^1,\tilde{C}_{\vec{k}}^0,
    \vec{d}_0,\sigma_\vec{k})} &= -
    \sum_{\vec{k}\,\epsilon\,\vec{K}}  
    \frac{(\vec{k}\cdot\vec{d_1}-\theta_\vec{k}^{\vec{d}_0})^2}
    {2\sigma_{\vec{k}}^2}
    + \frac{1}{2}\log{2\pi} + \log{\sigma_{\vec{k}}}.
    \end{align}
    
    At the maximum likelihood the gradient of \ref{eqn:dloglikelihoodFull} 
    will be zero which  we can solve to find,
    
    \begin{align}
    \vec{d}_1 &= 
    \left(
    \sum_{\vec{k}\,\epsilon\,\vec{K}}  \frac{\vec{k}\otimes\vec{k}}{\sigma_{\vec{k}}^2} 
    \right)^{-1}
    \sum_{\vec{k}\,\epsilon\,\vec{K}}  \vec{k}
    \frac{\theta_\vec{k}^{\vec{d}_0}}
    {\sigma_{\vec{k}}^2},
    \end{align}
    """
    def __init__(self, Natoms, box, permlist=None, dim=3, origin=None, 
                 maxd=None, tol=None, maxdisps=10, likeThresh=2.):
        self.Natoms = Natoms
        self.boxvec = np.asanyarray(box, dtype=float)
        self.dim = dim
        self.likeThresh = likeThresh
        if permlist is None:
            self.perm = [np.arange(self.Natoms)]
        else:
            self.perm = permlist
        self.permlist = [[p] for p in self.perm] + [self.perm] + [[np.arange(self.Natoms)]]
        if tol is None:
            self.tol = (np.prod(self.boxvec)/self.Natoms)**(1./self.dim) * 1e-2
        else:
            self.tol = tol
        self.maxdisps = maxdisps
        ####
        if maxd is None:
            maxd = np.prod(self.boxvec)**(1./self.dim)/5.
        self.setMaxDisplacment(maxd)
        ###
        if origin is None:
            self.origin=np.zeros_like(self.boxvec)
        else:
            self.origin = np.asanyarray(origin)
        self.minVals = self.origin - 0.5 * self.boxvec
        self.maxVals = self.origin + 0.5 * self.boxvec
    @classmethod
    def calcWaveVectors(cls, n, p):
        def setindex(inds):
            out = [0]*n
            for i in inds:
                out[i] += 1
            return out
        def negIndices(out):
            nonzero = [i for i, a in enumerate(out) if a!=0]
            ret = [out]
            for i in nonzero[1:]:
                extend = []
                for inds in ret:
                    extend.append([-a if i==j else a for j, a in enumerate(inds)])
                ret.extend(extend)
            return ret
        indices =  list(chain(*(negIndices(setindex(combin)) for r in range(1,p+1)
                                for combin in combinations_with_replacement(range(n), r))))
        return indices
    @classmethod
    def removeDuplicates(cls, disps, tol, *args):
        n = len(disps)
        triu = np.triu_indices(n,1)
        dists = norm(disps[triu[0]]-disps[triu[1]], axis=1)
        remove = (dists<tol).nonzero()[0]
        tokeep = np.setdiff1d(np.arange(n), triu[1][remove])
        if args:
            return [disps[tokeep]] + [a[tokeep] for a in args]
        else:
            return disps[tokeep]
    ##
    def setMaxDisplacment(self, maxd):
        # This routines calculates, for a given threshold distance move the set of
        # starting displacements to calculate the MLE of the displacements.
        self.n = int(np.ceil(self.boxvec.max()/maxd))
        ps = np.array(self.calcWaveVectors(self.dim, self.n))
        ks = 2.*np.pi/self.boxvec[None,:]*ps
        absks = norm(ks, axis=1)
        argsort = absks.argsort()
        maxI = np.searchsorted(absks[argsort], 2.*np.pi/maxd)
        self.ps = ps[argsort[:maxI]]
        self.permps = np.r_[(self.ps,)*len(self.perm)]
        self.absps = norm(self.ps, axis=1)
        self.ks = ks[argsort[:maxI]]
        self.permks = np.r_[(self.ks,)*len(self.perm)]
        self.splits = np.diff(self.absps).nonzero()[0]+1
        self.psplit = np.split(self.ps, self.splits)
        maxp = int(np.ceil(self.absps.max()))
        self.startDisps = np.array(list(product(*[xrange(-maxp, maxp+1)]*self.dim)),dtype=float)
        self.startDisps *= 0.5 * self.boxvec[None,:] / self.absps.max()
    ##
    def FourierShift(self, pos1, pos2):
        disps, logLikes = [np.concatenate(a) for a in 
                           zip(*[self.removeDuplicates(disp[:self.maxdisps], self.tol, likes[:self.maxdisps])
                                 for disp, likes in ((d[like-like[0]<self.likeThresh],
                                                      like[like-like[0]<self.likeThresh])
                                 for d, like in (self.calcDispLogLike(pos1, pos2, self.startDisps, p) 
                                                 for p in self.permlist))])]
        disps, logLikes = self.removeDuplicates(disps, self.tol, logLikes)
        return disps, logLikes
    ##
    def calcDispLogLike(self, pos1, pos2, startDisps, perm=None):
        startDisps = np.atleast_2d(startDisps)
        if perm is None:
            perm = [np.arange(len(pos1))]
        nperm = len(perm)
        C1s = [np.exp(-1j * (pos1[p,None,:] * self.ks[None,...]).sum(2))
               for p in perm]
        C2s = [np.exp(-1j * (pos2[p,None,:] * self.ks[None,...]).sum(2))
               for p in perm]
        C1 = np.concatenate([C.sum(0) for C in C1s])
        C2 = np.concatenate([C.sum(0) for C in C2s])
        absC1, absC2 = np.abs(C1), np.abs(C2)
        angC1, angC2 = np.angle(C1), np.angle(C2)
        theta = (angC1-angC2)
        ## Estimating distribution of displacement vector
        nk = len(self.ks)
        Fs = [2*(np.abs(absC1[i*nk:(i+1)*nk] - absC2[i*nk:(i+1)*nk])).mean()+1e-6
             for i in xrange(nperm)]
        std = np.concatenate([F*(norm(self.ps,axis=1)**-3.+1.)
                              for F in Fs])
        absC1C2 = 2*(absC1**-1+absC2**-1)**-1
        permks = np.r_[(self.ks,)*nperm]
        stdk = std * norm(permks,axis=1) / absC1C2
        fs = stdk**-2
        ##
        A = (permks[:,None,:]*permks[:,:,None]*fs[:,None,None]).sum(0)
        invA = np.linalg.inv(A)
        dtheta = theta[:,None]-permks.dot(startDisps.T)
        dtheta -= np.round(dtheta / (2*np.pi)) * 2 * np.pi
        disps = np.einsum("jk,li,l,lk->ij", invA, dtheta, fs, permks)
        logLikes = (fs[:,None]*(permks.dot(disps.T)-dtheta)**2).sum(0) - np.log(stdk[:,None]).sum(0)
        argsort = logLikes.argsort()
        disps[:] = (startDisps+disps)[argsort]
        logLikes[:] = logLikes[argsort]/nperm
        disps -= np.round((disps-self.origin)/self.boxvec) * self.boxvec[None,:]
        return disps, logLikes
    ##
    def findDisps(self, pos1, pos2):
        return self.FourierShift(pos1, pos2)[0]
    
    
class PeriodicAlign(BasePeriodicAlignment):
    """
    Finds the best alignment between two configurations of a periodic system
    
    Parameters
    ----------
    Natoms : int 
    boxvec : array like floats
        defines periodicity of system being aligned
    permlist : optional
        list of allowed permutations. If nothing is given, all atoms will be
        considered as permutable. For no permutations give an empty list []
    scale : optional float
        determins the size of the Gaussian kernels automatically set to be
        1/3 of interatomic separation
    maxk : optional float
        the value of wavevector at which the calculation is truncated
        the higher this is set the more accurate the overlap calculation
    dim : optional int
        dimensionality of the system TODO: TEST
    """
    def __init__(self, Natoms, boxvec, permlist=None, dim=3, 
                 scale=None, maxk=None):       
        self.Natoms = Natoms
        self.boxvec = np.array(boxvec, dtype=float)
        self.dim = dim
        if permlist is None:
            self.perm = [np.arange(self.Natoms)]
        else:
            self.perm = map(np.array, permlist)
        self.pos1 = np.zeros((self.Natoms, self.dim))
        self.pos2 = np.zeros((self.Natoms, self.dim))
        if scale is None:
            scale = (np.prod(self.boxvec)/self.Natoms)**(1./self.dim)/3.
        if maxk is None:
            maxk = 1.5 / scale
        self.scale = scale
        self.factor = 2 * (np.pi*self.scale**2)**(-self.dim*0.5)*self.scale**2 / np.prod(self.boxvec)
        self.setks(maxk)
        self.setPos(self.pos1, self.pos2)
        #self.setPotentials()
    ##
    def setks(self, maxk):
        self.n = int(np.ceil(2*np.pi/self.boxvec.min()*maxk))
        ps = np.indices((self.n*2+1,)*self.dim) - self.n
        self.ks = 2.*np.pi/self.boxvec[(slice(None),)+(None,)*self.dim]*ps
        self.absks = norm(self.ks, axis=0)
        self.C1 = np.empty((len(self.perm),)+self.absks.shape, dtype=complex)
        self.C2 = np.empty((len(self.perm),)+self.absks.shape, dtype=complex)
        self.C = np.empty_like(self.absks, dtype=np.complex128)
        shape = np.array(self.C.shape)*2+1
        self.fshape = tuple(_next_fast_len(int(d)) for d in shape)
        self.f = np.empty(self.fshape, dtype=complex)
        self.fabs = np.empty(self.fshape, dtype=float)
    ##
    def setScale(self, scale):
        self.scale = scale
        self.factor = 2 * (np.pi*self.scale**2)**(-self.dim*0.5)*self.scale**2 / np.prod(self.boxvec)
    ##
    def calcFourierCoeff(self, pos, out=None):
        if out is None:
            out = np.empty((len(self.perm),)+self.absks.shape, dtype=complex)
        for p, C in izip(self.perm, out):
            np.exp(-1j * np.tensordot(pos[p,:], self.ks, [1,0])).sum(0, out=C)
            C[(self.n,)*self.dim] = 0
        return out
    ##
    def setPos(self, pos1=None, pos2=None, Cs=None):
        if pos1 is not None:
            self.pos1 = pos1
        if pos2 is not None:
            self.pos2 = pos2
        if Cs is None:
            self.calcFourierCoeff(pos1, self.C1)
            self.calcFourierCoeff(pos2, self.C2)
        else:
            self.C1, self.C2 = Cs
        self.Csum = (((self.C1*self.C1.conj()).real + (self.C2*self.C2.conj()).real) * 
                     np.exp(-self.absks[None,:]**2 * (self.scale**2))).real.sum() * 0.5 * self.factor
        (self.C1 * self.C2.conj() * np.exp(-self.absks**2 * (self.scale**2))).sum(0, out=self.C)
        self.f[:] = np.fft.fftn(self.C, self.fshape)
        np.abs(self.f, out=self.fabs)
    ##
    def findDisps(self, pos1, pos2, Cs=None):
        self.setPos(pos1, pos2, Cs)
        disp = findMax(self.fabs) * self.boxvec / self.fabs.shape
        return disp[None,:]
    ##
    def align(self, pos1, pos2, Cs=None):
        disps = self.findDisps(pos1, pos2, Cs)
        return self.refine(pos1, pos2, disps)
    ##
    def alignGroup(self, coords, keepCoords=False):
        n = len(coords)
        if keepCoords:
            aligned = np.empty((2, n, n, self.Natoms, self.dim))
        coeffs = [self.calcFourierCoeff(p) for p in coords]
        dists = np.zeros((n, n))
        for i, (pos1, c1) in enumerate(izip(coords, coeffs)):
            for j, (pos2, c2) in enumerate(izip(coords, coeffs)):
                dist, x1, x2 = self.align(pos1, pos2, [c1,c2])[:3]
                if keepCoords:
                    aligned[0,i,j,...] = x1
                    aligned[1,i,j,...] = x2
                dists[i,j] = dist
        if keepCoords:
            return dists, aligned
        else:
            return dists
            
class PeriodicAlignFortran(BasePeriodicAlignment):
    ##
    def __init__(self, Natoms, boxVec, permlist=None, dim=3, 
                 scale=None, maxk=None):       
        self.Natoms = Natoms
        self.boxVec = np.array(boxVec, dtype=float)
        self.dim = dim
        if permlist is None:
            self.setPerm([np.arange(self.Natoms)])
        else:
            self.setPerm(permlist)
    ##
    def setPerm(self, perm):
        self.perm = perm
        self.nperm = len(permlist)
        self.npermsize = map(len, permlist)
        self.permgroup = np.concatenate([np.asanyarray(p)+1 for p in permlist])
        fast.bulkfastoverlap.initialise(self.Natoms, self.permgroup, self.npermsize)
    ##
    def align(self, pos1, pos2, ndisps=10, perm=None):
        if perm is not None:
            self.setPerm(perm)
        coordsb = pos1.flatten()
        coordsa = pos2.flatten()
        args = (coordsb, coordsa,True,
                self.boxVec[0],self.boxVec[1],self.boxVec[2],False,False,ndisps)
        dist, dist2, disp, disps = fast.bulkfastoverlap.align(*args)
        return dist, coordsb.reshape(pos1.shape), coordsa.reshape(pos2.shape), disp
    
if __name__ == "__main__":
    print 'testing alignment on example data, distance should = 1.559'
    import os
    import csv
    datafolder = "../examples/BLJ256"
    def readFile(filename):
        with open(filename, 'rb') as f:
            reader = csv.reader(f, delimiter=' ')
            dist = [map(float, row) for row in reader]
        return np.array(dist)
    
    pos1 = readFile(os.path.join(datafolder, 'coords'))
    pos2 = readFile(os.path.join(datafolder, 'finish'))
    
    natoms = 256
    ntypeA = 204
    shape = (natoms, 3)
    boxSize = np.ones(3)*5.975206329 
    permlist = [np.arange(ntypeA), np.arange(ntypeA, natoms)]
    overlap = PeriodicAlign(256, boxSize, permlist)
    align = FourierAlign(256, boxSize, permlist)
    if have_fortran:
        overlap_f = PeriodicAlignFortran(256, boxSize, permlist)
        
    print 'PeriodicAlign aligment:', overlap(pos1, pos2)[0]
    if have_fortran:
        print 'PeriodicAlignFortran aligment:', overlap_f(pos1, pos2)[0]
    

    """
    import timeit
    print 'Timing python implementation'
    print  timeit.timeit("overlap(pos1, pos2)", setup="from __main__ import overlap, pos1, pos2")    
    if  have_fortran:
        print 'Timing Fortran implementation'
        print  timeit.timeit("overlap_f.align(pos1, pos2,1)", setup="from __main__ import overlap_f, pos1, pos2")    
    """
    """
    Testing Run time in Ipython
    %prun dists, aligned = overlap.alignGroup(coords13, True)
    """
        
        
        