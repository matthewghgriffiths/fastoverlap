# -*- coding: utf-8 -*-
"""
Created on Wed Feb 15 19:06:01 2017

@author: mg542
"""

from heapq import heappop, heappush

import numpy as np
from numpy.linalg import norm
from numpy import sin, cos, sqrt, pi, exp
from scipy.special import eval_jacobi, gamma
from scipy.special.orthogonal import genlaguerre, eval_genlaguerre

try:
    from scipy.special.basic import factorial, comb
except ImportError:
    from scipy.misc import factorial, comb


from scipy.optimize import curve_fit
from scipy.spatial.distance import cdist

"""
Tries to import code from pele:
Documentation: http://pele-python.github.io/pele/

Source Code: https://github.com/pele-python/pele

If it can't import the code from pele, manually redefine needed functions
these functions are direct copies from pele
"""
try:
    import munkres
    hungarian = munkres.Munkres().compute

    def lap(cost):
        """Solve Linear Assignment problem
        """
        return [pair[1] for pair in sorted(hungarian(cost))]
except:
    def lap(cost):
        """Solve Linear Assignment problem
        """
        raise NotImplementedError
        

try:
    from pele.mindist.rmsfit import findrotation
    from pele.mindist import find_best_permutation
    
except ImportError:
    
    import munkres
    hungarian = munkres.Munkres().compute
    
    def _make_cost_matrix(X1, X2):
        """
        return the cost matrix for use in the hungarian algorithm.
        
        the cost matrix is the distance matrix (squared) for all atoms 
        in atomlist
        """
        return cdist(X1, X2, 'sqeuclidean')
    
    def find_permutations_hungarian( X1, X2, 
                                    make_cost_matrix=_make_cost_matrix ):
        cost = make_cost_matrix(X1, X2)
        #########################################
        # run the hungarian algorithm
        #########################################
        perm = lap(cost)
        #########################################
        # apply the permutation
        #########################################
        # TODO: how to get new distance?
        dist = -1
        return dist, perm
        
    def find_best_permutation(X1, X2, permlist=None, user_algorithm=None, 
                                 reshape=True, 
                                 user_cost_matrix=_make_cost_matrix,
                                 **kwargs):
        """
        find the permutation of the atoms which minimizes the distance |X1-X2|
        
        With all the default parameters, findBestPermutation assumes that 
        X1, X2 are arrays of atoms in 3d space and performs reshaping on the 
        coordinates. However, if you want to pass a 2D system or a custom 
        array with own cost function, you can turn automatic reshaping off. 
        
        Parameters
        ----------
        X1, X2 : 
            the structures to align
        permlist : a list of lists
            A list of lists of atoms which are interchangable.
            e.g. for a 50/50 binary mixture::
            
                permlist = [range(1,natoms/2), range(natoms/2,natoms)]
            
            If permlist is None all atoms are assumed to be permutable.
    
        user_algoriithm : None or callable
            you can optionally pass which algorithm to use.
        gen_cost_matrix : None or callable
            user function to generate the cost matrix
        reshape : boolean
            shall coordinate reshaping be performed.
        box_lengths : float array
            array of floats giving the box lengths for periodic boundary conditions.
            Set to None for no periodic boundary conditions.
        
        Returns
        -------
        dist : float
            the minimum distance WARNING: THIS IS NOT NECESSARILY CORRECT, IT SHOULD BE 
            RECALCULATED.  THIS WILL BE REMOVED IN THE FUTURE.
        perm:
            a permutation which will best align coords2 with coords1
        
        Notes
        -----
        For each list of interchangeable atoms in permlist the permutation
        which minimizes the distance between the two structures is found.  This minimimization
        is done by mapping the problem onto the linear assignment problem which can then be solved
        using graph theoretic techniques.  
        
        http://en.wikipedia.org/wiki/Linear_assignment_problem
        http://en.wikipedia.org/wiki/Hungarian_algorithm
    
        there are several packages in pypi which solve the linear assignment problem
        
        hungarian : c++ code wrapped in python.  scales roughly like natoms**2.5
        
        munkres : completely in python. scales roughly like natoms**3.  very slow for natoms > 10
        
        in addition we have wrapped the OPTIM version for use in pele.  It uses the sparse 
        version of the Jonker-Volgenant algorithm.  Furthermore the cost matrix calculated in 
        a compiled language for an additional speed boost. It scales roughly like natoms**2
    
        """
        if reshape:
            X1 = X1.reshape([-1,3])
            X2 = X2.reshape([-1,3])
        
        if permlist is None:
            permlist = [list(range(len(X1)))]
        
        newperm = list(range(len(X1)))
        disttot = 0.
        
        for atomlist in permlist:
            if len(atomlist) == 0:
                continue
            if user_algorithm is not None:
                dist, perm = user_algorithm(X1[atomlist], X2[atomlist], make_cost_matrix=user_cost_matrix, **kwargs)
            else:
                dist, perm = find_permutations_hungarian(X1[atomlist], X2[atomlist], make_cost_matrix=user_cost_matrix, **kwargs)

            disttot += dist**2
            for atom, i in zip(atomlist,list(range(len(atomlist)))):
                newperm[atom] = atomlist[perm[i]]
        dist = sqrt(disttot)
        return dist, newperm

    def findrotation(x1, x2, align_com=True):
        """
        Return the rotation matrix which aligns XB with XA
        #
        Return the matrix which
        aligns structure XB to be as similar as possible to structure XA.
        To be precise, rotate XB, so as to minimize the distance |XA - XB|.
        #
        Rotations will be done around the origin, not the center of mass
        #
        Rotational alignment follows the prescription of
        Kearsley, Acta Cryst. A, 45, 208-210, 1989
        http://dx.doi.org/10.1107/S0108767388010128
        """
        if x1.size != x2.size:
            raise ValueError("dimension of arrays does not match")
        # reshape the arrays
        x1 = x1.reshape([-1,3]).copy()
        x2 = x2.reshape([-1,3]).copy()
        # determine number of atoms
        natoms = x1.shape[0]
        # set both com to zero
        if align_com:
            com1 = np.sum(x1,axis=0) / float(natoms)
            com2 = np.sum(x2,axis=0) / float(natoms)
            x1 -= com1
            x2 -= com2
        x1 = x1.ravel() 
        x2 = x2.ravel()
        # TODO: this is very dirty!
        #########################################
        # Create matrix QMAT
        #########################################
        QMAT = np.zeros([4,4], np.float64)
        for J1 in range(natoms):
            J2 = 3* J1 -1
            XM = x1[J2+1] - x2[J2+1]
            YM = x1[J2+2] - x2[J2+2]
            ZM = x1[J2+3] - x2[J2+3]
            XP = x1[J2+1] + x2[J2+1]
            YP = x1[J2+2] + x2[J2+2]
            ZP = x1[J2+3] + x2[J2+3]
            QMAT[0,0] = QMAT[0,0] + XM**2 + YM**2 + ZM**2
            QMAT[0,1] = QMAT[0,1] - YP*ZM + YM*ZP
            QMAT[0,2] = QMAT[0,2] - XM*ZP + XP*ZM
            QMAT[0,3] = QMAT[0,3] - XP*YM + XM*YP
            QMAT[1,1] = QMAT[1,1] + YP**2 + ZP**2 + XM**2
            QMAT[1,2] = QMAT[1,2] + XM*YM - XP*YP
            QMAT[1,3] = QMAT[1,3] + XM*ZM - XP*ZP
            QMAT[2,2] = QMAT[2,2] + XP**2 + ZP**2 + YM**2
            QMAT[2,3] = QMAT[2,3] + YM*ZM - YP*ZP
            QMAT[3,3] = QMAT[3,3] + XP**2 + YP**2 + ZM**2
        #
        QMAT[1,0] = QMAT[0,1]
        QMAT[2,0] = QMAT[0,2]
        QMAT[2,1] = QMAT[1,2]
        QMAT[3,0] = QMAT[0,3]
        QMAT[3,1] = QMAT[1,3]
        QMAT[3,2] = QMAT[2,3]
        ###########################################
        """
        Find eigenvalues and eigenvectors of QMAT.  The eigenvector corresponding
        to the smallest eigenvalue is the quaternion which rotates XB into best
        alignment with XA.  The smallest eigenvalue is the squared distance between
        the resulting structures.
        """
        ###########################################
        (eigs, vecs) = np.linalg.eig(QMAT)
    
        imin = np.argmin(eigs)
        eigmin = eigs[imin] # the minimum eigenvector
        Q2 = vecs[:,imin]  # the eigenvector corresponding to the minimum eigenvalue
        if eigmin < 0.:
            if abs(eigmin) < 1e-6:
                eigmin = 0.
            else:
                print(('minDist> WARNING minimum eigenvalue is ',eigmin,' change to absolute value'))
                eigmin = -eigmin
        #
        dist = sqrt(eigmin) # this is the minimized distance between the two structures
        #
        Q2 = np.real_if_close(Q2, 1e-10)
        if np.iscomplexobj(Q2):
            raise ValueError("Q2 is complex")
        return dist, q2mx(Q2)
        
    def q2mx(qin):
        """quaternion to rotation matrix"""
        Q = qin / np.linalg.norm(qin)
        RMX = np.zeros([3, 3], np.float64)
        Q2Q3 = Q[1] * Q[2]
        Q1Q4 = Q[0] * Q[3]
        Q2Q4 = Q[1] * Q[3]
        Q1Q3 = Q[0] * Q[2]
        Q3Q4 = Q[2] * Q[3]
        Q1Q2 = Q[0] * Q[1]
        #
        RMX[0, 0] = 2. * (0.5 - Q[2] * Q[2] - Q[3] * Q[3])
        RMX[1, 1] = 2. * (0.5 - Q[1] * Q[1] - Q[3] * Q[3])
        RMX[2, 2] = 2. * (0.5 - Q[1] * Q[1] - Q[2] * Q[2])
        RMX[0, 1] = 2. * (Q2Q3 - Q1Q4)
        RMX[1, 0] = 2. * (Q2Q3 + Q1Q4)
        RMX[0, 2] = 2. * (Q2Q4 + Q1Q3)
        RMX[2, 0] = 2. * (Q2Q4 - Q1Q3)
        RMX[1, 2] = 2. * (Q3Q4 - Q1Q2)
        RMX[2, 1] = 2. * (Q3Q4 + Q1Q2)
        return RMX
    

def _next_fast_len(target):
    if target <= 6:
        return target
    # Quickly check if it's already a power of 2
    if not (target & (target-1)):
        return target
    match = float('inf')  # Anything found will be smaller
    p5 = 1
    while p5 < target:
        p35 = p5
        while p35 < target:
            # Ceiling integer division, avoiding conversion to float
            # (quotient = ceil(target / p35))
            quotient = -(-target // p35)
            # Quickly find next power of 2 >= quotient
            try:
                p2 = 2**((quotient - 1).bit_length())
            except AttributeError:
                # Fallback for Python <2.7
                p2 = 2**(len(bin(quotient - 1)) - 2)
            N = p2 * p35
            if N == target:
                return N
            elif N < match:
                match = N
            p35 *= 3
            if p35 == target:
                return p35
        if p35 < match:
            match = p35
        p5 *= 5
        if p5 == target:
            return p5
    if p5 < match:
        match = p5
    return match

def multiargmax(a):
    a = np.asanyarray(a)
    return np.unravel_index(a.argmax(), a.shape)
    
def findMax(a):
    """
    This finds the interpolated max of multidimensional array
    result given as fractional coordinates of array size
    """
    a = np.asanyarray(a)
    shape = a.shape
    dim = len(shape)
    ind = np.unravel_index(a.argmax(), a.shape)
    i1 = tuple(tuple((i0+1)%shape[i] if i==j else i0 for j in range(dim))
               for i, i0 in enumerate(ind))
    i2 = tuple(tuple(i0 for j in range(dim))
               for i, i0 in enumerate(ind))
    i3 = tuple(tuple(i0-1 if i==j else i0 for j in range(dim))
               for i, i0 in enumerate(ind))
    y1 = np.abs(a[i1])
    y2 = np.abs(a[i2])
    y3 = np.abs(a[i3])
    d = (y3 - y1) / (2 * (2 * y2 - y1 - y3))
    return (np.array(ind) - d)
    
def indtoEuler(ind, n):
    ind = np.atleast_2d(ind)
    indFactor = np.array([2*pi/n, pi/n, 2*pi/n])
    rot = indFactor*ind
    rot[:,1] += 0.5*pi/n
    return rot.squeeze()

def _gaussian(x, A, mu, *alphax0):
    x = np.atleast_2d(x)
    dim = len(x)
    sigma = np.zeros((dim,dim), float)
    sigma[np.triu_indices(dim)] = alphax0[:-dim]
    x0 = x - np.array(alphax0[-dim:])[:,None]
    return A * np.exp(-np.einsum("ik,jk,ij->k",x0,x0,sigma)) + mu

def fitPeak(f, ind, n=2):
    dim = len(f.shape)
    peak = f[np.ix_(*[np.arange(i-n,i+n+1)%s for i,s in zip(ind,f.shape)])]
    flatindices = np.indices((2*n+1,)*dim).reshape((dim,-1))
    flatpeak = peak.flatten()
    p0 = ([peak[(n,)*dim], 0.] + 
          [1. if i==j else 0. for i in range(dim) for j in range(i, dim)] + 
          [0.]*dim)
    popt, pcov = curve_fit(_gaussian, flatindices-n, flatpeak, p0=p0)
    return popt, pcov                   
    
def findPeaks(a, npeaks=10, width=2):
    f = np.array(a)
    f -= f.min()
    dim = len(f.shape)
    indices = np.indices(f.shape).reshape((dim,-1))
    tupinds = tuple(indices)
    peaks = []
    amplitude = []
    mean = []
    sigma = []
    for i in range(npeaks):
        ind = multiargmax(f)
        try:
            popt = fitPeak(f, ind, width)[0]
            peaks.append((popt[-dim:] + ind))
            amplitude.append(popt[0])
            mean.append(popt[1])
            sigma.append((2*popt[2:-dim])**-0.5)
            popt[-dim:] += ind 
            f[tupinds] -= _gaussian(indices, *popt)
        except RuntimeError:
            break
        except ValueError:
            break
    peaks = np.array(peaks)
    if len(peaks) == 0:
        peaks = findMax(a)[None,:]
        amplitude.append(a.max())
        mean.append(0)
        sigma.append(np.NaN)
    return peaks, amplitude, mean, sigma, f
    
def eval_grad_jacobi(n, alpha, beta, x, out=None):
    fact = gamma(alpha + beta + n + 2)/2/gamma(alpha+beta+n+1)
    return eval_jacobi(n-1, alpha+1, beta+1, x, out) * fact
    
def BruteOverlap(pos1, pos2, scale):
    pos2, pos1 = list(map(np.atleast_2d, sorted([pos2, pos1], key=len)))
    inds = np.indices((len(pos1),len(pos2))).reshape((2,-1))
    rs = norm(pos1[inds[0]]-pos2[inds[1]], axis=1)
    return exp(-rs**2 / 4 / scale**2).sum() * (pi * scale**2)**(1.5)
    
def norm_harmonicBasis(n, l, r0):
    """
    returns the normalisaton constant of the harmonic basis function
    """
    return (2 * factorial(n) * r0**(-2*l-3) / gamma(1.5+n+l))**0.5
    
def coeffs_harmonicBasis(n, l, r0):
    """
    Calculates the coefficients of the polynomial that defines the
    radial harmonic basis function of degree n and angular momentum order l
    r0 is the scale of the harmonic basis function
    
    returns \tilde{g}^{nl}_s given the following    
    r^l L^{l+1/2}_n = \sum_{s=l}^{n+l} \tilde{g}^{nl}_s r^s,
    """
    N = norm_harmonicBasis(n, l, r0)
    coeffs = np.zeros(2*n+l+1)
    coeffs[l::2] = genlaguerre(n, l+0.5).coeffs[::-1]
    coeffs *= N
    return coeffs
    
def eval_harmonicBasis(n, l, r, r0=1.):
    """ 
    For testing purposes 
    
    Evaluates the harmonic basis function at positions r
    """
    N = norm_harmonicBasis(n, l, r0)
    return (N*r**l* exp(-0.5*r**2/r0**2) * 
            eval_genlaguerre(n, l+0.5, r**2/r0**2))

def calcThetaPhiR(pos):
    pos = np.atleast_2d(pos)
    X, Y, Z = pos.T
    R = norm(pos, axis=1)
    phi = np.arctan2(Y, X)
    theta = np.arccos(Z/R)
    return theta, phi, R
    
def EulerM(a, b, y):
    sina, cosa = sin(a), cos(a)
    sinb, cosb = sin(b), cos(b)
    siny, cosy = sin(y), cos(y)
    Ma = np.array(((cosa, -sina, 0),
                   (sina, cosa, 0),
                   (0,0,1)))
    Mb = np.array(((cosb, 0, -sinb),
                   (0, 1, 0),
                   (sinb, 0, cosb)))
    My = np.array(((cosy, -siny, 0),
                   (siny, cosy, 0),
                   (0,0,1)))
    return My.dot(Mb).dot(Ma)
    
def angle_axis2mat(vector):
    ''' Rotation matrix of angle `theta` around `vector`
    Parameters
    ----------
    vector : 3 element sequence
       vector specifying axis for rotation. Norm of vector gives angle of
       rotation.
    Returns
    -------
    mat : array shape (3,3)
       rotation matrix for specified rotation
    Notes
    -----
    From: https://en.wikipedia.org/wiki/Rotation_matrix#Axis_and_angle
    '''
    vector = np.asanyarray(vector)
    theta = norm(vector)
    if theta==0.:
        return np.eye(3)
    x, y, z = vector/theta
    c, s = cos(theta), sin(theta)
    C = 1 - c
    xs, ys, zs = x * s, y * s, z * s
    xC, yC, zC = x * C, y * C, z * C
    xyC, yzC, zxC = x * yC, y * zC, z * xC
    return np.array([[x * xC + c, xyC - zs, zxC + ys],
                     [xyC + zs, y * yC + c, yzC - xs],
                     [zxC - ys, yzC + xs, z * zC + c]])
                     
def mat2angle_axis(M):
    ''' Calculates rotation vector where the norm of the rotation vector
    indicates the angle of rotation from a rotation matrix M
    
    Parameters
    ----------
    M : (3,3) array like
        matrix encoding rotation matrix
        
    Returns
    -------
    v: array shape (3)
        rotation vector
    '''
    M = np.asanyarray(M)
    theta = np.arccos(0.5*np.trace(M)-0.5)
    v = np.array([M[2,1]-M[1,2],M[0,2]-M[2,0],M[1,0]-M[0,1]])
    v *= 0.5*theta/sin(theta)
    return v
    
class PriorityQueueHeap(object):
    def __init__(self):
        self.heap = []
        
    def get(self):
        return heappop(self.heap)
        
    def put(self, item):
        heappush(self.heap, item)
        
    def empty(self):
        return len(self.heap) == 0
        
    def qsize(self):
        return len(self.heap)