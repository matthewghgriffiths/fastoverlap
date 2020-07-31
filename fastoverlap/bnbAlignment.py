# -*- coding: utf-8 -*-


from itertools import product
from heapq import heappop, heappush
from functools import total_ordering

import numpy as np
from numpy import sin, sqrt, cos, pi
from numpy.linalg import norm

from scipy.spatial import cKDTree

from .utils import (
    angle_axis2mat, mat2angle_axis, PriorityQueueHeap, lap, _make_cost_matrix)

from . import f90
if f90.have_libbnb:
    from .f90 import libbnb
    


class BranchandBoundMaster(object):
    def __init__(self, initialNodes=[], rtol=5e-2, atol=1e-3):
        self.queue = PriorityQueueHeap()
        self.nodes = 0
        self.ncalc = 0
        self.quickcount = 0
        self.rtol = rtol
        self.atol = atol
        self.record = []
        if len(initialNodes):
            self.best = initialNodes[0]
            initialNodes[0].addtoQueue(self.queue)
        for node in initialNodes[1:]:
            if node.upperbound < self.best.upperbound:
                self.best = node
            if node.lowerbound < self.best.upperbound:
                node.addtoQueue(self.queue)

    def run(self, niter=10000, force=False, iprint=0):
        while(not self.queue.empty()):
            node = self.queue.get()[-1]
            self.nodes += 1
            if (node.lowerbound <
                self.best.upperbound - self.rtol*self.best.upperbound - self.atol):
                self.record.append((
                    self.nodes, self.ncalc, self.queue.qsize(), self.quickcount, node.rot, node.width,
                    node.lowerbound, node.upperbound, self.best.upperbound))
                if iprint:
                    if self.nodes%iprint == 0:
                        print((self.record[-1]))
                newnodes = node.branch()
                for new in newnodes:
                    self.ncalc += 1
                    new.calcbounds(self.best.upperbound, force=force)
                    if new.quick:
                        self.quickcount += 1
                    if new.lowerbound < self.best.upperbound:
                        if new.upperbound < self.best.upperbound:
                            self.best = new
                            self.best.iter = self.nodes
                        new.addtoQueue(self.queue)
                if self.nodes > niter:
                    break
        return self.best

    def __call__(self, pos1, pos2, invert=True, niter=np.inf,
                 extra_inputs=None, force=False, iprint=0):

        pos1 = np.array(pos1)
        pos2 = np.array(pos2)
        cm1 = pos1.mean(0)[None,:]
        cm2 = pos2.mean(0)[None,:]
        pos1 -= cm1
        pos2 -= cm2

        rot = np.zeros(3)

        r1 = norm(pos1, axis=1)
        r2 = norm(pos2, axis=1)

        r1tr2 = (r1[:,None]*r2[None,:])
        r1pr2 = r1[:,None]**2 + r2[None,:]**2
        pos1tree = cKDTree(pos1, 10)

        startNode = BranchNodeCluster(rot, 2*np.pi, pos1, pos2,
                               r1, r2, r1tr2, r1pr2, pos1tree)
        startNode.calcbounds(force=True)

        inputs = [startNode]
        if invert:
            inputs.append(BranchNodeCluster(rot, 2*np.pi, pos1, -pos2,
                                     r1, r2, r1tr2, r1pr2, pos1tree))
            inputs[-1].calcbounds(force=True)
        if extra_inputs is not None:
            inputs.extend(extra_inputs)

        self.nodes=0
        self.__init__(inputs, self.rtol, self.atol)
        self.run(niter, force=force, iprint=iprint)

        bestNode = self.best
        return bestNode.upperbound, pos1+cm1, bestNode.rpos2+cm1

    def plot(self):
        import matplotlib.pyplot as plt
        n, ncalc, qsize, quick, vs, widths, lower, upper, best = \
            list(map(np.array, list(zip(*self.record))))

        fig, axes = plt.subplots(3,1, sharex=True)

        axes[0].plot(n, best)
        axes[0].plot(n, lower)
        axes[0].scatter(n, upper, marker='+')

        bestNode = self.best

        axes[1].plot(n, -np.ones_like(n)*pi, lw=1, c='k')
        axes[1].scatter(n, vs[:,0], marker='+', c='r')
        axes[1].plot(n, np.ones_like(n)*bestNode.rot[0], lw=2, c='r')
        axes[1].plot(n, np.ones_like(n)*pi, lw=1, c='k')
        axes[1].scatter(n, vs[:,1] + 2*np.pi, marker='+', c='b')
        axes[1].plot(n, np.ones_like(n)*bestNode.rot[1] + 2*np.pi, lw=2, c='b')
        axes[1].plot(n, np.ones_like(n)*3*pi, lw=1, c='k')
        axes[1].scatter(n, vs[:,2] + 4*np.pi, marker='+', c='g')
        axes[1].plot(n, np.ones_like(n)*bestNode.rot[2] + 4*np.pi, lw=2, c='g')
        axes[1].plot(n, np.ones_like(n)*5*pi, lw=1, c='k')
        axes[1].plot(n, widths + 5*np.pi)
        axes[1].plot(n, np.ones_like(n)*5*np.pi)

        axes[2].plot(n, ncalc)
        axes[2].plot(n, qsize)
        axes[2].plot(n, quick)

        axes[0].set_xlim(min(n)-10, max(n)+10)
        fig.tight_layout(h_pad=0.1)
        return fig, axes

@total_ordering
class BranchNode(object):
    def addtoQueue(self, queue):
        queue.put((self.lowerbound, self))
    def branch(self):
        raise NotImplementedError
    def calcbounds(self, bestupper=np.inf):
        raise NotImplementedError
        
    def __eq__(self, other):
        return self.lowerbound == other.lowerbound 
    
    def __gt__(self, other):
        return self.lowerbound > other.lowerbound


class BranchNodeBulk(BranchNode):
    branchvects = np.array(list(product([1,-1],repeat=3)))
    def __init__(self, disp, width, pos1, pos2, boxvec, perm):
        self.disp = disp
        self.width = width
        self.pos1, self.pos2 = pos1, pos2
        self.boxvec = boxvec
        self.perm = perm

    def branch(self):
        newdisps = self.branchvects*self.width*0.25
        newnodes = [BranchNodeBulk(disp, self.width*0.5, self.pos1, self.pos2,
                                   self.boxvec, self.perm)
                    for disp in newdisps]
        return newnodes

    def calcbounds(self, bestupper=np.inf, force=False):
        tpos2 = self.pos2 - np.atleast_2d(self.disp)
        disps = [self.pos1[p[:,None],:] - tpos2[p[None,:],:]
                 for p in self.perm]
        for d in disps:
            d -= np.round(d/self.boxvec) * self.boxvec
        dists = [np.norm(d, axis=2) for d in disps]
        ldists = [np.clip(d-0.5*self.width*sqrt(3),0,np.inf)**2 for d in dists]
        dists = [d**2 for d in dists]
        if not force:
            upperbound = sqrt(sum(d.min(0).sum() for d in dists))
        else:
            upperbound = 0.

        if upperbound < bestupper or force:
            self.quick = False
            lowerperm = [lap(ld) for ld in ldists]
            self.lowerbound = sqrt(
                sum(ldists[p,lp] for p,lp in zip(self.perm, lowerperm)))
            if self.lowerbound < bestupper:
                upperperm = [lap(d) for d in dists]
                self.upperbound = sqrt(
                    sum(dists[p,lp] for p,lp in zip(self.perm, upperperm)))
            else:
                upperperm = lowerperm
                self.upperbound = np.inf
        else:
            self.quick = True
            self.lowerbound = sqrt(sum(d.min(0).sum() for d in ldists))
            upperperm = None
        return self.lowerbound, self.upperbound, upperperm

class BranchNodeCluster(BranchNode):
    branchvects = np.array(list(product([1,-1],repeat=3)))
    def __init__(self, rot, width, pos1, pos2,
                 r1, r2, r1tr2, r1pr2, pos1tree):
        self.rot = rot
        self.width = width
        self.pos1, self.pos2 = pos1, pos2
        #self.rpos2 = self.pos2.dot(angle_axis2mat(self.rot))
        self.r1, self.r2 = r1, r2
        self.r1tr2 = r1tr2
        self.r1pr2 = r1pr2
        self.pos1tree = pos1tree
        self.n = 2

    def branch(self):
        n = self.n
        n = 2
        newrots = np.indices((n,n,n), dtype=float).reshape(3,-1).T
        newrots *= 2
        newrots -= (n-1.)
        newrots *= 0.5 * self.width / n
        newrots += self.rot[None,:]
        newnodes = [BranchNodeCluster(rot, self.width/n,
                               self.pos1, self.pos2, self.r1, self.r2,
                               self.r1tr2, self.r1pr2, self.pos1tree)
                    for rot in newrots]
        return newnodes

    def calcbounds(self, bestupper=np.inf, force=False):
        rpos2 = self.pos2.dot(angle_axis2mat(self.rot))

        if not force:
            d, p = self.pos1tree.query(rpos2)
            upperbound = norm(d)
        else:
            upperbound = 0

        if upperbound < bestupper or force:
            self.quick = False
            dists = _make_cost_matrix(self.pos1, rpos2)

            cosa = 0.5*(self.r1pr2 - dists)/self.r1tr2
            sina = sqrt(np.clip(1.-cosa**2,0.,1.))
            theta = min((sqrt(3) * self.width / 2),np.pi)
            cosd = cos(theta)
            sind = abs(sin(theta))
            cosdm = np.where(cosa>cosd, 1., cosa*cosd + sina*sind)

            lowerbounddists = np.clip(self.r1pr2 - 2*self.r1tr2*cosdm, 0., np.inf)
            lowerperm = lap(lowerbounddists.copy())
            self.lowerbound = sqrt(sum(lowerbounddists[i,j] for i,j in enumerate(lowerperm)))

            if self.lowerbound < bestupper:
                upperperm = lap(dists)
                self.upperbound = norm(self.pos1 - rpos2[upperperm])
            else:
                upperperm = lowerperm
                self.upperbound = np.inf
        else:
            self.quick = True
            cosa = 0.5*(self.r1[p]**2 + self.r2**2 - d**2)/self.r1[p]/self.r2
            sina = sqrt(np.clip(1.-cosa**2,0.,1.))
            theta = min((sqrt(3) * self.width / 2),np.pi)
            cosd = cos(theta)
            sind = abs(sin(theta))
            cosdm = np.where(cosa>cosd, 1., cosa*cosd + sina*sind)
            lowerbounddists = np.clip(self.r1[p]**2 + self.r2**2 - 2*self.r1[p]*self.r2*cosdm, 0., np.inf)
            lowerbound = lowerbounddists.sum()**0.5
            self.upperbound = upperbound
            self.lowerbound = lowerbound
            upperperm = p

        self.rpos2 = rpos2
        return self.lowerbound, self.upperbound, upperperm

class BranchnBoundAlignment(object):

    def __init__(self,invert=True,boxSize=None,perm=None):
        self.invert = invert
        self.libbnb = libbnb
        self.gopermdist = libbnb.gopermdist
        self.qlen = self.gopermdist.queuelen
        self.commons = libbnb.commons
        if boxSize is None:
            self.bulk = False
            self.boxvec = np.zeros(3)
            self.gopermdist.setcluster(invert)
        else:
            self.bulk = True
            self.boxvec = np.array(boxSize, dtype=float)
            self.gopermdist.setbulk(invert)
        if perm is not None:
            self.setPerm(perm)
        self.Natoms = None

    def setPerm(self, perm):
        self.Natoms = sum(map(len,perm))
        self.perm = perm
        self.nperm = len(perm)
        self.npermsize = list(map(len, perm))
        self.permgroup = np.concatenate([np.asanyarray(p)+1 for p in perm])
        self.gopermdist.setperm(self.Natoms, self.permgroup, self.npermsize)

    def initialise(self, pos1, pos2, perm=None, debug=False):
        pos1 = np.asanyarray(pos1).reshape(-1,3)
        pos2 = np.asanyarray(pos2).reshape(-1,3)

        if perm is not None:
            self.setPerm(perm)
        elif len(pos1) != self.Natoms:
            self.Natoms = len(pos1)
            self.setPerm([np.arange(self.Natoms)])

        self.coordsb = np.array(pos1).flatten()
        self.coordsa = np.array(pos2).flatten()
        self.gopermdist.initialise(
            self.coordsb, self.coordsa, self.boxvec[0], self.boxvec[1],
            self.boxvec[2], self.bulk)
        self.libbnb.alignutils.debug = debug

    def __call__(self, pos1, pos2, perm=None, invert=None, debug=False,
                 force=False, niter=1000, iprint=1):
        if invert is None:
            invert = self.invert
        if invert:
            self.invert = invert
            if self.bulk:
                self.gopermdist.setbulk(invert)
            else:
                self.gopermdist.setcluster(invert)
        self.initialise(pos1, pos2, perm=perm, debug=debug)
        self.libbnb.alignutils.nstored = 0
        bestupper = np.array(np.inf)
        if self.bulk:
            width = max(self.boxvec)
        else:
            width = 2*pi
        self.gopermdist.addnode(np.zeros(3),width,1,bestupper,True)
        if self.bulk and self.invert:
            for i in range(2,49):
                self.gopermdist.addnode(np.zeros(3),width,i,bestupper,force)
        elif self.invert:
            self.gopermdist.addnode(np.zeros(3),width,2,bestupper,force)
        self.gopermdist.run(niter,force,iprint,bestupper)
        bestid = self.gopermdist.bestid.item()-1
        coordsb = self.gopermdist.savecoordsb.reshape(pos1.shape)
        coordsa = self.gopermdist.bestcoordsa[:,bestid].reshape(pos2.shape)
        if self.bulk:
            return bestupper.item(), coordsb, coordsa
        else:
            rmat = self.gopermdist.bestrmat[:,:,bestid]
            return bestupper.item(), coordsb, coordsa, rmat


if __name__ == "__main__":
    import os
    import csv

    # Turn debug on if you want status messages
    # You will get A LOT of print statements!
    debug=False
    # Turn periodic alignment on, much slower...
    periodic = True
    # Turn octahedral symmetry searching on, 48 times slower!
    ohcell = False #True


    def readFile(filename):
        with open(filename, 'rb') as f:
            reader = csv.reader(f, delimiter=' ')
            dist = [list(map(float, row)) for row in reader]
        return np.array(dist)

    datafolder = "../examples/LJ38"
    natoms = 38
    pos1 = readFile(os.path.join(datafolder, 'coords'))
    pos2 = readFile(os.path.join(datafolder, 'finish'))

    bnbpy = BranchandBoundMaster()


    if not f90.have_libbnb:
        # Much slower!
        dpyclus, coordsb, coordsa = bnbpy(pos1, pos2)
        fig, axes = bnbpy.plot()

    if f90.have_libbnb:
        bnbcluster = BranchnBoundAlignment()
        dcluster, coordsb, coordsa, rmat = bnbcluster(pos1, pos2, debug=debug, niter=100)

    print('Cluster alignment:')
    print('On example LJ38 data, distance should = 1.4767')
    if not f90.have_libbnb:
        print(('Branch and bound python alignment: ', dpyclus))
    if f90.have_libbnb:
        print(('Branch and bound Fortran alignment: ', dcluster))

    if f90.have_libbnb and periodic:
        datafolder = "../examples/BLJ256"

        pos1 = readFile(os.path.join(datafolder, 'coords'))
        pos2 = readFile(os.path.join(datafolder, 'finish'))

        natoms = 256
        ntypeA = 204
        shape = (natoms, 3)
        boxSize = np.ones(3)*5.975206329
        permlist = [np.arange(ntypeA), np.arange(ntypeA, natoms)]
        bnbbulk = BranchnBoundAlignment(invert=False, boxSize=boxSize)

        dbulk, coordsab, coordsa = bnbbulk(pos1, pos2, debug=debug, niter=50)

        print('\nPeriodic alignment:')
        print('On example BLJ256 data, distance should = 1.559')
        print(('Branch and bound alignment: ', dbulk))

        # Testing for octahderal symetries will take ~48 times longer!

        if ohcell:
            bnbbulk = BranchnBoundAlignment(boxSize=boxSize)
            dbulk, coordsab, coordsa = bnbbulk(pos1, pos2.dot([[0,1,0],[1,0,0],[0,0,1]]),
                                               debug=False, niter=2000, invert=True)
            print(("Octahedral alignment: ", dbulk))
