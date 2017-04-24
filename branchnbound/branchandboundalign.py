# -*- coding: utf-8 -*-


import numpy as np
from numpy import pi

import libbnb

class BranchnBoundAlignment(object):
    
    def __init__(self,invert=True,boxSize=None):
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
        self.Natoms = None

    def setPerm(self, perm):
        self.Natoms = sum(map(len,perm))
        self.perm = perm
        self.nperm = len(perm)
        self.npermsize = map(len, perm)
        self.permgroup = np.concatenate([np.asanyarray(p)+1 for p in perm])
        self.gopermdist.setperm(self.Natoms, self.permgroup, self.npermsize)

    def initialise(self, pos1, pos2, perm=None, debug=False):
        if perm is not None:
            self.setPerm(perm)
        elif len(pos1) != self.Natoms:
            self.Natoms = len(pos1)
            self.setPerm([np.arange(self.Natoms)])
            
        self.coordsb = np.asanyarray(pos1).flatten()
        self.coordsa = np.asanyarray(pos2).flatten()
        self.gopermdist.initialise(
            self.coordsb, self.coordsa, self.boxvec[0], self.boxvec[1], 
            self.boxvec[2], self.bulk)
        self.gopermdist.debug=debug
            
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
        bestupper = np.array(np.inf)
        if self.bulk:
            width = max(self.boxvec)
        else:
            width = 2*pi
        self.gopermdist.addnode(np.zeros(3),width,1,bestupper,True)
        if self.bulk and self.invert:
            for i in xrange(2,49):
                self.gopermdist.addnode(np.zeros(3),width,i,bestupper,force)
        elif self.invert:
            self.gopermdist.addnode(np.zeros(3),width,2,bestupper,force)
        self.gopermdist.run(niter,force,iprint,bestupper)
        bestid = self.gopermdist.bestid.item()-1
        coordsb = self.gopermdist.savecoordsb.reshape(pos1.shape)
        coordsa = self.gopermdist.savecoordsa[:,bestid].reshape(pos2.shape)
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
        
    datafolder = "../examples/LJ38"
    def readFile(filename):
        with open(filename, 'rb') as f:
            reader = csv.reader(f, delimiter=' ')
            dist = [map(float, row) for row in reader]
        return np.array(dist)
    
    pos1 = readFile(os.path.join(datafolder, 'coords'))
    pos2 = readFile(os.path.join(datafolder, 'finish'))
    
    bnbcluster = BranchnBoundAlignment()
    
    natoms = 38
    
    dcluster, coordsb, coordsa, rmat = bnbcluster(pos1, pos2, debug=False, niter=1e6)
                                    
    datafolder = "../examples/BLJ256"

    pos1 = readFile(os.path.join(datafolder, 'coords'))
    pos2 = readFile(os.path.join(datafolder, 'finish'))

    natoms = 256
    ntypeA = 204
    shape = (natoms, 3)
    boxSize = np.ones(3)*5.975206329
    permlist = [np.arange(ntypeA), np.arange(ntypeA, natoms)]

    bnbbulk = BranchnBoundAlignment(invert=False, boxSize=boxSize)
    
    # Testing for octahderal symetries will take ~48 times longer!
    dbulk, coordsab, coordsa = bnbbulk(pos1, pos2, debug=False, niter=1e6)
                              
    print 'Summary:'
    print 'Cluster alignment:'
    print 'On example LJ38 data, distance should = 1.4767'                            
    print 'Branch and bound alignment: ', dcluster
           
    print '\nPeriodic alignment:'
    print 'On example BLJ256 data, distance should = 1.559'                            
    print 'Branch and bound alignment: ', dbulk