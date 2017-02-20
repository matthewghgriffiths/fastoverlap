# -*- coding: utf-8 -*-
import os
import csv
import numpy as np
from fastoverlap import SphericalAlign, SphericalHarmonicAlign

datafolder = "BLJ256/"

def readFile(filename):
    with open(filename, 'rb') as f:
        reader = csv.reader(f, delimiter=' ')
        dist = [map(float, row) for row in reader]
    return np.array(dist)
    
pos1 = readFile(os.path.join(datafolder, 'coords'))
pos2 = readFile(os.path.join(datafolder, 'finish'))

natoms = 38
scale = 0.5 # Set this to be ~ half interatomic separation
maxl = 15 # Max degree of spherical harmonic
maxn = 7 # Number of harmonic basis functions to use.

align = SphericalAlign(scale, maxl)
harm = SphericalHarmonicAlign(scale, maxn, maxl)

permRMS = align.Hungarian(pos1, pos2)[0]*natoms**-0.5

if __name__ == "__main__":
    print 'Performing fastoverlap alignment'
    fastdist = align(pos1, pos2)[0]
    fastRMS = fastdist*natoms**-0.5
    print 'RMSD = {:0.4f}'.format(fastRMS)
    
    print 'Performing fastoverlap alignmentin harmonic basis'
    harmdist = harm(pos1, pos2)[0]
    harmRMS = fastdist*natoms**-0.5
    print 'RMSD = {:0.4f}'.format(harmRMS)
    
    import timeit
    print 'Timing fastoverlap alignment:'
    alignTimer = timeit.Timer(stmt="salign.align(salign.pos1, salign.pos2)", 
                     setup="import alignSpherical as salign")  
    aligntime =  alignTimer.timeit(10)/10.
    harmTimer = timeit.Timer(stmt="salign.harm(salign.pos1, salign.pos2)", 
                     setup="import alignSpherical as salign")  
    harmtime =  harmTimer.timeit(10)/10.
        
    print 'Average time to align for fastoverlap {:0.3} s'.format(aligntime)
    print 'Average time to align for harmonics basis {:0.3} s'.format(harmtime)
     