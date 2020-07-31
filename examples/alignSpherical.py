# -*- coding: utf-8 -*-
import os
import csv
import numpy as np
from fastoverlap import SphericalAlign

try:
    from fastoverlap import BranchnBoundAlignment
    alignbnb = True
except ImportError:
    alignbnb = False

datafolder = "LJ38/"

def readFile(filename):
    with open(filename, 'rb') as f:
        reader = csv.reader(f, delimiter=' ')
        dist = [list(map(float, row)) for row in reader]
    return np.array(dist)

pos1 = np.loadtxt(os.path.join(datafolder, 'coords'))
pos2 = np.loadtxt(os.path.join(datafolder, 'finish'))

natoms = 38
scale = 0.3 # Set this to be ~ half interatomic separation
maxl = 15 # Max degree of spherical harmonic

align = SphericalAlign(scale, maxl)

if alignbnb:
    bnb = BranchnBoundAlignment()



if __name__ == "__main__":

    print('Performing permutational alignment with Hungarian algorithm')
    permRMS = align.Hungarian(pos1, pos2)[0]*natoms**-0.5
    print('RMSD = {:0.4f}'.format(permRMS))

    print('Performing fastoverlap alignment')
    fastdist = align(pos1, pos2)[0]
    fastRMS = fastdist*natoms**-0.5
    print('RMSD = {:0.4f}'.format(fastRMS))
    if alignbnb:
        print('Performing branch and bound alignment')
        bnbdist = bnb(pos1, pos2)[0]
        bnbRMS = bnbdist*natoms**-0.5
        print('RMSD = {:0.4f}'.format(bnbRMS))

    import timeit
    print('Timing fastoverlap alignment:')
    alignTimer = timeit.Timer(stmt="salign.align(salign.pos1, salign.pos2)",
                     setup="import alignSpherical as salign")
    aligntime =  alignTimer.timeit(10)/10.

    print('Average time to align for fastoverlap {:0.3} s'.format(aligntime))

    if alignbnb:
        bnbTimer = timeit.Timer(stmt="salign.bnb(salign.pos1, salign.pos2)",
                         setup="import alignSpherical as salign")
        bnbtime =  bnbTimer.timeit(10)/10.
        print('Average time to align for branch and bound alignment {:0.3} s'.format(bnbtime))


