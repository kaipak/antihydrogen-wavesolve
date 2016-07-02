# Math and science libraries
import wavesolve
import mpmath as mpm
import numpy as np
from scipy.optimize import fmin_cobyla
# System libraries
import argparse
import sys
import datetime
import itertools
import timeit

iterations = 0

# Default physical parameters modifiable by argparse
A1       = 0.67820882
A2       = 0.14058198
B1       = 0.63572127
B2       = 1.14989192
G1       = -0.17638067
G2       = 0.39904493
NSIZE    = 32 # number of terms
ETA      = 1 - 8.439*(10**-6)
Z_PROTON = 1

def main():
    runtime = str(datetime.datetime.now())
    time_start = timeit.default_timer()
    # filename = "cobyla_wavesolve_run" + str(NSIZE) + runtime
    np.finfo(np.float64)
    print np.float128(.9918237412347219341284219)
    wavesolve.solve([A1,A2,B1,B2,G1,G2], Z_PROTON, ETA, NSIZE)

if __name__ == '__main__':
    main()
