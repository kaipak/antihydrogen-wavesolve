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
A1       = np.longdouble(0.22277)
A2       = np.longdouble(1.58047)
B1       = np.longdouble(0.98603)
B2       = np.longdouble(1.33237)
G1       = np.longdouble(-0.16261)
G2       = np.longdouble(0.76359)
NSIZE    = 8 # number of terms
ETA      = np.longdouble(1 + 1.42*(10**-9))
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
