# Math libraries
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

# Default physical parameters modifiable by getopts.
PREC       = 128
A1         = mpm.mpf(.2480)
A2         = mpm.mpf(.8270)
B1         = mpm.mpf(.852)
B2         = mpm.mpf(1.1260)
G1         = mpm.mpf(-.0520)
G2         = mpm.mpf(.1050)
NSIZE      = 1 # number of terms
ETA        = 1 - 8.439*(10**-6)
Z_PROTON   = 1

# Options related to COBYLA function.
MAXFUN     = 1000
RHOBEG     = .1
RHOEND     = '1e-8'

def objective(args):
    return wavesolve.solve(args, Z_PROTON, ETA, NSIZE)

# Calculate a coefficient
def n_subK(L1, L2, k, rooted):
    frac = mpm.mpf(((k*(k+1)*np.sqrt(rooted))/2)%1)
    num = mpm.mpf(((L2 - L1)*frac) + L1)
    return num

def get_constraints(args):
    global iterations
    iterations = iterations + 1
    print 'iterations:' + str(iterations)
    constraints = []
    alpha_ks = []
    beta_ks  = []
    gamma_ks = []
    for i in xrange(1, NSIZE + 1):
        # alpha_k + beta_k
        constraints.append(n_subK(args[0], args[1], i, 2) + \
                           n_subK(args[2], args[3], i, 3))
        # alpha_k + gamma_k
        constraints.append(n_subK(args[0], args[1], i, 2) + \
                           n_subK(args[4], args[5], i, 5))
        # beta_k + gamma_k
        constraints.append(n_subK(args[2], args[3], i, 3) + \
                           n_subK(args[4], args[5], i, 5))
    return constraints

def parse_args():
    """Parse command line arguments
    """
    global Z_PROTON
    global NSIZE
    global PREC
    global MAXFUN
    global RHOBEG
    global RHOEND
    parser = argparse.ArgumentParser(description='Iteratively find \
                                                  parallelotope parameters \
                                                  through application of \
                                                  COBYLA algorithm on \
                                                  wavesolve.')
    parser.add_argument("-z", "--protons", type=int,
                        help="Z number of atom to evaluate. Defaults to \
                              hydrogen (Z=1).")
    parser.add_argument("-n", "--nsize", type=int,
                        help="Number of terms wave function will have. \
                              Defaults to 1 term.")
    parser.add_argument("-p", "--precision", type=int,
                        help="Precision of numbers used in calculations. \
                              Defaults to quadruple precision or 128 bits.")
    parser.add_argument("-m", "--maxfun", type=int,
                        help="Maximum number of tries COBYLA will try to \
                              minimize function before giving up.  Default \
                              is 1000.")
    parser.add_argument("-b", "--rhobeg", type=float,
                        help="Initial changes in parameters COBYLA will try \
                              when minimizing.  Default is increments of .1.")
    parser.add_argument("-e", "--rhoend", type=str,
                        help="Lower bound accuracy in final solution COBYLA \
                              will produce.  Defaults to 1e-8.")
    args = parser.parse_args()

    if args.protons:
        Z_PROTON = args.protons
    if args.nsize:
        NSIZE = args.nsize
    if args.precision:
        PREC = args.precision
    if args.maxfun:
        MAXFUN = args.maxfun
    if args.rhobeg:
        RHOBEG = args.rhobeg
    if args.rhoend:
        RHOEND = args.rhoend
    # print args.accumulate(args.integers)

def main():
    parse_args()
    mpm.mp.prec = PREC
    runtime = str(datetime.datetime.now())
    time_start = timeit.default_timer()
    # filename = "cobyla_wavesolve_run" + str(NSIZE) + runtime
    fmin_cobyla(objective, [A1,A2,B1,B2,G1,G2], get_constraints, maxfun=MAXFUN,
                rhobeg=RHOBEG, rhoend=RHOEND)

if __name__ == '__main__':
    main()
