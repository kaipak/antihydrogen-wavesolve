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
PREC       = 128
A1=0.67820882
A2=0.14058198
B1=0.63572127
B2=1.14989192
G1=-0.17638067
G2=0.39904493
NSIZE      = 1 # number of terms
ETA        = 1 - 8.439*(10**-6)
Z_PROTON   = 1

# Options related to COBYLA function.
MAXFUN     = 10000
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
    print('iterations:' + str(iterations))
    constraints = []
    alpha_ks = []
    beta_ks  = []
    gamma_ks = []
    for i in range(1, NSIZE + 1):
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
    global A1
    global A2
    global B1
    global B2
    global G1
    global G2

    description="""Iteratively optimize parallelotope parameters through \
    application of COBYLA algorithm that minimizes energy of wavesolve.\
    Defaults are generally provided for each option as described below--\
    including for the paralleltope params: \n\
    A1=0.67820882,  A2=0.14058198\n \
    B1=0.63572127,  B2=1.14989192\n \
    G1=-0.17638067, G2=0.39904493\n \
    \n\
    Example usage:\
    \n\
    calculation on helium, 10000 iterations, A1-G2 defined):\n\
    python main_cobyla.py -z 2 -m 10000 .7 .15 .65 1.2 .2 .4\n
    """

    parser = argparse.ArgumentParser(description=description)
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
                        help="Maximum number of tries COBYLA will attempt to \
                              minimize function before giving up.  Default \
                              is 10000.")
    parser.add_argument("-b", "--rhobeg", type=float,
                        help="Initial changes in parameters COBYLA will try \
                              when minimizing.  Default is increments of .1.")
    parser.add_argument("-e", "--rhoend", type=str,
                        help="Lower bound accuracy in final solution COBYLA \
                              will produce.  Defaults to 1e-8.")
    parser.add_argument("A1", nargs='?', type=float, \
                        help="A1 parallelotope parameter", default = A1)
    parser.add_argument("A2", nargs='?', type=float, \
                        help="A2 parallelotope parameter", default = A2)
    parser.add_argument("B1", nargs='?', type=float, \
                        help="B1 parallelotope parameter", default = B1)
    parser.add_argument("B2", nargs='?', type=float,
                        help="B2 parallelotope parameter", default = B2)
    parser.add_argument("G1", nargs='?', type=float,
                        help="G1 parallelotope parameter", default = G1)
    parser.add_argument("G2", nargs='?', type=float,
                        help="G2 parallelotope parameter", default = G2)
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
    if args.A1:
        A1 = args.A1
    if args.A2:
        A2 = args.A2
    if args.B1:
        B1 = args.B1
    if args.B2:
        B2 = args.B2
    if args.G1:
        G1 = args.G1
    if args.G2:
        G2 = args.G2
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
