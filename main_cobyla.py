# Math libraries
import wavesolve
import mpmath as mpm
import numpy as np
from scipy.optimize import fmin_cobyla
# System libraries
import sys
import datetime
import itertools
import timeit
# -0.5277304493274347`+0.527730449
iterations = 0
PREC = 128
A1 = mpm.mpf(.2480)
A2 = mpm.mpf(.8270)
B1 = mpm.mpf(.852)
B2 = mpm.mpf(1.1260)
G1 = mpm.mpf(-.0520)
G2 = mpm.mpf(.1050)
NSIZE = 100 # number of terms

def objective(args):
    return wavesolve.solve(args, NSIZE)

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
        constraints.append(n_subK(args[0], args[1], i, 2) + n_subK(args[2], args[3], i, 3))
        # alpha_k + gamma_k
        constraints.append(n_subK(args[0], args[1], i, 2) + n_subK(args[4], args[5], i, 5))
        # beta_k + gamma_k
        constraints.append(n_subK(args[2], args[3], i, 3) + n_subK(args[4], args[5], i, 5))
    return constraints

def main():
    mpm.mp.prec = PREC
    runtime = str(datetime.datetime.now())
    time_start = timeit.default_timer()
    # mpm.mp.prec = PREC
    # filename = "cobyla_wavesolve_run" + str(NSIZE) + runtime
    fmin_cobyla(objective, [A1,A2,B1,B2,G1,G2], get_constraints, maxfun=100000, rhobeg=.1, rhoend=1e-16)

if __name__ == '__main__':
    main()
