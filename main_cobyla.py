import wavesolve
import mpmath as mpm
import numpy as np
import sys
from scipy.optimize import fmin_cobyla

A1 = .1
A2 = .5
B1 = .5
B2 = 1.0
G1 = -.1
G2 = .1
NSIZE = 10 # number of terms

def objective(args):
    return wavesolve.solve(args, NSIZE)

# Calculate a coefficient
def n_subK(L1, L2, k, rooted):
    frac = mpm.mpf(((k*(k+1)*np.sqrt(rooted))/2)%1)
    num = mpm.mpf(((L2 - L1)*frac) + L1)
    return num

def get_constraints(args):
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

def main():
    #objective([.1490,1.199,.899,1.18,-.077,.2250])
    #get_constraints(A1, A2, B1, B2, G1, G2)
    fmin_cobyla(objective, [A1,A2,B1,B2,G1,G2], get_constraints, rhobeg=.1, rhoend=.0000000000000001)



if __name__ == '__main__':
    main()
