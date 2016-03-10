import wavesolve
import mpmath as mpm
import numpy as np
import sys
from scipy.optimize import fmin_cobyla

A1 = .2480
A2 = .8270
B1 = .852
B2 = 1.1260
G1 = -.0520
G2 = .1050
NSIZE = 10 # number of terms

def objective(args):
    wavesolve.solve(args, NSIZE)

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
        alpha_ks.append(n_subK(A1, A2, i, 2))
        beta_ks.append(n_subK(B1, B2, i, 3))
        gamma_ks.append(n_subK(G1, G2, i, 5))

    constraints.extend(alpha_ks)
    constraints.extend(beta_ks)
    constraints.extend(gamma_ks)
    for cons in constraints:
        print cons
    return constraints


def constr1(args):
    return args[0] + args[2]

def constr2(args):
    return args[0] + args[4]

def constr3(args):
    return args[0] + args[4]


def main():
    #objective([.1490,1.199,.899,1.18,-.077,.2250])
    #get_constraints(A1, A2, B1, B2, G1, G2)
    sys.exit
    fmin_cobyla(objective, [.300,.8,.8,1.0,-.1,.1], get_constraints, rhobeg=.01, rhoend=.000001)



if __name__ == '__main__':
    main()
