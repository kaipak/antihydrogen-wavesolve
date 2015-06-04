"""
  Modified main for Jack solving explicitly correlated wavefunction.
  200 terms
  
  --author - Kai Pak
  --start_date - May 22, 2015
  
"""

# Standard libraries
import itertools
import timeit
import multiprocessing as mp
# from numpy import array, matrix, linalg, pi, set_printoptions, modf
from IPython.display import display
import sympy as sym
import mpmath as mpm # for matrix math in arbitrary precision
import sys

# Custom libraries
import ws_maths
import ws_physics

# Physical parameters

A1 = .2480
A2 = .8270
B1 = .852
B2 = 1.1260
G1 = -.0520
G2 = .1050
ETA = 1 - 8.439*(10**-6)

#A1 = .22277
#A2 = 1.58047
#B1 = .98603
#B2 = 1.33237
#G1 = -.16261
#G2 = .76359
#ETA = 1 + 1.42*(10**-9)

Z  = 1

# Application attributes
NSIZE = 10
PREC = 128 # in bits
DPS = 256 # decimal places

# set_printoptions(precision=PREC)

def main():
    time_start = timeit.default_timer()
    mpm.mp.prec = PREC
    ws_physics.static_params(A1, A2, B1, B2, G1, G2, ETA, Z, NSIZE, PREC, DPS)
    alphas = ws_physics.thakar_smith_param(A1, A2, 2)
    betas = ws_physics.thakar_smith_param(B1, B2, 3)
    gammas = ws_physics.thakar_smith_param(G1, G2, 5) 

    print alphas
    print betas
    print gammas


    # sys.exit()

    """Create A Matrix"""
    time_matA_start = timeit.default_timer()
    a_mat = mpm.matrix(NSIZE) 
    # Iterate through matrix starting at index 0
    for m in xrange(0, NSIZE):
        for n in range(0, NSIZE):
            a_mat[m,n] = ws_physics.psif_H_psii(alphas[n], betas[n], gammas[n], \
                                         alphas[m], betas[m], gammas[m]) +\
                  ws_physics.psif_H_psii(betas[n], alphas[n], gammas[n], \
                                         alphas[m], betas[m], gammas[m]) +\
                  ws_physics.psif_H_psii(alphas[n], betas[n], gammas[n],\
                                         betas[m], alphas[m], gammas[m]) +\
                  ws_physics.psif_H_psii(betas[n], alphas[n], gammas[n], \
                                         betas[m], alphas[m], gammas[m])
        print "Done with a_mat, row", m
    time_matA_end = timeit.default_timer()
    print "A Matrix:"
    print a_mat


    """Create B Matrix"""
    time_matB_start = timeit.default_timer()
    b_mat = mpm.matrix(NSIZE)
    for m in xrange(0, NSIZE):
        for n in range(0, NSIZE):
            b_mat[m,n] = mpm.mpf(8 * mpm.pi**2 * 
                                 (ws_physics.hfs_gamma(1, 1, 1, alphas[n] + alphas[m], 
                                                       betas[n] + betas[m], 
                                                       gammas[n] + gammas[m]) +
                                  ws_physics.hfs_gamma(1, 1, 1, alphas[n] + betas[m], 
                                                       alphas[m] + betas[n], 
                                                       gammas[n] + gammas[m]) +
                                  ws_physics.hfs_gamma(1, 1, 1, alphas[m] + betas[n], 
                                                       alphas[n] + betas[m], 
                                                       gammas[n] + gammas[m]) +
                                  ws_physics.hfs_gamma(1, 1, 1, betas[n] + betas[m], 
                                                       alphas[n] + alphas[m], 
                                                       gammas[n] + gammas[m])))
        print "Done with b_mat, row", m
    time_matB_end = timeit.default_timer()
    print "B Matrix:"
    print b_mat


    ws_maths.eigensolve2(a_mat, b_mat)

    time_end = timeit.default_timer()

    print "\n\nmat_A generation time: ", time_matA_end - time_matA_start
    print "mat_B generation time: ", time_matB_end - time_matB_start
    print "Total execution time: ", time_end - time_start


def build_matrix(psis_i, psis_j, bracket_notation):
    """Generate nxn matrix by determining inner product of psi_i and psi_j
    
    Keyword arguments:
    psi_i -- an array of bras (according to Dirac bra-ket notation)
    psi_i -- an array of kets (according to Dirac bra-ket notation)
    bracket_notation -- some string describing what bra-ket value you're
                        trying to construct, e.g., <psi_i|H|psi_j>
                        
    Returns:
    matrix -- an nxn matrix with all calculated inner products
    
     Build matrix in following manner:
    [ n11, n12, n13]
    [      n22, n23]
    [           n33]
   
    As these are symmetric matrices, we save roughly n/2 processing time by
    just mirroring from upper triangular
    
    Multiprocessing version
    """
    matrix = [[0.0] * NSIZE for i in xrange(NSIZE)]
    
    
    for i in xrange(0, NSIZE):
        print 'Constructing row', i, 'of', bracket_notation
        
        pool = mp.Pool(processes = None)
        psirowobj = [pool.apply_async(ws_physics.get_qstate,
                                      args=(psis_i[i],psis_j[j],
                                            ALPHA, BETA, GAMMA))
                     for j in xrange(i, NSIZE)]
        psirow = [p.get() for p in psirowobj]
        pool.terminate()
        
        # Now populate row and its transpose column since this is a symmetric matrix
        for j in xrange(0, len(psirow)):
            matrix[i][i+j] = psirow[j]
            matrix[i+j][i] = psirow[j]
        
    return matrix
    
if __name__ == '__main__':
    main()
