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
from numpy import array, matrix, linalg, pi, set_printoptions, modf
from IPython.display import display
import sympy as sym

# Custom libraries
import ws_maths
import ws_physics

# Physical parameters
A1 = .1280
A2 = 1.2460
B1 = .8890
B2 = 1.2030
G1 = -.1100
G2 = .3420
Z  = 1
ETA = 1 + 2.314*(10**-9)

# Application attributes
NSIZE = 60
PREC = 16

set_printoptions(precision=PREC)

def main():
    start_time = timeit.default_timer()
    ws_physics.static_params(A1, A2, B1, B2, G1, G2, ETA, Z, NSIZE, PREC)
    alphas = ws_physics.thakar_smith_param(A1, A2, 2)
    betas = ws_physics.thakar_smith_param(B1, B2, 3)
    gammas = ws_physics.thakar_smith_param(G1, G2, 5) 
    print alphas


    """Create A Matrix"""
    a_mat = [[0.0] * NSIZE for i in xrange(NSIZE)]
    # Iterate through matrix starting at index 0
    for m in xrange(0, NSIZE):
        for n in range(0, NSIZE):
            a_mat[m][n] = float(ws_physics.psif_H_psii(alphas[n], betas[n], gammas[n], \
                                         alphas[m], betas[m], gammas[m]) +\
                  ws_physics.psif_H_psii(betas[n], alphas[n], gammas[n], \
                                         alphas[m], betas[m], gammas[m]) +\
                  ws_physics.psif_H_psii(alphas[n], betas[n], gammas[n],\
                                         betas[m], alphas[m], gammas[m]) +\
                  ws_physics.psif_H_psii(betas[n], alphas[n], gammas[n], \
                                         betas[m], alphas[m], gammas[m]))
        print "Done with a_mat, row", m

    """Create B Matrix"""
    b_mat = [[0.0] * NSIZE for i in xrange(NSIZE)]
    for m in xrange(0, NSIZE):
        for n in range(0, NSIZE):
            alphanm = alphas[n] + alphas[m]
            betanm = betas[n] + betas[m]
            gammanm = gammas[n] + gammas[m]
            b_mat[m][n] = float(8 * sym.pi**2 * (ws_physics.hfs_gamma(1, 1, 1, alphanm, betanm, gammanm) +
                                                 ws_physics.hfs_gamma(1, 1, 1, alphas[n] + betas[m], alphas[m] + betas[n], gammas[n] + gammas[m]) +
                                                 ws_physics.hfs_gamma(1, 1, 1, alphas[m] + betas[n], alphas[n] + betas[m], gammas[n] + gammas[m]) +
                                                 ws_physics.hfs_gamma(1, 1, 1, betas[n] + betas[m], alphas[n] + alphas[m], gammas[n] + gammas[m])))
        print "Done with b_mat, row", m

    ws_maths.eigensolve(a_mat, b_mat)

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
