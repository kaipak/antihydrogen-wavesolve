"""
  Modified main for Jack solving explicitly correlated wavefunction.
  200 terms
  
  --author - Kai Pak
  --start_date - May 22, 2015
  
"""

# Standard libraries
import datetime
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
# P10
#a1 = .2480
#a2 = .8270
#b1 = .852
#b2 = 1.1260
#g1 = -.0520
#g2 = .1050
#eTA = 1 - 8.439*(10**-6)

## P60
#A1 = .1280
#A2 = 1.2460
#B1 = .8890
#B2 = 1.2030
#G1 = -.1100
#G2 = .3420
#ETA = 1 + 2.314*(10**-9)

# P200
A1 = .22277
A2 = 1.58047
B1 = .98603
B2 = 1.33237
G1 = -.16261
G2 = .76359
ETA = 1 + 1.42*(10**-9)

Z  = 1

# Application attributes
NSIZE = 200
PREC = 128 # in bits
DPS = 256 # decimal places

# set_printoptions(precision=PREC)

def main():
    # Set up support stuff
    runtime = str(datetime.datetime.now())
    time_start = timeit.default_timer()
    mpm.mp.prec = PREC
    filename = "runwavesolve_n" + str(NSIZE) + "_prec" + str(PREC) + "_" + runtime
    f = open(filename, 'w')

    # Physical parameters
    ws_physics.static_params(A1, A2, B1, B2, G1, G2, ETA, Z, NSIZE, PREC, DPS)
    alphas = ws_physics.thakar_smith_param(A1, A2, 2)
    betas = ws_physics.thakar_smith_param(B1, B2, 3)
    gammas = ws_physics.thakar_smith_param(G1, G2, 5) 

    # Create A Matrix
    time_matA_start = timeit.default_timer()
    a_mat = mpm.matrix(NSIZE) 
    # Iterate through matrix starting at index 0
    for m in xrange(0, NSIZE):
        pool = mp.Pool(processes = None)
        rowobj = [pool.apply_async(ws_physics.thakar_smith_amat,
                                   args=(alphas[n], betas[n], gammas[n],
                                         alphas[m], betas[m], gammas[m]))
                  for n in xrange(m, NSIZE)]
        row = [p.get() for p in rowobj]
        pool.terminate()
        
        # Now populate row and its transpose column since this is a symmetric matrix
        for n in xrange(0, len(row)):
            a_mat[m,m+n] = row[n]
            a_mat[m+n,m] = row[n]
        print "Done with a_mat, row", m
    time_matA_end = timeit.default_timer()


    # Create B Matrix
    time_matB_start = timeit.default_timer()
    b_mat = mpm.matrix(NSIZE)
    for m in xrange(0, NSIZE):
        pool = mp.Pool(processes = None)
        rowobj = [pool.apply_async(ws_physics.thakar_smith_bmat,
                                   args=(alphas[n], betas[n], gammas[n],
                                         alphas[m], betas[m], gammas[m]))
                  for n in xrange(m, NSIZE)]
        row = [p.get() for p in rowobj]
        pool.terminate()
        
        # Now populate row and its transpose column since this is a symmetric matrix
        for n in xrange(0, len(row)):
            b_mat[m,m+n] = row[n]
            b_mat[m+n,m] = row[n]
        print "Done with b_mat, row", m
    time_matB_end = timeit.default_timer()

    ui_mat, eigvals, eigvecs = ws_maths.eigensolve(a_mat, b_mat)
    zn_mat, energy, coeff    = ws_maths.normalize_Z(ui_mat, eigvecs, eigvals)

    time_end = timeit.default_timer()

    # Output to file and display times
    f.write("Matrix UI:================================================\n")
    f.write(str(ui_mat))
    f.write("\n\nEigenvalues:==========================================\n")
    f.write(str(eigvals))
    f.write("\n\nEigenvectors:=========================================\n")
    f.write(str(eigvecs))
    f.write("\n\nMatrix Zn:============================================\n")
    f.write(str(zn_mat))
    f.write("\n\nLowest Eigenvalue:====================================\n")
    f.write(str(energy))
    f.write("\n\nCorresponding Coeffs:=================================\n")
    f.write(str(coeff))

    matA_time = time_matA_end - time_matA_start
    matB_time = time_matB_end - time_matB_start
    total_time = time_end - time_start
    
    print "\n\nmat_A generation time: ", matA_time
    f.write("\n\nmat_A generation time: " + str(matA_time))
    print "mat_B generation time: ", matB_time
    f.write("\n\nmat_B generation time: " + str(matB_time))
    print "Total execution time: ", total_time
    f.write("\n\nTotal execution time: " + str(total_time))
    f.close()

    
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
