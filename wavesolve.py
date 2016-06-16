"""
  Generalized wavesolve program for use with Powell COBYLA methodology of
  iteratively guessing correct parallelotrope parameters

  --author - Kai Pak
  --start_date - February 01, 2016

"""

# Standard libraries
import datetime
import itertools
import timeit
import multiprocessing as mp
import re
import mpmath as mpm # for matrix math in arbitrary precision
import numpy as np
import sys

# Custom libraries
import ws_maths
import ws_physics


def solve(args, z_proton, eta, nsize):
    '''
    A1-G2 parameters passed as an array due the data type accepted by the python
    COBYLA wrapper.
    '''
    A1 = args[0]
    A2 = args[1]
    B1 = args[2]
    B2 = args[3]
    G1 = args[4]
    G2 = args[5]
    # Set up support stuff
    runtime = str(datetime.datetime.now())
    time_start = timeit.default_timer()
    # filename = "runwavesolve_n" + str(NSIZE) + "_prec" + str(PREC) + "_" + runtime
    # f = open(filename, 'w')

    # Physical parameters
    ws_physics.static_params(A1, A2, B1, B2, G1, G2, eta, z_proton, nsize)
    alphas = ws_physics.thakar_smith_param(A1, A2, 2)
    betas = ws_physics.thakar_smith_param(B1, B2, 3)
    gammas = ws_physics.thakar_smith_param(G1, G2, 5)

    # Create A Matrix
    time_matA_start = timeit.default_timer()
    a_mat = mpm.matrix(nsize)
    a_mat_np = np.zeros(shape=(nsize,nsize))
    print a_mat_np
    # Iterate through matrix starting at index 0
    for m in xrange(0, nsize):
        pool = mp.Pool(processes = None)
        rowobj = [pool.apply_async(ws_physics.thakar_smith_amat,
                                   args=(alphas[n], betas[n], gammas[n],
                                         alphas[m], betas[m], gammas[m]))
                  for n in xrange(m, nsize)]
        row = [p.get() for p in rowobj]
        pool.terminate()

        # Now populate row and its transpose column since this is a symmetric matrix
        for n in xrange(0, len(row)):
            a_mat_np[m,m+n] = row[n]
            a_mat_np[m+n,m] = row[n]
        #print "Done with a_mat, row", m
    time_matA_end = timeit.default_timer()
    print a_mat_np


    # Create B Matrix
    time_matB_start = timeit.default_timer()
    b_mat = mpm.matrix(nsize)
    for m in xrange(0, nsize):
        pool = mp.Pool(processes = None)
        rowobj = [pool.apply_async(ws_physics.thakar_smith_bmat,
                                   args=(alphas[n], betas[n], gammas[n],
                                         alphas[m], betas[m], gammas[m]))
                  for n in xrange(m, nsize)]
        row = [p.get() for p in rowobj]
        pool.terminate()

        # Now populate row and its transpose column since this is a symmetric matrix
        for n in xrange(0, len(row)):
            b_mat[m,m+n] = row[n]
            b_mat[m+n,m] = row[n]
        #print "Done with b_mat, row", m
    time_matB_end = timeit.default_timer()

    #ui_mat, eigvals, eigvecs = ws_maths.eigensolve(a_mat, b_mat)
    #zn_mat, energy, coeff    = ws_maths.normalize_Z(ui_mat, eigvecs, eigvals)

    time_end = timeit.default_timer()

    # Output to file and display times
    #f.write("Matrix UI:================================================\n")
    #f.write(matrix_format(ui_mat))
    #f.write("\n\nEigenvalues:==========================================\n")
    #f.write(matrix_format(eigvals))
    #f.write("\n\nEigenvectors:=========================================\n")
    #f.write(matrix_format(eigvecs))
    #f.write("\n\nMatrix Zn:============================================\n")
    #f.write(matrix_format(zn_mat))
    #f.write("\n\nLowest Eigenvalue:====================================\n")
    #f.write(str(energy))
    #f.write("\n\nCorresponding Coeffs:=================================\n")
    #f.write(matrix_format(coeff))

    matA_time = time_matA_end - time_matA_start
    matB_time = time_matB_end - time_matB_start
    total_time = time_end - time_start

    #print "\n\nmat_A generation time: ", matA_time
    #f.write("\n\nmat_A generation time: " + str(matA_time))
    #print "mat_B generation time: ", matB_time
    #f.write("\n\nmat_B generation time: " + str(matB_time))
    print "Total execution time: ", total_time
    #f.write("\n\nTotal execution time: " + str(total_time))
    #f.close()
    # return energy


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

def matrix_format(matrix):
    """Change formatting to be compatible with Mathematica natively"""
    matrix_text = str(matrix).split('\n')
    matrix_return = ''
    for line in matrix_text:
        line = re.sub('\[', '{', line)
        line = re.sub('\]', '},', line)
        line = re.sub(r'(\.\d+)(\s+)', r'\1, ', line)
        line = re.sub('mpf\(\'', '', line)
        line = re.sub(r'\'\)', '', line)
        line = re.sub(r'\,\}', '}', line)
        line = re.sub(r'\,\'\)', ' ', line)
        matrix_return += line + '\n'
    return matrix_return
