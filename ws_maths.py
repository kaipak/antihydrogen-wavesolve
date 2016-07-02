"""
    Set of math functions that'll handle linear algebra calculations for
    wavesolve.py program

    @author Kai Pak
    @start_date January 1, 2014
    @current_version 1.1


    History

    July 22 2014 - Separating physics wave function generation into own set of
    functions.

"""

import mpmath as mpm
import sympy as sym
import numpy as np
import scipy as sc
import sys

from IPython.display import display
import timeit
perftime = True

def eigentest(mat_A, mat_B):
    print "B Matrix\n", mat_B
    print "B Matrix\n", mat_B.dtype
    mat_BI = sc.linalg.inv(mat_B)
    mat_H  = mat_A.dot(mat_BI)
    print "H Matrix\n", mat_H.dtype
    la, v  = sc.linalg.eig(mat_H.astype(np.float64))
    print v
    for eigenvalue in la:
        print eigenvalue
    exit

def eigensolve_numpy(mat_A, mat_B):
    ident, lower, upper = sc.linalg.lu(mat_B)
    print "ALU Matrix \n", lower, upper
    matrix_UT  = sc.linalg.cholesky(mat_B, lower=True)
    #matrix_UT = lower
    #print "Chol Matrix \n", matrix_chol
    print "UT Matrix \n", matrix_UT
    matrix_U   = matrix_UT.transpose()
    print "U Matrix \n", matrix_U
    matrix_UI  = sc.linalg.inv(matrix_U)
    print "UI Matrix \n", matrix_UI
    matrix_UIT = matrix_UI.transpose()
    print "UIT Matrix \n", matrix_UIT
    matrix_C = matrix_UIT.dot(mat_A)
    print "C Matrix \n", matrix_C
    matrix_C = matrix_C.dot(matrix_UI)
    print "C Matrix \n", matrix_C

    eigsolve_stime = timeit.default_timer()
    la, v      = sc.linalg.eig(matrix_C)
    print "eigsolve: ", timeit.default_timer() - eigsolve_stime
    print v
    for eigenvalue in la:
        print eigenvalue

def eigensolve(mat_A, mat_B):
    """Compute eigenvalues and eigenvectors: A.Z = E*B.Z
    Refactoring of eigensolve routine.

    """
    matrix_dim = len(mat_A)

    chol_stime = timeit.default_timer()
    matrix_UT = mpm.cholesky(mat_B) # lower triangular T
    print "matrix_UT (cholesky): ", timeit.default_timer() - chol_stime
    matU_stime = timeit.default_timer()
    matrix_U  = matrix_UT.transpose() # upper triangular
    print "matrix_U (transpose): ", timeit.default_timer() - matU_stime

    # Now get inverse of U
    # matrix_UI = mpm.inverse(matrix_U)
    matUI_stime = timeit.default_timer()
    matrix_UI = matrix_U**-1
    print "matrix_UI (inverse): ", timeit.default_timer() - matUI_stime

    # UI Transpose
    matUIT_stime = timeit.default_timer()
    matrix_UIT = matrix_UI.transpose()
    print "matrix_UIT (transpose): ", timeit.default_timer() - matUIT_stime

    # mat_C = UIT*A*UI
    matrixmult_stime = timeit.default_timer()
    matrix_C = matrix_UIT * mat_A
    matrix_C *= matrix_UI
    print "matrix_C (2 multiplications): ", timeit.default_timer() - matrixmult_stime
    print matrix_C



    eigsolve_stime = timeit.default_timer()
    eigval_C, eigvec_C = mpm.eig(matrix_C)
    print "eigsolve: ", timeit.default_timer() - eigsolve_stime

    eigsolve2_stime = timeit.default_timer()
    eigval_C, eigvec_C = LA.eig(np.array(matrix_C.tolist(),dtype=np.float64))
    print "eigsolve2: ", timeit.default_timer() - eigsolve2_stime

    # Normalization routine to make sure all leading column elements are
    # positive
    for i in xrange(0, matrix_dim):
        if eigvec_C[0, i] < 0:
            eigvec_C[:,i] = -eigvec_C[:,i]

    return (matrix_UI, eigval_C, eigvec_C)

def normalize_Z(matrix_UI, eigvec_C, eigval_C):
    matrix_dim = len(matrix_UI)
    matrix_Z = mpm.matrix(matrix_dim)
    for i in xrange(0, matrix_dim):
        vec_UIdotEVEC = matrix_UI * eigvec_C[:,i]
        for j in xrange(0, matrix_dim):
            matrix_Z[i,j] = vec_UIdotEVEC[j]

    matrix_Zn = mpm.matrix(matrix_dim)

    for i in xrange(0, matrix_dim):
        absmax = 0.0
        for j in matrix_Z[i,:]:
            if mpm.fabs(j) > mpm.fabs(absmax):
                absmax = j

        for j in xrange(0, matrix_dim):
            matrix_Zn[i, j] = mpm.fdiv(matrix_Z[i, j], absmax)

    low_E = mpm.mpf(eigval_C[0] )
    loc  = 0
    curr = 0
    for eigenvalue in eigval_C:
        if eigenvalue < low_E:
            loc = curr
            low_E = eigenvalue
        curr += 1

    print "Low energy state: ", low_E
    # print "Corresponding with Matrix Z[", loc, "]: ", matrix_Zn[loc,:]

    return (matrix_Zn, low_E, matrix_Zn[loc,:])
