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
# import sympy as sym
import scipy as sc
import sys
import timeit
from IPython.display import display

def eigensolve(mat_A, mat_B):
    """Compute eigenvalues and eigenvectors: A.Z = E*B.Z
    Refactoring of eigensolve routine.

    """
    matrix_dim = len(mat_A)

    matrix_UT = mpm.cholesky(mat_B) # lower triangular T
    # print "UT Matrix \n", matrix_UT
    matrix_U  = matrix_UT.transpose() # upper triangular
    #print "U Matrix \n", matrix_U
    # Now get inverse of U
    # matrix_UI = mpm.inverse(matrix_U)
    matrix_UI = matrix_U**-1
    #print "UI Matrix \n", matrix_UI
    # return matrix_UI

    # UI Transpose
    matrix_UIT = matrix_UI.transpose()
    #print "UIT Matrix \n", matrix_UIT

    # mat_C = UIT*A*UI
    matrix_C = matrix_UIT * mat_A
    #print "C Matrix \n", matrix_C
    matrix_C *= matrix_UI
    #print "C Matrix \n", matrix_C

    # eigval_C, eigvec_C = mpm.eig(matrix_C)

    # Convert to numpy matrix for solving
    print '===========================NP CONVERSION==========================='
    np_matrix = sc.zeros(dtype=sc.longdouble, shape=(matrix_dim,matrix_dim))

    for i in range(0,matrix_C.rows):
        for j in range(0,matrix_C.cols):
            np_matrix[i,j] = matrix_C[i,j]

    eigsolve_stime = timeit.default_timer()
    eigval, eigvec = sc.linalg.eig(np_matrix)
    print "eigsolve: ", timeit.default_timer() - eigsolve_stime
    print eigval, eigvec

    low_E = eigval[0]
    for value in eigval:
        if value < low_E:
            low_E = value

    print '\n\n Lowest Eigenvalue'
    print low_E





    # Normalization routine to make sure all leading column elements are
    # positive
    #for i in xrange(0, matrix_dim):
    #    if eigvec_C[0, i] < 0:
    #        eigvec_C[:,i] = -eigvec_C[:,i]

    #return (matrix_UI, eigval_C, eigvec_C)

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
