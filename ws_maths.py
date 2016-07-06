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

def matrix_mpmath_2_numpy(matrix_mpm):
    """Transform mpmath matrix into numpy equivelent with longdouble as datatype
    """
    matrix_dim = len(matrix_mpm)
    matrix_np   = sc.zeros(dtype=sc.longdouble, shape=(matrix_dim,matrix_dim))

    for i in range(0, matrix_mpm.rows):
        for j in range(0, matrix_mpm.cols):
            matrix_np[i,j] = matrix_mpm[i,j]

    return matrix_np

def fixmul(A, B, prec=128):
    m = A.rows; p = B.rows; n = B.cols;
    A = [[A[i,j].to_fixed(prec) for j in range(p)] for i in range(m)]
    B = [[B[i,j].to_fixed(prec) for j in range(n)] for i in range(p)]
    C = [([0] * n) for r in range(m)]
    for i in range(m):
        for j in range(n):
            s = 0
            for k in range(p):
                s += A[i][k] * B[k][j]
            C[i][j] = s
    return mpm.matrix(C) * mpm.mpf(2)**(-2*prec)


def eigensolve_numpy(mat_A, mat_B):
    mat_A_np =  matrix_mpmath_2_numpy(mat_A)

    # matrix_LT  = sc.linalg.cholesky(mat_B_np, lower=True)
    # Need to continue using mpmath values for positive definiteness in
    # matrix for Cholesky step
    print "yiss?"
    matrix_LT = mpm.cholesky(mat_B)
    matrix_LT_np = matrix_mpmath_2_numpy(matrix_LT)
    print "UT Matrix \n", matrix_LT_np
    matrix_U   = matrix_LT_np.transpose()
    print "U Matrix \n", matrix_U
    matrix_UI  = sc.linalg.inv(matrix_U)
    print "UI Matrix \n", matrix_UI
    matrix_UIT = matrix_UI.transpose()
    print "UIT Matrix \n", matrix_UIT
    matrix_C = matrix_UIT.dot(mat_A_np)
    print "C Matrix \n", matrix_C
    matrix_C = matrix_C.dot(matrix_UI)
    print "C Matrix \n", matrix_C

    eigsolve_stime = timeit.default_timer()
    eigval, v      = sc.linalg.eig(matrix_C)
    print "eigsolve: ", timeit.default_timer() - eigsolve_stime

    low_E = eigval[0]
    for value in eigval:
        if value < low_E:
            low_E = value

    print '\n\n Lowest Eigenvalue'
    print low_E

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
    matUI_stime = timeit.default_timer()
    matrix_UI = matrix_U**-1
    print "matrix_UI (inverse): ", timeit.default_timer() - matUI_stime
    #print "UI Matrix \n", matrix_UI
    # return matrix_UI

    # UI Transpose
    matrix_UIT = matrix_UI.transpose()
    #print "UIT Matrix \n", matrix_UIT

    # mat_C = UIT*A*UI
    matrixmult_stime = timeit.default_timer()
    #matrix_C = matrix_UIT * mat_A
    matrix_C = fixmul(matrix_UIT, mat_A)
    #print "C Matrix \n", matrix_C
    #matrix_C *= matrix_UI
    matrix_C = fixmul(matrix_C, matrix_UI)
    print "matrix_C (2 multiplications): ", timeit.default_timer() - matrixmult_stime
    #print "C Matrix \n", matrix_C

    # eigval_C, eigvec_C = mpm.eig(matrix_C)

    # Convert to numpy matrix for solving
    print '===========================NP CONVERSION==========================='
    np_matrix = sc.zeros(dtype=sc.longdouble, shape=(matrix_dim,matrix_dim))

    copy_stime = timeit.default_timer()
    for i in range(0,matrix_C.rows):
        for j in range(0,matrix_C.cols):
            np_matrix[i,j] = matrix_C[i,j]
    print "copy to numpy matrix: ", timeit.default_timer() - copy_stime

    eigsolve_stime = timeit.default_timer()
    eigval, eigvec = sc.linalg.eig(np_matrix)
    print "eigsolve: ", timeit.default_timer() - eigsolve_stime
    # print eigval, eigvec

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
