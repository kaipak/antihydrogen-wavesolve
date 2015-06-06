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

import numpy as np
import mpmath as mpm
import sympy as sym
import sys
from math import sqrt
from pprint import pprint
from scipy import linalg, array
from IPython.display import display

def comb(n,r):
    """
    Return combination of r choices from n elements.
    nCr = n!/r!(n-r)!
    
    """
    return (np.math.factorial(n))/(np.math.factorial(r)*np.math.factorial(n-r))

#setcontext(quadruple_precision)

""" 
    Cholesky from linalg may not work as it doesn't have precision required.

"""

def cholesky_L(matrix):
    return linalg.cholesky(matrix, lower=True)
       
def cholesky_U(matrix):
    return linalg.cholesky(matrix, lower=False)

def b_cholesky_L(matrix):
    """
    Performs a Cholesky decomposition of A, which must 
    be a symmetric and positive definite matrix. The function
    returns the lower variant triangular matrix, L.
    
    Modified source code from www.quantstart.com.
    
    06/06/14 - Initial results not great. Complexity appears to be polynomial:
               O(f(n^2...).
    
    """
    
    n = len(matrix)

    # Create zero matrix for L
    L = [[0.0] * n for i in xrange(n)]

    # Perform the Cholesky decomposition
    for i in xrange(n):
        for k in xrange(i+1):
            tmp_sum = sum(L[i][j] * L[k][j] for j in xrange(k))
            
            if (i == k): # Diagonal elements
                # LaTeX: l_{kk} = \sqrt{ a_{kk} - \sum^{k-1}_{j=1} l^2_{kj}}
                L[i][k] = sqrt(matrix[i][i] - tmp_sum)
            else:
                """
                LaTeX: l_{ik} = \frac{1}{l_{kk}} \left( a_{ik} - 
                \sum^{k-1}_{j=1} l_{ij} l_{kj} \right)
                """
                L[i][k] = (1.0 / L[k][k] * (matrix[i][k] - tmp_sum))
    return L

def eigensolve1(mat_A, mat_B):
    """
    Compute eigenvalues and eigenvectors of generalized problem A.Z = E*B.Z
    where A and B are symmetric nxn matrices and B must be positive-definite
    due to use of Cholesky decomposition.
    
    Routine is borrowed from IMSL
    """
    print '\nEigensolve1: '
    matrix_dim = len(mat_A)
    
    # B = UT.U
    # First step is to use Cholesky decomposition to break B into 
    # U-transpose and U where U is the upper triangular component.
    matrix_U = cholesky_U(mat_B)
    matrix_UT = matrix_U.transpose()
    
    
    # Now get inverse of U
    matrix_UI = linalg.inv(matrix_U)

    # U inverse, transpose
    matrix_UIT = matrix_UI.transpose()
    
    # Create a matrix C = UIT.A.UI
    matrix_C = matrix_UIT.dot(mat_A).dot(matrix_UI)

    eigval_C, eigvec_C = linalg.eigh(matrix_C)

    matrix_Z =  [[0.0] * matrix_dim for i in xrange(matrix_dim)]
    for i in xrange(0, matrix_dim):
        vec_UIdotEVEC = matrix_UI.dot(eigvec_C[:,i])
        for j in xrange(0, matrix_dim):
            matrix_Z[i][j] = vec_UIdotEVEC[j]
    
    for i in xrange(0, matrix_dim):
        absmax = 0.0
        for j in matrix_Z[i]:
            if np.abs(j) > np.abs(absmax):
                absmax = j
        
        for j in xrange(0, matrix_dim):
            matrix_Z[i][j] = np.divide(matrix_Z[i][j], absmax)

    print matrix_Z
         
def eigensolve(mat_A, mat_B):
    """Compute eigenvalues and eigenvectors: A.Z = E*B.Z
    Refactoring of eigensolve routine.

    """
    matrix_dim = len(mat_A)

    matrix_UT = mpm.cholesky(mat_B) # lower triangular T
    matrix_U  = matrix_UT.transpose() # upper triangular

    # Now get inverse of U
    # matrix_UI = mpm.inverse(matrix_U)
    matrix_UI = matrix_U**-1
    # return matrix_UI

    # UI Transpose
    matrix_UIT = matrix_UI.transpose()

    # mat_C = UIT*A*UI 
    matrix_C = matrix_UIT * mat_A 
    matrix_C *= matrix_UI

    eigval_C, eigvec_C = mpm.eig(matrix_C)
    # eigval_C, eigvec_C = mpm.eig_sort(eigval_C, eigvec_C)
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
            matrix_Zn[i, j] = np.divide(matrix_Z[i, j], absmax)

    low_E = eigval_C[0] 
    loc  = 0
    curr = 0
    for eigenvalue in eigval_C:
        if eigenvalue < low_E:
            loc = curr
            low_E = eigenvalue
        curr += 1

    print "Low energy state: ", low_E
    print "Corresponding with Matrix Z[", loc, "]: ", matrix_Zn[loc,:]

    return (matrix_Zn, low_E, matrix_Zn[loc,:])


def rand_matrix(msize):

    """Generate and return a random matrix

    Matrix that's *almost* always positive definite
    and Hermitian.  Primarily for testing purposes
    
    """    

    A = np.random.rand(msize,msize)
    
    # Generate a Hermitian Matrix that's semi-definite
    matrix = A + A.transpose()
    
    # Ensure positive definiteness.  Hackish, but just generates an nxn identity
    # matrix to add to the randomly generated Hermitian.
    eigenadd = 100 * np.eye(msize)
    
    matrix = matrix + eigenadd
    
    return matrix