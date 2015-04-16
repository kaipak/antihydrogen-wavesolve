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
import ws_maths
from math import sqrt
from pprint import pprint
from scipy import linalg, array
from IPython.display import display
from sympy import *

PREC = 16

np.set_printoptions(precision=PREC)

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

def eigensolve(mat_A, mat_B):
    """
    Compute eigenvalues and eigenvectors of generalized problem A.Z = E*B.Z
    where A and B are symmetric nxn matrices and B must be positive-definite
    due to use of Cholesky decomposition.
    
    Routine is borrowed from IMSL
    """
    matrix_dim = len(mat_A)
    
    # B = UT.U
    # First step is to use Cholesky decomposition to break B into U-transpose and U
    # where U is the upper triangular component.
    matrix_U = cholesky_U(mat_B)
    matrix_UT = matrix_U.transpose()
    
    print matrix_U
    
    # Now get inverse of U
    matrix_UI = linalg.inv(matrix_U)
    
    # U inverse, transpose
    matrix_UIT = matrix_UI.transpose()
    
    # Create a matrix C = UIT.A.UI
    matrix_C = matrix_UIT.dot(mat_A).dot(matrix_UI)
    print '\nMATRIX C'
    display(matrix_C)

    # Note: this doesn't give an array of Eigenvectors as may be expected.  Instead
    # gives matrix where eigenvectors are columns. Inverting so easier to do dot product
    # below.
    eigval_C, eigvec_C = linalg.eigh(matrix_C)
    eigvec_C = eigvec_C.transpose()
    
    print 'EIGENVALUES'
    display(eigval_C)
    print '\nEIGENVECTORS'
    display(eigvec_C)
    
    matrix_Z =  [[0.0] * matrix_dim for i in xrange(matrix_dim)]
    print '\nUI.EigVecC'
    for i in xrange(0, matrix_dim):
        vec_UIdotEVEC = matrix_UI.dot(eigvec_C[i])
        for j in xrange(0, matrix_dim):
            matrix_Z[i][j] = vec_UIdotEVEC[j]
    print matrix_Z
    
    for i in xrange(0, matrix_dim):
        absmax = 0.0
        for j in matrix_Z[i]:
            if np.abs(j) > np.abs(absmax):
                absmax = j
        
        for j in xrange(0, matrix_dim):
            matrix_Z[i][j] = np.divide(matrix_Z[i][j], absmax)
         
    print '\n'
    print 'Coefficients: '
    print matrix_Z[0]
    
    
def rand_matrix(msize):

    """
    Generate and return a random matrix that's *almost* always positive definite
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