"""

  Wavesolve program main program.  Runs main routine and calls calculations 
  from functions from math and physics libraries.
  
  
  Built using Canopy Enthought IDE and frameworks.
  
  @required_packages
  bigfloat - Python wrapper for MPFR for handling arbitrarily large precision 
             floating point numbers.
  MPFR - GNU C library for multiple-precision floating-point computation.  Note,
         GMP is required for installation.
         http://www.mpfr.org
  GMP - GNU Multiple Precision library for arbitrary precision 
        arithmetic, operating on signed integers, rational numbers, and 
        floating-point numbers.
        https://gmplib.org/
  
  @author Kai Pak
  @start_date January 1, 2014
  @current_version 1.0
  
  History
  
"""

# Standard libraries
import sympy
import itertools
import timeit
#from bigfloat import *
from multiprocessing import Pool
from numpy import matrix, linalg
from IPython.display import display

# Custom libraries
import ws_maths
import ws_physics

# Bits of precison that will be used throughout application
#setcontext(quadruple_precision)

def main():
    start_time = timeit.default_timer()
    
    # TEST CHOLESKY
    # bigmatrix = ws_maths.rand_matrix(10)
    # for row in bigmatrix:
    #    print row
    #print linalg.eig(bigmatrix)
    #L_alt = ws_maths.b_cholesky_L(bigmatrix)
    #L = ws_maths.cholesky_L(bigmatrix)
    #print "\n\nSome Cholesky Decomposition of a matrix:\n"
    #for row in L_alt:
    #    print row
    
    # Generate wave equations.  Note comments on make_waves function as this creates (n+1)^3
    # Wave equations.
    wave_equations = ws_physics.make_waves(10)
    testbed = []
    
    # Pare down list to desire number of equations
    # This is a n*n matrix, so variable to determine the ultimate size.
    n_size = 2
    
    for i in xrange(0, n_size):
        testbed.append(wave_equations[i])
        i += 1
        
    # Generate <Psi_i|Psi_j> over matrix
    psi_ij = [[0.0] * n_size for i in xrange(n_size)] 
    
    for i in xrange(0, n_size):
        for j in xrange(0, n_size):
            # Get coefficient, L, M, N values.
            print '<psi', i, '|', 'psi', j, '>'
            clmns = ws_physics.extract_clmn(testbed[i], testbed[j])
            print '\n'
            
            testvalue = 0
                
            for k in clmns:
                c,l,m,n = k
                tempval = c * ws_physics.hfs_gamma(l, m, n)
                testvalue += tempval
            psi_ij[i][j] = testvalue
            j += 1
        i += 1
    
    for row in psi_ij:
        print row
    
    neo = matrix(psi_ij)
    
    print '\n\nCholesky decomp below:\n'
    print linalg.cholesky(neo)
    
    stop_time = timeit.default_timer()
    print "\n\nTime elapsed in seconds: "
    print stop_time - start_time

if __name__ == '__main__':
    main()