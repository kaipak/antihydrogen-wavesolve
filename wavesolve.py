"""

  Wavesolve program main program.  Runs main routine and calls calculations 
  from functions from math and physics libraries.
  
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
import itertools

import timeit
from sympy import *
#from bigfloat import *
from multiprocessing import Pool
from numpy import matrix, linalg, pi
from IPython.display import display

# Custom libraries
import ws_maths
import ws_physics

NSIZE = 100
Z = 1

# Bits of precison that will be used throughout application
#setcontext(quadruple_precision)

def main():
    start_time = timeit.default_timer()
    
    init_printing()
    
    # Generate wave equations.  Note comments on make_waves function as this creates (n+1)^3
    # Wave equations.
    wave_equations = ws_physics.make_waves(10)
    psis = []
    
    # Pare down list to desire number of equations
    for i in xrange(0, NSIZE):
        psis.append(wave_equations[i])
        i += 1
    """
    # Generate <Psi_i|Psi_j> over matrix
    # Chris Keating wave functions for test
    
    psis.append(ws_physics.gen_wavefunction(0, 0, 0))
    psis.append(ws_physics.gen_wavefunction(0, 0, 1))
    psis.append(ws_physics.gen_wavefunction(0, 1, 0))
    psis.append(ws_physics.gen_wavefunction(1, 0, 0))
    psis.append(ws_physics.gen_wavefunction(2, 0, 0))
    psis.append(ws_physics.gen_wavefunction(0, 0, 2))    
    """    
    
    psi_ij = [[0.0] * NSIZE for i in xrange(NSIZE)] 
    
    # Build matrix in following manner:
    # [ n11, n12, n13]
    # [      n22, n23]
    # [           n33]
    #
    # As this is a symmetric matrix (<psi_i|psi_j> <psi_j|psi_i>, we save roughly n/2
    # processing time by just mirroring from upper triangular
    for i in xrange(0, NSIZE):
        for j in xrange(i, NSIZE):
            psi_ij_val = ws_physics.get_qstate(psis[i], psis[j])
            psi_ij[i][j] = psi_ij_val
            # As this is a symmetric matrix, save some calculation time by populating
            # opposite side of matrix with same value.
            if i != j:
                psi_ij[j][i] = psi_ij_val
    
    # Now build <Psi_i|H|Psi_j> over matrix in similar manner as previous step
    hamiltonians = [] # Where we'll store wave equations with applied H operator
    for i in psis:
        hamiltonians.append(ws_physics.hamiltonian_r(i, Z))
    
    # Generate Summation of <Psi_i|H|Psi_j>
    psi_iHj = [[0.0] * NSIZE for i in xrange(NSIZE)]
        
    for i in xrange(0, NSIZE):
        for j in xrange(0, NSIZE):
            psi_iHj_val = ws_physics.get_qstate(psis[i], hamiltonians[j])
            psi_iHj[i][j] = psi_iHj_val
            # As this is a symmetric matrix, save some calculation time by populating
            # opposite side of matrix with same value.
                        
            if i != j:
                psi_iHj[j][i] = psi_iHj_val
            
    E = Symbol('E')
    matrix_ij = E * Matrix(psi_ij)
    matrix_iHj = Matrix(psi_iHj)
    print '\n'
    display(matrix_ij)
    print '\n'
    display(matrix_iHj)
     
    print '\n'
    submatrix = matrix_iHj - matrix_ij
    display(submatrix)
    print '\n'
    
    print '\n'
    #display(submatrix.det().normal().expand())
    #print solve(submatrix.det(), rational=False)
    
    stop_time = timeit.default_timer()
    print "\n\nTime elapsed in seconds: "
    print stop_time - start_time
    
if __name__ == '__main__':
    main()