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
from multiprocessing import Pool
from numpy import matrix, linalg, pi, set_printoptions
from IPython.display import display

# Custom libraries
import ws_maths
import ws_physics

NSIZE = 6
Z = 1
PREC = 32

set_printoptions(precision=PREC)

def main():
    start_time = timeit.default_timer()
    
    init_printing()
    
    # Generate wave equations.  Note comments on make_waves function as this creates (n+1)^3
    # Wave equations.
    # wave_equations = ws_physics.make_waves(2)
    psis = []
    
    
    # Chris Wave funcs
    psis.append(ws_physics.gen_wavefunction(0, 0, 0))
    psis.append(ws_physics.gen_wavefunction(0, 0, 1))
    psis.append(ws_physics.gen_wavefunction(0, 1, 0))
    psis.append(ws_physics.gen_wavefunction(1, 0, 0))
    psis.append(ws_physics.gen_wavefunction(2, 0, 0))
    psis.append(ws_physics.gen_wavefunction(0, 0, 2))  
    
    """
    # Pare down list to desire number of equations
    for i in xrange(0, NSIZE):
        psis.append(wave_equations[i])
        i += 1
    """
    
    matrix_ij = [[0.0] * NSIZE for i in xrange(NSIZE)] 
    
    # Build matrix in following manner:
    # [ n11, n12, n13]
    # [      n22, n23]
    # [           n33]
    #
    # As this is a symmetric matrix (<psi_i|psi_j> <psi_j|psi_i>, we save roughly n/2
    # processing time by just mirroring from upper triangular
    for i in xrange(0, NSIZE):
        print 'Constructing row', i, 'of <i|j>'
        for j in xrange(i, NSIZE):
            psi_ij = ws_physics.get_qstate(psis[i], psis[j])
            matrix_ij[i][j] = psi_ij
            # As this is a symmetric matrix, save some calculation time by populating
            # opposite side of matrix with same value.
            if i != j:
                matrix_ij[j][i] = psi_ij
    
    matrix_ij_time = timeit.default_timer()
    print "\n\nTime elapsed in seconds: "
    print matrix_ij_time - start_time
    
    # Now build <Psi_i|H|Psi_j> over matrix in similar manner as previous step
    hamiltonians = [] # Where we'll store wave equations with applied H operator
    for i in psis:
        hamiltonians.append(ws_physics.hamiltonian_r(i, Z))
    
    # Generate Summation of <Psi_i|H|Psi_j>
    matrix_iHj = [[0.0] * NSIZE for i in xrange(NSIZE)]
        
    for i in xrange(0, NSIZE):
        print 'Constructing row', i, 'of <i|H|j>'
        for j in xrange(0, NSIZE):
            psi_iHj = ws_physics.get_qstate(psis[i], hamiltonians[j])
            matrix_iHj[i][j] = psi_iHj
            # As this is a symmetric matrix, save some calculation time by populating
            # opposite side of matrix with same value.
            if i != j:
                matrix_iHj[j][i] = psi_iHj
            
    matrix_iHj_time = timeit.default_timer()
    print "\n\nTime elapsed in seconds: "
    print matrix_iHj_time - start_time
    
    print '\nEigensolving...'
    ws_maths.eigensolve(matrix_iHj, matrix_ij)
    
    print '\n'
    stop_time = timeit.default_timer()
    print "\n\nTime elapsed in seconds: "
    print stop_time - start_time
    
if __name__ == '__main__':
    main()