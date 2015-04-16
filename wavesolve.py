"""

  Wavesolve program main program.  Runs main routine and calls calculations 
  from functions from math and physics libraries.
  
  @author Kai Pak
  @start_date January 1, 2014
  @current_version 1.0
  
  History
  
"""

# Standard libraries
import itertools
import timeit
import multiprocessing as mp
from numpy import array, matrix, linalg, pi, set_printoptions
from IPython.display import display

# Custom libraries
import ws_maths
import ws_physics

NSIZE = 4
Z = 1
PREC = 16

set_printoptions(precision=PREC)

def main():
    start_time = timeit.default_timer()
    
    # Generate wave equations.  Note comments on make_waves function as this creates (n+1)^3
    # Wave equations.
    wave_equations = ws_physics.make_waves(4)
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
    matrix_ij = matrix_ij
    print matrix_ij
    
    """
    Build matrix in following manner:
    [ n11, n12, n13]
    [      n22, n23]
    [           n33]
   
    As this is a symmetric matrix (<psi_i|psi_j> <psi_j|psi_i>, we save roughly n/2
    processing time by just mirroring from upper triangular
    
    Multiprocessing version
    
    """
    for i in xrange(0, NSIZE):
        print 'Constructing row', i, 'of <i|j>'
        pool = mp.Pool(processes=None)
        psirowobj = [pool.apply_async(ws_physics.get_qstate, args=(psis[i],psis[j])) \
                  for j in xrange(0, NSIZE)]
        psirow = [p.get() for p in psirowobj]
        matrix_ij[i] = psirow
        
        """
        # Now populate row and it's transpose column since this is a symmetric matrix
        for j in xrange(0, len(psirow)):
            matrix_ij[i][i+j] = psirow[j]
            matrix_ij[i+j][i] = psirow[j]
        """    
        
    for row in matrix_ij:
        print row
        
    matrix_ij_time = timeit.default_timer()
    print "\n\nTime elapsed in seconds: "
    print matrix_ij_time - start_time
    
    # Now build <Psi_i|H|Psi_j> over matrix in similar manner as previous step
    # Multiprocessor version
    hamiltonians = [] # Where we'll store wave equations with applied H operator
    for i in psis:
        hamiltonians.append(ws_physics.hamiltonian_r(i, Z))
    
    # Generate Summation of <Psi_i|H|Psi_j>
    matrix_iHj = [[0.0] * NSIZE for i in xrange(NSIZE)]
        
    for i in xrange(0, NSIZE):
        print 'Constructing row', i, 'of <i|H|j>'
        pool = mp.Pool(processes=None)
        psi_iHj_rowobj = [pool.apply_async(ws_physics.get_qstate, args=(psis[i],\
                          hamiltonians[j])) for j in xrange(i, NSIZE)]
        psi_iHj_row = [p.get() for p in psi_iHj_rowobj]
        matrix_iHj[i] = psi_iHj_row
            
    print matrix_iHj
    
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