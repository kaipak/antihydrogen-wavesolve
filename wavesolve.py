"""
  Wavesolve program main program.  Runs main routine and calls calculations 
  from functions from math and physics libraries.
  
  --author - Kai Pak
  --start_date - January 1, 2014
  --relase date - April 30, 2015
  --current_version - 1.0
  
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

# Physical parameters
ALPHA = .70120
BETA = .70120
GAMMA = 0
ZED = 1

# Application attributes
NSIZE = 20
PREC = 32

set_printoptions(precision=PREC)

def main():
    start_time = timeit.default_timer()
    
    
    ws_physics.set_params(ALPHA, BETA, GAMMA, ZED, PREC)
    # Generate wave equations.  Note comments on make_waves function as this
    # creates (n+1)^3 wave equations.
    wave_equations = ws_physics.make_waves(4)
    psis = []
    
    
    # Chris/Jack Wave funcs
    psis.append(ws_physics.gen_wavefunction(0, 0, 0))
    psis.append(ws_physics.gen_wavefunction(0, 0, 1))
    psis.append(ws_physics.gen_wavefunction(0, 2, 0))
    psis.append(ws_physics.gen_wavefunction(1, 0, 0))
    psis.append(ws_physics.gen_wavefunction(2, 0, 0))
    psis.append(ws_physics.gen_wavefunction(0, 0, 2))  
    psis.append(ws_physics.gen_wavefunction(1, 0, 1))
    psis.append(ws_physics.gen_wavefunction(0, 2, 1))
    psis.append(ws_physics.gen_wavefunction(0, 0, 3))
    psis.append(ws_physics.gen_wavefunction(0, 2, 2))
    psis.append(ws_physics.gen_wavefunction(1, 2, 0))
    psis.append(ws_physics.gen_wavefunction(3, 0, 0))  
    psis.append(ws_physics.gen_wavefunction(0, 2, 4))  
    psis.append(ws_physics.gen_wavefunction(0, 0, 4))
    psis.append(ws_physics.gen_wavefunction(0, 0, 5))  
    psis.append(ws_physics.gen_wavefunction(0, 2, 3))  
    psis.append(ws_physics.gen_wavefunction(2, 2, 0))
    psis.append(ws_physics.gen_wavefunction(4, 0, 0))  
    psis.append(ws_physics.gen_wavefunction(1, 2, 1))  
    psis.append(ws_physics.gen_wavefunction(0, 4, 0))  
    
    """    
    # Pare down list to desire number of equations
    for i in xrange(0, NSIZE):
        psis.append(wave_equations[i])
        i += 1
    """

    matrix_ij =  build_matrix(psis, psis, '<i|j>')
    
    matrix_ij_time = timeit.default_timer()
    print "\n\nTime elapsed in seconds: "
    print matrix_ij_time - start_time
    
    # Now build <Psi_i|H|Psi_j> 
    hamiltonians = [] # wave equations with applied H operator
    for i in psis:
        hamiltonians.append(ws_physics.hamiltonian_r(i))
    
    matrix_iHj = build_matrix(psis, hamiltonians, '<i|H|j>')
        
    matrix_iHj_time = timeit.default_timer()
    print "\n\nTime elapsed in seconds: "
    print matrix_iHj_time - start_time
    
    print '\nEigensolving...'
    ws_maths.eigensolve(matrix_iHj, matrix_ij)
    
    print '\n'
    stop_time = timeit.default_timer()
    print "\n\nTime elapsed in seconds: "
    print stop_time - start_time

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
                                      args=(psis_i[i],psis_j[j]))
                     for j in xrange(i, NSIZE)]
        psirow = [p.get() for p in psirowobj]
        pool.terminate()
        
        # Now populate row and its transpose column since this is a symmetric matrix
        for j in xrange(0, len(psirow)):
            matrix[i][i+j] = psirow[j]
            matrix[i+j][i] = psirow[j]
        
    return matrix
    
if __name__ == '__main__':
    main()
