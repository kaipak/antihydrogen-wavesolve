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
import sys
import mpmath as mpm
from numpy import array, matrix, linalg, pi, set_printoptions
from IPython.display import display
from scipy import linalg


# Custom libraries
import ws_maths
import ws_physics

# Physical parameters
ALPHA = .70120
BETA = .70120
GAMMA = 0
ZED = 1

# Application attributes
NSIZE = 4
PREC = 128

set_printoptions(precision=PREC)

def main():
    start_time = timeit.default_timer()
    mpm.mp.prec = PREC
    
    
    ws_physics.static_params(0, 0, 0, 0, 0, 0, 0, 1, NSIZE, 16)
    psis = []
    
    
    # Chris Wave funcs. Base case, don't change!!
    psis.append(ws_physics.gen_wavefunction(0, 0, 0, ALPHA, BETA, GAMMA))
    psis.append(ws_physics.gen_wavefunction(0, 0, 1, ALPHA, BETA, GAMMA))
    psis.append(ws_physics.gen_wavefunction(0, 2, 0, ALPHA, BETA, GAMMA))
    psis.append(ws_physics.gen_wavefunction(1, 0, 0, ALPHA, BETA, GAMMA))
    psis.append(ws_physics.gen_wavefunction(2, 0, 0, ALPHA, BETA, GAMMA))
    psis.append(ws_physics.gen_wavefunction(0, 0, 2, ALPHA, BETA, GAMMA))

    """    
    # Pare down list to desire number of equations
    for i in xrange(0, NSIZE):
        psis.append(wave_equations[i])
        i += 1
    """

    matrix_ij1 =  build_matrix1(psis, psis, '<i|j>')
    matrix_ij2 =  build_matrix2(psis, psis, '<i|j>')
    
    matrix_ij_time = timeit.default_timer()
    print "\n\nTime elapsed in seconds: "
    print matrix_ij_time - start_time
    
    # Now build <Psi_i|H|Psi_j> 
    hamiltonians = [] # wave equations with applied H operator
    for i in psis:
        hamiltonians.append(ws_physics.hamiltonian_r(i, ALPHA,
                                                      BETA, GAMMA))
    
    matrix_iHj1 = build_matrix1(psis, hamiltonians, '<i|H|j>')
    matrix_iHj2 = build_matrix2(psis, hamiltonians, '<i|H|j>')
    print matrix_ij1
    print '\n'
    print matrix_iHj1
        
    matrix_iHj_time = timeit.default_timer()
    print "\n\nTime elapsed in seconds: "
    print matrix_iHj_time - start_time
    
    print '\nEigensolving...'
    ws_maths.eigensolve1(matrix_iHj1, matrix_ij1)
    ws_maths.eigensolve2(matrix_iHj2, matrix_ij2)

    print '\n'
    #print matZ1
    print '\nFrom eigensolve2'
    #print matZ2
    
    print '\n'
    stop_time = timeit.default_timer()
    print "\n\nTime elapsed in seconds: "
    print stop_time - start_time

def build_matrix1(psis_i, psis_j, bracket_notation):
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
    #matrix = mpm.matrix(NSIZE)
    
    
    for i in xrange(0, NSIZE):
        print 'Constructing row', i, 'of', bracket_notation
        
        pool = mp.Pool(processes = None)
        psirowobj = [pool.apply_async(ws_physics.get_qstate,
                                      args=(psis_i[i],psis_j[j], 
                                            ALPHA, BETA, GAMMA))
                     for j in xrange(i, NSIZE)]
        psirow = [p.get() for p in psirowobj]
        pool.terminate()
        
        # Now populate row and its transpose column since this is a symmetric matrix
        for j in xrange(0, len(psirow)):
            matrix[i][i+j] = psirow[j]
            matrix[i+j][i] = psirow[j]
        
    return matrix

def build_matrix2(psis_i, psis_j, bracket_notation):
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
    #matrix = [[0.0] * NSIZE for i in xrange(NSIZE)]
    matrix = mpm.matrix(NSIZE)
    
    
    for i in xrange(0, NSIZE):
        print 'Constructing row', i, 'of', bracket_notation
        
        pool = mp.Pool(processes = None)
        psirowobj = [pool.apply_async(ws_physics.get_qstate,
                                      args=(psis_i[i],psis_j[j], 
                                            ALPHA, BETA, GAMMA))
                     for j in xrange(i, NSIZE)]
        psirow = [p.get() for p in psirowobj]
        pool.terminate()
        
        # Now populate row and its transpose column since this is a symmetric matrix
        for j in xrange(0, len(psirow)):
            matrix[i,i+j] = psirow[j]
            matrix[i+j,i] = psirow[j]
        
    return matrix
    
if __name__ == '__main__':
    main()
