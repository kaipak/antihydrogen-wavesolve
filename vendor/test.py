"""
    Wavesolve main for testing functions
    
"""
# Standard libraries
import itertools
import timeit
import numpy as np
import multiprocessing as mp
from sympy import *
from multiprocessing import Pool
from IPython.display import display

# Custom libraries
import ws_maths
import ws_physics

NSIZE = 200
Z = 1
PREC = 32

np.set_printoptions(precision=PREC)

def main():
    start_time = timeit.default_timer()
    init_printing()
    
    wave_equations = ws_physics.make_waves(4)
    psis = []
     
    # Pare down list to desire number of equations
    for i in xrange(0, NSIZE):
        psis.append(wave_equations[i])
        i += 1
    pool = mp.Pool(processes=None)
    hamiltonians = [pool.apply_async(ws_physics.hamiltonian_r,args=(psis[x], Z)) for x in xrange(0, NSIZE)]
    print hamiltonians
    output = [p.get() for p in hamiltonians]
    print output
    
    
    print '\n'
    stop_time = timeit.default_timer()
    print "\n\nTime elapsed in seconds: "
    print stop_time - start_time

if __name__ == '__main__':
    main()