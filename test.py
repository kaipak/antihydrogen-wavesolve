"""
    Wavesolve main for testing functions
    
"""
# Standard libraries
import itertools
import timeit
import numpy as np
from sympy import *
from multiprocessing import Pool
from IPython.display import display

# Custom libraries
import ws_maths
import ws_physics

NSIZE = 20
Z = 1
PREC = 32

np.set_printoptions(precision=PREC)

@vectorize(['float32(float32, float32)'], target='gpu')
def Add(a, b):
    return a + b

def main():
    start_time = timeit.default_timer()
    init_printing()
    
    wave_equations = ws_physics.make_waves(2)
    psis = []
     
    # Pare down list to desire number of equations
    for i in xrange(0, NSIZE):
        psis.append(wave_equations[i])
        i += 1
    
    print psis
    
    
    print '\n'
    stop_time = timeit.default_timer()
    print "\n\nTime elapsed in seconds: "
    print stop_time - start_time

if __name__ == '__main__':
    main()