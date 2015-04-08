"""

  Pseudostates program main program.  Runs most basic physical calculations 
  and calls mathematical functions from libraries.
  
  Using RMPS, calculates set of orthonormal hydrogenic pseudostates for a given
  angular momentum.
  
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
	meh
  
"""

import cProfile
import numpy as np
import pstates_maths
from bigfloat import *
from enthought.traits.api import *
from enthought.traits.ui.api \
    import View, ArrayEditor, Item, Group, VSplit
from enthought.traits.ui.menu import NoButtons
from scipy import linalg, random, misc
from sympy import *
import timeit

# Bits of precison that will be used throughout application
setcontext(quadruple_precision)



def rand_matrix(msize):
    
    #tarray = []
    """
    for column in range(0, msize):
        # Ranges for random ints
        lower = msize
        upper = msize + 10
        arow = []
        for row in range(0, msize):
            arow.append(random.randint(lower,upper))
            lower -= 1
            upper -= 1
        tarray.append(arow)
    """
    
    A = random.rand(msize,msize)
    
    # Generate a Hermitian Matrix that's semi-definite
    tmatrix = A + A.transpose()
    
    # Ensure positive definiteness.  Hackish, but just generates an nxn identity
    # matrix to add to the randomly generated Hermitian.
    eigenadd = 100 * np.eye(msize)
    
    tmatrix = tmatrix + eigenadd
    
    return tmatrix
    
def comb(n,r):
    """
    Return combination of r choices from n elements
    
    """
    return (np.math.factorial(n))/(np.math.factorial(r)*np.math.factorial(n-r))
    
def gen_gamma(l, m, n):
    """
    Using Harris, Frolov, and Smith paper (J. Chem Phys., Vol. 121, No. 13
    Oct. 2004), generate radial integral or gamma function 
    (alpha=a, beta=b, gamma=g) using l,m,n quantum states.
    
    """
    fact_coef = 2 * np.math.factorial(l) * np.math.factorial(m) \
                * np.math.factorial(n)
    
    
    a = Symbol('a')
    b = Symbol('b')
    g = Symbol('g')
    
    x = a + b
    y = a + g
    z = b + g
    
    # Prime our function with the f, g, h functions and without the coefficient
    # so that doesn't also get summed
    
    l_prime = 0
    m_prime = 0
    n_prime = 0
    
    
    big_gamma = 1/(x**2*y**2*z**2)
    
    #big_gamma = 1/(x**(m)*y**(l)*z**(n))
    #print big_gamma
   
    test = 0
    
    for l_prime in range(0,l+1):
        for m_prime in range(0,m+1):     
            for n_prime in range(1,n+1):
                big_gamma = big_gamma + (comb((m - m_prime + l_prime), l_prime) * \
                 comb((l - l_prime + n_prime), n_prime) * \
                 comb((n - n_prime + m_prime), m_prime) / \
                (x**(m-m_prime+l_prime+1)*y**(l-l_prime+n_prime+1)*z**(n-n_prime+m_prime+1)))
    
    big_gamma = fact_coef*big_gamma
    
    #return str(big_gamma.simplify()).replace("**","^")
    return big_gamma

def main():
    #mw = MainWindow()
    #mw.configure_traits()
    
    start_time = timeit.default_timer()
    print gen_gamma(1,1,1)
    
    """
    bigmatrix = rand_matrix(32)
    for row in bigmatrix:
        print row

    # print linalg.eig(bigmatrix)
    L_alt = pstates_maths.b_cholesky_L(bigmatrix)
    #L = pstates_maths.cholesky_L(bigmatrix)
    
    #print "\n\nDisplaying some maths.\n"
    #print x
    
    
    print "\n\nSome Cholesky Decomposition of a matrix:\n"
    
    for row in L_alt:
        print row
    """
    stop_time = timeit.default_timer()
    
    
    print "\n\nTime elapsed in seconds: "
    print stop_time - start_time


if __name__ == '__main__':
    main()
