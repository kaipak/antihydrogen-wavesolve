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
  
"""

import cProfile
import numpy as np
import pstates_maths
from bigfloat import *
from enthought.traits.api import *
from enthought.traits.ui.api \
    import View, ArrayEditor, Item, Group, VSplit
from enthought.traits.ui.menu import NoButtons
from scipy import linalg, random
from sympy import *
import timeit

# Bits of precison that will be used throughout application
setcontext(quadruple_precision)

class InputBox(HasTraits):
    """
       Class to handle input of quantum variables.
    
    """
    angMo = Float(0, label="Angular Momentum (L)")
    zCharge = Int(0, label="Nuclear Charge (Z)")
    calc = Button("Calculate")
    view = View(Item('angMo'), 
                Item('zCharge'), 
                Item('calc', show_label=False))
         
class Pstates(HasTraits):
    """
        Class to handle some output stuff.  So you can see what application is
        doing.  Most calculations will be methods in this class.
    """
    string = String()
    view= View( Item('string',show_label=False, springy=True, style='readonly' ))

    def test_calc(Float):
        return sqrt(Float)
    

class MainWindow(HasTraits):
    """
        Main interface of app that takes input from user.  Output is rather 
        large so it'll go to an output file.
        
    """ 
    inbox = Instance(InputBox, ())
    pstates = Instance(Pstates, ())
    
    view = View(VSplit(Item('inbox', style='custom',
                             springy=True,
                             label="Input"),
                       Item('pstates', style='custom',
                             springy=True,
                             label="pstates")),
                resizable=True)
    
    
    

    def _get_some_calc(self): 
        return const_pi()

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
    
def gen_gamma(l, m, n):
    """
    Using Harris, Frolov, and Smith paper (J. Chem Phys., Vol. 121, No. 13
    Oct. 2004), generate radial integral or gamma function 
    (alpha=a, beta=b, gamma=g) using l,m,n quantum states.
    """
    a = Symbol('a')
    b = Symbol('b')
    g = Symbol('g')
    
    f = (a + b)(a + g)

def main():
    #mw = MainWindow()
    #mw.configure_traits()
    
    start_time = timeit.default_timer()
    gen_gamma(1,2,3)
    
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