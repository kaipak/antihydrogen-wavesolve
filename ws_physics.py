from IPython.display import display
import numpy as np
import scipy.misc as sc
from sympy import *
import itertools
import ws_maths



A = Symbol('A')
a = Symbol('alpha')
b = Symbol('beta')
g = Symbol('gamma')
m = Symbol('m')
n = Symbol('n')
r1 = Symbol('r1')
r2 = Symbol('r2')
r12 = Symbol('r12')

# Hylleras coordinate stuff (transform later?)
s = r1 + r2
t = r1 - r2
u = r12

Z = Symbol('Z')

# Constants
LMN_LENGTH = 3
    
def hfs_gamma(l, m, n):
    """
    Using Harris, Frolov, and Smith paper (J. Chem Phys., Vol. 121, No. 13
    Oct. 2004), generate radial integral or gamma function 
    (alpha=a, beta=b, gamma=g) using l,m,n quantum states.
    
    """
    fact_coef = 2 * np.math.factorial(l) * np.math.factorial(m) \
                * np.math.factorial(n)
        
    x = a + b
    y = a + g
    z = b + g
        
    big_gamma = Mul(0)
    
    for l_prime in range(0,l+1):
        for m_prime in range(0,m+1):     
            for n_prime in range(0,n+1):
                big_gamma = big_gamma + (sc.comb((m - m_prime + l_prime), l_prime) * \
                 sc.comb((l - l_prime + n_prime), n_prime) * \
                 sc.comb((n - n_prime + m_prime), m_prime) / \
                (x**(m-m_prime+l_prime+1)*y**(l-l_prime+n_prime+1)*z**(n-n_prime+m_prime+1)))
    
    big_gamma = fact_coef*big_gamma
    
    print str(big_gamma.simplify()).replace("**","^")
    print '\n'
    print str(big_gamma).replace("**","^")
    
    return big_gamma
    
def gen_wavefunction(l, m, n):
    """
    Generate wave function with arbitrary number of terms.  Need to later convert
    to Hylleraas coordinate system.  
    
    """
    exponential = -((a * r1) + (b * r2) + (g * r12)) 
    psi = exp(exponential)
    # Generate coefficients
    coef = Symbol('c' + str(l) + str(m) + str(n))
    wave_equation = coef * (s)**l * (t)**(2*m) * (u)**n * psi
    return wave_equation

def make_waves(iterables):
    """
    Make a bunch of wave functions over Cartesian product from 0 to iterables over length of LMN (3).  
    Somewhat of a binary counting system that increments l, m, n.
    
    Note on iterables variable.  To make it more clear to user, up to what cardinal number he/she
    is using. 1, up to '1', 2, up to '2' and so on.  Otherwise, it 2 would be 0 and 1, and that gets
    confusing.
    
    Produces (iterables + 1)^3 permutations of wave functions
    
    """
    iterables += 2
    my_list = []
    lower_bound = 0
    for i in xrange(0, iterables):
        for j in list(itertools.product(range(0,i), repeat = LMN_LENGTH)):
            if j not in my_list:
                my_list.append((j))
        
    for tup in my_list:
        l,m,n = tup
        display(gen_wavefunction(l,m,n))
    
    print len(my_list)    
    
def hamiltonian_r(wfunc):
    """
    Solve Hamiltonian of wave function in r1, r2, r12 coordinate system
    
    """
    # hamiltonian = Mul(0)
    hamiltonian = (-diff(wfunc, r1,2)/2 - diff(wfunc, r2, 2)/2 - diff(wfunc, r12, 2) - \
                  ((1/r1) * diff(wfunc, r1, 1)) - \
                  ((1/r2) * diff(wfunc, r2, 1)) - \
                  ((2/r12) * diff(wfunc, r12, 1)) - \
                  (((r1**2 - r2**2 + r12**2)/(r1 * r12)) * diff(diff(wfunc, r12,1),r1,1)) - \
                  (((r2**2 - r1**2 + r12**2)/(r2 * r12)) * diff(diff(wfunc, r12,1),r2,1)) - \
                  ((Z/r1) + (Z/r2) - (1/r12)) * wfunc) * \
                  (r1 * r2 * r12)
            
    return hamiltonian
    
    
    
def quantum_state(wf1, wf_H2):
    """
    Solve <psy1|H|psy2> in r1, r2, r12 coordinate system.  Second wf should have
    Hamiltonian already applied to it. 
    
    """

    # Combine <psy1|H|psy2> with Volume element for integration: 
    # 8pi^2 psy1(conj)*H*psy2(r1)(r2)(r12)
    equation = wf1 * wf_H2 * r1 * r2 * r12
    display(equation)
    # equation = 8 * pi**2 * wf1 * r1 * r2 * r12
    
    piece1 = integrate(equation, (r12, (r1-r2), (r1+r2)), conds='none')
    piece1 = integrate(piece1, (r1, 0, r2), conds='none' )
    display(piece1)
    
    piece2 = integrate(equation, (r12, (r2-r1), (r1+r2)), conds='none').simplify()
    piece2 = integrate(piece2, (r1, r2, oo))
    
    #display(integrate(exp(-r1), (r1, r2, oo))) 
    
    final = 8*pi**2 * (piece1 + piece2)
    #display(final)
    
    display(integrate(final, (r2, 0, oo), conds='none'))
    

    """
    # Separate absolute value parts into piece-wise integrals
    q_state1 = integrate(equation, (r1, r1-r2, r1+r2), conds='none')
    display(q_state1)
    q_state1 = integrate(q_state1, (r12, 0, r2), conds='none')
    display(q_state1)
    
    print "q_state2"
    q_state2 = integrate(integrate(equation, (r1, -(r1-r2), r1+r2), conds='none'),(r12, r2, oo), conds='none')
    display(q_state2)
    
    q_state = q_state1 + q_state2
    display(q_state)
    
    #q_statefinal = (8*pi**2) * integrate(q_state, (r2, 0, oo), conds='none')
    #display(q_statefinal)
    
    
    #return q_statefinal
    """
    
    
    
    
    