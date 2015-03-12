from IPython.display import display
import numpy as np
import scipy.misc as sc
from sympy import *
from re import *
import itertools
import ws_maths

A = Symbol('A')
a = .70120
b = .70120
g = 0
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
EXPONENTIAL = -((a * r1) + (b * r2) + (g * r12))
# Base wave equation we'll build on
PSI = exp(EXPONENTIAL)
LMN_LENGTH = 3
    
def hfs_gamma(l, m, n):
    """
    Using Harris, Frolov, and Smith paper (J. Chem Phys., Vol. 121, No. 13
    Oct. 2004), generate radial integral or gamma function 
    (alpha=a, beta=b, gamma=g) using l,m,n quantum states.
    
    """
    fact_coef = 2 * np.math.factorial(l) * np.math.factorial(m) \
                * np.math.factorial(n)
        
    x = 2*a + 2*b
    y = 2*a + 2*g
    z = 2*b + 2*g
        
    big_gamma = Mul(0)
    
    for l_prime in range(0,l+1):
        for m_prime in range(0,m+1):     
            for n_prime in range(0,n+1):
                big_gamma = big_gamma + (sc.comb((m - m_prime + l_prime), l_prime) * \
                 sc.comb((l - l_prime + n_prime), n_prime) * \
                 sc.comb((n - n_prime + m_prime), m_prime) / \
                (x**(m-m_prime+l_prime+1)*y**(l-l_prime+n_prime+1)*z**(n-n_prime+m_prime+1)))
    
    big_gamma = fact_coef*big_gamma
    
    return big_gamma
    
def test_mp(wf1):
    print wf1 * wf1

def extract_clmn(func1, func2):
    """
    Generalized function to extract coefficient, r1, r2, and r12 (correspending to L, M, and N)
    values from <Psi_i|Psi_j> or <Psi_i|Operator|Psi_j> for primary purpose of being applied
    to Harris, Frolov, and Smith Equation.
    
    <Psi_i|Psi_j> = integral(Psy_i*Psy_j*r1*r2*r12dr1dr2dr12).  Pull out powers ot r1, r2, r12 from resultant
    polynomial to generate lmn numbers and gather coefficients.
    
    From earlier comment, should be used for <Psi_i|H|Psi_j> as well.
    
    Return: list of tuples representing Coefficient,L,M,N values for use in Harris, Frolov, Smith
    
    """
    integrand = ((func1 * func2 * r1 * r2 * r12)/PSI**2).expand()
    
    #display(func1, func2)
    print integrand
    
    clmn = []
    
    # Extract coefficients and powers from integrand by turning integrand into a string and using
    # regular expressions.
    integrand_list = str(integrand).split(' ')
    sign = 'pos'
    coefficient = 1
            
    for k in integrand_list:
        print k
        r1_pow = 1
        r2_pow = 1
        r12_pow = 1
        r1_search = search('r1\*\*(\d+)', k)
        r2_search = search('r2\*\*(\d+)', k)
        r12_search = search('r12\*\*(\d+)', k)
        if r1_search:
            r1_pow = int(r1_search.group(1))
        if r2_search:
            r2_pow = int(r2_search.group(1))
        if r12_search:
            r12_pow = int(r12_search.group(1))
        # Check to see if coefficient should be negative
        if k == '-':
            sign = 'neg'
            continue
        if k == '+':
            sign = 'pos'
            continue
        if k[0] != 'r':
            coefficient = float(search(r'([-+]?\d*\.\d+|\d+)', k).group(0))
        if k[0] == 'r':
            coefficient = 1
        if sign == 'neg':
            coefficient = float(coefficient * -1)
        #print 'coeff=', coefficient, 'L=', r1_pow, 'M=' , r2_pow, 'N=', r12_pow
        clmn.append((coefficient, r1_pow, r2_pow, r12_pow))
            
    return clmn
    
def gen_wavefunction(l, m, n):
    """
    Generate wave function with arbitrary number of terms.  Need to later convert
    to Hylleraas coordinate system.  
    
    """
    # Generate coefficients
    coef = Symbol('c' + str(l) + str(m) + str(n))
    wave_equation = (s)**l * (t)**(2*m) * (u)**n * PSI
    return wave_equation

def make_waves(iterables):
    """
    Generate array of wave functions over Cartesian product from 0 to iterables over length of LMN (3).  
    Somewhat of a binary counting system that increments l, m, n.
    
    Note on iterables variable.  To make it more clear to user, up to what cardinal number he/she
    is using. 1, up to '1', 2, up to '2' and so on.  Otherwise, it 2 would be 0 and 1, and that gets
    confusing.
    
    Produces (iterables + 1)^3 permutations of wave functions
    
    """
    iterables += 2
    lmn_values = []
    wave_equations = []
    lower_bound = 0
    for i in xrange(0, iterables):
        for j in list(itertools.product(range(0,i), repeat = LMN_LENGTH)):
            if j not in lmn_values:
                lmn_values.append((j))
        
    for i in lmn_values:
        l,m,n = i
        wave_equations.append(gen_wavefunction(l,m,n))
    
    return wave_equations
    
def hamiltonian_r(wfunc, z_value):
    """
    Apply Hamiltonian to wave function in r1, r2, r12 coordinate system.  z_value is atomic number, Z.
    
    """
    Z = z_value
    
    hamiltonian = (-diff(wfunc, r1,2)/2 - diff(wfunc, r2, 2)/2 - diff(wfunc, r12, 2) - \
                  ((1/r1) * diff(wfunc, r1, 1)) - \
                  ((1/r2) * diff(wfunc, r2, 1)) - \
                  ((2/r12) * diff(wfunc, r12, 1)) - \
                  (((r1**2 - r2**2 + r12**2)/(r1 * r12)) * diff(diff(wfunc, r12,1),r1,1)) - \
                  (((r2**2 - r1**2 + r12**2)/(r2 * r12)) * diff(diff(wfunc, r12,1),r2,1)) - \
                  ((Z/r1) + (Z/r2) - (1/r12)) * wfunc) * \
                  (r1 * r2 * r12)
            
    return hamiltonian
