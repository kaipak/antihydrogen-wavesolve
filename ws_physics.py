from IPython.display import display
import numpy as np
import scipy.misc as sc
import sympy as sym
import mpmath as mpm
import itertools
import re

r1 = sym.Symbol('r1')
r2 = sym.Symbol('r2')
r12 = sym.Symbol('r12')

alpha = sym.Symbol('alpha')
beta = sym.Symbol('beta')
gamma = sym.Symbol('gamma')

# Hylleras coordinate stuff (transform later?)
s = r1 + r2
t = r1 - r2
u = r12

LMN_LENGTH = 3 # determine Cartesian product for generating wave funcs.

def static_params(a1, a2, b1, b2, g1, g2, eta, zed, nsize, prec=16, dps=8): 
    """
    Globally set STATIC parameters for later function calls.  Something of
    a poor man's OOP and an effort to curtail the dreaded creep of
    spaghetti code.  This is generally frowned upon practice, but done to
    to make setting of global parameters more user friendly from main.  
    Dynamic parameters such as alpha, beta should NOT be defined here.

    Keyword arguments:
    a1, a2, b1, etc. - parameters required for particular problem that do 
                       not change.

    Z - atomic number

    """
    global A1
    global A2
    global B1
    global B2
    global ETA
    global Z
    Z = zed
    
    # Specific to p10diag problem from JS
    A1 = a1
    A2 = a2
    B1 = b1
    B2 = b2
    ETA = eta
    
    global PREC
    global DPS
    global NSIZE
    
    NSIZE = nsize
    PREC = prec
    DPS = dps

def test():
    return mpm.pi

def hfs_gamma(l, m, n, alpha, beta, gamma):
    """
    Using Harris, Frolov, and Smith paper (J. Chem Phys., Vol. 121, No. 13
    Oct. 2004), generate radial integral or gamma function 
    (alpha=a, beta=b, gamma=g) using l,m,n quantum states.
    
    Keyword arguments:
    l, m, n -- parameters corresponding to powers of s, t, u respectively.
    a, b, g -- independently derived parameters unique to problem being solved.
    
    Returns:
    float(16,32,64,etc) -- A numerical value of the gamma function
    
    """
    a = alpha
    b = beta
    g = gamma

    fact_coef = mpm.mpf(2 * mpm.factorial(l) * mpm.factorial(m) \
                * mpm.factorial(n))
        
    x = mpm.mpf(a + b)
    y = mpm.mpf(a + g)
    z = mpm.mpf(b + g)
        
    #big_gamma = sym.Mul(0)
    big_gamma = mpm.mpf(0)
    
    for l_prime in range(0,l+1):
        for m_prime in range(0,m+1):     
            for n_prime in range(0,n+1):
                big_gamma = big_gamma + (mpm.binomial((m - m_prime + l_prime), l_prime) * \
                 mpm.binomial((l - l_prime + n_prime), n_prime) * \
                 mpm.binomial((n - n_prime + m_prime), m_prime) / \
                (x**(m-m_prime+l_prime+1)*y**(l-l_prime+n_prime+1)*z**(n-n_prime+m_prime+1)))
    
    big_gamma = mpm.mpf(fact_coef*big_gamma)
    
    return big_gamma

def extract_clmn(func1, func2, alpha, beta, gamma):
    """ Generalized function to extract coefficient, r1, r2, and r12 

    Values from <Psi_i|Psi_j> or <Psi_i|Operator|Psi_j> for
    primary purpose of being applied to Harris, Frolov, and Smith Equation.
    
    <Psi_i|Psi_j> = integral(Psy_i*Psy_j*r1*r2*r12dr1dr2dr12).  Pull out powers
    of r1, r2, r12 from resultant polynomial to generate lmn numbers and gather
    coefficients.
    
    From earlier comment, should be used for <Psi_i|H|Psi_j> as well.
    
    Keyword Arguments:
    func1, func2 -- wave functions (bra according to Dirac bra-ket notation)
    
    Returns:
    clmn -- list of tuples representing Coefficient,L,M,N values for use in
            Harris, Frolov, Smith
    
    """
    exponential = -((alpha * r1) + (beta * r2) + (gamma * r12))
    psi = sym.exp(exponential)

    integrand = ((func1 * func2 * r1 * r2 * r12)/psi**2).expand()
    
    # display(func1, func2)
    # display(integrand)
    
    clmn = []
    
    # Extract coefficients and powers from integrand by turning 
    # integrand into a string and using regular expressions.
    integrand_list = str(integrand).split(' ')
    sign = 'pos'
    coefficient = 1
            
    for k in integrand_list:
        
        # Cascade through possible r values to determine power.  Should possibly
        # make more elegant later        
        r1_pow = 0
        r2_pow = 0
        r12_pow = 0
        
        # check to see if r values exist
        r1_search = re.search('(r1\*|r1\s)', k)
        r2_search = re.search('r2', k)
        r12_search = re.search('r12', k)
        
        if r1_search:
            r1_pow = 1
        if r2_search:
            r2_pow = 1    
        if r12_search:
            r12_pow = 1        
        
        # check to see if r values have powers
        r1_powsearch = re.search('r1\*\*(\d+)', k)
        r2_powsearch = re.search('r2\*\*(\d+)', k)
        r12_powsearch = re.search('r12\*\*(\d+)', k)
        
        if r1_powsearch:
            r1_pow = int(r1_powsearch.group(1))
        if r2_powsearch:
            r2_pow = int(r2_powsearch.group(1))
        if r12_powsearch:
            r12_pow = int(r12_powsearch.group(1))
        # Check to see if coefficient should be negative
        if k == '-':
            sign = 'neg'
            continue
        if k == '+':
            sign = 'pos'
            continue
        if k[0] != 'r':
            coefficient = sym.N(re.search(r'([-+]?\d*\.\d+|\d+)', k).group(0), DPS)
        if k[0] == 'r':
            coefficient = 1
        if sign == 'neg':
            coefficient = sym.N(coefficient * -1, DPS)
        # print 'coeff =', coefficient, 'L =', r1_pow, 'M =' , r2_pow, 'N =', r12_pow
        clmn.append((coefficient, r1_pow, r2_pow, r12_pow))
            
    return clmn

def get_qstate(bra, ket, alpha, beta, gamma):
    """Calculate inner-product of wavefunctions

    Entry point to get quantum state of bra and ket utilizing Frolov, 
    Smith, and Harris analytic equation for radial integral.  If 
    operator applied (e.g.,, Hamiltonian, that should happen to ket 
    before this function is called

    """
    clmns = extract_clmn(bra, ket, alpha, beta, gamma)
    innerprod = 0
    
    for i in clmns:
        c,l,m,n = i
        innerprod += sym.N(c * hfs_gamma(l, m, n, alpha, beta, gamma), DPS)
    
    innerprod = sym.N(8 * sym.pi**2 * innerprod, DPS)
    return innerprod
    
def gen_wavefunction(l, m, n, alpha, beta, gamma):
    """
    Generate wave function with arbitrary number of terms.  Need to 
    later convert to Hylleraas coordinate system.  
    
    """
    
    exponential = -((alpha * r1) + (beta * r2) + (gamma * r12))
    psi = sym.exp(exponential)
    # Generate coefficients
    # coef = sym.Symbol('c' + str(l) + str(m) + str(n))


    wave_equation = (s)**l * (t)**(m) * (u)**n * psi
    return wave_equation

def make_waves(iterables):
    """
    Generate array of wave functions over Cartesian product from 0 to 
    iterables over length of LMN (3). Somewhat of a binary counting 
    system that increments l, m, n.
    
    Note on iterables variable.  To make it more clear to user, up to 
    what cardinal number he/she is using. 1, up to '1', 2, up to '2' 
    and so on.  Otherwise, it 2 would be 0 and 1, and that gets
    confusing.
    
    Produces (iterables + 1)^3 permutations of wave functions
    
    """
    iterables += 2
    lmn_values = []
    wave_equations = []
    lower_bound = 0
    
    """
    for i in xrange(0, iterables):
        for j in list(itertools.product(range(0,i), repeat = LMN_LENGTH)):
            if j not in lmn_values:
                lmn_values.append(j)
    """
    for i in list(itertools.product(range(0,iterables), repeat = LMN_LENGTH)):
        l,m,n = i
        print i
        wave_equations.append(gen_wavefunction(l,m,n))

    """
    for i in lmn_values:
        l,m,n = i
        print i
        wave_equations.append(gen_wavefunction(l,m,n))
    """
    
    return wave_equations

def hamiltonian_r(wfunc, alpha, beta, gamma):
    """Calcualte Hamiltonian to wave function in r1, r2, r12
    
    """
    hamiltonian = (-sym.diff(wfunc, r1,2)/2 - sym.diff(wfunc, r2, 2)/2 - sym.diff(wfunc, r12, 2) - \
                  ((1/r1) * sym.diff(wfunc, r1, 1)) - \
                  ((1/r2) * sym.diff(wfunc, r2, 1)) - \
                  ((2/r12) * sym.diff(wfunc, r12, 1)) - \
                  (((r1**2 - r2**2 + r12**2)/(2 * r1 * r12)) * sym.diff(sym.diff(wfunc, r1,1),r12,1)) - \
                  (((r2**2 - r1**2 + r12**2)/(2 * r2 * r12)) * sym.diff(sym.diff(wfunc, r2,1),r12,1)) - \
                  ((Z/r1) + (Z/r2) - (1/r12)) * wfunc)
    
    return hamiltonian

def psif_H_psii(alpha_n, beta_n, gamma_n, alpha_m, beta_m, gamma_m):
    """Return innner product of <psi_f| H | psi_i>"""

    alphanm = alpha_n + alpha_m
    betanm = beta_n + beta_m
    gammanm = gamma_n + gamma_m

    innerprod = sym.N(8 * sym.pi**2 * (-.5 * (alpha_m**2 + beta_m**2 + (2*gamma_m**2)) * hfs_gamma(1, 1, 1, alphanm, betanm, gammanm) +
                                   ((alpha_m - Z) * hfs_gamma(0, 1, 1, alphanm, betanm, gammanm)) +
                                   ((beta_m - Z) * hfs_gamma(1, 0, 1, alphanm, betanm, gammanm)) +
                                   ((1 + (2 * gamma_m)) * hfs_gamma(1, 1, 0, alphanm, betanm, gammanm)) -
                                   (((gamma_m * alpha_m)/2) * hfs_gamma(2, 1, 0, alphanm, betanm, gammanm)) -
                                   (((gamma_m * beta_m)/2) * hfs_gamma(1, 2, 0, alphanm, betanm, gammanm)) -
                                   (((gamma_m * alpha_m)/2) * hfs_gamma(0, 1, 2, alphanm, betanm, gammanm)) -
                                   (((gamma_m * beta_m)/2) * hfs_gamma(1, 0, 2, alphanm, betanm, gammanm)) +
                                   (((gamma_m * alpha_m)/2) * hfs_gamma(0, 3, 0, alphanm, betanm, gammanm)) +
                                   (((gamma_m * beta_m)/2) * hfs_gamma(3, 0, 0, alphanm, betanm, gammanm))), DPS)

    return innerprod

def thakar_smith_param(L1, L2, num_rad):
    """ Parameter values as defined by Thakkar and Vedene H.Smith, Jr.
        Phys Rev A 15, 1 (1977). Defined for alpha, beta, and gamma in
        paper.
    
    Returns:
    list -- array of parameter values
    
    """
    params = []
    
    for i in xrange(1, NSIZE + 1):
        # decompose number extract fractional part
        # num_comps = np.modf((i*(i+1)*np.sqrt(num_rad))/2)
        # frac = num_comps[0]

        frac = mpm.mpf(((i*(i+1)*np.sqrt(num_rad))/2)%1)

        num = mpm.mpf(ETA * (((L2 - L1)*frac) + L1))
        params.append(num)
    
    return params