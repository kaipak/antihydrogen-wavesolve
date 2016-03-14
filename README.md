##  Synopsis

Wavesolve solves eigenvalue/eigenfunction on matrices describing atomic physics
characteristics including energy states and coefficients for resultant
exponential function. Integral for radial portion of wavefunction is calculated
utilizing Harris, Frolov, and Smith analytic methodology in their paper
(J. Chem Phys., Vol. 121, No. 13 Oct. 2004), generate radial integral (known as
gamma function throughout this program) with l,m,n quantum states.  

A methodology of optimizing by iteratively finding low energy states through
Powell's COBYLA algorithm that has a wrapper in mpmath library is used on
parallelotope parameters to arrive at correct alpha, beta, gamma parameters.
Further detail can be found from work from Cann and Thakkar.

## Code Example

Workflow is as follows:
* Seed parallelope parameters with a guess for A1, A2, B1, B2, G1, G2.  
* Run through wavesolve program to return an energy value.
* Check contraint conditions and ensure the following:
** a_k + b_k >= 0
** a_k + g_k >= 0
** b_k + g_k >= 0
* If not, discard guess and try something else.
* Keep going until a minimum E is found.

Eigenvectors and Eigenvalues are outputs of the wavesolve program.

## Installation
Libraries required:
mpmath
sympy
numpy
scipy

## Interface
For modifying physical evaluation parameters, you generally want to be modifying
variables in main_cobyla.py.  Here, you can change NSIZE (effectively, number of
terms in the wave function), Z (atomic number), A1-G2 (parallelotope parameters),
and PREC (precision of float).

## TO DO
* Create better interface (API).  
* Logging functionality.
* Improve matrix performance.

## Contributors
Jack Straton
Chris Keating
