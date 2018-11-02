from sympy import *
from sympy import init_session
#init_session(quiet=True)

# dirac equation
x, m, e = symbols('x m e')
gm = symbols('gamma__mu')
psi1 = Function('psi')(x)
psi2 = Function('psibar')(x)
A = Function('A_mu')(x)

field=[psi1, psi2]

L = I*psi2*gm*diff(psi1,x) -m*psi1*psi2 -e*psi1*gm*psi2*A
deltPhi=[I*psi1,-I*psi2]
deltX=0
pprint(L)
