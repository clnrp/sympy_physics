from sympy import *
#from sympy import init_session
#init_session(quiet=True)

t, q, qdot = symbols('t q \dot{q}')
L = Function('L')(q,qdot,t)
eq=diff(L,q)-diff(diff(L,t),qdot)
pprint(eq)
