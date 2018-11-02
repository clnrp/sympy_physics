from sympy import *
#from sympy import init_session
#init_session(quiet=True)

#taylor series f(x+h)

x, h = symbols('x, h')
f = Function('f')
taylor = sum(h**i/factorial(i) * f(x).diff(x, i) for i in range(4))
#Taylor = series(f(x), x, x0=h, n=3)

pprint(taylor)
