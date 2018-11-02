from sympy import *
#from sympy import init_session
#init_session(quiet=True)

x = Symbol('x', real=true)
#f = Function('phi', real=true)(x)
f = Function('phi')(x)
v=Symbol('v', real=true)
phi=exp(I*v)*f

conjugate(phi)*phi
diff(conjugate(phi),x)*diff(phi,x)
r = simplify(diff(conjugate(phi),x)*diff(phi,x))
pprint(r)

