from sympy import *
from sympy import init_session
#init_session(quiet=True)

xm, xn = symbols('x_mu x_nu')
x, m = symbols('x m')
n = symbols('n_mu_nu')
phi1 = Function('phi')(x)
phi2 = Function('phi__*')(x)
L = diff(phi1,x)*diff(phi2,x)+m**2*phi1*phi2 # lagrangiana do campo escalar
pprint(L)

field=[phi1, phi2]

eqs=[diff(L,field[i])-diff(diff(L,diff(field[i],x)),x) for i in range(2)] # equacoes de campo
pprint(eqs) 

deltPhi=[I*phi1,-I*phi2]
deltX=0

J=sum([diff(L,diff(field[i],x))*(diff(field[i],x)*deltX - deltPhi[i]) for i in range(2)])-L*deltX # corrente
pprint(J) 

T=sum([diff(L,diff(field[i],x))*diff(field[i],x) for i in range(2)])-n*L # tensor energia-momento
pprint(T) 









