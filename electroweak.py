from sympy import *
from sympy import LeviCivita
#from sympy import init_session
#init_session(quiet=True)

v, H = symbols('v, H', real=true)
g1, g2 = symbols("g, g'")
x, gamma = symbols("x, gamma__mu")

w_ma, w_me, C = symbols("{w__+}, {w__-}, C_mu")
B1, B2, B3 = symbols("{B_mu_1}, {B_mu_2}, {B_mu_3}")
A, Z, theta = symbols("A_mu, Z_mu, theta")

# setor leptons
ll, lr, vl = symbols("l_L, l_R, nu_L")

l=Matrix([[vl], [ll]])
Dl=Matrix([[diff(vl(x),x)], [diff(ll(x),x)]])+I*g1/2*Matrix([[B3, sqrt(2)*w_ma], [sqrt(2)*w_me, -B3]])*l-I*g2*C*l
Dlr=diff(lr(x),x)-2*I*g2*C*lr

L=(l.adjoint()*I*gamma*Dl)[0]+lr.conjugate()*I*gamma*Dlr
L=L.simplify()
L=L.expand(basic=true)
L.factor(vl.conjugate())
pprint(L)

# setor de higgs
DphiDag=Matrix([[0, diff(v+H(x),x)]])-I*g1/2*Matrix([[0, v+H]])*Matrix([[B3, sqrt(2)*w_ma], [sqrt(2)*w_me, -B3]])-I*g2*Matrix([[0, v*C+H*C]])

Dphi=Matrix([[0], [diff(v+H(x),x)]])+I*g1/2*Matrix([[B3, sqrt(2)*w_ma], [sqrt(2)*w_me, -B3]])*Matrix([[0], [v+H]])+I*g2*Matrix([[0], [v*C+H*C]])


L=DphiDag*Dphi
L=L[0].simplify()
L=L.expand(basic=true)
L.factor(v)
pprint(L)

# massa do boson Z
M=Matrix([[g1**2/4, -g1*g2/2], [-g1*g2/2, g2**2]])
M.det()
valores=M.eigenvals()
a_vet=Matrix([[Z], [A]])

a_va=Matrix([[valores.items()[0][0], 0], [0, valores.items()[1][0]]])
v**2/2*a_vet.transpose()*a_va*a_vet

L=L.subs(B3, -Z*g1/2-A*g2).expand(basic=true)
L=L.subs(C, Z*g2-A*g1/2).expand(basic=true)

MO=1/sqrt(g1**2/4+g2**2)*Matrix([[-g1/2, g2], [-g2, -g1/2]])
pprint(L)

#setor gauge FF
Bnu1, Bnu2, Bnu3 = symbols("{B_nu_1}, {B_nu_2}, {B_nu_3}")
Bmu=[B1,B2,B3]
Bnu=[Bnu1,Bnu2,Bnu3]

a=1
G=sum([sum([LeviCivita(a, b, c)*Bmu[b]*Bnu[c] for b in range(3)]) for c in range(3)])
pprint(G)



