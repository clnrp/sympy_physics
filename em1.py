from sympy import *
from sympy import init_session
#init_session(quiet=True)

t,c = symbols('t c')
Ex,Ey,Ez = symbols('E_x E_y E_z')
Bx,By,Bz = symbols('B_x B_y B_z')

t,x,y,z = symbols('t x y z')
n_mu_nu=Matrix([[1,0,0,0],[0,-1,0,0],[0,0,-1,0],[0,0,0,-1]])
x_mu=Matrix([t,x,y,z])
x_nu=n_mu_nu*x_mu

F=Matrix([[0,-Ex/c,-Ey/c,-Ez/c],[Ex/c,0,-Bz,By],[Ey/c,Bz,0,-Bx],[Ez/c,-By,Bx,0]])
sum([sum([n_mu_nu[0,i]*n_mu_nu[3,j]*F[i,j] for i in range(4)]) for j in range(4)])

F__mn=Matrix([[sum([sum([n_mu_nu[i,k]*n_mu_nu[j,l]*F[k,l] for k in range(4)]) for l in range(4)]) for j in range(4)] for i in range(4)])

# para o eixo x
th=Symbol('theta')
L=Matrix([[cosh(th),-sinh(th),0,0],[-sinh(th),cosh(th),0,0],[0,0,1,0],[0,0,0,1]])

Ft=Matrix([[sum([sum([L[i,k]*L[j,l]*F[k,l] for k in range(4)]) for l in range(4)]) for j in range(4)] for i in range(4)])

# apenas campo eletrico no eixo z
F=Matrix([[0,0,0,Ez/c],[0,0,0,0],[0,0,0,0],[0,0,0,0]])
Fl=Matrix([[sum([sum([L[i,k]*L[j,l]*F[k,l] for k in range(4)]) for l in range(4)]) for j in range(4)] for i in range(4)])
pprint(Fl)
