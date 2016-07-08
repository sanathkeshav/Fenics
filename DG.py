from dolfin import *

r = 1

tau = 1.0

gamma = 10.0*r*r

f = Expression("-exp(3*x[0])*(-9*pow(x[0],2)-3*x[0]+4)*sin(2*pi*x[1]) \
		+(x[0]-pow(x[0],2))*exp(3*x[0])*4*pow(pi,2)*sin(2*pi*x[1])")
g = Constant(0.0)

Th = UnitSquareMesh(30,30)

Vh = FunctionSpace(Th, "DG", r)

u = TrialFunction(Vh)

v = TestFunction(Vh)

n = FacetNormal(Th)

h = CellSize(Th)

a = dot(grad(u), grad(v))*dx - dot(jump(v, n), avg(grad(u)))*dS \
	- tau*dot(jump(u, n), avg(grad(v)))*dS \
	+ gamma/avg(h)*dot(jump(u, n), jump(v, n))*dS \
	- v*dot(grad(u), n)*ds - tau*u*dot(grad(v), n)*ds + gamma/h*u*v*ds

L = f*v*dx - tau*g*v*ds + gamma/h*g*v*ds

u = Function(Vh)

solve(a==L, u)

plot(u, interactive=True)
