"""This demo program solves the Biharmonic equation,

    nabla^4 u(x, y) = f(x, y)

on the unit square with source f given by

    f(x, y) = 4 pi^4 sin(pi*x)*sin(pi*y)

and boundary conditions given by

    u(x, y)         = 0
    nabla^2 u(x, y) = 0
"""
from dolfin import *
#from ufl import *
#i = Index()
#j = Index()

mesh = UnitSquare(64, 64)
V = FunctionSpace(mesh, "CG", 4)
#element = TensorElement("CG", triangle, 2)

class DirichletBoundary(SubDomain):
    def inside(self, x, on_boundary):
        return on_boundary
"""
class Source(Expression):
    def eval(self, values, x):
        #values[0] = 4.0*pi**4*sin(pi*x[0])*sin(pi*x[1])
	values[0] = 10;
f = Source(degree=2)
"""
f = Constant("10.0")


# Define trial and test functions
u = TrialFunction(V)
v = TestFunction(V)



L = f*v*dx
a = inner(div(grad(u)), div(grad(v)))*dx #\
#  - inner(D(D(u,i),j), D(D(v,i),j))*dx

# Define boundary condition
u0 = Constant(0.0)
bc = DirichletBC(V, u0, DirichletBoundary())

u = Function(V)
solve(a == L, u, bc)
# Plot solution
plot(u, interactive=True)
