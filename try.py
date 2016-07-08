from dolfin import *
import sys

def compute(nx, ny, degree):
	mesh = UnitSquareMesh(nx, ny)
	V = FunctionSpace(mesh, "Lagrange", degree)
	
	def u0_boundary(x, on_boundary):
        	return on_boundary
	
	bc = DirichletBC(V, Constant(0.0), u0_boundary)
	
	# Exact solution
    	omega = 1.0
    	u_e = Expression('sin(omega*pi*x[0])*sin(omega*pi*x[1])',omega=omega)
	
	# Define variational problem
	u = TrialFunction(V)
	v = TestFunction(V)
	f = 2*pi**2*omega**2*u_e
	a = inner(nabla_grad(u), nabla_grad(v))*dx
	L = f*v*dx

	# Compute solution
	u = Function(V)
	solve(a == L, u, bc)
	#plot(u)
	#interactive()
	
	u_e_V = interpolate(u_e, V)
    	E5 = abs(u_e_V.vector().array() - u.vector().array()).max()
	errors = {'infinity norm (of dofs)': E5}
	return errors

degree = 2 #int(sys.argv[1])
h = []  # element sizes
E = []  # errors
for nx in [10, 20, 40, 80, 160]:
    h.append(1.0/nx)
    E.append(compute(nx, nx, degree))  # list of dicts

# Convergence rates
from math import log as ln  # log is a dolfin name too
error_types = E[0].keys()
for error_type in sorted(error_types):
    print '\nError norm based on', error_type
    for i in range(1, len(E)):
        Ei   = E[i][error_type]  # E is a list of dicts
        Eim1 = E[i-1][error_type]
        r = ln(Ei/Eim1)/ln(h[i]/h[i-1])
        print 'h=%8.2E E=%8.2E r=%.5f' % (h[i], Ei, r)
