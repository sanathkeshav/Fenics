from dolfin import *
import numpy
p=0
n = 2
maxdiff = []
order = []
for p in range(0,5):
	# Create mesh and define function space
	nx = ny = n
	mesh = UnitSquare(nx, ny)
	V = FunctionSpace(mesh, 'Lagrange', 1)

	# Define boundary conditions
	
	u0 = Expression('exp(-t)*sin(pi*x[0])*sin(pi*x[1])',t=0)

	class Boundary(SubDomain):  # define the Dirichlet boundary
	    def inside(self, x, on_boundary):
	        return on_boundary

	boundary = Boundary()
	bc = DirichletBC(V, u0, boundary)

	# Initial condition
	u_1 = interpolate(u0, V)
	#u_1 = project(u0, V)  # will not result in exact solution!

	dt = 0.1      # time step

	# Define variational problem
	u = TrialFunction(V)
	v = TestFunction(V)
	
	a = u*v*dx + dt*inner(nabla_grad(u), nabla_grad(v))*dx
	

	A = assemble(a)   # assemble only once, before the time stepping
	b = None          # necessary for memory saving assemeble call

	# Compute solution
	u = Function(V)   # the unknown at a new time level
	T = 2.0           # total simulation time
	t = dt
	
	while t <= T:
	    #print 'time =', t
	    
	    f = Expression('(2*pi*pi - 1)*exp(-t)*sin(pi*x[0])*sin(pi*x[1])', t=t)
	    L = (u_1 + dt*f)*v*dx
	    b = assemble(L, tensor=b)	
	    u0.t = t
	    bc.apply(A, b)
	    solve(A, u.vector(), b)
	
	    # Verify
	    u_e = interpolate(u0, V)
	    diff = numpy.abs(u_e.vector().array() - u.vector().array()).max()
	    
	    #print 'Max error, t=%.2f: %-10.6f' % (t, diff)	

    	    t += dt
	    u_1.assign(u)
	    
	    
	n = n*2
	maxdiff.append(diff)
	plot(u_1,interactive=True)
	
k=0
for k in range(0,5):
	print 'Max error in iteration %d is %-10.20f' % (k+1, maxdiff[k])
j=0
for j in range(0,4):
	#order[j] = ln(maxdiff[j]/maxdiff[j+1])/ln(2)
	#print 'order is %-10.6f' % (order[j])
	print ln(maxdiff[j]/maxdiff[j+1])/ln(2)

	 
