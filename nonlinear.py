from dolfin import *

Th = UnitSquareMesh(20,20)

Vh = FunctionSpace(Th, "CG", 1)

mu = Constant(0.01)
rho = Constant(1.0)

# Define boundary conditions
# Define Dirichlet boundary (x = 0 or x = 1)


def on_boundary(x):
    return x[0] < DOLFIN_EPS or x[0] > 1.0 - DOLFIN_EPS

bc = DirichletBC(Vh, Constant("0.0"), on_boundary)


# Define variational problem

du = TrialFunction(Vh)
v = TestFunction(Vh)
u = Function(Vh)
F = (mu*inner(grad(u), grad(v))+rho*(u*u*v)-v)*dx


J = derivative(F, u, du)
# automatic computation of the Jacobian

# Compute solution


problem = NonlinearVariationalProblem(F, u, bc, J) # set nonlinear problem

solver = NonlinearVariationalSolver(problem) # set nonlinear solver

solver.solve() # solve

plot(u, interactve=True)
