"""Smoke test: Poisson problem on the unit cube with a known exact solution."""
import numpy as np
import dolfinx
import ufl
from mpi4py import MPI
from dolfinx import fem, mesh
from dolfinx.fem.petsc import LinearProblem

print("dolfinx version:", dolfinx.__version__)
msh = mesh.create_unit_cube(MPI.COMM_WORLD, 8, 8, 8)
V = fem.functionspace(msh, ("Lagrange", 2))
x = ufl.SpatialCoordinate(msh)
u_exact = 1 + x[0] ** 2 + 2 * x[1] ** 2 + 3 * x[2] ** 2
f = fem.Constant(msh, -12.0)
u_bc = fem.Function(V)
u_bc.interpolate(lambda X: 1 + X[0] ** 2 + 2 * X[1] ** 2 + 3 * X[2] ** 2)
msh.topology.create_connectivity(2, 3)
facets = mesh.exterior_facet_indices(msh.topology)
bc = fem.dirichletbc(u_bc, fem.locate_dofs_topological(V, 2, facets))
u, v = ufl.TrialFunction(V), ufl.TestFunction(V)
problem = LinearProblem(ufl.inner(ufl.grad(u), ufl.grad(v)) * ufl.dx, f * v * ufl.dx, bcs=[bc],
                        petsc_options={"ksp_type": "preonly", "pc_type": "lu"}, petsc_options_prefix="smoke_")
uh = problem.solve()
err = np.sqrt(msh.comm.allreduce(fem.assemble_scalar(fem.form((uh - u_exact) ** 2 * ufl.dx)), op=MPI.SUM))
print(f"L2 error vs exact quadratic solution: {err:.3e}")
assert err < 1e-8
print("Smoke test passed")
