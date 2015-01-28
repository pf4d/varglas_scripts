from pylab  import *
from fenics import *

out_dir  = 'dump/bed/balance_velocity/'
in_dir   = 'dump/bed/01/'

mesh   = Mesh(in_dir + 'submesh.xdmf')
Q      = FunctionSpace(mesh, 'CG', 1)
QB     = FunctionSpace(mesh, 'B',  3)
W      = Q
V      = MixedFunctionSpace([Q,Q,Q])
Q2     = MixedFunctionSpace([Q,Q])

S      = Function(Q)
B      = Function(Q)
adot   = Function(Q)

File(in_dir + 'S_s.xml')    >> S
File(in_dir + 'B_s.xml')    >> B
File(in_dir + 'adot_s.xml') >> adot

#parameters['form_compiler']['quadrature_degree'] = 3
params = {"newton_solver":
         {"linear_solver"        : 'mumps',
          "preconditioner"       : 'default',
          "maximum_iterations"   : 100,
          "relaxation_parameter" : 1.0,
          "relative_tolerance"   : 1e-3,
          "absolute_tolerance"   : 1e-16}}

class Gamma(SubDomain):
  def inside(self, x, on_boundary):
    return on_boundary
gamma = Gamma()
ff    = FacetFunction('uint', mesh)
gamma.mark(ff,1)
ds    = ds[ff]

#===============================================================================
# calculate direction of flow (down gradient) :

U    = TrialFunction(W)
phi  = TestFunction(W)
Nx   = TrialFunction(W)
Ny   = TrialFunction(W)

rho   = 917.0                             # density of ice
g     = 9.8                               # gravitational acceleration
H     = S - B                             # thickness
kappa = Constant(20.0)                    # smoothing radius

# calculate surface slope and norm :
dSdx   = project(S.dx(0), W)
dSdy   = project(S.dx(1), W)

a_dSdx = + Nx * phi * dx \
         + (kappa*H)**2 * (phi.dx(0)*Nx.dx(0) + phi.dx(1)*Nx.dx(1)) * dx
L_dSdx = rho * g * H * dSdx * phi * dx \

a_dSdy = + Ny * phi * dx \
         + (kappa*H)**2 * (phi.dx(0)*Ny.dx(0) + phi.dx(1)*Ny.dx(1)) * dx
L_dSdy = rho * g * H * dSdy * phi*dx \

Nx = Function(W)
Ny = Function(W)
solve(a_dSdx == L_dSdx, Nx)
solve(a_dSdy == L_dSdy, Ny)

# normalize the direction vector :
dSdx_v = Nx.vector().array()
dSdy_v = Ny.vector().array()
dSn_v  = np.sqrt(dSdx_v**2 + dSdy_v**2 + 1e-16)
dSdx.vector().set_local(-dSdx_v / dSn_v)
dSdy.vector().set_local(-dSdy_v / dSn_v)
dSdx.vector().apply('insert')
dSdy.vector().apply('insert')
dS     = as_vector([dSdx, dSdy, 0.0]) # unit normal surface slope


#===============================================================================
# calculate balance-velocity :

# SUPG method :
cellh   = CellSize(mesh)
#phihat  = phi + cellh/2 * dot(dS, grad(phi))  # reduce where slope is low 
phihat  = phi + cellh/(2*H) * ((H*dS[0]*phi).dx(0) + (H*dS[1]*phi).dx(1))

def L(u, uhat):
  return div(uhat)*u + dot(grad(u), uhat)

B = L(U*H, dS) * phihat * dx
a = adot * phihat * dx

U = Function(W, name='$\Vert U \Vert$')
solve(B == a, U)


U_v = U.vector().array()

print 'U <min,max>:', U_v.min(), U_v.max()

File(out_dir + 'Ubar_s.xml') << U
File(out_dir + 'Ubar_s.pvd') << U



