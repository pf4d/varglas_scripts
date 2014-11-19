from pylab  import *
from fenics import *

thklim  = 1.0
out_dir = 'output/'

mesh = Mesh('meshes/bed_mesh.xml')
Q    = FunctionSpace(mesh, 'CG', 1)
QB   = FunctionSpace(mesh, 'B',  4)
W    = Q
V    = MixedFunctionSpace([Q,Q,Q])
Q2   = MixedFunctionSpace([Q,Q])

beta  = Function(Q)
adot  = Function(Q)
Mb    = Function(Q)
T     = Function(Q)
S     = Function(Q)
B     = Function(Q)
H     = Function(Q)
u     = Function(Q)
v     = Function(Q)
w     = Function(Q)

File('test/bed/beta_s.xml') >> beta
File('test/bed/adot_s.xml') >> adot
File('test/bed/Mb_s.xml')   >> Mb
File('test/bed/T_s.xml')    >> T
File('test/bed/S_s.xml')    >> S
File('test/bed/B_s.xml')    >> B
File('test/bed/H_s.xml')    >> H
File('test/bed/u_s.xml')    >> u
File('test/bed/v_s.xml')    >> v
File('test/bed/w_s.xml')    >> w

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
z     = SpatialCoordinate(mesh)[2]        # z-coordinate of bed
H     = S - B                             # thickness
kappa = 0.0                               # smoothing radius

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
#solve(a_dSdx == L_dSdx, Nx)
#solve(a_dSdy == L_dSdy, Ny)

#dSdx_v = Nx.vector().array()
#dSdy_v = Ny.vector().array()
dSdx_v = dSdx.vector().array()
dSdy_v = dSdy.vector().array()
dSn_v  = np.sqrt(dSdx_v**2 + dSdy_v**2 + 1e-16)
dSdx.vector().set_local(-dSdx_v / dSn_v)
dSdy.vector().set_local(-dSdy_v / dSn_v)
dSdx.vector().apply('insert')
dSdy.vector().apply('insert')
dS     = as_vector([dSdx, dSdy, 0.0]) # unit normal surface slope

adot_v = adot.vector().array()
adot_v[adot_v < 0] = 0
adot.vector().set_local(adot_v)
adot.vector().apply('insert')
    
# SUPG method :
cellh   = CellSize(mesh)
phihat  = phi + cellh/2 * dot(dS, grad(phi))
phihat  = phi + cellh/(2*H) * ((H*dS[0]*phi).dx(0) + (H*dS[1]*phi).dx(1))

def L(u, uhat):
  return div(uhat)*u + dot(grad(u), uhat)

B = L(U*H, dS) * phihat * dx
a = adot * phihat * dx

U = Function(W)
solve(B == a, U)


U_v = U.vector().array()

print 'U <min,max>:', U_v.min(), U_v.max()

File(out_dir + 'U.pvd')    << U
File(out_dir + 'H.pvd')    << project(H,W)
File(out_dir + 'adot.pvd') << adot
#File(out_dir + 'uhat.pvd') << project(uhat, V)



