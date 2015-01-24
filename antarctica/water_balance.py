from pylab  import *
from fenics import *

out_dir  = 'dump/bed/balance_water/'
in_dir   = 'dump/bed/02/'

mesh   = Mesh(in_dir + 'submesh.xdmf')
Q      = FunctionSpace(mesh, 'CG', 1)
QB     = FunctionSpace(mesh, 'B',  3)
W      = Q
V      = MixedFunctionSpace([Q,Q,Q])
Q2     = MixedFunctionSpace([Q,Q])

S      = Function(Q)
B      = Function(Q)
Mb     = Function(Q)

File(in_dir + 'S_s.xml')    >> S
File(in_dir + 'B_s.xml')    >> B
File(in_dir + 'Mb_s.xml')   >> Mb

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
# calculate direction of basal water flow (down pressure gradient) :

phi  = TestFunction(W)
Nx   = TrialFunction(W)
Ny   = TrialFunction(W)

rhoi  = 917.0                             # density of ice
rhow  = 1000.0                            # density of water
g     = 9.8                               # gravitational acceleration
H     = S - B                             # thickness
z     = SpatialCoordinate(mesh)[2]        # z-coordinate of bed
kappa = Constant(20.0)                    # smoothing radius
 
Pw    = rhoi*g*H + rhow*g*z               # basal water pressure
      
gPx   = project(Pw.dx(0), W)
gPy   = project(Pw.dx(1), W)

a_dSdx = + Nx * phi * dx \
         + (kappa*H)**2 * (phi.dx(0)*Nx.dx(0) + phi.dx(1)*Nx.dx(1)) * dx
L_dSdx = gPx * phi * dx \

a_dSdy = + Ny * phi * dx \
         + (kappa*H)**2 * (phi.dx(0)*Ny.dx(0) + phi.dx(1)*Ny.dx(1)) * dx
L_dSdy = gPy * phi * dx \

Nx = Function(W)
Ny = Function(W)
solve(a_dSdx == L_dSdx, Nx)
solve(a_dSdy == L_dSdy, Ny)

# normalize the direction vector :
gPx_v = Nx.vector().array()
gPy_v = Ny.vector().array()
gPn_v = np.sqrt(gPx_v**2 + gPy_v**2 + 1e-16)

gPx.vector().set_local(-dPx_v / gPn_v)
gPy.vector().set_local(-dPy_v / gPn_v)
gPx.vector().apply('insert')
gPy.vector().apply('insert')

gPw  = as_vector([gPx, gPy, 0.0])        # 2D pressure gradient

#===============================================================================

q    = TrialFunction(W)
phi  = TestFunction(W)
dq   = TrialFunction(W)

def L(u, uhat):
  return div(uhat)*u + dot(grad(u), uhat)

# SUPG method phihat :
h       = 1.0
np      = 0.02
U       = 1/np * (h/2)**2 * (gPw / (rhow * g))
Unorm   = sqrt(dot(U, U) + 1e-16)
cellh   = CellSize(mesh)
#phihat  = phi + cellh/(2*Unorm)*dot(U, grad(phi))
phihat  = phi + cellh/(2*Unorm)*((h*gPw[0]*phi).dx(0) + (h*gPw[1]*phi).dx(1))

w = 1/np * (q/2)**2 / (rhow * g)
B = L(q, gPw) * phihat * dx
a = Mb * phihat * dx

q = Function(W)
solve(B == a, q)

q_v = q.vector().array()

print 'q <min,max>:', q_v.min(), q_v.max()

File(out_dir + 'q.pvd')    << q
#File(out_dir + 'Pw.pvd')   << project(Pw, Q)
File(out_dir + 'U.pvd')    << project(U, V)



