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
# calculate direction of basal water flow (down pressure gradient) :

rhoi  = 917.0                             # density of ice
rhow  = 1000.0                            # density of water
g     = 9.8                               # gravitational acceleration
z     = SpatialCoordinate(mesh)[2]        # z-coordinate of bed
      
Pw    = rhoi*g*H + rhow*g*z               # basal water pressure
      
gPx   = project(Pw.dx(0), W)
gPy   = project(Pw.dx(1), W)
gPz   = project(Pw.dx(2), W)
gPw   = as_vector([gPx, gPy, gPz])        # 2D pressure gradient
gPx_v = gPw[0].vector().array()
gPy_v = gPw[1].vector().array()
gPz_v = gPw[2].vector().array()
gPn_v = np.sqrt(gPx_v**2 + gPy_v**2 + gPz_v**2 + 1e-16)

gPn   = Function(W)                       # norm of pressure
gPn.vector().set_local(gPn_v)
gPn.vector().apply('insert')

uhat  = -gPw / gPn                        # flow direction unit vector

#===============================================================================

q    = TrialFunction(W)
phi  = TestFunction(W)
dq   = TrialFunction(W)

def L(u, uhat):
  return (uhat[0].dx(0) + uhat[1].dx(1))*u + dot(grad(u), uhat)

# SUPG method phihat :
h       = 0.0000001
U       = gPw * h
Unorm   = sqrt(dot(U, U) + 1e-16)
cellh   = CellSize(mesh)
phihat  = phi + cellh/(2*Unorm)*dot(U, grad(phi))
phihat  = phi + cellh/(2*Unorm)*((h*gPw[0]*phi).dx(0) + (h*gPw[1]*phi).dx(1))

B = L(q,uhat) * phihat * dx
a = Mb * phihat * dx
q = Function(W)

solve(B == a, q)



q_v = q.vector().array()

print 'q <min,max>:', q_v.min(), q_v.max()

File(out_dir + 'q.pvd')    << q
File(out_dir + 'Pw.pvd')   << project(Pw, Q)
File(out_dir + 'uhat.pvd') << project(uhat, V)



