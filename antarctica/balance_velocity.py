import sys
src_directory = '../../statistical_modeling'
sys.path.append(src_directory)

import varglas.physical_constants as pc
import varglas.model              as model
from src.regstats                 import linRegstats
from varglas.mesh.mesh_factory    import MeshFactory
from varglas.data.data_factory    import DataFactory
from varglas.utilities            import DataInput, DataOutput
from pylab                        import *
from fenics                       import *

thklim  = 1.0
out_dir = 'test/bed/'

mesh = Mesh('meshes/bed_mesh.xml')
Q    = FunctionSpace(mesh, 'CG', 1)
QB   = FunctionSpace(mesh, 'B',  4)
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

q    = Function(Q+QB)
q_1  = Function(Q+QB)
phi  = TestFunction(Q+QB)
dq   = TrialFunction(Q+QB)
z    = SpatialCoordinate(mesh)[2]

rhoi = 917.0
rhow = 1000.0
g    = 9.8
dt   = 1.0

Pw   = rhoi*g*H + rhow*g*z

gPx  = project(Pw.dx(0), Q)
gPy  = project(Pw.dx(1), Q)
gPz  = project(Pw.dx(2), Q)

gPw  = as_vector([gPx, gPy, 0.0])

gPx_v = gPx.vector().array()
gPy_v = gPy.vector().array()
#gPz_v = gPz.vector().array()
gPn_v = np.sqrt(gPx_v**2 + gPy_v**2)# + gPz_v**2)

gPn   = Function(Q)
gPn.vector().set_local(gPn_v)
gPn.vector().apply('insert')

uhat  = -gPw / gPn

class Gamma(SubDomain):
  def inside(self, x, on_boundary):
    return on_boundary
gamma = Gamma()
ff    = FacetFunction('uint', mesh)
gamma.mark(ff,1)
ds    = ds[ff]


def L(u):
  return div(uhat)*u + dot(grad(u), uhat)

F     = Mb
qb    = Constant(0.0)
    
## SUPG method psihat :
#vnorm   = sqrt(dot(w, w) + 1e-10)
#cellh   = CellSize(mesh)
#phihat  = phi + cellh/(2*vnorm)*dot(w, phi.dx(0))

delta = + (uhat[0].dx(0) + uhat[1].dx(1)) * q * phi * dx \
        + (q.dx(0) * uhat[0] + q.dx(1) * uhat[1]) * phi * dx \
        + F * phi * dx
J     = derivative(delta, q, dq)

solve(delta == 0, q, J=J, solver_parameters=params)

q_v = q.vector().array()

print 'q <min,max>:', q_v.min(), q_v.max()




File(out_dir + 'q.pvd')    << q
File(out_dir + 'Pw.pvd')   << project(Pw, Q)
File(out_dir + 'uhat.pvd') << project(uhat, V)



