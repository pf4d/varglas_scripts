from pylab                     import *
from fenics                    import *
from varglas.mesh.mesh_factory import MeshFactory
from varglas.data.data_factory import DataFactory
from varglas.utilities         import DataInput, DataOutput

mesh  = MeshFactory.get_antarctica_3D_10k()

bmesh   = BoundaryMesh(mesh, 'exterior')
cellmap = bmesh.entity_map(2)
pb      = CellFunction("size_t", bmesh, 0)
for c in cells(bmesh):
  if Facet(mesh, cellmap[c.index()]).normal().z() < 0:
    pb[c] = 1
submesh = SubMesh(bmesh, pb, 1)           # subset of surface mesh

out_dir  = 'output/balance_velocity/'
thklim   = 1.0

measures = DataFactory.get_ant_measures(res=900)
bedmap1  = DataFactory.get_bedmap1(thklim=thklim)
bedmap2  = DataFactory.get_bedmap2(thklim=thklim)

dm       = DataInput(measures, mesh=submesh)
db1      = DataInput(bedmap1,  mesh=submesh)
db2      = DataInput(bedmap2,  mesh=submesh)
    
db2.data['B'] = db2.data['S'] - db2.data['H']
db2.set_data_val('H', 32767, thklim)
db2.data['S'] = db2.data['B'] + db2.data['H']

S      = db2.get_projection("S",     near=True)
B      = db2.get_projection("B",     near=True)
#M      = db2.get_projection("mask",  near=True)
#Ts     = db1.get_projection("temp",  near=True)
#q_geo  = db1.get_projection("ghfsr", near=True)
adot   = db1.get_projection("acca",  near=True)
#U_ob   = dm.get_projection("U_ob",   near=True)

Q      = FunctionSpace(submesh, 'CG', 1)
QB     = FunctionSpace(submesh, 'B',  3)
W      = Q
V      = MixedFunctionSpace([Q,Q,Q])
Q2     = MixedFunctionSpace([Q,Q])

beta   = Function(Q)
Mb     = Function(Q)
Tb     = Function(Q)
u      = Function(Q)
v      = Function(Q)
w      = Function(Q)

File('test/bed/beta_s.xml') >> beta
File('test/bed/Mb_s.xml')   >> Mb
File('test/bed/T_s.xml')    >> Tb
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
H     = S - B                             # thickness
kappa = Constant(10.0)                    # smoothing radius

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

# remove areas with negative accumulation :
#adot_v = adot.vector().array()
#adot_v[adot_v < 0] = 0
#adot.vector().set_local(adot_v)
#adot.vector().apply('insert')
    
#===============================================================================
# calculate balance-velocity :

# SUPG method :
cellh   = CellSize(submesh)
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

File('test/bed/Ubar_s.xml') << U
File(out_dir + 'U.pvd')    << U
File(out_dir + 'H.pvd')    << project(H,W)
File(out_dir + 'adot.pvd') << adot
#File(out_dir + 'uhat.pvd') << project(uhat, V)



