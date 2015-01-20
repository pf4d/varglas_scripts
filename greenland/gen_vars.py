import varglas.model              as model
from varglas.mesh.mesh_factory    import MeshFactory
from varglas.data.data_factory    import DataFactory
from varglas.io                   import DataInput
from fenics                       import *
from time                         import time
from termcolor                    import colored, cprint

t0 = time()

out_dir  = 'dump/vars/'
thklim   = 1.0

# collect the raw data :
searise  = DataFactory.get_searise(thklim = thklim)
bamber   = DataFactory.get_bamber(thklim = thklim)
fm_qgeo  = DataFactory.get_gre_qgeo_fox_maule()
rignot   = DataFactory.get_gre_rignot_updated()

# define the mesh :
mesh = MeshFactory.get_greenland_detailed()

# create data objects to use with varglas :
dsr     = DataInput(searise,  mesh=mesh)
dbm     = DataInput(bamber,   mesh=mesh)
dfm     = DataInput(fm_qgeo,  mesh=mesh)
drg     = DataInput(rignot,   mesh=mesh)
    
# change the projection of the measures data to fit with other data :
drg.change_projection(dsr)

# get the expressions used by varglas :
S     = dbm.get_expression('S',        near=True)
B     = dbm.get_expression('B',        near=True)
M     = dbm.get_expression('mask',     near=True)
adot  = dsr.get_expression('adot',     near=True)
T_s   = dsr.get_interpolation('T',     near=True)
q_geo = dfm.get_interpolation('q_geo', near=True)
u     = drg.get_interpolation('vx',    near=True)
v     = drg.get_interpolation('vy',    near=True)
U_ob  = drg.get_interpolation('U_ob',  near=True)

model = model.Model()
model.set_mesh(mesh)
model.set_geometry(S, B, deform=True)
model.calculate_boundaries(M, adot=adot)

# constraints on optimization for beta :
class Beta_max(Expression):
  def eval(self, values, x):
    if M(x[0], x[1], x[2]) > 0:
      values[0] = DOLFIN_EPS
    else:
      values[0] = 4000

# constraints on optimization for b :
class B_max(Expression):
  def eval(self, values, x):
    if M(x[0], x[1], x[2]) > 0:
      values[0] = 5e6
    else:
      values[0] = DOLFIN_EPS

bedMesh  = model.get_bed_mesh()

beta_min = interpolate(Constant(0.0), model.Q)
beta_max = interpolate(Beta_max(element = model.Q.ufl_element()), model.Q)

b_min    = interpolate(Constant(0.0), model.Q)
b_max    = interpolate(B_max(element = model.Q.ufl_element()), model.Q)

adot     = interpolate(adot, model.Q)

XDMFFile(mesh.mpi_comm(),    out_dir + 'mesh.xdmf')    << model.mesh
XDMFFile(bedMesh.mpi_comm(), out_dir + 'bedMesh.xdmf') << bedMesh

# save the state of the model :
f = HDF5File(mesh.mpi_comm(), out_dir + 'vars.h5', 'w')
f.write(model.ff,     'ff')
f.write(model.cf,     'cf')
f.write(model.ff_acc, 'ff_acc')
f.write(model.S,      'S')
f.write(model.B,      'B')
f.write(T_s,          'T_s')
f.write(q_geo,        'q_geo')
f.write(adot,         'adot')
f.write(u,            'u')
f.write(v,            'v')
f.write(U_ob,         'U_ob')
f.write(beta_min,     'beta_min')
f.write(beta_max,     'beta_max')
f.write(b_min,        'b_min')
f.write(b_max,        'b_max')

tf = time()

# calculate total time to compute
s = tf - t0
m = s / 60.0
h = m / 60.0
s = s % 60
m = m % 60
if model.MPI_rank == 0:
  s    = "Total time to compute: %02d:%02d:%02d" % (h,m,s)
  text = colored(s, 'red', attrs=['bold'])
  print text

