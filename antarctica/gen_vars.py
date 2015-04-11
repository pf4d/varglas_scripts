import varglas.physical_constants as pc
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
measures = DataFactory.get_ant_measures(res=900)
bedmap1  = DataFactory.get_bedmap1(thklim=thklim)
bedmap2  = DataFactory.get_bedmap2(thklim=thklim)

mesh = MeshFactory.get_antarctica_3D_gradS_detailed()
#mesh = MeshFactory.get_antarctica_3D_gradS_crude()
#mesh = MeshFactory.get_antarctica_3D_10k()
#mesh = Mesh('meshes/ant_ultraa_gradS.xml.gz')

dm = DataInput(measures, mesh=mesh)
d1 = DataInput(bedmap1,  mesh=mesh)
d2 = DataInput(bedmap2,  mesh=mesh)

S     = d2.get_expression("S",        near=True)
B     = d2.get_expression("B",        near=True)
M     = d2.get_expression("mask",     near=True)
adot  = d1.get_expression("acca",     near=True)
T_s   = d1.get_interpolation("temp",  near=True)
q_geo = d1.get_interpolation("ghfsr", near=True)
u     = dm.get_interpolation("vx",    near=True)
v     = dm.get_interpolation("vy",    near=True)
U_ob  = dm.get_interpolation("U_ob",  near=True)

model = model.Model()
model.set_mesh(mesh)
model.calculate_boundaries(mask=M, adot=adot)
model.set_geometry(S, B, deform=True)
model.set_parameters(pc.IceParameters())

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

beta_min = interpolate(Constant(0.0), model.Q)
beta_max = interpolate(Beta_max(element = model.Q.ufl_element()), model.Q)

b_min    = interpolate(Constant(0.0), model.Q)
b_max    = interpolate(B_max(element = model.Q.ufl_element()), model.Q)

adot     = interpolate(adot, model.Q)

XDMFFile(mesh.mpi_comm(),    out_dir + 'mesh.xdmf')    << model.mesh

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

