import varglas.physics            as physics
import varglas.physical_constants as pc
import varglas.model              as model
from varglas.helper               import default_config
from fenics                       import *
from time                         import time
from termcolor                    import colored


# get the input args :
out_dir = 'dump/antarctica/bed/balance_velocity/'
in_dir  = 'dump/antarctica/vars/'

mesh   = Mesh(in_dir + 'mesh.xdmf')
Q      = FunctionSpace(mesh, 'CG', 1)
ff     = MeshFunction('size_t', mesh)
cf     = MeshFunction('size_t', mesh)
ff_acc = MeshFunction('size_t', mesh)

S      = Function(Q)
B      = Function(Q)
adot   = Function(Q)

f = HDF5File(mesh.mpi_comm(), in_dir + 'vars.h5', 'r')

f.read(S,     'S')
f.read(B,     'B')
f.read(adot,  'adot')
f.read(ff,    'ff')
f.read(cf,    'cf')
f.read(ff_acc,'ff_acc')

model = model.Model()
model.set_mesh(mesh)
model.set_surface_and_bed(S, B)
model.set_subdomains(ff, cf, ff_acc)
model.set_parameters(pc.IceParameters())
model.initialize_variables()

model.adot = adot

config = default_config()
config['output_path']               = out_dir
config['balance_velocity']['kappa'] = Constant(10.0)
config['balance_velocity']['adot']  = adot

F = physics.VelocityBalance(model, config)

t0 = time()
F.solve()
tf = time()

File(out_dir + 'Ubar.pvd') << model.Ubar
File(out_dir + 'Ubar.xml') << model.Ubar

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



