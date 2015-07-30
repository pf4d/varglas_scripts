import varglas.solvers            as solvers
import varglas.model              as model
from varglas.helper               import default_config
from fenics                       import *

# get the input args :
out_dir = 'dump/ant_spacing/balance_velocity/'
in_dir  = 'dump/vars_ant_spacing/'

mesh   = Mesh(in_dir + 'mesh.xdmf')
Q      = FunctionSpace(mesh, 'CG', 1)

S      = Function(Q)
B      = Function(Q)
adot   = Function(Q)

f = HDF5File(mesh.mpi_comm(), in_dir + 'vars.h5', 'r')

f.read(S,     'S')
f.read(B,     'B')
f.read(adot,  'adot')

config = default_config()
config['output_path']               = out_dir
config['balance_velocity']['kappa'] = 10.0

model = model.Model(config)
model.set_mesh(mesh)
model.set_surface_and_bed(S, B)
model.initialize_variables()

model.init_adot(adot)

F = solvers.BalanceVelocitySolver(model, config)

F.solve()

model.save_xml(model.Ubar, 'Ubar')



