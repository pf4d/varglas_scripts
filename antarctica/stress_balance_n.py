import varglas.solvers            as solvers
import varglas.model              as model
from varglas.mesh.mesh_factory    import MeshFactory
from varglas.data.data_factory    import DataFactory
from varglas.helper               import default_nonlin_solver_params, \
                                         default_config
from varglas.io                   import DataInput, DataOutput
from fenics                       import *

set_log_active(False)

thklim  = 1.0
in_dir  = 'dump/test/03/'
out_dir = 'dump/stress_n/'
var_dir = 'dump/vars/'

mesh   = Mesh(var_dir + 'mesh.xdmf')
Q      = FunctionSpace(mesh, 'CG', 1)
ff     = MeshFunction('size_t', mesh)
cf     = MeshFunction('size_t', mesh)
ff_acc = MeshFunction('size_t', mesh)

S      = Function(Q)
B      = Function(Q)

f = HDF5File(mesh.mpi_comm(), var_dir + 'vars.h5', 'r')

f.read(S,        'S')
f.read(B,        'B')
f.read(ff,       'ff')
f.read(cf,       'cf')
f.read(ff_acc,   'ff_acc')

model = model.Model()
model.set_mesh(mesh)
model.set_surface_and_bed(S, B)
model.set_subdomains(ff, cf, ff_acc)
model.initialize_variables()

model.init_beta(in_dir + 'beta.xml')
model.init_U(in_dir + 'u.xml',
             in_dir + 'v.xml',
             in_dir + 'w.xml')
model.init_T(in_dir + 'T.xml')
model.init_W(in_dir + 'W.xml')
model.init_E(1.0)

model.calc_eta()

config = default_config()
config['output_path']                      = out_dir
config['stokes_balance']['viscosity_mode'] = 'linear'
config['stokes_balance']['eta']            = model.eta
config['stokes_balance']['vert_integrate'] = False

T = solvers.StokesBalanceSolver(model, config)
T.solve()



