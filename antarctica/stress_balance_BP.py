import varglas.solvers            as solvers
import varglas.model              as model
from varglas.mesh.mesh_factory    import MeshFactory
from varglas.data.data_factory    import DataFactory
from varglas.helper               import default_config
from fenics                       import *

set_log_active(False)

thklim  = 1.0
in_dir  = 'dump/high/07/'
out_dir = 'dump/stress/BP/'
var_dir = 'dump/vars_high/'

mesh   = Mesh(var_dir + 'mesh.xdmf')
Q      = FunctionSpace(mesh, 'CG', 1)
ff     = MeshFunction('size_t', mesh)
cf     = MeshFunction('size_t', mesh)
ff_acc = MeshFunction('size_t', mesh)

S      = Function(Q)
B      = Function(Q)
mask   = Function(Q)

f = HDF5File(mesh.mpi_comm(), var_dir + 'vars.h5', 'r')

f.read(S,       'S')
f.read(B,       'B')
f.read(ff,      'ff')
f.read(cf,      'cf')
f.read(ff_acc,  'ff_acc')
f.read(mask,    'mask')

config = default_config()
config['output_path']                      = out_dir
config['model_order']                      = 'BP'
config['use_dukowicz']                     = False
config['stokes_balance']['vert_integrate'] = True

model = model.Model(config)
model.set_mesh(mesh)
model.set_surface_and_bed(S, B)
model.set_subdomains(ff, cf, ff_acc)
model.initialize_variables()
model.init_viscosity_mode('full')

model.init_mask(mask)
model.init_beta(in_dir + 'beta.xml')
model.init_U(in_dir + 'u.xml',
             in_dir + 'v.xml',
             in_dir + 'w.xml')
model.init_T(in_dir + 'T.xml')
model.init_W(in_dir + 'W.xml')
model.init_E_shf(in_dir + 'E_shf.xml')
model.init_E_gnd(in_dir + 'E_gnd.xml')

E_shf_v = model.E_shf.vector().array()
E_shf_v[E_shf_v < 1e-2] = 1e-2
model.assign_variable(model.E_shf, E_shf_v)

T = solvers.StokesBalanceSolver(model, config)
T.solve()

model.save_pvd(model.eta,    'eta')
model.save_xml(model.tau_id, 'tau_id')
model.save_xml(model.tau_jd, 'tau_jd')
model.save_xml(model.tau_ib, 'tau_ib')
model.save_xml(model.tau_jb, 'tau_jb')
model.save_xml(model.tau_ip, 'tau_ip')
model.save_xml(model.tau_jp, 'tau_jp')
model.save_xml(model.tau_ii, 'tau_ii')
model.save_xml(model.tau_ij, 'tau_ij')
model.save_xml(model.tau_iz, 'tau_iz')
model.save_xml(model.tau_ji, 'tau_ji')
model.save_xml(model.tau_jj, 'tau_jj')
model.save_xml(model.tau_jz, 'tau_jz')



