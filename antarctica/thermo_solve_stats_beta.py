import varglas.solvers            as solvers
import varglas.physics            as physics
import varglas.model              as model
from varglas.data.data_factory    import DataFactory
from varglas.io                   import DataInput
from varglas.mesh.mesh_factory    import MeshFactory
from varglas.helper               import default_nonlin_solver_params, \
                                         default_config
from fenics                       import *
from time                         import time
from termcolor                    import colored, cprint
from scipy.io                     import loadmat


t0 = time()

# get the input args :
out_dir = 'dump/linear_model_Ubar/'
var_dir = 'dump/vars_high/'
in_dir  = 'dump/high/07/'
bv_dir  = 'dump/bed/07/bv/'

mesh   = Mesh(var_dir + 'mesh.xdmf')
Q      = FunctionSpace(mesh, 'CG', 1)
ff     = MeshFunction('size_t', mesh)
cf     = MeshFunction('size_t', mesh)
ff_acc = MeshFunction('size_t', mesh)

bedmap1  = DataFactory.get_bedmap1(thklim=1.0)
d1       = DataInput(bedmap1,  mesh=mesh)

Ubar_v   = loadmat(bv_dir + 'Ubar_5.mat')

d1.data['Ubar'] = Ubar_v['map_data']

Ubar   = d1.get_expression("Ubar", near=True)

S      = Function(Q)
B      = Function(Q)
T_s    = Function(Q)
U_ob   = Function(Q)
adot   = Function(Q)
q_geo  = Function(Q)
mask   = Function(Q)

f = HDF5File(mesh.mpi_comm(), var_dir + 'vars.h5', 'r')

f.read(S,     'S')
f.read(B,     'B')
f.read(T_s,   'T_s')
f.read(U_ob,  'U_ob')
f.read(q_geo, 'q_geo')
f.read(adot,  'adot')
f.read(mask,  'mask')
f.read(ff,    'ff')
f.read(cf,    'cf')
f.read(ff_acc,'ff_acc')

# specify non-linear solver parameters :
params = default_nonlin_solver_params()
params['nonlinear_solver']                          = 'newton'
params['newton_solver']['relaxation_parameter']     = 0.7
params['newton_solver']['relative_tolerance']       = 1e-9
params['newton_solver']['maximum_iterations']       = 25
params['newton_solver']['error_on_nonconvergence']  = False
params['newton_solver']['linear_solver']            = 'cg'
params['newton_solver']['preconditioner']           = 'hypre_amg'
parameters['form_compiler']['quadrature_degree']    = 2

config = default_config()
config['output_path']                     = out_dir
config['model_order']                     = 'BP'
config['use_dukowicz']                    = False
config['coupled']['on']                   = True
config['coupled']['max_iter']             = 5
config['velocity']['newton_params']       = params
config['velocity']['vert_solve_method']   = 'mumps'#'superlu_dist'
config['velocity']['calc_pressure']       = False
config['velocity']['transient_beta']      = 'stats'
config['enthalpy']['on']                  = True
config['enthalpy']['solve_method']        = 'mumps'#'superlu_dist'
config['age']['on']                       = True
config['age']['use_smb_for_ela']          = True
config['balance_velocity']['kappa']       = 5.0

model = model.Model(config)
model.set_mesh(mesh)
model.set_surface_and_bed(S, B)
model.set_subdomains(ff, cf, ff_acc)
model.initialize_variables()

model.init_viscosity_mode('full')
model.init_mask(mask)
model.init_q_geo(model.ghf)
model.init_T_surface(T_s)
model.init_adot(adot)
model.init_E(1.0)

model.init_T(in_dir + 'T.xml')             # temp
#model.init_beta(in_dir + 'beta.xml')      # friction
model.init_W(in_dir + 'W.xml')             # water
model.init_E_shf(in_dir + 'E_shf.xml')     # enhancement
model.init_U(in_dir + 'u.xml',
             in_dir + 'v.xml',
             in_dir + 'w.xml')
model.init_Ubar(Ubar)

model.save_pvd(model.Ubar, 'Ubar')

## get the balance velocity :
#F = solvers.BalanceVelocitySolver(model, config)
#F.solve()

model.init_beta_stats('Ubar')

model.save_pvd(model.beta, 'beta')
model.save_xml(model.beta, 'beta')

# solve the BP model :
F = solvers.SteadySolver(model, config)
F.solve()

model.save_xml(model.T,  'T')
model.save_xml(model.W,  'W')
model.save_xml(model.u,  'u')
model.save_xml(model.v,  'v')
model.save_xml(model.w,  'w')
model.save_xml(model.Mb, 'Mb')

# calculate total time to compute
tf = time()
s  = tf - t0
m  = s / 60.0
h  = m / 60.0
s  = s % 60
m  = m % 60
if model.MPI_rank == 0:
  s    = "Total time to compute: %02d:%02d:%02d" % (h,m,s)
  text = colored(s, 'red', attrs=['bold'])
  print text



