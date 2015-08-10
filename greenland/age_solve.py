import varglas.solvers            as solvers
import varglas.physics            as physics
import varglas.model              as model
from varglas.data.data_factory    import DataFactory
from varglas.io                   import DataInput
from varglas.mesh.mesh_factory    import MeshFactory
from varglas.helper               import default_nonlin_solver_params, \
                                         default_config
from fenics                       import *
from scipy.io                     import loadmat

# get the input args :
out_dir = 'dump/age_run/'
var_dir = 'dump/vars_ant_spacing/'
in_dir  = 'dump/ant_spacing/08/'

mesh   = Mesh(var_dir + 'mesh.xdmf')
Q      = FunctionSpace(mesh, 'CG', 1)
ff     = MeshFunction('size_t', mesh)
cf     = MeshFunction('size_t', mesh)
ff_acc = MeshFunction('size_t', mesh)

S      = Function(Q)
B      = Function(Q)

f = HDF5File(mesh.mpi_comm(), var_dir + 'vars.h5', 'r')

f.read(S,     'S')
f.read(B,     'B')
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
config['model_order']                     = 'stokes'#'BP'
config['use_dukowicz']                    = False
config['coupled']['on']                   = False
config['coupled']['max_iter']             = 5
config['velocity']['newton_params']       = params
config['velocity']['vert_solve_method']   = 'mumps'#'superlu_dist'
config['velocity']['calc_pressure']       = False
config['velocity']['transient_beta']      = None #'stats'
config['enthalpy']['on']                  = True
config['enthalpy']['solve_method']        = 'mumps'#'superlu_dist'
config['age']['on']                       = True
config['age']['use_smb_for_ela']          = True
config['balance_velocity']['kappa']       = 5.0
config['velocity']['on']       = False
config['enthalpy']['on']       = False

model = model.Model(config)
model.set_mesh(mesh)
model.set_surface_and_bed(S, B)
model.set_subdomains(ff, cf, ff_acc)
model.initialize_variables()

model.init_U(in_dir + 'u.xml',
             in_dir + 'v.xml',
             in_dir + 'w.xml')

# solve the BP model :
F = solvers.SteadySolver(model, config)
F.solve()



