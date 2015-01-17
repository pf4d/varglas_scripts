import sys
import varglas.solvers            as solvers
import varglas.physical_constants as pc
import varglas.model              as model
from varglas.mesh.mesh_factory    import MeshFactory
from varglas.data.data_factory    import DataFactory
from varglas.helper               import default_nonlin_solver_params, \
                                         default_config
from varglas.utilities            import DataInput, DataOutput
from fenics                       import *
from time                         import time
from termcolor                    import colored, cprint

t0 = time()

# get the input args :
i = int(sys.argv[2])           # assimilation number
dir_b = sys.argv[1] + '/0'     # directory to save

# set the output directory :
out_dir = dir_b + str(i) + '/'
in_dir  = 'dump/vars/'

mesh   = Mesh(in_dir + 'mesh.xdmf')
Q      = FunctionSpace(mesh, 'CG', 1)
ff     = MeshFunction('size_t', mesh)
cf     = MeshFunction('size_t', mesh)
ff_acc = MeshFunction('size_t', mesh)

S        = Function(Q)
B        = Function(Q)
T_s      = Function(Q)
U_ob     = Function(Q)
adot     = Function(Q)
q_geo    = Function(Q)
beta_min = Function(Q)
beta_max = Function(Q)
b_min    = Function(Q)
b_max    = Function(Q)
u        = Function(Q)
v        = Function(Q)

f = HDF5File(mesh.mpi_comm(), in_dir + 'vars.h5', 'r')

f.read(S,        'S')
f.read(B,        'B')
f.read(T_s,      'T_s')
f.read(U_ob,     'U_ob')
f.read(q_geo,    'q_geo')
f.read(adot,     'adot')
f.read(ff,       'ff')
f.read(cf,       'cf')
f.read(ff_acc,   'ff_acc')
f.read(beta_min, 'beta_min') 
f.read(beta_max, 'beta_max')
f.read(b_min,    'b_min')
f.read(b_max,    'b_max')
f.read(u,        'u')
f.read(v,        'v')

model = model.Model()
model.set_mesh(mesh)
model.set_surface_and_bed(S, B)
model.set_subdomains(ff, cf, ff_acc)
model.initialize_variables()

# specifify non-linear solver parameters :
nonlin_solver_params = default_nonlin_solver_params()
nonlin_solver_params['newton_solver']['relaxation_parameter']    = 0.7
nonlin_solver_params['newton_solver']['relative_tolerance']      = 1e-3
nonlin_solver_params['newton_solver']['maximum_iterations']      = 16
nonlin_solver_params['newton_solver']['error_on_nonconvergence'] = False
nonlin_solver_params['newton_solver']['linear_solver']           = 'cg'
nonlin_solver_params['newton_solver']['preconditioner']          = 'hypre_amg'
parameters['form_compiler']['quadrature_degree']                 = 2

config = default_config()
config['output_path']                     = out_dir
config['coupled']['on']                   = True
config['coupled']['max_iter']             = 2
config['velocity']['newton_params']       = params
config['velocity']['approximation']       = 'fo'#'stokes'
config['velocity']['viscosity_mode']      = 'full'
config['velocity']['use_T0']              = True
config['velocity']['use_U0']              = False
config['velocity']['use_beta0']           = False
config['velocity']['T0']                  = model.T_w - 30.0
config['velocity']['init_beta_from_U_ob'] = True
config['velocity']['U_ob']                = U_ob
config['enthalpy']['on']                  = True
config['enthalpy']['T_surface']           = T_s
config['enthalpy']['q_geo']               = q_geo
config['age']['on']                       = False
config['age']['use_smb_for_ela']          = True
config['adjoint']['max_fun']              = 75

# use T0 and beta0 from the previous run :
if i > 0:
  config['velocity']['init_beta_from_U_ob'] = False
  config['velocity']['use_beta0']           = True
  config['velocity']['use_T0']              = True
  config['velocity']['use_U0']              = True
  config['velocity']['beta0']               = dir_b + str(i-1) + '/beta.xml'
  config['velocity']['T0']                  = dir_b + str(i-1) + '/T.xml'
  config['velocity']['u0']                  = dir_b + str(i-1) + '/u.xml'
  config['velocity']['v0']                  = dir_b + str(i-1) + '/v.xml'
  config['velocity']['w0']                  = dir_b + str(i-1) + '/w.xml'

F = solvers.SteadySolver(model, config)
File(out_dir + 'beta0.pvd') << model.beta

t0 = time()
F.solve()
t1 = time()

params['newton_solver']['maximum_iterations'] = 25
config['velocity']['init_beta_from_U_ob']     = False
config['velocity']['use_T0']                  = False
config['velocity']['use_U0']                  = False
config['velocity']['use_beta0']               = False
config['enthalpy']['on']                      = False
config['coupled']['on']                       = False

params['newton_solver']['relaxation_parameter']  = 1.0
config['velocity']['viscosity_mode']             = 'linear'
config['velocity']['eta_shf']                    = model.eta_shf
config['velocity']['eta_gnd']                    = model.eta_gnd
config['adjoint']['surface_integral']            = 'grounded'
config['adjoint']['alpha']                       = 0
config['adjoint']['bounds']                      = (beta_min, beta_max)
config['adjoint']['control_variable']            = model.beta

A = solvers.AdjointSolver(model, config)
A.set_target_velocity(u=u, v=v)
#uf = dir_b + str(i-1) + '/u.xml'
#vf = dir_b + str(i-1) + '/v.xml'
#wf = dir_b + str(i-1) + '/w.xml'
#A.set_velocity(uf, vf, wf)
t2 = time()
A.solve()
t3 = time()

eta   = project(model.eta, model.Q)

File(out_dir + 'T.xml')       << model.T
File(out_dir + 'S.xml')       << model.S
File(out_dir + 'B.xml')       << model.B
File(out_dir + 'u.xml')       << model.u 
File(out_dir + 'v.xml')       << model.v 
File(out_dir + 'w.xml')       << model.w 
File(out_dir + 'beta.xml')    << model.beta
File(out_dir + 'Mb.xml')      << model.Mb
File(out_dir + 'eta.xml')     << eta

#XDMFFile(mesh.mpi_comm(), out_dir + 'mesh.xdmf')   << model.mesh
#
## save the state of the model :
#if i !=0: rw = 'a'
#else:     rw = 'w'
#f = HDF5File(mesh.mpi_comm(), out_dir + 'floating_shelves_0'+str(i)+'.h5', rw)
#f.write(model.mesh,  'mesh')
#f.write(model.beta,  'beta')
#f.write(model.Mb,    'Mb')
#f.write(model.T,     'T')
#f.write(model.S,     'S')
#f.write(model.B,     'B')
#f.write(model.U,     'U')
#f.write(model.eta,   'eta')

# calculate total time to compute
s = (t1 - t0) + (t3 - t2)
m = s / 60.0
h = m / 60.0
s = s % 60
m = m % 60
if model.MPI_rank == 0:
  s    = "Total time to compute: %02d:%02d:%02d" % (h,m,s)
  text = colored(s, 'red', attrs=['bold'])
  print text



