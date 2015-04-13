# Total time to compute SteadySolver: 00:04:03
# Total time to compute SteadySolver: 00:08:09

# Total time to compute SteadySolver: 00:04:58
 

import sys
import varglas.solvers            as solvers
import varglas.physical_constants as pc
import varglas.model              as model
from varglas.helper               import default_nonlin_solver_params, \
                                         default_config
from varglas.io                   import print_min_max, print_text
from fenics                       import *
from time                         import time

t0 = time()

#set_log_active(False)
#set_log_level(PROGRESS)

# get the input args :
i       = int(sys.argv[2])       # assimilation number
dir_b   = sys.argv[1] + '/0'     # directory to save

# set the output directory :
out_dir = dir_b + str(i) + '/'
in_dir  = 'dump/basin_vars_crude/'

mesh   = Mesh(in_dir + 'mesh.xdmf')
Q      = FunctionSpace(mesh, 'CG', 1)
ff     = MeshFunction('size_t', mesh)
cf     = MeshFunction('size_t', mesh)
ff_acc = MeshFunction('size_t', mesh)

S        = Function(Q)
B        = Function(Q)
T_s      = Function(Q)
adot     = Function(Q)
q_geo    = Function(Q)
beta_min = Function(Q)
beta_max = Function(Q)
b_min    = Function(Q)
b_max    = Function(Q)
u_ob     = Function(Q)
v_ob     = Function(Q)

f = HDF5File(mesh.mpi_comm(), in_dir + 'vars.h5', 'r')

f.read(S,        'S')
f.read(B,        'B')
f.read(T_s,      'T_s')
f.read(q_geo,    'q_geo')
f.read(adot,     'adot')
f.read(ff,       'ff')
f.read(cf,       'cf')
f.read(ff_acc,   'ff_acc')
f.read(beta_min, 'beta_min') 
f.read(beta_max, 'beta_max')
f.read(b_min,    'b_min')
f.read(b_max,    'b_max')
f.read(u_ob,     'u')
f.read(v_ob,     'v')
    
# specify non-linear solver parameters :
params = default_nonlin_solver_params()
#params['nonlinear_solver']                          = 'snes'
#params['snes_solver']['method']                     = 'newtonls'
#params['snes_solver']['line_search']                = 'basic'
#params['snes_solver']['error_on_nonconvergence']    = False
#params['snes_solver']['absolute_tolerance']         = 1.0
#params['snes_solver']['relative_tolerance']         = 1e-3
#params['snes_solver']['maximum_iterations']         = 20
#params['snes_solver']['linear_solver']              = 'cg'
#params['snes_solver']['preconditioner']             = 'hypre_amg'
params['nonlinear_solver']                          = 'newton'
params['newton_solver']['relaxation_parameter']     = 0.7
params['newton_solver']['relative_tolerance']       = 1e-3
params['newton_solver']['maximum_iterations']       = 30
params['newton_solver']['error_on_nonconvergence']  = False
params['newton_solver']['linear_solver']            = 'cg'
params['newton_solver']['preconditioner']           = 'hypre_amg'
parameters['form_compiler']['quadrature_degree'] = 2
#parameters['form_compiler']['optimize']          = True
parameters['form_compiler']['cpp_optimize']      = True


config = default_config()
config['output_path']                     = out_dir
config['model_order']                     = 'BP'
config['use_dukowicz']                    = True
config['coupled']['on']                   = True
config['coupled']['max_iter']             = 2
config['velocity']['newton_params']       = params
config['velocity']['viscosity_mode']      = 'full'
config['velocity']['vert_solve_method']   = 'mumps'#'superlu_dist'
config['velocity']['calc_pressure']       = False
config['enthalpy']['on']                  = True
config['enthalpy']['solve_method']        = 'mumps'#'superlu_dist'
config['age']['on']                       = False
config['age']['use_smb_for_ela']          = True
config['adjoint']['max_fun']              = 2000

model = model.Model(config)
model.set_mesh(mesh)
model.set_surface_and_bed(S, B)
model.set_subdomains(ff, cf, ff_acc)
model.initialize_variables()

model.init_viscosity_mode('full')
model.init_q_geo(model.ghf)
model.init_T_surface(T_s)
model.init_adot(adot)
model.init_U_ob(u_ob, v_ob)
model.init_E(1.0)

model.T_surface.vector()[model.shf_dofs] = 0.0
File(out_dir + 'T_test.pvd')       << model.T_surface


# use T0 and beta0 from the previous run :
if i > 0:
  model.init_T(dir_b + str(i-1) + '/T.xml')
  model.init_beta(dir_b + str(i-1) + '/beta.xml')
else:
  model.init_T(model.T_w - 30.0)
  #model.init_beta_SIA_new_slide()
  model.init_beta_SIA()

File(out_dir + 'beta0.pvd') << model.beta
File(out_dir + 'U_ob.pvd')  << model.U_ob

F = solvers.SteadySolver(model, config)
F.solve()

#===============================================================================
# print the time to compute :
tf = time()

# calculate total time to compute
s   = tf - t0
m   = s / 60.0
h   = m / 60.0
s   = s % 60
m   = m % 60
txt = "Total time to compute SteadySolver: %02d:%02d:%02d" % (h,m,s)
print_text(txt, 'red', 1)
#===============================================================================

params['newton_solver']['maximum_iterations'] = 25
config['velocity']['use_U0']                  = False
config['enthalpy']['on']                      = False
config['coupled']['on']                       = False
config['adjoint']['objective_function']       = 'linear'
config['log_history']                         = True

# invert for basal friction over grounded ice :
if i % 2 == 0:
  params['newton_solver']['relaxation_parameter']  = 1.0
  params['newton_solver']['relative_tolerance']    = 1e-10
  params['newton_solver']['maximum_iterations']    = 3
  config['adjoint']['surface_integral']            = 'grounded'
  config['adjoint']['alpha']                       = 0.0
  config['adjoint']['bounds']                      = (beta_min, beta_max)
  config['adjoint']['control_variable']            = model.beta
  model.init_viscosity_mode('linear')

# invert for enhancement over shelves :
else:
  params['newton_solver']['relaxation_parameter'] = 0.6
  params['newton_solver']['relative_tolerance']   = 1e-3
  params['newton_solver']['maximum_iterations']   = 18
  config['adjoint']['surface_integral'] = 'shelves'
  config['adjoint']['alpha']            = 0
  config['adjoint']['bounds']           = (1e-4, 1.0)
  config['adjoint']['control_variable'] = model.E_shf

A = solvers.AdjointSolver(model, config)
A.solve()

b_shf = project(model.b_shf, model.Q)
b_gnd = project(model.b_gnd, model.Q)

File(out_dir + 'T.xml')       << model.T
File(out_dir + 'W.xml')       << model.W 
File(out_dir + 'S.xml')       << model.S
File(out_dir + 'B.xml')       << model.B
File(out_dir + 'u.xml')       << model.u 
File(out_dir + 'v.xml')       << model.v 
File(out_dir + 'w.xml')       << model.w 
File(out_dir + 'beta.xml')    << model.beta
File(out_dir + 'Mb.xml')      << model.Mb
File(out_dir + 'b_shf.xml')   << b_shf
File(out_dir + 'b_gnd.xml')   << b_gnd
File(out_dir + 'E_shf.xml')   << model.E_shf

tf = time()

# calculate total time to compute
s   = tf - t0
m   = s / 60.0
h   = m / 60.0
s   = s % 60
m   = m % 60
txt = "Total time to compute: %02d:%02d:%02d" % (h,m,s)
print_text(txt, 'red', 1)



