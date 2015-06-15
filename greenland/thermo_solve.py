import sys
import varglas.solvers            as solvers
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
out_dir = 'dump/GLM/'
var_dir = 'dump/vars_ant_spacing/'
in_dir  = 'dump/ant_spacing/08/'

mesh   = Mesh(var_dir + 'mesh.xdmf')
Q      = FunctionSpace(mesh, 'CG', 1)
ff     = MeshFunction('size_t', mesh)
cf     = MeshFunction('size_t', mesh)
ff_acc = MeshFunction('size_t', mesh)

S        = Function(Q)
B        = Function(Q)
T_s      = Function(Q)
adot     = Function(Q)
mask     = Function(Q)
q_geo    = Function(Q)
u_ob     = Function(Q)
v_ob     = Function(Q)

f = HDF5File(mesh.mpi_comm(), var_dir + 'vars.h5', 'r')

f.read(S,        'S')
f.read(B,        'B')
f.read(T_s,      'T_s')
f.read(q_geo,    'q_geo')
f.read(adot,     'adot')
f.read(mask,     'mask')
f.read(ff,       'ff')
f.read(cf,       'cf')
f.read(ff_acc,   'ff_acc')
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
params['newton_solver']['relative_tolerance']       = 1e-9
params['newton_solver']['maximum_iterations']       = 30
params['newton_solver']['error_on_nonconvergence']  = False
params['newton_solver']['linear_solver']            = 'cg'
params['newton_solver']['preconditioner']           = 'hypre_amg'
parameters['form_compiler']['quadrature_degree']    = 2
parameters['form_compiler']['cpp_optimize']         = True


config = default_config()
config['log_history']                     = False
config['output_path']                     = out_dir
config['model_order']                     = 'BP'
config['use_dukowicz']                    = False
config['coupled']['on']                   = True
config['coupled']['max_iter']             = 20
config['velocity']['newton_params']       = params
config['velocity']['vert_solve_method']   = 'mumps'#'superlu_dist'
config['velocity']['calc_pressure']       = False
config['velocity']['transient_beta']      = 'stats'
config['balance_velocity']['kappa']       = 10
config['enthalpy']['on']                  = True
config['enthalpy']['solve_method']        = 'mumps'#'superlu_dist'
config['age']['on']                       = False
config['age']['use_smb_for_ela']          = True

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
model.init_U_ob(u_ob, v_ob)
model.init_E(1.0)

## solve the balance velocity for stats beta :
#F = solvers.BalanceVelocitySolver(model, config)
#F.solve()
#
#model.save_xml(model.Ubar, 'Ubar')

# use T0 and beta0 from the previous run :
model.init_T(in_dir + 'T.xml')           # temp
model.init_W(in_dir + 'W.xml')           # water
#model.init_Mb(in_dir + 'Mb.xml')         # basal melt rate
model.init_Mb('dump/test/Mb.xml')         # basal melt rate
model.init_U(in_dir + 'u.xml',
             in_dir + 'v.xml',
             in_dir + 'w.xml')           # velocity
model.init_Ubar('dump/BV/Ubar.xml')      # balance velocity
model.init_beta_stats()                  # friction

model.save_pvd(model.beta, 'beta0')

# solve the steady state solver :
F = solvers.SteadySolver(model, config)
F.solve()

File(out_dir + 'T.xml')       << model.T
File(out_dir + 'W.xml')       << model.W 
File(out_dir + 'S.xml')       << model.S
File(out_dir + 'B.xml')       << model.B
File(out_dir + 'u.xml')       << model.u 
File(out_dir + 'v.xml')       << model.v 
File(out_dir + 'w.xml')       << model.w 
File(out_dir + 'beta.xml')    << model.beta
File(out_dir + 'Mb.xml')      << model.Mb

#===============================================================================
# calculate total time to compute
s   = time() - t0
m   = s / 60.0
h   = m / 60.0
s   = s % 60
m   = m % 60
txt = "Total time to compute SteadySolver: %02d:%02d:%02d" % (h,m,s)
print_text(txt, 'red', 1)
#===============================================================================


