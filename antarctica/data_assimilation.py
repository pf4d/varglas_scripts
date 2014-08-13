# beta:         F =   1292801664727.8921,  Total time to compute: 01:03:21
# beta^2:       F =   1714728897787.7820,  Total time to compute: 01:16:09
# r=1:          F =   2646114341918.4277,  Total time to compute: 01:15:44
# r=1, beta^2:  F =   5170367358780.6270,  Total time to compute: 01:17:37


import sys
import varglas.solvers            as solvers
import varglas.physical_constants as pc
import varglas.model              as model
from varglas.mesh.mesh_factory    import MeshFactory
from varglas.data.data_factory    import DataFactory
from varglas.helper               import default_nonlin_solver_params
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
in_dir  = 'vars/'

set_log_active(True)

thklim = 1.0

measures  = DataFactory.get_ant_measures(res=900)
bedmap1   = DataFactory.get_bedmap1(thklim=thklim)
bedmap2   = DataFactory.get_bedmap2(thklim=thklim)

mesh = MeshFactory.get_antarctica_3D_gradS_detailed()
#mesh = MeshFactory.get_antarctica_3D_gradS_crude()

dm  = DataInput(measures, mesh=mesh)
db1 = DataInput(bedmap1,  mesh=mesh)
db2 = DataInput(bedmap2,  mesh=mesh)

db2.data['B'] = db2.data['S'] - db2.data['H']
db2.set_data_val('H', 32767, thklim)
db2.data['S'] = db2.data['B'] + db2.data['H']

H      = db2.get_nearest_expression("H")
S      = db2.get_nearest_expression("S")
B      = db2.get_nearest_expression("B")
M      = db2.get_nearest_expression("mask")
T_s    = db1.get_nearest_expression("srfTemp")
q_geo  = db1.get_nearest_expression("q_geo")
adot   = db1.get_nearest_expression("adot")
U_ob   = dm.get_projection("U_ob", near=True)
u      = dm.get_nearest_expression("vx")
v      = dm.get_nearest_expression("vy")

model = model.Model()
model.set_mesh(mesh)
model.set_geometry(S, B,deform=True)
model.set_parameters(pc.IceParameters())
model.calculate_boundaries(mask=M, adot=adot)
model.initialize_variables()

# constraints on optimization for beta :
class Beta_max(Expression):
  def eval(self, values, x):
    if M(x[0], x[1], x[2]) > 0:
      values[0] = 0.0
    else:
      values[0] = 4000

beta_min = interpolate(Constant(0.0), model.Q)
beta_max = interpolate(Beta_max(element = model.Q.ufl_element()), model.Q)

# specifify non-linear solver parameters :
nonlin_solver_params = default_nonlin_solver_params()
nonlin_solver_params['newton_solver']['relaxation_parameter']    = 0.7
nonlin_solver_params['newton_solver']['relative_tolerance']      = 1e-3
nonlin_solver_params['newton_solver']['maximum_iterations']      = 16
nonlin_solver_params['newton_solver']['error_on_nonconvergence'] = False
nonlin_solver_params['newton_solver']['linear_solver']           = 'mumps'
nonlin_solver_params['newton_solver']['preconditioner']          = 'default'
parameters['form_compiler']['quadrature_degree']                 = 2


config = { 'mode'                         : 'steady',
           't_start'                      : None,
           't_end'                        : None,
           'time_step'                    : None,
           'output_path'                  : out_dir,
           'wall_markers'                 : [],
           'periodic_boundary_conditions' : False,
           'log'                          : True,
           'coupled' : 
           { 
             'on'                  : True,
             'inner_tol'           : 0.0,
             'max_iter'            : 1
           },
           'velocity' : 
           { 
             'on'                  : True,
             'newton_params'       : nonlin_solver_params,
             'viscosity_mode'      : 'full',
             'b_linear'            : None,
             'use_T0'              : True,
             'T0'                  : model.T_w - 30.0,
             'A0'                  : 1e-16,
             'beta'                : None,
             'init_beta_from_U_ob' : True,
             'U_ob'                : U_ob,
             'r'                   : 0.0,
             'E'                   : 1.0,
             'approximation'       : 'fo',
             'boundaries'          : None,#'user_defined',
             'u_lat_boundary'      : u,
             'v_lat_boundary'      : v,
             'log'                 : True
           },
           'enthalpy' : 
           { 
             'on'                  : True,
             'use_surface_climate' : False,
             'T_surface'           : T_s,
             'q_geo'               : q_geo,
             'lateral_boundaries'  : None,
             'log'                 : True 
           },
           'free_surface' :
           { 
             'on'                  : False,
             'lump_mass_matrix'    : True,
             'thklim'              : thklim,
             'use_pdd'             : False,
             'observed_smb'        : adot,
           },  
           'age' : 
           { 
             'on'                  : True,
             'use_smb_for_ela'     : True,
             'ela'                 : None,
           },
           'surface_climate' : 
           { 
             'on'                  : False,
             'T_ma'                : None,
             'T_ju'                : None,
             'beta_w'              : None,
             'sigma'               : None,
             'precip'              : None
           },
           'adjoint' :
           { 
             'alpha'               : 0.0,
             'gamma1'              : 1.0,
             'gamma2'              : 100.0,
             'max_fun'             : 20,
             'objective_function'  : 'logarithmic',
             'bounds'              : (beta_min, beta_max),
             'control_variable'    : model.beta,
             'regularization_type' : 'Tikhonov'
           }}

if i != 0:
  #config['velocity']['approximation']       = 'stokes'
  config['velocity']['init_beta_from_U_ob'] = False
  config['velocity']['beta']                = dir_b + str(i-1) + '/beta.xml'
  config['velocity']['T0']                  = dir_b + str(i-1) + '/T.xml'

F = solvers.SteadySolver(model, config)
File(out_dir + 'beta_0.pvd') << model.beta
F.solve()

params = config['velocity']['newton_params']['newton_solver']
params['maximum_iterations']              = 25
config['velocity']['init_beta_from_U_ob'] = False
config['enthalpy']['on']                  = False
config['surface_climate']['on']           = False
config['coupled']['on']                   = False
config['velocity']['use_T0']              = False
config['velocity']['beta']                = model.beta

#if i == 0:
#  params['relaxation_parameter']         = 1.0
#  config['velocity']['viscosity_mode']   = 'linear'
#  config['velocity']['eta']              = model.eta
#  config['adjoint']['alpha']             = 0
#  config['adjoint']['bounds']            = (beta_min, beta_max)
#  config['adjoint']['control_variable']  = model.beta
#
#elif i == 1:
#  b_shf   = project(model.b, model.Q)
#  b_gnd   = b_shf.copy()
#  #b_min  = b_shf.vector().min()/10.0
#  #b_max  = b_shf.vector().max()*10.0
#  b_max   = 1e16
#  b_min   = 0.0
#  model.b = b_shf
#  model.print_min_max(model.b, 'b')
#  File(out_dir + 'b_0.pvd') << model.b
#  
#  params['relaxation_parameter']         = 0.4
#  config['velocity']['viscosity_mode']   = 'b_control'
#  config['velocity']['b_shf']            = b_shf
#  config['velocity']['b_gnd']            = b_gnd
#  config['adjoint']['alpha']             = 0
#  config['adjoint']['bounds']            = (b_min, b_max)
#  config['adjoint']['control_variable']  = b_shf
#
#elif i % 2 == 1:
#  params['relaxation_parameter']         = 0.4
#  config['velocity']['viscosity_mode']   = 'b_control'
#  config['velocity']['b_shf']            = b_shf
#  config['velocity']['b_gnd']            = b_gnd
#  config['adjoint']['alpha']             = 0
#  config['adjoint']['bounds']            = (b_min, b_max)
#  config['adjoint']['control_variable']  = b_shf
#
#else:
#  params['relaxation_parameter']         = 1.0
#  config['velocity']['viscosity_mode']   = 'constant_b'
#  config['velocity']['b']                = File(dir_b + str(i-1) + 'b.xml')
#  config['adjoint']['alpha']             = 0
#  config['adjoint']['bounds']            = (beta_min, beta_max)
#  config['adjoint']['control_variable']  = model.beta

b_shf   = project(model.b, model.Q)
b_gnd   = b_shf.copy()
b_max   = 1e16
b_min   = 0.0
model.b = b_shf
model.print_min_max(model.b, 'b')
File(out_dir + 'b_0.pvd') << model.b
params['relaxation_parameter']         = 0.4
config['velocity']['viscosity_mode']   = 'b_control'
config['velocity']['b_shf']            = b_shf
config['velocity']['b_gnd']            = b_gnd
config['adjoint']['alpha']             = [1e-7, 0]
config['adjoint']['bounds']            = [(beta_min, beta_max), (b_min, b_max)]
config['adjoint']['control_variable']  = [model.beta, b_shf]

A = solvers.AdjointSolver(model, config)
A.set_target_velocity(u=u, v=v)
A.solve()

File(out_dir + 'T.xml')       << model.T
File(out_dir + 'S.xml')       << model.S
File(out_dir + 'B.xml')       << model.B
File(out_dir + 'u.xml')       << project(model.u, model.Q) 
File(out_dir + 'v.xml')       << project(model.v, model.Q) 
File(out_dir + 'w.xml')       << model.w 
File(out_dir + 'beta.xml')    << model.beta
File(out_dir + 'eta.xml')     << project(model.eta, model.Q)
File(out_dir + 'age.xml')     << model.age
File(out_dir + 'b.xml')       << model.b
File(out_dir + 'b.pvd')       << project(model.b,   model.Q)

#XDMFFile(mesh.mpi_comm(), out_dir + 'mesh.xdmf')   << model.mesh
#
## save the state of the model :
#if i !=0: rw = 'a'
#else:     rw = 'w'
#f = HDF5File(mesh.mpi_comm(), out_dir + '3D_5H_stokes.h5', rw)
#f.write(model.mesh,  'mesh')
#f.write(model.beta,  'beta')
#f.write(model.Mb,    'Mb')
#f.write(model.T,     'T')
#f.write(model.S,     'S')
#f.write(model.B,     'B')
#f.write(model.U,     'U')
#f.write(model.eta,   'eta')

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



