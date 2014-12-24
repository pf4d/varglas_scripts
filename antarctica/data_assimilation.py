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
#mesh  = MeshFactory.get_antarctica_3D_10k()

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
T_s    = db1.get_nearest_expression("temp")
q_geo  = db1.get_nearest_expression("ghfsr")
adot   = db1.get_nearest_expression("acca")
U_ob   = dm.get_projection("U_ob", near=True)
u      = dm.get_nearest_expression("vx")
v      = dm.get_nearest_expression("vy")

model = model.Model()
model.set_mesh(mesh)
model.set_geometry(S, B,deform=True)
model.set_parameters(pc.IceParameters())
model.calculate_boundaries(mask=M, adot=adot)
model.initialize_variables()

File(out_dir + 'U_ob.pvd') << U_ob

# constraints on optimization for beta :
class Beta_max(Expression):
  def eval(self, values, x):
    if M(x[0], x[1], x[2]) > 0:
      values[0] = 0.0
    else:
      values[0] = 4000

# constraints on optimization for b :
class B_max(Expression):
  def eval(self, values, x):
    if M(x[0], x[1], x[2]) > 0:
      values[0] = 1e10
    else:
      values[0] = 0.0

beta_min = interpolate(Constant(0.0), model.Q)
beta_max = interpolate(Beta_max(element = model.Q.ufl_element()), model.Q)

b_min    = interpolate(Constant(0.0), model.Q)
b_max    = interpolate(B_max(element = model.Q.ufl_element()), model.Q)

# specifify non-linear solver parameters :
params = default_nonlin_solver_params()
params['nonlinear_solver']                          = 'snes'
params['snes_solver']['method']                     = 'newtonls'
params['snes_solver']['line_search']                = 'bt'
params['snes_solver']['error_on_nonconvergence']    = False
params['snes_solver']['absolute_tolerance']         = 1.0
params['snes_solver']['relative_tolerance']         = 1e-3
params['snes_solver']['maximum_iterations']         = 20
params['snes_solver']['linear_solver']              = 'gmres'
params['snes_solver']['preconditioner']             = 'hypre_amg'
#params['nonlinear_solver']                          = 'newton'
#params['newton_solver']['relaxation_parameter']     = 1.0
#params['newton_solver']['relative_tolerance']       = 1e-3
#params['newton_solver']['maximum_iterations']       = 16
#params['newton_solver']['error_on_nonconvergence']  = False
#params['newton_solver']['linear_solver']            = 'bicgstab'
#params['newton_solver']['preconditioner']           = 'hypre_amg'
params['newton_solver']['krylov_solver']['monitor_convergence']  = True
parameters['form_compiler']['quadrature_degree']    = 2
#parameters['krylov_solver']['monitor_convergence']  = True
#parameters['lu_solver']['verbose']                  = True


config = default_config()
config['output_path']                     = out_dir
config['coupled']['on']                   = True
config['coupled']['max_iter']             = 5
config['velocity']['newton_params']       = params
config['velocity']['approximation']       = 'fo'#'stokes'
config['velocity']['viscosity_mode']      = 'full'
config['velocity']['use_T0']              = True
config['velocity']['use_U0']              = False
config['velocity']['use_beta0']           = False
config['velocity']['T0']                  = model.T_w - 30.0
config['velocity']['init_beta_from_U_ob'] = True
config['velocity']['init_b_from_U_ob']    = False
config['velocity']['U_ob']                = U_ob
config['velocity']['boundaries']          = None#'user_defined',
config['velocity']['u_lat_boundary']      = u
config['velocity']['v_lat_boundary']      = v
config['enthalpy']['on']                  = True
config['enthalpy']['T_surface']           = T_s
config['enthalpy']['q_geo']               = model.q_geo
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
tf = time()

#params['newton_solver']['maximum_iterations'] = 25
#config['velocity']['init_beta_from_U_ob']     = False
#config['velocity']['use_T0']                  = False
#config['velocity']['use_U0']                  = False
#config['velocity']['use_beta0']               = False
#config['velocity']['use_b_shf0']              = False
#config['enthalpy']['on']                      = False
#config['coupled']['on']                       = False
#
#if i % 2 == 0:
#  params['newton_solver']['relaxation_parameter']  = 1.0
#  config['velocity']['viscosity_mode']             = 'linear'
#  config['velocity']['eta_shf']                    = model.eta_shf
#  config['velocity']['eta_gnd']                    = model.eta_gnd
#  config['adjoint']['surface_integral']            = 'grounded'
#  config['adjoint']['alpha']                       = 0
#  config['adjoint']['bounds']                      = (beta_min, beta_max)
#  config['adjoint']['control_variable']            = model.beta
#
#else:
#  if i > 2:
#    config['velocity']['use_b_shf0']       = True
#    config['velocity']['b_shf']            = dir_b + str(i-2) + '/b_shf.xml'
#  params['relaxation_parameter']         = 0.6
#  b = project(model.b_shf)
#  model.print_min_max(b, 'b')
#  config['velocity']['viscosity_mode']   = 'b_control'
#  config['velocity']['b_shf']            = b
#  config['velocity']['b_gnd']            = b.copy()
#  b_min, b_max = (0.0, 1e10)
#  config['adjoint']['surface_integral']  = 'shelves'
#  config['adjoint']['alpha']             = 0
#  config['adjoint']['bounds']            = (b_min, b_max)
#  config['adjoint']['control_variable']  = b
#  #params['relaxation_parameter']         = 0.6
#  #E = model.E
#  #model.print_min_max(E, 'E')
#  #config['velocity']['viscosity_mode']   = 'E_control'
#  #config['velocity']['E_shf']            = E
#  #config['velocity']['E_gnd']            = E.copy()
#  #E_min, E_max = (1e-16, 100.0)
#  #config['adjoint']['surface_integral']  = 'shelves'
#  #config['adjoint']['alpha']             = 0
#  #config['adjoint']['bounds']            = (E_min, E_max)
#  #config['adjoint']['control_variable']  = E
#
#A = solvers.AdjointSolver(model, config)
#A.set_target_velocity(u=u, v=v)
##uf = dir_b + str(i-1) + '/u.xml'
##vf = dir_b + str(i-1) + '/v.xml'
##wf = dir_b + str(i-1) + '/w.xml'
##A.set_velocity(uf, vf, wf)
#A.solve()
#
#eta   = project(model.eta, model.Q)
#b_shf = project(model.b_shf, model.Q)
#b_gnd = project(model.b_gnd, model.Q)
#
#File(out_dir + 'T.xml')       << model.T
#File(out_dir + 'S.xml')       << model.S
#File(out_dir + 'B.xml')       << model.B
#File(out_dir + 'u.xml')       << model.u 
#File(out_dir + 'v.xml')       << model.v 
#File(out_dir + 'w.xml')       << model.w 
#File(out_dir + 'beta.xml')    << model.beta
#File(out_dir + 'Mb.xml')      << model.Mb
#File(out_dir + 'eta.xml')     << eta
#File(out_dir + 'b_shf.xml')   << b_shf
#File(out_dir + 'b_shf.pvd')   << b_shf
#File(out_dir + 'b_gnd.xml')   << b_gnd
#File(out_dir + 'E_shf.xml')   << model.E_shf
#File(out_dir + 'E_shf.pvd')   << model.E_shf

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
#f.write(model.b_gnd, 'b_gnd')
#f.write(model.b_shf, 'b_shf')

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



