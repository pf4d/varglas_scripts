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

set_log_active(True)
#set_log_level(PROGRESS)

thklim = 200.0

measures  = DataFactory.get_ant_measures(res=900)
bedmap1   = DataFactory.get_bedmap1(thklim=thklim)
bedmap2   = DataFactory.get_bedmap2(thklim=thklim)

#mesh = MeshFactory.get_antarctica_coarse()
mesh = MeshFactory.get_antarctica_3D_50H()

dm  = DataInput(None, measures, mesh=mesh)
db1 = DataInput(None, bedmap1,  mesh=mesh)
db2 = DataInput(None, bedmap2,  mesh=mesh)

db2.set_data_val('H',   32767, thklim)
db2.set_data_val('S',   32767, 0.0)
dm.set_data_min('U_ob', 0.0,   0.0)

db2.data['B'] = db2.data['S'] - db2.data['H']

H     = db2.get_spline_expression("H")
S     = db2.get_spline_expression("S")
B     = db2.get_spline_expression("B")
M     = db2.get_nearest_expression("mask")
T_s   = db1.get_spline_expression("srfTemp")
SMB   = db1.get_spline_expression("q_geo")
adot  = db1.get_spline_expression("adot")
U_ob  = dm.get_spline_expression("U_ob")

model = model.Model()
model.set_mesh(mesh)
model.set_geometry(S, B, mask=M, deform=True)
model.set_parameters(pc.IceParameters())
model.initialize_variables()

File(out_dir + 'ff.pvd') << model.ff

# specifify non-linear solver parameters :
nonlin_solver_params = default_nonlin_solver_params()
nonlin_solver_params['newton_solver']['relaxation_parameter']    = 0.7
nonlin_solver_params['newton_solver']['relative_tolerance']      = 1e-3
nonlin_solver_params['newton_solver']['maximum_iterations']      = 25
nonlin_solver_params['newton_solver']['error_on_nonconvergence'] = False
#nonlin_solver_params['newton_solver']['linear_solver']           = 'gmres'
#nonlin_solver_params['newton_solver']['preconditioner']          = 'hypre_amg'
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
             'T0'                  : 268.0,
             'A0'                  : 1e-16,
             'beta2'               : 0.5,
             'r'                   : 1.0,
             'E'                   : 1.0,
             'approximation'       : 'fo',
             'boundaries'          : None
           },
           'enthalpy' : 
           { 
             'on'                  : True,
             'use_surface_climate' : False,
             'T_surface'           : T_s,
             'q_geo'               : SMB,# 0.042*60*60*24*365,
             'lateral_boundaries'  : None
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
             'on'                  : False,
             'use_smb_for_ela'     : False,
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
             'alpha'               : [H**2],
             'beta'                : 0.0,
             'max_fun'             : 20,
             'objective_function'  : 'logarithmic',
             'bounds'              : [(0,20)],
             'control_variable'    : None,
             'regularization_type' : 'Tikhonov'
           }}

F = solvers.SteadySolver(model,config)
if i != 0: 
  File(dir_b + str(i-1) + '/beta2.xml') >> model.beta2
  config['velocity']['approximation'] = 'stokes'
F.solve()

params = config['velocity']['newton_params']['newton_solver']
params['relaxation_parameter']         = 1.0
config['velocity']['viscosity_mode']   = 'linear'
config['velocity']['b_linear']         = model.eta
config['enthalpy']['on']               = False
config['surface_climate']['on']        = False
config['coupled']['on']                = False
config['velocity']['use_T0']           = False
config['adjoint']['control_variable']  = [model.beta2]

A = solvers.AdjointSolver(model,config)
A.set_target_velocity(U = U_ob)
if i != 0: File(dir_b + str(i-1) + '/beta2.xml') >> model.beta2
A.solve()
    
File(out_dir + 'S.xml')       << model.S
File(out_dir + 'B.xml')       << model.B
File(out_dir + 'u.xml')       << project(model.u, model.Q)
File(out_dir + 'v.xml')       << project(model.v, model.Q)
File(out_dir + 'w.xml')       << project(model.w, model.Q)
File(out_dir + 'beta2.xml')   << model.beta2
File(out_dir + 'eta.xml')     << project(model.eta, model.Q)

#XDMFFile(mesh.mpi_comm(), out_dir + 'mesh.xdmf')   << model.mesh
#
## save the state of the model :
#if i !=0: rw = 'a'
#else:     rw = 'w'
#f = HDF5File(out_dir + '3D_5H_stokes.h5', rw)
#f.write(model.mesh,  'mesh')
#f.write(model.beta2, 'beta2')
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
if MPI.rank(mpi_comm_world()) == 0:
  s    = "Total time to compute: %02d:%02d:%02d" % (h,m,s)
  text = colored(s, 'red', attrs=['bold'])
  print text



