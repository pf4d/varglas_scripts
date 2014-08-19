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

thklim = 200.0

# collect the raw data :
searise  = DataFactory.get_searise(thklim = thklim)
bamber   = DataFactory.get_bamber(thklim = thklim)
fm_qgeo  = DataFactory.get_gre_qgeo_fox_maule()
rignot   = DataFactory.get_gre_rignot_updated()

# define the mesh :
mesh = MeshFactory.get_greenland_detailed()

# create data objects to use with varglas :
dsr     = DataInput(searise,  mesh=mesh)
dbm     = DataInput(bamber,   mesh=mesh)
dfm     = DataInput(fm_qgeo,  mesh=mesh)
drg     = DataInput(rignot,   mesh=mesh)

# change the projection of the measures data to fit with other data :
drg.change_projection(dsr)

# get the expressions used by varglas :
H     = dbm.get_nearest_expression('H')
S     = dbm.get_nearest_expression('S')
B     = dbm.get_nearest_expression('B')
T_s   = dsr.get_nearest_expression('T')
q_geo = dfm.get_nearest_expression('q_geo')
adot  = dsr.get_nearest_expression('adot')
U_ob  = drg.get_nearest_expression('U_ob')
u     = drg.get_nearest_expression("vx")
v     = drg.get_nearest_expression("vy")

# inspect the data values :
#do    = DataOutput('results_pre/')
#do.write_one_file('vmag',           drg.get_projection('U_ob'))
#do.write_one_file('h',              dbm.get_projection('H'))
#do.write_one_file('Ubmag_measures', dbv.get_projection('Ubmag_measures'))
#do.write_one_file('sr_qgeo',        dsr.get_projection('q_geo'))
#exit(0)

model = model.Model()
model.set_mesh(mesh)
model.set_geometry(S, B, deform=True)
model.set_parameters(pc.IceParameters())
model.calculate_boundaries(adot=adot)
model.initialize_variables()

# specifify non-linear solver parameters :
nonlin_solver_params = default_nonlin_solver_params()
nonlin_solver_params['newton_solver']['relaxation_parameter']    = 0.7
nonlin_solver_params['newton_solver']['relative_tolerance']      = 1e-3
nonlin_solver_params['newton_solver']['absolute_tolerance']      = 1e2
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
             'max_iter'            : 2
           },                      
           'velocity' :            
           {                       
             'on'                  : True,
             'newton_params'       : nonlin_solver_params,
             'init_beta_from_U_ob' : True,
             'U_ob'                : U_ob,
             'viscosity_mode'      : 'full',
             'eta'                 : None,
             'use_T0'              : True,
             'use_beta0'           : False,
             'T0'                  : 268.0,
             'A0'                  : None,
             'beta0'               : 0.5,
             'r'                   : 1.0,
             'E'                   : 1.0,
             'approximation'       : 'fo',
             'boundaries'          : None,
             'log'                 : True
           },
           'enthalpy' : 
           { 
             'on'                  : True,
             'use_surface_climate' : False,
             'T_surface'           : T_s,
             'q_geo'               : q_geo,
             'lateral_boundaries'  : None,
             'log'                 : True,
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
             'alpha'               : 0.0,#[H**2],
             'max_fun'             : 20,
             'objective_function'  : 'logarithmic',
             'bounds'              : (0,100),
             'control_variable'    : model.beta,
             'regularization_type' : 'Tikhonov'
           }}

if i != 0: 
  config['velocity']['approximation']       = 'stokes'
  config['velocity']['init_beta_from_U_ob'] = False
  config['velocity']['beta']                = dir_b + str(i-1) + '/beta.xml'
  config['velocity']['T0']                  = dir_b + str(i-1) + '/T.xml'

F = solvers.SteadySolver(model,config)
File(out_dir + 'beta_0.pvd') << model.beta
F.solve()

params = config['velocity']['newton_params']['newton_solver']
params['relaxation_parameter']         = 1.0
params['relative_tolerance']           = 1e-3
params['absolute_tolerance']           = 0.0
params['maximum_iterations']           = 12
config['velocity']['viscosity_mode']   = 'linear'
config['velocity']['eta']              = model.eta
config['enthalpy']['on']               = False
config['surface_climate']['on']        = False
config['coupled']['on']                = False
config['velocity']['use_T0']           = False
config['velocity']['use_beta0']        = False

A = solvers.AdjointSolver(model,config)
A.set_target_velocity(u=u, v=v)
A.solve()
    
File(out_dir + 'S.xml')      << model.S
File(out_dir + 'B.xml')      << model.B
File(out_dir + 'u.xml')      << project(model.u, model.Q)
File(out_dir + 'v.xml')      << project(model.v, model.Q)
File(out_dir + 'w.xml')      << project(model.w, model.Q)
File(out_dir + 'P.xml')      << project(model.P, model.Q)
File(out_dir + 'T.xml')      << model.T
File(out_dir + 'beta.xml')   << model.beta
File(out_dir + 'eta.xml')    << project(model.eta, model.Q)

#XDMFFile(mesh.mpi_comm(), out_dir + 'mesh.xdmf')   << model.mesh
#
## save the state of the model :
#if i !=0: rw = 'a'
#else:     rw = 'w'
#f = HDF5File(out_dir + '3D_5H_stokes.h5', rw)
#f.write(model.mesh, 'mesh')
#f.write(model.beta, 'beta')
#f.write(model.Mb,   'Mb')
#f.write(model.T,    'T')
#f.write(model.S,    'S')
#f.write(model.B,    'B')
#f.write(model.U,    'U')
#f.write(model.eta,  'eta')

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



