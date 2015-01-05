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
out_dir = 'test_all_procs/'     # directory to save
thklim  = 1.0

measures  = DataFactory.get_ant_measures(res=900)
bedmap1   = DataFactory.get_bedmap1(thklim=thklim)
bedmap2   = DataFactory.get_bedmap2(thklim=thklim)

#mesh = MeshFactory.get_antarctica_3D_gradS_detailed()
#mesh = MeshFactory.get_antarctica_3D_10k()
mesh = Mesh('meshes/ant_high.xml')

dm  = DataInput(measures, mesh=mesh)
db1 = DataInput(bedmap1,  mesh=mesh)
db2 = DataInput(bedmap2,  mesh=mesh)
    
db2.data['B'] = db2.data['S'] - db2.data['H']
db2.set_data_val('H', 32767, thklim)
db2.data['S'] = db2.data['B'] + db2.data['H']

S    = db2.get_nearest_expression("S")
B    = db2.get_nearest_expression("B")
M    = db2.get_nearest_expression("mask")
T_s  = db1.get_nearest_expression("temp")
U_ob = dm.get_projection("U_ob", near=True)

model = model.Model()
model.set_mesh(mesh)
model.set_geometry(S, B, deform=True)
model.set_parameters(pc.IceParameters())
model.calculate_boundaries(mask=M)
model.initialize_variables()

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
params['newton_solver']['maximum_iterations']       = 100
params['newton_solver']['error_on_nonconvergence']  = False
params['newton_solver']['linear_solver']            = 'cg'
params['newton_solver']['preconditioner']           = 'hypre_amg'
#params['newton_solver']['krylov_solver']['monitor_convergence']  = True
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
config['velocity']['T0']                  = model.T_w - 30.0
config['velocity']['init_beta_from_U_ob'] = True
config['velocity']['U_ob']                = U_ob
config['enthalpy']['on']                  = True
config['enthalpy']['T_surface']           = T_s
config['enthalpy']['q_geo']               = model.q_geo

F = solvers.SteadySolver(model, config)

t0 = time()
F.solve()
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



