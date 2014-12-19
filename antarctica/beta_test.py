import varglas.solvers            as solvers
import varglas.physical_constants as pc
import varglas.model              as model
from varglas.mesh.mesh_factory    import MeshFactory
from varglas.data.data_factory    import DataFactory
from varglas.helper               import default_nonlin_solver_params, \
                                         default_config
from varglas.utilities            import DataInput
from fenics                       import *
from time                         import time
from termcolor                    import colored, cprint

t0      = time()
out_dir = 'output/beta_test/'
thklim  = 1.0

measures  = DataFactory.get_ant_measures(res=900)
bedmap1   = DataFactory.get_bedmap1(thklim=thklim)
bedmap2   = DataFactory.get_bedmap2(thklim=thklim)

mesh = MeshFactory.get_antarctica_3D_10k()

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

# specifify non-linear solver parameters :
nonlin_solver_params = default_nonlin_solver_params()
nonlin_solver_params['newton_solver']['relaxation_parameter']    = 0.7
nonlin_solver_params['newton_solver']['relative_tolerance']      = 1e-3
nonlin_solver_params['newton_solver']['maximum_iterations']      = 16
nonlin_solver_params['newton_solver']['error_on_nonconvergence'] = False
nonlin_solver_params['newton_solver']['linear_solver']           = 'mumps'
nonlin_solver_params['newton_solver']['preconditioner']          = 'default'
parameters['form_compiler']['quadrature_degree']                 = 2

config = default_config()
config['output_path']                     = out_dir
config['coupled']['on']                   = True
config['coupled']['max_iter']             = 10
config['velocity']['newton_params']       = nonlin_solver_params
config['velocity']['approximation']       = 'fo'#'stokes'
config['velocity']['viscosity_mode']      = 'full'
config['velocity']['use_T0']              = True
config['velocity']['use_beta0']           = False
config['velocity']['T0']                  = model.T_w - 30.0
config['velocity']['use_stats_beta']      = True
config['velocity']['U_ob']                = U_ob
config['enthalpy']['on']                  = True
config['enthalpy']['T_surface']           = T_s
config['enthalpy']['q_geo']               = model.q_geo

F = solvers.SteadySolver(model, config)
F.solve()

File(out_dir + 'beta.pvd') << project(model.beta, model.Q)

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



