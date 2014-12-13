import varglas.solvers            as solvers
import varglas.physical_constants as pc
import varglas.model              as model
from varglas.helper               import default_nonlin_solver_params, \
                                         default_config
from fenics                       import *


alpha = 0.5 * pi / 180 
L     = 100000
nx    = 20
ny    = 20 
nz    = 10

model = model.Model()
model.generate_uniform_mesh(nx, ny, nz, xmin=0, xmax=L, ymin=0, ymax=L) 

surface = Expression('1300 - x[0] * tan(alpha)', alpha=alpha, 
                     element=model.Q.ufl_element())
bed     = Expression(  '1300 - x[0] * tan(alpha) - 1000.0 + 500.0 * ' \
                     + ' sin(2*pi*x[0]/L) * sin(2*pi*x[1]/L)',
                     alpha=alpha, L=L, element=model.Q.ufl_element())
beta    = Expression(  '10 + 10 * sin(2*pi*x[0]/L) * sin(2*pi*x[1]/L)',
                     alpha=alpha, L=L, element=model.Q.ufl_element())

model.set_geometry(surface, bed, deform=True)
model.set_parameters(pc.IceParameters())
model.calculate_boundaries()
model.initialize_variables()

File('output/beta.pvd') << interpolate(beta, model.Q)

# specifify non-linear solver parameters :
nonlin_solver_params = default_nonlin_solver_params()
nonlin_solver_params['newton_solver']['relaxation_parameter']    = 0.9
nonlin_solver_params['newton_solver']['relative_tolerance']      = 1e-14
nonlin_solver_params['newton_solver']['maximum_iterations']      = 25
nonlin_solver_params['newton_solver']['error_on_nonconvergence'] = False
nonlin_solver_params['newton_solver']['linear_solver']           = 'mumps'
nonlin_solver_params['newton_solver']['preconditioner']          = 'default'
parameters['form_compiler']['quadrature_degree']                 = 2

config = default_config()
config['output_path']                  = 'output/'
config['velocity']['newton_params']    = nonlin_solver_params
config['velocity']['approximation']    = 'stokes'
config['velocity']['beta0']            = beta

F = solvers.SteadySolver(model, config)
F.solve()

T = solvers.StokesBalanceSolver(model, config)
T.solve()


