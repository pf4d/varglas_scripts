import varglas.solvers            as solvers
import varglas.physical_constants as pc
import varglas.model              as model
from varglas.helper               import default_nonlin_solver_params, \
                                         default_config
from fenics                       import *


t    = 1000000.0 / 2
S0   = 200.0
bm   = 1000.0
xmin = -t
xmax = t
ymin = -t
ymax = t

mesh = Mesh('../meshes/unit_cyl_mesh.xml')

# width and origin of the domain for deforming x coord :
width_x  = xmax - xmin
offset_x = xmin

# width and origin of the domain for deforming y coord :
width_y  = ymax - ymin
offset_y = ymin

# Deform the square to the defined geometry :
for x in mesh.coordinates():
  # transform x :
  x[0]  = x[0]  * width_x

  # transform y :
  x[1]  = x[1]  * width_y

def gauss(x, y, sigx, sigy):
  return exp(-((x/(2*sigx))**2 + (y/(2*sigy))**2))

class Surface(Expression):
  def eval(self,values,x):
    values[0] = S0 + 1200*gauss(x[0], x[1], t/2, t/2)
S = Surface(element = Q.ufl_element())

class Bed(Expression):
  def eval(self,values,x):
    values[0] = S0 - 500.0 \
                - 1000.0 * gauss(x[0], x[1], t/2, t/2)
B = Bed(element = Q.ufl_element())

class Beta(Expression):
  def eval(self, values, x):
    values[0] = bm * gauss(x[0], x[1], t/2, t/2)
fric = Beta(element = Q.ufl_element())

model = model.Model()
model.set_mesh(mesh)
model.set_geometry(S, B, deform=True)
model.set_parameters(pc.IceParameters())
model.calculate_boundaries()
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
config['output_path']                     = 'output'
config['velocity']['newton_params']       = nonlin_solver_params

F = solvers.SteadySolver(model, config)
F.solve()

T = solvers.StokesBalanceSolver(model, config)
T.solve()


