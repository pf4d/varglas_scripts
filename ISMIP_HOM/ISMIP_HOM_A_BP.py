from varglas.model              import Model
from varglas.solvers            import SteadySolver
from varglas.physical_constants import IceParameters
from varglas.helper             import default_nonlin_solver_params
from fenics                     import set_log_active, File, Expression, pi, \
                                       sin, tan, interpolate

set_log_active(True)

nonlin_solver_params = default_nonlin_solver_params()
nonlin_solver_params['newton_solver']['linear_solver']   = 'gmres'
nonlin_solver_params['newton_solver']['preconditioner']  = 'hypre_amg'

config = { 'mode'                         : 'steady',
           'output_path'                  : './results_A_BP/',
           'wall_markers'                 : [],
           'periodic_boundary_conditions' : True,
           't_start'                      : None,
           't_end'                        : None,
           'time_step'                    : None,
           'log'                          : True,
           'coupled' : 
           { 
             'on'                  : False,
             'inner_tol'           : 0.0,
             'max_iter'            : 1
           },
           'velocity' : 
           { 
             'on'                  : True,
             'newton_params'       : nonlin_solver_params,
             'viscosity_mode'      : 'isothermal',
             'b_linear'            : None,
             'use_T0'              : False,
             'T0'                  : None,
             'A0'                  : 1e-16,
             'beta2'               : 1e3,
             'r'                   : 1,
             'E'                   : 1,
             'approximation'       : 'fo',
             'boundaries'          : None
           },
           'enthalpy' : 
           { 
             'on'                  : False,
             'use_surface_climate' : False,
             'T_surface'           : None,
             
           },
           'free_surface' :
           { 
             'on'                  : False,
             'thklim'              : None,
             'use_pdd'             : False,
             'observed_smb'        : None,
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
             'alpha'               : None,
             'beta'                : None,
             'max_fun'             : None,
             'objective_function'  : 'logarithmic',
             'animate'             : False
           }}


alpha = 0.5 * pi / 180 
L     = 40000
nx    = 50
ny    = 50 
nz    = 10


model = Model()
model.generate_uniform_mesh(nx, ny, nz, xmin=0, xmax=L, ymin=0, ymax=L, 
                            generate_pbcs=True)

Surface = Expression('- x[0] * tan(alpha)', alpha=alpha, 
                     element=model.Q.ufl_element())
Bed     = Expression(  '- x[0] * tan(alpha) - 1000.0 + 500.0 * ' \
                     + ' sin(2*pi*x[0]/L) * sin(2*pi*x[1]/L)',
                     alpha=alpha, L=L, element=model.Q.ufl_element())

File('./results_A_BP/' + 'B.pvd') << interpolate(Bed, model.Q)

model.set_geometry(Surface, Bed, deform=True)
model.set_parameters(IceParameters())
model.calculate_boundaries()
model.initialize_variables()
 
F = SteadySolver(model, config)
F.solve()



