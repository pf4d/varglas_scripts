from varglas.model              import Model
from varglas.solvers            import SteadySolver
from varglas.physical_constants import IceParameters
from varglas.helper             import default_nonlin_solver_params
from fenics                     import set_log_active, File, Expression, pi, \
                                       sin, tan, interpolate, project, \
                                       parameters

set_log_active(True)

alpha = 0.1 * pi / 180
L     = 100000
nx    = 50
ny    = 10 
nz    = 10


model = Model()
model.generate_uniform_mesh(nx, ny, nz, xmin=0, xmax=500000, 
                            ymin=0, ymax=L, 
                            generate_pbcs=True)

Surface = Expression('- x[1] * tan(alpha)', alpha=alpha, 
                     element=model.Q.ufl_element())
Bed     = Expression('- x[1] * tan(alpha) - 2500.0', alpha=alpha, 
                     element=model.Q.ufl_element())
Beta2   = Expression(  '1000 + 1000 * sin(pi*x[0]/L)',
                     alpha=alpha, L=L, element=model.Q.ufl_element())
T_s     = Expression(  '240 +  5 * sin(pi*x[0]/L)',
                     alpha=alpha, L=L, element=model.Q.ufl_element())

model.set_geometry(Surface, Bed, deform=True)
model.set_parameters(IceParameters())
model.calculate_boundaries()
model.initialize_variables()

nonlin_solver_params = default_nonlin_solver_params()
#nonlin_solver_params['newton_solver']['linear_solver']  = 'gmres'
#nonlin_solver_params['newton_solver']['preconditioner'] = 'hypre_amg'
nonlin_solver_params['newton_solver']['relaxation_parameter']    = 1.0
nonlin_solver_params['newton_solver']['error_on_nonconvergence'] = False
nonlin_solver_params['newton_solver']['maximum_iterations']      = 2
nonlin_solver_params['newton_solver']['absolute_tolerance']      = 1.0
nonlin_solver_params['newton_solver']['linear_solver']           = 'mumps'
nonlin_solver_params['newton_solver']['preconditioner']          = 'default'
parameters['form_compiler']['quadrature_degree']        = 2

config = { 'mode'                         : 'steady',
           'output_path'                  : './results_C_BP/',
           'wall_markers'                 : [],
           'periodic_boundary_conditions' : True,
           't_start'                      : None,
           't_end'                        : None,
           'time_step'                    : None,
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
             'beta2'               : Beta2,
             'r'                   : 0.0,
             'E'                   : 1,
             'approximation'       : 'fo',
             'boundaries'          : None
           },
           'enthalpy' : 
           { 
             'on'                  : True,
             'use_surface_climate' : False,
             'T_surface'           : T_s,
             'q_geo'               : 0.00042*60**2*24*365,
             'lateral_boundaries'  : None,
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
             'on'                  : True,
             'use_smb_for_ela'     : False,
             'ela'                 : -10,
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
 
File('./results_C_BP/' + 'B.pvd')     << interpolate(Bed, model.Q)
File('./results_C_BP/' + 'S.pvd')     << project(Surface, model.Q)
File('./results_C_BP/' + 'beta2.pvd') << interpolate(Beta2, model.Q)

F = SteadySolver(model, config)
F.solve()



