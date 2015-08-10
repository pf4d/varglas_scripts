from fenics                    import *
from varglas.solvers           import HybridTransientSolver
from varglas.data.data_factory import DataFactory
from varglas.mesh.mesh_factory import MeshFactory
from varglas.io                import DataInput, print_min_max
from varglas.helper            import default_config
from varglas.model             import Model

set_log_active(False)

thklim = 1.0

# collect the raw data :
searise = DataFactory.get_searise(thklim = thklim)
bamber  = DataFactory.get_bamber(thklim = thklim)

# load a mesh :
mesh  = MeshFactory.get_greenland_2D_5H()

# create data objects to use with varglas :
dsr   = DataInput(searise, mesh=mesh)
dbm   = DataInput(bamber,  mesh=mesh)

dbm.data['S'] = dbm.data['B'] + thklim

B     = dbm.get_expression('B',    near=True)
S     = dbm.get_expression('S',    near=True)
adot  = dsr.get_expression('adot', near=True)
T_s   = dsr.get_interpolation('T', near=True)

out_dir = 'dump/transient_H/'

parameters['form_compiler']['quadrature_degree'] = 2
parameters['form_compiler']['precision']         = 30
#parameters['form_compiler']['optimize']          = True
parameters['form_compiler']['cpp_optimize']      = True
parameters['form_compiler']['representation']    = 'quadrature'

config = default_config()
config['log']                          = True
config['log_history']                  = False
config['mode']                         = 'transient'
config['model_order']                  = 'L1L2'
config['output_path']                  = out_dir
config['t_start']                      = 0.0
config['t_end']                        = 35000.0
config['time_step']                    = 10.0
config['periodic_boundary_conditions'] = False
config['velocity']['poly_degree']      = 2
config['enthalpy']['on']               = True
config['enthalpy']['N_T']              = 8
config['free_surface']['on']           = True
config['free_surface']['thklim']       = thklim
config['velocity']['transient_beta']   = 'eismint_H'

model = Model(config)
model.set_mesh(mesh)

model.set_geometry(S, B, deform=False)
model.initialize_variables()

model.init_adot(adot)
model.init_beta(sqrt(1e9))
model.init_T_surface(T_s)
model.init_H(thklim)
model.init_H_bounds(thklim, 1e4)
model.init_q_geo(model.ghf)

model.eps_reg = 1e-10

T = HybridTransientSolver(model, config)
T.solve()

File(out_dir + 'Ts.xml')      << model.Ts
File(out_dir + 'Tb.xml')      << model.Tb
File(out_dir + 'Mb.xml')      << model.Mb
File(out_dir + 'H.xml')       << model.H
File(out_dir + 'S.xml')       << model.S
File(out_dir + 'B.xml')       << model.B
File(out_dir + 'u_s.xml')     << model.u_s
File(out_dir + 'v_s.xml')     << model.v_s
File(out_dir + 'w_s.xml')     << model.w_s
File(out_dir + 'u_b.xml')     << model.u_b
File(out_dir + 'v_b.xml')     << model.v_b
File(out_dir + 'w_b.xml')     << model.w_b
File(out_dir + 'Ubar.xml')    << model.Ubar
File(out_dir + 'beta.xml')    << model.beta



