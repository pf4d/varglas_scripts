import varglas.solvers            as solvers
import varglas.model              as model
from varglas.mesh.mesh_factory    import MeshFactory
from varglas.data.data_factory    import DataFactory
from varglas.helper               import default_config
from fenics                       import *

set_log_active(False)

in_dir  = 'dump/bed/08/'
out_dir = 'dump/stress/SSA/'

mesh  = Mesh(in_dir + 'submesh.xdmf')
Q     = FunctionSpace(mesh, 'CG', 1)

S     = Function(Q)
B     = Function(Q)
 
File(in_dir + 'S_s.xml') >> S
File(in_dir + 'B_s.xml') >> B

config = default_config()
config['output_path']     = out_dir
config['model_order']     = 'SSA'
config['use_dukowicz']    = False

model = model.Model(config)
model.set_mesh(mesh)
model.set_surface_and_bed(S, B)
model.initialize_variables()

model.init_beta(in_dir + 'beta_s.xml')
model.init_component_Ubar(in_dir + 'ubar_s.xml',
                          in_dir + 'vbar_s.xml')
model.init_etabar(in_dir + 'etabar_s.xml')

T = solvers.StokesBalanceSolver(model, config)
T.solve()

Ubar = project(as_vector([model.ubar, model.vbar]))

model.save_pvd(model.etabar, 'etabar')
model.save_pvd(Ubar,         'Ubar')
model.save_xml(model.tau_id, 'tau_id')
model.save_xml(model.tau_jd, 'tau_jd')
model.save_xml(model.tau_ib, 'tau_ib')
model.save_xml(model.tau_jb, 'tau_jb')
model.save_xml(model.tau_ii, 'tau_ii')
model.save_xml(model.tau_ij, 'tau_ij')
model.save_xml(model.tau_ji, 'tau_ji')
model.save_xml(model.tau_jj, 'tau_jj')



