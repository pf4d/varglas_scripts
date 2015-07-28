import varglas.model           as model
from fenics                    import *
from pylab                     import *
from varglas.helper            import default_config
from varglas.data.data_factory import DataFactory

out_dir  = 'dump/linear_model/'
in_dir   = 'dump/linear_model/'
var_dir  = 'dump/vars_high/'

mesh   = Mesh('dump/meshes/ant_mesh_high.xml')
Q      = FunctionSpace(mesh, 'CG', 1)
ff     = MeshFunction('size_t', mesh)
cf     = MeshFunction('size_t', mesh)
ff_acc = MeshFunction('size_t', mesh)

S      = Function(Q)
B      = Function(Q)

f = HDF5File(mesh.mpi_comm(), var_dir + 'vars.h5', 'r')

f.read(S,     'S')
f.read(B,     'B')
f.read(ff,    'ff')
f.read(cf,    'cf')
f.read(ff_acc,'ff_acc')

config = default_config()
config['output_path']     = out_dir
config['model_order']     = 'BP'

model = model.Model(config)
model.set_mesh(mesh)
model.set_surface_and_bed(S, B)
model.set_subdomains(ff, cf, ff_acc)
model.initialize_variables()

model.init_beta(in_dir + 'beta.xml')
model.init_U(in_dir + 'u.xml',
             in_dir + 'v.xml',
             in_dir + 'w.xml')

submesh = model.get_bed_mesh()
#submesh = Mesh('dump/bed/07/submesh.xdmf')

vertex_indices = submesh.data().array("parent_vertex_indices", 0)

Q_b     = FunctionSpace(submesh, 'CG', 1)
beta_s  = Function(Q_b)

File(out_dir + 'beta_s.xml') >> beta_s

model.beta[vertex_indices] = beta_s.vector().array()

model.save_pvd(model.beta, 'beta_s_test')

