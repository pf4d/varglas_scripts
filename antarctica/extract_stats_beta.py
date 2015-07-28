import varglas.model           as model
from fenics                    import *
from pylab                     import *
from varglas.helper            import default_config
from varglas.data.data_factory import DataFactory

out_dir  = 'dump/linear_model/'
in_dir   = 'dump/linear_model/'
var_dir  = 'dump/vars_high/'

#mesh   = Mesh(var_dir + 'mesh.xdmf')
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
config['model_order']     = 'SSA'

model = model.Model(config)
model.set_mesh(mesh)
model.set_surface_and_bed(S, B)
model.set_subdomains(ff, cf, ff_acc)
model.initialize_variables()

model.init_beta(in_dir + 'beta.xml')
model.init_U(in_dir + 'u.xml',
             in_dir + 'v.xml',
             in_dir + 'w.xml')

us = model.vert_extrude(model.u, d='down')
vs = model.vert_extrude(model.v, d='down')
ws = model.vert_extrude(model.w, d='down')

#submesh = model.get_bed_mesh()
submesh = Mesh('dump/bed/07/submesh.xdmf')

Q_b     = FunctionSpace(submesh, 'CG', 1)
#beta_s  = Function(Q_b)
us_s    = Function(Q_b)
vs_s    = Function(Q_b)
ws_s    = Function(Q_b)

lg      = LagrangeInterpolator()

#lg.interpolate(beta_s, model.beta)
lg.interpolate(us_s,   us)
lg.interpolate(vs_s,   vs)
lg.interpolate(ws_s,   ws)

#model.save_xml(beta_s, 'beta_s')
model.save_xml(us_s,   'us_s')
model.save_xml(vs_s,   'vs_s')
model.save_xml(ws_s,   'ws_s')

#measures   = DataFactory.get_ant_measures(res=900)
#dm         = DataInput(measures, gen_space=False)
#
#plotIce(dm, beta,  name='beta', direc=out_dir, cmap='gist_yarg',
#        title=r'$\beta$', scale='lin', umin=1.0, umax=200,
#        numLvls=12, tp=False, tpAlpha=0.5, extend='max', show=False)
