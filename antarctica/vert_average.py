import varglas.model              as model
from varglas.mesh.mesh_factory    import MeshFactory
from fenics                       import *
from varglas.helper               import default_config

# set the in/out directory :
out_dir = 'dump/high/vert_average/'
in_dir  = 'dump/high/07/'
var_dir = 'dump/vars_high/'

mesh   = Mesh(var_dir + 'mesh.xdmf')
Q      = FunctionSpace(mesh, 'CG', 1)
ff     = MeshFunction('size_t', mesh)
cf     = MeshFunction('size_t', mesh)
ff_acc = MeshFunction('size_t', mesh)
S      = Function(Q)
B      = Function(Q)
mask   = Function(Q)

f = HDF5File(mesh.mpi_comm(), var_dir + 'vars.h5', 'r')

f.read(ff,     'ff')
f.read(cf,     'cf')
f.read(ff_acc, 'ff_acc')
f.read(S,      'S')
f.read(B,      'B')
f.read(mask,   'mask')

config = default_config()
config['output_path']    = out_dir
config['model_order']    = 'SSA'

model = model.Model(config)
model.set_mesh(mesh)
model.set_surface_and_bed(S, B)
model.set_subdomains(ff, cf, ff_acc)
model.initialize_variables()

model.init_mask(mask)
model.init_U(in_dir + 'u.xml',
             in_dir + 'v.xml',
             in_dir + 'w.xml')
model.init_T(in_dir + 'T.xml')
model.init_W(in_dir + 'W.xml')
model.init_E_shf(in_dir + 'E_shf.xml')
model.init_E_gnd(in_dir + 'E_gnd.xml')

E_shf_v = model.E_shf.vector().array()
E_shf_v[E_shf_v < 1e-2] = 1e-2
model.assign_variable(model.E_shf, E_shf_v)

model.calc_eta()

us = model.vert_extrude(model.u, d='down')
vs = model.vert_extrude(model.v, d='down')
ws = model.vert_extrude(model.w, d='down')

etabar = model.calc_vert_average(model.eta)
ubar   = model.calc_vert_average(model.u)
vbar   = model.calc_vert_average(model.v)
wbar   = model.calc_vert_average(model.w)

model.save_pvd(etabar, 'etabar')
model.save_pvd(ubar,   'ubar')
model.save_pvd(vbar,   'vbar')
model.save_pvd(wbar,   'wbar')

model.save_xml(etabar, 'etabar')
model.save_xml(ubar,   'ubar')
model.save_xml(vbar,   'vbar')
model.save_xml(wbar,   'wbar')
model.save_xml(us,     'us')
model.save_xml(vs,     'vs')
model.save_xml(ws,     'ws')



