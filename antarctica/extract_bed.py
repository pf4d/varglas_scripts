import varglas.model as model
from fenics          import *
from pylab           import *
from varglas.helper  import default_config

out_dir  = 'dump/bed/07/'
in_dir   = 'dump/high/07/'
bv_dir   = 'dump/high/balance_velocity/'
var_dir  = 'dump/vars_high/'
avg_dir  = 'dump/high/vert_average/'
str_dir  = 'dump/stress/BP/'

#mesh   = Mesh(var_dir + 'mesh.xdmf')
mesh = Mesh('dump/meshes/ant_mesh_high.xml')
Q      = FunctionSpace(mesh, 'CG', 1)
ff     = MeshFunction('size_t', mesh)
cf     = MeshFunction('size_t', mesh)
ff_acc = MeshFunction('size_t', mesh)

S      = Function(Q)
B      = Function(Q)
T_s    = Function(Q)
U_ob   = Function(Q)
adot   = Function(Q)
q_geo  = Function(Q)
mask   = Function(Q)

tau_id = Function(Q)
tau_jd = Function(Q)
tau_ii = Function(Q)
tau_ij = Function(Q)
tau_iz = Function(Q)
tau_ji = Function(Q)
tau_jj = Function(Q)
tau_jz = Function(Q)

us     = Function(Q)
vs     = Function(Q)
ws     = Function(Q)

f = HDF5File(mesh.mpi_comm(), var_dir + 'vars.h5', 'r')

f.read(S,     'S')
f.read(B,     'B')
f.read(T_s,   'T_s')
f.read(U_ob,  'U_ob')
f.read(q_geo, 'q_geo')
f.read(adot,  'adot')
f.read(ff,    'ff')
f.read(cf,    'cf')
f.read(ff_acc,'ff_acc')
f.read(mask,  'mask')

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
model.init_T(in_dir + 'T.xml')
model.init_Mb(in_dir + 'Mb.xml')
model.init_W(in_dir + 'W.xml')
model.init_Ubar(bv_dir + 'Ubar.xml')
model.init_etabar(avg_dir + 'etabar.xml')
model.init_component_Ubar(avg_dir + 'ubar.xml',
                          avg_dir + 'vbar.xml')

model.assign_variable(tau_id, str_dir + 'tau_id.xml')
model.assign_variable(tau_jd, str_dir + 'tau_jd.xml')
model.assign_variable(tau_ii, str_dir + 'tau_ii.xml')
model.assign_variable(tau_ij, str_dir + 'tau_ij.xml')
model.assign_variable(tau_iz, str_dir + 'tau_iz.xml')
model.assign_variable(tau_ji, str_dir + 'tau_ji.xml')
model.assign_variable(tau_jj, str_dir + 'tau_jj.xml')
model.assign_variable(tau_jz, str_dir + 'tau_jz.xml')

model.assign_variable(us, avg_dir + 'us.xml')
model.assign_variable(vs, avg_dir + 'vs.xml')
model.assign_variable(ws, avg_dir + 'ws.xml')

submesh = model.get_bed_mesh()

Q_b       = FunctionSpace(submesh, 'CG', 1)
beta_s    = Function(Q_b)
S_s       = Function(Q_b)
B_s       = Function(Q_b)
H_s       = Function(Q_b)
Tb_s      = Function(Q_b)
Ts_s      = Function(Q_b)
ub_s      = Function(Q_b)
vb_s      = Function(Q_b)
wb_s      = Function(Q_b)
us_s      = Function(Q_b)
vs_s      = Function(Q_b)
ws_s      = Function(Q_b)
Mb_s      = Function(Q_b)
W_s       = Function(Q_b)
adot_s    = Function(Q_b)
qgeo_s    = Function(Q_b)
U_ob_s    = Function(Q_b)
Ubar_s    = Function(Q_b)
etabar_s  = Function(Q_b)
ubar_s    = Function(Q_b)
vbar_s    = Function(Q_b)
tau_id_s  = Function(Q_b)
tau_jd_s  = Function(Q_b)
tau_ii_s  = Function(Q_b)
tau_ij_s  = Function(Q_b)
tau_iz_s  = Function(Q_b)
tau_ji_s  = Function(Q_b)
tau_jj_s  = Function(Q_b)
tau_jz_s  = Function(Q_b)
mask_s    = Function(Q_b)


lg      = LagrangeInterpolator()

lg.interpolate(beta_s,   model.beta)
lg.interpolate(Tb_s,     model.T)
lg.interpolate(Ts_s,     T_s)
lg.interpolate(S_s,      S)
lg.interpolate(B_s,      B)
lg.interpolate(ub_s,     model.u)
lg.interpolate(vb_s,     model.v)
lg.interpolate(wb_s,     model.w)
lg.interpolate(us_s,     us)
lg.interpolate(vs_s,     vs)
lg.interpolate(ws_s,     ws)
lg.interpolate(Mb_s,     model.Mb)
lg.interpolate(W_s,      model.W)
lg.interpolate(adot_s,   adot)
lg.interpolate(qgeo_s,   q_geo)
lg.interpolate(U_ob_s,   U_ob)
lg.interpolate(Ubar_s,   model.Ubar)
lg.interpolate(etabar_s, model.etabar)
lg.interpolate(ubar_s,   model.ubar)
lg.interpolate(vbar_s,   model.vbar)
lg.interpolate(tau_id_s, tau_id)
lg.interpolate(tau_jd_s, tau_jd)
lg.interpolate(tau_ii_s, tau_ii)
lg.interpolate(tau_ij_s, tau_ij)
lg.interpolate(tau_iz_s, tau_iz)
lg.interpolate(tau_ji_s, tau_ji)
lg.interpolate(tau_jj_s, tau_jj)
lg.interpolate(tau_jz_s, tau_jz)
lg.interpolate(mask_s,   mask)

XDMFFile(submesh.mpi_comm(), out_dir + 'submesh.xdmf') << submesh

model.save_xml(beta_s,   'beta_s')
model.save_xml(Tb_s,     'Tb_s')
model.save_xml(Ts_s,     'Ts_s')
model.save_xml(S_s,      'S_s')
model.save_xml(B_s,      'B_s')
model.save_xml(ub_s,     'ub_s')
model.save_xml(vb_s,     'vb_s')
model.save_xml(wb_s,     'wb_s')
model.save_xml(us_s,     'us_s')
model.save_xml(vs_s,     'vs_s')
model.save_xml(ws_s,     'ws_s')
model.save_xml(Mb_s,     'Mb_s')
model.save_xml(W_s,      'W_s')
model.save_xml(adot_s,   'adot_s')
model.save_xml(qgeo_s,   'qgeo_s')
model.save_xml(U_ob_s,   'U_ob_s')
model.save_xml(Ubar_s,   'Ubar_s')
model.save_xml(etabar_s, 'etabar_s')
model.save_xml(ubar_s,   'ubar_s')
model.save_xml(vbar_s,   'vbar_s')
model.save_xml(tau_id_s, 'tau_id_s')
model.save_xml(tau_jd_s, 'tau_jd_s')
model.save_xml(tau_ii_s, 'tau_ii_s')
model.save_xml(tau_ij_s, 'tau_ij_s')
model.save_xml(tau_iz_s, 'tau_iz_s')
model.save_xml(tau_ji_s, 'tau_ji_s')
model.save_xml(tau_jj_s, 'tau_jj_s')
model.save_xml(tau_jz_s, 'tau_jz_s')
model.save_xml(mask_s,   'mask_s')



