import varglas.model as model
from fenics          import *
from pylab           import *

out_dir  = 'dump/bed/08/'
in_dir   = 'dump/ant_spacing/08/'
bv_dir   = 'dump/ant_spacing/balance_velocity/'
var_dir  = 'dump/vars_ant_spacing/'
str_dir  = 'dump/stress/'

mesh   = Mesh(var_dir + 'mesh.xdmf')
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

tau_dn = Function(Q)
tau_dt = Function(Q)
tau_bn = Function(Q)
tau_bt = Function(Q)
tau_nn = Function(Q)
tau_nt = Function(Q)
tau_nz = Function(Q)
tau_tn = Function(Q)
tau_tt = Function(Q)
tau_tz = Function(Q)

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

File(str_dir + 'tau_dn.xml') >> tau_dn
File(str_dir + 'tau_dt.xml') >> tau_dt
File(str_dir + 'tau_bn.xml') >> tau_bn
File(str_dir + 'tau_bt.xml') >> tau_bt
File(str_dir + 'tau_nn.xml') >> tau_nn
File(str_dir + 'tau_nt.xml') >> tau_nt
File(str_dir + 'tau_nz.xml') >> tau_nz
File(str_dir + 'tau_tn.xml') >> tau_tn
File(str_dir + 'tau_tt.xml') >> tau_tt
File(str_dir + 'tau_tz.xml') >> tau_tz

model = model.Model()
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

submesh = model.get_bed_mesh()

Q_b       = FunctionSpace(submesh, 'CG', 1)
beta_s    = Function(Q_b)
S_s       = Function(Q_b)
B_s       = Function(Q_b)
H_s       = Function(Q_b)
Tb_s      = Function(Q_b)
Ts_s      = Function(Q_b)
u_s       = Function(Q_b)
v_s       = Function(Q_b)
w_s       = Function(Q_b)
Mb_s      = Function(Q_b)
W_s       = Function(Q_b)
adot_s    = Function(Q_b)
qgeo_s    = Function(Q_b)
U_ob_s    = Function(Q_b)
Ubar_s    = Function(Q_b)
tau_dn_s  = Function(Q_b)
tau_dt_s  = Function(Q_b)
tau_nn_s  = Function(Q_b)
tau_nt_s  = Function(Q_b)
tau_nz_s  = Function(Q_b)
tau_tn_s  = Function(Q_b)
tau_tt_s  = Function(Q_b)
tau_tz_s  = Function(Q_b)

lg      = LagrangeInterpolator()

#lg.interpolate(beta_s,   model.beta)
#lg.interpolate(Tb_s,     model.T)
#lg.interpolate(Ts_s,     T_s)
#lg.interpolate(S_s,      S)
#lg.interpolate(B_s,      B)
#lg.interpolate(u_s,      model.u)
#lg.interpolate(v_s,      model.v)
#lg.interpolate(w_s,      model.w)
#lg.interpolate(Mb_s,     model.Mb)
#lg.interpolate(W_s,      model.W)
#lg.interpolate(adot_s,   adot)
#lg.interpolate(qgeo_s,   q_geo)
#lg.interpolate(U_ob_s,   U_ob)
#lg.interpolate(Ubar_s,   model.Ubar)
#lg.interpolate(tau_dn_s, tau_dn)
#lg.interpolate(tau_dt_s, tau_dt)
#lg.interpolate(tau_bn_s, tau_bn)
#lg.interpolate(tau_bt_s, tau_bt)
#lg.interpolate(tau_nn_s, tau_nn)
#lg.interpolate(tau_nt_s, tau_nt)
lg.interpolate(tau_nz_s, tau_nz)
#lg.interpolate(tau_tn_s, tau_tn)
#lg.interpolate(tau_tt_s, tau_tt)
lg.interpolate(tau_tz_s, tau_tz)

XDMFFile(submesh.mpi_comm(),  out_dir + 'submesh.xdmf') << submesh

#File(out_dir + 'beta_s.xml')    << beta_s
#File(out_dir + 'Mb_s.xml')      << Mb_s
#File(out_dir + 'W_s.xml')       << W_s
#File(out_dir + 'Tb_s.xml')      << Tb_s
#File(out_dir + 'Ts_s.xml')      << Ts_s
#File(out_dir + 'S_s.xml')       << S_s
#File(out_dir + 'B_s.xml')       << B_s
#File(out_dir + 'u_s.xml')       << u_s
#File(out_dir + 'v_s.xml')       << v_s
#File(out_dir + 'w_s.xml')       << w_s
#File(out_dir + 'adot_s.xml')    << adot_s
#File(out_dir + 'qgeo_s.xml')    << qgeo_s
#File(out_dir + 'U_ob_s.xml')    << U_ob_s
#File(out_dir + 'Ubar_s.xml')    << Ubar_s
#File(out_dir + 'tau_dn_s.xml')  << tau_dn_s
#File(out_dir + 'tau_dt_s.xml')  << tau_dt_s
#File(out_dir + 'tau_bn_s.xml')  << tau_bn_s
#File(out_dir + 'tau_bt_s.xml')  << tau_bt_s
#File(out_dir + 'tau_nn_s.xml')  << tau_nn_s
#File(out_dir + 'tau_nt_s.xml')  << tau_nt_s
File(out_dir + 'tau_nz_s.xml')  << tau_nz_s
#File(out_dir + 'tau_tn_s.xml')  << tau_tn_s
#File(out_dir + 'tau_tt_s.xml')  << tau_tt_s
File(out_dir + 'tau_tz_s.xml')  << tau_tz_s



