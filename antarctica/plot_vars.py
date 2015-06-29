from fenics                    import *
from pylab                     import *
from varglas.helper            import plotIce
from varglas.data.data_factory import DataFactory
from varglas.io                import DataInput

in_dir  = 'dump/bed/07/'

mesh  = Mesh(in_dir + 'submesh.xdmf')
Q     = FunctionSpace(mesh, 'CG', 1)

S      = Function(Q)
B      = Function(Q)
adot   = Function(Q)
qgeo   = Function(Q)
beta   = Function(Q)
Mb     = Function(Q)
Tb     = Function(Q)
Ts     = Function(Q)
us     = Function(Q)
vs     = Function(Q)
ws     = Function(Q)
ub     = Function(Q)
vb     = Function(Q)
wb     = Function(Q)
Ubar1  = Function(Q)
Ubar5  = Function(Q)
Ubar10 = Function(Q)
Ubar20 = Function(Q)
Ubar50 = Function(Q)
U_ob   = Function(Q)
etabar = Function(Q)
ubar   = Function(Q)
vbar   = Function(Q)
tau_id = Function(Q)
tau_jd = Function(Q)
tau_ii = Function(Q)
tau_ij = Function(Q)
tau_iz = Function(Q)
tau_ji = Function(Q)
tau_jj = Function(Q)
tau_jz = Function(Q)
mask   = Function(Q)

File(in_dir + 'S_s.xml')      >> S
File(in_dir + 'B_s.xml')      >> B
File(in_dir + 'adot_s.xml')   >> adot
File(in_dir + 'qgeo_s.xml')   >> qgeo
File(in_dir + 'beta_s.xml')   >> beta
File(in_dir + 'Mb_s.xml')     >> Mb
File(in_dir + 'Tb_s.xml')     >> Tb
File(in_dir + 'Ts_s.xml')     >> Ts
File(in_dir + 'us_s.xml')     >> us
File(in_dir + 'vs_s.xml')     >> vs
File(in_dir + 'ws_s.xml')     >> ws
File(in_dir + 'ub_s.xml')     >> ub
File(in_dir + 'vb_s.xml')     >> vb
File(in_dir + 'wb_s.xml')     >> wb
File(in_dir + 'bv/Ubar_1.xml')  >> Ubar1
File(in_dir + 'bv/Ubar_5.xml')  >> Ubar5
File(in_dir + 'bv/Ubar_10.xml') >> Ubar10
File(in_dir + 'bv/Ubar_20.xml') >> Ubar20
File(in_dir + 'bv/Ubar_50.xml') >> Ubar50
File(in_dir + 'etabar_s.xml') >> etabar
File(in_dir + 'ubar_s.xml')   >> ubar
File(in_dir + 'vbar_s.xml')   >> vbar
File(in_dir + 'U_ob_s.xml')   >> U_ob
File(in_dir + 'U_ob_s.xml')   >> U_ob

File(in_dir + 'tau_id_s.xml') >> tau_id
File(in_dir + 'tau_jd_s.xml') >> tau_jd
File(in_dir + 'tau_ii_s.xml') >> tau_ii
File(in_dir + 'tau_ij_s.xml') >> tau_ij
File(in_dir + 'tau_iz_s.xml') >> tau_iz
File(in_dir + 'tau_ji_s.xml') >> tau_ji
File(in_dir + 'tau_jj_s.xml') >> tau_jj
File(in_dir + 'tau_jz_s.xml') >> tau_jz
File(in_dir + 'mask_s.xml')   >> mask

H          = Function(Q)
S_v        = S.vector().array()
B_v        = B.vector().array()
H_v        = S_v - B_v
H.vector()[:] = H_v

dSdx    = project(S.dx(0), Q)
dSdy    = project(S.dx(1), Q)

dBdx    = project(B.dx(0), Q)
dBdy    = project(B.dx(1), Q)

gradS   = Function(Q)
dSdx_v  = dSdx.vector().array()
dSdy_v  = dSdy.vector().array()
gradS_v = sqrt(dSdx_v**2 + dSdy_v**2 + 1e-16)
gradS.vector()[:] = gradS_v

gradB   = Function(Q)
dBdx_v  = dBdx.vector().array()
dBdy_v  = dBdy.vector().array()
gradB_v = sqrt(dBdx_v**2 + dBdy_v**2 + 1e-16)
gradB.vector()[:] = gradB_v

D      = Function(Q)
D_v    = zeros(len(B_v))
D_v[B_v < 0] = B_v[B_v < 0]
D.vector()[:] = D_v

us_v       = us.vector().array()
vs_v       = vs.vector().array()
ws_v       = ws.vector().array()
Us_mag     = sqrt(us_v**2 + vs_v**2 + ws_v**2 + 1e-16)
Us_mag_f   = Function(Q)
Us_mag_f.vector()[:] = Us_mag

ub_v       = ub.vector().array()
vb_v       = vb.vector().array()
wb_v       = wb.vector().array()
Ub_mag     = sqrt(ub_v**2 + vb_v**2 + wb_v**2 + 1e-16)
Ub_mag_f   = Function(Q)
Ub_mag_f.vector()[:] = Ub_mag

ubar_v     = ubar.vector().array()
vbar_v     = vbar.vector().array()
Ubar_mag   = sqrt(ubar_v**2 + vbar_v**2 + 1e-16)
Ubar_mag_f = Function(Q)
Ubar_mag_f.vector()[:] = Ubar_mag

etabar_v   = etabar.vector().array()
etabar_v  /= 1e6
etabar.vector()[:] = etabar_v

measures   = DataFactory.get_ant_measures(res=900)
dm         = DataInput(measures, gen_space=False)

plotIce(dm, mask, name='mask', direc='images/model/', cmap='gist_yarg',
        title=r'', scale='bool', umin=0.0, umax=1.0,
        numLvls=12, tp=False, tpAlpha=0.5, extend='neither', show=False)
plotIce(dm, S, name='S', direc='images/model/', cmap='gist_yarg',
        title=r'$S$', scale='lin', umin=None, umax=None,
        numLvls=12, tp=False, tpAlpha=0.5, extend='neither', show=False)
plotIce(dm, B, name='B', direc='images/model/', cmap='gist_yarg',
        title=r'$B$', scale='lin', umin=None, umax=None,
        numLvls=12, tp=False, tpAlpha=0.5, extend='neither', show=False)
plotIce(dm, gradS, name='gradS', direc='images/model/', cmap='gist_yarg',
        title=r'$\Vert \nabla S \Vert$', scale='log', umin=1e-4, umax=None,
        numLvls=12, tp=False, tpAlpha=0.5, extend='min', show=False)
plotIce(dm, gradB, name='gradB', direc='images/model/', cmap='gist_yarg',
        title=r'$\Vert \nabla B \Vert$', scale='log', umin=1e-4, umax=None,
        numLvls=12, tp=False, tpAlpha=0.5, extend='min', show=False)
plotIce(dm, H, name='H', direc='images/model/', cmap='gist_yarg',
        title=r'$H$', scale='lin', umin=None, umax=None,
        numLvls=12, tp=False, tpAlpha=0.5, extend='neither', show=False)
plotIce(dm, D, name='D', direc='images/model/', cmap='gist_yarg',
        title=r'$D$', scale='lin', umin=None, umax=-1e-12,
        numLvls=12, tp=False, tpAlpha=0.5, extend='neither', show=False)
plotIce(dm, Tb, name='Tb', direc='images/model/', cmap='gist_yarg',
        title=r'$T_B$', scale='lin', umin=240, umax=273.15,
        numLvls=12, tp=False, tpAlpha=0.5, extend='min', show=False)
plotIce(dm, Ts, name='Ts', direc='images/model/', cmap='gist_yarg',
        title=r'$T_S$', scale='lin', umin=210, umax=273.15,
        numLvls=12, tp=False, tpAlpha=0.5, extend='min', show=False)
plotIce(dm, adot, name='adot', direc='images/model/', cmap='gist_yarg',
        title=r'$\dot{a}$', scale='lin', umin=0.0, umax=1.0,
        numLvls=12, tp=False, tpAlpha=0.5, extend='both', show=False)
plotIce(dm, etabar, name='etabar', direc='images/model/', cmap='gist_yarg',
        title=r'$\bar{\eta}$', scale='log', umin=1.0, umax=5e3,
        numLvls=12, tp=False, tpAlpha=0.5, extend='both', show=False)
plotIce(dm, Ub_mag_f, name='Ub', direc='images/model/', cmap='gist_yarg',
        title=r'$\Vert\mathbf{u}_B\Vert$', scale='log', umin=1.0, umax=4000,
        numLvls=12, tp=False, tpAlpha=0.5, extend='both', show=False)
plotIce(dm, Us_mag_f, name='Us', direc='images/model/', cmap='gist_yarg',
        title=r'$\Vert\mathbf{u}_S\Vert$', scale='log', umin=1.0, umax=4000,
        numLvls=12, tp=False, tpAlpha=0.5, extend='both', show=False)
plotIce(dm, Ubar_mag_f, name='Ubar', direc='images/model/', cmap='gist_yarg',
        title=r'$\Vert\mathbf{\bar{u}}\Vert$', scale='log', umin=1.0, umax=4000,
        numLvls=12, tp=False, tpAlpha=0.5, extend='both', show=False)
plotIce(dm, Ubar1, name='Ubar_bv_1', direc='images/model/', cmap='gist_yarg',
        title=r'$\Vert\mathbf{\bar{u}}_{bv1}\Vert$', scale='log',
        umin=1.0, umax=4000, numLvls=12, tp=False, tpAlpha=0.5,
        extend='both', show=False)
plotIce(dm, Ubar5, name='Ubar_bv_5', direc='images/model/', cmap='gist_yarg',
        title=r'$\Vert\mathbf{\bar{u}}_{bv5}\Vert$', scale='log',
        umin=1.0, umax=4000, numLvls=12, tp=False, tpAlpha=0.5,
        extend='both', show=False)
plotIce(dm, Ubar10, name='Ubar_bv_10', direc='images/model/', cmap='gist_yarg',
        title=r'$\Vert\mathbf{\bar{u}}_{bv10}\Vert$', scale='log',
        umin=1.0, umax=4000, numLvls=12, tp=False, tpAlpha=0.5,
        extend='both', show=False)
plotIce(dm, Ubar20, name='Ubar_bv_20', direc='images/model/', cmap='gist_yarg',
        title=r'$\Vert\mathbf{\bar{u}}_{bv20}\Vert$', scale='log',
        umin=1.0, umax=4000, numLvls=12, tp=False, tpAlpha=0.5,
        extend='both', show=False)
plotIce(dm, Ubar50, name='Ubar_bv_50', direc='images/model/', cmap='gist_yarg',
        title=r'$\Vert\mathbf{\bar{u}}_{bv50}\Vert$', scale='log',
        umin=1.0, umax=4000, numLvls=12, tp=False, tpAlpha=0.5,
        extend='both', show=False)
plotIce(dm, U_ob,   name='U_ob', direc='images/model/', cmap='gist_yarg',
        title=r'$\Vert \mathbf{u}_{ob} \Vert$', scale='log',
        umin=1.0, umax=4000, numLvls=12, tp=False, tpAlpha=0.5,
        extend='both', show=False)
plotIce(dm, beta,  name='beta', direc='images/model/', cmap='gist_yarg',
        title=r'$\beta$', scale='log', umin=1.0, umax=200,
        numLvls=12, tp=False, tpAlpha=0.5, extend='max', show=False)
plotIce(dm, Mb,  name='Mb', direc='images/model/', cmap='gist_yarg',
        title=r'$M_B$', scale='log', umin=0.003, umax=1.0,
        numLvls=13, tp=False, tpAlpha=0.5, extend='both', show=False)
plotIce(dm, tau_id,  name='tau_id', direc='images/model/', cmap='RdGy',
        title=r'$\tau_{id}$', scale='lin', umin=-1e5, umax=1e5,
        numLvls=13, tp=False, tpAlpha=0.5, extend='both', show=False)
plotIce(dm, tau_jd,  name='tau_jd', direc='images/model/', cmap='RdGy',
        title=r'$\tau_{jd}$', scale='lin', umin=-1e5, umax=1e5,
        numLvls=13, tp=False, tpAlpha=0.5, extend='both', show=False)
plotIce(dm, tau_ii,  name='tau_ii', direc='images/model/', cmap='RdGy',
        title=r'$\tau_{ii}$', scale='lin', umin=-1e4, umax=1e4,
        numLvls=13, tp=False, tpAlpha=0.5, extend='both', show=False)
plotIce(dm, tau_ij,  name='tau_ij', direc='images/model/', cmap='RdGy',
        title=r'$\tau_{ij}$', scale='lin', umin=-1e4, umax=1e4,
        numLvls=13, tp=False, tpAlpha=0.5, extend='both', show=False)
plotIce(dm, tau_iz,  name='tau_iz', direc='images/model/', cmap='RdGy',
        title=r'$\tau_{iz}$', scale='lin', umin=-1e5, umax=1e5,
        numLvls=13, tp=False, tpAlpha=0.5, extend='both', show=False)
plotIce(dm, tau_ji,  name='tau_ji', direc='images/model/', cmap='RdGy',
        title=r'$\tau_{ji}$', scale='lin', umin=-1e4, umax=1e4,
        numLvls=13, tp=False, tpAlpha=0.5, extend='both', show=False)
plotIce(dm, tau_jj,  name='tau_jj', direc='images/model/', cmap='RdGy',
        title=r'$\tau_{jj}$', scale='lin', umin=-1e4, umax=1e4,
        numLvls=13, tp=False, tpAlpha=0.5, extend='both', show=False)
plotIce(dm, tau_jz,  name='tau_jz', direc='images/model/', cmap='RdGy',
        title=r'$\tau_{jz}$', scale='lin', umin=-1e5, umax=1e5,
        numLvls=13, tp=False, tpAlpha=0.5, extend='both', show=False)



     



