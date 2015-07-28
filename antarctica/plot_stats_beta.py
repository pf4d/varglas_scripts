from fenics                    import *
from pylab                     import *
from varglas.helper            import default_config, plotIce
from varglas.data.data_factory import DataFactory
from varglas.io                import DataInput

out_dir  = 'dump/linear_model/'
in_dir   = 'dump/linear_model/'

submesh = Mesh('dump/bed/07/submesh.xdmf')

Q     = FunctionSpace(submesh, 'CG', 1)
beta  = Function(Q)
us    = Function(Q)
vs    = Function(Q)
ws    = Function(Q)

File(in_dir + 'beta_s.xml') >> beta
File(in_dir + 'us_s.xml')   >> us
File(in_dir + 'vs_s.xml')   >> vs
File(in_dir + 'ws_s.xml')   >> ws

us_v       = us.vector().array()
vs_v       = vs.vector().array()
ws_v       = ws.vector().array()
Us_mag     = sqrt(us_v**2 + vs_v**2 + ws_v**2 + 1e-16)
Us_mag_f   = Function(Q)
Us_mag_f.vector()[:] = Us_mag


measures   = DataFactory.get_ant_measures(res=900)
dm         = DataInput(measures, gen_space=False)

plotIce(dm, beta,  name='beta', direc=out_dir, cmap='gist_yarg',
        title=r'$\beta$', scale='lin', umin=0.0, umax=200,
        numLvls=12, tp=False, tpAlpha=0.5, extend='max', show=False)

plotIce(dm, Us_mag_f, name='Us', direc=out_dir, cmap='gist_yarg',
        title=r'$\Vert\mathbf{u}_S\Vert$', scale='log', umin=1.0, umax=4000,
        numLvls=12, tp=False, tpAlpha=0.5, extend='both', show=False)
