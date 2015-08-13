import varglas.physical_constants as pc
import varglas.model              as model
from varglas.mesh.mesh_factory    import MeshFactory
from varglas.data.data_factory    import DataFactory
from varglas.io                   import DataInput
from varglas.helper               import plotIce
from fenics                       import *
from pylab                        import *

mesh = Mesh('dump/meshes/jakobshavn_2D_5H_mesh_block.xml.gz')

# collect the raw data :
searise  = DataFactory.get_searise()
bamber   = DataFactory.get_bamber()
rignot   = DataFactory.get_rignot()

# create data objects to use with varglas :
dsr     = DataInput(searise,  mesh = mesh)
dbm     = DataInput(bamber,   mesh = mesh)
drg     = DataInput(rignot,   mesh = mesh)

# change the projection of all data to Rignot projection :
drg.change_projection(dbm)

# get the expressions used by varglas :
#S     = dbm.get_interpolation('S',     near=True)
#B     = dbm.get_interpolation('B',     near=True)
#adot  = dsr.get_interpolation('adot',  near=True)
#T_s   = dsr.get_interpolation('T',     near=True)
u     = drg.get_interpolation('vx',    near=True)
v     = drg.get_interpolation('vy',    near=True)

model = model.Model()
model.config['model_order'] = 'SSA'
model.set_mesh(mesh)
#model.set_surface_and_bed(S, B)
model.initialize_variables()

model.init_U_ob(u, v)

#plotIce(dbm, 'mask', name='mask', direc='images/data/', cmap='gist_yarg',
#        title=r'', scale='lin', umin=None, umax=None, numLvls=num_mask)
#
#plotIce(dbm, 'ref', name='ref', direc='images/data/', cmap='gist_yarg',
#        title=r'', scale='lin', umin=None, umax=None, numLvls=3)
#
#plotIce(dbm, 'H', name='H', direc='images/data/', cmap='gist_yarg',
#        title=r'$H$', scale='lin', umin=None, umax=None, numLvls=25)
#
#plotIce(dbm, 'S', name='S', direc='images/data/', cmap='gist_yarg',
#        title=r'$S$', scale='lin', umin=0.0, umax=None, numLvls=25)
#
plotIce(dbm, model.U_ob, name='U_ob', direc='images/data/', cmap='gist_yarg',
        title=r'$\Vert \mathbf{u}_{ob} \Vert$', scale='log',
        umin=1.0, umax=4000, numLvls=25, tp=True)



