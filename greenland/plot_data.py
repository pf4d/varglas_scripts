import varglas.physical_constants as pc
import varglas.model              as model
from varglas.mesh.mesh_factory    import MeshFactory
from varglas.data.data_factory    import DataFactory
from varglas.io                   import DataInput
from varglas.helper               import plotIce
from fenics                       import *
from pylab                        import *

# collect the raw data :
searise  = DataFactory.get_searise()
bamber   = DataFactory.get_bamber()
rignot   = DataFactory.get_rignot()

# create data objects to use with varglas :
dsr     = DataInput(searise,  gen_space=False)
dbm     = DataInput(bamber,   gen_space=False)
drg     = DataInput(rignot,   gen_space=False)

# change the projection of all data to Rignot projection :
#drg.change_projection(dbm)

## get the expressions used by varglas :
#S     = dbm.get_expression('S',        near=True)
#B     = dbm.get_expression('B',        near=True)
#M     = dbm.get_expression('mask',     near=True)
#adot  = dsr.get_expression('adot',     near=True)
#T_s   = dsr.get_interpolation('T',     near=True)
#q_geo = dfm.get_interpolation('q_geo', near=True)
#u     = drg.get_interpolation('vx',    near=True)
#v     = drg.get_interpolation('vy',    near=True)

#model = model.Model()
#model.set_mesh(mesh)
#model.calculate_boundaries(mask=M, adot=adot)
#model.set_geometry(S, B, deform=True)
#model.set_parameters(pc.IceParameters())

#m = dbm.data['mask']
#
## calculate surface gradient :
#gradM = gradient(m)
#gM    = sqrt(gradM[0]**2 + gradM[1]**2 + 1e-16)
#dbm.data['gM'] = gM
#
#gM[gM > 0.1] = 100.0
#gM[gM < 100] = 0.0
#
#ref = m - gM
#
#ref[ref > 1]    = 100
#ref[ref < 100]  = 1
#ref[ref == 100] = 0
#
#dbm.data['ref'] = ref
#
#num_mask = len(unique(m))

#x1 = -500000; y1 = -2190000
#x2 = -270000; y2 = -2320000
#
#x = dbm.x
#y = dbm.y
#
#x_valid  = where(x > x1)[0]
#x_valid  = intersect1d(x_valid, where(x < x2)[0])
#
#y_valid  = where(y < y1)[0]
#y_valid  = intersect1d(y_valid, where(y > y2)[0])
#
#for i in y_valid:
#  for j in x_valid:
#    dbm.data['H'][i,j] = 1.0

#plotIce(dbm, 'mask', name='mask', direc='images/data/', cmap='gist_yarg',
#        title=r'', scale='lin', umin=None, umax=None, numLvls=num_mask)
#
#plotIce(dbm, 'gM', name='gM', direc='images/data/', cmap='gist_yarg',
#        title=r'', scale='lin', umin=None, umax=None, numLvls=25)
#
#plotIce(dbm, 'ref', name='ref', direc='images/data/', cmap='gist_yarg',
#        title=r'', scale='lin', umin=None, umax=None, numLvls=3)
#
plotIce(dbm, 'H', name='H', direc='images/data/', cmap='gist_yarg',
        title=r'$H$', scale='lin', umin=None, umax=None, numLvls=25)

#plotIce(dbm, 'S', name='S', direc='images/data/', cmap='gist_yarg',
#        title=r'$S$', scale='lin', umin=0.0, umax=None, numLvls=25)
#
#plotIce(dsr, 'U_sar', name='U_ob', direc='images/data/', cmap='gist_yarg',
#        title=r'$\Vert \mathbf{u}_{ob} \Vert$', scale='log',
#        umin=1.0, umax=4000, numLvls=25)



