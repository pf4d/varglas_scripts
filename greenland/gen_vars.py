import varglas.physical_constants as pc
import varglas.model              as model
from varglas.mesh.mesh_factory    import MeshFactory
from varglas.data.data_factory    import DataFactory
from varglas.io                   import DataInput
from fenics                       import *

out_dir  = 'dump/vars_ant_spacing/'
thklim   = 1.0

# collect the raw data :
searise  = DataFactory.get_searise(thklim = thklim)
bamber   = DataFactory.get_bamber(thklim = thklim)
fm_qgeo  = DataFactory.get_gre_qgeo_fox_maule()
rignot   = DataFactory.get_gre_rignot()

# define the mesh :
mesh = Mesh('dump/meshes/gre_mesh_ant_spacing.xml')

# create data objects to use with varglas :
dsr     = DataInput(searise,  mesh=mesh)
dbm     = DataInput(bamber,   mesh=mesh)
dfm     = DataInput(fm_qgeo,  mesh=mesh)
drg     = DataInput(rignot,   mesh=mesh)
    
# change the projection of all data to Rignot projection :
dsr.change_projection(drg)
dbm.change_projection(drg)
dfm.change_projection(drg)

# get the expressions used by varglas :
S     = dbm.get_expression('S',        near=True)
B     = dbm.get_expression('B',        near=True)
M     = dbm.get_expression('mask',     near=True)
adot  = dsr.get_expression('adot',     near=True)
T_s   = dsr.get_interpolation('T',     near=True)
q_geo = dfm.get_interpolation('q_geo', near=True)
u     = drg.get_interpolation('vx',    near=True)
v     = drg.get_interpolation('vy',    near=True)
U_ob  = drg.get_interpolation("U_ob",  near=True)

model = model.Model()
model.set_mesh(mesh)
model.calculate_boundaries(mask=M, adot=adot)
model.set_geometry(S, B, deform=True)
model.set_parameters(pc.IceParameters())

adot     = interpolate(adot, model.Q)
mask     = interpolate(M,    model.Q)

XDMFFile(mesh.mpi_comm(),    out_dir + 'mesh.xdmf')    << model.mesh

# save the state of the model :
f = HDF5File(mesh.mpi_comm(), out_dir + 'vars.h5', 'w')
f.write(model.ff,     'ff')
f.write(model.cf,     'cf')
f.write(model.ff_acc, 'ff_acc')
f.write(model.S,      'S')
f.write(model.B,      'B')
f.write(T_s,          'T_s')
f.write(q_geo,        'q_geo')
f.write(adot,         'adot')
f.write(mask,         'mask')
f.write(u,            'u')
f.write(v,            'v')
f.write(U_ob,         'U_ob')



