import varglas.model              as model
from varglas.mesh.mesh_factory    import MeshFactory
from varglas.data.data_factory    import DataFactory
from varglas.io                   import DataInput
from varglas.helper               import default_config
from fenics                       import *
from pylab                        import *

out_dir  = 'dump/vars_jakobshavn/'
thklim   = 1.0

# collect the raw data :
searise  = DataFactory.get_searise(thklim = thklim)
bamber   = DataFactory.get_bamber(thklim = thklim)
rignot   = DataFactory.get_rignot()

# define the mesh :
mesh = Mesh('dump/meshes/jakobshavn_3D_5H_mesh_block.xml.gz')

# create data objects to use with varglas :
dsr     = DataInput(searise,  mesh=mesh)
dbm     = DataInput(bamber,   mesh=mesh)
drg     = DataInput(rignot,   mesh=mesh)

m = dbm.data['mask']

# calculate surface gradient :
gradM = gradient(m)
gM    = sqrt(gradM[0]**2 + gradM[1]**2 + 1e-16)

gM[gM > 0.1] = 100.0
gM[gM < 100] = 0.0

ref = m - gM

ref[ref > 1]    = 100
ref[ref < 100]  = 1
ref[ref == 100] = 0

dbm.data['ref'] = ref

# change the projection of all data to Rignot projection :
drg.change_projection(dbm)

# get the expressions used by varglas :
S     = dbm.get_expression('S',        near=False, kx=1, ky=1)
B     = dbm.get_expression('B',        near=False, kx=1, ky=1)
M     = dbm.get_expression('ref',      near=True)
adot  = dsr.get_expression('adot',     near=False)
T_s   = dsr.get_interpolation('T',     near=False)
u     = drg.get_interpolation('vx',    near=False)
v     = drg.get_interpolation('vy',    near=False)

config = default_config()
config['output_path'] = out_dir

model = model.Model(config)
model.set_mesh(mesh)
model.calculate_boundaries(mask=M, adot=adot)
model.set_geometry(S, B, deform=True)

model.save_pvd(model.ff, 'ff')

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
f.write(adot,         'adot')
f.write(mask,         'mask')
f.write(u,            'u')
f.write(v,            'v')


