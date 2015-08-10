import varglas.model              as model
from varglas.mesh.mesh_factory    import MeshFactory
from varglas.data.data_factory    import DataFactory
from varglas.io                   import DataInput
from varglas.helper               import default_config
from fenics                       import *

out_dir  = 'dump/vars_thwaites/'
thklim   = 1.0
measures = DataFactory.get_ant_measures(res=900)
bedmap1  = DataFactory.get_bedmap1(thklim=thklim)
bedmap2  = DataFactory.get_bedmap2(thklim=thklim)

#mesh = MeshFactory.get_antarctica_3D_gradS_detailed()
#mesh = MeshFactory.get_antarctica_3D_gradS_crude()
#mesh = MeshFactory.get_antarctica_3D_10k()
mesh = Mesh('dump/meshes/thwaites_2D_1H_mesh.xml.gz')

dm = DataInput(measures, mesh=mesh)
d1 = DataInput(bedmap1,  mesh=mesh)
d2 = DataInput(bedmap2,  mesh=mesh)

S     = d2.get_expression("S",        near=True)
B     = d2.get_expression("B",        near=True)
M     = d2.get_expression("mask",     near=True)
adot  = d1.get_expression("acca",     near=True)
T_s   = d1.get_interpolation("temp",  near=True)
q_geo = d1.get_interpolation("ghfsr", near=True)
u     = dm.get_interpolation("vx",    near=True)
v     = dm.get_interpolation("vy",    near=True)
U_ob  = dm.get_interpolation("U_ob",  near=True)

config = default_config()
config['output_path']    = out_dir
config['model_order']    = 'SSA'

model = model.Model()
model.set_mesh(mesh)
#model.calculate_boundaries(mask=M, adot=adot)
#model.set_geometry(S, B, deform=True)
model.set_geometry(S, B, deform=False)

adot     = interpolate(adot, model.Q)
mask     = interpolate(M, model.Q)

# do not need this in 3D (not deformed) :
#XDMFFile(mesh.mpi_comm(),    out_dir + 'mesh.xdmf')    << model.mesh

# save the state of the model :
f = HDF5File(mesh.mpi_comm(), out_dir + 'vars.h5', 'w')
#f.write(model.ff,     'ff')
#f.write(model.cf,     'cf')
#f.write(model.ff_acc, 'ff_acc')
f.write(model.S,      'S')
f.write(model.B,      'B')
f.write(T_s,          'T_s')
f.write(q_geo,        'q_geo')
f.write(adot,         'adot')
f.write(mask,         'mask')
f.write(u,            'u')
f.write(v,            'v')
f.write(U_ob,         'U_ob')



