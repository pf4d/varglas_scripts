import sys
import varglas.physical_constants as pc
import varglas.model              as model
from varglas.mesh.mesh_factory    import MeshFactory
from varglas.data.data_factory    import DataFactory
from varglas.utilities            import DataInput
from fenics                       import *

# get the input args :
out_dir = sys.argv[1] + '/'    # directory to save

set_log_active(True)

thklim = 200.0

measures  = DataFactory.get_ant_measures(res=900)
bedmap1   = DataFactory.get_bedmap1(thklim=thklim)
bedmap2   = DataFactory.get_bedmap2(thklim=thklim)

mesh = MeshFactory.get_antarctica_3D_50H()

dm  = DataInput(None, measures, mesh=mesh)
db1 = DataInput(None, bedmap1,  mesh=mesh)
db2 = DataInput(None, bedmap2,  mesh=mesh)

db2.set_data_val('H',   32767, thklim)
db2.set_data_val('S',   32767, 0.0)
dm.set_data_min('U_ob', 0.0,   0.0)

db2.data['B'] = db2.data['S'] - db2.data['H']

H     = db2.get_spline_expression("H")
S     = db2.get_spline_expression("S")
B     = db2.get_spline_expression("B")
M     = db2.get_nearest_expression("mask")
T_s   = db1.get_spline_expression("srfTemp")
q_geo = db1.get_spline_expression("q_geo")
adot  = db1.get_spline_expression("adot")
U_ob  = dm.get_spline_expression("U_ob")

model = model.Model()
model.set_mesh(mesh)
model.set_geometry(S, B, mask=M, deform=True)
model.set_parameters(pc.IceParameters())
model.initialize_variables()

XDMFFile(mesh.mpi_comm(), out_dir + 'mesh.xdmf')   << model.mesh

# save the state of the model :
f = HDF5File(mesh.mpi_comm(), out_dir + '3D_5H_stokes.h5', 'w')
f.write(model.mesh,                 'mesh')
f.write(model.H,                    'H')
f.write(model.S,                    'S')
f.write(model.B,                    'B')
f.write(model.ff,                   'ff')
f.write(interpolate(T_s, model.Q),  'T_s')
f.write(interpolate(q_geo, model.Q),'q_geo')
f.write(interpolate(adot, model.Q), 'adot')
f.write(interpolate(U_ob, model.Q), 'U_ob')



