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

thklim = 1.0

measures  = DataFactory.get_ant_measures(res=900)  # res is 900 or 450
bedmap1   = DataFactory.get_bedmap1(thklim=thklim)
bedmap2   = DataFactory.get_bedmap2(thklim=thklim)

mesh = MeshFactory.get_antarctica_3D_gradS_detailed()

dm  = DataInput(None, measures, mesh=mesh)
db1 = DataInput(None, bedmap1,  mesh=mesh)
db2 = DataInput(None, bedmap2,  mesh=mesh)

db2.data['B'] = db2.data['S'] - db2.data['H']
db2.set_data_val('H', 32767, thklim)
db2.data['S'] = db2.data['B'] + db2.data['H']

S     = db2.get_spline_expression("S")
B     = db2.get_spline_expression("B")
M     = db2.get_nearest_expression("mask")
#H     = db2.get_interpolation("H")
#T_s   = db1.get_interpolation("srfTemp")
#q_geo = db1.get_interpolation("q_geo")
#adot  = db1.get_interpolation("adot")
#u     = dm.get_interpolation("vx")
#v     = dm.get_interpolation("vy")

model = model.Model()
model.set_mesh(mesh)
model.set_geometry(S, B, mask=M, deform=True)
model.set_parameters(pc.IceParameters())
model.calculate_boundaries()
model.initialize_variables()

# save the mesh and subdomains :
File(out_dir + 'mesh.xml')       << model.mesh
File(out_dir + 'flat_mesh.xml')  << model.flat_mesh
File(out_dir + 'ff.xml')         << model.ff
File(out_dir + 'ff_flat.xml')    << model.ff_flat

## save the state of the model :
#File(out_dir + 'H.xml')     << H
#File(out_dir + 'S.xml')     << model.S
#File(out_dir + 'B.xml')     << model.B
#File(out_dir + 'T_s.xml')   << T_s
#File(out_dir + 'q_geo.xml') << q_geo
#File(out_dir + 'adot.xml')  << adot
#File(out_dir + 'u_obs.xml') << u
#File(out_dir + 'v_obs.xml') << v


