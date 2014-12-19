import varglas.physical_constants as pc
import varglas.model              as model
from varglas.mesh.mesh_factory    import MeshFactory
from varglas.data.data_factory    import DataFactory
from varglas.utilities            import DataInput, DataOutput
from fenics                       import *
from pylab                        import *

thklim  = 1.0
out_dir = 'test/01/'

measures  = DataFactory.get_ant_measures(res=900)
bedmap1   = DataFactory.get_bedmap1(thklim=thklim)
bedmap2   = DataFactory.get_bedmap2(thklim=thklim)

#mesh = MeshFactory.get_antarctica_3D_gradS_detailed()
mesh  = Mesh('meshes/mesh.xml')

dm  = DataInput(measures, mesh=mesh)
db1 = DataInput(bedmap1,  mesh=mesh)
db2 = DataInput(bedmap2,  mesh=mesh)

db2.data['B'] = db2.data['S'] - db2.data['H']
db2.set_data_val('H', 32767, thklim)
db2.data['S'] = db2.data['B'] + db2.data['H']

H      = db2.get_nearest_expression("H")
S      = db2.get_nearest_expression("S")
B      = db2.get_nearest_expression("B")
M      = db2.get_nearest_expression("mask")
T_s    = db1.get_nearest_expression("temp")
q_geo  = db1.get_nearest_expression("ghfsr")
adot   = db1.get_nearest_expression("acca")
U_ob   = dm.get_projection("U_ob", near=True)
u      = dm.get_nearest_expression("vx")
v      = dm.get_nearest_expression("vy")

model = model.Model()
model.set_mesh(mesh)
model.calculate_boundaries(mask=M, adot=adot)
model.set_geometry(S, B,deform=True)
model.set_parameters(pc.IceParameters())
model.initialize_variables()

File(out_dir + 'beta.xml') >> model.beta
File(out_dir + 'u.xml')    >> model.u
File(out_dir + 'v.xml')    >> model.v
File(out_dir + 'w.xml')    >> model.w
File(out_dir + 'T.xml')    >> model.T
File(out_dir + 'Mb.xml')   >> model.Mb

bmesh   = BoundaryMesh(model.mesh, 'exterior')
cellmap = bmesh.entity_map(2)
pb      = CellFunction("size_t", bmesh, 0)
for c in cells(bmesh):
  if Facet(mesh, cellmap[c.index()]).normal().z() < 0:
    pb[c] = 1
submesh = SubMesh(bmesh, pb, 1)           # subset of surface mesh

Q_b     = FunctionSpace(submesh, 'CG', 1)
beta_s  = Function(Q_b)
S_s     = Function(Q_b)
B_s     = Function(Q_b)
H_s     = Function(Q_b)
T_s     = Function(Q_b)
u_s     = Function(Q_b)
v_s     = Function(Q_b)
w_s     = Function(Q_b)
Mb_s    = Function(Q_b)
adot_s  = Function(Q_b)

lg      = LagrangeInterpolator()

lg.interpolate(beta_s, model.beta)
lg.interpolate(T_s,    model.T)
lg.interpolate(S_s,    model.S)
lg.interpolate(B_s,    model.B)
lg.interpolate(H_s,    interpolate(H, model.Q))
lg.interpolate(u_s,    model.u)
lg.interpolate(v_s,    model.v)
lg.interpolate(w_s,    model.w)
lg.interpolate(Mb_s,   model.Mb)
lg.interpolate(adot_s, interpolate(model.adot, model.Q))

File('test/bed/beta_s.xml') << beta_s
File('test/bed/Mb_s.xml')   << Mb_s
File('test/bed/T_s.xml')    << T_s
File('test/bed/S_s.xml')    << S_s
File('test/bed/B_s.xml')    << B_s
File('test/bed/H_s.xml')    << H_s
File('test/bed/u_s.xml')    << u_s
File('test/bed/v_s.xml')    << v_s
File('test/bed/w_s.xml')    << w_s
File('test/bed/adot_s.xml') << adot_s



