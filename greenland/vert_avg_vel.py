import varglas.model              as model
from varglas.mesh.mesh_factory    import MeshFactory
from varglas.data.data_factory    import DataFactory
from varglas.utilities            import DataInput, DataOutput
from fenics                       import *

# set the in/out directory :
out_dir = 'dump/vert_avg_vel/'
in_dir  = 'dump/test/04/'

thklim  = 200.0

# collect the raw data :
bamber  = DataFactory.get_bamber(thklim = thklim)

# define the mesh :
mesh    = MeshFactory.get_greenland_detailed()

# create data objects to use with varglas :
dbm     = DataInput(bamber,   mesh=mesh)

# get the expressions used by varglas :
S = dbm.get_expression('S',    near=True)
B = dbm.get_expression('B',    near=True)
M = dbm.get_expression('mask', near=True)

model = model.Model()
model.set_mesh(mesh)
model.set_geometry(S, B, deform=True)
model.calculate_boundaries(M)
model.initialize_variables()

File(in_dir + 'beta.xml') >> model.beta
File(in_dir + 'u.xml')    >> model.u
File(in_dir + 'v.xml')    >> model.v
File(in_dir + 'w.xml')    >> model.w
  
ubar = model.calc_vert_average(model.u)
vbar = model.calc_vert_average(model.v)
beta = model.extrude(model.beta, [3,5], 2)

File(out_dir + 'ubar.pvd') << ubar
File(out_dir + 'vbar.pvd') << vbar
File(out_dir + 'beta.pvd') << beta

# write data to matlab matrix file :
do = DataOutput(out_dir)
do.write_matlab(dbm, ubar, 'ubar')
do.write_matlab(dbm, vbar, 'ubar')
do.write_matlab(dbm, beta, 'ubar')



