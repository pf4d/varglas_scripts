import varglas.model              as model
from varglas.mesh.mesh_factory    import MeshFactory
from varglas.data.data_factory    import DataFactory
from varglas.utilities            import DataInput, DataOutput
from fenics                       import *

# set the in/out directory :
out_dir = 'dump/vert_avg_vel/'
in_dir  = 'dump/test/04/'
var_dir = 'dump/vars/'

mesh   = MeshFactory.get_greenland_detailed()
#mesh   = Mesh(var_dir + 'mesh.xdmf')
Q      = FunctionSpace(mesh, 'CG', 1)
ff     = MeshFunction('size_t', mesh)
cf     = MeshFunction('size_t', mesh)
ff_acc = MeshFunction('size_t', mesh)
S      = Function(Q)
B      = Function(Q)

f = HDF5File(mesh.mpi_comm(), var_dir + 'vars.h5', 'r')

f.read(ff,     'ff')
f.read(cf,     'cf')
f.read(ff_acc, 'ff_acc')
f.read(S,      'S')
f.read(B,      'B')

model = model.Model()
model.set_mesh(mesh)
model.set_surface_and_bed(S, B)
model.set_subdomains(ff, cf, ff_acc)
model.initialize_variables()

#File(in_dir + 'beta.xml') >> model.beta
#File(in_dir + 'u.xml')    >> model.u
#File(in_dir + 'v.xml')    >> model.v
#File(in_dir + 'w.xml')    >> model.w
#  
#ubar = model.calc_vert_average(model.u)
#vbar = model.calc_vert_average(model.v)
#beta = model.extrude(model.beta, [3,5], 2)
#
#File(out_dir + 'ubar.pvd') << ubar
#File(out_dir + 'vbar.pvd') << vbar
#File(out_dir + 'beta.pvd') << beta
#
#File(out_dir + 'ubar.xml') << ubar
#File(out_dir + 'vbar.xml') << vbar
#File(out_dir + 'beta.xml') << beta

ubar = Function(Q)
vbar = Function(Q)
beta = Function(Q)

File(out_dir + 'ubar.xml') >> ubar
File(out_dir + 'vbar.xml') >> vbar
File(out_dir + 'beta.xml') >> beta

# write data to matlab matrix file :
bamber = DataFactory.get_bamber(thklim = 1.0)
dbm    = DataInput(bamber, mesh=mesh)
do     = DataOutput(out_dir)
do.write_matlab(dbm, ubar, 'ubar')
do.write_matlab(dbm, vbar, 'vbar')
do.write_matlab(dbm, beta, 'beta')



