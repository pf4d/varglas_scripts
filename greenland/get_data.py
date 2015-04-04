from fenics import *

in_dir  = 'dump/vars/'
out_dir = 'dump/data/'

mesh    = Mesh(in_dir + 'mesh.xdmf')
Q       = FunctionSpace(mesh, 'CG', 1)

S       = Function(Q)
B       = Function(Q)
adot    = Function(Q)
q_geo   = Function(Q)
u_ob    = Function(Q)
v_ob    = Function(Q)

f = HDF5File(mesh.mpi_comm(), in_dir + 'vars.h5', 'r')

f.read(S,     'S')
f.read(B,     'B')
f.read(q_geo, 'q_geo')
f.read(adot,  'adot')
f.read(u_ob,  'u')
f.read(v_ob,  'v')

File(out_dir + 'S.pvd')     << S
File(out_dir + 'B.pvd')     << B
File(out_dir + 'q_geo.pvd') << q_geo
File(out_dir + 'adot.pvd')  << adot
File(out_dir + 'u_ob.pvd')  << u_ob
File(out_dir + 'v_ob.pvd')  << v_ob



