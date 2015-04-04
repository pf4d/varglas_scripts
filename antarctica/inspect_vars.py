import varglas.model              as model
from varglas.data.data_factory    import DataFactory
from varglas.io                   import DataInput, print_min_max
from varglas.helper               import plotIce
from fenics                       import *

out_dir  = 'dump/basin_vars_low/'
thklim   = 1.0
bedmap2  = DataFactory.get_bedmap2(thklim=thklim)

print_min_max(bedmap2['S'], 'S')

mesh = Mesh('dump/meshes/byrd_basin_mesh_low.xml')

#d2 = DataInput(bedmap2, gen_space=False)
d2 = DataInput(bedmap2, mesh=mesh)

#plotIce(d2, 'S', '.', 'gist_yarg', scale='lin', name=r'$S$', 
#        numLvls=12, tp=False, tpAlpha=0.5)
#
#plotIce(d2, 'B', '.', 'gist_yarg', scale='lin', name=r'$B$', 
#        numLvls=12, tp=False, tpAlpha=0.5)
#
#plotIce(d2, 'H', '.', 'gist_yarg', scale='lin', name=r'$H$', 
#        numLvls=12, tp=False, tpAlpha=0.5)
#
#plotIce(d2, 'mask', '.', 'gist_yarg', scale='lin', name='mask', 
#        numLvls=12, tp=False, tpAlpha=0.5)


S     = d2.get_expression("S",     )#     near=True)
B     = d2.get_expression("B",     )#     near=True)
M     = d2.get_expression("mask",  near=True)

model = model.Model()
model.set_mesh(mesh)
model.calculate_boundaries(mask=M)
model.set_geometry(S, B, deform=True)

File('crap/S.pvd')  << model.S
File('crap/ff.pvd') << model.ff

