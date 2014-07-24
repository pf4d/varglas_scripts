from varglas.utilities         import DataInput, MeshGenerator, MeshRefiner
from varglas.data.data_factory import DataFactory
from pylab                     import *


#===============================================================================
# data preparation :
thklim = 0.0

# create meshgrid for contour :
bedmap2 = DataFactory.get_bedmap2()

# process the data :
db2 = DataInput(bedmap2, gen_space=False)
#db2.set_data_val("H", 32767, thklim)
#db2.set_data_val('S', 32767, 0.0)
#
#db2.data['B'] = db2.data['S'] - db2.data['H']
#
#gradS = gradient(db2.data['S'])
#gS_n  = sqrt(gradS[0]**2 + gradS[1]**2) + 1
#
#db2.data['ref'] = db2.data['H'] / gS_n

db2.data['mask'][db2.data['mask'] == 1]   = 100
db2.data['mask'][db2.data['mask'] == 127] = 0
db2.data['mask'][db2.data['mask'] == 1]   = 0

# plot to check :
imshow(db2.data['mask'][::-1,:])
colorbar()
tight_layout()
show()


#===============================================================================
# generate the contour :
m = MeshGenerator(db2, 'mesh', '')

m.create_contour('mask', zero_cntr=1, skip_pts=2)
m.eliminate_intersections(dist=40)
m.plot_contour()
m.write_gmsh_contour(boundary_extend=False)
m.extrude(h=100000, n_layers=10)
m.close_file()


#===============================================================================
# refine :
thklim = 100.0
db2.set_data_min('H', boundary=thklim, val=thklim)

ref_b2 = MeshRefiner(db2, 'H', gmsh_file_name='mesh')


#===============================================================================
# refine on thickness :
a,aid = ref_b2.add_static_attractor(10)
ref_b2.set_background_field(aid)


#===============================================================================
# finish stuff up :
ref_b2.finish(gui=False)



