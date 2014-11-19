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

#db2.data['mask'][db2.data['mask'] == 1]   = 100
#db2.data['mask'][db2.data['mask'] == 127] = 0
#db2.data['mask'][db2.data['mask'] == 1]   = 0

## plot to check :
#imshow(db2.data['mask'][::-1,:])
#colorbar()
#tight_layout()
#show()

#===============================================================================
# data preparation :
thklim = 0.0

# create meshgrid for contour :
measure = DataFactory.get_ant_measures()

# process the data :
dbm = DataInput(measure, gen_space=False)
#dbm.set_data_max('U_ob', boundary=100.0, val=100.0)
dbm.set_data_min('U_ob', boundary=0.0,   val=0.0)

dbm.data['U_ob'] = log(dbm.data['U_ob'] + e)

# plot to check :
imshow(dbm.data['U_ob'][::-1,:])
colorbar()
tight_layout()
show()


#===============================================================================
# generate the contour :
m = MeshGenerator(db2, 'mesh', 'meshes/')

m.create_contour('mask', zero_cntr=1, skip_pts=20)
m.eliminate_intersections(dist=40)
#m.plot_contour()
m.write_gmsh_contour(boundary_extend=False)
m.extrude(h=100000, n_layers=10)
m.close_file()


#===============================================================================
# refine :
#thklim = 100.0
#db2.set_data_min('ref', boundary=thklim, val=thklim)

ref_bm = MeshRefiner(dbm, 'U_ob', gmsh_file_name='meshes/mesh')
#ref_bm = MeshRefiner(db2, 'ref', gmsh_file_name='mesh')


#===============================================================================
## refine on velocity and divide using inverse option :
#lmax = 70000
#lmin = 2000
#
#a1,a1id = ref_sr.add_linear_attractor(log(1.0), lmin, lmax, inv=True)
#a2,a2id = ref_sr.add_linear_attractor(log(1.0), lmin, lmax, inv=False)
#
#m1  = ref_sr.add_min_field([a1.op, a2.op])
#ref_sr.set_background_field(mid)


#===============================================================================
# refine on thickness :
#a,aid = ref_bm.add_static_attractor()
a,aid = ref_bm.add_linear_attractor(1.0, 100000, 500000, inv=True) 
ref_bm.set_background_field(aid)


#===============================================================================
# finish stuff up :
ref_bm.finish(gui=False)



