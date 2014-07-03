from varglas.utilities         import DataInput, MeshGenerator, MeshRefiner
from varglas.data.data_factory import DataFactory
from pylab                     import *


#===============================================================================
# data preparation :
thklim = 0.0

# create meshgrid for contour :
bedmap2 = DataFactory.get_bedmap2()

# process the data :
db2 = DataInput(None, bedmap2, gen_space=False)
db2.set_data_val("H", 32767, thklim)
db2.set_data_val('S', 32767, 0.0)

db2.data['B'] = db2.data['S'] - db2.data['H']

gradS = gradient(db2.data['S'])
gS_n  = sqrt(gradS[0]**2 + gradS[1]**2) + 1

db2.data['ref'] = db2.data['H'] / gS_n

#===============================================================================
# generate the contour :
m = MeshGenerator(db2, 'mesh', '')

m.create_contour('mask', zero_cntr=1, skip_pts=2)
m.eliminate_intersections(dist=40)
#m.plot_contour()
m.write_gmsh_contour(boundary_extend=False)
m.extrude(h=100000, n_layers=10)
m.close_file()


#===============================================================================
# refine :
thklim = 100.0
db2.set_data_min('ref', boundary=thklim, val=thklim)

## plot to check :
#imshow(db2.data['mask'][::-1,:])
#colorbar()
#tight_layout()
#show()

ref_bm = MeshRefiner(db2, 'ref', gmsh_file_name='mesh')


#===============================================================================
## refine on velocity and divide using inverse option :
#lmax = 70000
#lmin = 2000
#
#a1,a1id = ref_sr.add_linear_attractor(log(1.0), lmin, lmax, inv=True, 
#                                      hard_cut=False)
#a2,a2id = ref_sr.add_linear_attractor(log(1.0), lmin, lmax, inv=False, 
#                                      hard_cut=False)
#
#m1  = ref_sr.add_min_field([a1.op, a2.op])
#ref_sr.set_background_field(mid)


#===============================================================================
# refine on thickness :
a,aid = ref_bm.add_static_attractor(400)
#H     = dbm.data['H']
#a,aid = ref_bm.add_linear_attractor(0, H.min(), H.max(), inv=False, 
#                                    hard_cut=False)
ref_bm.set_background_field(aid)


#===============================================================================
# finish stuff up :
ref_bm.finish(gui=False)



