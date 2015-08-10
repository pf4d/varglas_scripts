from varglas.meshing           import MeshGenerator, MeshRefiner, GetBasin
from varglas.data.data_factory import DataFactory
from varglas.io                import DataInput, print_min_max
from pylab                     import *
from scipy.interpolate         import interp2d

kappa = 1.0  # ice thickness to refine

#===============================================================================
# data preparation :
out_dir   = 'dump/meshes/'
mesh_name = 'thwaites_3D_%iH' % int(kappa)

# get the data :
measure = DataFactory.get_ant_measures()
bedmap2 = DataFactory.get_bedmap2()

# process the data :
dbm = DataInput(measure, gen_space=False)
db2 = DataInput(bedmap2, gen_space=False)

dbm.set_data_min('U_ob', boundary=0.0, val=0.0)
db2.set_data_val("H",    32767,        0.0)
db2.set_data_val('S',    32767,        0.0)

# calculate surface gradient :
gradS = gradient(db2.data['S'])
gS_n  = sqrt(gradS[0]**2 + gradS[1]**2 + 1e-16)

# create interpolation object to convert measures to bedmap2 coordinates :
interp = interp2d(db2.x, db2.y, gS_n)
gS_n   = interp(dbm.x, dbm.y)
interp = interp2d(db2.x, db2.y, db2.data['mask'])
mask   = interp(dbm.x, dbm.y)

# get dofs for shelves where we restrict the element size :
slp = gS_n >  30
shf = mask >= 0.9
gnd = mask <  0.9
nan = mask >  10


#===============================================================================
# form field from which to refine :
#dbm.data['ref'] = (0.15 + 1/(1 + dbm.data['U_ob'])) * 40000

# restrict element size on the shelves and outside the domain of the data :
#dbm.data['ref'][slp] = 2000.0
#dbm.data['ref'][shf] = 10000.0
#dbm.data['ref'][nan] = 10000.0

db2.data['ref'] = kappa*db2.data['H'].copy()
db2.data['ref'][db2.data['ref'] < kappa*1000.0] = kappa*1000.0

print_min_max(db2.data['ref'], 'ref')

## plot to check :
#imshow(dbm.data['ref'][::-1,:])
#colorbar()
#tight_layout()
#show()


#===============================================================================
# generate the contour :
m = MeshGenerator(db2, mesh_name, out_dir)

m.create_contour('mask', zero_cntr=1, skip_pts=1)

gb = GetBasin(dbm, basin='21', edge_resolution=500)
#gb.extend_edge(20000)
gb.intersection(m.longest_cont)

gb.plot_xycoords_buf(Show=True, other=m.longest_cont)

m.set_contour(gb.get_xy_contour())
m.eliminate_intersections(dist=200)
#m.plot_contour()
m.write_gmsh_contour(boundary_extend=False)
m.extrude(h=100000, n_layers=10)
m.close_file()


#===============================================================================
# refine :
ref_bm = MeshRefiner(db2, 'ref', gmsh_file_name = out_dir + mesh_name)

a,aid = ref_bm.add_static_attractor()
ref_bm.set_background_field(aid)


#===============================================================================
# finish stuff up :
ref_bm.finish(gui = False, out_file_name = out_dir + mesh_name)



