import varglas.solvers            as solvers
import varglas.physics            as physics
import varglas.model              as model
from varglas.helper               import default_config
from fenics                       import *
from varglas.data.data_factory    import DataFactory
from varglas.mesh.mesh_factory    import MeshFactory
from varglas.io                   import DataInput, DataOutput


# get the input args :
out_dir = 'dump/bed/09/smb/'
in_dir  = 'dump/bed/09/'

mesh   = Mesh(in_dir + 'submesh.xdmf')
Q      = FunctionSpace(mesh, 'CG', 1)

S      = Function(Q)
B      = Function(Q)

File(in_dir + 'S_s.xml') >> S
File(in_dir + 'B_s.xml') >> B

config = default_config()
config['output_path']  = out_dir
config['model_order']  = 'SSA'

model = model.Model(config)
model.set_mesh(mesh)
model.set_surface_and_bed(S, B)
model.initialize_variables()

model.init_adot(in_dir + 'adot_s.xml')
model.init_Mb(in_dir + 'Mb_s.xml')
model.init_component_Ubar(in_dir + 'ubar_s.xml',
                          in_dir + 'vbar_s.xml',
                          in_dir + 'wbar_s.xml')

model.save_pvd(model.adot, 'adot')

F = physics.SurfaceMassBalance(model, config)

F.solve()

model.save_pvd(model.adot, 'smb')
model.save_xml(model.adot, 'smb')

bamber  = DataFactory.get_bamber(thklim = 0.0)
dbm     = DataInput(bamber,  mesh=mesh)

do = DataOutput(out_dir)
do.write_matlab(dbm, model.adot, 'smb', val=0.0)



