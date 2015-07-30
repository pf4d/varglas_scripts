import varglas.solvers            as solvers
import varglas.physics            as physics
import varglas.model              as model
from varglas.data.data_factory    import DataFactory
from varglas.mesh.mesh_factory    import MeshFactory
from varglas.io                   import DataInput, DataOutput
from varglas.helper               import default_config
from fenics                       import *
from scipy.io                     import loadmat

# get the input args :
out_dir = 'dump/bed/09/bv_smb/'
in_dir  = 'dump/bed/09/'

mesh   = Mesh(in_dir + 'submesh.xdmf')
Q      = FunctionSpace(mesh, 'CG', 1)

bamber = DataFactory.get_bamber()
dbm    = DataInput(bamber,  mesh=mesh)

S      = Function(Q)
B      = Function(Q)

File(in_dir + 'S_s.xml') >> S
File(in_dir + 'B_s.xml') >> B

adot_v = loadmat(in_dir + 'smb/smb.mat')
dbm.data['smb'] = adot_v['map_data']

adot = dbm.get_expression('smb', near=True)

config = default_config()
config['output_path']               = out_dir
config['balance_velocity']['kappa'] = 5.0
config['model_order']               = 'SSA'

model = model.Model(config)
model.set_mesh(mesh)
model.set_surface_and_bed(S, B)
model.initialize_variables()

#model.init_adot(in_dir + 'adot_s.xml')
model.init_adot(adot)

F = solvers.BalanceVelocitySolver(model, config)

F.solve()

#model.save_xml(model.Ubar, 'Ubar_5')

