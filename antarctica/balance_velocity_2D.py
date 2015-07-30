import varglas.solvers            as solvers
import varglas.model              as model
from varglas.helper               import default_config
from varglas.data.data_factory    import DataFactory
from varglas.io                   import DataInput, DataOutput
from fenics                       import *

# get the input args :
out_dir = 'dump/bed/07/bv/'
in_dir  = 'dump/bed/07/'

mesh   = Mesh(in_dir + 'submesh.xdmf')
Q      = FunctionSpace(mesh, 'CG', 1)

S      = Function(Q)
B      = Function(Q)

File(in_dir + 'S_s.xml') >> S
File(in_dir + 'B_s.xml') >> B

config = default_config()
config['output_path']               = out_dir
config['balance_velocity']['kappa'] = 5.0
config['model_order']               = 'SSA'

model = model.Model(config)
model.set_mesh(mesh)
model.set_surface_and_bed(S, B)
model.initialize_variables()

model.init_adot(in_dir + 'adot_s.xml')

F = solvers.BalanceVelocitySolver(model, config)

F.solve()

model.save_xml(model.Ubar, 'Ubar_5')

bedmap   = DataFactory.get_bedmap1()
bm       = DataInput(bedmap, gen_space=False)

do = DataOutput(out_dir)
do.write_matlab(bm, model.Ubar, 'Ubar_5', val=0.0)



