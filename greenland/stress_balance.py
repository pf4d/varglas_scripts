import sys
import varglas.solvers            as solvers
import varglas.physical_constants as pc
import varglas.model              as model
from varglas.mesh.mesh_factory    import MeshFactory
from varglas.data.data_factory    import DataFactory
from varglas.utilities            import DataInput, DataOutput
from fenics                       import parameters, set_log_active, File
from time                         import time
from termcolor                    import colored, cprint

t0 = time()

out_dir = './stress_balance/'
in_dir  = './test/03/'

set_log_active(True)

thklim = 200.0

# collect the raw data :
bamber = DataFactory.get_bamber(thklim = thklim)

# define the meshes :
mesh = MeshFactory.get_greenland_detailed()

# create data objects to use with varglas :
dbm  = DataInput(bamber, mesh=mesh)

# get the expressions used by varglas :
S    = dbm.get_spline_expression('S')
B    = dbm.get_spline_expression('B')

model = model.Model(out_dir = out_dir)
model.set_mesh(mesh)
model.set_geometry(S, B, deform=True)
model.set_parameters(pc.IceParameters())
model.calculate_boundaries()
model.initialize_variables()
parameters['form_compiler']['quadrature_degree'] = 2

File(in_dir + 'u.xml')    >>  model.u
File(in_dir + 'v.xml')    >>  model.v
File(in_dir + 'w.xml')    >>  model.w
File(in_dir + 'beta.xml') >>  model.beta
File(in_dir + 'eta.xml')  >>  model.eta

config = {'output_path'     : out_dir,
          'log'             : True,
          'solver_params'   : 
          {
            'linear_solver' : 'mumps'
          }}

F = solvers.StokesBalanceSolver(model, config)
F.solve()

tf = time()

# calculate total time to compute
s = tf - t0
m = s / 60.0
h = m / 60.0
s = s % 60
m = m % 60
if model.MPI_rank == 0:
  s    = "Total time to compute: %02d:%02d:%02d" % (h,m,s)
  text = colored(s, 'red', attrs=['bold'])
  print text
        


