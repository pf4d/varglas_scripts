import sys
src_directory = '../../statistical_modeling'
sys.path.append(src_directory)

from src.regstats                 import linRegstats
from varglas.mesh.mesh_factory    import MeshFactory
from fenics                       import *
from pylab                        import *


#===============================================================================
# get the data from the model output on the bed :

mesh  = MeshFactory.get_antarctica_3D_10k()

bmesh   = BoundaryMesh(mesh, 'exterior')
cellmap = bmesh.entity_map(2)
pb      = CellFunction("size_t", bmesh, 0)
for c in cells(bmesh):
  if Facet(mesh, cellmap[c.index()]).normal().z() < 0:
    pb[c] = 1
submesh = SubMesh(bmesh, pb, 1)           # subset of surface mesh

Q_b   = FunctionSpace(submesh, 'CG', 1)
beta  = Function(Q_b)
S     = Function(Q_b)
B     = Function(Q_b)
H     = Function(Q_b)
T     = Function(Q_b)
u     = Function(Q_b)
v     = Function(Q_b)
w     = Function(Q_b)
Mb    = Function(Q_b)
adot  = Function(Q_b)

File('test/bed/beta_s.xml') >> beta
File('test/bed/Mb_s.xml')   >> Mb
File('test/bed/T_s.xml')    >> T
File('test/bed/S_s.xml')    >> S
File('test/bed/B_s.xml')    >> B
File('test/bed/H_s.xml')    >> H
File('test/bed/u_s.xml')    >> u
File('test/bed/v_s.xml')    >> v
File('test/bed/w_s.xml')    >> w
File('test/bed/adot_s.xml') >> adot

dSdx   = project(S.dx(0), Q_b)
dSdy   = project(S.dx(1), Q_b)

dBdx   = project(B.dx(0), Q_b)
dBdy   = project(B.dx(1), Q_b)

# vectors :
beta_v = beta.vector().array()
S_v    = S.vector().array()
B_v    = B.vector().array()
H_v    = S_v - B_v
T_v    = T.vector().array()
u_v    = u.vector().array()
v_v    = v.vector().array()
w_v    = w.vector().array()
Mb_v   = Mb.vector().array()
U_mag  = sqrt(u_v**2 + v_v**2 + w_v**2 + 1e-16)
dSdx_v = dSdx.vector().array()
dSdy_v = dSdy.vector().array()
gradS  = sqrt(dSdx_v**2 + dSdy_v**2 + 1e-16)
dBdx_v = dBdx.vector().array()
dBdy_v = dBdy.vector().array()
gradB  = sqrt(dBdx_v**2 + dBdy_v**2 + 1e-16)

valid  = intersect1d(where(beta_v > 1e-7)[0], where(U_mag < 800)[0])

x1   = Mb_v[valid]
x2   = S_v[valid]
x3   = 273.15 - T_v[valid]
x4   = gradS[valid]
x5   = abs(B_v[valid])
x6   = gradB[valid]
x7   = H_v[valid]
x8   = u_v[valid]
x9   = v_v[valid]
x10  = w_v[valid]
x11  = U_mag[valid]

#===============================================================================
# formulte design matrix and do some EDA :
names = [r'$M_b$', 
         r'$S$', 
         r'$-T$', 
         r'$\Vert \nabla S \Vert$', 
         r'$|B|$',
         r'$\Vert \nabla B \Vert$', 
         r'$H$',
         r'$u$', 
         r'$v$', 
         r'$w$', 
         r'$\Vert \mathbf{U} \Vert$']

X    = [x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11]

ii  = [3,4,5,6,[3,4],[3,6],[4,6]]
fig = figure()
Xt  = []

for k,i in enumerate(ii):
  
  ax = fig.add_subplot(330 + k + 1)
  
  if type(i) == list:
    n = ''
    x = 1.0
    for j in i:
      m = j - 1
      x *= X[m]
      n += names[m]
  else:
    m = i - 1
    x = X[m]
    n = names[m]

  Xt.append(x)
  
  j = argsort(x)
  ax.plot(x[j], beta_v[valid][j], 'ko', alpha=0.1)
  ax.set_xlabel(n)
  ax.set_ylabel(r'$\beta$')
  ax.grid()
show()

# cell declustering :
#h_v    = project(CellSize(submesh), Q_b).vector().array()
#A      = sum(h_v)
#wt     = h_v / A
#beta_bar = sum(beta_v * h_v) / A

Xt = array(Xt)                        # design matrix
yt = sqrt(log(beta_v[valid] + 1))     # lhs


#===============================================================================
# perform regression :

out  = linRegstats(Xt, yt, 0.95)

yhat = out['yhat']
bhat = out['bhat']
cibl = out['CIB'][0]
cibh = out['CIB'][1]

print "<F_pval, pval>:", out['F_pval'], out['pval']

#===============================================================================
# plot y, yhat :
fig  = figure()
ax   = fig.add_subplot(111)

j    = argsort(Xt[0])

ax.plot(Xt[0][j], yt[j],   'ko', alpha=0.1)
ax.plot(Xt[0][j], yhat[j], 'r-', lw=2.0)
ax.set_ylabel(r'$\beta$')
grid()
tight_layout()
show()

#===============================================================================
# plot parameter values with confidence intervals:
fig  = figure()
ax   = fig.add_subplot(111)

j    = argsort(Xt[0])

x    = range(len(ii) + 1)

ax.plot(x, cibh, 'r--', lw=2.0)
ax.plot(x, bhat, 'k-',  lw=2.0)
ax.plot(x, cibl, 'r--', lw=2.0)
ax.set_ylabel(r'$\hat{\beta}_i$')
ax.set_xlabel(r'$i$')
grid()
tight_layout()
show()



