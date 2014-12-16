import sys
src_directory = '../../statistical_modeling'
sys.path.append(src_directory)

from src.regstats                 import linRegstats
from varglas.mesh.mesh_factory    import MeshFactory
from fenics                       import *
from pylab                        import *

mesh  = MeshFactory.get_antarctica_3D_10k()

bmesh   = BoundaryMesh(mesh, 'exterior')
cellmap = bmesh.entity_map(2)
pb      = CellFunction("size_t", bmesh, 0)
for c in cells(bmesh):
  if Facet(mesh, cellmap[c.index()]).normal().z() < 0:
    pb[c] = 1
submesh = SubMesh(bmesh, pb, 1)           # subset of surface mesh

Q_b     = FunctionSpace(submesh, 'CG', 1)
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
x2   = 273.15 - T_v[valid]
x3   = gradS[valid]
x4   = abs(B_v[valid])
x5   = gradB[valid]
x6   = H_v[valid]
x7   = u_v[valid]
x8   = v_v[valid]
x9   = w_v[valid]
x10  = U_mag[valid]
x11  = x5 * x9

i    = argsort(x2)
x1   = x1[i]
x2   = x2[i]
x3   = x3[i]
x4   = x4[i]
x5   = x5[i]
x6   = x6[i]
x7   = x7[i]
x8   = x8[i]
x9   = x9[i]
x10  = x10[i]
x11  = x11[i]

# cell declustering :
#h_v    = project(CellSize(submesh), Q_b).vector().array()
#A      = sum(h_v)
#wt     = h_v / A
#beta_bar = sum(beta_v * h_v) / A

X    = array([x2,x3,x4,x5,x2*x3,x2*x5,x3*x5])
yt   = sqrt(log(beta_v[valid][i] + 1))
#yt   = beta_v[valid][i]

names = [r'$-T$', r'$\Vert \nabla S \Vert$', 
         r'$|B|$', r'$\Vert \nabla B \Vert$', r'$\Vert \mathbf{U} \Vert$']

fig  = figure()
for i, (x,n) in enumerate(zip(X,names)):
  ax = fig.add_subplot(230 + i + 1)
  j  = argsort(x)
  ax.plot(x[j], beta_v[valid][j], 'ko', alpha=0.1)
  ax.set_title(n)
  ax.set_xlabel(n)
  ax.set_ylabel(r'$\beta$')
  ax.grid()
show()

out  = linRegstats(X, yt, 0.95)

bhat = out['bhat']
yhat = out['yhat']
ciy  = out['CIY']

print "<F_pval, pval>:", out['F_pval'], out['pval']

fig  = figure()
ax   = fig.add_subplot(111)

ax.plot(X[0], yt,     'ko', alpha=0.1)
ax.plot(X[0], yhat,   'r-', lw=2.0)
#ax.plot(u_v, ciy[0], 'k:', lw=2.0)
#ax.plot(u_v, ciy[1], 'k:', lw=2.0)
#ax.set_xlabel(r'$\Vert \mathbf{U} \Vert$')
#ax.set_ylabel(r'$\beta$')
grid()
tight_layout()
show()



