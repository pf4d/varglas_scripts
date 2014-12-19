import sys
src_directory = '../../statistical_modeling'
sys.path.append(src_directory)

from scipy.stats               import probplot 
from src.regstats              import *
from fenics                    import *
from pylab                     import *
from varglas.mesh.mesh_factory import MeshFactory
from varglas.data.data_factory import DataFactory
from varglas.utilities         import DataInput, DataOutput


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

thklim   = 1.0

measures = DataFactory.get_ant_measures(res=900)
bedmap1  = DataFactory.get_bedmap1(thklim=thklim)
bedmap2  = DataFactory.get_bedmap2(thklim=thklim)

dm       = DataInput(measures, mesh=submesh)
db1      = DataInput(bedmap1,  mesh=submesh)
db2      = DataInput(bedmap2,  mesh=submesh)
    
db2.data['B'] = db2.data['S'] - db2.data['H']
db2.set_data_val('H', 32767, thklim)
db2.data['S'] = db2.data['B'] + db2.data['H']

S     = db2.get_projection("S",     near=True)
B     = db2.get_projection("B",     near=True)
M     = db2.get_projection("mask",  near=True)
Ts    = db1.get_projection("temp",  near=True)
q_geo = db1.get_projection("ghfsr", near=True)
adot  = db1.get_projection("acca",  near=True)
U_ob  = dm.get_projection("U_ob",   near=True)

Q     = FunctionSpace(submesh, 'CG', 1)
beta  = Function(Q)
Tb    = Function(Q)
u     = Function(Q)
v     = Function(Q)
w     = Function(Q)
Mb    = Function(Q)
Ubar  = Function(Q)

File('test/bed/beta_s.xml') >> beta
File('test/bed/Mb_s.xml')   >> Mb
File('test/bed/T_s.xml')    >> Tb
File('test/bed/u_s.xml')    >> u
File('test/bed/v_s.xml')    >> v
File('test/bed/w_s.xml')    >> w
File('test/bed/Ubar_s.xml') >> Ubar

dSdx   = project(S.dx(0), Q)
dSdy   = project(S.dx(1), Q)

dBdx   = project(B.dx(0), Q)
dBdy   = project(B.dx(1), Q)

# vectors :
beta_v = beta.vector().array()
S_v    = S.vector().array()
B_v    = B.vector().array()
H_v    = S_v - B_v
Tb_v   = Tb.vector().array()
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
Ubar_v = Ubar.vector().array()
M_v    = M.vector().array()
qgeo_v = q_geo.vector().array()
adot_v = adot.vector().array()
U_ob_v = U_ob.vector().array()
Ts_v   = Ts.vector().array()


valid  = where(beta_v > 1e-7)[0]
valid  = intersect1d(valid, where(U_mag < 800)[0])
valid  = intersect1d(valid, where(adot_v > 0.0)[0])
valid  = intersect1d(valid, where(log(1 + U_ob_v) > 1e-3)[0])

x1   = Mb_v[valid]
x2   = S_v[valid]
x3   = 273.15 - Tb_v[valid]
x4   = 273.15 - Ts_v[valid]
x5   = gradS[valid]
x6   = abs(B_v[valid])
x7   = gradB[valid]
x8   = H_v[valid]
x9   = u_v[valid]
x10  = v_v[valid]
x11  = w_v[valid]
x12  = U_mag[valid]
x13  = abs(Ubar_v[valid])
x14  = U_ob_v[valid]
x15  = qgeo_v[valid]
x16  = adot_v[valid]

#===============================================================================
# formulte design matrix and do some EDA :
names = [r'$M_b$', 
         r'$S$', 
         r'$-T_B$', 
         r'$-T_S$', 
         r'$\Vert \nabla S \Vert$', 
         r'$|B|$',
         r'$\Vert \nabla B \Vert$', 
         r'$H$',
         r'$u$', 
         r'$v$', 
         r'$w$', 
         r'$\Vert \mathbf{U} \Vert$',
         r'$\Vert \bar{\mathbf{U}} \Vert$',
         r'$\Vert \mathbf{U}_{ob} \Vert$',
         r'$Q_{geo}$',
         r'$\dot{a}$']

X   = [x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16]
y   = log(beta_v[valid] + 1)

ii     = [2,3,4,5,6,11,13,14,15]
ii_int = []
ii_int.extend(ii)

for i,m in enumerate(ii):
  for j,n in enumerate(ii[i+1:]):
    ii_int.append([m,n])

fig = figure()
Xt  = []

for k,i in enumerate(ii_int):
  
  ax = fig.add_subplot(8,7,k+1)
  
  if type(i) == list:
    n = ''
    x = 1.0
    for j in i:
      x *= log(X[j])
      n += names[j]
  else:
    x = log(X[i] + 1)
    n = names[i]

  Xt.append(x)
  
  ax.plot(x, y, 'ko', alpha=0.1)
  ax.set_xlabel(n)
  ax.set_ylabel(r'$\beta$')
  ax.grid()
show()


#===============================================================================
# perform regression :

# cell declustering :
#h_v    = project(CellSize(submesh), Q).vector().array()
#A      = sum(h_v)
#wt     = h_v / A
#beta_bar = sum(beta_v * h_v) / A

out  = linRegstats(array(Xt), y, 0.95)

yhat = out['yhat']
bhat = out['bhat']
cibl = out['CIB'][0]
cibh = out['CIB'][1]

print "<F_pval, pval>:", out['F_pval'], out['pval']


# ==============================================================================
# residual plot and normal quantile plot for residuals :
fig = figure(figsize=(12,5))
ax1 = fig.add_subplot(121)
ax2 = fig.add_subplot(122)

ax1.plot(out['yhat'], out['resid'], 'ko', alpha=0.1)
ax1.set_xlabel('Predicted Values') 
ax1.set_ylabel('Residuals') 
ax1.set_title('Residual Plot') 
ax1.grid()

# Normal quantile plot of residuals
p = prbplotObj(ax2)
probplot(out['resid'], plot=p)
ax2.set_xlabel('Standard Normal Quantiles') 
ax2.set_ylabel('Residuals') 
ax2.set_title('Normal Quantile Plot')
ax2.grid()
#savefig('images/resid-NQ.png', dpi=300)
show()


#===============================================================================
# plot y, yhat :
fig = figure()

for k,i in enumerate(ii):
  
  ax = fig.add_subplot(3,4,k+1)
  
  if type(i) == list:
    n = ''
    x = 1.0
    for j in i:
      x *= log(X[j])
      n += names[j]
  else:
    x = log(X[i] + 1)
    n = names[i]

  ax.plot(x, y,    'ko', alpha=0.1)
  ax.plot(x, yhat, 'ro', alpha=0.1)
  ax.set_xlabel(n)
  ax.set_ylabel(r'$\beta$')
  ax.grid()
show()

#===============================================================================
# plot parameter values with confidence intervals:
fig  = figure()
ax   = fig.add_subplot(111)

xb   = range(len(ii_int) + 1)

ax.plot(xb, cibh, 'r--', lw=2.0)
ax.plot(xb, bhat, 'k-',  lw=2.0)
ax.plot(xb, cibl, 'r--', lw=2.0)
ax.set_ylabel(r'$\hat{\beta}_i$')
ax.set_xlabel(r'$i$')
grid()
tight_layout()
show()



