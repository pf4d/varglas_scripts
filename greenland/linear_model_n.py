from pylab             import *
from scipy.stats       import t, distributions, scoreatpercentile, \
                              distributions, probplot, chisquare
from scipy.special     import fdtrc
from scipy.sparse      import diags
from scipy.interpolate import interp1d

import os
import sys

from fenics                    import *
from pylab                     import *
from varglas.helper            import plotIce
from varglas.data.data_factory import DataFactory
from varglas.io                import DataInput

lognorm  = distributions.lognorm

def iqr(arr):
  arr           = sort(arr.copy())
  upperQuartile = scoreatpercentile(arr,.75)
  lowerQuartile = scoreatpercentile(arr,.25)
  iqr           = upperQuartile - lowerQuartile
  return iqr

def normal(x, mu, sigma):
  """ 
  Function which returns the normal value at <x> with mean <mu> and 
  standard deviation <sigma>.
  """
  return 1.0/(sigma * sqrt(2.0 * pi)) * exp(-(x - mu)**2 / (2.0 * sigma**2))

def glm(x,y,w=1.0):

  p,n    = shape(x)                    # sample size
  p     += 1                           # add one for intercept
  dof    = n - p                       # degrees of freedom
  
  sig    = var(y)                      # variance
  mu     = (y + mean(y))/2.0           # initial mean estimate
  eta    = log(mu)                     # initial predictor
  X      = vstack((ones(n), x)).T      # observed x-variable matrix

  # Newton-Raphson :
  converged = False
  rtol      = 1e-15
  dtol      = 1e-15
  lmbda     = 1.0
  nIter     = 0
  deviance  = 1
  D         = 1
  ahat      = zeros(p)   # initial parameters
  rel_res   = zeros(p)   # initial relative residual
  maxIter   = 65

  rel_a = []
  dev_a = []

  while not converged and nIter < maxIter:
    W       = diags(w*mu**2/sig, 0)         # compute weights
    z       = eta + (y - mu)/mu             # adjusted dependent variable

    WX      = W.dot(X)
    XTWX    = dot(X.T, WX)
    iXTWX   = inv(XTWX)
    Wz      = W.dot(z)

    ahat_n  = dot(iXTWX, dot(X.T, Wz))
    
    eta     = dot(X, ahat_n)               # compute estimates
    mu      = exp(eta)                     # linear predictor

    # calculate residual :
    rel_res  = norm(ahat - ahat_n, inf)
    rel_a.append(rel_res)
    ahat     = ahat_n

    D_n      = sum((y - mu)**2)
    deviance = abs(D_n - D)
    D        = D_n
    dev_a.append(deviance)
    
    if rel_res < rtol or deviance < dtol: converged = True
    nIter +=  1

    string = "Newton iteration %d: d (abs) = %.2e, (tol = %.2e) r (rel) = %.2e (tol = %.2e)"
    print string % (nIter, deviance, dtol, rel_res, rtol)
  
  # calculate statistics :
  varA   = diag(iXTWX)            # variance of alpha hat
  sea    = sqrt(varA)             # vector of standard errors for alpha hat
  t_a    = ahat / sea
  pval   = t.sf(abs(t_a), dof) * 2
  conf   = 0.95                      # 95% confidence interval
  tbonf  = t.ppf((1 - conf/p), dof)  # bonferroni corrected t-value
  ci     = tbonf*sea                 # confidence interval for ahat
  resid  = (y - mu)                  # 'working' residual
                                       
  RSS    = sum((y - mu)**2)            # residual sum of squares
  TSS    = sum((y - mean(y))**2)       # total sum of squares
  R2     = (TSS-RSS)/TSS               # R2
  F      = (TSS-RSS)/(p-1) * (n-p)/RSS # F-statistic
  F_p    = fdtrc(p-1, dof, F)          # F-Stat. p-value

  # log-likelihood :
  L      = sum((y*mu - mu**2/2)/(2*sig) - y**2/(2*sig) - 0.5*log(2*pi*sig))
  AIC    = (-2*L + 2*p)/n              # AIC statistic

  # estimated error variance :
  sighat = 1/(n-p) * RSS
                                        
  vara = { 'ahat'  : ahat,              
           'yhat'  : mu,                
           'sea'   : sea,               
           'ci'    : ci,                
           'dof'   : dof,               
           'resid' : resid,             
           'rel_a' : rel_a,
           'dev_a' : dev_a,
           'R2'    : R2,
           'F'     : F,
           'AIC'   : AIC,
           'sighat': sighat}
  return vara


#===============================================================================
# create directories and such :

file_n = sys.argv[1]

if sys.argv[2] == 'weighted':
  file_n = file_n + '/weighted/'

elif sys.argv[2] == 'limited':
  file_n = file_n + '/limited_weighted/'

else:
  file_n = file_n + '/normal/'

print file_n

dirs = ['images/stats/' + file_n, 'dat/' + file_n]

for di in dirs:
  if not os.path.exists(di):
    os.makedirs(di)


#===============================================================================
# get the data from the model output on the bed :

in_dir  = 'dump/bed/09/'
out_dir = 'dump/linear_model/09/'

mesh  = Mesh(in_dir + 'submesh.xdmf')
Q     = FunctionSpace(mesh, 'CG', 1)

S      = Function(Q)
B      = Function(Q)
adot   = Function(Q)
qgeo   = Function(Q)
beta   = Function(Q)
Mb     = Function(Q)
Tb     = Function(Q)
Ts     = Function(Q)
u      = Function(Q)
v      = Function(Q)
w      = Function(Q)
Ubar5  = Function(Q)
Ubar10 = Function(Q)
Ubar20 = Function(Q)
U_ob   = Function(Q)
tau_id = Function(Q)
tau_jd = Function(Q)
tau_ii = Function(Q)
tau_ij = Function(Q)
tau_iz = Function(Q)
tau_ji = Function(Q)
tau_jj = Function(Q)
tau_jz = Function(Q)
mask   = Function(Q)

File(in_dir + 'S_s.xml')    >> S
File(in_dir + 'B_s.xml')    >> B
File(in_dir + 'adot_s.xml') >> adot
File(in_dir + 'qgeo_s.xml') >> qgeo
File(in_dir + 'beta_s.xml') >> beta
File(in_dir + 'Mb_s.xml')   >> Mb
File(in_dir + 'Tb_s.xml')   >> Tb
File(in_dir + 'Ts_s.xml')   >> Ts
File(in_dir + 'ub_s.xml')   >> u
File(in_dir + 'vb_s.xml')   >> v
File(in_dir + 'wb_s.xml')   >> w
File(in_dir + 'bv/Ubar_5.xml')  >> Ubar5
File(in_dir + 'bv/Ubar_10.xml') >> Ubar10
File(in_dir + 'bv/Ubar_20.xml') >> Ubar20
File(in_dir + 'U_ob_s.xml') >> U_ob

File(in_dir + 'tau_id_s.xml') >> tau_id
File(in_dir + 'tau_jd_s.xml') >> tau_jd
File(in_dir + 'tau_ii_s.xml') >> tau_ii
File(in_dir + 'tau_ij_s.xml') >> tau_ij
File(in_dir + 'tau_iz_s.xml') >> tau_iz
File(in_dir + 'tau_ji_s.xml') >> tau_ji
File(in_dir + 'tau_jj_s.xml') >> tau_jj
File(in_dir + 'tau_jz_s.xml') >> tau_jz
File(in_dir + 'mask_s.xml')   >> mask

dSdx   = project(S.dx(0), Q)
dSdy   = project(S.dx(1), Q)
                          
dBdx   = project(B.dx(0), Q)
dBdy   = project(B.dx(1), Q)

# vectors :
beta_v = beta.vector().array()
S_v    = S.vector().array()
B_v    = B.vector().array()
adot_v = adot.vector().array()
qgeo_v = qgeo.vector().array()
Mb_v   = Mb.vector().array()
Tb_v   = Tb.vector().array()
Ts_v   = Ts.vector().array()
u_v    = u.vector().array()
v_v    = v.vector().array()
w_v    = w.vector().array()
Ubar5_v  = Ubar5.vector().array()
Ubar10_v = Ubar10.vector().array()
Ubar20_v = Ubar20.vector().array()
U_ob_v = U_ob.vector().array()
tau_id_v = tau_id.vector().array()
tau_jd_v = tau_jd.vector().array()
tau_ii_v = tau_ii.vector().array()
tau_ij_v = tau_ij.vector().array()
tau_iz_v = tau_iz.vector().array()
tau_ji_v = tau_ji.vector().array()
tau_jj_v = tau_jj.vector().array()
tau_jz_v = tau_jz.vector().array()
mask_v   = mask.vector().array()

H_v    = S_v - B_v
U_mag  = sqrt(u_v**2 + v_v**2 + w_v**2 + 1e-16)
dSdx_v = dSdx.vector().array()
dSdy_v = dSdy.vector().array()
gradS  = sqrt(dSdx_v**2 + dSdy_v**2 + 1e-16)
dBdx_v = dBdx.vector().array()
dBdy_v = dBdy.vector().array()
gradB  = sqrt(dBdx_v**2 + dBdy_v**2 + 1e-16)
D      = zeros(len(B_v))
D[B_v < 0] = B_v[B_v < 0]

Ubar_avg = (Ubar5_v + Ubar10_v + Ubar20_v) / 3.0
ini      = sqrt(917.0 * 9.8 * H_v * gradS / (Ubar_avg + 0.1))

# areas of cells for weighting :
h_v  = project(CellSize(mesh), Q).vector().array()

#===============================================================================
# remove areas with garbage data :
valid  = where(mask_v < 1.0)[0]
valid  = intersect1d(valid, where(beta_v > 1e-14)[0])
valid  = intersect1d(valid, where(beta_v < 1000)[0])
if sys.argv[2] == 'limited':
  valid  = intersect1d(valid, where(U_mag > 20)[0])
else:
  valid  = intersect1d(valid, where(U_mag > 0)[0])
valid  = intersect1d(valid, where(U_ob_v > 1e-9)[0])
valid  = intersect1d(valid, where(abs(Mb_v) < 200)[0])
valid  = intersect1d(valid, where(S_v > -100)[0])
valid  = intersect1d(valid, where(Ts_v > 100)[0])
valid  = intersect1d(valid, where(h_v > 0)[0])
valid  = intersect1d(valid, where(S_v - B_v > 60)[0])
valid  = intersect1d(valid, where(adot_v > -100)[0])

valid_f          = Function(Q)
valid_f_v        = valid_f.vector().array()
valid_f_v[valid] = 1.0
valid_f.vector().set_local(valid_f_v)
valid_f.vector().apply('insert')

rignot   = DataFactory.get_gre_rignot()
drg      = DataInput(rignot, gen_space=False)

betaMax = 200.0

ini_f = Function(Q)
ini_f.vector()[:] = ini
plotIce(drg, ini_f, name='ini', direc='images/stats/' + file_n,
        cmap='gist_yarg', scale='log', numLvls=12, tp=False,
        tpAlpha=0.5, show=False, umin=1.0, umax=betaMax)

Ubar_avg_f = Function(Q)
Ubar_avg_f.vector()[:] = Ubar_avg
plotIce(drg, Ubar_avg_f, name='Ubar_avg', direc='images/stats/' + file_n,
        title=r'$\Vert \mathbf{\bar{u}}_{\bar{bv}} \Vert$', cmap='gist_yarg', 
        scale='log', numLvls=12, tp=False,
        tpAlpha=0.5, show=False, umin=1.0, umax=4000.0)

plotIce(drg, valid_f, name='valid', direc='images/stats/' + file_n,
        cmap='gist_yarg', scale='bool', numLvls=12, tp=False,
        tpAlpha=0.5, show=False)

#===============================================================================
# cell declustering :
n        = len(valid)
h_v      = h_v[valid]
A        = sum(h_v)
wt       = n * h_v / A
beta_bar = 1.0/n * sum(beta_v[valid] * wt)

#===============================================================================

#g_data = [log(g_beta_v + 1), g_S_v, g_B_v, g_H_v, g_adot_v, g_Mb_v, 
#          g_Tb_v, g_Ts_v, log(g_Ubar_v + 1), g_qbar_v, g_u_v, 
#          g_v_v, g_w_v, log(g_U_mag + 1), sqrt(g_gradS + 1), sqrt(g_gradB + 1)]
#data = [beta_v,  S_v,     B_v,    gradS,   gradB, 
#        H_v,     adot_v,  Ts_v,   Tb_v,    Mb_v,
#        Ubar_v,  u_v,     v_v,    w_v,     U_mag]
#names = [r'$\beta$',
#         r'$S$',
#         r'$D$',
#         r'$\Vert \nabla S \Vert$', 
#         r'$\Vert \nabla B \Vert$', 
#         r'$H$',
#         r'$\dot{a}$',
#         r'$T_S$', 
#         r'$T_B$', 
#         r'$M_B$', 
#         r'$\Vert \bar{\mathbf{u}} \Vert$',
#         r'$u$', 
#         r'$v$', 
#         r'$w$', 
#         r'$\Vert \mathbf{u}_B \Vert$']
#
#fig = figure(figsize=(25,15))
#for k,(n,d) in enumerate(zip(names, data)):
#  ax = fig.add_subplot(4,4,k+1)
#  m, bins, pat = hist(d[valid], 1000, normed=1, histtype='stepfilled')
#  setp(pat, 'facecolor', 'b', 'alpha', 0.75)
#  ax.set_xlabel(n)
#  ax.set_ylabel(r'$n$')
#  ax.grid()
#fn = 'images/greenland_data.png'
#savefig(fn, dpi=100)
#show()

#===============================================================================
# data analysis :
y     = beta_v[valid]

#===============================================================================
# do the glm :

v0   = S_v
v1   = Ts_v
v2   = gradS
v3   = D
v4   = gradB
v5   = H_v
v6   = qgeo_v
v7   = adot_v
v8   = Tb_v
v9   = Mb_v
v10  = u_v
v11  = v_v
v12  = w_v
v13  = log(Ubar5_v+1)
v14  = log(Ubar10_v+1)
v15  = log(Ubar20_v+1)
v16  = log(U_mag+1)
v17  = tau_id_v
v18  = tau_jd_v
v19  = tau_ii_v
v20  = tau_ij_v
v21  = tau_iz_v
v22  = tau_ji_v
v23  = tau_jj_v
v24  = tau_jz_v
v25  = ini

x0   = v0[valid]
x1   = v1[valid]
x2   = v2[valid]
x3   = v3[valid]
x4   = v4[valid]
x5   = v5[valid]
x6   = v6[valid]
x7   = v7[valid]
x8   = v8[valid]
x9   = v9[valid]
x10  = v10[valid]
x11  = v11[valid]
x12  = v12[valid]
x13  = v13[valid]
x14  = v14[valid]
x15  = v15[valid]
x16  = v16[valid]
x17  = v17[valid]
x18  = v18[valid]
x19  = v19[valid]
x20  = v20[valid]
x21  = v21[valid]
x22  = v22[valid]
x23  = v23[valid]
x24  = v24[valid]
x25  = v25[valid]

#===============================================================================
# formulte design matrix and do some EDA :
names = [r'$S$', 
         r'$T_S$', 
         r'$\Vert \nabla S \Vert$', 
         r'$D$',
         r'$\Vert \nabla B \Vert$', 
         r'$H$',
         r'$q_{geo}$',
         r'$\dot{a}$',
         r'$T_B$', 
         r'$M_B$', 
         r'$u$', 
         r'$v$', 
         r'$w$', 
         r'$\ln\left(\Vert \bar{\mathbf{u}}_{5} \Vert + 1\right)$',
         r'$\ln\left(\Vert \bar{\mathbf{u}}_{10} \Vert + 1\right)$',
         r'$\ln\left(\Vert \bar{\mathbf{u}}_{20} \Vert + 1\right)$',
         r'$\ln\left(\Vert \mathbf{u}_B \Vert + 1\right)$',
         r'$\tau_{id}$',
         r'$\tau_{jd}$',
         r'$\tau_{ii}$',
         r'$\tau_{ij}$',
         r'$\tau_{iz}$',
         r'$\tau_{ji}$',
         r'$\tau_{jj}$',
         r'$\tau_{jz}$',
         r'ini']

X      = [x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16,x17,x18,
          x19,x20,x21,x22,x23,x24,x25]
V      = [v0,v1,v2,v3,v4,v5,v6,v7,v8,v9,v10,v11,v12,v13,v14,v15,v16,v17,v18,
          v19,v20,v21,v22,v23,v24,v25]

# no u,v,w :
if sys.argv[1] == 'full':
  index  = [0,1,2,4,5,7,8,9,13,14,15,16,19,20,25]

# no stresses :
elif sys.argv[1] == 'no_stress':
  index  = [0,1,2,4,5,7,8,9,13,14,15,16,25]

# no basal velocity :
elif sys.argv[1] == 'no_U':
  index  = [0,1,2,4,5,7,8,9,13,14,15,25]

# no basal melt :
elif sys.argv[1] == 'no_Mb':
  index  = [0,1,2,4,5,7,8,13,14,15,25]

# independent only :
elif sys.argv[1] == 'ind_only':
  index  = [0,1,2,4,5,7,13,14,15,25]

ii     = index
ii_int = []
ii_int.extend(ii)

for i,m in enumerate(ii):
  for j,n in enumerate(ii[i:]):
    if not(m==n and (m == 17 or m==18 or m==19 or m==20 or m==22 or m==23)):
      ii_int.append([m,n])

#fig = figure(figsize=(25,15))
Xt   = []
Vt   = []
ex_n = []

for k,i in enumerate(ii_int):
  
  if type(i) == list:
    n = ''
    x = 1.0
    v = 1.0
    for jj,j in enumerate(i):
      x *= X[j]
      v *= V[j]
      n += names[j]
      if jj < len(i) - 1:
        n += r'$ \star $'
  else:
    x = X[i]
    v = V[i]
    n = names[i]
    #ax = fig.add_subplot(3,4,k+1)
    #ax.plot(x, y, 'ko', alpha=0.1)
    #ax.set_xlabel(n)
    #ax.set_ylabel(r'$\beta$')
    #ax.grid()

  ex_n.append(n)
  Xt.append(x)
  Vt.append(v)

ex_n.insert(0, '$\mathbf{1}$')
 
#show()

#===============================================================================
# plot beta distribution and lognorm fit :

ln_fit   = lognorm.fit(y)
g_x      = linspace(y.min(), y.max(), 1000)
ln_freq  = lognorm.pdf(g_x, *ln_fit)

fig      = figure()
ax       = fig.add_subplot(111)

ax.hist(y, 300, histtype='stepfilled', color='k', alpha=0.5, normed=True,
        label=r'$\beta$')
ax.plot(g_x, ln_freq, lw=2.0, color='r', label=r'$\mathrm{LogNorm}$')
ax.set_xlim([0,200])
#ax.set_ylim([0,0.025])
ax.set_xlabel(r'$\beta$')
ax.set_ylabel('Frequency')
ax.legend(loc='upper right')
ax.grid()
tight_layout()
fn = 'images/stats/' + file_n + 'beta_distribution.png'
savefig(fn, dpi=100)
#show()
close(fig)

#===============================================================================
# fit the glm :

if sys.argv[2] == 'weighted' or sys.argv[2] == 'limited':
  out   = glm(array(Xt), y, wt)
else:
  out   = glm(array(Xt), y)
yhat  = out['yhat']
resid = out['resid']
ahat  = out['ahat']
ci    = out['ci']
  
yhat_f = Function(Q)
y_f    = Function(Q)
resi_f = Function(Q)

yhat_f.vector()[valid] = yhat
y_f.vector()[valid] = y
resi_f.vector()[valid] = resid

plotIce(drg, yhat_f, name='GLM_beta', direc='images/stats/' + file_n, 
        title=r'$\hat{\beta}$', cmap='gist_yarg', scale='log', 
        umin=1.0, umax=betaMax, numLvls=12, tp=False, tpAlpha=0.5, show=False)

plotIce(drg, y_f, name='beta', direc='images/stats/' + file_n, 
        title=r'$\beta$', cmap='gist_yarg', scale='log', 
        umin=1.0, umax=betaMax, numLvls=12, tp=False, tpAlpha=0.5, show=False)

plotIce(drg, resi_f, name='GLM_resid', direc='images/stats/' + file_n, 
        title=r'$d$', cmap='RdGy', scale='lin', 
        umin=-50, umax=50, numLvls=13, tp=False, tpAlpha=0.5, show=False)

#===============================================================================
# data analysis :

fig      = figure()
ax       = fig.add_subplot(111)

ax.hist(y,    300, histtype='step', color='k', lw=1.5, alpha=1.0, normed=True,
        label=r'$\beta$')
ax.hist(yhat, 300, histtype='step', color='r', lw=1.5, alpha=1.0, normed=True,
        label=r'$\hat{\beta}$')
ax.set_xlim([0,200])
#ax.set_ylim([0,0.03])
ax.set_xlabel(r'$\hat{\beta}$')
ax.set_ylabel('Frequency')
ax.legend(loc='upper right')
ax.grid()
tight_layout()
fn = 'images/stats/' + file_n + 'GLM_beta_distributions.png'
savefig(fn, dpi=100)
#show()
close(fig)
  
#=============================================================================
# residual plot and normal quantile plot for residuals :
fig = figure(figsize=(12,5))
ax1 = fig.add_subplot(121)
ax2 = fig.add_subplot(122)

# Normal quantile plot of residuals
((osm,osr), (m, b, r)) = probplot(resid)
interp = interp1d(osm, osr)
yl     = interp(-2.5)
yh     = interp(2.5)
ax1.plot(osm, m*osm + b, 'r-', lw=2.0, label=r'LS fit')
ax1.plot(osm, osr,       'k.', alpha=1.0, label='$\mathbf{d}$')
ax1.set_xlabel('Standard Normal Quantiles')
ax1.set_ylabel('Residuals')
#ax1.set_title('Normal Quantile Plot')
ax1.set_xlim([-2.5, 2.5])
ax1.set_ylim([yl,   yh])
ax1.legend(loc='lower right')
ax1.grid()

ax2.plot(yhat, resid, 'k.', alpha=0.10)
ax2.set_xlabel(r'$\hat{\beta}$')
ax2.set_ylabel('Residuals')
#ax2.set_title('Residual Plot')
ax2.set_xlim([0,  150])
ax2.set_ylim([yl, yh])
ax2.grid()

tight_layout()
fn = 'images/stats/' + file_n + 'GLM_resid_NQ.png'
savefig(fn, dpi=100)
#show()
close(fig)


#=============================================================================
# plot newton residuals :
fig = figure()
ax  = fig.add_subplot(111)

ax.plot(out['rel_a'], 'k-', lw=2.0,
        label=r'$\Vert \alpha - \alpha_n \Vert^2$')
ax.plot(out['dev_a'], 'r-', lw=2.0,
        label=r'$\Vert \mathbf{d} - \mathbf{d}_n \Vert^2$')
ax.set_xlabel(r'Iteration')
ax.set_yscale('log')
ax.set_xlim([0, len(out['dev_a'])-1])
ax.grid()
ax.legend()
tight_layout()
fn = 'images/stats/' + file_n + 'GLM_newton_resid.png'
savefig(fn, dpi=100)
#show()
close(fig)

##=============================================================================
## create partial-residual plot :
#fig  = figure(figsize=(25,15))
#s1   = int(ceil(sqrt(len(ii_int))))
#s2   = int(floor(sqrt(len(ii_int))))
#
#for k,i in enumerate(ii_int):
#  
#  if type(i) == list:
#    n = ''
#    x = 1.0
#    v = 1.0
#    for jj,j in enumerate(i):
#      x *= X[j]
#      v *= V[j]
#      n += names[j]
#      if jj < len(i) - 1:
#        n += r'$ \star $'
#  else:
#    x = X[i]
#    v = V[i]
#    n = names[i]
#   
#  alpha = ahat[k+1]
#  eta   = alpha*x
#  
#  ax = fig.add_subplot(s1,s1,k+1)
#  ax.plot(x, resid + eta, 'ko', alpha=0.1)
#  ax.set_xlabel(n)
#  ax.set_ylabel(r'Residual')
#  ax.grid()
#
#tight_layout()
#
#fn = 'images/stats/' + file_n + 'GLM_partial_residual.png'
#savefig(fn, dpi=100)
#show()
#close(fig)

#===============================================================================
# create tables :
n        = len(valid)

mu       = mean(y)                      # mean
med      = median(y)                    # median
sigma    = std(y)                       # standard deviation
fe_iqr   = iqr(y)                       # IQR
v_m_rat  = sigma**2 / mu                # variance-to-mean ratio
stats_y  = [mu, med, sigma**2, fe_iqr, v_m_rat]

mu       = mean(yhat)                  # mean
med      = median(yhat)                # median
sigma    = std(yhat)                   # standard deviation
fe_iqr   = iqr(yhat)                   # IQR
v_m_rat  = sigma**2 / mu               # variance-to-mean ratio
stats_yh = [mu, med, sigma**2, fe_iqr, v_m_rat, 
            out['R2'], out['F'], out['AIC'], out['sighat']]

f = open('dat/' + file_n + 'alpha.dat', 'w')
for n, a, c in zip(ex_n, ahat, ci):
  al = a-c
  ah = a+c
  if sign(al) != sign(ah):
    strng = '\\color{red}%s & \\color{red}%.1e & \\color{red}%.1e & ' + \
            '\\color{red}%.1e \\\\\n'
  else:
    strng = '%s & %.1e & %.1e & %.1e \\\\\n'
  f.write(strng % (n, al, a, ah))
f.write('\n')
f.close()

names = ['$\mu$', 'median', '$\sigma^2$', 'IQR',   '$\sigma^2 / \mu$',
         '$R^2$', 'F',      'AIC',        '$\hat{\sigma}^2$']

f = open('dat/' + file_n + 'stats.dat', 'w')
for n, s_yh in zip(names, stats_yh):
  strng = '%s & %g \\\\\n' % (n, s_yh)
  f.write(strng)
f.write('\n')
f.close()

#===============================================================================
# reduce the model to explanitory variables with meaning :

ex_a   = array(ex_n)
ahat_n = ahat.copy()
ci_n   = ci.copy()
X_i    = array(Xt)

# find out how many to eliminate first:
v = []
for i,(a,c) in enumerate(zip(ahat_n, ci_n)):
  al = a-c
  ah = a+c
  if sign(al) == sign(ah):
    v.append(i)
v = array(v)
exterminated = len(ahat_n) - len(v)
print "eliminated %i fields" % exterminated

while exterminated > 0:

  ex_a  = ex_a[v]
  X_i   = X_i[v[1:]-1]
  
  out_n = glm(X_i, y)
  
  yhat_n  = out_n['yhat']
  resid_n = out_n['resid']
  ahat_n  = out_n['ahat']
  ci_n    = out_n['ci']
    
  yhat_f = Function(Q)
  resi_f = Function(Q)
  
  yhat_f.vector()[valid] = yhat_n
  resi_f.vector()[valid] = resid_n
  
  plotIce(drg, yhat_f, name='GLM_beta_reduced',
          direc='images/stats/' + file_n, title=r'$\hat{\beta}$',
          cmap='gist_yarg', scale='log', umin=1.0, umax=betaMax, 
          numLvls=12, tp=False, tpAlpha=0.5, show=False)
  
  plotIce(drg, resi_f, name='GLM_resid_reduced',
          direc='images/stats/' + file_n, title=r'$d$',
          cmap='RdGy', scale='lin', umin=-50, umax=50,
          numLvls=13, tp=False, tpAlpha=0.5, show=False)
  
  #=============================================================================
  # data analysis :
  
  fig      = figure()
  ax       = fig.add_subplot(111)
  
  ax.hist(y,      300, histtype='step', color='k', lw=1.5, alpha=1.0,
          normed=True, label=r'$\beta$')
  ax.hist(yhat_n, 300, histtype='step', color='r', lw=1.5, alpha=1.0, 
          normed=True, label=r'$\hat{\beta}$')
  ax.set_xlim([0,200])
  #ax.set_ylim([0,0.03])
  ax.set_xlabel(r'$\hat{\beta}$')
  ax.set_ylabel('Frequency')
  ax.legend(loc='upper right')
  ax.grid()
  tight_layout()
  fn = 'images/stats/'+file_n+'GLM_beta_distributions_reduced.png'
  savefig(fn, dpi=100)
  #show()
  close(fig)
    
  #=============================================================================
  # residual plot and normal quantile plot for residuals :
  fig = figure(figsize=(12,5))
  ax1 = fig.add_subplot(121)
  ax2 = fig.add_subplot(122)
  
  #rtol    = 30
  #Xr      = X_i[:, resid_n < rtol]
  #yhat_n  = yhat_n[resid_n < rtol]
  #resid_n = resid_n[resid_n < rtol]
  
  # Normal quantile plot of residuals
  ((osm,osr), (m, b, r)) = probplot(resid_n)
  interp = interp1d(osm, osr)
  yl     = interp(-2.5)
  yh     = interp(2.5)
  ax1.plot(osm, m*osm + b, 'r-', lw=2.0, label=r'LS fit')
  ax1.plot(osm, osr,       'k.', alpha=1.0, label='$\mathbf{d}$')
  ax1.set_xlabel('Standard Normal Quantiles')
  ax1.set_ylabel('Residuals')
  #ax1.set_title('Normal Quantile Plot')
  ax1.set_xlim([-2.5, 2.5])
  ax1.set_ylim([yl,   yh])
  ax1.legend(loc='lower right')
  ax1.grid()
  
  ax2.plot(yhat_n, resid_n, 'k.', alpha=0.10)
  ax2.set_xlabel(r'$\hat{\beta}$')
  ax2.set_ylabel('Residuals')
  #ax2.set_title('Residual Plot')
  ax2.set_xlim([0,  150])
  ax2.set_ylim([yl, yh])
  ax2.grid()
  
  tight_layout()
  fn = 'images/stats/'+file_n+'GLM_resid-NQ_reduced.png'
  savefig(fn, dpi=100)
  #show()
  close(fig)

  #=============================================================================
  # plot newton residuals :
  fig = figure()
  ax  = fig.add_subplot(111)
  
  ax.plot(out_n['rel_a'], 'k-', lw=2.0,
          label=r'$\Vert \alpha - \alpha_n \Vert^2$')
  ax.plot(out_n['dev_a'], 'r-', lw=2.0,
          label=r'$\Vert \mathbf{d} - \mathbf{d}_n \Vert^2$')
  ax.set_xlabel(r'Iteration')
  ax.set_yscale('log')
  ax.set_xlim([0, len(out_n['dev_a'])-1])
  ax.grid()
  ax.legend()
  fn = 'images/stats/' + file_n + 'GLM_newton_resid_reduced.png'
  tight_layout()
  savefig(fn, dpi=100)
  #show()
  close(fig)

  #=============================================================================
  # save tables :
  fn = open('dat/'+file_n+'alpha_reduced.dat', 'w')
  for n, a, c in zip(ex_a, ahat_n, ci_n):
    al = a-c
    ah = a+c
    if sign(al) != sign(ah):
      strng = '\\color{red}%s & \\color{red}%.1e & \\color{red}%.1e & ' + \
              '\\color{red}%.1e \\\\\n'
    else:
      strng = '%s & %.1e & %.1e & %.1e \\\\\n'
    fn.write(strng % (n, al, a, ah))
  fn.write('\n')
  fn.close()
  
  mu       = mean(yhat_n)                  # mean
  med      = median(yhat_n)                # median
  sigma    = std(yhat_n)                   # standard deviation
  fe_iqr   = iqr(yhat_n)                   # IQR
  v_m_rat  = sigma**2 / mu                 # variance-to-mean ratio
  stats_yh = [mu, med, sigma**2, fe_iqr, v_m_rat, 
              out_n['R2'], out_n['F'], out_n['AIC'], out_n['sighat']]

  names = ['$\mu$', 'median', '$\sigma^2$', 'IQR',   '$\sigma^2 / \mu$',
           '$R^2$', 'F',      'AIC',        '$\hat{\sigma}^2$']
  
  fn = open('dat/' + file_n + 'stats_reduced.dat', 'w')
  for n, s_yh in zip(names, stats_yh):
    strng = '%s & %g \\\\\n' % (n, s_yh)
    fn.write(strng)
  fn.write('\n')
  fn.close()

  v = []
  for i,(a,c) in enumerate(zip(ahat_n, ci_n)):
    al = a-c
    ah = a+c
    if sign(al) == sign(ah):
      v.append(i)
  v = array(v)
  exterminated = len(ahat_n) - len(v)
  print "eliminated %i fields" % exterminated
  





