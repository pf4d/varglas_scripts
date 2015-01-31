from pylab          import *
from scipy.stats    import t, distributions
from scipy.special  import fdtrc

chi2cdf = distributions.chi2.cdf

def linRegstats(x, y, conf):

  if size(shape(x)) == 1: 
    p    = 2
    n    = len(x)                 # Sample size
  else:
    p,n  = shape(x)               # (samples size, # of parameters)
    p    = p + 1                  # add one for intercept
  dof    = n - p                  # Degrees of freedom
  nov    = p - 1                  # number of variables
  ym     = mean(y)                # Mean of log recovery
  X      = vstack((ones(n), x)).T # observed x-variable matrix
  XTX    = dot(X.T, X)
  iXTX   = inv(XTX)
  bhat   = dot(dot(iXTX, X.T), y)
  yhat   = dot(X, bhat)           # Linear fit line
  SSE    = sum((y    - yhat)**2)  # Sum of Squared Errors
  SST    = sum((y    - ym  )**2)  # Sum of Squared Total
  SSR    = sum((yhat - ym  )**2)  # Sum of Squared Residuals (SSR = SST - SSE)
  R2     = SSR / SST              # R^2 Statistic (rval**2)
  MSE    = SSE / dof              # Mean Squared Error (MSE)
  MSR    = SSR / nov              # Mean Squared Residual (MSR)
  F      = MSR / MSE              # F-Statistic
  F_p    = fdtrc(nov, dof, F)     # F-Stat. p-value
  
  # variance of beta estimates :
  VARb   = MSE * iXTX             # diag(VARb) == varB
  varB   = diag(VARb)             # variance of beta hat
  seb    = sqrt(varB)             # vector of standard errors for beta hat
 
  ## variance of y estimates : 
  #VARy   = MSE * dot(dot(X, iXTX), X.T)
  #varY   = diag(VARy)             # variance of y hat
  #sey    = sqrt(varY)             # standard errors for yhat
 
  # calculate t-statistic :
  t_b    = bhat / seb
  #t_y    = yhat / sey
 
  # calculate p-values :
  pval   = t.sf(abs(t_b), dof) * 2
  
  tbonf  = t.ppf((1+conf)/2.0, dof)  # uncorrected t*-value
  
  ci_b   = [bhat - tbonf*seb,     # Confidence intervals for betas
            bhat + tbonf*seb]     #   in 2 columns (lower,upper)
  
  #ci_y   = [yhat - tbonf*sey,     # Confidence intervals for estimates
  #          yhat + tbonf*sey]     #   in 2 columns (lower,upper)

  resid  = y - yhat

  #vara = { 'SSR'   : SSR,
  #         'SSE'   : SSE,
  #         'SST'   : SST,
  #         'df'    : (nov, dof, n-1),
  #         'MSR'   : MSR,
  #         'MSE'   : MSE,
  #         'F'     : F,
  #         'F_pval': F_p,
  #         'varY'  : varY,
  #         'varB'  : varB,
  #         'SEB'   : seb,
  #         'SEY'   : sey,
  #         'tbonf' : tbonf,
  #         't_beta': t_b,
  #         't_y'   : t_y,
  #         'pval'  : pval,
  #         'CIB'   : array(ci_b),
  #         'CIY'   : array(ci_y),
  #         'bhat'  : bhat,
  #         'yhat'  : yhat,
  #         'R2'    : R2,
  #         'resid' : resid}
  vara = { 'SSR'   : SSR,
           'SSE'   : SSE,
           'SST'   : SST,
           'df'    : (nov, dof, n-1),
           'MSR'   : MSR,
           'MSE'   : MSE,
           'F'     : F,
           'F_pval': F_p,
           'varB'  : varB,
           'SEB'   : seb,
           'tbonf' : tbonf,
           't_beta': t_b,
           'pval'  : pval,
           'CIB'   : array(ci_b),
           'bhat'  : bhat,
           'yhat'  : yhat,
           'R2'    : R2,
           'resid' : resid}
  return vara


import sys
src_directory = '../../statistical_modeling'
sys.path.append(src_directory)

from scipy.stats               import probplot 
from src.regstats              import prbplotObj
from fenics                    import *
from pylab                     import *

#===============================================================================
# get the data from the model output on the bed :

out_dir  = 'dump/bed/linear_model/'
in_dir   = 'dump/bed/03/'

mesh  = Mesh(in_dir + 'submesh.xdmf')
Q     = FunctionSpace(mesh, 'CG', 1)

S     = Function(Q)
B     = Function(Q)
adot  = Function(Q)
qgeo  = Function(Q)
beta  = Function(Q)
Mb    = Function(Q)
Tb    = Function(Q)
Ts    = Function(Q)
u     = Function(Q)
v     = Function(Q)
w     = Function(Q)
Ubar  = Function(Q)
qbar  = Function(Q)

File(in_dir + 'S_s.xml')    >> S
File(in_dir + 'B_s.xml')    >> B
File(in_dir + 'adot_s.xml') >> adot
File(in_dir + 'qgeo_s.xml') >> qgeo
File(in_dir + 'beta_s.xml') >> beta
File(in_dir + 'Mb_s.xml')   >> Mb
File(in_dir + 'Tb_s.xml')   >> Tb
File(in_dir + 'Ts_s.xml')   >> Ts
File(in_dir + 'u_s.xml')    >> u
File(in_dir + 'v_s.xml')    >> v
File(in_dir + 'w_s.xml')    >> w
File('dump/bed/balance_velocity/Ubar_s.xml') >> Ubar
File('dump/bed/balance_water/q.xml')         >> qbar

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
Ubar_v = Ubar.vector().array()
qbar_v = qbar.vector().array()

H_v    = S_v - B_v
U_mag  = sqrt(u_v**2 + v_v**2 + w_v**2 + 1e-16)
dSdx_v = dSdx.vector().array()
dSdy_v = dSdy.vector().array()
gradS  = sqrt(dSdx_v**2 + dSdy_v**2 + 1e-16)
dBdx_v = dBdx.vector().array()
dBdy_v = dBdy.vector().array()
gradB  = sqrt(dBdx_v**2 + dBdy_v**2 + 1e-16)

#B_v -= B_v.min()

valid  = where(beta_v > 1e-10)[0]
valid  = intersect1d(valid, where(Ubar_v > 0)[0])
valid  = intersect1d(valid, where(abs(Mb_v) < 40)[0])
valid  = intersect1d(valid, where(abs(adot_v) < 2)[0])

print "sample size:", len(valid)

x0   = S_v[valid]
x1   = Ts_v[valid]
x2   = gradS[valid]
x3   = abs(B_v[valid])
x4   = gradB[valid]
x5   = H_v[valid]
x6   = u_v[valid]
x7   = v_v[valid]
x8   = w_v[valid]
x9   = qgeo_v[valid]
x10  = adot_v[valid]
x11  = log(Ubar_v[valid] + 1)
x12  = Tb_v[valid]
x13  = abs(Mb_v[valid])
x14  = log(U_mag[valid] + 1)

#===============================================================================
# formulte design matrix and do some EDA :
names = [r'$S$', 
         r'$T_S$', 
         r'$\Vert \nabla S \Vert$', 
         r'$|B|$',
         r'$\Vert \nabla B \Vert$', 
         r'$H$',
         r'$u$', 
         r'$v$', 
         r'$w$', 
         r'$Q_{geo}$',
         r'$\dot{a}$',
         r'$\Vert \bar{\mathbf{U}} \Vert$',
         r'$T_B$', 
         r'$|M_B|$', 
         r'$\Vert \mathbf{U}_B \Vert$']

X   = [x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14]
y   = log(beta_v[valid] + 100)

#ii     = [0,1,2,3,4,5,10,11,12,13,14]
ii     = [0,1,2,3,4,5,10,11,12,13,14]
ii_int = []
ii_int.extend(ii)

for i,m in enumerate(ii):
  for j,n in enumerate(ii[i+1:]):
    ii_int.append([m,n])

#fig = figure(figsize=(25,10))
Xt  = []

for k,i in enumerate(ii_int):
  
  if type(i) == list:
    n = ''
    x = 1.0
    for j in i:
      x *= X[j]
      n += names[j]
  else:
    x = X[i]
    n = names[i]
    #ax = fig.add_subplot(3,4,k+1)
    #ax.plot(x, y, 'ko', alpha=0.1)
    #ax.set_xlabel(n)
    #ax.set_ylabel(r'$\beta$')
    #ax.grid()

  Xt.append(x)
  
#show()


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


#===============================================================================
# residual plot and normal quantile plot for residuals :
fig = figure(figsize=(12,5))
ax1 = fig.add_subplot(121)
ax2 = fig.add_subplot(122)

ax1.plot(out['yhat'], out['resid'], 'ko', alpha=0.1)
ax1.set_xlabel('Predicted Values') 
ax1.set_ylabel('Residuals') 
ax1.set_title('Residual Plot') 
#ax1.set_xlim([4, 7.5])
#ax1.set_ylim([-3,3])
ax1.grid()

# Normal quantile plot of residuals
p = prbplotObj(ax2)
probplot(out['resid'], plot=p)
ax2.set_xlabel('Standard Normal Quantiles') 
ax2.set_ylabel('Residuals') 
ax2.set_title('Normal Quantile Plot')
ax2.set_xlim([-6,6])
ax2.set_ylim([-3,3])
ax2.grid()
savefig('../images/linear_model/resid-NQ_' + str(len(ii_int)) + '.png', dpi=100)
#show()


#===============================================================================
# plot y, yhat :
fig = figure(figsize=(25,15))

for k,i in enumerate(ii):
  
  ax = fig.add_subplot(4,3,k+1)
  
  if type(i) == list:
    n = ''
    x = 1.0
    for j in i:
      x *= X[j]
      n += names[j]
  else:
    x = X[i]
    n = names[i]

  ax.plot(x, y,    'ko', alpha=0.1)
  ax.plot(x, yhat, 'ro', alpha=0.1)
  ax.set_ylim([4,8])
  ax.set_xlabel(n)
  ax.set_ylabel(r'$\beta$')
  ax.grid()
savefig('../images/linear_model/approx_' + str(len(ii_int)) + '.png', dpi=100)
#show()

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
#show()



