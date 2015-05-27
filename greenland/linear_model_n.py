from pylab          import *
from scipy.stats    import t, distributions, scoreatpercentile, distributions, \
                           probplot, chisquare
from scipy.special  import fdtrc
from scipy.sparse   import diags

import sys

from fenics                    import *
from pylab                     import *
from varglas.helper            import plotIce
from varglas.data.data_factory import DataFactory
from varglas.io                import DataInput

poisson  = distributions.poisson
gamma    = distributions.gamma
chi2     = distributions.chi2
cauchy   = distributions.cauchy
expon    = distributions.expon
lognorm  = distributions.lognorm
invgauss = distributions.invgauss
invgamma = distributions.invgamma
weibull  = distributions.frechet_l
betad    = distributions.beta


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
  
  sig    = var(y)                      # variance
  mu     = (y + mean(y))/2.0           # initial mean estimate
  eta    = log(mu)                     # initial predictor
  X      = vstack((ones(n), x)).T      # observed x-variable matrix

  # Newton-Raphson :
  converged = False
  rtol      = 1e-6
  dtol      = 1e-7
  lmbda     = 1.0
  nIter     = 0
  deviance  = 1
  D         = 1
  ahat      = zeros(p)   # initial parameters
  rel_res   = zeros(p)   # initial relative residual
  maxIter   = 100

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
    ahat     = ahat_n

    D_n      = sum((y - mu)**2)
    deviance = abs(D_n - D)
    D        = D_n
    
    residual = norm(y - mu)**2 / sig
    
    if rel_res < rtol or deviance < dtol: converged = True
    nIter +=  1

    string = "Newton iteration %d: d (abs) = %.2e, (tol = %.2e) r (rel) = %.2e (tol = %.2e)"
    print string % (nIter, deviance, dtol, rel_res, rtol)

  vara = { 'ahat'  : ahat,
           'yhat'  : mu,
           'resid' : y - mu}
  return vara


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

#===============================================================================
# get the data from the model output on the bed :

in_dir  = 'dump/bed/08/'
out_dir = 'dump/linear_model/08/'

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
U_ob  = Function(Q)

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
File(in_dir + 'Ubar_s.xml') >> Ubar
File(in_dir + 'U_ob_s.xml') >> U_ob

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
U_ob_v = U_ob.vector().array()

H_v    = S_v - B_v
U_mag  = sqrt(u_v**2 + v_v**2 + 1e-16)
dSdx_v = dSdx.vector().array()
dSdy_v = dSdy.vector().array()
gradS  = sqrt(dSdx_v**2 + dSdy_v**2 + 1e-16)
dBdx_v = dBdx.vector().array()
dBdy_v = dBdy.vector().array()
gradB  = sqrt(dBdx_v**2 + dBdy_v**2 + 1e-16)

# areas of cells for weighting :
h_v  = project(CellSize(mesh), Q).vector().array()

#===============================================================================
# remove areas with garbage data :
valid  = where(beta_v > 1e-4)[0]
valid  = intersect1d(valid, where(beta_v < 1000)[0])
valid  = intersect1d(valid, where(log(Ubar_v + 1) > 0)[0])
valid  = intersect1d(valid, where(U_mag > 2)[0])
valid  = intersect1d(valid, where(U_ob_v > 1e-2)[0])
valid  = intersect1d(valid, where(abs(Mb_v) < 200)[0])
valid  = intersect1d(valid, where(S_v > -100)[0])
valid  = intersect1d(valid, where(Ts_v > 100)[0])
valid  = intersect1d(valid, where(h_v > 0)[0])
valid  = intersect1d(valid, where(S_v - B_v > 10)[0])

valid_f = Function(Q)
valid_f_v = valid_f.vector().array()
valid_f_v[valid] = 1.0
valid_f.vector().set_local(valid_f_v)
valid_f.vector().apply('insert')

#beta_v[beta_v < 1e-9] = 1e-9
#beta.vector().set_local(beta_v)
#beta.vector().apply('insert')

rignot   = DataFactory.get_gre_rignot()
drg      = DataInput(rignot, gen_space=False)

betaMax = 500.0

#plotIce(drg, valid_f, name='valid', direc='images/', cmap='gist_yarg', 
#        scale='bool', numLvls=12, tp=False, tpAlpha=0.5)
#plotIce(drg, beta,  name='beta', direc='images/', cmap='gist_yarg',
#        scale='log', umin=1.0, umax=betaMax, numLvls=12, tp=False, tpAlpha=0.5)

#===============================================================================
# cell declustering :
h_v      = h_v[valid]
A        = sum(h_v)
wt       = h_v / A
beta_bar = sum(beta_v[valid] * wt)

#===============================================================================

#g_data = [log(g_beta_v + 1), g_S_v, g_B_v, g_H_v, g_adot_v, g_Mb_v, 
#          g_Tb_v, g_Ts_v, log(g_Ubar_v + 1), g_qbar_v, g_u_v, 
#          g_v_v, g_w_v, log(g_U_mag + 1), sqrt(g_gradS + 1), sqrt(g_gradB + 1)]
data = [beta_v,  S_v,     B_v,    gradS,   gradB, 
        H_v,     adot_v,  Ts_v,   Tb_v,    Mb_v,
        Ubar_v,  u_v,     v_v,    w_v,     U_mag]
names = [r'$\beta$',
         r'$S$',
         r'$B$',
         r'$\Vert \nabla S \Vert$', 
         r'$\Vert \nabla B \Vert$', 
         r'$H$',
         r'$\dot{a}$',
         r'$T_S$', 
         r'$T_B$', 
         r'$M_B$', 
         r'$\Vert \bar{\mathbf{U}} \Vert$',
         r'$u$', 
         r'$v$', 
         r'$w$', 
         r'$\Vert \mathbf{U}_B \Vert$']

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
n     = len(valid)

num  = 1000
vals = arange(0, 8)

mu       = mean(y)                      # mean
med      = median(y)                    # median
sigma    = std(y)                       # standard deviation
fe_iqr   = iqr(y)                       # IQR
v_m_rat  = sigma**2 / mu                # variance-to-mean ratio
px       = poisson.pmf(vals, mu)        # compute probabilities
expfreq  = px * len(y)                  # computes expected frequencies

##===============================================================================
## fitting to distributions :
#gam_fit  = gamma.fit(y)
#chi_fit  = chi2.fit(y)
##cau_fit  = cauchy.fit(y)
#exp_fit  = expon.fit(y)
#ln_fit   = lognorm.fit(y)
#ing_fit  = invgauss.fit(y)
##wei_fit  = weibull.fit(y)
#ingam_fit = invgamma.fit(y)
##beta_fit  = betad.fit(y)
#g_x      = linspace(y.min(), y.max(), num)
#
##===============================================================================
## plotting :
#
#fig      = figure()
#ax       = fig.add_subplot(111)
#
#obs_freq, bins, patches = ax.hist(y, num, histtype='stepfilled', normed=True)
#nml_freq    = normal(g_x, mu, sigma)
##cau_freq    = cauchy.pdf(g_x, *cau_fit)
#ln_freq     = lognorm.pdf(g_x, *ln_fit)
#ing_freq    = invgauss.pdf(g_x, *ing_fit)
##wei_freq    = weibull.pdf(g_x, *wei_fit)
#ingam_freq  = invgamma.pdf(g_x, *ingam_fit)
##betad_freq  = betad.pdf(g_x, *beta_fit)
##chi_freq  = chi2.pdf(g_x, *chi_fit)
##exp_freq  = expon.pdf(g_x, *exp_fit)
##gam_freq  = gamma.pdf(g_x, *gam_fit)
##chi2tst   = chisquare(obs_freq, f_exp=exp_freq, ddof=1)[1]
#ax.plot(g_x, nml_freq,    'r-',  lw=2.0,  label=r'$\mathcal{N}$')
##ax.plot(g_x, cau_freq,    'm-',  lw=2.0,  label=r'$\mathrm{Cauchy}$')
#ax.plot(g_x, ln_freq,     'g-',  lw=2.0,  label=r'$\mathrm{LogNorm}$')
#ax.plot(g_x, ing_freq,    'k-',  lw=2.0,  label=r'$\mathrm{InvGauss}$')
##ax.plot(g_x, wei_freq,    'y-',  lw=2.0,  label=r'$\mathrm{Weibull}$')
#ax.plot(g_x, ingam_freq,  'k:',  lw=2.0,  label=r'$\mathrm{InvGamma}$')
##ax.plot(g_x, betad_freq,  'k--', lw=2.0,  label=r'$\mathrm{Beta}$')
##ax.plot(g_x, chi_freq, 'g-', lw=2.0, label=r'$\chi$')
##ax.plot(g_x, exp_freq, 'k-', lw=2.0, label=r'$\exp$')
##ax.plot(g_x, gam_freq, 'y-', lw=2.0, label=r'$\Gamma$')
#ax.set_xlim([0,200])
#ax.set_ylim([0,0.025])
#ax.set_xlabel(r'$\beta$')
#ax.set_ylabel('Frequency')
#ax.set_title(r'Histogram of $\beta$ values')
#ax.grid()
#leg = ax.legend(loc='lower right')
#leg.get_frame().set_alpha(0.5)
#textstr =   '$n = %i$\n' \
#          + '$\mu = %.2f$\n' \
#          + '$\mathrm{median} = %.2f$\n' \
#          + '$\sigma = %.2f$\n' \
#          + '$\sigma^2 = %.2f$\n' \
#          + '$\mathrm{IQR} = %.3f$\n' \
#          + '$\sigma^2 / \mu = %.2f$\n'
##          + '$\chi^2 \mathrm{test} = %.2f$'
#textstr = textstr % (n, mu, med, sigma, sigma**2, fe_iqr, v_m_rat)
#
## these are matplotlib.patch.Patch properties
#props   = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
#
## place a text box in upper left in axes coords
#ax.text(0.65, 0.95, textstr, transform=ax.transAxes, fontsize=14,
#        verticalalignment='top', bbox=props)
#savefig('images/beta_distribution.png', dpi=200)
#show()

#===============================================================================
# do the glm :

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
x11  = log(Ubar_v[valid]+1)
x12  = Tb_v[valid]
x13  = abs(Mb_v[valid])
x14  = log(U_mag[valid])

v0   = S_v
v1   = Ts_v
v2   = gradS
v3   = abs(B_v)
v4   = gradB
v5   = H_v
v6   = u_v
v7   = v_v
v8   = w_v
v9   = qgeo_v
v10  = adot_v
v11  = log(Ubar_v+1)
v12  = Tb_v
v13  = abs(Mb_v)
v14  = log(U_mag)

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
         r'$\ln\left(\Vert \bar{\mathbf{U}} \Vert\right)$',
         r'$T_B$', 
         r'$|M_B|$', 
         r'$\ln\left(\Vert \mathbf{U}_B \Vert\right)$']

X      = [x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14]
V      = [v0,v1,v2,v3,v4,v5,v6,v7,v8,v9,v10,v11,v12,v13,v14]
#index  = [0,1,2,3,4,5,10,11]                # dependent only
index  = [0,1,2,3,4,5,6,7,8,10,11,12,13,14]  # full
#index  = [0,1,2,3,4,5,10,11,12,13,14]
#index  = [1,2,4,10,11,12,13,14]
ii     = index
ii_int = []
ii_int.extend(ii)

for i,m in enumerate(ii):
  for j,n in enumerate(ii[i+1:]):
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

ex_n.insert(0, '$\hat{a}_0$')
 
#show()

#===============================================================================
# fit the glm :

out_glm   = glm(array(Xt), y)#, wt)
glm_yhat  = out_glm['yhat']
resid     = out_glm['resid']
ahat      = out_glm['ahat']
  
ahat_f = Function(Q)

ahat_f.vector()[valid] = glm_yhat

#plotIce(drg, ahat_f, name='GLM_beta', direc='images/', cmap='gist_yarg',
#        scale='log', umin=1.0, umax=betaMax, numLvls=12, tp=False, tpAlpha=0.5)

#===============================================================================
# data analysis :
n        = len(valid)
mu       = mean(glm_yhat)                  # mean
med      = median(glm_yhat)                # median
sigma    = std(glm_yhat)                   # standard deviation
fe_iqr   = iqr(glm_yhat)                   # IQR
v_m_rat  = sigma**2 / mu                   # variance-to-mean ratio

#===============================================================================
# plotting :
fig      = figure()
ax       = fig.add_subplot(111)

ax.hist(y,    800, histtype='step', alpha=1.0, normed=True,
        label=r'$\beta$')
ax.hist(glm_yhat, 800, histtype='step', alpha=1.0, normed=True,
        label=r'$\hat{\mu}$')
ax.set_xlim([0,200])
ax.set_ylim([0,0.025])
ax.set_xlabel(r'$\hat{\mu}$')
ax.set_ylabel('Frequency')
ax.set_title(r'Histogram of $\hat{\mu}$ values')
ax.legend(loc='center right')
ax.grid()
textstr =   '$n = %i$\n' \
          + '$\mu = %.2f$\n' \
          + '$\mathrm{median} = %.2f$\n' \
          + '$\sigma = %.2f$\n' \
          + '$\sigma^2 = %.2f$\n' \
          + '$\mathrm{IQR} = %.3f$\n' \
          + '$\sigma^2 / \mu = %.2f$'
textstr = textstr % (n, mu, med, sigma, sigma**2, fe_iqr, v_m_rat)

# these are matplotlib.patch.Patch properties
props   = dict(boxstyle='round', facecolor='wheat', alpha=0.5)

# place a text box in upper left in axes coords
ax.text(0.65, 0.95, textstr, transform=ax.transAxes, fontsize=14,
        verticalalignment='top', bbox=props)
fn = 'images/GLM_beta_distributions.png'
savefig(fn, dpi=200)
show()
  
#=============================================================================
# residual plot and normal quantile plot for residuals :
fig = figure(figsize=(12,5))
ax1 = fig.add_subplot(121)
ax2 = fig.add_subplot(122)

ax1.plot(glm_yhat, resid, 'ko', alpha=0.1)
ax1.set_xlabel('Predicted Values') 
ax1.set_ylabel('Residuals') 
ax1.set_title('Residual Plot') 
#ax1.set_xlim([4, 7.5])
#ax1.set_ylim([-3,3])
ax1.grid()

# Normal quantile plot of residuals
probplot(resid, plot=ax2)
ax2.set_xlabel('Standard Normal Quantiles') 
ax2.set_ylabel('Residuals') 
ax2.set_title('Normal Quantile Plot')
#ax2.set_xlim([-6,6])
#ax2.set_ylim([-3,3])
ax2.grid()
fn = 'images/GLM_resid-NQ.png'
savefig(fn, dpi=200)
show()

#=============================================================================
# plot y, yhat :
#fig = figure(figsize=(25,15))
#
#for k,i in enumerate(ii):
#  
#  ax = fig.add_subplot(4,3,k+1)
#  
#  if type(i) == list:
#    n = ''
#    x = 1.0
#    for j in i:
#      x *= X[j]
#      n += names[j]
#  else:
#    x = X[i]
#    n = names[i]
#
#  ax.plot(x, y,        'ko', alpha=0.1)
#  ax.plot(x, glm_yhat, 'ro', alpha=0.1)
#  #ax.set_ylim([4,8])
#  ax.set_xlabel(n)
#  ax.set_ylabel(r'$\beta$')
#  ax.grid()
#fn = 'images/glm_result.png'
#savefig(fn, dpi=100)
#show()
#sys.exit(0)

#===============================================================================
# perform multiple linear regression :
y   = log(beta_v[valid])

#for idx in range(1,len(index)+1):
idx = len(index)+1
ii     = index[:idx]
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

#=============================================================================
# perform regression :

out  = linRegstats(array(Xt), y, 0.95)

print "bhat:", out['bhat']

yhat  = out['yhat']
bhat  = out['bhat']
resid = out['resid']
cibl  = out['CIB'][0]
cibh  = out['CIB'][1]

f = open('table.dat', 'w')
for n, a, b in zip(ex_n, ahat, bhat):
  strng = '%s & %.1e & %.1e \\\\\n' % (n, a, b)
  f.write(strng)
f.write('\n')
f.close()


#=============================================================================
# plot the fields :

#ahat_f = Function(Q)
#bhat_f = Function(Q)
#
#ahat_v = ahat[0]*ones(len(beta_v))
#bhat_v = bhat[0]*ones(len(beta_v))
#
#for k,i in enumerate(ii_int):
#  
#  if type(i) == list:
#    v = 1.0
#    for j in i:
#      v *= Vt[j]
#  else:
#    v = Vt[i]
#
#  print ahat[k+1]
#  ahat_v += ahat[k+1] * v
#  bhat_v += bhat[k+1] * v
#
#ahat_v = exp(ahat_v)
#bhat_v = exp(bhat_v)
#
##ahat_v[ahat_v < 1e-9] = 1e-9
##bhat_v[bhat_v < 1e-9] = 1e-9
#
#ahat_f.vector().set_local(ahat_v)
#bhat_f.vector().set_local(bhat_v)
#
#ahat_f.vector().apply('insert')
#bhat_f.vector().apply('insert')

bhat_f = Function(Q)

bhat_f.vector()[valid] = exp(yhat)

plotIce(drg, bhat_f, name='OLM_beta', direc='images/', cmap='gist_yarg',
        scale='log', umin=1.0, umax=betaMax, numLvls=12, tp=False, tpAlpha=0.5)

#=============================================================================
# remove any data above 1000 for plotting purposes :
lt1e3   = where(exp(yhat) < 1e3)[0]
expYhat = exp(yhat[lt1e3])
expY    = exp(y[lt1e3])

#=============================================================================
# plotting :

n        = len(valid)
mu       = mean(yhat)                      # mean
med      = median(yhat)                    # median
sigma    = std(yhat)                       # standard deviation
fe_iqr   = iqr(yhat)                       # IQR
v_m_rat  = sigma**2 / mu                   # variance-to-mean ratio

fig      = figure()
ax       = fig.add_subplot(111)

y_s = [expY, glm_yhat, expYhat]
l_s = [r'$\beta$', r'GLM - $\hat{\mu}$', r'OLM - $\hat{\mu}$']

ax.hist(y_s, 800, histtype='step', alpha=1.0, normed=True,
        label=l_s, stacked=False)
ax.set_xlim([0,200])
ax.set_ylim([0,0.025])
ax.set_xlabel(r'$\hat{\mu}$')
ax.set_ylabel('Frequency')
ax.set_title(r'Histogram of $\hat{\mu}$ values')
ax.legend(loc='center right')
ax.grid()
textstr =   '$n = %i$\n' \
          + '$\mu = %.2f$\n' \
          + '$\mathrm{median} = %.2f$\n' \
          + '$\sigma = %.2f$\n' \
          + '$\sigma^2 = %.2f$\n' \
          + '$\mathrm{IQR} = %.3f$\n' \
          + '$\sigma^2 / \mu = %.2f$'
textstr = textstr % (n, mu, med, sigma, sigma**2, fe_iqr, v_m_rat)

# these are matplotlib.patch.Patch properties
props   = dict(boxstyle='round', facecolor='wheat', alpha=0.5)

# place a text box in upper left in axes coords
ax.text(0.65, 0.95, textstr, transform=ax.transAxes, fontsize=14,
        verticalalignment='top', bbox=props)
fn = 'images/beta_distributions.png'
savefig(fn, dpi=200)
show()


#=============================================================================
# residual plot and normal quantile plot for residuals :
fig = figure(figsize=(12,5))
ax1 = fig.add_subplot(121)
ax2 = fig.add_subplot(122)

ax1.plot(yhat, resid, 'ko', alpha=0.1)
ax1.set_xlabel('Predicted Values') 
ax1.set_ylabel('Residuals') 
ax1.set_title('Residual Plot') 
#ax1.set_xlim([4, 7.5])
#ax1.set_ylim([-3,3])
ax1.grid()

# Normal quantile plot of residuals
probplot(resid, plot=ax2)
ax2.set_xlabel('Standard Normal Quantiles') 
ax2.set_ylabel('Residuals') 
ax2.set_title('Normal Quantile Plot')
#ax2.set_xlim([-6,6])
#ax2.set_ylim([-3,3])
ax2.grid()
fn = 'images/OLM_resid-NQ_' + str(idx) + '.png'
savefig(fn, dpi=200)
show()


##=============================================================================
## plot y, yhat :
#fig = figure(figsize=(25,25))
#
#for k,i in enumerate(ii):
#  
#  ax = fig.add_subplot(4,4,k+1)
#  
#  if type(i) == list:
#    n = ''
#    x = 1.0
#    for j in i:
#      x *= X[j][lt1e3]
#      n += names[j]
#  else:
#    x = X[i][lt1e3]
#    n = names[i]
#
#  ax.plot(x, expY,    'ko', alpha=0.1)
#  ax.plot(x, expYhat, 'ro', alpha=0.1)
#  #ax.set_ylim([4,8])
#  ax.set_xlabel(n)
#  ax.set_ylabel(r'$\beta$')
#  ax.grid()
#fn = 'images/olm_approx_' + str(idx) + '.png'
#savefig(fn, dpi=100)
#show()
#
##=============================================================================
## plot parameter values with confidence intervals:
#fig  = figure()
#ax   = fig.add_subplot(111)
#
#xb   = range(len(ii_int) + 1)
#
#ax.plot(xb, cibh, 'r--', lw=2.0)
#ax.plot(xb, bhat, 'k-',  lw=2.0)
#ax.plot(xb, cibl, 'r--', lw=2.0)
#ax.set_ylabel(r'$\hat{\beta}_i$')
#ax.set_xlabel(r'$i$')
#grid()
#tight_layout()
#show()



