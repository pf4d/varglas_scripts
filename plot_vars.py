from pylab          import *
from scipy.stats    import t, distributions, scoreatpercentile, distributions, \
                           probplot, chisquare
from scipy.special  import fdtrc
from scipy.sparse   import diags

from varglas.data.data_factory    import DataFactory
from varglas.io                   import DataInput
from varglas.helper               import plotIce
#from plot           import plotIce

import sys
src_directory = '../statistical_modeling'
sys.path.append(src_directory)

from src.regstats              import prbplotObj
from fenics                    import *
from pylab                     import *

poisson = distributions.poisson
gamma   = distributions.gamma
chi2    = distributions.chi2
cauchy  = distributions.cauchy
expon   = distributions.expon

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
  mu     = y                           # initial estimate of mu is data
  X      = vstack((ones(n), x)).T      # observed x-variable matrix
  W      = diags(w*mu**2/sig**2, 0)    # sparse diagonal weight matrix
  
  XTW    = dot(X.T, W)
  XTWX   = dot(XTW, X)
  iXTWX  = inv(XTWX)
  eta    = log(mu)
  z      = eta - 1.0
  
  bhat   = dot(dot(iXTWX, XTW), z)   # parameter initial estimates

  # Newton-Raphson :
  converged = False
  rtol      = 1e-10
  lmbda     = 1.0
  nIter     = 0
  residual  = 1
  rel_res   = bhat
  maxIter   = 13

  while not converged and nIter < maxIter:
    nIter +=  1
    W      = diags(mu**2/sig**2, 0)       # compute weights
    
    XTW    = dot(X.T, W)
    XTWX   = dot(XTW, X)
    iXTWX  = inv(XTWX)
    eta    = dot(X, bhat)                 # compute estimates
    mu     = exp(eta)                     # linear predictor
    z      = eta + y/mu - 1.0

    bhat_n = dot(dot(iXTWX, XTW), z)

    # calculate residual :
    rel_res = norm(bhat - bhat_n, inf)
    bhat    = bhat_n
    if rel_res < rtol: converged = True
    residual = norm(y - mu)**2 / sig**2

    string = "Newton iteration %d: r (abs) = %.3e, r (rel) = %.3e (tol = %.3e)"
    print string % (nIter, residual, rel_res, rtol)

  resid  = y - mu
  vara = { 'bhat'  : bhat,
           'yhat'  : mu,
           'resid' : resid}
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

ant_in_dir = 'antarctica/dump/bed/03/'
gre_in_dir = 'greenland/dump/bed/01/'

a_mesh  = Mesh(ant_in_dir + 'submesh.xdmf')
a_Q     = FunctionSpace(a_mesh, 'CG', 1)

a_S     = Function(a_Q)
a_B     = Function(a_Q)
a_adot  = Function(a_Q)
a_qgeo  = Function(a_Q)
a_beta  = Function(a_Q)
a_Mb    = Function(a_Q)
a_Tb    = Function(a_Q)
a_Ts    = Function(a_Q)
a_u     = Function(a_Q)
a_v     = Function(a_Q)
a_w     = Function(a_Q)
a_Ubar  = Function(a_Q)
a_qbar  = Function(a_Q)

File(ant_in_dir + 'S_s.xml')    >> a_S
File(ant_in_dir + 'B_s.xml')    >> a_B
File(ant_in_dir + 'adot_s.xml') >> a_adot
File(ant_in_dir + 'qgeo_s.xml') >> a_qgeo
File(ant_in_dir + 'beta_s.xml') >> a_beta
File(ant_in_dir + 'Mb_s.xml')   >> a_Mb
File(ant_in_dir + 'Tb_s.xml')   >> a_Tb
File(ant_in_dir + 'Ts_s.xml')   >> a_Ts
File(ant_in_dir + 'u_s.xml')    >> a_u
File(ant_in_dir + 'v_s.xml')    >> a_v
File(ant_in_dir + 'w_s.xml')    >> a_w
File('antarctica/dump/bed/balance_velocity/Ubar_s.xml') >> a_Ubar
File('antarctica/dump/bed/balance_water/q.xml')         >> a_qbar

#===============================================================================
# plotting :
bedmap1  = DataFactory.get_bedmap1(thklim=1.0)
d1       = DataInput(bedmap1, mesh=a_mesh)
  
cmap = 'gist_earth'
cmap = 'RdYlGn'
cmap = 'RdGy'
cmap = 'gist_gray'
cmap = 'Purples'
cmap = 'Reds'
cmap = 'Oranges'

plotIce(d1, 'B', '.', 'gist_yarg', scale='lin', name=r'$B$', 
        numLvls=12, tp=False, tpAlpha=0.5)

#===============================================================================

a_dSdx   = project(a_S.dx(0), a_Q)
a_dSdy   = project(a_S.dx(1), a_Q)
                                
a_dBdx   = project(a_B.dx(0), a_Q)
a_dBdy   = project(a_B.dx(1), a_Q)

# vectors :
a_beta_v = a_beta.vector().array()
a_S_v    = a_S.vector().array()
a_B_v    = a_B.vector().array()
a_adot_v = a_adot.vector().array()
a_qgeo_v = a_qgeo.vector().array()
a_Mb_v   = a_Mb.vector().array()
a_Tb_v   = a_Tb.vector().array()
a_Ts_v   = a_Ts.vector().array()
a_u_v    = a_u.vector().array()
a_v_v    = a_v.vector().array()
a_w_v    = a_w.vector().array()
a_Ubar_v = a_Ubar.vector().array()
a_qbar_v = a_qbar.vector().array()

a_H_v    = a_S_v - a_B_v
a_U_mag  = sqrt(a_u_v**2 + a_v_v**2 + 1e-16)
a_dSdx_v = a_dSdx.vector().array()
a_dSdy_v = a_dSdy.vector().array()
a_gradS  = sqrt(a_dSdx_v**2 + a_dSdy_v**2 + 1e-16)
a_dBdx_v = a_dBdx.vector().array()
a_dBdy_v = a_dBdy.vector().array()
a_gradB  = sqrt(a_dBdx_v**2 + a_dBdy_v**2 + 1e-16)

g_mesh  = Mesh(gre_in_dir + 'submesh.xdmf')
g_Q     = FunctionSpace(g_mesh, 'CG', 1)

g_S     = Function(g_Q)
g_B     = Function(g_Q)
g_adot  = Function(g_Q)
g_qgeo  = Function(g_Q)
g_beta  = Function(g_Q)
g_Mb    = Function(g_Q)
g_Tb    = Function(g_Q)
g_Ts    = Function(g_Q)
g_u     = Function(g_Q)
g_v     = Function(g_Q)
g_w     = Function(g_Q)
g_Ubar  = Function(g_Q)
g_qbar  = Function(g_Q)

File(gre_in_dir + 'S_s.xml')    >> g_S
File(gre_in_dir + 'B_s.xml')    >> g_B
File(gre_in_dir + 'adot_s.xml') >> g_adot
File(gre_in_dir + 'qgeo_s.xml') >> g_qgeo
File(gre_in_dir + 'beta_s.xml') >> g_beta
File(gre_in_dir + 'Mb_s.xml')   >> g_Mb
File(gre_in_dir + 'Tb_s.xml')   >> g_Tb
File(gre_in_dir + 'Ts_s.xml')   >> g_Ts
File(gre_in_dir + 'u_s.xml')    >> g_u
File(gre_in_dir + 'v_s.xml')    >> g_v
File(gre_in_dir + 'w_s.xml')    >> g_w
File('greenland/dump/bed/balance_velocity/Ubar_s.xml') >> g_Ubar
File('greenland/dump/bed/balance_water/q.xml')         >> g_qbar

#===============================================================================
# plotting :
rignot   = DataFactory.get_gre_rignot()
bamber   = DataFactory.get_bamber()
drg      = DataInput(rignot, mesh=g_mesh)
dbm      = DataInput(bamber, mesh=g_mesh)

plotIce(drg, 'U_ob', '.', 'gist_yarg', scale='log', name=r'$B$', 
        numLvls=12, tp=False, tpAlpha=0.5)

sys.exit(0)

#===============================================================================

g_dSdx   = project(g_S.dx(0), g_Q)
g_dSdy   = project(g_S.dx(1), g_Q)
                                
g_dBdx   = project(g_B.dx(0), g_Q)
g_dBdy   = project(g_B.dx(1), g_Q)

# vectors :
g_beta_v = g_beta.vector().array()
g_S_v    = g_S.vector().array()
g_B_v    = g_B.vector().array()
g_adot_v = g_adot.vector().array()
g_qgeo_v = g_qgeo.vector().array()
g_Mb_v   = g_Mb.vector().array()
g_Tb_v   = g_Tb.vector().array()
g_Ts_v   = g_Ts.vector().array()
g_u_v    = g_u.vector().array()
g_v_v    = g_v.vector().array()
g_w_v    = g_w.vector().array()
g_Ubar_v = g_Ubar.vector().array()
g_qbar_v = g_qbar.vector().array()

g_H_v    = g_S_v - g_B_v
g_U_mag  = sqrt(g_u_v**2 + g_v_v**2 + 1e-16)
g_dSdx_v = g_dSdx.vector().array()
g_dSdy_v = g_dSdy.vector().array()
g_gradS  = sqrt(g_dSdx_v**2 + g_dSdy_v**2 + 1e-16)
g_dBdx_v = g_dBdx.vector().array()
g_dBdy_v = g_dBdy.vector().array()
g_gradB  = sqrt(g_dBdx_v**2 + g_dBdy_v**2 + 1e-16)

g_valid  = where(g_beta_v > 1e-10)[0]
g_valid  = intersect1d(g_valid, where(g_Ubar_v > 0)[0])
g_valid  = intersect1d(g_valid, where(g_Ubar_v < 1000)[0])
g_valid  = intersect1d(g_valid, where(abs(g_Mb_v) < 0.1)[0])
g_valid  = intersect1d(g_valid, where(g_S_v >= 0)[0])
g_valid  = intersect1d(g_valid, where(g_H_v >  1)[0])
g_valid  = intersect1d(g_valid, where(g_B_v > -1000)[0])
g_valid  = intersect1d(g_valid, where(abs(g_qbar_v) < 2000)[0])
g_valid  = intersect1d(g_valid, where(abs(g_w_v) < 50)[0])
g_valid  = intersect1d(g_valid, where(abs(g_v_v) < 500)[0])
g_valid  = intersect1d(g_valid, where(abs(g_u_v) < 500)[0])
g_valid  = intersect1d(g_valid, where(g_gradS < 0.5)[0])
g_valid  = intersect1d(g_valid, where(g_gradB < 0.5)[0])

beta_v = hstack((a_beta_v, g_beta_v))
S_v    = hstack((a_S_v,    g_S_v))
Ts_v   = hstack((g_Ts_v,   a_Ts_v))
gradS  = hstack((a_gradS,  g_gradS))
B_v    = hstack((a_B_v,    g_B_v))
gradB  = hstack((a_gradB,  g_gradB))
H_v    = hstack((a_H_v,    g_H_v))
u_v    = hstack((a_u_v,    g_u_v))
v_v    = hstack((a_v_v,    g_v_v))
w_v    = hstack((a_w_v,    g_w_v))
qgeo_v = hstack((a_qgeo_v, g_qgeo_v))
adot_v = hstack((a_adot_v, g_adot_v))
Ubar_v = hstack((a_Ubar_v, g_Ubar_v))
Tb_v   = hstack((a_Tb_v,   g_Tb_v))
Mb_v   = hstack((a_Mb_v,   g_Mb_v))
U_mag  = hstack((a_U_mag,  g_U_mag))

# areas of cells for weighting :
a_h_v  = project(CellSize(a_mesh), a_Q).vector().array()
g_h_v  = project(CellSize(g_mesh), g_Q).vector().array()
h_v    = hstack((a_h_v, g_h_v))

#===============================================================================
# remove areas with garbage data :
valid  = where(beta_v > 1e-9)[0]
valid  = intersect1d(valid, where(beta_v < 1e3)[0])
valid  = intersect1d(valid, where(Ubar_v > 0)[0])
valid  = intersect1d(valid, where(U_mag < 25)[0])
valid  = intersect1d(valid, where(abs(Mb_v) < 200)[0])
valid  = intersect1d(valid, where(S_v > -100)[0])
valid  = intersect1d(valid, where(Ts_v > 100)[0])
valid  = intersect1d(valid, where(h_v > 0)[0])

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
g_data = [g_beta_v, g_S_v, g_B_v, g_H_v, g_adot_v, g_Mb_v, 
          g_Tb_v, g_Ts_v, g_Ubar_v, g_qbar_v, g_u_v, 
          g_v_v, g_w_v, g_U_mag, g_gradS, g_gradB]
a_data = [a_beta_v, a_S_v, a_B_v, a_H_v, a_adot_v, a_Mb_v,
          a_Tb_v, a_Ts_v, a_Ubar_v, a_qbar_v, a_u_v,
          a_v_v, a_w_v, a_U_mag, a_gradS, a_gradB]
names = [r'$\beta$',
         r'$S$',
         r'$B$',
         r'$H$',
         r'$\dot{a}$',
         r'$M_B$', 
         r'$T_B$', 
         r'$T_S$', 
         r'$\Vert \bar{\mathbf{U}} \Vert$',
         r'$\Vert \bar{\mathbf{q}} \Vert$',
         r'$u$', 
         r'$v$', 
         r'$w$', 
         r'$\Vert \mathbf{U}_B \Vert$',
         r'$\Vert \nabla S \Vert$', 
         r'$\Vert \nabla B \Vert$'] 

#fig = figure(figsize=(25,15))
#for k,(n,d) in enumerate(zip(names, g_data)):
#  ax = fig.add_subplot(4,4,k+1)
#  m, bins, pat = hist(d[g_valid], 1000, normed=1, histtype='stepfilled')
#  setp(pat, 'facecolor', 'b', 'alpha', 0.75)
#  ax.set_xlabel(n)
#  ax.set_ylabel(r'$n$')
#  ax.grid()
#fn = 'images/linear_model/combined/greenland_data.png'
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

#valid = intersect1d(valid, where(y > 0)[0])
#y  = y[valid]

#===============================================================================
# fitting to distributions :
gam_fit  = gamma.fit(y)
chi_fit  = chi2.fit(y)
cau_fit  = cauchy.fit(y)
exp_fit  = expon.fit(y)
g_x      = linspace(y.min(), y.max(), num)

#===============================================================================
# plotting :

fig      = figure()
ax       = fig.add_subplot(111)

obs_freq, bins, patches = ax.hist(y, num, histtype='stepfilled', normed=True)
nml_freq  = normal(g_x, mu, sigma)
chi_freq  = chi2.pdf(g_x, chi_fit[0], chi_fit[1], chi_fit[2])
cau_freq  = cauchy.pdf(g_x, cau_fit[0], cau_fit[1])
exp_freq  = expon.pdf(g_x, exp_fit[0], exp_fit[1])
gam_freq  = gamma.pdf(g_x, gam_fit[0], gam_fit[1], gam_fit[2])
#chi2tst   = chisquare(obs_freq, f_exp=exp_freq, ddof=1)[1]
ax.plot(g_x, nml_freq, 'r-', lw=2.0, label=r'$\mathcal{N}$')
ax.plot(g_x, chi_freq, 'g-', lw=2.0, label=r'$\chi$')
ax.plot(g_x, cau_freq, 'm-', lw=2.0, label=r'$\mathrm{Cauchy}$')
ax.plot(g_x, exp_freq, 'k-', lw=2.0, label=r'$\exp$')
ax.plot(g_x, gam_freq, 'y-', lw=2.0, label=r'$\Gamma$')
ax.set_xlim([0,400])
ax.set_ylim([0,0.02])
ax.set_xlabel(r'$\beta$')
ax.set_ylabel('Frequency')
ax.set_title(r'Histogram of $\beta$ values')
ax.grid()
leg = ax.legend(loc='lower right')
leg.get_frame().set_alpha(0.5)
textstr =   '$n = %i$\n' \
          + '$\mu = %.2f$\n' \
          + '$\mathrm{median} = %.2f$\n' \
          + '$\sigma = %.2f$\n' \
          + '$\sigma^2 = %.2f$\n' \
          + '$\mathrm{IQR} = %.3f$\n' \
          + '$\sigma^2 / \mu = %.2f$\n'
#          + '$\chi^2 \mathrm{test} = %.2f$'
textstr = textstr % (n, mu, med, sigma, sigma**2, fe_iqr, v_m_rat)

# these are matplotlib.patch.Patch properties
props   = dict(boxstyle='round', facecolor='wheat', alpha=0.5)

# place a text box in upper left in axes coords
ax.text(0.65, 0.95, textstr, transform=ax.transAxes, fontsize=14,
        verticalalignment='top', bbox=props)
show()

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
         r'$\ln\left(\Vert \bar{\mathbf{U}} + 1 \Vert\right)$',
         r'$T_B$', 
         r'$|M_B|$', 
         r'$\ln\left(\Vert \mathbf{U}_B + 1 \Vert\right)$']

X      = [x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14]
index  = [0,1,2,3,4,5,10,11,12,13,14]
ii     = index
ii_int = []
ii_int.extend(ii)

for i,m in enumerate(ii):
  for j,n in enumerate(ii[i+1:]):
    ii_int.append([m,n])

#fig = figure(figsize=(25,15))
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
    #n = names[i]
    #ax = fig.add_subplot(3,4,k+1)
    #ax.plot(x, y, 'ko', alpha=0.1)
    #ax.set_xlabel(n)
    #ax.set_ylabel(r'$\beta$')
    #ax.grid()

  Xt.append(x)
 
#show()

#===============================================================================
# fit the glm :

out_glm   = glm(array(Xt), y)
yhat      = out_glm['yhat']
resid     = out_glm['resid']
bhat_glm  = out_glm['bhat']

# remove any data above 1000 for plotting purposes :
#lt1e3     = where(yhat < 1000)[0]
#yhat      = yhat[lt1e3]
#y         = y[y < 1000]
#resid     = resid[lt1e3]

#===============================================================================
# data analysis :
n        = len(valid)
mu       = mean(yhat)                      # mean
med      = median(yhat)                    # median
sigma    = std(yhat)                       # standard deviation
fe_iqr   = iqr(yhat)                       # IQR
v_m_rat  = sigma**2 / mu                   # variance-to-mean ratio

#===============================================================================
# plotting :
fig      = figure()
ax       = fig.add_subplot(111)

ax.hist(y,    800, histtype='stepfilled', alpha=0.5, normed=True)
ax.hist(yhat, 800, histtype='stepfilled', alpha=0.5, normed=True)
ax.set_xlim([0,400])
ax.set_ylim([0,0.02])
ax.set_xlabel(r'$\hat{\beta}$')
ax.set_ylabel('Frequency')
ax.set_title(r'Histogram of $\hat{\beta}$ values')
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
fn = 'images/linear_model/combined/glm_beta_distributions.png'
savefig(fn, dpi=100)
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
fn = 'images/linear_model/combined/resid-glm-NQ.png'
savefig(fn, dpi=100)
show()

#=============================================================================
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
  #ax.set_ylim([4,8])
  ax.set_xlabel(n)
  ax.set_ylabel(r'$\beta$')
  ax.grid()
fn = 'images/linear_model/combined/glm_result.png'
savefig(fn, dpi=100)
show()

#===============================================================================
# perform multiple linear regression :
y   = log(beta_v[valid] + 1)

#for idx in range(1,len(index)+1):
for idx in [len(index)+1]:
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
  
  print "<F_pval, pval>:", out['F_pval'], out['pval']
  
  yhat  = out['yhat']
  bhat  = out['bhat']
  resid = out['resid']
  cibl  = out['CIB'][0]
  cibh  = out['CIB'][1]

  # remove any data above 1000 for plotting purposes :
  lt1e3   = where(exp(yhat)-1 < 1e3)[0]
  expYhat = exp(yhat[lt1e3]) - 1
  expY    = exp(y[lt1e3]) - 1

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
  
  ax.hist(expY,    800, histtype='stepfilled', alpha=0.5, normed=True)
  ax.hist(expYhat, 800, histtype='stepfilled', alpha=0.5, normed=True)
  ax.set_xlim([0,400])
  ax.set_ylim([0,0.02])
  ax.set_xlabel(r'$\hat{\beta}$')
  ax.set_ylabel('Frequency')
  ax.set_title(r'Histogram of $\hat{\beta}$ values')
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
  fn = 'images/linear_model/combined/beta_distributions.png'
  savefig(fn, dpi=100)
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
  fn = 'images/linear_model/combined/resid-NQ_' + str(idx) + '.png'
  savefig(fn, dpi=100)
  show()
  
  
  #=============================================================================
  # plot y, yhat :
  fig = figure(figsize=(25,15))
  
  for k,i in enumerate(ii):
    
    ax = fig.add_subplot(4,3,k+1)
    
    if type(i) == list:
      n = ''
      x = 1.0
      for j in i:
        x *= X[j][lt1e3]
        n += names[j]
    else:
      x = X[i][lt1e3]
      n = names[i]
  
    ax.plot(x, expY,    'ko', alpha=0.1)
    ax.plot(x, expYhat, 'ro', alpha=0.1)
    #ax.set_ylim([4,8])
    ax.set_xlabel(n)
    ax.set_ylabel(r'$\beta$')
    ax.grid()
  fn = 'images/linear_model/combined/approx_' + str(idx) + '.png'
  savefig(fn, dpi=100)
  show()
  
  #=============================================================================
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



