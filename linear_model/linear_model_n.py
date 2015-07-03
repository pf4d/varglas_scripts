from pylab             import *
from scipy.stats       import t, distributions, scoreatpercentile, \
                              distributions, probplot, chisquare
from scipy.special     import fdtrc
from scipy.sparse      import diags
from scipy.interpolate import interp1d

import sys
import os

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
  rtol      = 1e-4
  dtol      = 1e-4
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
  conf   = 0.95                        # 95% confidence interval
  tbonf  = t.ppf((1 - conf/p), dof)    # bonferroni corrected t-value
  ci     = tbonf*sea                   # confidence interval for ahat
  resid  = (y - mu)                    # 'working' residual
                                       
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

a_fn = 'images/stats/' + file_n + 'antarctica/'
g_fn = 'images/stats/' + file_n + 'greenland/'

dirs = [a_fn, g_fn, 'dat/' + file_n]

for di in dirs:
  if not os.path.exists(di):
    os.makedirs(di)


#===============================================================================
# get the data from the model output on the bed :

a_in_dir  = '../antarctica/dump/bed/07/'
g_in_dir  = '../greenland/dump/bed/09/'

#===============================================================================
# antarctica :

a_mesh  = Mesh(a_in_dir + 'submesh.xdmf')
a_Q     = FunctionSpace(a_mesh, 'CG', 1)

a_S      = Function(a_Q)
a_B      = Function(a_Q)
a_adot   = Function(a_Q)
a_qgeo   = Function(a_Q)
a_beta   = Function(a_Q)
a_Mb     = Function(a_Q)
a_Tb     = Function(a_Q)
a_Ts     = Function(a_Q)
a_u      = Function(a_Q)
a_v      = Function(a_Q)
a_w      = Function(a_Q)
a_Ubar5  = Function(a_Q)
a_Ubar10 = Function(a_Q)
a_Ubar20 = Function(a_Q)
a_U_ob   = Function(a_Q)
a_tau_id = Function(a_Q)
a_tau_jd = Function(a_Q)
a_tau_ii = Function(a_Q)
a_tau_ij = Function(a_Q)
a_tau_iz = Function(a_Q)
a_tau_ji = Function(a_Q)
a_tau_jj = Function(a_Q)
a_tau_jz = Function(a_Q)
a_mask   = Function(a_Q)

File(a_in_dir + 'S_s.xml')    >> a_S
File(a_in_dir + 'B_s.xml')    >> a_B
File(a_in_dir + 'adot_s.xml') >> a_adot
File(a_in_dir + 'qgeo_s.xml') >> a_qgeo
File(a_in_dir + 'beta_s.xml') >> a_beta
File(a_in_dir + 'Mb_s.xml')   >> a_Mb
File(a_in_dir + 'Tb_s.xml')   >> a_Tb
File(a_in_dir + 'Ts_s.xml')   >> a_Ts
File(a_in_dir + 'ub_s.xml')   >> a_u
File(a_in_dir + 'vb_s.xml')   >> a_v
File(a_in_dir + 'wb_s.xml')   >> a_w
File(a_in_dir + 'bv/Ubar_5.xml')  >> a_Ubar5
File(a_in_dir + 'bv/Ubar_10.xml') >> a_Ubar10
File(a_in_dir + 'bv/Ubar_20.xml') >> a_Ubar20
File(a_in_dir + 'U_ob_s.xml') >> a_U_ob

File(a_in_dir + 'tau_id_s.xml') >> a_tau_id
File(a_in_dir + 'tau_jd_s.xml') >> a_tau_jd
File(a_in_dir + 'tau_ii_s.xml') >> a_tau_ii
File(a_in_dir + 'tau_ij_s.xml') >> a_tau_ij
File(a_in_dir + 'tau_iz_s.xml') >> a_tau_iz
File(a_in_dir + 'tau_ji_s.xml') >> a_tau_ji
File(a_in_dir + 'tau_jj_s.xml') >> a_tau_jj
File(a_in_dir + 'tau_jz_s.xml') >> a_tau_jz
File(a_in_dir + 'mask_s.xml')   >> a_mask

a_dSdx   = project(a_S.dx(0), a_Q)
a_dSdy   = project(a_S.dx(1), a_Q)
                          
a_dBdx   = project(a_B.dx(0), a_Q)
a_dBdy   = project(a_B.dx(1), a_Q)
                          
a_dHdx   = project((a_S - a_B).dx(0), a_Q)
a_dHdy   = project((a_S - a_B).dx(1), a_Q)

# vectors :
a_beta_v   = a_beta.vector().array()
a_S_v      = a_S.vector().array()
a_B_v      = a_B.vector().array()
a_adot_v   = a_adot.vector().array()
a_qgeo_v   = a_qgeo.vector().array()
a_Mb_v     = a_Mb.vector().array()
a_Tb_v     = a_Tb.vector().array()
a_Ts_v     = a_Ts.vector().array()
a_u_v      = a_u.vector().array()
a_v_v      = a_v.vector().array()
a_w_v      = a_w.vector().array()
a_Ubar5_v  = a_Ubar5.vector().array()
a_Ubar10_v = a_Ubar10.vector().array()
a_Ubar20_v = a_Ubar20.vector().array()
a_U_ob_v   = a_U_ob.vector().array()
a_tau_id_v = a_tau_id.vector().array()
a_tau_jd_v = a_tau_jd.vector().array()
a_tau_ii_v = a_tau_ii.vector().array()
a_tau_ij_v = a_tau_ij.vector().array()
a_tau_iz_v = a_tau_iz.vector().array()
a_tau_ji_v = a_tau_ji.vector().array()
a_tau_jj_v = a_tau_jj.vector().array()
a_tau_jz_v = a_tau_jz.vector().array()
a_mask_v   = a_mask.vector().array()

a_H_v    = a_S_v - a_B_v
a_dHdx_v = a_dHdx.vector().array()
a_dHdy_v = a_dHdy.vector().array()
a_gradH  = sqrt(a_dHdx_v**2 + a_dHdy_v**2 + 1e-16)
a_U_mag  = sqrt(a_u_v**2 + a_v_v**2 + 1e-16)
a_dSdx_v = a_dSdx.vector().array()
a_dSdy_v = a_dSdy.vector().array()
a_gradS  = sqrt(a_dSdx_v**2 + a_dSdy_v**2 + 1e-16)
a_dBdx_v = a_dBdx.vector().array()
a_dBdy_v = a_dBdy.vector().array()
a_gradB  = sqrt(a_dBdx_v**2 + a_dBdy_v**2 + 1e-16)
a_D      = zeros(len(a_B_v))
a_D[a_B_v < 0] = a_B_v[a_B_v < 0]

a_taux = -917.0 * 9.8 * a_H_v * a_dSdx_v
a_tauy = -917.0 * 9.8 * a_H_v * a_dSdy_v
a_tau_mag = sqrt(a_taux**2 + a_tauy**2 + 1e-16)

if sys.argv[1] == 'Ubar':
  a_uhat = a_taux / a_tau_mag
  a_vhat = a_tauy / a_tau_mag

elif sys.argv[1] == 'U' or sys.argv[1] == 'stress':
  a_uhat = a_u_v / a_U_mag
  a_vhat = a_v_v / a_U_mag

a_dBdi = a_dBdx_v * a_uhat + a_dBdy_v * a_vhat
a_dBdj = a_dBdx_v * a_vhat - a_dBdy_v * a_uhat

a_dSdi = a_dSdx_v * a_uhat + a_dSdy_v * a_vhat
a_dSdj = a_dSdx_v * a_vhat - a_dSdy_v * a_uhat

a_dHdi = a_dHdx_v * a_uhat + a_dHdy_v * a_vhat
a_dHdj = a_dHdx_v * a_vhat - a_dHdy_v * a_uhat

a_Ubar_avg = (a_Ubar5_v + a_Ubar10_v + a_Ubar20_v) / 3.0

if sys.argv[1] == 'Ubar':
  a_ini_i    = 917.0 * 9.8 * a_H_v * a_dSdi / (a_Ubar5_v + 0.1)
  a_ini_j    = 917.0 * 9.8 * a_H_v * a_dSdj / (a_Ubar5_v + 0.1)

elif sys.argv[1] == 'U' or sys.argv[1] == 'stress':
  a_ini_i    = 917.0 * 9.8 * a_H_v * a_dSdi / (a_U_mag + 0.1)
  a_ini_j    = 917.0 * 9.8 * a_H_v * a_dSdj / (a_U_mag + 0.1)

# areas of cells for weighting :
a_h_v  = project(CellSize(a_mesh), a_Q).vector().array()

# number of dofs :
a_n = len(a_beta_v)

#===============================================================================
# greenland :

g_mesh  = Mesh(g_in_dir + 'submesh.xdmf')
g_Q     = FunctionSpace(g_mesh, 'CG', 1)

g_S      = Function(g_Q)
g_B      = Function(g_Q)
g_adot   = Function(g_Q)
g_qgeo   = Function(g_Q)
g_beta   = Function(g_Q)
g_Mb     = Function(g_Q)
g_Tb     = Function(g_Q)
g_Ts     = Function(g_Q)
g_u      = Function(g_Q)
g_v      = Function(g_Q)
g_w      = Function(g_Q)
g_Ubar5  = Function(g_Q)
g_Ubar10 = Function(g_Q)
g_Ubar20 = Function(g_Q)
g_U_ob   = Function(g_Q)
g_tau_id = Function(g_Q)
g_tau_jd = Function(g_Q)
g_tau_ii = Function(g_Q)
g_tau_ij = Function(g_Q)
g_tau_iz = Function(g_Q)
g_tau_ji = Function(g_Q)
g_tau_jj = Function(g_Q)
g_tau_jz = Function(g_Q)
g_mask   = Function(g_Q)

File(g_in_dir + 'S_s.xml')    >> g_S
File(g_in_dir + 'B_s.xml')    >> g_B
File(g_in_dir + 'adot_s.xml') >> g_adot
File(g_in_dir + 'qgeo_s.xml') >> g_qgeo
File(g_in_dir + 'beta_s.xml') >> g_beta
File(g_in_dir + 'Mb_s.xml')   >> g_Mb
File(g_in_dir + 'Tb_s.xml')   >> g_Tb
File(g_in_dir + 'Ts_s.xml')   >> g_Ts
File(g_in_dir + 'ub_s.xml')   >> g_u
File(g_in_dir + 'vb_s.xml')   >> g_v
File(g_in_dir + 'wb_s.xml')   >> g_w
File(g_in_dir + 'bv/Ubar_5.xml')  >> g_Ubar5
File(g_in_dir + 'bv/Ubar_10.xml') >> g_Ubar10
File(g_in_dir + 'bv/Ubar_20.xml') >> g_Ubar20
File(g_in_dir + 'U_ob_s.xml') >> g_U_ob

File(g_in_dir + 'tau_id_s.xml') >> g_tau_id
File(g_in_dir + 'tau_jd_s.xml') >> g_tau_jd
File(g_in_dir + 'tau_ii_s.xml') >> g_tau_ii
File(g_in_dir + 'tau_ij_s.xml') >> g_tau_ij
File(g_in_dir + 'tau_iz_s.xml') >> g_tau_iz
File(g_in_dir + 'tau_ji_s.xml') >> g_tau_ji
File(g_in_dir + 'tau_jj_s.xml') >> g_tau_jj
File(g_in_dir + 'tau_jz_s.xml') >> g_tau_jz
File(g_in_dir + 'mask_s.xml')   >> g_mask

g_dSdx   = project(g_S.dx(0), g_Q)
g_dSdy   = project(g_S.dx(1), g_Q)

g_dBdx   = project(g_B.dx(0), g_Q)
g_dBdy   = project(g_B.dx(1), g_Q)
                          
g_dHdx   = project((g_S - g_B).dx(0), g_Q)
g_dHdy   = project((g_S - g_B).dx(1), g_Q)

# vectors :
g_beta_v   = g_beta.vector().array()
g_S_v      = g_S.vector().array()
g_B_v      = g_B.vector().array()
g_adot_v   = g_adot.vector().array()
g_qgeo_v   = g_qgeo.vector().array()
g_Mb_v     = g_Mb.vector().array()
g_Tb_v     = g_Tb.vector().array()
g_Ts_v     = g_Ts.vector().array()
g_u_v      = g_u.vector().array()
g_v_v      = g_v.vector().array()
g_w_v      = g_w.vector().array()
g_Ubar5_v  = g_Ubar5.vector().array()
g_Ubar10_v = g_Ubar10.vector().array()
g_Ubar20_v = g_Ubar20.vector().array()
g_U_ob_v   = g_U_ob.vector().array()
g_tau_id_v = g_tau_id.vector().array()
g_tau_jd_v = g_tau_jd.vector().array()
g_tau_ii_v = g_tau_ii.vector().array()
g_tau_ij_v = g_tau_ij.vector().array()
g_tau_iz_v = g_tau_iz.vector().array()
g_tau_ji_v = g_tau_ji.vector().array()
g_tau_jj_v = g_tau_jj.vector().array()
g_tau_jz_v = g_tau_jz.vector().array()
g_mask_v   = g_mask.vector().array()

g_H_v    = g_S_v - g_B_v
g_dHdx_v = g_dHdx.vector().array()
g_dHdy_v = g_dHdy.vector().array()
g_gradH  = sqrt(g_dHdx_v**2 + g_dHdy_v**2 + 1e-16)
g_U_mag  = sqrt(g_u_v**2 + g_v_v**2 + 1e-16)
g_dSdx_v = g_dSdx.vector().array()
g_dSdy_v = g_dSdy.vector().array()
g_gradS  = sqrt(g_dSdx_v**2 + g_dSdy_v**2 + 1e-16)
g_dBdx_v = g_dBdx.vector().array()
g_dBdy_v = g_dBdy.vector().array()
g_gradB  = sqrt(g_dBdx_v**2 + g_dBdy_v**2 + 1e-16)
g_D      = zeros(len(g_B_v))
g_D[g_B_v < 0] = g_B_v[g_B_v < 0]

g_taux = -917.0 * 9.8 * g_H_v * g_dSdx_v
g_tauy = -917.0 * 9.8 * g_H_v * g_dSdy_v
g_tau_mag = sqrt(g_taux**2 + g_tauy**2 + 1e-16)

if sys.argv[1] == 'Ubar':
  g_uhat = g_taux / g_tau_mag
  g_vhat = g_tauy / g_tau_mag

elif sys.argv[1] == 'U' or sys.argv[1] == 'stress':
  g_uhat = g_u_v / g_U_mag
  g_vhat = g_v_v / g_U_mag

g_dBdi = g_dBdx_v * g_uhat + g_dBdy_v * g_vhat
g_dBdj = g_dBdx_v * g_vhat - g_dBdy_v * g_uhat

g_dSdi = g_dSdx_v * g_uhat + g_dSdy_v * g_vhat
g_dSdj = g_dSdx_v * g_vhat - g_dSdy_v * g_uhat

g_dHdi = g_dHdx_v * g_uhat + g_dHdy_v * g_vhat
g_dHdj = g_dHdx_v * g_vhat - g_dHdy_v * g_uhat

g_Ubar_avg = (g_Ubar5_v + g_Ubar10_v + g_Ubar20_v) / 3.0

if sys.argv[1] == 'Ubar':
  g_ini_i    = 917.0 * 9.8 * g_H_v * g_dSdi / (g_Ubar5_v + 0.1)
  g_ini_j    = 917.0 * 9.8 * g_H_v * g_dSdj / (g_Ubar5_v + 0.1)

elif sys.argv[1] == 'U' or sys.argv[1] == 'stress':
  g_ini_i    = 917.0 * 9.8 * g_H_v * g_dSdi / (g_U_mag + 0.1)
  g_ini_j    = 917.0 * 9.8 * g_H_v * g_dSdj / (g_U_mag + 0.1)

# areas of cells for weighting :
g_h_v  = project(CellSize(g_mesh), g_Q).vector().array()

# number of dofs :
g_n = len(g_beta_v)

#===============================================================================
# indicies for each ice sheet :

a_indicies = arange(a_n)
g_indicies = arange(a_n, g_n + a_n)

#===============================================================================
# combined :

beta_v    = hstack((a_beta_v,   g_beta_v  ))
S_v       = hstack((a_S_v,      g_S_v     )) 
B_v       = hstack((a_B_v,      g_B_v     )) 
Ts_v      = hstack((a_Ts_v,     g_Ts_v    ))
gradH     = hstack((a_gradH,    g_gradH   ))
dHdi      = hstack((a_dHdi,     g_dHdi    ))
dHdj      = hstack((a_dHdj,     g_dHdj    ))
gradS     = hstack((a_gradS,    g_gradS   ))
dSdi      = hstack((a_dSdi,     g_dSdi    ))
dSdj      = hstack((a_dSdj,     g_dSdj    ))
D         = hstack((a_D,        g_D       ))
gradB     = hstack((a_gradB,    g_gradB   ))
dBdi      = hstack((a_dBdi,     g_dBdi    ))
dBdj      = hstack((a_dBdj,     g_dBdj    ))
H_v       = hstack((a_H_v,      g_H_v     ))
qgeo_v    = hstack((a_qgeo_v,   g_qgeo_v  ))
adot_v    = hstack((a_adot_v,   g_adot_v  ))
Tb_v      = hstack((a_Tb_v,     g_Tb_v    ))
Mb_v      = hstack((a_Mb_v,     g_Mb_v    ))
u_v       = hstack((a_u_v,      g_u_v     ))
v_v       = hstack((a_v_v,      g_v_v     ))
w_v       = hstack((a_w_v,      g_w_v     ))
Ubar5_v   = hstack((a_Ubar5_v,  g_Ubar5_v ))
Ubar10_v  = hstack((a_Ubar10_v, g_Ubar10_v))
Ubar20_v  = hstack((a_Ubar20_v, g_Ubar20_v))
U_ob_v    = hstack((a_U_ob_v,   g_U_ob_v  ))
U_mag     = hstack((a_U_mag,    g_U_mag   ))
tau_id_v  = hstack((a_tau_id_v, g_tau_id_v))
tau_jd_v  = hstack((a_tau_jd_v, g_tau_jd_v))
tau_ii_v  = hstack((a_tau_ii_v, g_tau_ii_v))
tau_ij_v  = hstack((a_tau_ij_v, g_tau_ij_v))
tau_iz_v  = hstack((a_tau_iz_v, g_tau_iz_v))
tau_ji_v  = hstack((a_tau_ji_v, g_tau_ji_v))
tau_jj_v  = hstack((a_tau_jj_v, g_tau_jj_v))
tau_jz_v  = hstack((a_tau_jz_v, g_tau_jz_v))
ini_i     = hstack((a_ini_i,    g_ini_i   ))
ini_j     = hstack((a_ini_j,    g_ini_j   ))
mask_v    = hstack((a_mask_v,   g_mask_v  ))
h_v       = hstack((a_h_v,      g_h_v     ))


#===============================================================================
# remove areas with garbage data :
valid  = where(mask_v < 1.0)[0]
valid  = intersect1d(valid, where(S_v > 0.0)[0])
valid  = intersect1d(valid, where(beta_v < 1000)[0])
valid  = intersect1d(valid, where(beta_v > 1e-14)[0])
if sys.argv[2] == 'limited':
  valid  = intersect1d(valid, where(U_mag > 20)[0])
else:
  valid  = intersect1d(valid, where(U_mag > 0)[0])
valid  = intersect1d(valid, where(U_ob_v > 1e-9)[0])
valid  = intersect1d(valid, where(Ts_v > 100)[0])
valid  = intersect1d(valid, where(h_v > 0)[0])
valid  = intersect1d(valid, where(S_v - B_v > 60)[0])
valid  = intersect1d(valid, where(adot_v > -100)[0])
#valid  = intersect1d(valid, where(gradS < 0.05)[0])
#valid  = intersect1d(valid, where(gradB < 0.2)[0])
#valid  = intersect1d(valid, where(Mb_v < 0.04)[0])
#valid  = intersect1d(valid, where(Mb_v > 0.0)[0])
#valid  = intersect1d(valid, where(adot_v < 1.2)[0])
#valid  = intersect1d(valid, where(adot_v > -1.0)[0])

#===============================================================================
# individual regions for plotting :

a_valid = intersect1d(a_indicies, valid)
g_valid = intersect1d(g_indicies, valid) - a_n

# to convert to fenics functions for plotting : 
a_conv  = arange(len(a_valid))
g_conv  = arange(len(a_valid), len(g_valid) + len(a_valid))

a_valid_f            = Function(a_Q)
a_valid_f_v          = a_valid_f.vector().array()
a_valid_f_v[a_valid] = 1.0
a_valid_f.vector().set_local(a_valid_f_v)
a_valid_f.vector().apply('insert')

g_valid_f            = Function(g_Q)
g_valid_f_v          = g_valid_f.vector().array()
g_valid_f_v[g_valid] = 1.0
g_valid_f.vector().set_local(g_valid_f_v)
g_valid_f.vector().apply('insert')

measures = DataFactory.get_ant_measures(res=900)
dm       = DataInput(measures, gen_space=False)

rignot   = DataFactory.get_gre_rignot()
drg      = DataInput(rignot, gen_space=False)

betaMax = 200.0

#===============================================================================

#plotIce(dm, a_valid_f, name='valid', direc=a_fn,
#        cmap='gist_yarg', scale='bool', numLvls=12, tp=False,
#        tpAlpha=0.5, show=False)
#
#plotIce(drg, g_valid_f, name='valid', direc=g_fn,
#        cmap='gist_yarg', scale='bool', numLvls=12, tp=False,
#        tpAlpha=0.5, show=False)
#
#a_dBdi_f = Function(a_Q)
#a_dBdj_f = Function(a_Q)
#a_dBdi_f.vector()[a_valid] = dBdi[a_conv]
#a_dBdj_f.vector()[a_valid] = dBdj[a_conv]
#  
#g_dBdi_f = Function(g_Q)
#g_dBdj_f = Function(g_Q)
#g_dBdi_f.vector()[g_valid] = dBdi[g_conv]
#g_dBdj_f.vector()[g_valid] = dBdj[g_conv]
#
#plotIce(dm, a_dBdi_f, name='dBdi', direc=a_fn, 
#        title=r'$\partial_i B$', cmap='RdGy', scale='lin', extend='max',
#        umin=-0.1, umax=0.1, numLvls=12, tp=False, tpAlpha=0.5, show=False)
#
#plotIce(dm, a_dBdj_f, name='dBdj', direc=a_fn, 
#        title=r'$\partial_j B$', cmap='RdGy', scale='lin', extend='max',
#        umin=-0.1, umax=0.1, numLvls=12, tp=False, tpAlpha=0.5, show=False)
#
#plotIce(drg, g_dBdi_f, name='dBdi', direc=g_fn, 
#        title=r'$\partial_i B$', cmap='RdGy', scale='lin', extend='max',
#        umin=-0.1, umax=0.1, numLvls=12, tp=False, tpAlpha=0.5, show=False)
#
#plotIce(drg, g_dBdj_f, name='dBdj', direc=g_fn, 
#        title=r'$\partial_j B$', cmap='RdGy', scale='lin', extend='max',
#        umin=-0.1, umax=0.1, numLvls=12, tp=False, tpAlpha=0.5, show=False)
#
#a_dSdi_f = Function(a_Q)
#a_dSdj_f = Function(a_Q)
#a_dSdi_f.vector()[a_valid] = dSdi[a_conv]
#a_dSdj_f.vector()[a_valid] = dSdj[a_conv]
#  
#g_dSdi_f = Function(g_Q)
#g_dSdj_f = Function(g_Q)
#g_dSdi_f.vector()[g_valid] = dSdi[g_conv]
#g_dSdj_f.vector()[g_valid] = dSdj[g_conv]
#
#plotIce(dm, a_dSdi_f, name='dSdi', direc=a_fn, 
#        title=r'$\partial_i S$', cmap='RdGy', scale='lin', extend='max',
#        umin=-0.1, umax=0.1, numLvls=12, tp=False, tpAlpha=0.5, show=False)
#
#plotIce(dm, a_dSdj_f, name='dSdj', direc=a_fn, 
#        title=r'$\partial_j S$', cmap='RdGy', scale='lin', extend='max',
#        umin=-0.1, umax=0.1, numLvls=12, tp=False, tpAlpha=0.5, show=False)
#
#plotIce(drg, g_dSdi_f, name='dSdi', direc=g_fn, 
#        title=r'$\partial_i S$', cmap='RdGy', scale='lin', extend='max',
#        umin=-0.1, umax=0.1, numLvls=12, tp=False, tpAlpha=0.5, show=False)
#
#plotIce(drg, g_dSdj_f, name='dSdj', direc=g_fn, 
#        title=r'$\partial_j S$', cmap='RdGy', scale='lin', extend='max',
#        umin=-0.1, umax=0.1, numLvls=12, tp=False, tpAlpha=0.5, show=False)

#===============================================================================
# cell declustering :
a_n_v    = len(a_valid)
g_n_v    = len(g_valid)
n        = len(valid)

a_h_v    = a_h_v[a_valid]
g_h_v    = g_h_v[g_valid]
a_A      = sum(a_h_v)
g_A      = sum(g_h_v)

gam      = float(g_n_v) / float(a_n_v)

a_wt     = n**2 / float(a_n_v) * a_h_v / a_A
g_wt     = n**2 / float(g_n_v) * g_h_v / g_A

wt       = hstack((a_wt, g_wt))

#h_v      = h_v[valid]
#A        = sum(h_v)
#wt       = n * h_v / A
#beta_bar = 1.0/n * sum(beta_v[valid] * wt)

#===============================================================================
#data = [beta_v,  S_v,     B_v,    gradS,   gradB, 
#        H_v,     adot_v,  Ts_v,   Tb_v,    Mb_v,
#        Ubar5_v, u_v,     v_v,    w_v,     U_mag]
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
#         r'$\Vert \bar{\mathbf{u}}_{bv} \Vert$',
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
#fn = 'images/data.png'
##savefig(fn, dpi=100)
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
v25  = ini_i
v26  = ini_j
v27  = dBdi
v28  = dBdj
v29  = dSdi
v30  = dSdj
v31  = gradH
v32  = dHdi
v33  = dHdj

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
x26  = v26[valid]
x27  = v27[valid]
x28  = v28[valid]
x29  = v29[valid]
x30  = v30[valid]
x31  = v31[valid]
x32  = v32[valid]
x33  = v33[valid]

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
         r'ini$_i$',
         r'ini$_j$',
         r'$\partial_i B$',
         r'$\partial_j B$',
         r'$\partial_i S$',
         r'$\partial_j S$',
         r'$\Vert \nabla H \Vert$',
         r'$\partial_i H$',
         r'$\partial_j H$']

X      = [x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16,x17,x18,
          x19,x20,x21,x22,x23,x24,x25,x26,x27,x28,x29,x30,x31,x32,x33]
V      = [v0,v1,v2,v3,v4,v5,v6,v7,v8,v9,v10,v11,v12,v13,v14,v15,v16,v17,v18,
          v19,v20,v21,v22,v23,v24,v25,v26,v27,v28,v29,v30,v31,v32,v33]

# with stress terms :
if sys.argv[1] == 'stress':
  index  = [0,1,5,7,16,27,29,17,19,20,22,23]

# U instead of Ubar :
elif sys.argv[1] == 'U':
  index  = [0,1,5,7,16,27,29]

# independent only :
elif sys.argv[1] == 'Ubar':
  index  = [0,1,5,7,13,27,29]

ii     = index
ii_int = []
ii_int.extend(ii)

for i,m in enumerate(ii):
  if sys.argv[1] == 'U' or sys.argv[1] == 'Ubar':
    k = i
  else:
    k = i+1
  for j,n in enumerate(ii[k:]):
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
ex_n = array(ex_n)
 
#show()

#==============================================================================
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
#ax.set_ylim([0,0.020])
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

a_yhat_f = Function(a_Q)
a_resi_f = Function(a_Q)
a_yhat_f.vector()[a_valid] = yhat[a_conv]
a_resi_f.vector()[a_valid] = resid[a_conv]
  
g_yhat_f = Function(g_Q)
g_resi_f = Function(g_Q)
g_yhat_f.vector()[g_valid] = yhat[g_conv]
g_resi_f.vector()[g_valid] = resid[g_conv]

plotIce(dm, a_yhat_f, name='GLM_beta', direc=a_fn, 
        title=r'$\hat{\beta}$', cmap='gist_yarg', scale='log', extend='max',
        umin=1.0, umax=betaMax, numLvls=12, tp=False, tpAlpha=0.5, show=False)

plotIce(dm, a_resi_f, name='GLM_resid', direc=a_fn, 
        title=r'$d$', cmap='RdGy', scale='lin', extend='both', 
        umin=-50, umax=50, numLvls=13, tp=False, tpAlpha=0.5, show=False)

plotIce(drg, g_yhat_f, name='GLM_beta', direc=g_fn, 
        title=r'$\hat{\beta}$', cmap='gist_yarg', scale='log', extend='max',
        umin=1.0, umax=betaMax, numLvls=12, tp=False, tpAlpha=0.5, show=False)

plotIce(drg, g_resi_f, name='GLM_resid', direc=g_fn, 
        title=r'$d$', cmap='RdGy', scale='lin', extend='both',
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

#rtol  = 30
#Xr    = array(X)
#Xr    = Xr[:, resid < rtol]
#yhat  = yhat[resid < rtol]
#resid = resid[resid < rtol]

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
fn = 'images/stats/' + file_n + 'GLM_newton_resid.png'
tight_layout()
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

srt = argsort(abs(ahat))[::-1]
f   = open('dat/' + file_n + 'alpha.dat', 'w')
for n, a, c in zip(ex_n[srt], ahat[srt], ci[srt]):
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

  a_yhat_f = Function(a_Q)
  a_resi_f = Function(a_Q)
  a_yhat_f.vector()[a_valid] = yhat_n[a_conv]
  a_resi_f.vector()[a_valid] = resid_n[a_conv]
    
  g_yhat_f = Function(g_Q)
  g_resi_f = Function(g_Q)
  g_yhat_f.vector()[g_valid] = yhat_n[g_conv]
  g_resi_f.vector()[g_valid] = resid_n[g_conv]
  
  plotIce(dm, a_yhat_f, name='GLM_beta_reduced', direc=a_fn, 
          title=r'$\hat{\beta}$', cmap='gist_yarg', scale='log', 
          umin=1.0, umax=betaMax, numLvls=12, tp=False, tpAlpha=0.5, show=False)
  
  plotIce(dm, a_resi_f, name='GLM_resid_reduced', direc=a_fn, 
          title=r'$d$', cmap='RdGy', scale='lin', 
          umin=-50, umax=50, numLvls=13, tp=False, tpAlpha=0.5, show=False)
  
  plotIce(drg, g_yhat_f, name='GLM_beta_reduced', direc=g_fn, 
          title=r'$\hat{\beta}$', cmap='gist_yarg', scale='log', 
          umin=1.0, umax=betaMax, numLvls=12, tp=False, tpAlpha=0.5, show=False)
  
  plotIce(drg, g_resi_f, name='GLM_resid_reduced', direc=g_fn, 
          title=r'$d$', cmap='RdGy', scale='lin', 
          umin=-50, umax=50, numLvls=13, tp=False, tpAlpha=0.5, show=False)
    
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
  srt = argsort(abs(ahat_n))[::-1]
  fn  = open('dat/'+file_n+'alpha_reduced.dat', 'w')
  for n, a, c in zip(ex_a[srt], ahat_n[srt], ci_n[srt]):
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



