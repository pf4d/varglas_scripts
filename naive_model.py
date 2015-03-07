from pylab          import *
from scipy.stats    import t, distributions
from scipy.special  import fdtrc

chi2cdf = distributions.chi2.cdf

import sys
src_directory = '../statistical_modeling'
sys.path.append(src_directory)

from scipy.stats               import probplot 
from src.regstats              import prbplotObj
from fenics                    import *
from pylab                     import *

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
qbar_v = hstack((a_qbar_v, g_qbar_v))

valid  = where(beta_v > 1e-10)[0]
valid  = intersect1d(valid, where(Ubar_v > 0)[0])
valid  = intersect1d(valid, where(abs(Mb_v) < 1500)[0])
valid  = intersect1d(valid, where(S_v > -100)[0])
valid  = intersect1d(valid, where(Ts_v > 100)[0])

#g_data = [log(g_beta_v + 1), g_S_v, g_B_v, g_H_v, g_adot_v, g_Mb_v, 
#          g_Tb_v, g_Ts_v, log(g_Ubar_v + 1), g_qbar_v, g_u_v, 
#          g_v_v, g_w_v, log(g_U_mag + 1), sqrt(g_gradS + 1), sqrt(g_gradB + 1)]
g_data = [g_beta_v, g_S_v, g_B_v, g_H_v, g_adot_v, g_Mb_v, 
          g_Tb_v, g_Ts_v, g_Ubar_v, g_qbar_v, g_u_v, 
          g_v_v, g_w_v, g_U_mag, g_gradS, g_gradB]
a_data = [a_beta_v, a_S_v, a_B_v, a_H_v, a_adot_v, a_Mb_v,
          a_Tb_v, a_Ts_v, a_Ubar_v, a_qbar_v, a_u_v,
          a_v_v, a_w_v, a_U_mag, a_gradS, a_gradB]
data   = array([beta_v, S_v, B_v, H_v, adot_v, Mb_v, 
                Tb_v, Ts_v, Ubar_v, qbar_v, u_v, 
                v_v, w_v, U_mag, gradS, gradB])
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
##savefig(fn, dpi=100)
#show()
#
#print "sample size:", len(g_valid)

n       = len(beta_v)
idx     = arange(n)
test_i  = randint(n, size = 1000)
train_i = setdiff1d(idx, test_i)

test    = data[1:,test_i].T
train   = data[1:,train_i].T

test_class  = data[0,test_i]
train_class = data[0, train_i]

def plot_dist(train, train_class, nbins, name, norm=True):
  """
  Function which plots a histogram of data <train> with associated training 
  class <train_class> with number of bins <nbins>, title <name>.  If <norm> 
  is True, the bins are normalized.

  Returns the means, standard deviations, bins, and bin counts as a tuple.
  
  Image is saved in directory ../doc/images/.
  """
  vs    = unique(train_class)              # different values
  m     = len(vs)                          # number of distinct values
  n     = shape(train)[1]                  # number of different attributes
  fig   = figure(figsize=(10,8))           # figure instance
  bi_a  = []                               # array of bins
  ct_a  = []                               # array of bin counts
  # iterate over each attribute i :
  for i in range(n):
    ax    = fig.add_subplot(220 + i+1)     # create a subplot
    mini  = train[:,i].min()               # min of attribute i
    maxi  = train[:,i].max()               # max of attribute i
    bi_i  = []                             # bin array for attribute i
    ct_i  = []                             # bin count array for attribute i
    # iterate over each class j :
    for j in range(m):
      wj   = where(train_class == j+1)[0]  # indicies for class j
      xj   = train[wj,i]                   # array of training values j
      rngi = linspace(mini, maxi, 1000)    # range of plotting normal curve
      lblh = classes[j+1] + ' (%i)' % (j+1)
      lbln = r'$\mathcal{N}(\mu_%i, \sigma_%i^2)$' % (j+1, j+1)
      # function returns the bin counts and bins
      ct, bins, ign = ax.hist(train[wj, i], nbins, label=lblh, alpha=0.7, 
                              normed=norm)
      bi_i.append(bins)                    # add the bin to the list
      ct_i.append(ct)                      # add the bin counts to the list
    bi_a.append(array(bi_i))               # convert to numpy array
    ct_a.append(array(ct_i))               # convert to numpy array
  
    ax.set_title(attrib[i])                # set the title
    if i == 2:
      leg = ax.legend()                    # add the legend
      leg.get_frame().set_alpha(0.5)       # transparent legend
    ax.grid()                              # gridlines
  tight_layout() 
  savefig('../doc/images/' + name, dpi=300)
  plt.close(fig)
  #show()
  return array(bi_a), array(ct_a)

out = plot_dist(train[1:], train[0], nbins=10, 'training', norm=False)
