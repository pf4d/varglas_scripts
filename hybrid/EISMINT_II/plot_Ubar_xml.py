from varglas.helper  import plot_variable
from fenics          import *
from pylab           import *
from matplotlib.ticker import ScalarFormatter, LogFormatter

in_dir  = 'dump2/Ubar_xml/'
img_dir = 'images_stats_Ubar/'

set_log_active(False)

alpha = 1.0 * pi / 180
L     = 2

mesh = Mesh('meshes/circle_mesh.xml')

L = 800.0

for x in mesh.coordinates():
  # transform x :
  x[0]  = x[0]  * L
  # transform y :
  x[1]  = x[1]  * L

Q    = FunctionSpace(mesh, 'CG', 1)

beta  = Function(Q)
U     = Function(Q)
S     = Function(Q)
H     = Function(Q)
u     = Function(Q)
v     = Function(Q)
w     = Function(Q)

Bmin = 30.0
Bmax = 400.0

Smin = -2233
Smax = 3500

Umin = 1.0
Umax = 200.0

Hmin = 1.0
Hmax = 3500

#===============================================================================
# save S colorbar :
fig = figure(figsize=(8,1))
ax  = fig.add_subplot(111)

# Set the colormap and norm to correspond to the data for which
# the colorbar will be used.
cmap = mpl.cm.gist_yarg

levels    = linspace(Smin, Smax, 12)
formatter = ScalarFormatter()
norm      = mpl.colors.BoundaryNorm(levels, cmap.N)

# ColorbarBase derives from ScalarMappable and puts a colorbar
# in a specified axes, so it has everything needed for a
# standalone colorbar.  There are many more kwargs, but the
# following gives a basic continuous colorbar with ticks
# and labels.
cb1 = mpl.colorbar.ColorbarBase(ax, cmap=cmap,
                                norm=norm,
                                orientation='horizontal',
                                format=formatter,
                                spacing='uniform',
                                ticks=levels)
tight_layout()
savefig(img_dir + 'S_cb.png', dpi=300)

#===============================================================================
# save H colorbar :
fig = figure(figsize=(8,1))
ax  = fig.add_subplot(111)

# Set the colormap and norm to correspond to the data for which
# the colorbar will be used.
cmap = mpl.cm.gist_yarg

levels    = linspace(Hmin, Hmax, 12)
formatter = ScalarFormatter()
norm      = mpl.colors.BoundaryNorm(levels, cmap.N)

# ColorbarBase derives from ScalarMappable and puts a colorbar
# in a specified axes, so it has everything needed for a
# standalone colorbar.  There are many more kwargs, but the
# following gives a basic continuous colorbar with ticks
# and labels.
cb1 = mpl.colorbar.ColorbarBase(ax, cmap=cmap,
                                norm=norm,
                                orientation='horizontal',
                                format=formatter,
                                spacing='uniform',
                                ticks=levels)
tight_layout()
savefig(img_dir + 'H_cb.png', dpi=300)

#===============================================================================
# save U colorbar :
fig = figure(figsize=(8,1))
ax  = fig.add_subplot(111)

# Set the colormap and norm to correspond to the data for which
# the colorbar will be used.
cmap = mpl.cm.gist_yarg

levels    = np.logspace(np.log10(Umin), np.log10(Umax), 12)
formatter = ScalarFormatter()
norm      = mpl.colors.BoundaryNorm(levels, cmap.N)

# ColorbarBase derives from ScalarMappable and puts a colorbar
# in a specified axes, so it has everything needed for a
# standalone colorbar.  There are many more kwargs, but the
# following gives a basic continuous colorbar with ticks
# and labels.
cb1 = mpl.colorbar.ColorbarBase(ax, cmap=cmap,
                                norm=norm,
                                orientation='horizontal',
                                format=formatter,
                                spacing='uniform',
                                ticks=levels)
tight_layout()
savefig(img_dir + 'U_mag_cb.png', dpi=300)

#===============================================================================
# save beta colorbar :
fig = figure(figsize=(8,1))
ax  = fig.add_subplot(111)

# Set the colormap and norm to correspond to the data for which
# the colorbar will be used.
cmap = mpl.cm.gist_yarg

levels    = np.logspace(np.log10(Bmin), np.log10(Bmax), 12)
formatter = ScalarFormatter()
norm      = mpl.colors.BoundaryNorm(levels, cmap.N)

# ColorbarBase derives from ScalarMappable and puts a colorbar
# in a specified axes, so it has everything needed for a
# standalone colorbar.  There are many more kwargs, but the
# following gives a basic continuous colorbar with ticks
# and labels.
cb1 = mpl.colorbar.ColorbarBase(ax, cmap=cmap,
                                norm=norm,
                                orientation='horizontal',
                                format=formatter,
                                spacing='uniform',
                                ticks=levels)
tight_layout()
savefig(img_dir + 'beta_cb.png', dpi=300)

#===============================================================================
# plot the observed values :
for i in range(0,35000,100):
  File(in_dir + 'beta_' + str(i) + '.xml') >> beta
  File(in_dir + 'S_'    + str(i) + '.xml') >> S
  File(in_dir + 'H_'    + str(i) + '.xml') >> H
  File(in_dir + 'u_s_'  + str(i) + '.xml') >> u
  File(in_dir + 'v_s_'  + str(i) + '.xml') >> v
  File(in_dir + 'w_s_'  + str(i) + '.xml') >> w

  u_v      = u.vector().array()
  v_v      = v.vector().array()
  w_v      = w.vector().array()
  U_mag_v  = sqrt(u_v**2 + v_v**2 + w_v**2 + 1e-16)
  U.vector()[:] = U_mag_v

  cb = False
  
  plot_variable(beta, 'beta_'+str(i),  img_dir, tp=False,
                umin=Bmin, umax=Bmax, scale='log',
                title='$\\beta$', use_colorbar=cb, hide_axis=True,
                show=False, hide_ax_tick_labels=True, label_axes=False,
                colorbar_loc='bottom')
  plot_variable(U,    'U_mag_'+str(i), img_dir, tp=False,
                umin=Umin, umax=Umax, scale='log', hide_axis=True,
                title='$\\Vert \\mathbf{u}_S \\Vert$',
                show=False, hide_ax_tick_labels=cb, label_axes=False,
                use_colorbar=cb, colorbar_loc='bottom')
  plot_variable(S,    'S_'+str(i),     img_dir, tp=False,
                umin=Smin, umax=Smax, scale='lin',
                title='$S$', use_colorbar=cb, hide_axis=True,
                show=False, hide_ax_tick_labels=True, label_axes=False,
                colorbar_loc='bottom')
  plot_variable(H,    'H_'+str(i),     img_dir, tp=False,
                umin=Hmin, umax=Hmax, scale='lin',
                title='$H$', use_colorbar=cb, hide_axis=True,
                show=False, hide_ax_tick_labels=True, label_axes=False,
                colorbar_loc='bottom')



