from varglas.helper  import plot_variable
from fenics          import *
from pylab           import *
from matplotlib.ticker import ScalarFormatter, LogFormatter

in_dir  = 'dump/EISMINT_vars/'
img_dir = 'images_EISMINT_vars/'

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

adot  = Function(Q)
T     = Function(Q)
B     = Function(Q)

File(in_dir + 'adot.xml') >> adot
File(in_dir + 'T_s.xml')  >> T
File(in_dir + 'B.xml')    >> B

adotMin = adot.vector().min()
adotMax = adot.vector().max() + 1e-10

Tmin = T.vector().min()
Tmax = T.vector().max()

Bmin = B.vector().min()
Bmax = B.vector().max()

#===============================================================================
# save adot colorbar :
fig = figure(figsize=(8,1))
ax  = fig.add_subplot(111)

# Set the colormap and norm to correspond to the data for which
# the colorbar will be used.
cmap = mpl.cm.gist_yarg

levels    = linspace(adotMin, adotMax, 12)
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
savefig(img_dir + 'adot_cb.png', dpi=300)

#===============================================================================
# save T colorbar :
fig = figure(figsize=(8,1))
ax  = fig.add_subplot(111)

# Set the colormap and norm to correspond to the data for which
# the colorbar will be used.
cmap = mpl.cm.gist_yarg

levels    = linspace(Tmin, Tmax, 12)
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
savefig(img_dir + 'T_cb.png', dpi=300)

#===============================================================================
# save B colorbar :
fig = figure(figsize=(8,1))
ax  = fig.add_subplot(111)

# Set the colormap and norm to correspond to the data for which
# the colorbar will be used.
cmap = mpl.cm.gist_yarg

levels    = linspace(Bmin, Bmax, 12)
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
savefig(img_dir + 'B_cb.png', dpi=300)


#===============================================================================
# plot the observed values :

cb = False

plot_variable(adot,    'adot',     img_dir, tp=False,
              umin=adotMin, umax=adotMax, scale='lin',
              title='$\dot{a}$', use_colorbar=cb, hide_axis=True,
              show=False, hide_ax_tick_labels=True, label_axes=False,
              colorbar_loc='bottom')

plot_variable(B,       'B',        img_dir, tp=False,
              umin=Bmin,    umax=Bmax,    scale='lin',
              title='$B$',       use_colorbar=cb, hide_axis=True,
              show=False, hide_ax_tick_labels=True, label_axes=False,
              colorbar_loc='bottom')

plot_variable(T,       'T_s',      img_dir, tp=False,
              umin=Tmin,    umax=Tmax,    scale='lin',
              title='$T_S$',     use_colorbar=cb, hide_axis=True,
              show=False, hide_ax_tick_labels=True, label_axes=False,
              colorbar_loc='bottom')



