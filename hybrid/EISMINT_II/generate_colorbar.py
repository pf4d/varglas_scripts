'''
Make a colorbar as a separate figure.
'''
  
#  # countour levels :
#  if scale == 'log':
#    v[v < vmin] = vmin + 1e-12
#    v[v > vmax] = vmax - 1e-12
#    from matplotlib.ticker import LogFormatter
#    levels      = np.logspace(np.log10(vmin), np.log10(vmax), numLvls)
#    formatter   = LogFormatter(10, labelOnlyBase=False)
#    norm        = colors.LogNorm()
#  
#  elif scale == 'lin':
#    v[v < vmin] = vmin + 1e-12
#    v[v > vmax] = vmax - 1e-12
#    from matplotlib.ticker import ScalarFormatter
#    levels    = np.linspace(vmin, vmax, numLvls)
#    formatter = ScalarFormatter()
#    norm      = None
#  
#  elif scale == 'bool':
#    from matplotlib.ticker import ScalarFormatter
#    levels    = [0, 1, 2]
#    formatter = ScalarFormatter()
#    norm      = None
#
#  fig = plt.figure(figsize=(11,10))
#  ax  = fig.add_subplot(111)
#  
#  c = ax.tricontourf(x, y, t, v, levels=levels, norm=norm, 
#                     cmap=pl.get_cmap(cmap))
#  plt.axis('equal')
#  
#  if tp == True:
#    p = ax.triplot(x, y, t, '-', lw=0.2, alpha=tpAlpha)
#  ax.set_xlim([x.min(), x.max()])
#  ax.set_ylim([y.min(), y.max()])
#  if label_axes:
#    ax.set_xlabel(r'$x$')
#    ax.set_ylabel(r'$y$')
#  if hide_ax_tick_labels:
#    ax.set_xticklabels([])
#    ax.set_yticklabels([])
#  if hide_axis:
#    plt.axis('off')
#  
#  # include colorbar :
#  if scale != 'bool' and use_colorbar:
#    divider = make_axes_locatable(plt.gca())
#    cax  = divider.append_axes(colorbar_loc, "5%", pad="3%")
#    cbar = plt.colorbar(c, cax=cax, format=formatter, 
#                        ticks=levels) 
#    tit = plt.title(title)

from matplotlib import pyplot
from matplotlib.ticker import ScalarFormatter
import matplotlib as mpl
import numpy      as np

mpl.rcParams['font.family']     = 'serif'
mpl.rcParams['legend.fontsize'] = 'medium'

Bmin = 30.0
Bmax = 400.0

Smin = -2233
Smax = 3500

Umin = 1.0
Umax = 200.0

Hmin = 1.0
Hmax = 3500

# Make a figure and axes with dimensions as desired.
fig = pyplot.figure(figsize=(8,1))
ax1 = fig.add_subplot(111)

# Set the colormap and norm to correspond to the data for which
# the colorbar will be used.
cmap = mpl.cm.gist_yarg

levels    = np.linspace(Smin, Smax, 12)
formatter = ScalarFormatter()
norm      = mpl.colors.BoundaryNorm(levels, cmap.N)

# ColorbarBase derives from ScalarMappable and puts a colorbar
# in a specified axes, so it has everything needed for a
# standalone colorbar.  There are many more kwargs, but the
# following gives a basic continuous colorbar with ticks
# and labels.
cb1 = mpl.colorbar.ColorbarBase(ax1, cmap=cmap,
                                norm=norm,
                                orientation='horizontal',
                                format=formatter,
                                spacing='uniform',
                                ticks=levels)
pyplot.tight_layout()
savefig(img_dir + 'beta_cb.png', dpi=300)
pyplot.show()
