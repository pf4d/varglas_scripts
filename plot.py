from mpl_toolkits.basemap import Basemap
from pylab                import *
from matplotlib           import colors
from pyproj               import *
from varglas.io           import DataInput
from fenics               import Function

mpl.rcParams['font.family']     = 'serif'
mpl.rcParams['legend.fontsize'] = 'medium'

def plotIce(di, u, cmap='gist_gray',  scale='lin', name='', 
            numLvls=12, tp=False, tpAlpha=0.5):
  """
  INPUTS :

    di :
      DataInput object with desired projection
    u :
      solution to plot; can be either a function on a 2D mesh, or a string 
      key to matrix variable in <di>.data.
    cmap :
      colormap to use - see images directory for sample and name
    scale :
      scale to plot, either 'log' or 'lin'
    name :
      title of the plot, latex accepted
    numLvls :
      number of levels for field values
    tp :
      boolean determins plotting of triangle overlay
    tpAlpha :
      alpha level of triangles 0.0 (transparent) - 1.0 (opaque)
  
  OUTPUT :
 
    A sigle 250 dpi .png in the source directory.
  
  """
  #=============================================================================
  # data gathering :
  cmap = 'gist_earth'
  cmap = 'RdYlGn'
  cmap = 'RdGy'
  cmap = 'gist_gray'
  cmap = 'Purples'
  cmap = 'Reds'
  cmap = 'Oranges'
  
  if isinstance(u, str):
    vx,vy  = meshgrid(di.x, di.y)
    v      = di.data[u]

  elif isinstance(u, Function):
    mesh  = u.function_space().mesh()
    coord = mesh.coordinates()
    fi    = mesh.cells()
    v     = u.compute_vertex_values(mesh)
    vx    = coord[:,0]
    vy    = coord[:,1]


  #=============================================================================
  # Get lon,lat from mesh coordinates
  lon,lat = di.proj(vx, vy, inverse=True)
  
  # the width/height numbers were calculated from vertices from a mesh :
  #w = 1.05 * (vx.max() - vx.min())
  #h = 1.05 * (vy.max() - vy.min())

  # Antarctica :
  if di.cont == 'antarctica':
    w   = 5513335.22665
    h   = 4602848.6605
    fig = plt.figure(figsize=(14,10))
    ax  = fig.add_axes()
    
    # new projection :
    m = Basemap(ax=ax, width=w, height=h, resolution='h', 
                projection='stere', lat_ts=-71, 
                lon_0=0, lat_0=-90)
   
    # draw lat/lon grid lines every 5 degrees.
    # labels = [left,right,top,bottom]
    m.drawmeridians(np.arange(0, 360, 20.0),
                    color = 'black',
                    labels = [True, False, True, True])
    m.drawparallels(np.arange(-90, 90, 5.0), 
                    color = 'black', 
                    labels = [True, False, True, True])
    m.drawmapscale(-130, -68, 0, -90, 400, 
                   yoffset  = 0.01 * (m.ymax - m.ymin), 
                   barstyle = 'fancy')
 
  # Greenland : 
  elif di.cont == 'greenland':
    w   = 1532453.49654
    h   = 2644074.78236
    fig = plt.figure(figsize=(8,10))
    ax  = fig.add_axes()
    
    # new projection :
    m = Basemap(ax=ax, width=w, height=h, resolution='h', 
                projection='stere', lat_ts=71, 
                lon_0=-41.5, lat_0=71)
    
    # draw lat/lon grid lines every 5 degrees.
    # labels = [left,right,top,bottom]
    m.drawmeridians(np.arange(0, 360, 5.0),
                    color = 'black',
                    labels = [False, False, False, True])
    m.drawparallels(np.arange(-90, 90, 5.0), 
                    color = 'black', 
                    labels = [True, False, True, False])
    m.drawmapscale(-34, 60.5, -41.5, 71, 400, 
                   yoffset  = 0.01 * (m.ymax - m.ymin), 
                   barstyle = 'fancy')

  # conversion to projection coordinates from lon, lat :
  x, y  = m(lon, lat)
 
  m.drawcoastlines(linewidth=0.25, color = 'black')
  #m.shadedrelief()
  #m.bluemarble()
  #m.etopo()
  

  #=============================================================================
  # plotting :
  
  # countour levels :
  if scale == 'log':
    from matplotlib.ticker import LogFormatter
    vmax             = ceil(v.max())
    vmin             = floor(v.min())
    v[where(v<=1.0)] = 1.0
    levels           = logspace(0.0, log10(vmax+1), numLvls)
    formatter        = LogFormatter(10, labelOnlyBase=False)
    norm             = colors.LogNorm()
  
  elif scale == 'lin':
    from matplotlib.ticker import ScalarFormatter
    vmin      = floor(v.min())
    vmax      = ceil(v.max())
    levels    = linspace(vmin, vmax+1, numLvls)
    formatter = ScalarFormatter()
    norm      = None
  
  elif scale == 'bool':
    from matplotlib.ticker import ScalarFormatter
    levels    = [0, 1, 2]
    formatter = ScalarFormatter()
    norm      = None
  

  if isinstance(u, str):
    #cs = pcolor(x, y, v, cmap=get_cmap(cmap), norm=norm)
    cs = contourf(x, y, v, levels=levels, 
                  cmap=get_cmap(cmap), norm=norm)
  
  elif isinstance(u, Function):
    #cs = tripcolor(x, y, fi, v, shading='gouraud', 
    #               cmap=get_cmap(cmap), norm=norm)
    cs = tricontourf(x, y, fi, v, levels=levels, 
                     cmap=get_cmap(cmap), norm=norm)
  
  # plot triangles :
  if tp == True:
    tp = triplot(x, y, fi, '-', lw=0.2, alpha=tpAlpha)

  # include colorbar :
  cbar = m.colorbar(cs, format=formatter, 
                    ticks=around(levels,decimals=1), 
                    location='right', pad="5%")
  
  # title :
  #tit = title(name)
  #tit.set_fontsize(40)
  
  #savefig(name + '.png', dpi=250)
  show()



