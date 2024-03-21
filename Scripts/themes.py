import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import numpy as np

colors = {'ISO' :'red',
          'CUBE':'orange',
          'XISO':'teal',
          'TRIG':'cyan',
          'TET' :'blue',
          'ORTH':'purple',
          'MONO':'black'}

symmetry_classes = ['ISO',
                    'CUBE',
                    'XISO',
                    'TRIG',
                    'TET',
                    'ORTH',
                    'MONO']

def global_plot(lons, lats, values, title, write=False, figname=None,
                central_longitude=0, cmap='jet', reverse_cmap=False, vmin=None,
                vmax=None, extend='neither', ticks=None):

    fig = plt.figure(figsize=(10, 5))
    ax = fig.add_subplot(1, 1, 1, projection=ccrs.Robinson(
        central_longitude=central_longitude))
    ax.set_global()
    ax.coastlines()

    cmap = plt.get_cmap(cmap)

    if reverse_cmap:
        cmap = cmap.reversed()

    if vmin == None and vmax == None:
        sc = ax.scatter(lons, lats, c=values,
                        cmap=cmap, alpha=0.6, edgecolors='w', linewidth=0.5,
                        transform=ccrs.PlateCarree())
    else:
        sc = ax.scatter(lons, lats, c=values, vmin=vmin, vmax=vmax,
                        cmap=cmap, alpha=0.6, edgecolors='w', linewidth=0.5,
                        transform=ccrs.PlateCarree())

    clb = plt.colorbar(sc, orientation='vertical', extend=extend,
                       fraction=0.0235, pad=0.05)
    if ticks is not None:
        clb.set_ticks(ticks)
    clb.ax.tick_params(labelsize=10)
    ax.set_title(title, fontsize=20)

    if write:
        if figname:
            plt.savefig(f'{figname}.png', bbox_inches='tight')
    plt.show()

def correlation_plot(x, y, xlabel, ylabel, write=False, figname=None):
    number_of_bins = 20
    bins = np.linspace(0, max(x), number_of_bins + 1)
    bin_width = bins[1] - bins[0]
    index = np.digitize(x, bins[:-1]) - 1

    bin_mids = bins[:-1] + bin_width / 2

    running_median = np.array(
        [np.median(y[index == i]) for i in range(number_of_bins)])
    running_mad = np.array(
        [np.median(np.abs(y[index == i] - running_median[i])) for i in
         range(number_of_bins)])

    plt.figure(figsize=(6, 6))
    plt.scatter(x, y)
    plt.plot(bin_mids, running_median, c='r')
    plt.plot(bin_mids, running_median + running_mad, 'r--')
    plt.plot(bin_mids, running_median - running_mad, 'r--')
    plt.xticks(np.round(np.linspace(0, max(x), 6), 3))

    plt.xlabel(xlabel, fontsize=15)
    plt.ylabel(ylabel, fontsize=15)

    if write: plt.savefig(figname, bbox_inches='tight')
    plt.show()