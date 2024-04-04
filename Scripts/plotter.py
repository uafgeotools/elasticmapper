import cartopy.crs as ccrs
from cartopy.feature import ShapelyFeature
import geopandas as gpd
import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt

class Plotter:

    def __init__(self, lons=None, lats=None, central_longitude=0, write=False,
                 highlight=None):
        self.lons = lons
        self.lats = lats
        self.central_longitude = central_longitude
        self.write = write
        self.highlight = highlight

    def correlation_plots(self, x, y, xlabel, ylabel,
                      cmaps1=('jet', False, None, None), cmaps2=None,
                      residual=False, rcmaps=('jet', False, None, None),
                      trend=False, aspect=None, filename=None):

        if cmaps2 is None: cmaps2 = cmaps1

        if residual:
            fig = plt.figure(figsize=(18, 10))
            gs = gridspec.GridSpec(2, 2, width_ratios=[1, 1])
            z = y - x
            zlabel = f"{ylabel} - {xlabel}"
            self.global_plot(z, zlabel, cmaps=rcmaps,
                             fig=fig, figloc=gs[0,1], fig_label='c)')
            self.correlation_plot(x, y, xlabel, ylabel,
                                  fig=fig, figloc=gs[1,1], fig_label='d)',
                                  residual=residual, trend=trend, aspect=aspect)
        else:
            fig = plt.figure(figsize=(14, 10))
            gs = gridspec.GridSpec(2, 3, width_ratios=[1, 0.05, 0.5])
            self.correlation_plot(x, y, xlabel, ylabel,
                                  fig=fig, figloc=gs[:,2], fig_label='c)',
                                  trend=trend, aspect=aspect)

        self.global_plot(y, ylabel, cmaps=cmaps2,
                         fig=fig, figloc=gs[0,0], fig_label='a)')
        self.global_plot(x, xlabel, cmaps=cmaps1,
                         fig=fig, figloc=gs[1,0], fig_label='b)')

        plt.tight_layout()
        if self.write: plt.savefig(f"{filename}.png", bbox_inches='tight')
        plt.show()

    def multi_correlation_plots(self, x, y, z, xlabel, ylabel, zlabel,
                                cmaps1=('jet', False, None, None),
                                cmaps2=None, cmaps3=None,
                                trend=False, aspect=None, filename=None):

        if cmaps2 is None: cmaps2 = cmaps1
        if cmaps3 is None: cmaps3 = cmaps1

        fig = plt.figure(figsize=(14, 15))
        gs = gridspec.GridSpec(3, 3, width_ratios=[1, 0.05, 0.5])

        self.global_plot(z, zlabel, cmaps=cmaps3,
                         fig=fig, figloc=gs[0, 0], fig_label='a)')
        self.global_plot(y, ylabel, cmaps=cmaps2,
                         fig=fig, figloc=gs[1, 0], fig_label='b)')
        self.global_plot(x, xlabel, cmaps=cmaps1,
                         fig=fig, figloc=gs[2, 0], fig_label='c)')
        self.correlation_plot(z, x, zlabel, xlabel,
                              fig=fig, figloc=gs[0,2], fig_label='d)',
                              trend=trend, aspect=aspect)
        self.correlation_plot(y, z, ylabel, zlabel,
                              fig=fig, figloc=gs[1,2], fig_label='e)',
                              trend=trend, aspect=aspect)
        self.correlation_plot(x, y, xlabel, ylabel,
                              fig=fig, figloc=gs[2,2], fig_label='f)',
                              trend=trend, aspect=aspect)

        plt.tight_layout()
        if self.write: plt.savefig(f"{filename}.png", bbox_inches='tight')
        plt.show()

    def global_plot(self, values, title, cmaps=('jet', False, None, None),
                    fig=None, figloc=111, fig_label=None, filename=None):

        show = False
        if fig is None:
            show = True
            fig = plt.figure(figsize=(9, 5))

        cmap, reverse_cmap, vmin, vmax = cmaps

        cmap = plt.get_cmap(cmap)
        if reverse_cmap: cmap = cmap.reversed()
        extend = 'neither'
        if vmin == None:
            vmin = min(values)
        elif vmin > min(values):
            extend = 'min'
        if vmax == None:
            vmax = values.max()
        elif vmax < max(values):
            if extend == 'min':
                extend = 'both'
            else:
                extend = 'max'

        ax = fig.add_subplot(figloc,
        projection=ccrs.Robinson(central_longitude=self.central_longitude))
        ax.set_title(title, fontsize=20)
        ax.set_global()
        ax.coastlines(color='grey')
        self.plot_boundaries(ax)

        scr = ax.scatter(self.lons, self.lats, c=values, alpha=0.6,
                         edgecolors='w', linewidth=0.5, cmap=cmap,
                         vmin=vmin, vmax=vmax, transform=ccrs.PlateCarree())
        clb = plt.colorbar(scr, extend=extend, orientation='vertical',
                           fraction=0.0235, pad=0.05)
        clb.ax.tick_params(labelsize=10)

        if self.highlight:
            index = self.highlight - 1
            ax.scatter(self.lons[index], self.lats[index], s=200,
                       facecolors=None, edgecolors='k', linewidth=3, alpha=0.5,
                       transform=ccrs.PlateCarree())

        if fig_label:
            ax.text(-0.02, 1.05, fig_label,
                    transform=ax.transAxes, fontsize='xx-large')

        if self.write:
            if filename:
                plt.savefig(f'{filename}.png', bbox_inches='tight')

        if show: plt.show()

    def correlation_plot(self, x, y, xlabel, ylabel,
                         residual=False, trend=False, aspect=None,
                         fig=None, figloc=111, fig_label=None, filename=None):

        show = False
        if fig is None:
            show = True
            fig = plt.figure(figsize=(5,5))

        ax = fig.add_subplot(figloc)
        ax.scatter(x, y, s=0.4)
        ax.set_xlabel(xlabel, fontsize=15)
        ax.set_ylabel(ylabel, fontsize=15)

        if trend:
            ax.plot([-100, 100], [-100, 100],
                    c='r', linestyle='--', linewidth=0.5)

        if self.highlight:
            index = self.highlight - 1
            ax.scatter(x[index], y[index], s=200,
                       facecolors=None, edgecolors='k', linewidth=3, alpha=0.5)
            ax.plot([x[index], x[index]], [-1, y[index]], c='k', linestyle='--')
            ax.plot([-1, x[index]], [y[index], y[index]], c='k', linestyle='--')

        side_x = max(x)
        side_y = max(y)
        if aspect == "equal":
            side = max(side_x, side_y)
            axis_limits = [-0.05, 1.05, -0.05, 1.05]
            ax.axis([item * side for item in axis_limits])
        else:
            ax.axis([-0.05*side_x, 1.05*side_x, -0.05*side_y, 1.05*side_y])

        ax.set_box_aspect(1)

        if residual:
            ax.text(-0.45, 0.95, fig_label,
                    transform=ax.transAxes, fontsize='xx-large')
        else:
            ax.text(-0.2, 1.05, fig_label,
                    transform=ax.transAxes, fontsize='xx-large')

        if self.write:
            if filename:
                plt.savefig(f'{filename}.png', bbox_inches='tight')

        if show: plt.show()

    def plot_boundaries(self, ax):

        geo_df = gpd.read_file('Data/PB2002_boundaries.json')
        for geometry in geo_df['geometry']:
            shape_feature = ShapelyFeature([geometry], ccrs.PlateCarree(),
                                           edgecolor='black', facecolor='none')
            ax.add_feature(shape_feature)
