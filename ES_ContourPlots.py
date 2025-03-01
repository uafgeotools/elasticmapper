import concurrent.futures
import os

import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np
import pyvista as pv
from functools import partial
from matplotlib.colors import ListedColormap

from Scripts.themes import symmetry_classes
from find_symmetry_groups import norm_matrix
from find_symmetry_groups import proj_to_vsig_of_u_new
from find_symmetry_groups import ZRot
from find_symmetry_groups import YRot
from safe_module import closest

def Alpha(t_mat, MONOorXISO, theta, phi):
    U = np.dot(ZRot(theta), YRot(phi))
    proj_matrix = proj_to_vsig_of_u_new(t_mat, U, MONOorXISO)
    return np.arccos(norm_matrix(proj_matrix) / norm_matrix(t_mat)) / np.pi * 180

def ES_ContourPlots(t_mat, sigma, vmin=None, vmax=None):

    assert sigma in ['MONO', 'XISO'], "sigma should either be 'MONO' or 'XISO'"

    # Generating contour plot similar to cpMONOprelim in Mathematica
    theta_values = np.linspace(0, 2 * np.pi, 100)
    phi_values = np.linspace(0, np.pi, 100)
    Theta, Phi = np.meshgrid(theta_values, phi_values)

    Z = np.array([[Alpha(t_mat, sigma, theta, phi) for theta in
                   theta_values] for phi
                  in phi_values])

    if vmin is not None and vmax is not None:
        fig, ax = plt.subplots(figsize=(20, 10))
        levels = np.linspace(vmin, vmax, 9)
        cmap = plt.cm.get_cmap('viridis', len(levels)-1)

        contourf = ax.contourf(Theta, Phi, Z, levels=levels, cmap=cmap)
        ax.contour(Theta, Phi, Z, levels=levels, colors='black')

        cbar = fig.colorbar(contourf, ax=ax, ticks=levels, boundaries=levels)
        cbar.formatter = ticker.FormatStrFormatter("%.2f")
        cbar.update_ticks()

        # Adjusting tick labels to be in terms of π
        ax.set_xticks([0, np.pi / 2, np.pi, 3 * np.pi / 2, 2 * np.pi])
        ax.set_xticklabels(['0', '90°', '180°', '270°', '360°'], fontsize=20)
        ax.set_yticks([0, np.pi / 2, np.pi])
        ax.set_yticklabels(['0', '90°', '180°'], fontsize=20)
        ax.set_xlabel('azimuthal angle', fontsize=20)
        ax.set_ylabel('polar angle', fontsize=20)
        # plt.show()

        return ax

    else:
        plt.figure(figsize=(12, 6))
        levels = 6
        plt.contourf(Theta, Phi, Z, levels=levels, alpha=0.7)
        plt.colorbar()
        plt.contour(Theta, Phi, Z, levels=levels, colors='black')
        plt.xlabel('Theta')
        plt.ylabel('Phi')
        plt.title('Contour Plot')
        plt.show()

        return plt

def alpha_mono_spheres(t_mats, lattice=False, pov=None, up_direction=None):

    workers = os.cpu_count()

    n_intervals = 10
    cmap = ListedColormap(plt.get_cmap("RdYlBu", n_intervals)
                          (np.linspace(0, 1, n_intervals)))
    cmap = cmap.reversed()

    # Generate spherical coordinates
    theta_values = np.linspace(0, 2 * np.pi, 500)
    phi_values = np.linspace(0, np.pi, 500)
    theta, phi = np.meshgrid(theta_values, phi_values)

    # Convert spherical coordinates to Cartesian coordinates
    x = np.sin(phi) * np.cos(theta)
    y = np.sin(phi) * np.sin(theta)
    z = np.cos(phi)

    x = x.flatten()
    y = y.flatten()
    z = z.flatten()

    # Create a PyVista mesh
    points = np.column_stack((x, y, z))

    if not isinstance(t_mats, list):
        t_mats = [t_mats]

    if lattice:
        t_mat = t_mats[0]
        t_mats = [closest(t_mat, symmetry_class) for symmetry_class in
                  symmetry_classes]
        t_mats.append(t_mat)

    n_t_mat = len(t_mats)
    if n_t_mat == 1:
        plotter = pv.Plotter(shape=(1, 1))
    else:
        plotter = pv.Plotter(shape=((n_t_mat + 1) // 2, 2))

    for i, t_mat in enumerate(t_mats):
        print("t_mat index = ", i+1)
        print("running ....")

        Alpha_new = partial(Alpha, t_mat, "MONO")
        with concurrent.futures.ProcessPoolExecutor(max_workers=workers) as executor:
            data = np.array(list(executor.map(Alpha_new, theta.flatten(), phi.flatten())))

        mesh = pv.PolyData(points)
        mesh['data'] = data

        plotter.subplot(i // 2, i % 2)
        plotter.add_mesh(mesh, scalars='data', cmap=cmap, opacity=1.0,
                         show_scalar_bar=False)
        plotter.add_axes()
        plotter.add_scalar_bar(title=f"alpha_mono", vertical=True,
                               n_labels=n_intervals)
        if pov is not None and up_direction is not None:
            plotter.view_vector(pov, viewup=up_direction)

    plotter.show()