import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np

from find_symmetry_groups import norm_matrix
from find_symmetry_groups import proj_to_vsig_of_u_new
from find_symmetry_groups import xyzTP
from find_symmetry_groups import ZRot
from find_symmetry_groups import YRot

def Alpha(t_mat, theta, phi, MONOorXISO):

    U = np.dot(ZRot(theta), YRot(phi))
    proj_matrix = proj_to_vsig_of_u_new(t_mat, U, MONOorXISO)
    return np.arccos(norm_matrix(proj_matrix) / norm_matrix(t_mat)) / np.pi * 180

def ES_ContourPlots(t_mat, sigma, vmin=None, vmax=None):

    assert sigma in ['MONO', 'XISO'], "sigma should either be 'MONO' or 'XISO'"

    # Generating contour plot similar to cpMONOprelim in Mathematica
    theta_values = np.linspace(0, 2 * np.pi, 100)
    phi_values = np.linspace(0, np.pi, 100)
    Theta, Phi = np.meshgrid(theta_values, phi_values)

    Z = np.array([[Alpha(t_mat, theta, phi, sigma) for theta in theta_values] for phi in phi_values])

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
        #plt.show()

        return plt

# SHOULD THE cpMONO FUNCTIONS BE RENAMED?
def cpMONOprelim(t_mat):
    theta, phi = np.mgrid[0:2 * np.pi:100j, 0:np.pi:50j]
    alpha = np.vectorize(lambda t, p: Alpha(t_mat, t, p, 'MONO'))(theta, phi)  # Adjust as necessary

    x, y, z = xyzTP([theta, phi])
    return x, y, z, alpha


# Function to create the 3D plot with Poly3DCollection
def cpMONO(t_mat):
    x, y, z, alpha = cpMONOprelim(t_mat)

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    # Choose a colormap suitable to your preference and the nature of the data
    color_map = plt.get_cmap('jet')
    ax.plot_surface(x, y, z, facecolors=color_map(alpha/ np.max(alpha)), rstride=1, cstride=1,
                    antialiased=True)

    ax.view_init(elev=30, azim=75)  # You can adjust the elevation and azimuth for the desired viewpoint
    ax.set_axis_off()  # Turn off the axis
    plt.show()
