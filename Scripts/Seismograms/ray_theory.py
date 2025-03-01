import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import minimize

def voigt_index(i, j):
    """
    Maps tensor indices (i, j) to Voigt notation index.

    Parameters:
        i (int): Row index (0-based).
        j (int): Column index (0-based).

    Returns:
        int: Voigt notation index (0-based).
    """
    if i > j:  # Ensure the indices are in upper triangular order (i <= j)
        i, j = j, i

    mapping = {
        (0, 0): 0,  # sigma_11
        (1, 1): 1,  # sigma_22
        (2, 2): 2,  # sigma_33
        (1, 2): 3,  # sigma_23
        (0, 2): 4,  # sigma_13
        (0, 1): 5  # sigma_12
    }

    return mapping.get((i, j), None)

def get_phase_group_velocity(a, n):
    """
    Params:
    =================================
    a: np.ndarray, shape(3,3,3,3)
        density normalized tensor: a = cijkl / rho
    n: np.ndarray, shape(3)
        direction

    Returns:
    ========================
    c: np.ndarray shape(3)
        phase velocity, ascending order
    u: np.ndarray(3,3)
        group velocity, shape(3,3) u[i] -> group veloc for c[i]
    """

    # christoffel matrix
    gamma = np.einsum("ijkl,j,l-> ik", a, n, n)  # aijkl n_j n_l

    # find eigenvalues
    eval, evec = np.linalg.eig(gamma)
    c = np.sqrt(eval)
    idx = np.argsort(abs(c))
    evec = evec[:, idx]
    c = c[idx]

    # group velocity
    u = np.zeros((3, 3))
    for i in range(3):
        u[i, :] = np.einsum("ijkl,l,j,k->i", a, n / c[i], evec[:, i],
                            evec[:, i])  # aijkl n_l/c g_j g_k

    return c, u


def get_spherical_angles(xs, xr):
    x, y, z = xr - xs

    # Compute the radius
    r = np.sqrt(x ** 2 + y ** 2 + z ** 2)

    # Compute the polar angle (theta)
    # theta = np.arctan2(z,np.sqrt(x**2+y**2))
    theta = np.arctan2(np.sqrt(x ** 2 + y ** 2), z)

    # Compute the azimuthal angle (phi)
    phi = np.arctan2(y, x)

    return theta, phi

def carl_c66():

    c66 = np.array([
        [119.288, 39.1133, 18.7898, 0.997705, -6.43441, 0.315116],
        [39.1133, 120.655, 16.9447, 2.65705, -2.7794, 0.328052],
        [18.7898, 16.9447, 178.161, -5.48014, 13.8157, -0.868537],
        [0.997705, 2.65705, -5.48014, 25.192, -0.911102, -1.91025],
        [-6.43441, -2.7794, 13.8157, -0.911102, 27.1276, 0.65295],
        [0.315116, 0.328052, -0.868537, -1.91025, 0.65295, 40.4283]
    ])

    cijkl = np.zeros((3, 3, 3, 3))
    for i in range(3):
        for j in range(3):
            i1 = voigt_index(i, j)
            for p in range(3):
                for q in range(3):
                    j1 = voigt_index(p, q)
                    cijkl[i, j, p, q] = c66[i1, j1]

    rho = 2.623
    return rho, cijkl

def plot_veloc_sheet_in_plane(phi_fix, rho, cijkl):

    # find group velocity for each theta
    theta = np.linspace(0,np.pi*2,200)
    phi = np.linspace(phi_fix,phi_fix,1)
    u = np.zeros((len(theta),3,3))
    slow = np.zeros((len(theta),3,3))
    for it in range(len(theta)):
        for ip in range(1):
            n = np.array([np.cos(phi[ip]) * np.sin(theta[it]),np.sin(phi[ip]) * np.sin(theta[it]),np.cos(theta[it])])
            c0,u0 = get_phase_group_velocity(cijkl / rho,n)
            for i in range(3):
                slow[it,i,:] =  1./ c0[i] * n
            u[it,...] = u0[...]

    # plot
    fig = plt.figure(2,figsize=(12,6))
    ax = fig.add_subplot(121)
    for i in range(3):
        unorm = np.linalg.norm(u[:,i,:],axis=-1)
        ax.plot(u[:,i,0],u[:,i,2])
    ax.set_title("Group Velocity Surface")
    ax.set_xlabel("U_1")
    ax.set_ylabel("U_3")

    ax1 = fig.add_subplot(122)
    for i in range(3):
        # v = 1. / np.sqrt(slow[:,i,0]**2 + slow[:,i,2]**2)
        # ang = np.arctan2(slow[:,i,0], slow[:,i,2])
        # vx = v * np.sin(ang)
        # vy = v * np.cos(ang)
        # ax1.plot(-vx,vy)
        ax1.plot(slow[:,i,0],slow[:,i,2])
    ax1.set_title("Slowness Surface")
    ax1.set_xlabel("p_1")
    ax1.set_ylabel("p_3")

def myfunc(x,aijkl,theta0,phi0,i_egn):
    """
    Parameters:
    ================
    x: np.ndarray, shape(2)
        initial guess
    theta0/phi0: float
        target theta0/phi0
    i_egn: int
        which eigenvalue is used, 0-2
    """
    th,ph = x

    n0 = np.array([np.cos(ph) * np.sin(th),np.sin(ph) * np.sin(th),np.cos(th)])
    _,u = get_phase_group_velocity(aijkl,n0)
    theta_u,phi_u = get_spherical_angles(u[i_egn,:] * 0,u[i_egn,:])

    return (theta_u - theta0)**2 + (phi_u - phi0)**2,theta_u,phi_u,np.linalg.norm(u[i_egn,:])

###############################################################################

def group_arrival_times(Cij, rho, xs, xr):

    # Aakash's cijkl/rho/source/station

    # Cij = Cij / 10**9
    # rho = rho / 1000

    cijkl = np.zeros((3, 3, 3, 3))
    for i in range(3):
        for j in range(3):
            i1 = voigt_index(i, j)
            for p in range(3):
                for q in range(3):
                    j1 = voigt_index(p, q)
                    cijkl[i, j, p, q] = Cij[i1, j1]

###############################################################################

    #plot_veloc_sheet_in_plane(1.1341691669813554, rho, cijkl)

###############################################################################

    r = np.linalg.norm(xr-xs)
    theta0,phi0 = get_spherical_angles(xs,xr)
    print("theta = {} phi = {} distance = {}".format(np.rad2deg(theta0),np.rad2deg(phi0),r))

    # compute phase/group velocity
    n = np.array([np.cos(phi0) * np.sin(theta0),np.sin(phi0) * np.sin(theta0),np.cos(theta0)])
    c,u = get_phase_group_velocity(cijkl / rho,n)
    c_time = np.zeros((3))
    for i in range(3):
        u0 = np.linalg.norm(u[i,:])
        theta_u,phi_u = get_spherical_angles(xs * 0,u[i,:])
        theta_u_deg = np.rad2deg(theta_u)
        phi_u_deg = np.rad2deg(phi_u)
        print(f'{i}-th phase time = {r / c[i]}')
        print(f'{i}-th phase/group_veloc = {c[i]} {u0}, theta_u, phi_u = {theta_u_deg} {phi_u_deg}')

        # time predicted by phase velocity
        c_time[i] = r / c[i]

    ###############################################################################

    np.random.seed(100)
    n = 200
    out_u = []
    out_ang = []
    out_iegn = []
    theta = np.random.rand(n) * np.pi
    phi = np.random.rand(n) * np.pi * 2

    for it in range(n): # loop every intial guess
        for i in range(3): # loop each eigenvector
            f = lambda x: myfunc(x,cijkl/rho,theta0,phi0,i)[0]
            msg = minimize(f,[theta[it],phi[it]],method='L-BFGS-B',bounds=[(0,np.pi),(0,np.pi*2)])

            if msg.success and msg.fun < 0.017:
                out_u.append(myfunc(msg.x,cijkl/rho,theta0,phi0,i)[-1])
                out_ang.append([msg.x[0],msg.x[1]])
                out_iegn.append(i)

    u_time_all = r / np.array(out_u)
    out_ang = np.array(out_ang)
    out_iegn = np.array(out_iegn)
    idx = np.unique(u_time_all.round(decimals=1),return_index=True)[-1]

    # select unique time/angs/eigenvalue index
    u_time = u_time_all[idx]
    angs = out_ang[idx]
    iegn_idx = out_iegn[idx]
    print('travel time by group velocity = ',u_time)
    print('travel time prediced by Christoffel = ',c_time[::-1])

    ###############################################################################
    # %%
    # make sure the group velocity close to theta0,phi0
    print(f"target group velocity direction theta/phi = {np.rad2deg(theta0)} {np.rad2deg(phi0)} ")
    for i in range(len(u_time)):
        _, theta_u, phi_u, _ = myfunc(angs[i, :], cijkl / rho, theta0, phi0, iegn_idx[i])
        print(f"from theta,phi = {np.rad2deg(angs[i, 0])} {np.rad2deg(angs[i, 0])}, the group veloc theta/phi = , {np.rad2deg(theta_u)} {np.rad2deg(phi_u)}")

    return u_time