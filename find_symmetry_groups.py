import numpy as np

from bases import bases

basis_BB = bases("TapeTape")

ChopT = 0.01 * np.pi / 180  # 0.01 Degree converted to radians

def chop_beta(beta):
    return np.where(abs(beta) < ChopT, 0, beta)

def xyzTP(theta_phi):
    theta, phi = theta_phi
    return [np.cos(theta) * np.sin(phi), np.sin(theta) * np.sin(phi), np.cos(phi)]

def inner_product_matrix(M, N):
    return np.dot(M.flatten(), N.flatten())

def norm_matrix(M):
    return np.sqrt(inner_product_matrix(M, M))

def YRot(t):
    return np.array([
        [np.cos(t), 0, np.sin(t)],
        [0, 1, 0],
        [-np.sin(t), 0, np.cos(t)]
    ])

def ZRot(t):
    return np.array([
        [np.cos(t), -np.sin(t), 0],
        [np.sin(t), np.cos(t), 0],
        [0, 0, 1]
    ])

def UsHat(theta_sigma_phi):
    theta, sigma, phi = theta_sigma_phi
    return np.dot(ZRot(theta), np.dot(YRot(phi), ZRot(sigma)))

id = np.eye(3)

def RotMatrixQ(w_xyz):
    w, xyz = w_xyz
    x, y, z = xyz
    return np.array([
        [w**2 + x**2 - y**2 - z**2, 2 * (x*y - w*z), 2 * (x*z + w*y)],
        [2 * (x*y + w*z), w**2 - x**2 + y**2 - z**2, 2 * (y*z - w*x)],
        [2 * (x*z - w*y), 2 * (y*z + w*x), w**2 - x**2 - y**2 + z**2]
    ])

def random_rotation():
    ttt = np.random.uniform(0, 2 * np.pi)
    ppp = np.arccos(np.random.uniform(-1, 1))
    sss = np.random.uniform(0, 2 * np.pi)
    rot_quat = [np.cos(sss / 2), np.sin(sss / 2) * np.dot(ZRot(ttt), np.dot(YRot(ppp), [0, 0, 1]))]
    return np.dot(RotMatrixQ(rot_quat), np.dot(ZRot(ttt), YRot(ppp)))

def matrixL(L):
    return [[inner_product_matrix(L(basis_BB[j]), basis_BB[i]) for j in range(6)] for i in range(6)]

def MatrixUbar(U):
    def func(x):
        return np.dot(U, np.dot(x, U.T))
    return matrixL(func)

def amat(Tmat, U):
    return np.dot(np.dot(np.array(MatrixUbar(U)).T, Tmat), MatrixUbar(U))

def DistToVSigmaofU(Tmat, U, Sigma):
    # This is d(T,Subscript[V, Sigma](U)) in the TapeTape2022 paper
    diff_matrix = amat(Tmat, U) - proj_to_vsig_of_u(amat(Tmat, U), id, Sigma)
    return norm_matrix(diff_matrix)

def proj_to_vsig_of_u_new(Tmat, U, Sigma):
    return np.dot(np.dot(MatrixUbar(U), proj_to_vsig_of_u(amat(Tmat, U), id, Sigma)), np.array(MatrixUbar(U)).T)

def proj_to_vsig_of_u(Tmat, id, mode):
    # Check if Tmat is a valid 6x6 matrix
    if not (isinstance(Tmat, np.ndarray) and Tmat.shape == (6, 6)):
        raise ValueError("Tmat must be a 6x6 numpy ndarray")

    # Extract values from Tmat
    a, g, m, q, t, v = Tmat[0, :]
    b, h, n, r, u = Tmat[1, 1:]
    c, i, o, s = Tmat[2, 2:]
    d, j, p = Tmat[3, 3:]
    e, k = Tmat[4, 4:]
    f = Tmat[5, 5]

    if mode == "MONO":
        return np.array([
            [a, g, 0, 0, 0, 0],
            [g, b, 0, 0, 0, 0],
            [0, 0, c, i, o, s],
            [0, 0, i, d, j, p],
            [0, 0, o, j, e, k],
            [0, 0, s, p, k, f]
        ])
    elif mode == "ORTH":
        return np.array([
            [a, 0, 0, 0, 0, 0],
            [0, b, 0, 0, 0, 0],
            [0, 0, c, 0, 0, 0],
            [0, 0, 0, d, j, p],
            [0, 0, 0, j, e, k],
            [0, 0, 0, p, k, f]
        ])
    elif mode == "TRIG":
        return np.array([
            [(a+b)/2, 0, (m+n)/2, 0, 0, 0],
            [0, (a+b)/2, 0, (m+n)/2, 0, 0],
            [(m+n)/2, 0, (c+d)/2, 0, 0, 0],
            [0, (m+n)/2, 0, (c+d)/2, 0, 0],
            [0, 0, 0, 0, e, k],
            [0, 0, 0, 0, k, f]
        ])
    elif mode == "TET":
        return np.array([
            [(a+b)/2, 0, 0, 0, 0, 0],
            [0, (a+b)/2, 0, 0, 0, 0],
            [0, 0, c, 0, 0, 0],
            [0, 0, 0, d, 0, 0],
            [0, 0, 0, 0, e, k],
            [0, 0, 0, 0, k, f]
        ])
    elif mode == "CUBE":
        avg = (a+b+c)/3
        avg_de = (d+e)/2
        return np.array([
            [avg, 0, 0, 0, 0, 0],
            [0, avg, 0, 0, 0, 0],
            [0, 0, avg, 0, 0, 0],
            [0, 0, 0, avg_de, 0, 0],
            [0, 0, 0, 0, avg_de, 0],
            [0, 0, 0, 0, 0, f]
        ])
    elif mode == "XISO":
        return np.array([
            [(a+b)/2, 0, 0, 0, 0, 0],
            [0, (a+b)/2, 0, 0, 0, 0],
            [0, 0, (c+d)/2, 0, 0, 0],
            [0, 0, 0, (c+d)/2, 0, 0],
            [0, 0, 0, 0, e, k],
            [0, 0, 0, 0, k, f]
        ])
    elif mode == "ISO":
        avg_all = (a+b+c+d+e)/5
        return np.array([
            [avg_all, 0, 0, 0, 0, 0],
            [0, avg_all, 0, 0, 0, 0],
            [0, 0, avg_all, 0, 0, 0],
            [0, 0, 0, avg_all, 0, 0],
            [0, 0, 0, 0, avg_all, 0],
            [0, 0, 0, 0, 0, f]
        ])
    else:
        raise ValueError(f"Unsupported mode: {mode}")

def dT(Tmat, sigma):
    return 999

def betaT(Tmat, sigma):
    if sigma == "TRIV":
        return 0
    else:
        return 999  # Assuming this value is in radians as per our previous conversations.

def UT(Tmat, sigma):
    id = np.eye(3)  # 3x3 identity matrix

    if sigma in ["TRIV", "ISO"]:
        return id
    # If you have more cases for sigma in the future, add them here.
    else:
        raise ValueError(f"Unsupported sigma value: {sigma}")

def betaT(Tmat, x):
    return np.arcsin(x / np.linalg.norm(Tmat))*(180/np.pi)

def objective(x, Tmat=None, Sigma=None):
    theta, sigma, phi = x
    return DistToVSigmaofU(Tmat, UsHat([theta, sigma, phi]), Sigma)