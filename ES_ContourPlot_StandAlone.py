import numpy as np
import matplotlib.pyplot as plt

def ES_ContourPlot(Tmat, Sigma):

    # Define BB matrices
    BB1 = np.array([[0, 0, 0], [0, 0, 1/np.sqrt(2)], [0, 1/np.sqrt(2), 0]])
    BB2 = np.array([[0, 0, 1/np.sqrt(2)], [0, 0, 0], [1/np.sqrt(2), 0, 0]])
    BB3 = np.array([[0, 1/np.sqrt(2), 0], [1/np.sqrt(2), 0, 0], [0, 0, 0]])
    BB4 = np.array([[-1/np.sqrt(2), 0, 0], [0, 1/np.sqrt(2), 0], [0, 0, 0]])
    BB5 = 1/np.sqrt(6) * np.array([[1, 0, 0], [0, 1, 0], [0, 0, -2]])
    BB6 = 1/np.sqrt(3) * np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1]])

    # List of BB matrices for easier access
    BB = [BB1, BB2, BB3, BB4, BB5, BB6]

    # Function Definitions
    def xyzTP(theta, phi):
        return [np.cos(theta) * np.sin(phi), np.sin(theta) * np.sin(phi), np.cos(phi)]

    def inner_product_matrix(M, N):
        return np.sum(M * N)

    def norm_matrix(M):
        return np.sqrt(inner_product_matrix(M, M))

    def YRot(t):
        return np.array([[np.cos(t), 0, np.sin(t)], [0, 1, 0], [-np.sin(t), 0, np.cos(t)]])

    def ZRot(t):
        return np.array([[np.cos(t), -np.sin(t), 0], [np.sin(t), np.cos(t), 0], [0, 0, 1]])

    def matrixL(L_func):
        return np.array([[inner_product_matrix(L_func(BB[j]), BB[i]) for j in range(6)] for i in range(6)])

    def MatrixUbar(U):
        def inner_func(matrix):
            return np.dot(U, np.dot(matrix, U.T))
        return matrixL(inner_func)

    def ProjToVSigOfU(Tmat, U, Sigma):
        Ubar = MatrixUbar(U)
        proj_matrix = ProjToVSigOfU_initial(np.dot(np.dot(Ubar.T, Tmat), Ubar))
        return np.dot(np.dot(Ubar, proj_matrix), Ubar.T)

    def ProjToVSigOfU_initial(Tmat):
        new_mat = np.zeros((6, 6))
        new_mat[0, 0] = Tmat[0, 0]
        new_mat[0, 1] = Tmat[0, 1]
        new_mat[1, 0] = Tmat[1, 0]
        new_mat[1, 1] = Tmat[1, 1]
        new_mat[2:, 2:] = Tmat[2:, 2:]
        return new_mat


    def Alpha(Tmat, theta, phi, MONOorXISO):
        U = np.dot(ZRot(theta), YRot(phi))
        proj_matrix = ProjToVSigOfU(Tmat, U, MONOorXISO)
        return np.arccos(norm_matrix(proj_matrix) / norm_matrix(Tmat))


    # Generating contour plot similar to cpMONOprelim in Mathematica
    theta_values = np.linspace(0, 2 * np.pi, 100)
    phi_values = np.linspace(0, np.pi, 100)
    Theta, Phi = np.meshgrid(theta_values, phi_values)

    Z = np.array([[Alpha(Tmat, theta, phi, 'MONO') for theta in theta_values] for phi in phi_values])

    # Specify levels explicitly for discrete contour fields. For example:
    levels = 6

    plt.figure(figsize=(12,6))

    plt.contourf(Theta, Phi, Z, levels=levels, alpha=0.7)
    plt.colorbar()

    plt.contour(Theta, Phi, Z, levels=levels, colors='black')

    plt.xlabel('Theta')
    plt.ylabel('Phi')
    plt.title('Contour Plot')
    plt.show()