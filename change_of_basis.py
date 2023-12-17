import numpy as np

from bases import bases

basis_BB = bases("TapeTape")
basis_PHI = bases("Phi")

# Define the inner product function
def inner_product_matrix(M, N):
    return np.dot(M.flatten(), N.flatten())

# Define Dmat
Dmat = np.diag([1, 1, 1, np.sqrt(2), np.sqrt(2), np.sqrt(2)])

# Define MatIsubCAbasisCisON function
def mat_isub_CA_basis_C_is_ON(C, A):
    return np.array([[inner_product_matrix(A[jj], C[ii]) for jj in range(6)] for ii in range(6)])

# Define MatIsubBA function
def mat_isub_BA(A):
    return mat_isub_CA_basis_C_is_ON(basis_BB, A)

# Define MatIsubBΦ and MatIsubΦB
MatIsubBΦ = mat_isub_BA(basis_PHI)
MatIsubΦB = np.linalg.inv(MatIsubBΦ)

# Define PHImatOfCmat function
def phi_mat_of_cmat(Cijkl):
    return Dmat @ Cijkl @ Dmat

# Define CmatOfPHImat function
def cmat_of_phi_mat(PHImat):
    return np.linalg.inv(Dmat) @ PHImat @ np.linalg.inv(Dmat)

# Define TmatOfPHImat function
def tmat_of_phi_mat(PHImat):
    return MatIsubBΦ @ PHImat @ MatIsubΦB

# Define PHImatOfTmat function
def phi_mat_of_tmat(Tmat):
    return MatIsubΦB @ Tmat @ MatIsubBΦ

# Define TmatOfCmat function
def tmat_of_cmat(Cijkl):
    return tmat_of_phi_mat(phi_mat_of_cmat(Cijkl))

# Define CmatOfTmat function
def cmat_of_tmat(Tmat):
    return cmat_of_phi_mat(phi_mat_of_tmat(Tmat))