import numpy as np

def Cvec2Tvec():

    r2 = np.sqrt(2)
    r3 = np.sqrt(3)

    P2 = np.diag([2, 2, 2, 1, 1/r3, r2/r3, 2, 2, 1, 1/r3, r2/r3, 2, 1, 1/r3, r2/r3, 1/2, 1/2/r3, 1/r2/r3, 1/6,
                   1/r2/3, 1/3])

    M2 = np.array([[   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   1,   0,   0,   0,   0,   0],
                   [   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   1,   0,   0,   0,   0],
                   [   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   1,   0,   0,   0],
                   [   0,   0,   0,  -1,   0,   0,   0,   0,   1,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0],
                   [   0,   0,   0,   1,   0,   0,   0,   0,   1,   0,   0,   0,  -2,   0,   0,   0,   0,   0,   0,   0,   0],
                   [   0,   0,   0,   1,   0,   0,   0,   0,   1,   0,   0,   0,   1,   0,   0,   0,   0,   0,   0,   0,   0],
                   [   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   1,   0,   0],
                   [   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   1,   0],
                   [   0,   0,   0,   0,  -1,   0,   0,   0,   0,   1,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0],
                   [   0,   0,   0,   0,   1,   0,   0,   0,   0,   1,   0,   0,   0,  -2,   0,   0,   0,   0,   0,   0,   0],
                   [   0,   0,   0,   0,   1,   0,   0,   0,   0,   1,   0,   0,   0,   1,   0,   0,   0,   0,   0,   0,   0],
                   [   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   1],
                   [   0,   0,   0,   0,   0,  -1,   0,   0,   0,   0,   1,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0],
                   [   0,   0,   0,   0,   0,   1,   0,   0,   0,   0,   1,   0,   0,   0,  -2,   0,   0,   0,   0,   0,   0],
                   [   0,   0,   0,   0,   0,   1,   0,   0,   0,   0,   1,   0,   0,   0,   1,   0,   0,   0,   0,   0,   0],
                   [   1,  -2,   0,   0,   0,   0,   1,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0],
                   [  -1,   0,   2,   0,   0,   0,   1,  -2,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0],
                   [  -1,   0,  -1,   0,   0,   0,   1,   1,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0],
                   [   1,   2,  -4,   0,   0,   0,   1,  -4,   0,   0,   0,   4,   0,   0,   0,   0,   0,   0,   0,   0,   0],
                   [   1,   2,  -1,   0,   0,   0,   1,  -1,   0,   0,   0,  -2,   0,   0,   0,   0,   0,   0,   0,   0,   0],
                   [   1,   2,   2,   0,   0,   0,   1,   2,   0,   0,   0,   1,   0,   0,   0,   0,   0,   0,   0,   0,   0]])

    return P2 @ M2

def TmatOfCmat(Cmat):
    return v2sm(Cvec2Tvec() @ sm2v(Cmat))

def CmatOfTmat(Tmat):
    return v2sm(np.inv(Cvec2Tvec()) @ sm2v(Tmat))

def v2sm(input_vec):

    '''
    vector to symmetric matrix

    :param input_vec:
    :return:
    '''

    assert len(input_vec) == 21, "The vector length is not 21"

    output_mat = np.zeros((6,6))

    output_mat[0,0:6] = input_vec[0:6]
    output_mat[1,1:6] = input_vec[6:11]
    output_mat[2,2:6] = input_vec[11:15]
    output_mat[3,3:6] = input_vec[15:18]
    output_mat[4,4:6] = input_vec[18:20]
    output_mat[5,5:6] = input_vec[20:21]

    output_mat = output_mat + output_mat.T - np.diag(np.diag(output_mat))

    return output_mat

def sm2v(input_mat):

    '''
    symmetric matrix to vector

    :param input_mat:
    :return:
    '''

    assert input_mat.shape == (6,6), "The matrix shape is not 6x6"

    output_vec = np.zeros(21)

    output_vec[0:6]   = input_mat[0,0:6]
    output_vec[6:11]  = input_mat[1,1:6]
    output_vec[11:15] = input_mat[2,2:6]
    output_vec[15:18] = input_mat[3,3:6]
    output_vec[18:20] = input_mat[4,4:6]
    output_vec[20:21] = input_mat[5,5:6]

    return output_vec