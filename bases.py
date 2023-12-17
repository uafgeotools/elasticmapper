import numpy as np
import sys

def bases(basis):

    '''
    Function to return basis for strain and stress

    input:
    basis (int) - Index to select a basis
    1 - Voigt's strain basis
    2 - Voigt's stress basis
    3 - Voigt's normalized stress basis
    4 - Tape & Tape basis

    returns:
    B1, B2, B3, B4, B5, B6 (numpy array) - Six basis elements for the chosen basis
    '''

    s2 = np.sqrt(2)
    s3 = np.sqrt(3)
    s6 = np.sqrt(6)

    if basis == "Voigt_Strain":

        # Voigt's strain basis

        B1 =        np.array([[ 1, 0, 0],
                              [ 0, 0, 0],
                              [ 0, 0, 0]])

        B2 =        np.array([[ 0, 0, 0],
                              [ 0, 1, 0],
                              [ 0, 0, 0]])

        B3 =        np.array([[ 0, 0, 0],
                              [ 0, 0, 0],
                              [ 0, 0, 1]])

        B4 = 1/2  * np.array([[ 0, 0, 0],
                              [ 0, 0, 1],
                              [ 0, 1, 0]])

        B5 = 1/2  * np.array([[ 0, 0, 1],
                              [ 0, 0, 0],
                              [ 1, 0, 0]])

        B6 = 1/2  * np.array([[ 0, 1, 0],
                              [ 1, 0, 0],
                              [ 0, 0, 0]])

    elif basis == "Voigt_Stress":

        # Voigt's stress basis

        B1 =        np.array([[ 1, 0, 0],
                              [ 0, 0, 0],
                              [ 0, 0, 0]])

        B2 =        np.array([[ 0, 0, 0],
                              [ 0, 1, 0],
                              [ 0, 0, 0]])

        B3 =        np.array([[ 0, 0, 0],
                              [ 0, 0, 0],
                              [ 0, 0, 1]])

        B4 =        np.array([[ 0, 0, 0],
                              [ 0, 0, 1],
                              [ 0, 1, 0]])

        B5 =        np.array([[ 0, 0, 1],
                              [ 0, 0, 0],
                              [ 1, 0, 0]])

        B6 =        np.array([[ 0, 1, 0],
                              [ 1, 0, 0],
                              [ 0, 0, 0]])

    elif basis == "Phi":

        # Voigt's normalized stress basis

        B1 =        np.array([[ 1, 0, 0],
                              [ 0, 0, 0],
                              [ 0, 0, 0]])

        B2 =        np.array([[ 0, 0, 0],
                              [ 0, 1, 0],
                              [ 0, 0, 0]])

        B3 =        np.array([[ 0, 0, 0],
                              [ 0, 0, 0],
                              [ 0, 0, 1]])

        B4 = 1/s2 * np.array([[ 0, 0, 0],
                              [ 0, 0, 1],
                              [ 0, 1, 0]])

        B5 = 1/s2 * np.array([[ 0, 0, 1],
                              [ 0, 0, 0],
                              [ 1, 0, 0]])

        B6 = 1/s2 * np.array([[ 0, 1, 0],
                              [ 1, 0, 0],
                              [ 0, 0, 0]])

    elif basis == "TapeTape":

        # Tape & Tape basis

        B1 = 1/s2 * np.array([[ 0, 0, 0],
                              [ 0, 0, 1],
                              [ 0, 1, 0]])

        B2 = 1/s2 * np.array([[ 0, 0, 1],
                              [ 0, 0, 0],
                              [ 1, 0, 0]])

        B3 = 1/s2 * np.array([[ 0, 1, 0],
                              [ 1, 0, 0],
                              [ 0, 0, 0]])

        B4 = 1/s2 * np.array([[-1, 0, 0],
                              [ 0, 1, 0],
                              [ 0, 0, 0]])

        B5 = 1/s6 * np.array([[ 1, 0, 0],
                              [ 0, 1, 0],
                              [ 0, 0,-2]])

        B6 = 1/s3 * np.array([[ 1, 0, 0],
                              [ 0, 1, 0],
                              [ 0, 0, 1]])

    else:

        sys.exit("Error: requested basis index invalid")

    return B1,B2,B3,B4,B5,B6