# construct a Christoffel matrix

import numpy as np

def Gamma_11(n, C):
    n1, n2, n3 = n
    return (C["C11"] * n1**2 + C["C66"] * n2**2 + C["C55"] * n3**2 +
            2 * (C["C16"] * n1 * n2 + C["C56"] * n2 * n3 + C["C15"] * n1 * n3))

def Gamma_22(n, C):
    n1, n2, n3 = n
    return (C["C66"] * n1**2 + C["C22"] * n2**2 + C["C44"] * n3**2 +
            2 * (C["C26"] * n1 * n2 + C["C24"] * n2 * n3 + C["C46"] * n1 * n3))

def Gamma_33(n, C):
    n1, n2, n3 = n
    return (C["C55"] * n1**2 + C["C44"] * n2**2 + C["C33"] * n3**2 +
            2 * (C["C45"] * n1 * n2 + C["C34"] * n2 * n3 + C["C35"] * n1 * n3))

def Gamma_12(n, C):
    n1, n2, n3 = n
    return (C["C16"] * n1**2 + C["C26"] * n2**2 + C["C45"] * n3**2 +
            (C["C12"] + C["C66"]) * n1 * n2 +
            (C["C25"] + C["C46"]) * n2 * n3 +
            (C["C14"] + C["C56"]) * n1 * n3)

def Gamma_13(n, C):
    n1, n2, n3 = n
    return (C["C15"] * n1**2 + C["C46"] * n2**2 + C["C35"] * n3**2 +
            (C["C14"] + C["C56"]) * n1 * n2 +
            (C["C36"] + C["C45"]) * n2 * n3 +
            (C["C13"] + C["C55"]) * n1 * n3)

def Gamma_23(n, C):
    n1, n2, n3 = n
    return (C["C56"] * n1**2 + C["C24"] * n2**2 + C["C34"] * n3**2 +
            (C["C25"] + C["C46"]) * n1 * n2 +
            (C["C23"] + C["C44"]) * n2 * n3 +
            (C["C36"] + C["C45"]) * n1 * n3)

def ChristoffelMatrix(n, C):
    return np.array([
        [Gamma_11(n, C), Gamma_12(n, C), Gamma_13(n, C)],
        [Gamma_12(n, C), Gamma_22(n, C), Gamma_23(n, C)],
        [Gamma_13(n, C), Gamma_23(n, C), Gamma_33(n, C)],
    ])