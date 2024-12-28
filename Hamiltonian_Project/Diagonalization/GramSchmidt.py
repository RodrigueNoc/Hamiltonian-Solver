import numpy as np


def matrix_to_familly(M: list) -> list:
    """Transforme une matrice en une liste de colone de la matrice"""
    L = []
    for i in range(len(M[0])):
        L.append(M[:, i].reshape(len(M), 1))
    return np.array(L)


def ps(a: list, b: list) -> float:
    """Fait le produit scalaire de Rn"""
    sum = 0
    for i in range(len(a)):
        sum += a[i]*b[i]
    return sum


def norm(u: list) -> float:
    """Norme associé au produit scalaire de Rn"""
    return ps(u, u)**(1/2)


def orthonorm_vec(F: list, v: list) -> list:
    """Orthonomalise un vecteur par rapport à une famille de vecteurs"""
    u = v
    for vk in F:
        u = u - ps(vk, v)*vk
    return u/norm(u)


def gram_schmidt_matrix(M: list) -> list:
    """Aplique l'algorithme de gramschmidt pour orthonomaliser la matrice M."""
    F = matrix_to_familly(M)
    return gram_schmidt(F)


def gram_schmidt(F: list) -> list:
    """Aplique l'algorithme de gramschmidt pour orthonomaliser la liste L"""
    ONF = []
    for e in F:
        ONF.append(orthonorm_vec(ONF, e))
    return np.transpose(np.array(ONF))
