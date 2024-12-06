import numpy as np

import GramSchmidt as GS
from sort import triFusion


###########################################################################
# Documentation : https://www.irisa.fr/sage/bernard/publis/DAVIDSON94.pdf #
#                 https://en.wikipedia.org/wiki/Lanczos_algorithm         #
###########################################################################

def davidson(M, m=1, l=1, seuil=1e-8, MIt=600) -> tuple:
    """
    M : la matrice à diagonaliser
    m : la précision de l'agorithme Km = {V0;MV0;M²V0;...}
    l : le nombre de vecteur/valeur propre à trouver
    seuil : précision de l'aproxiamtion des valeurs propres
    MIt : maximum d'itiération

    Fonction qui retourne le coupe ([valeur propre], [vecteur propre]) de
    la matrice M
    """
    np.set_printoptions(precision=3, suppress=True)

    # Initialisation
    M = np.array(M)
    N = len(M)
    V1 = np.eye(N, l)
    I = np.identity(N)

    if not(l<=m<=N):
        raise ValueError(f"value must be {l} <= {m} <= {N}")

    # Concstruction du sous espace réduit de (Kernal)
    v = [V1]

    # itération
    for k in range(MIt):
        # Résolution du sous-espace réduit
        T = np.dot(v[k].T, np.dot(M, v[k]))
        val, vec = np.linalg.eig(T)
        Eigencouple = list(zip(val,vec))
        #Eigencouple = triFusion(Eigencouple)[-l:] # Plus grande VP
        Eigencouple = triFusion(Eigencouple)[:l] # Plus petite VP
        val, vec = list(zip(*Eigencouple))
        vec = np.array(vec).T
        val = np.array(val)

        x = np.array([[] for _ in range(N)])
        r = np.array([[] for _ in range(N)])
        t = np.array([[] for _ in range(N)])

        # Expansion du sous-espace
        for i in range(l):
            yi = vec[:, i].reshape((T.shape[0], 1))

            xi = np.dot(v[k], yi)
            x = np.concatenate((x, xi), axis=1)

            ri = val[i]*xi-np.dot(M, np.dot(v[k], yi))
            r = np.concatenate((r, ri), axis=1)

            # Seuil à 1e-300 pour eviter les problèmes
            # de division par 0
            mu_i = [1/max(np.abs(val[i]-M[j, j]), 1e-200)
                    for j in range(N)]
            Ci = np.diag(mu_i)
            ti = np.dot(Ci, ri)
            t = np.concatenate((t, ti), axis=1)
        
        #print(f"t : {t}")        
        #print(f"x : {x}\nr : {r}")
                
        # Stop si l'algo converge        
        if np.linalg.norm(t) < seuil:
            return val, x

        # Sinon on regénère une nouvelle matrice de passage
        if v[k].shape[0] <= m - l :
            v.append(GS.gram_schmidt_matrix(
                np.concatenate((v[k], t), axis=1)
                ).reshape(N, v[k].shape[1]+l))
            #print(f"V : {V}\nT : {T}")
        else:
            v.append(GS.gram_schmidt_matrix(
                np.concatenate((x, t), axis=1)
                ).reshape(N, 2*l))
            #print(f"V : {v[k]}\nT : {T}")
    return val, x

n = 1600
print('Dimension of the matrix',n,'*',n)
sparsity = 0.001
A = np.zeros((n,n))
for i in range(0,n) : 
    A[i,i] = i-9
A = A + sparsity*np.random.randn(n,n)
A = (A.T + A)/2

print(davidson(A))
print(np.linalg.eig(A))