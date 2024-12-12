import numpy as np

from base import B_sin
from Diagonalization.diagonalisation import davidson


### Constantes et fonctions
V = input("Entre la fonction du potentiel")
N = int(input("Entre la dimension de ta base"))
L = float(input("Entre la valeur de L"))


### BP
base = B_sin(N, L)
base_de_sinus = base.ret_base()

H = np.zeros((N,N))

for i in range(N):
    for j in range(i + 1):
        H[i,j] = base.prdt_scalaire_Ec(i, j) +  base.prdt_scalaire_Ep(i, j, V)

H = H + H.T - np.diag(np.diag(H))

E, phi = davidson(H)

print(E)