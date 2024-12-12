import numpy as np

from base import B_sin
from Diagonalization.diagonalisation import davidson
import DotProduct as DP

### BP
def CalculHamiltonian(base):

    H = np.zeros((N,N))

    for i in range(N):
        for j in range(i + 1):
            H[i,j] = base.prdt_scalaire_Ec(i, j) +  base.prdt_scalaire_Ep(i, j)

    H = H + H.T - np.diag(np.diag(H))

    return H

### Constantes et fonctions
m = 1
w = 1
h = 1
k = m*(w**2)
V = f"({k}*x**2)/2"
N = 5
L = 1
#B = B_sin(N, L, V)
B = DP.B(N,f'np.sqrt(2/{L})*np.sin( (a*np.pi*x) /{L})','a', V)

H = CalculHamiltonian(B)
print(H)
E, phi = davidson(H, N)
#E = 0
eig, val = np.linalg.eig(H)
print(np.sort(eig))

print("Energie numerique : ", E)
print("Energie analitique : ", h*w/2)