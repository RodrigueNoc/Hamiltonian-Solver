### Ce code permet de creer les bases pour calculer l'hamiltonien
from math import pi, sqrt, sin, factorial

class B_sin():
    """ Cette classe permet de calculer la base des sinus pour le calcul
    de l'hamiltonien """
    def __init__(self, N, L, V):
        """ Methode qui initialise la classe de l'objet B_sin() """
        self.N = N
        self.L = L
        self.V = V
        x = sqrt(2)
        self.k = eval(self.V)
    
    def ret_base(self):
        """ Methode qui calcule la liste des sinus qui formera la base """
        return [f'sqrt(2/{self.L}) * sin(({n}*pi*x)/{self.L})' for n in 
                range(self.N)]

    def prdt_scalaire_Ec(self, i, j):
        """ Methode qui donne le produit scalaire des sinus (calcul a la 
        main) """
        if i == j:
            return - (j*pi)**2 / self.L
        else :
            return 0
    
    def prdt_scalaire_Ep(self, i, j):
        """ Method that gives the dot product of sines (hand calculation) 
        for a harmonic potential"""
        if i == j :
            return (self.k * (self.L)**2) / 3
        else :
            return ((2 * self.k * (self.L)**2) / pi**2) * ( (((-1)**(i+j))/(i+j)**2) - (((-1)**(i-j))/(i-j)**2) )




##### /!\ EN TRAVAUX /!\ #####

## Permet de calculer le polynôme d'Hermite grâce à une somme (source : wikipédia de "Polunôme d'Hermite")
# def P_Hermite(n,x):
#     S = 0
#     for k in range(n//2 + 1):
#         S += (-1)**k * (factorial(n)/(2**k * factorial(k) * factorial(n-2*k))) * x**(n-2*k)
#     return S

## Construction de la classe de la base de l'oscillateur harmonique
# # class B_OH():
#     """ Cette classe permet de calculer la base des sinus pour le calcul
#     de l'hamiltonien """
#     def __init__(self, N, m, omega, h_bar):
#         """ Methode qui initialise la classe de l'objet B_sin() """
#         self.N = N
#         self.m = m
#         self.omega = omega
#         self.h_bar = h_bar
    
#     def ret_base(self):
#         """ Methode qui calcule la liste des Phin(x) qui formera la base """
#         return [f'(1/sqrt((2**{n})*factorial({n}))) * (({self.m}*{self.omega})/
#                 (pi*{self.h_bar}))**(1/4) * P_Hermite({n}, sqrt({self.m}*{self.omega}/
#                 {self.h_bar})*x) * exp((-{self.m}*{self.omega}*(x**2))/(2*{self.h_bar}))' 
#                 for n in range(self.N)]

## Ci-dessus 'ret_base' permet de calculer la base grâce à la formule trouvé sur :
## wikipédia de "Oscillateur harmonique quantique" dans la partie "États propres de l'opérateur N"




### Des idées pour la suite ....
# oscillateur harmonique (en partie fait mais probablement à revoir)
# polynome de Legendre (rotateur rigide)

