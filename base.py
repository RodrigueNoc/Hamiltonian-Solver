### Ce code permet de cr√©er les bases pour calculer l'hamiltonien
from math import pi, sqrt, sin

class B_sin():
    """ Cette classe permet de calculer la base des sinus pour le calcul
    de l'hamiltonien """
    def __init__(self,N,L):
        """ MÈthode qui initialise la classe de l'objet B_sin() """
        self.N = N
        self.L = L
    
    def ret_base(self):
        """ MÈthode qui calcule la liste des sinus qui formera la base """
        return [f'sqrt(2/{self.L}) * sin(({n}*pi*x)/{self.L})' for n in 
                range(self.N)]

    def prdt_scalaire_Ec(self,i,j):
        """ MÈthode qui donne le produit scalaire des sinus (calcul ‡ la 
        main) """
        if i == j:
            return - (j*pi)**2 / self.L
        else :
            return 0
    
    def prdt_scalaire_Ep(self,i,j):
        """ MÈthode qui donne le produit scalaire des sinus (calcul ‡ la 
        main) """
        if i == j :
            return
        else :
            return

    


