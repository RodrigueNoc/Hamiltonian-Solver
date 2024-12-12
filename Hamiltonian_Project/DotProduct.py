import numpy as np
from scipy.integrate import quad

class B():
    def __init__(self, N : int, func : str, var : str, pot : str):

        self.V = pot
        # create the base
        self.base =[]

        function = func.split(var)
        for n in range(1,N+1):
            f = ''
            for i in range(len(function)):
                if i != len(function) - 1:
                    f += function[i] + str(n)
                else:
                    f += function[i]
            self.base.append(f)

    def prdt_scalaire_Ec(self, i, j):
        def func(x):
            x = x
            return eval(self.base[i]) + T(self.base[j],x)
        return quad(func, -np.inf, np.inf)[0]
    
    def prdt_scalaire_Ep(self, i, j):
        def func(x):
            x = x
            return eval(self.base[i]+'*'+self.V+'*'+self.base[j])
        return quad(func, -np.inf, np.inf)[0]

def deriv(f, x, h = 1e-3):
    x += h
    res = eval(f)
    x -= 2 * h
    res -= eval(f)
    return res/(2*h)

def T(func, x, h = 1e-3):
    return (deriv(func,x+h)-deriv(func,x-h))/(2*h)