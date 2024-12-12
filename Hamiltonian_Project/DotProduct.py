import numpy as np
from scipy.integrate import quad

class B():
    def __init__(self, N : int, func : str, var : str):

        # create the base
        self.base =[]

        function = func.split(var)
        for n in range(N):
            f = ''
            for _ in function:
                f += _ + str(n)
            self.base.append(f[:-1])

    def prdt_scalaire_Ec(self, i, j):
        def func(x):
            x = x
            return eval(self.base[i]) + T(self.base[j],x)
        return quad(func, -np.inf, np.inf)
    
    def prdt_scalaire_Ep(self,i,j,V):
        def func(x):
            x = x
            return eval(self.base[i]+'*'+V+'*'+self.base[j])

def deriv(f, x, h = 1e-3):
    x += h
    res = eval(f)
    x -= 2 * h
    res -= eval(f)
    return res/(2*h)

def T(func, x, h = 1e-3):
    return (deriv(func,x+h)-deriv(func,x-h))/(2*h)