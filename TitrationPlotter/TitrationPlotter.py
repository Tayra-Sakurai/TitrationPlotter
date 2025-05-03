from math import nan
import numpy as np
import scipy.optimize as opt

ca = 0.0
cb = 0.0
va = 0.0
vb = 0.0
kw = 1e-14

def titration_pH (vb: np.ndarray, ka: float) -> np.ndarray:
    '''This returns the array of the roots.

    This function solves the equation about the concentration of proton.

    The meaningd of the paramaters are:

    ka: Ka value of the titrated acid.
    ca: initial concentration of the concentration of the titrated acid.
    va: Initial volume of the acid.
    cb: initial concentration of the concentration of the titrating base.
    vb: volume of added base.
    kw: Kw value at the tempreture.
    '''
    result = list()
    # For each vb, make and solve the equation
    for arr in np.nditer(vb):
        equation = np.poly1d([1, ka + (vb * cb) / (va + vb), ((ka / (va + vb)) * (cb * vb - ca * va)) - kw, - ka * kw])
        roots = equation.roots
        # Find suitable solution
        for root in np.nditer(roots):
            if type(root) == float:
                if root > 0:
                    result.append(root)
                    break
                else:
                    continue

        
    
