# TitrationPlotter

This program is made to plot pH during the titration with considering the data you provide.

This is built with Python 3.13.1. To use this program, you have to activate it from [Python Official Release](https://www.python.org/downloads/).

The csv file to hand to the program must be in the following format and not have any heading rows. Otherwise the program throws an error and no plot.

```csv
[titrated volume of the base],[pH value]
[titrated volume of the base],[pH value]
...
```

You are free to use this program for any purpose allowed under the laws.

## How It Works

This program uses least-square method on the following function:

```pyhton
from csv import reader
from math import ceil, nan
import numpy as np
import scipy.optimize as opt
import matplotlib.pyplot as plt
from tkinter import VERTICAL, Button, Listbox, Widget, PhotoImage, Tk, StringVar, BooleanVar
from tkinter.ttk import Frame, Entry, Label, Scrollbar, Spinbox, Button
def titration_pH (vb: np.ndarray, ka: float, delta1: float, delta2: float) -> np.ndarray:
    '''This returns the array of the roots.

    This function solves the equation about the concentration of proton.

    The meaningd of the paramaters are:

    ka: Ka value of the titrated acid.
    vb: volume of added base.
    delta1: correction of concentration of proton
    delta2: correction
    '''
    result = list()
    print(ka)
    # For each vb, make and solve the equation
    for arr in np.nditer(vb):
        delta = delta1 * (arr ** 2) + delta2
        polyn = [1, ka + (arr * cb) / (va + arr), ((ka / (va + arr)) * (cb * arr - ca * va)) - kw - delta * k1, - (ka * kw + 2 * delta * k2 + delta * ka * k1), -2 * k2 * ka]
        polyn.reverse()
        f = lambda x: np.polynomial.polynomial.polyval(x, polyn)
        rf = opt.fsolve(f, [0.1])
        result.append(rf[0])
    print('result=',result)
    rarray = -np.log10(np.array(result))
    print(rarray)
    return rarray
```

The second and latter arguments represent fitted variables.
