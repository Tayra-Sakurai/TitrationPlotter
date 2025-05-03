from math import nan
import numpy as np
import scipy.optimize as opt
import matplotlib.pyplot as plt
from tkinter import Button, Widget, PhotoImage, Tk, StringVar, BooleanVar
from tkinter.ttk import Frame, Entry, Label, Spinbox, Button

ca = 0.0
cb = 0.0
va = 110.0
vb = 0.0
kw = 1e-14
root = Tk()
root.title('設定')
csvfilename = StringVar()
caS = StringVar()
kwS = StringVar()
cbS = StringVar()

def cmd ():
    '''Button command
    '''
    pass

def settings () -> None:
    '''This sets the options of the plot

    This has no argumanets.
    '''
    global ca,cb,kw,root,caS,cbS,kwS,csvfilename
    content = Frame(root)
    # CSV file of the data
    datainput = Entry(content, textvariable=csvfilename)
    datalabel = Label(content, text='滴定データを含むcsvファイルのパス')
    # concentration of the acid
    caEntry = Spinbox(content, from_=0.00, to=10.000, increment=0.0001, textvariable=caS)
    caLabel = Label(content, text='被滴定溶液の容量モル濃度 / M')
    # Kw
    kwEntry = Entry(content, textvariable=kwS)
    kwLabel = Label(content, text='Kwの値')
    # Cb
    cbEntry = Entry(content, textvariable=cbS)
    cbLabel = Label(content, text='滴定時に加えた塩基の濃度/ mol L-1')
    # Button
    btn = Button(content, text='プロットを作成')
    content.grid(column=0, row=0)
    datainput.grid(column=1, row=0, columnspan=2)
    datalabel.grid(column=0, row=0)
    caEntry.grid(column=1, row=1, columnspan=2)
    caLabel.grid(column=0, row=1)
    cbEntry.grid(column=1, row=2, columnspan=2)
    cbLabel.grid(column=0, row=2)
    kwEntry.grid(column=1, row=3, columnspan=2)
    kwLabel.grid(column=0, row=3)
    btn.grid(column=2, row=4)
    root.mainloop()
    return None

settings()

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
    return np.array(result)
