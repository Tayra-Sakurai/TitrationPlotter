#    This program is a titration plotter. This is TitrationPlotter.py written in Python 3.13.2.
#    Copyright (C) 2025 Tayra Sakurai
#
#    This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or any later version.
#
#    This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License along with this program. If not, see <https://www.gnu.org/licenses/>.
#
from csv import reader
from math import ceil, log10, nan
from sys import version
import numpy as np
import scipy
import scipy.optimize as opt
from scipy import differentiate
import matplotlib.pyplot as plt
from tkinter import VERTICAL, Button, Listbox, Widget, PhotoImage, Tk, StringVar, BooleanVar
from tkinter.ttk import Frame, Entry, Label, Scrollbar, Spinbox, Button
import scipy.version

ca = 0.0
cb = 0.0
va = 110.0
kw = 1e-14
k1 = 4.335e-7
k2 = 1.982e-17
root = Tk()
root.title('設定')
content = Frame(root)
csvfilename = StringVar()
caS = StringVar()
kwS = StringVar()
cbS = StringVar()
ylabel = StringVar()
xlabel = StringVar()
originallabel = StringVar()
curvelabel = StringVar()
# Result
resultbox = Listbox(content, height=3, width=50)
def titration_pH (vb: np.ndarray | float | np.float64, ka: float, delta1: float | np.float64, delta2: float | np.float64, delta3: float | np.float64) -> np.ndarray | float:
    '''This returns the array of the roots.

    This function solves the equation about the concentration of proton.

    The meaningd of the paramaters are:

    ka: Ka value of the titrated acid.
    vb: volume of added base.
    delta1: correction of concentration of proton
    delta2: correction
    delta3: correction
    '''
    result = list()
    print('ka=', ka)
    # For each vb, make and solve the equation
    if type(vb) == np.ndarray:
        print('size=', vb.size)
        # Counter
        count = 0
        for arr in np.nditer(vb):
            count += 1
            delta = delta1 * (arr**2) + delta2 * arr + delta3
            polyn = [1, ka + (arr * cb) / (va + arr), ((ka / (va + arr)) * (cb * arr - ca * va)) - kw - delta * k1, - (ka * kw + 2 * delta * k2 + delta * ka * k1), -2 * k2 * ka]
            polyn.reverse()
            f = lambda x: np.polynomial.polynomial.polyval(x, polyn)
            rf = opt.fsolve(f, [0.1])
            countmsg = f'Processing...    {str(count).rjust(len(str(vb.size)))} / {vb.size}'
            delete = '\b' * len(countmsg)
            print(delete + countmsg, end='')
            result.append(rf[0])
        print()
    elif type(vb) == float or type(vb) == np.float64:
        # If it is a float, do the same
        arr = vb
        delta = delta1 * (arr**2) + delta2 * arr + delta3
        polyn = [1, ka + (arr * cb) / (va + arr), ((ka / (va + arr)) * (cb * arr - ca * va)) - kw - delta * k1, - (ka * kw + 2 * delta * k2 + delta * ka * k1), -2 * k2 * ka]
        polyn.reverse()
        f = lambda x: np.polynomial.polynomial.polyval(x, polyn)
        rf = opt.fsolve(f, [0.1])
        result.append(rf[0])
        return -log10(rf[0])
    else:
        raise TypeError('Invalid type of vb')
    print('result=',result)
    rarray = -np.log10(np.array(result))
    rarray = rarray.reshape(vb.shape)
    print('rarray=', rarray)
    return rarray

def cmd () -> bool:
    '''Button command
    '''
    global root,ca,cb,kw
    # Get values
    # read csv file
    pH = list()
    vb = list()
    with open(csvfilename.get(), 'r') as file:
        r = reader(file)
        for line in r:
            vb.append(float(line[0]))
            pH.append(float(line[1]))
    vb = np.array(vb)
    pH = np.array(pH)
    ca = float(caS.get())
    cb = float(cbS.get())
    kw = float(kwS.get())
    (popt, pcov) = opt.curve_fit(titration_pH, vb, pH, (8.77e-6, 0, 0, 1.31e-5))
    for p in np.nditer(popt):
        resultbox.insert('end', p)
    xrange = np.linspace(int(min(*vb)),ceil(max(*vb)),10000)
    # Find the equivalent point
    # Calculate the derivative
    # Find the points where the derivative is 1
    derivobj = differentiate.derivative(lambda x: titration_pH(x, *popt), np.arange(vb.min(),ceil(vb.max()), 0.01))
    deriv = derivobj.df
    # Find the maximum
    maxderiv = np.argmax(deriv)
    equivalent = np.arange(vb.min(), ceil(vb.max()), 0.01)[maxderiv]
    # Find the pH at the equivalent point
    pH_equivalent = titration_pH(equivalent, *popt)
    # Output the equivalent point
    resultbox.insert('end',
                     f'Equivalent point: {equivalent} mL, pH: {pH_equivalent}')
    # Plot style setting
    plt.style.use('grayscale')
    # Plot the point
    plt.plot(equivalent, pH_equivalent, 'o', label="Equivalent Point")
    # Plot horizontal line
    plt.axhline(pH_equivalent, linestyle='--')
    # Plot vertical line
    plt.axvline(equivalent, linestyle='--')
    plt.plot(xrange, titration_pH(xrange, *popt), linestyle='-', label=curvelabel.get())
    plt.plot(vb,pH,'.', label=originallabel.get())
    plt.grid(True, 'both')
    plt.minorticks_on()
    plt.xticks(np.arange(int(min(vb)), ceil(max(vb)), 0.02), minor=True)
    plt.yticks(np.arange(int(min(pH)), ceil(max(pH)), 0.05), minor=False)
    plt.xticks(np.arange(int(min(vb)), ceil(max(vb))), minor=True)
    plt.yticks(np.arange(int(min(pH)), ceil(max(pH))), minor=False)
    plt.ylabel(ylabel.get())
    plt.xlabel(xlabel.get())
    # Show the legend
    plt.legend()
    plt.show()
    return False

def settings () -> None:
    '''This sets the options of the plot

    This has no argumanets.
    '''
    global root,caS,cbS,kwS,csvfilename,xlabel,ylabel,graphlabel
    # CSV file of the data
    datainput = Entry(content, textvariable=csvfilename)
    datalabel = Label(content, text='滴定データを含むcsvファイルのパス')
    # concentration of the acid
    caEntry = Spinbox(content, from_=0.00, to=10.000, increment=0.000001, textvariable=caS)
    caLabel = Label(content, text='被滴定溶液の容量モル濃度 / M')
    # Kw
    kwEntry = Entry(content, textvariable=kwS)
    kwLabel = Label(content, text='Kwの値')
    # Cb
    cbEntry = Entry(content, textvariable=cbS)
    cbLabel = Label(content, text='滴定時に加えた塩基の濃度/ mol L-1')
    # X label
    xLabelEntry = Entry(content, textvariable=xlabel)
    xLabelLabel = Label(content, text='x軸のラベル')
    # Y label
    yLabelEntry = Entry(content, textvariable=ylabel)
    yLabelLabel = Label(content, text='y軸のラベル')
    #Scroll bar
    sBar = Scrollbar(content, orient=VERTICAL, command=resultbox.yview)
    # Label for points
    pointLabelEntry = Entry(content, textvariable=originallabel)
    pointLabelLabel = Label(content, text='csvでインポートしたデータの凡例')
    # Label for curve
    curveLabelEntry = Entry(content, textvariable=curvelabel)
    curveLabelLabel = Label(content, text='フィッティングしたデータの凡例')
    # Button
    btn = Button(content, text='プロットを作成', command=cmd)
    content.grid(column=0, row=0)
    datainput.grid(column=1, row=0, columnspan=3)
    datalabel.grid(column=0, row=0)
    caEntry.grid(column=1, row=1, columnspan=3)
    caLabel.grid(column=0, row=1)
    cbEntry.grid(column=1, row=2, columnspan=3)
    cbLabel.grid(column=0, row=2)
    kwEntry.grid(column=1, row=3, columnspan=3)
    kwLabel.grid(column=0, row=3)
    xLabelEntry.grid(column=1, row=4, columnspan=3)
    xLabelLabel.grid(column=0, row=4)
    yLabelEntry.grid(column=1, row=5, columnspan=3)
    yLabelLabel.grid(column=0, row=5)
    pointLabelEntry.grid(column=1, row=6, columnspan=3)
    pointLabelLabel.grid(column=0, row=6)
    curveLabelEntry.grid(column=1, row=7, columnspan=3)
    curveLabelLabel.grid(column=0, row=7)
    resultbox.grid(column=0,row=8,columnspan=3)
    resultbox['yscrollcommand'] = sBar.set
    sBar.grid(column=3, row=8)
    btn.grid(column=2, row=9, columnspan=2)
    root.mainloop()
    return None

print('This is Python',version)
print('with NumPy',np.version.full_version)
print('and SciPy',scipy.__version__)
settings()