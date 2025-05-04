from csv import reader
from math import ceil, log10, nan
import numpy as np
import scipy.optimize as opt
from scipy import differentiate
import matplotlib.pyplot as plt
from tkinter import VERTICAL, Button, Listbox, Widget, PhotoImage, Tk, StringVar, BooleanVar
from tkinter.ttk import Frame, Entry, Label, Scrollbar, Spinbox, Button

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
# Result
resultbox = Listbox(content, height=3, width=50)
def titration_pH (vb: np.ndarray | float | np.float64, ka: float, delta1: float, delta2: float) -> np.ndarray | float:
    '''This returns the array of the roots.

    This function solves the equation about the concentration of proton.

    The meaningd of the paramaters are:

    ka: Ka value of the titrated acid.
    vb: volume of added base.
    delta1: correction of concentration of proton
    delta2: correction
    '''
    result = list()
    print('ka=', ka)
    # For each vb, make and solve the equation
    if type(vb) == np.ndarray:
        print('size=', vb.size)
        for arr in np.nditer(vb):
            delta = delta1 * (arr ** 2) + delta2
            polyn = [1, ka + (arr * cb) / (va + arr), ((ka / (va + arr)) * (cb * arr - ca * va)) - kw - delta * k1, - (ka * kw + 2 * delta * k2 + delta * ka * k1), -2 * k2 * ka]
            polyn.reverse()
            f = lambda x: np.polynomial.polynomial.polyval(x, polyn)
            rf = opt.fsolve(f, [0.1])
            print('rf=', rf)
            result.append(rf[0])
    elif type(vb) == float or type(vb) == np.float64:
        # If it is a float, do the same
        arr = vb
        delta = delta1 * (arr ** 2) + delta2
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
    (popt, pcov) = opt.curve_fit(titration_pH, vb, pH, (1.754e-5, 0, 1.31e-5))
    for p in np.nditer(popt):
        resultbox.insert('end', p)
    xrange = np.arange(min(*vb),max(*vb),0.00005)
    # Find the equivalent point
    equivlist = [
        titration_pH(vb[i+1], *popt) - titration_pH(vb[i], *popt) for i in range(len(vb)-1)
        ]
    equiv = vb[np.argmax(equivlist) + 1]
    # Output the equivalent point
    resultbox.insert('end', f'equivalent point: {equiv} mL')
    # Output the pH at the equivalent point
    resultbox.insert(
        'end', f'pH at the equivalent point: {titration_pH(equiv, *popt)}')
    # Find the pH at the equivalent point
    equivpH = titration_pH(equiv, *popt)
    plt.axhline(equivpH, 0, equiv, linestyle='--')
    plt.axvline(equiv, 0, 1, linestyle='--')
    plt.plot(xrange, titration_pH(xrange, *popt))
    plt.plot(vb,pH,'.')
    plt.grid(True, 'both')
    plt.minorticks_on()
    plt.xticks(np.arange(int(min(vb)), ceil(max(vb)), 0.02), minor=True)
    plt.yticks(np.arange(int(min(pH)), ceil(max(pH)), 0.05), minor=False)
    plt.xticks(np.arange(int(min(vb)), ceil(max(vb))), minor=True)
    plt.yticks(np.arange(int(min(pH)), ceil(max(pH))), minor=False)
    plt.ylabel(ylabel.get())
    plt.xlabel(xlabel.get())
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
    resultbox.grid(column=0,row=6,columnspan=3)
    resultbox['yscrollcommand'] = sBar.set
    sBar.grid(column=3, row=6)
    btn.grid(column=2, row=7, columnspan=2)
    root.mainloop()
    return None

settings()