# -*- coding: utf-8 -*-
"""
Created on Thu Jan 28 16:32:48 2016

@author: jclark

this code uses the Ghosh method to determine the apparent resistivities
for a layered earth model. Either schlumberger or Wenner configurations
can be used
"""
"""
GUI Implementation Additions
Created Winter 2021

@author: Rebecca DiCosola
"""

# Imports
import numpy as np
import random
import matplotlib.pyplot as plt
import sys
from tkinter import *
from tkinter import filedialog

# Constants
# Schlumberger filter
FLTR1 = [0., .00046256, -.0010907, .0017122, -.0020687,
         .0043048, -.0021236, .015995, .017065, .098105, .21918, .64722,
         1.1415, .47819, -3.515, 2.7743, -1.201, .4544, -.19427, .097364,
         -.054099, .031729, -.019109, .011656, -.0071544, .0044042,
         -.002715, .0016749, -.0010335, .00040124]

# Wenner Filter
FLTR2 = [0., .000238935, .00011557, .00017034, .00024935,
         .00036665, .00053753, .0007896, .0011584, .0017008, .0024959,
         .003664, .0053773, .007893, .011583, .016998, .024934, .036558,
         .053507, .078121, .11319, .16192, .22363, .28821, .30276, .15523,
         -.32026, -.53557, .51787, -.196, .054394, -.015747, .0053941,
         -.0021446, .000665125]

# Array size
# 65 is completely arbitrary
ARRAYSIZE = 65

# I know there must be a better method to assign lists. And probably numpy
# arrays would be best. But my Python wasn't up to it. If the last letter
# is an 'l' that means it is a log10 of the value

p = [0] * 20
r = [0] * ARRAYSIZE
rl = [0] * ARRAYSIZE
t = [0] * 50
b = [0] * ARRAYSIZE
asav = [0] * ARRAYSIZE
asavl = [0] * ARRAYSIZE
adatl = [0] * ARRAYSIZE
rdatl = [0] * ARRAYSIZE
adat = [0] * ARRAYSIZE
rdat = [0] * ARRAYSIZE
pkeep = [0] * ARRAYSIZE
rkeep = [0] * ARRAYSIZE
rkeepl = [0] * ARRAYSIZE
pltanswer = [0] * ARRAYSIZE
pltanswerl = [0] * ARRAYSIZE
pltanswerkeep = [0] * ARRAYSIZE
pltanswerkeepl = [0] * ARRAYSIZE

small = [0] * ARRAYSIZE
xlarge = [0] * ARRAYSIZE

x = [0] * 100
y = [0] * 100
y2 = [0] * 100
u = [0] * 5000

# this variable is never used
# new_x = [0] * 1000
# this variable is never used
# new_y = [0] * 1000

# Input
algorithm_choice = 1
SHCHLUMBERGER = 1
WENNER = 2
layers_choice = 0  # number of layers
n = 2 * layers_choice - 1
iter = 10000  # number of iterations for the Monte Carlo guesses. to be input on GUI

# this is where the range in parameters should be input from a GUI
# I'm hard coding this in for now

# enter thickenss range for each layer and then resistivity range.
# for 3 layers small[1] and small[2] are low end of thickness range
# small[3], small[4] and small[5] are the low end of resistivities

# small[1] = 1.
# xlarge[1] = 5
# small[2] = 10.
# xlarge[2] = 75.
# small[3] = 20.
# xlarge[3] = 200.
# small[4] = 2.
# xlarge[4] = 100
# small[5] = 500.
# xlarge[5] = 3000.

# hard coded data input - spacing and apparent resistivities measured
# in the field
ndat = 13
# adat = [0., 0.55, 0.95, 1.5, 2.5, 3., 4.5, 5.5, 9., 12., 20., 30., 70.]
# rdat = [0., 125., 110., 95., 40., 24., 15., 10.5, 8., 6., 6.5, 11., 25.]
one30 = 1.e30
rms = one30
errmin = 1.e10

spac = 0.2  # smallest electrode spacing
m = 20  # number of points where resistivity is calculated

spac = np.log(spac)
delx = np.log(10.0) / 6.

# these lines apparently find the computer precision ep
ep = 1.0
ep = ep / 2.0
fctr = ep + 1.
while fctr > 1.:
    ep = ep / 2.0
    fctr = ep + 1.

# GUI initialization
mainwindow = Tk()
mainwindow.title('VES Inverse Monte Carlo')
mainwindow.configure(background="gainsboro")

frame = Frame(mainwindow)
frame.pack()

# variables used for GUI
algorithm_index = IntVar(mainwindow, 1)
num_layers = IntVar(mainwindow, 3)
num_iter = IntVar(mainwindow, 10000)
# num_datapoints = IntVar(mainwindow)
# resistivity_file = StringVar()

# function definitions
def openGUI():
    global algorithm_choice
    global layers_choice
    global n
    global iter
    global ndat
    global file_explore
    global file_view
    global resistivity_file
    global file_content
    global algorithm_index
    global rdat
    global adat
    global ndat

    # full window label
    # main_label = Label(mainwindow, bg="gainsboro", font=("TkDefaultFont", 15),
    #               text="Input Data for Resistivity Survey")
    # main_label.pack(side=TOP, anchor=NW)

    # file explore button
    preframe = Frame(mainwindow, background="gainsboro")
    preframe.pack(side=TOP, anchor=NW)
    selected_file = Label(preframe, bg="gainsboro", text = "Selected File Path")
    selected_file.grid(row=1, column=2)
    file_view = Label(preframe, bg="gainsboro", text = "No file", width=40, wraplength=220, justify="center")
    file_view.grid(row=2, column=2)
    file_explore = Button(preframe, text = "Select Resistivity Data File",
                        command = pickFile)
    file_explore.grid(row=1, column=1, rowspan=2)

    # drop down menu to pick number of layers
    dropdown_label = Label(preframe, bg="gainsboro",
                  text="Number of Layers", width=20)
    dropdown_label.grid(row=1, column=3)
    layerlist = [
        1, 2, 3, 4, 5
    ]
    layersmenu = OptionMenu(preframe, num_layers, *layerlist)
    layersmenu.config(bg="gainsboro")
    layersmenu.grid(row=2, column=3)

    # box to enter number of iterations
    iter_label = Label(preframe, bg="gainsboro",
                  text="Number of Iterations", width=20)
    iter_label.grid(row=1, column=4)
    iterentry = Entry(preframe, textvariable=num_iter, width=10)
    iterentry.grid(row=2, column=4, pady=5)

    # button to add details below for num layers
    layersframe = Frame(mainwindow, background="gainsboro")
    layersframe.pack(side=TOP, anchor=NW)
    layer_details = Button(layersframe, text = "Add layer details",
                        command = layerDetails)
    layer_details.grid(row=1, columnspan=7, pady=5)

    # execute VEScurves button
    executionframe = Frame(mainwindow, background="gainsboro")
    executionframe.pack(side=BOTTOM, anchor=SW)
    execute_VES = Button(executionframe, text = "Execute VEScurves",
                        command = executeVES)
    execute_VES.grid(row=1, column=1, pady=5)

    # view predicted model button
    view_model = Button(executionframe, text = "View Predicted Model",
                        command = executeVES)
    view_model.grid(row=1, column=2, pady=5)

    # plot curves button
    plot_curves = Button(executionframe, text = "Plot the Curves",
                        command = executeVES)
    plot_curves.grid(row=1, column=3, pady=5)

    return

def pickFile():
    global file_explore
    global file_view
    global resistivity_file
    global file_content
    global algorithm_index
    global file_list
    global algorithm_choice
    global adat
    global rdat
    global ndat

    # get file
    resistivity_file = filedialog.askopenfilename(initialdir="/",
                        title="Open File",
                        filetypes=(("Text Files", "*.txt"),
                        ("All Files", "*.*")))

    # set label by file button to resistivity file link
    file_view.config(text=resistivity_file)

    # file_content = open(resistivity_file, 'r')
    # file_list = file_content.split("\n")
    # print(file_list)

    # open the file
    with open(resistivity_file, "r") as file_content:
	    file_list = file_content.readlines()

    for i in range(2, len(file_list)):
        print(file_list[i])

    # split the file
    if file_list[1].strip() == "1":
        algorithm_choice = 1
    elif file_list[1].strip() == "2":
        algorithm_choice = 2

    for i in range(2, len(file_list)):
        fields = file_list[i].split()
        spacing_val = fields[0].strip()
        # print('-->' + spacing_val + '<--')
        resis_val = fields[1].strip()
        # print('-->' + resis_val + '<--')
        spacing_float = float(spacing_val)
        resis_float = float(resis_val)
        adat[i-2] = spacing_float
        rdat[i-2] = resis_float

    ndat = len(file_list) - 2

    # old readData()
    for i in range(1, ndat, 1):
        adatl[i] = np.log10(adat[i])
        rdatl[i] = np.log10(rdat[i])
    return

# def readFile():
#     global resistivity_file
#     global file_content

#     # read file
#     # file_content = resistivity_file.read()

#     file_content = open(resistivity_file, 'r')
#     return

# def readData():
#     # normally this is where the data would be read from the csv file
#     # but now I'm just hard coding it in as global lists
#     for i in range(1, ndat, 1):
#         adatl[i] = np.log10(adat[i])
#         print(adatl[i])
#         rdatl[i] = np.log10(rdat[i])
#         print(rdatl[i])
#     return

def layerDetails():
    # thick_max_layer = IntVar(mainwindow)
    # res_min_layer = IntVar(mainwindow)
    # res_max_layer = IntVar(mainwindow)

    global layers_choice
    global thick_min_layer
    global thick_max_layer
    global res_min_layer
    global res_max_layer

    # set number of layers
    layers_choice = num_layers.get()

    # initialize arrays
    thick_min_layer = [0]*(layers_choice - 1)
    thick_max_layer = [0]*(layers_choice - 1)
    res_min_layer = [0]*(layers_choice)
    res_max_layer = [0]*(layers_choice)

    # full frame
    layerinputframe = Frame(mainwindow, background="gainsboro")
    layerinputframe.pack(side=BOTTOM, anchor=SW)

    # note while working on functionality
    note_label = Label(layerinputframe, bg="gainsboro", font=("TkDefaultFont", 7),
                  text="  --> For predictable results, enter data commented out that was originally hard coded in.")
    note_label.grid(row=8, column=1, columnspan=3, pady=5)

    # thickness range
    thickness_label = Label(layerinputframe, bg="gainsboro", font=("TkDefaultFont", 13),
                  text="Model Range in Thickness (m)")
    thickness_label.grid(row=1, column=1, columnspan=3, pady=5)

    # resistivity range
    resistivity = Label(layerinputframe, bg="gainsboro", font=("TkDefaultFont", 13),
                  text="Model Range in Resistivity (m)")
    resistivity.grid(row=1, column=5, columnspan=3, pady=5)

    # thickness minimum values
    thick_min_label = Label(layerinputframe, bg="gainsboro",
                  text="Minimum\nValue", width=15)
    thick_min_label.grid(row=2, column=1)
    for i in range(0, layers_choice - 1, 1):
        thick_min_layer[i] = IntVar(mainwindow)
        thick_min_data = Entry(layerinputframe, textvariable=thick_min_layer[i], width=10)
        thick_min_data.grid(row=i+3, column=1)
    infinite_label = Label(layerinputframe, bg="gainsboro",
                  text="Infinite Thickness")
    infinite_label.grid(row=layers_choice+2, column=1, columnspan=2)

    # thickness maximum values
    thick_max_label = Label(layerinputframe, bg="gainsboro",
                  text="Maximum\nValue", width=15)
    thick_max_label.grid(row=2, column=2)
    for i in range(0, layers_choice - 1, 1):
        thick_max_layer[i] = IntVar(mainwindow)
        thick_max_data = Entry(layerinputframe, textvariable=thick_max_layer[i], width=10)
        thick_max_data.grid(row=i+3, column=2)

    # thickness predictions
    thick_pred_label = Label(layerinputframe, bg="gainsboro",
                  text="Thickness\nPrediction", width=15)
    thick_pred_label.grid(row=2, column=3)

    # resistivity predictions
    res_pred_label = Label(layerinputframe, bg="gainsboro",
                  text="Resistivity\nPrediction", width=15)
    res_pred_label.grid(row=2, column=5)

    # resistivity minimum values
    res_min_label = Label(layerinputframe, bg="gainsboro",
                  text="Minimum\nValue", width=15)
    res_min_label.grid(row=2, column=6)
    for i in range(0, layers_choice, 1):
        res_min_layer[i] = IntVar(mainwindow)
        res_min_data = Entry(layerinputframe, textvariable=res_min_layer[i], width=10)
        res_min_data.grid(row=i+3, column=6)

    # resistivity maximum values
    res_max_label = Label(layerinputframe, bg="gainsboro",
                  text="Maximum\nValue", width=15)
    res_max_label.grid(row=2, column=7)
    res_max_list = [0] * layers_choice
    for i in range(0, layers_choice, 1):
        res_max_layer[i] = IntVar(mainwindow)
        res_max_data = Entry(layerinputframe, textvariable=res_max_layer[i], width=10)
        res_max_data.grid(row=i+3, column=7)

    return

def executeVES():
    return



def error():
    sumerror = 0.
    #pltanswer = [0] * 64
    spline(m, one30, one30, asavl, rl, y2)
    for i in range(1, ndat, 1):
        ans = splint(m, adatl[i], asavl, rl, y2)
        sumerror = sumerror + (rdatl[i] - ans) * (rdatl[i] - ans)
        #print(i, sum1, rdat[i], rdatl[i], ans)
        pltanswerl[i] = ans
        pltanswer[i] = np.power(10, ans)
    rms = np.sqrt(sumerror / (ndat-1))

    # check the spline routine
#    for i in range(1, m+1, 1):
#        anstest = splint(m, asavl[i], asavl, rl, y2)
#        print( asavl[i], rl[i], anstest)
    #print(' rms  =  ', rms)
# if you erally want to get a good idea of all perdictions from Montecarlo
# perform the following plot (caution - change iter to a smaller number)
    #plt.loglog(adat[1:ndat], pltanswer[1:ndat])

    return rms


def transf(y, i):
    u = 1. / np.exp(y)
    t[1] = p[n]
    for j in range(2, layers_choice+1, 1):
        pwr = -2. * u * p[layers_choice+1-j]
        if pwr < np.log(2. * ep):
            pwr = np.log(2. * ep)
        a = np.exp(pwr)
        b = (1. - a)/(1. + a)
        rs = p[n + 1 - j]
        tpr = b * rs
        t[j] = (tpr + t[j-1]) / (1. + tpr * t[j-1] / (rs * rs))
    r[i] = t[layers_choice]
    return


def filters(b, k):
    for i in range(1, m+1, 1):
        re = 0.
        for j in range(1, k+1, 1):
            re = re + b[j] * r[i + k - j]
        r[i] = re
    return


def rmsfit():
    if algorithm_choice == SHCHLUMBERGER:
        y = spac - 19. * delx - 0.13069
        mum1 = m + 28
        for i in range(1, mum1 + 1, 1):
            transf(y, i)
            y = y + delx
        filters(FLTR1, 29)
    elif algorithm_choice == WENNER:
        s = np.log(2.)
        y = spac-10.8792495 * delx
        mum2 = m + 33
        for i in range(1, mum2 + 1, 1):
            transf(y, i)
            a = r[i]
            y1 = y + s
            transf(y1, i)
            r[i] = 2. * a - r[i]
            y = y + delx
        filters(FLTR2, 34)
    else:
        print(" type of survey not indicated")
        sys.exit()

    x = spac
    #print("A-Spacing   App. Resistivity")
    for i in range(1, m+1, 1):
        a = np.exp(x)
        asav[i] = a
        asavl[i] = np.log10(a)
        rl[i] = np.log10(r[i])
        x = x + delx
        #print("%7.2f   %9.3f " % ( asav[i], r[i]))

    rms = error()

    return rms

# my code to do a spline fit to predicted data at the nice spacing of Ghosh
# use splint to determine the spline interpolated prediction at the
# spacing where the measured resistivity was taken - to compare observation
# to prediction


def spline(n, yp1, ypn, x=[], y=[], y2=[]):
    u = [0] * 1000
    one29 = 0.99e30
    #print(x, y)
    if yp1 > one29:
        y2[0] = 0.
        u[0] = 0.
    else:
        y2[0] = -0.5
        u[0] = (3. / (x[1] - x[0])) * ((y[1] - y[0]) /
                                       (x[1] - x[0]) - yp1)

    for i in range(1, n):
        #print(i, x[i])
        sig = (x[i] - x[i - 1]) / (x[i + 1] - x[i - 1])
        p = sig * y2[i - 1] + 2.
        y2[i] = (sig - 1.) / p
        u[i] = ((6. * ((y[i + 1] - y[i]) / (x[i + 1] - x[i]) -
                       (y[i] - y[i - 1]) / (x[i] - x[i - 1])) /
                 (x[i + 1] - x[i - 1]) - sig * u[i - 1]) / p)

    if ypn > one29:
        qn = 0.
        un = 0.
    else:
        qn = 0.5
        un = (3. / (x[n] - x[n - 1])) * \
            (ypn - (y[n] - y[n-1]) / (x[n] - x[n - 1]))

    y2[n] = (un - qn * u[n - 1]) / (qn * y2[n - 1] + 1.)
    for k in range(n - 1, -1, -1):
        y2[k] = y2[k] * y2[k + 1] + u[k]

    return


def splint(n, x, xa=[], ya=[], y2a=[]):
    klo = 0
    khi = n
    while khi - klo > 1:
        k = int((khi + klo) // 2)
        if xa[k] > x:
            khi = k
        else:
            klo = k
    h = xa[khi] - xa[klo]
    if abs(h) < 1e-20:
        print(" bad xa input")
    #print(x, xa[khi], xa[klo])
    a = (xa[khi] - x) / h
    b = (x - xa[klo]) / h
    y = (a * ya[klo] + b * ya[khi] + ((a * a * a - a) * y2a[klo] +
                                      (b * b * b - b) * y2a[khi]) * (h * h) / 6.)
    #print("x=   ", x,"y=  ", y, "  ya=  ", ya[khi],"  y2a=  ", y2a[khi], "  h=  ",h)

    return y


# Main
if __name__ == '__main__':
    # Seed the RNG so we don't have randomness while testing
    random.seed(0)
    openGUI()
    mainwindow.mainloop()

    # set algorithm to Shchlumberger or Wenner
    # algorithm_choice = algorithm_index

    # set number of iterations
    iter = num_iter.get()

    # get number of layers
    layers_choice = num_layers.get()
    n = 2 * layers_choice - 1

    # set small[] and xlarge[]
    for i in range (layers_choice - 1):
        small[i + 1] = thick_min_layer[i].get()
    for i in range (layers_choice - 1):
        xlarge[i + 1] = thick_max_layer[i].get()
    for i in range (layers_choice):
        small[i + layers_choice] = res_min_layer[i].get()
    for i in range (layers_choice):
        xlarge[i + layers_choice] = res_max_layer[i].get()

    # set number of data points
    # ndat = num_datapoints.get()

    # readData()

    print(adat[1:ndat], rdat[1:ndat])
    for iloop in range(1, iter + 1, 1):
        #print( '  iloop is ', iloop)
        for i in range(1, layers_choice + 1, 1):
            randNumber = random.random()
            #print(randNumber, '  random')
            p[i] = (xlarge[i] - small[i]) * randNumber + small[i]
        for i in range(layers_choice + 1, n + 1, 1):
            randNumber = random.random()
            #print(randNumber, '  random')
            p[i] = (xlarge[i] - small[i]) * randNumber + small[i]

        rms = rmsfit()

        if rms < errmin:
            print('rms  ', rms, '   errmin ', errmin)
            for i in range(1, n+1, 1):
                pkeep[i] = p[i]
            for i in range(1, m+1, 1):
                rkeep[i] = r[i]
                rkeepl[i] = rl[i]
            for i in range(1, ndat + 1, 1):
                pltanswerkeepl[i] = pltanswerl[i]
                pltanswerkeep[i] = pltanswer[i]
            errmin = rms

# output the best fitting earth model
    print(' Layer ', '     Thickness  ', '   Res_ohm-m  ')
    for i in range(1, layers_choice, 1):
        print(i, pkeep[i], pkeep[layers_choice + i - 1])

    print(layers_choice, '  Infinite ', pkeep[n])
    for i in range(1, m + 1, 1):
        asavl[i] = np.log10(asav[i])

# output the error of fit
    print(' RMS error   ', errmin)
    print('  Spacing', '  Res_pred  ', ' Log10_spacing  ', ' Log10_Res_pred ')
    for i in range(1, m + 1, 1):
        #print(asav[i], rkeep[i], asavl[i], rkeepl[i])
        print("%7.2f   %9.3f  %9.3f  %9.3f" % (asav[i], rkeep[i],
                                               asavl[i], rkeepl[i]))

    plt.loglog(asav[1:m], rkeep[1:m], '-')  # resistivity prediction curve
    plt.loglog(adat[1:ndat], pltanswerkeep[1:ndat],
               'ro')  # predicted data red dots
    s = 7
    plt.loglog(adat[1:ndat], rdat[1:ndat], 'bo',
               markersize=s)  # original data blue dots
    plt.show()
    plt.grid(True)
    sys.exit(0)
