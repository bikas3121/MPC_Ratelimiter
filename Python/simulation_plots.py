
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""  Plot the  variance plots for the given methods 

@author: Bikash Adhikari
@date: 24.03.2024
@license: BSD 3-Clause
"""

 # %%
# imports
import numpy as np
import matplotlib.pyplot as plt



# %%  Bar plots for the Error Variance, SINAD or ENOB
def bar_plot(descr = '', **kwargs):
    # Plots the bar plots for the given values 
    x = []
    y = []
    y1 = []
    for key, value in kwargs.items() :
        x.append(key)
        y.append(np.round(value,5))
        y1.append(value)

    fig, ax = plt.subplots()
    x_pos = np.arange(0, len(x),1) 
    bar_labels = x
    ax.bar(x,y, label=bar_labels)
    plt.xticks([i  for i in range(len(y))],x[0:len(x)])
    for i in range(len(x)):
        # plt.text(x = x_pos[i]-0.1, y= y[i] + 0.015 , s= y[i], size = 10)
        plt.text(x = x_pos[i]-0.1, y= y[i] + 0.015*y[i] , s= y[i], size = 10)
    ax.set_ylabel(descr)
    ax.legend()
    plt.show()


    for i in range(len(x)):
         print(x[i] + ' ' + descr + ':{}'.format(y1[i]))

# %% LInear plots
def line_plot(lin_methods, values1, values2, SINAD_th, descr = ''):

    x = np.linspace(1,len(lin_methods),len(lin_methods))

    # theoretical sinad value 
    txt1 = "SINAD_TH = " 
    txt2 = str(SINAD_th)
    txt3 = " dB"
    txt = txt1 + txt2 + txt3

    # Plots the bar plots for the given values 
    fig, ax = plt.subplots()
    ax.plot(x, values1, linestyle = "-", marker = "o", label= "$H[z]$")
    ax.plot(x, values2, linestyle = "-", marker = "o", label= "$H^{*}[z]$")
    ax.axhline(y = 79.9, color = 'r', linestyle= '--')
    # ax.text(3, 82, "$\mathrm{SINAD}_{th}$ =  79.9 dB")
    ax.text(3, 82, txt)
    ax.set_xticks(x)
    ax.set_xticklabels(lin_methods)
    ax.set_xlabel("Linearisation methods")
    ax.set_ylabel(descr)
    ax.grid( visible = True, which ='both', axis = 'y', linestyle =':')
    ax.legend()
    ax.set_ylim(20, 90) 
