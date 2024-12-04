#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Gains with some moderate pains
# Plots the variance plots
@author: Bikash Adhikari
@date: 22.02.2024
@license: BSD 3-Clause
"""

# %%
import numpy as np
from matplotlib import pyplot as plt

# %% 
Nb = 12
# %% ENOB save (normal analog butterworth filter in the output)
# ENOBS non-optimal
ENOB10_NOPT_DIRECT = 11.069
ENOB10_NOPT_NSD = 11.069
ENOB10_NOPT_MPC_1 = 6.781 # (N = 1)
ENOB10_NOPT_MPC_2 = 16.060 # (N = 2)

# ENOBS optimal
ENOB10_OPT_DIRECT = 11.069
ENOB10_OPT_NSD = 16.654
ENOB10_OPT_MPC_1 = 16.654 # (N = 1)
ENOB10_OPT_MPC_2 = 16.730 # (N = 2)


# Butterworth filter with 100k cutoff frequency
# ENOBS Non- optimal
ENOB100_NOPT_DIRECT = 10.824
ENOB100_NOPT_NSD = 10.832
ENOB100_NOPT_MPC_1 = 4.166 # (N = 1)
ENOB100_NOPT_MPC_2 = 7.988 # (N = 2)


# ENOBS optimal
ENOB100_OPT_DIRECT = 10.824
ENOB100_OPT_NSD = 12.938
ENOB100_OPT_MPC_1 = 12.938 # (N = 1)
ENOB100_OPT_MPC_2 = 12.940 # (N = 2)



ENOB_10_NOPT = np.array([ENOB10_NOPT_DIRECT, ENOB10_NOPT_NSD, ENOB10_NOPT_MPC_1, ENOB10_NOPT_MPC_2])
ENOB_10_OPT = np.array([ENOB10_OPT_DIRECT, ENOB10_OPT_NSD, ENOB10_OPT_MPC_1, ENOB10_OPT_MPC_2])


ENOB_100_NOPT = np.array([ENOB100_NOPT_DIRECT, ENOB100_NOPT_NSD, ENOB100_NOPT_MPC_1, ENOB100_NOPT_MPC_2])
ENOB_100_OPT = np.array([ENOB100_OPT_DIRECT, ENOB100_OPT_NSD, ENOB100_OPT_MPC_1, ENOB100_OPT_MPC_2])

lin_methods = ['Direct', 'NSD', 'MHOQ (N=1)', 'MHOQ (N=2)']#, 'ILC']
barWidth = 0.25
# set position of the bar on X axis
bar1 = np.arange(ENOB_10_OPT.size)
bar2 = [x + barWidth for x in bar1] 

# % Draw plot
fig, ax = plt.subplots(figsize = (7,5))
plt.grid()
# plt.axhline(y = 0, color = 'black', linestyle = '-')
b1 = plt.bar(bar1, ENOB_100_NOPT, width = barWidth, color = 'tab:blue', edgecolor = 'white', label = 'Static')
b2 = plt.bar(bar2, ENOB_100_OPT, width = barWidth, color = 'tab:orange', edgecolor = 'white', label = 'Static')

plt.xlabel('Linearisation method', fontweight ='bold', fontsize = 15) 
plt.ylabel('ENOB', fontsize = 13) 

pos_xticks = np.array([r  for r in range(len(ENOB_10_OPT))]) + barWidth/2
plt.xticks(pos_xticks, lin_methods , fontsize = 10)
plt.legend (['Non-optimal(Fc = 100k)', 'Optimal(Fc = 100k)'])

ah = []
for rect in b1 + b2 :
    height = rect.get_height()
    ah.append(height)
    if height > 0 :
        plt.text(rect.get_x() + rect.get_width()/2.0 - 0.03, 0.3, f'{height:.2f} bits', rotation = 90, fontsize  = 10)        
    if height < 0 :
        plt.text(rect.get_x() + rect.get_width()/2.0 - 0.03, 0.3, f'{height:.2f} bits', rotation = 90, fontsize  = 10)        

ax.set_axisbelow(True)
ax.grid(zorder=0, axis = "y")
plt.show()
# %% DAC  and sampling frequency
DAC_16bit_Fs_1MHz = 1
DAC_6bit_Fs_1MHz = 2
DAC_6bit_Fs_250MHz = 3
DAC_16bit_Fs_32MHz = 4
DAC_6bit_Fs_32MHz = 5


# %%
match 1:
    case 1: # 16-bit 1 MHz Trond

        methods ="static-spice"
        Nb = 16
        Fs = 1.02
        static_baseline = 11.121
        static_physcal = 16.487
        static_phfd = 11.283
        static_shpd = 10.750 
        static_nsdcal = 14.416
        static_dem = 8.645
        #static_ilc= 15.551

        # Spice simulation results
        spice_baseline = 11.111
        spice_physcal = 15.189
        spice_phfd = 11.173
        spice_shpd = 10.672
        spice_nsdcal = 14.249
        spice_dem = 8.640
        #spice_ilc= 15.004

        # Calculate gains
        #static_gains = np.array([static_physcal,  static_phfd, static_shpd, static_dem, static_nsdcal,static_ilc ])- static_baseline
        #spice_gains = np.array([spice_physcal,  spice_phfd, spice_shpd, spice_dem, spice_nsdcal,spice_ilc])- spice_baseline
        
        static_gains = np.array([static_physcal,  static_phfd, static_shpd, static_nsdcal, static_dem]) - static_baseline
        spice_gains = np.array([spice_physcal,  spice_phfd, spice_shpd, spice_nsdcal, spice_dem]) - spice_baseline

        lin_methods = ['PHYSCAL', 'PHFD', 'SHPD', 'NSDCAL', 'DEM']#, 'ILC']


# %% Plots
barWidth = 0.25
# set position of the bar on X axis
bar1 = np.arange(static_gains.size)
bar2 = [x + barWidth for x in bar1] 
bar3 = [x + barWidth for x in bar2]

# % Draw plot
fig, ax = plt.subplots(figsize = (7,5))
plt.axhline(y = 0, color = 'black', linestyle = '-')
b1 = plt.bar(bar2, static_gains, width = barWidth, color = 'tab:blue', edgecolor = 'white', label = 'Static')
match methods:
    case "static-spice":
        b2 = plt.bar(bar3, spice_gains, width = barWidth,  color = 'tab:orange',edgecolor = 'white', label = 'SPICE')
    case "static-spectre":
        b2 = plt.bar(bar3, spectre_gains, width = barWidth,  color = 'tab:orange',edgecolor = 'white', label = 'Spectre')

plt.xlabel('Linearisation method', fontweight ='bold', fontsize = 15) 
plt.ylabel('ENOB gain', fontsize = 13) 

pos_xticks = np.array([r + barWidth for r in range(len(static_gains))]) + barWidth/2
plt.xticks(pos_xticks, lin_methods , fontsize = 13)

ah = []
for rect in b1 + b2 :
    height = rect.get_height()
    ah.append(height)
    if height > 0 :
        plt.text(rect.get_x() + rect.get_width()/2.0 - 0.03, 0.3, f'{height:.2f} bits', rotation = 90, fontsize  = 13)        
    if height < 0 :
        plt.text(rect.get_x() + rect.get_width()/2.0 - 0.03, 0.3, f'{height:.2f} bits', rotation = 90, fontsize  = 13)        

plt.title(f"{int(Nb)}-bit DAC with sampling rate Fs = {Fs} MHz", fontsize = "13")
plt.legend(fontsize="13", loc='upper right')
ax.set_axisbelow(True)
ax.grid(zorder=0, axis = "y")
fig.tight_layout()
# plt.savefig(f"Gainplot-{Nb}bits.pdf")

# %%
# fname = "figures/Gainplot-" + str(int(Nb)) + str("bits") + str(int(Fs)) + str("MHz-") + methods
# fname = str(fname) + ".pdf"
# fig.savefig(fname, format='pdf', bbox_inches='tight')
