
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Gains with some moderate pains

@author: Bikash Adhikari
@date: 22.02.2024
@license: BSD 3-Clause
"""

# %%
import numpy as np
from matplotlib import pyplot as plt

# %% 
Nb = 12

# ENOB Perception filter vs second order butter filter
ENOB_P_DIR = 12.344 
ENOB_P_NSD = 14.853
ENOB_P_MPC_1 = 14.853
ENOB_P_MPC_2 = 14.90


ENOB10_B_DIR = 11.069 
ENOB10_B_NSD = 11.069
ENOB10_B_MPC_1 = 6.781
ENOB10_B_MPC_2 = 16.060 

ENOB100_B_DIR = 10.824
ENOB100_B_NSD = 10.832
ENOB100_B_MPC_1 = 4.166 # (N = 1)
ENOB100_B_MPC_2 = 7.988 # (N = 2)

ENOB_P = np.array([ENOB_P_DIR , ENOB_P_NSD, ENOB_P_MPC_1 , ENOB_P_MPC_2 ])
ENOB_B10 = np.array([ENOB10_B_DIR , ENOB10_B_NSD, ENOB10_B_MPC_1 , ENOB10_B_MPC_2 ])
ENOB_B100 = np.array([ENOB100_B_DIR , ENOB100_B_NSD, ENOB100_B_MPC_1 , ENOB100_B_MPC_2 ])

# Butterworth filter with 10k cutoff frequency
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
bar3 = [x + barWidth for x in bar2]
# % Draw plot
fig, ax = plt.subplots(figsize = (7,5))
plt.grid()
# plt.axhline(y = 0, color = 'black', linestyle = '-')
b1 = plt.bar(bar1, ENOB_P, width = barWidth, color = 'tab:blue', edgecolor = 'white', label = 'Static')
b2 = plt.bar(bar2, ENOB_B10, width = barWidth, color = 'tab:orange', edgecolor = 'white', label = 'Static')
b3 = plt.bar(bar3, ENOB_B100, width = barWidth, color = 'tab:red', edgecolor = 'white', label = 'Static')

plt.xlabel('Linearisation method', fontweight ='bold', fontsize = 15) 
plt.ylabel('ENOB', fontsize = 13) 

pos_xticks = np.array([r  for r in range(len(ENOB_10_OPT))]) + barWidth/2 
plt.xticks(pos_xticks, lin_methods , fontsize = 10)
plt.legend (['Psychoacoustic Filter', 'Butterworth Filter (10k)', 'Butterworth Filter (100k)'])

ah = []
for rect in b1 + b2 + b3 :
    height = rect.get_height()
    ah.append(height)
    if height > 0 :
        plt.text(rect.get_x() + rect.get_width()/2.0 - 0.03, 0.3, f'{height:.2f} bits', rotation = 90, fontsize  = 10)        
    if height < 0 :
        plt.text(rect.get_x() + rect.get_width()/2.0 - 0.03, 0.3, f'{height:.2f} bits', rotation = 90, fontsize  = 10)        

ax.set_axisbelow(True)
ax.grid(zorder=0, axis = "y")
plt.show()

