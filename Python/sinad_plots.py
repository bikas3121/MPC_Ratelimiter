# %%
import numpy as np
from matplotlib import pyplot as plt
from simulation_plots import line_plot
# %% Non_opt filter
S_NO_I = [57.4, 57.8, 22.2, 43.7, 43.5, 42.8]
S_NO_INL = [43.5, 43.5, 22.2, 42.9, 43.3, 43.1]

S_O_I = [57.4, 60.9, 60.9, 60.8, 61.0, 61.2]
S_O_INL = [43.5, 56.2, 56.2, 56.3, 56.3, 56.5]

lin_methods = np.array(["DQ", "NSQ", "MHOQ  \n(N=1)", "MHOQ \n(N=2)", "MHOQ \n(N=3)", "MHOQ \n(N=4)"])
x = np.linspace(1,len(lin_methods),len(lin_methods))


# %% plot line plots
line_plot(lin_methods, S_NO_I, S_O_I, 79.9, descr='SINAD')
line_plot(lin_methods, S_NO_INL, S_O_INL, 79.9, descr='SINAD')

# %% SINAD plots for ideal case
# fig, ax = plt.subplots()
# ax.plot(x, S_NO_I, linestyle = "-", marker = "o", label= "$F_{B}[z]$")
# ax.plot(x, S_O_I, linestyle = "-", marker = "o", label= "$F_{B}^{*}[z]$")
# ax.axhline(y = 79.9, color = 'r', linestyle= '--')
# ax.text(3, 82, "$\mathrm{SINAD}_{th}$ = 79.9 dB")
# ax.set_xticks(x)
# ax.set_xticklabels(lin_methods)
# ax.set_xlabel("Linearisation methods")
# ax.set_ylabel("SINAD")
# ax.grid( visible = True, which ='both', axis = 'y', linestyle =':')
# ax.legend()
# ax.set_ylim(20, 90)
# # %% SINAD plots for non-ideal case